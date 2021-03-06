#' Data augmentation with a pseudo-marginal Markov chain Monte Carlo algorithm for multivariate generalized Pareto models
#'
#'
#' @param dat n by D matrix of observations
#' @param mthresh vector of marginal thresholds under which data are censored
#' @param thresh functional max threshold determining the risk region
#' @param likt string indicating the type of likelihood, with an additional contribution for the non-exceeding components: one of  \code{"mgp"}, \code{"binom"} and \code{"pois"}.
#' @param coord matrix of coordinates, with longitude and latitude in the first two columns and additional covariates for the latent Gaussian model
#' @param blockupsize size of block for updates of the scale parameter; \code{ncol(dat)} yields individual updates
#' @param start named list with starting values for the parameters, with arguments:
#' \itemize{
#' \item\code{scale}: a \code{D} vector of scale parameters, strictly positive.
#' \item\code{shape}: a scale containing the shape and satisfying support constraints for all sites.
#' \item\code{marg.pcov}: initial proposal covariance for the marginal parameters (scale and shape).
#' \item\code{dep}: initial values for the dependence function.
#' \item\code{dep}: dependence function, with distance matrix as first argument and \code{dep} as second argument.
#' \item\code{aniso}: initial values for the anisotropy parameters, scale and angle.
#' \item\code{dep.lb}: lower bounds for the dependence parameters
#' \item\code{dep.ub}: upper bounds for the dependence parameters
#' }
#' If any of \code{scale}, \code{shape} or \code{marg.pcov} are missing, the function will attempt to find starting values.
#' @param geoaniso logical; should geometric anisotropy be included? Default to \code{TRUE}.
#' @param lambdau probability of exceedance of the threshold for censored observations
#' @param numiter number of iterations to be returned
#' @param thin thining parameter; only every \code{thin} iteration is saved
#' @param burnin number of initial parameters for adaptation and discarded values.
#' @param verbose report current values via print every \code{verbose} iterations.
#' @param filename name of file for save.
#' @param keepburnin logical; should initial runs during \code{burnin} be kept for diagnostic. Default to \code{TRUE}.
#' @param saveinterm integer indicating when to save results. Default to \code{500L}.
#' @inheritDotParams clikmgp
#' @export
#' @keywords internal
#' @return a list with \code{res} containing the results of the chain
da.mgp <- function(dat, mthresh, thresh, lambdau = 1, coord, start,
                   numiter = 4e4L, burnin = 5e3L, thin = 1L, verbose = 50L, filename,
                   keepburnin = TRUE, geoaniso = TRUE, blockupsize = ncol(dat), likt = c("mgp", "pois", "binom"), transform = FALSE,
                   saveinterm = 500L, ...) {
  slurm_arrayid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
  slurm_jobid <- Sys.getenv("SLURM_JOB_ID")
  filename <- paste0(filename, ifelse(slurm_arrayid == "", "", "_"), slurm_arrayid, "_da", ifelse(slurm_jobid == "", "", "_"), slurm_jobid)

  likt <- match.arg(likt)
  ellips <- list(...)
  B <- numiter * thin + burnin
  if (isTRUE(is.finite(saveinterm))) {
    if (saveinterm < 0) {
      saveinterm <- B + 1L
    } else {
      saveinterm <- as.integer(saveinterm)
    }
  } else {
    saveinterm <- B + 1L
  }
  n <- nrow(dat)
  D <- ncol(dat)
  mthresh <- rep(mthresh, length.out = D)
  dat <- t(t(dat) - mthresh)
  mthresh <- rep(0, D)
  lambdau <- rep(lambdau, length.out = D)
  # Censoring
  censored <- t(t(dat) < mthresh)
  numAbovePerRow <- D - rowSums(censored)

  # Remove observations that are not extreme
  zeros <- which(numAbovePerRow == 0)
  if (length(zeros) > 0) {
    dat <- dat[-zeros, ]
    n <- nrow(dat)
    censored <- t(t(dat) < mthresh)
    numAbovePerRow <- D - rowSums(censored)
  }
  numAbovePerCol <- n - colSums(censored)


  # Coordinates for distance and regression covariates for spatial random effect model on scale
  loc <- as.matrix(coord[, 1:2])
  Xm <- cbind(1, scale(as.matrix(coord)))
  di <- distg(loc = coord[, 1:2], scale = 1, rho = 0)

  ellips <- list(...)
  # Precision and generating vectors for Monte Carlo routines
  genvec1 <- ellips$genvec1
  genvec2 <- ellips$genvec2
  B1 <- ifelse(is.null(ellips$B1), 499L, ellips$B1)
  B2 <- ifelse(is.null(ellips$B2), 307L, ellips$B2)
  if (is.null(genvec1)) {
    genvec1 <- mvPot::genVecQMC(B1, ncol(dat) - 1L)$genVec
  }
  if (is.null(genvec2)) {
    genvec2 <- mvPot::genVecQMC(B2, ncol(dat) - 1L)$genVec
  }
  ncores <- ifelse(is.null(ellips$ncores), 1L, ellips$ncores)
  mmin <- rep(0, D)
  mmax <- apply(dat, 2, max)
  # Independence likelihood
  indeplik <- function(par, dat, mthresh, ...) {
    stopifnot(ncol(dat) == length(par) - 1L)
    D <- length(par) - 1
    -sum(sapply(1:D, function(j) {
      mev::gpd.ll(par = c(par[j], par[D + 1]), dat = dat[dat[, j] > mthresh[j], j] - mthresh[j])
    }))
  }

  # Starting values for the chain and proposal covariance for the marginal parameters
  scale.c <- start$scale
  shape.c <- start$shape
  if (is.null(shape.c) || length(shape.c) == 1L) {
    cshape <- TRUE
  } else {
    cshape <- FALSE
  }
  marg.pcov <- start$marg.pcov
  if (is.null(marg.pcov) || is.null(scale.c) || is.null(shape.c)) {
    # Starting values for the parameters
    # Independence likelihood for marginal parameters with GP likelihood pointwise
    margpars <- suppressWarnings(sapply(1:D, function(j) {
      mev::fit.gpd(xdat = as.vector(dat[, j]), threshold = mthresh[j])$est
    }))
    mean0 <- c(margpars[1, ], mean(margpars[2, ]))
  } else {
    mean0 <- c(scale.c, shape.c)
  }
  # These may not be feasible, but auglag can take the argument nevertheless as starting value
  opt.ind <- try(alabama::auglag(
    par = mean0, fn = indeplik, hin = function(par, dat, mthresh, ...) {
      c(par[-length(par)], par[-length(par)] + par[length(par)] * (apply(dat, 2, max) - mthresh))
    },
    mthresh = mthresh, dat = dat, control.outer = list(trace = FALSE)
  ))

  if (!is.character(opt.ind) && isTRUE(opt.ind$kkt1 && opt.ind$kkt2)) {
    scale.c <- opt.ind$par[1:D]
    shape.c <- opt.ind$par[D + 1]
    marg.pcov <- solve(opt.ind$hessian)
    # Somehow the matrix returned by alabama is not symmetric...
    marg.pcov[lower.tri(marg.pcov)] <- t(marg.pcov)[lower.tri(marg.pcov)]
  } else {
    stop("Could not find default starting values for scale and shape parameters.")
  }

  # Copy dependence parameters
  dep.c <- start$dep
  # Constant for scaling of covariance matrix

  ndep <- length(dep.c)
  if (ndep == 0) {
    stop("Invalid dependence parameter `dep` in `start`")
  }
  dep.lb <- start$dep.lb
  dep.ub <- start$dep.ub
  stopifnot(length(dep.lb) == ndep, length(dep.ub) == ndep, dep.lb < dep.ub, dep.c >= dep.lb, dep.c <= dep.ub)
  dep.fun <- start$dep.fun
  stopifnot(is.function(dep.fun))
  dep.pcov <- start$dep.pcov
  if (is.null(dep.pcov)) {
    dep.pcov <- diag(ndep)
  }
  # log-prior function for the dependence parameters
  dep.lpriorfn <- start$dep.lpriorfn
  stopifnot(is.function(dep.lpriorfn))
  if (is.null(dep.pcov)) {
    stop("Missing `dep.prior` in `start`")
  }


  # geometric anisotropy
  if (geoaniso) {
    aniso.c <- start$aniso
    aniso.pcov <- diag(c(0.01, 0.01))
    if (is.null(aniso.c)) {
      aniso.c <- c(1.2, 0)
    }
    aniso.lb <- c(1, -Inf)
  }

  # set number of parameters, minus regression parameters (separate containers appended via cbind to res)
  npar <- length(scale.c) + length(shape.c) + ndep + ifelse(geoaniso, 2, 0)
  # Positions of parameters for saving
  scale.i <- 1:D
  shape.i <- D + 1
  dep.i <- (D + 2):(D + 1 + length(dep.c))
  if (geoaniso) {
    aniso.i <- (npar - 1):npar
  }

  lscalelm.i <- (npar + 1):(npar + ncol(Xm) + 2)
  npar <- npar + ncol(Xm) + 2
  if (transform) {
    transform.fn <- function(x) {
      transfo.pars(x, lbound = dep.lb, ubound = dep.ub)
    }
  } else {
    transform.fn <- identity
  }

  # Hyperpriors for Bayesian linear model
  ols <- lm(log(scale.c) ~ -1 + Xm)
  lscale.hyp.mean.c <- ols$coef
  if (!requireNamespace("geoR", quietly = TRUE)) {
    warning("`geoR` package is not installed.")
    lscale.hyp.tausq.c <- sum(resid(ols)^2) / ols$df.residual
    lscale.hyp.rho.c <- max(di[lower.tri(di)]) / 100
  } else {
    lan.vcov.o <- geoR::likfit(
      coords = loc, trend = geoR::trend.spatial(~ -1 + Xm),
      data = log(scale.c), cov.model = "exponential",
      ini.cov.pars = rep(1, 2), lik.method = "ML", messages = FALSE
    )
    lscale.hyp.rho.c <- lan.vcov.o$phi
    lscale.hyp.tausq.c <- lan.vcov.o$sigmasq
  }
  # Set hyperpriors for random effect on log scale
  lscale.hyp.precis.c <- solve(powerexp.cor(h = di, scale = lscale.hyp.rho.c))
  lscale.fhyp.mean.Vinv <- diag(ncol(Xm)) * 5
  lscale.fhyp.mean.b <- rep(0, ncol(Xm))
  lscale.fhyp.tausq.a <- 0.5
  lscale.fhyp.tausq.b <- 0.1
  lscale.fhyp.crossprod <- c(t(lscale.fhyp.mean.b) %*% lscale.fhyp.mean.Vinv %*% lscale.fhyp.mean.b)
  lscale.mu <- as.vector(Xm %*% lscale.hyp.mean.c)


  # Function to create list with parameter, input is model dependent
  makepar <- function(dep, distm, df = NULL) {
    Lambda <- dep.fun(distm, par = dep)
    return(list(Lambda = Lambda))
  }
  if (geoaniso) {
    distm.c <- distg(loc, scale = aniso.c[1], rho = aniso.c[2])
    aniso.lpriorfn <- function(aniso) {
      dnorm(aniso[1] - 1, sd = 1, mean = 0, log = TRUE)
      # dgamma(aniso[1] - 1, scale = 5, shape = 1, log = TRUE) #anisotropy
    }
    aniso.pcov <- diag(c(0.01, 0.01))
  } else { # pointer to distance matrix
    distm.c <- di
    aniso.pcov <- NULL
  }
  ntot <- ellips$ntot

  par.c <- makepar(dep = dep.c, distm = distm.c)



  imputefn <- function(scale, shape, par, ...) {
    impute(
      dat = dat, thresh = thresh, mthresh = mthresh, loc = rep(0, D),
      scale = scale, shape = shape, lambdau = lambdau, riskr = "max",
      par = par, map = FALSE, ...
    )
  }


  # Set block update size for scale parameters
  blockupsize <- min(max(1L, as.integer(blockupsize)), D)
  facs <- as.factor(rep(1:blockupsize, length = D))


  # Define containers
  lpost <- rep(0, B)
  accept <- rep(0L, npar)
  attempt <- rep(0L, npar)
  res <- matrix(data = 0, ncol = npar, nrow = B)

  # Start timer
  time.0 <- proc.time()[3] # Better would be to time every run and keep median


  ########################################################################
  ####################          LOOP         #############################
  ########################################################################
  thin <- as.integer(thin)
  for (b in 1:B) {
    idat <- imputefn(scale = scale.c, shape = shape.c, par = par.c)
    loglikfn <- function(scale, shape, par, ...) {
      # Initial evaluation of the log-likelihood
      likmgp(
        dat = idat, thresh = thresh, loc = rep(0, D), scale = scale,
        shape = shape, par = par, model = "br", likt = likt, lambdau = lambdau,
        B1 = B1, genvec1 = genvec1, ncores = ncores, ntot = ntot
      )
    }
    loglik.c <- loglikfn(scale = scale.c, shape = shape.c, par = par.c)

    attempt <- attempt + 1L

    ####          UPDATE MARGINAL SCALE PARAMETERS          ####
    # Metropolis within Gibbs with random scans
    # Update scale parameter
    order <- split(sample.int(D, size = D, replace = FALSE), f = facs)
    scale.lpriorfn <- function(margpar) {
      scale <- margpar[1:ncol(lscale.hyp.precis.c)]
      dmvnorm.precis(
        x = log(scale), mean = as.vector(Xm %*% lscale.hyp.mean.c),
        precis = lscale.hyp.precis.c / lscale.hyp.tausq.c, logd = TRUE
      ) - sum(log(scale))
    }
    scale.loglikfn <- function(scale) {
      loglikfn(scale = scale, shape = shape.c, par = par.c)
    }
    for (k in 1:length(order)) {
      scale.lb <- pmax(-shape.c * mmax[order[[k]]], 0)
      update <- mh.fun(
        cur = c(scale.c, shape.c), lb = scale.lb, ind = order[[k]], lik.fun = scale.loglikfn,
        ll = loglik.c, pcov = marg.pcov, cond = FALSE, prior.fun = scale.lpriorfn
      )
      # Increase acceptance count, update log
      if (update$accept) {
        accept[order[[k]]] <- accept[order[[k]]] + 1L
        loglik.c <- update$ll
        scale.c <- update$cur[scale.i]
      }
    }

    ####          UPDATE HYPERPRIORS ON SCALE               ####
    # Conjugate updates from NIG model b_sigma,
    # b_xi | tausq ~ No(0, tausq*Vb_xi)
    # tausq ~ IG(a, b)
    # Gibbs steps, Conjugate Priors in Bayesian Linear Model

    # Step 1: sample tausq from IG
    lscale.c <- log(scale.c)
    RinvX <- lscale.hyp.precis.c %*% Xm
    lscale.hyp.mean.precis <- crossprod(Xm, RinvX) + lscale.fhyp.mean.Vinv
    lscale.hyp.mean.mu <- c(crossprod(lscale.c, RinvX)) + c(lscale.fhyp.mean.Vinv %*% lscale.fhyp.mean.b)
    Mm <- solve(lscale.hyp.mean.precis) %*% lscale.hyp.mean.mu
    lscale.hyp.tausq.c <- 1 / rgamma(
      n = 1, shape = lscale.fhyp.tausq.a + 0.5 * (D + ncol(Xm)),
      rate = lscale.fhyp.tausq.b +
        0.5 * (lscale.fhyp.crossprod +
          c(t(lscale.c) %*% lscale.hyp.precis.c %*% lscale.c) -
          t(lscale.hyp.mean.mu) %*% Mm)
    )
    # Conditional on tausq, sample beta (mean) from Normal
    lscale.hyp.mean.c <- c(rmvnormprec(n = 1, mu = c(Mm), precis = lscale.hyp.mean.precis / lscale.hyp.tausq.c))

    # Laplace approximation for the range parameter
    lscale.hyp.rho.c <- updt.range(
      tau = lscale.c - c(Xm %*% lscale.hyp.mean.c), alpha = 1 / lscale.hyp.tausq.c,
      lambda = lscale.hyp.rho.c, di = di, a = 2, b = 2, discount = 0.4, lb = 1e-2, maxstep = 2
    )
    if (attributes(lscale.hyp.rho.c)$accept) { # only update is move is accepted
      lscale.hyp.precis.c <- solve(powerexp.cor(h = di, scale = lscale.hyp.rho.c))
    }

    ####          UPDATE MARGINAL SHAPE PARAMETER           ####

    shape.lb <- max(-0.49, -min(scale.c / mmax))
    shape.lpriorfn <- function(shape) {
      dnorm(x = shape, mean = 0, sd = 0.2, log = TRUE)
    }
    shape.loglikfn <- function(shape) {
      loglikfn(scale = scale.c, shape = shape, par = par.c)
    }
    update <- mh.fun(
      cur = shape.c, lb = shape.lb, ind = 1, lik.fun = shape.loglikfn, prior.fun = shape.lpriorfn,
      ll = loglik.c, pcov = as.matrix(marg.pcov[D + 1, D + 1]), cond = FALSE
    )
    # Increase acceptance count, update log
    if (update$accept) {
      accept[shape.i] <- accept[shape.i] + 1L
      loglik.c <- update$ll
      shape.c <- update$cur
    }


    ####           UPDATE DEPENDENCE PARAMETERS             ####
    # adding prior contribution for the anisotropy parameters


    dep.loglikfn <- function(dep) {
      par <- list(Lambda = dep.fun(distm.c, dep))
      if (isTRUE(any(is.nan(par$Lambda)))) { # If alpha too close to zero, invalid matrix
        ll <- -Inf
      } else {
        ll <- loglikfn(scale = scale.c, shape = shape.c, par = par)
      }
      attributes(ll)$par <- par
      ll
    }

    # Perform first updates parameter by parameter
    if (b < min(burnin, 2000L)) {
      for (i in 1:ndep) {
        update <- mh.fun(
          cur = dep.c, lb = dep.lb[i], ub = dep.ub[i], ind = i, lik.fun = dep.loglikfn,
          ll = loglik.c, pcov = dep.pcov, cond = TRUE, prior.fun = dep.lpriorfn, transform = transform
        )
        if (update$accept) {
          par.c <- attributes(update$ll)$par
          loglik.c <- update$ll
          dep.c <- update$cur
          accept[dep.i[i]] <- accept[dep.i[i]] + 1L
        }
      }
    } else {
      update <- mh.fun(
        cur = dep.c, lb = dep.lb, ub = dep.ub, ind = 1:ndep, lik.fun = dep.loglikfn,
        ll = loglik.c, pcov = dep.pcov, cond = FALSE, transform = transform, prior.fun = dep.lpriorfn
      )
      if (update$accept) {
        par.c <- attributes(update$ll)$par
        loglik.c <- update$ll
        dep.c <- update$cur
        accept[dep.i] <- accept[dep.i] + 1L
      }
    }

    if (geoaniso) {
      aniso.loglikfn <- function(aniso) {
        aniso[2] <- wrapAng(aniso[2])
        distm <- distg(loc = loc, scale = aniso[1], rho = aniso[2])
        Lambda <- dep.fun(distm, dep.c)
        par <- list(Lambda = Lambda)
        if (isTRUE(any(is.nan(par$Lambda)))) { # If alpha too close to zero, invalid matrix
          ll <- -Inf
        } else {
          ll <- loglikfn(scale = scale.c, shape = shape.c, par = par)
        }
        attributes(ll)$par <- par
        attributes(ll)$distm <- distm
        ll
      }

      update <- mh.fun(
        cur = aniso.c, lb = aniso.lb, ind = 1:2, lik.fun = aniso.loglikfn,
        ll = loglik.c, pcov = aniso.pcov, cond = FALSE, prior.fun = aniso.lpriorfn
      )
      if (update$accept) {
        accept[aniso.i] <- accept[aniso.i] + 1L
        par.c <- attributes(update$ll)$par
        distm.c <- attributes(update$ll)$distm
        loglik.c <- update$ll
        aniso.c <- update$cur
        aniso.c[2] <- wrapAng(aniso.c[2])
      }
    }

    # Save current value of the log-likelihood and of all parameters at the end of the iteration
    if (thin == 1 || b %% thin == 0) {
      i <- as.integer(b / thin)
      lpost[i] <- loglik.c
      res[i, scale.i] <- scale.c
      res[i, shape.i] <- shape.c
      res[i, dep.i] <- dep.c
      if (geoaniso) {
        res[i, aniso.i] <- aniso.c
      }
      res[i, lscalelm.i] <- c(lscale.hyp.mean.c, lscale.hyp.tausq.c, lscale.hyp.rho.c)
    }

    # Adapt covariance matrix, but using only previous iterations
    # Stop adapting after burnin
    if (b < burnin) {
      if (b %% 20 * thin == 0 && b > 200L && b <= 2000L) {
        # Update covariance matrix of the marginal proposals - only diagonal elements to begin with
        updiag <- sqrt(diag(marg.pcov))
        for (j in 1:(D + 1)) { # scale and shape parameters
          ada <- adaptive(attempts = attempt[j], acceptance = accept[j], sd.p = updiag[j])
          updiag[j] <- ada$sd / updiag[j]
          accept[j] <- ada$acc
          attempt[j] <- ada$att
        }
        marg.pcov <- diag(updiag) %*% marg.pcov %*% diag(updiag)
        # Update covariance matrix of the correlation function and degrees of freedom proposals
        updiag <- sqrt(diag(dep.pcov))
        for (j in 1:ndep) {
          ada <- adaptive(attempts = attempt[dep.i[j]], acceptance = accept[dep.i[j]], sd.p = updiag[j])
          updiag[j] <- ada$sd / updiag[j]
          accept[dep.i[j]] <- ada$acc
          attempt[dep.i[j]] <- ada$att
        }
        dep.pcov <- diag(updiag) %*% dep.pcov %*% diag(updiag)
      } else if (b > 2000L && b < burnin && (b %% 200) == 0) {
        mb <- max(b - 1000, 200)
        marg.pcov <- marg.pcov + 0.1 * (cov(res[mb:b, c(scale.i, shape.i)]) + 1e-4 * diag(D + 1) - marg.pcov)
        dep.pcov <- dep.pcov + 0.1 * (cov(transform.fn(res[mb:b, dep.i])) + 1e-4 * diag(ncol(dep.pcov)) - dep.pcov)
      }
      if (b > 200 && b < burnin && (b %% 200) == 0) {
        if (geoaniso) {
          mb <- max(b - 1000, 200)
          aniso.pcov <- aniso.pcov + 0.1 * (cov(res[mb:b, aniso.i]) + 1e-5 * diag(2) - aniso.pcov)
        }
      }
    }
    if (b %% verbose == 0) {
      cat("Current values: ", round(c(mean(scale.c), shape.c, dep.c), 2), "\n")
      elapsed.time <- round((proc.time()[3] - time.0) / 60 / 60, 2)
      remaining.time <- round((elapsed.time / b) * (B - b), 2)
      print(paste("Iteration", b, "out of", B, "completed at", Sys.time()))
      cat("  Elapsed time:", elapsed.time, "hours\n")
      cat("  Remaining time:", remaining.time, "hours\n\n")
    }
    if (b %% saveinterm == 0) {
      save(res, dat, Xm, lpost, dep.pcov, marg.pcov, aniso.pcov, file = paste0(filename, ".RData"))
    }
  }
  time <- round((proc.time()[3] - time.0) / 60 / 60, 2)
  save(res, time, dat, Xm, lpost, dep.pcov, marg.pcov, aniso.pcov, file = paste0(filename, ".RData"))
  invisible(return(list(
    res = res, time = time, dat = dat, Xm = Xm, coord = coord, lpost = lpost,
    dep.pcov = dep.pcov, marg.pcov = marg.pcov, aniso.pcov = aniso.pcov
  )))
}

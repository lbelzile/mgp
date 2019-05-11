#' Two-step Markov chain Monte Carlo algorithm for multivariate generalized Pareto models
#'
#' The algorithm estimates marginal parameters in a preliminary run, then fix the latter to the posterior median
#' and samples from the dependence parameters in a second time.
#'
#' @inheritParams mcmc.mgp
#' @inheritDotParams clikmgp
#' @export
#' @keywords internal
#' @return a list with \code{res} containing the results of the chain
twostep.mgp <- function(dat, mthresh, thresh, lambdau = 1, model = c("br", "xstud"), coord, start,
                        numiter = 4e4L, burnin = 5e3L, thin = 1L, verbose = 50L, transform = FALSE, censor = TRUE, filename,
                        keepburnin = TRUE, geoaniso = TRUE, likt = c("mgp", "pois", "binom"),
                        saveinterm = 500L, ...) {
  slurm_arrayid <- Sys.getenv("SLURM_ARRAY_TASK_ID")
  slurm_jobid <- Sys.getenv("SLURM_JOB_ID")
  filename <- paste0(
    filename, ifelse(slurm_arrayid == "", "", "_"), slurm_arrayid, "_",
    model, "_twostep", ifelse(slurm_jobid == "", "", "_"), slurm_jobid
  )

  model <- match.arg(model)
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
  if (censor) {
    mmin <- rep(0, D)
  } else {
    mmin <- apply(dat, 2, min)
  }
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
  if (length(shape.c) == 1L) {
    cshape <- TRUE
  } else {
    cshape <- FALSE
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


  if (is.null(scale.c) || is.null(shape.c)) {
    margpars <- sapply(1:D, function(j) {
      mev::fit.gpd(xdat = as.vector(dat[, j]), threshold = mthresh[j])$est
    })
    scale.c <- as.vector(margpars[1, ])
    shape.c <- as.vector(margpars[2, ])
    cshape <- FALSE
  } else if (!is.null(shape.c)) {
    cshape <- ifelse(length(shape.c) > 1, FALSE, TRUE)
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
  scale.i <- 1:length(scale.c)
  shape.i <- (length(scale.c) + 1:length(shape.c))
  dep.i <- (shape.i[length(shape.i)] + 1):(shape.i[length(shape.i)] + length(dep.c))
  if (geoaniso) {
    aniso.i <- (npar - 1):npar
  }
  if (model == "xstud") {
    df.c <- start$df
    df.pcov <- 0.2
    if (is.null(df.c)) {
      df.c <- 2
    }
    df.lb <- 1 + 1e-5
    npar <- npar + 1L
    df.i <- npar
    # Port to dep. to perform updates together with dep
    dep.c <- c(dep.c, df.c)
    dep.npcov <- matrix(0, ndep + 1, ndep + 1)
    dep.npcov[1:ndep, 1:ndep] <- dep.pcov
    dep.npcov[ndep + 1, ndep + 1] <- df.pcov
    dep.pcov <- dep.npcov
    dep.lb <- c(dep.lb, df.lb)
    dep.ub <- c(dep.ub, Inf)
    dep.i <- c(dep.i, df.i)
    ndep <- ndep + 1L
  } else {
    df.pcov <- NULL
  }

  if (transform) {
    transform.fn <- function(x) {
      transfo.pars(x, lbound = dep.lb, ubound = dep.ub)
    }
  } else {
    transform.fn <- identity
  }
  ckst <- 2.38 * 2.38 / ndep
  lscalelm.i <- (npar + 1):(npar + ncol(Xm) + 2)
  npar <- npar + ncol(Xm) + 2
  if (!cshape) {
    shapelm.i <- (npar + 1):(npar + ncol(Xm) + 2)
    npar <- npar + ncol(Xm) + 2
  }
  # Function to create list with parameter, input is model dependent
  makepar <- function(dep, model = c("br", "xstud"), distm) {
    model <- match.arg(model)
    if (model == "xstud") {
      Sigma <- dep.fun(distm, par = dep)
      return(list(Sigma = Sigma, df = dep[ndep]))
    } else if (model == "br") {
      Lambda <- dep.fun(distm, par = dep)
      return(par <- list(Lambda = Lambda))
    }
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

  par.c <- makepar(dep = dep.c, model = model, distm = distm.c)
  # keep only marginal exceedances
  ldat <- apply(dat, 2, function(col) {
    as.vector(col[col > 0])
  })
  # Define containers
  lpost <- rep(0, B)
  accept <- rep(0L, ndep)
  attempt <- rep(0L, ndep)
  res <- matrix(data = 0, ncol = npar, nrow = B)
  margpar <- mcmc.mgp(
    dat = dat, mthresh = mthresh, thresh = thresh,
    lambdau = lambdau, model = "lgm", coord = coord,
    start = list(scale = scale.c, shape = shape.c),
    numiter = numiter, burnin = burnin, thin = thin,
    verbose = B, saveinterm = Inf, filename = filename
  )
  # Copy results and transform margins
  res[, scale.i] <- margpar$res[, scale.i]
  res[, shape.i] <- margpar$res[, shape.i]
  res[, lscalelm.i] <- margpar$res[, (max(shape.i) + 1):(max(shape.i) + 6)]
  scale.c <- apply(margpar$res[-(1:min(burnin, 2000)), scale.i], 2, median)
  if (!cshape) {
    shape.c <- apply(margpar$res[-(1:min(burnin, 2000)), shape.i], 2, median)
    res[, shapelm.i] <- margpar$res[, (max(shape.i) + 7):(max(shape.i) + 12)]
  } else {
    shape.c <- median(margpar$res[-(1:min(burnin, 2000)), shape.i], na.rm = TRUE)
  }
  par.c <- makepar(dep = dep.c, model = model, distm = distm.c)
  if (!censor) {
    loglikfn <- function(scale, shape, par, ...) {
      # Initial evaluation of the log-likelihood
      likmgp(
        dat = dat, thresh = thresh, loc = rep(0, D), scale = scale,
        shape = shape, par = par, model = model, likt = likt, lambdau = lambdau,
        B1 = B1, genvec1 = genvec1, ncores = ncores, ntot = ntot
      )
    }
  } else if (censor) {
    loglikfn <- function(scale, shape, par, ...) {
      # Initial evaluation of the log-likelihood
      clikmgp(
        dat = dat, thresh = thresh, loc = rep(0, D), scale = scale,
        shape = shape, par = par, model = model, likt = likt, lambdau = lambdau,
        B1 = B1, B2 = B2, genvec1 = genvec1, genvec2 = genvec2, ncores = ncores,
        numAbovePerRow = numAbovePerRow, numAbovePerCol = numAbovePerCol, censored = censored, ntot = ntot
      )
    }
  }
  loglik.c <- loglikfn(scale = scale.c, shape = shape.c, par = par.c)
  # Start timer
  time.0 <- proc.time()[3] # Better would be to time every run and keep median
  ########################################################################
  ####################          LOOP         #############################
  ########################################################################
  thin <- as.integer(thin)
  browser()
  for (b in 1:B) {
    attempt <- attempt + 1L

    ####           UPDATE DEPENDENCE PARAMETERS             ####
    # adding prior contribution for the anisotropy parameters

    if (model == "xstud") {
      dep.loglikfn <- function(dep) {
        Sigma <- dep.fun(distm.c, dep)
        # par = list(Sigma = Sigma, df = par.c$df)
        par <- list(Sigma = Sigma, df = dep[ndep])
        ll <- loglikfn(scale = scale.c, shape = shape.c, par = par)
        attributes(ll)$par <- par
        ll
      }
    } else {
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
        accept[dep.i] <- accept[dep.i] + 1L
        par.c <- attributes(update$ll)$par
        loglik.c <- update$ll
        dep.c <- update$cur
      }
    }

    if (geoaniso) {
      if (model == "xstud") {
        aniso.loglikfn <- function(aniso) {
          aniso[2] <- wrapAng(aniso[2])
          distm <- distg(loc = loc, scale = aniso[1], rho = aniso[2])
          Sigma <- dep.fun(distm, dep.c)
          par <- list(Sigma = Sigma, df = dep.c[ndep])
          ll <- loglikfn(scale = scale.c, shape = shape.c, par = par)
          attributes(ll)$par <- par
          attributes(ll)$distm <- distm
          ll
        }
      } else {
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
      res[i, dep.i] <- dep.c
      if (geoaniso) {
        res[i, aniso.i] <- aniso.c
      }
    }

    # Adapt covariance matrix, but using only previous iterations
    # Stop adapting after burnin
    if (b < burnin) {
      if (b %% 20 * thin == 0 && b > 200L && b <= 2000L) {
        # Update covariance matrix of the correlation function and degrees of freedom proposals
        updiag <- sqrt(diag(dep.pcov))
        for (j in 1:ndep) {
          ada <- adaptive(attempts = attempt[j], acceptance = accept[j], sd.p = updiag[j])
          updiag[j] <- ada$sd / updiag[j]
          accept[j] <- ada$acc
          attempt[j] <- ada$att
        }
        dep.pcov <- diag(updiag) %*% dep.pcov %*% diag(updiag)
      }
      if (b > 2000L && b < burnin && (b %% 200) == 0) {
        mb <- max(b - 1000, 200)
        dep.pcov <- dep.pcov + 0.1 * (cov(transform.fn(res[mb:b, dep.i])) + 1e-4 * diag(ncol(dep.pcov)) - dep.pcov)
      }
      if (b > 200 && b < burnin && (b %% 200) == 0) {
        if (geoaniso) {
          mb <- max(b - 1000, 200)
          aniso.pcov <- aniso.pcov + 0.1 * (cov(res[mb:b, aniso.i]) + 1e-4 * diag(2) - aniso.pcov)
        }
      }
    }
    if (b %% verbose == 0) {
      elapsed.time <- round((proc.time()[3] - time.0) / 60 / 60, 2)
      remaining.time <- round((elapsed.time / b) * (B - b), 2)
      print(paste("Iteration", b, "out of", B, "completed at", Sys.time()))
      cat("  Elapsed time:", elapsed.time, "hours\n")
      cat("  Remaining time:", remaining.time, "hours\n\n")
    }
    if (b %% saveinterm == 0) {
      save(res, dat, Xm, lpost, dep.pcov, aniso.pcov, model, file = paste0(filename, ".RData"))
    }
  }
  time <- round((proc.time()[3] - time.0) / 60 / 60, 2)
  save(res, time, dat, Xm, lpost, dep.pcov, aniso.pcov, model, file = paste0(filename, ".RData"))
  invisible(return(list(
    res = res, time = time, dat = dat, Xm = Xm, coord = coord, lpost = lpost,
    dep.pcov = dep.pcov, aniso.pcov = aniso.pcov, model = model
  )))
}

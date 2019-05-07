#' Update range parameter via Metropolis-Hastings
#'
#' This function, adapted from \code{\link[spatial.gev.bma]{gev.update.lambda}},
#' uses the Laplace approximation to build a proposal for the range parameter of an
#' exponential covariance model with precision \eqn{\alpha} and range \eqn{\lambda}.
#'
#' The difference with the \code{spatial.gev.bma} package is that proposals
#' are truncated Gaussian, so may be accepted even if the Newton discount is centered at negative values.
#' If the current \eqn{\lambda} value is in a region where the Hessian matrix is not negative, the proposal
#' is automatically rejected. The function \code{discount} controls the size of the Newton steps.
#'
#' @param tau current value of centered random effects
#' @param alpha precision of the random effect
#' @param lambda range parameter to be updated
#' @param di matrix of pairwise distances between points
#' @param a shape parameter of the Gamma prior for \eqn{\lambda}
#' @param b rate parameter of the Gamma prior for \eqn{\lambda}
#' @param lb lower bound for admissible value of \eqn{\lambda}. Default to 0.01.
#' @param discount size of discount in Newton update; default to 0.2.
#' @param maxstep maximum discount size for the truncated Normal proposal
#' @author code from spatial.gev.bma by Alex Lenkoski
#' @return new value of \code{lambda}
#' @export
#' @seealso \code{\link[spatial.gev.bma]{gev.update.lambda}}
updt.range <- function (tau, alpha, lambda, di, a, b, lb = 1e-2, discount = 0.2, maxstep = Inf){
  attributes(lambda) <- list(accept = FALSE)
  logdet <- function(A) {return(2 * sum(log(diag(chol(A)))))}
  l.prime <-  function (tau, alpha, lambda, D, a, b){
    E.l <- exp(-1/lambda * D)
    diag(E.l) <- diag(E.l) + 1e-05
    E.inv <- solve(E.l)
    F.l <- 1/lambda^2 * D * E.l
    M.l <- E.inv %*% (-F.l) %*% E.inv
    res <- -0.5 * sum(diag(E.inv %*% F.l)) - 0.5 * alpha * t(tau) %*% M.l %*% tau - b + (a - 1)/lambda
    return(res[1])
  }
  l.double.prime <- function (tau, alpha, lambda, D, a, b) {
    E.l <- exp(-1/lambda * D)
    diag(E.l) <- diag(E.l) + 1e-05
    E.inv <- solve(E.l)
    F.l <- 1/lambda^2 * D * E.l
    M.l <- E.inv %*% (-F.l) %*% E.inv
    G.l <- -2/lambda^3 * (D * E.l) + 1/lambda^2 * (D * F.l)
    L.l <- M.l %*% F.l + E.inv %*% G.l
    N.l <- M.l %*% (-F.l) %*% E.inv + E.inv %*% (-G.l) %*% E.inv + E.inv %*% (-F.l) %*% M.l
    res <- -0.5 * sum(diag(L.l)) - 0.5 * alpha * t(tau) %*% N.l %*%  tau - (a - 1) * lambda^(-2)
    return(res[1])
  }
  stopifnot(maxstep > 0)
  l.curr <- l.prime(tau, alpha, lambda, di, a, b)
  l.double.curr <- l.double.prime(tau, alpha, lambda, di, a, b)
  d.curr <- -l.double.curr
  b.curr <- discount * l.curr / d.curr + lambda
if (isTRUE(d.curr > 0)) {
  lambda.new <- TruncatedNormal::mvrandn(pmax(lb, lambda - maxstep), lambda + maxstep, Sig = matrix(1/d.curr), n = 1, mu = b.curr)
  #rnorm(1, b.curr/d.curr, sd = sqrt(1/d.curr))
  l.new <- l.prime(tau, alpha, lambda.new, di, a, b)
  l.double.new <- l.double.prime(tau, alpha, lambda.new, di, a, b)
  d.new <- -l.double.new
  b.new <- discount * l.new / d.new + lambda.new
  if (isTRUE(d.new > 0)) {
    E.l.curr <- exp(-1/lambda * di)
    diag(E.l.curr) <- diag(E.l.curr) + 1e-05
    E.l.curr.inv <- solve(E.l.curr)
    E.l.new <- exp(-1/lambda.new * di)
    diag(E.l.new) <- diag(E.l.new) + 1e-05
    E.l.new.inv <- solve(E.l.new)
    L.curr <- -0.5 * alpha * t(tau) %*% E.l.curr.inv %*%
      tau - 0.5 * logdet(E.l.curr)
    L.new <- -0.5 * alpha * t(tau) %*% E.l.new.inv %*%
      tau - 0.5 * logdet(E.l.new)
    prior.curr <- dgamma(lambda, a, b, log = TRUE)
    prior.new <- dgamma(lambda.new, a, b, log = TRUE)
    prop.curr <- dnorm(lambda.new, b.curr, sd = sqrt(1/d.curr), log = TRUE) -
      log(suppressWarnings(TruncatedNormal::mvNcdf(l = pmax(lb, lambda - maxstep) - b.curr,
                                                   u = lambda + maxstep - b.curr, Sig = 1/d.curr, n = 1)$prob))
    prop.new <- dnorm(lambda, b.new, sd = sqrt(1/d.new), log = TRUE) -
      log(suppressWarnings(TruncatedNormal::mvNcdf(l = pmax(lb, lambda.new - maxstep) - b.new,
                                                   u = lambda.new + maxstep - b.new, Sig = 1/d.new, n = 1)$prob))
    mh <- L.new - L.curr + prior.new - prior.curr + prop.new - prop.curr
    if (isTRUE(log(runif(1)) < mh)) {
      lambda <- lambda.new
      attributes(lambda) <- list(accept = TRUE)
    }
  }
}
return(lambda)
}



#' Proposals for random walk Metropolis-Hastings
#'
#' This function transforms a vector \code{cur} to an unconstrained scale based on \code{lbound} and \code{ubound},
#' then samples draws from a multivariate Gaussian vector with covariance matrix \code{cov} centered at the current (unconstrained) value.
#' The function then returns a list containing the proposal on the unconstrained scale (\code{trprop}),
#' on the original scale (\code{prop}) along with the log of the jacobian of the transformation (\code{logjac}).
#'
#' @param cur vector of current parameter on the original (constrained) scale
#' @param trcur vector of current parameter values on the unconstrained scale
#' @param cov covariance matrix of the proposal for the \code{trcur} scale
#' @param lbound lower bounds for the parameters
#' @param ubound upper bounds for the parameter
#' @return a list with components \code{prop}, \code{trprop} and \code{logjac}
#' @export
#' @keywords internal
propRWMH <- function(cur = NULL, trcur = NULL, cov, lbound = rep(-Inf, ncol(cov)), ubound = rep(Inf, ncol(cov))){
  if(is.null(trcur) && is.null(cur)){
    stop("Either `trcur` or `cur` must be provided")
  } else if(!is.null(trcur)){
    d <- length(trcur)
    stopifnot(d == dim(cov), length(lbound) == d, length(ubound) == d)
    tr.cur <-  trcur
  } else{
    d <- length(cur)
    tr.cur <- rep(0, d);
    stopifnot(d == dim(cov), length(lbound) == d, length(ubound) == d)
    for(j in 1:d){
      if(lbound[j] == -Inf && ubound[j] == Inf){
        tr.cur[j] <- cur[j]
      } else if(lbound[j] > -Inf && ubound[j] == Inf){
        tr.cur[j] <- log(cur[j] - lbound[j])
      } else if(lbound[j] == -Inf && ubound[j] < Inf){
        tr.cur[j] <- log(ubound[j] - cur[j])
      } else{
        tr.cur[j] <-  logit((cur[j] - lbound[j]) / (ubound[j] - lbound[j]))
      }
    }
  }
  prop <- rep(0, d)
  logratio <- rep(0, d)
  expit <- function(x){1/(1+exp(-x))}
  logit <- function(x){log(x)-log(1-x)}
  if(d == 1L){
    tr.prop <- rnorm(n = 1, mean = tr.cur, sd =  sqrt(cov))
  } else{
    tr.prop <- mev::mvrnorm(n = 1, mu = tr.cur, Sigma = cov)
  }
  for(j in 1:d){
    if(lbound[j] == -Inf && ubound[j] == Inf){
      prop[j] <- tr.prop[j]
    } else if(lbound[j] > -Inf && ubound[j] == Inf){
      prop[j] <- lbound[j] + exp(tr.prop[j])
      logratio[j] <- tr.prop[j] - tr.cur[j]
    } else if(lbound[j] == -Inf && ubound[j] < Inf){
      prop[j] <- ubound[j] - exp(tr.prop[j])
      logratio[j] <- tr.prop[j] - tr.cur[j]
    } else{
      prop[j] <- lbound[j] + (ubound[j] - lbound[j]) * expit(tr.prop[j])
      logratio[j] <- log(1/(cur[j] - lbound[j]) + 1/(ubound[j] - cur[j])) - log(1/(prop[j] - lbound[j]) + 1/(ubound[j] - prop[j]))
    }
  }
  return(list(prop = prop, trprop = tr.prop, logjac = logratio))
}

#' Variance adaptation
#'
#' Adapt standard deviation of a proposal for a Markov chain Monte Carlo algorithm, based on the acceptance rate.
#' The function is targeting an acceptance probability of 0.44 for univariate Metropolis--Hastings proposals, which is
#' optimal for certain Gaussian regimes. If the number of attempts is large enough and the acceptance rate is not close
#' to the target, the standard deviation of the proposal will be increased or reduced and the number acceptance and attempts
#' reinitialized. If no change is made, the components \code{acc} will be the same as acceptance and similarly for attempts.
#'
#' @param attempts integer indicating number of attempts
#' @param acceptance integer giving the number of moves accepted by the algorithm
#' @param sd.p standard deviation of proposal
#' @return a list with components \code{sd}, \code{acc}, \code{att}
#' @export
#' @keywords internal
adaptive <- function(attempts, acceptance, sd.p){
  stopifnot(sd.p > 0)
  if(attempts < acceptance){
    stop("Invalid input: the number of attempts must be larger than the number of acceptance")
  }
  att <- attempts
  acc <- acceptance
  newsd <- sd.p
  if(att > 20){
    if(acc / att > 0.9){
      newsd <- sd.p * 1.5;
      att = 0L; acc <- 0L
    } else if(acc / att < 0.05){
      newsd <- sd.p * 0.5;
      att = 0L; acc <- 0L
    }
  }
  if(att > 30){
    if(acc / att > 0.5){
      newsd <- sd.p * 1.25;
      att = 0L; acc <- 0L
    } else if(acc / att < 0.10){
      newsd <- sd.p * 0.75;
      att = 0L; acc <- 0L
    }
  }
  if(att > 50){
    if(acc / att > 0.35){
      newsd <- sd.p * 1.1;
      att = 0L; acc <- 0L
    } else if(acc / att < 0.2){
      newsd <- sd.p * 0.95;
      att = 0L; acc <- 0L
    }
  }

  return(list(sd = newsd, acc = acc, att = att))
}

#' Wrap angular component in anisotropy
#' @param ang angle
#' @return angle between \eqn{(-\pi/2,\pi/2)}
#' @keywords internal
#' @export
wrapAng <- function(ang){
  stopifnot(length(ang) == 1L)
  if(ang > -pi/2 && ang < pi/2){
    return(ang)
  } else{
    stang <- ang %% (2*pi)
    if(stang > 1.5*pi){
      return(stang - 2*pi)
    } else if(stang < 0.5*pi){
      return(stang)
    } else if(stang > pi){
      return(-0.5*pi + (1.5*pi - stang))
    } else{
      return(pi/2 - (stang - pi/2))
    }
  }
}





#' Update for latent Gaussian model for the scale parameter of a generalized Pareto
#'
#' The scale has a log-Gaussian prior with variance \code{lscale.tausq} and precision (the inverse of the correlation matrix)
#' given by \code{lscale.precis}.
#'
#' @param scale vector of scale parameters for the generalized Pareto distribution
#' @param shape vector of shape parameters for the generalized Pareto distribution
#' @param ldat list with exceedances at each site
#' @param lscale.mu mean of the log-Gaussian process for the scale
#' @param lscale.precis precision matrix of the log-Gaussian process corresponding to the inverse of the correlation matrix
#' @param lscale.tausq variance of the log-Gaussian process
#' @param mmax vector of maximum of each series in \code{ldat}
#' @param discount numeric giving the discount factor for the Newton  Default to 1.
#' @param maxstep maximum step size for the MCMC (proposal will be at most \code{maxstep} units away from the current value).
#' @return a vector of scale parameters
#' @export
lscale.lgm <- function(scale, shape, ldat, lscale.mu, lscale.precis, lscale.tausq, mmax, discount = 1, maxstep = Inf){
  # Gradient and Hessian for the log-gaussian random effect model
  gradlgauss <- function(x, mu, Q){ c(-1/x- Q %*% (log(x)-mu)*(1/x))}
  hesslgauss <- function(x, mu, Q){ if(length(x) > 1L){
    diag(1/x^2) - Q * tcrossprod(1/x,1/x) + diag(c(Q %*% (log(x)-mu))/(x^2))
  } else{
    (1- Q) / x^2 + Q * (log(x)-mu)/(x^2)
  }}
  # likelihood function
  llfun <- function(scale, shape, ind = 1:D, ...){
    scale <- rep(scale, length.out = D)
    shape <- rep(shape, length.out = D)
    sum(sapply(ind, function(i){mev::gpd.ll(par = c(scale[i], shape[i]), dat = ldat[[i]])}))
  }
  if(length(shape) == 1L){
    cshape <- TRUE
  } else{
    stopifnot(length(scale) == length(shape))
    cshape <- FALSE
  }

  # Componentwise maximum at site to compute bounds for simulating scale
  if(missing(mmax)){
    mmax <- as.vector(unlist(lapply(ldat, max)))
  }


  #Start update
  D <- length(scale)
  for(k in sample.int(D, D)){
    #Conditional mean and precision for log-Gaussian - both do not depend on current value of "k"
    condmean <- c(lscale.mu[k] - solve(lscale.precis[k,k, drop = FALSE]) %*% lscale.precis[k,-k, drop = FALSE] %*% (log(scale[-k]) - lscale.mu[-k]))
    condprec <- lscale.precis[k,k, drop = FALSE]/lscale.tausq
    #Lower bound for the scale - simulations are from a truncated Gaussian
    if(cshape){
      lbscale <- pmax(-shape*mmax[k], 0, scale[k] - maxstep)
      ubscale <- pmin(Inf, scale[k] + maxstep)
    } else{
      lbscale <- pmax(-shape[k]*mmax[k], 0, scale[k] - maxstep)
      ubscale <- pmin(Inf, scale[k] + maxstep)
    }
    #Copy scale, perform Newton discount from current value
    scale.p <- scale
    #Gradient and Hessian at current value
    fp <- sapply(k, function(j){mev::gpd.score(par = c(scale[j], ifelse(cshape, shape, shape[j])), dat = ldat[[j]])[1]}) +
      gradlgauss(x = scale[k], mu = condmean, Q = condprec)
    if(length(k) > 1L){
      fpp <- diag(sapply(k, function(j){-mev::gpd.infomat(par = c(scale[j], ifelse(cshape, shape, shape[j])), dat = ldat[[j]], method = "obs")[1,1]})) +
        hesslgauss(x = scale[k], mu = condmean, Q = condprec)
    } else{
      fpp <- -mev::gpd.infomat(par = c(scale[k], ifelse(cshape, shape, shape[k])), dat = ldat[[k]], method = "obs")[1,1] +
        hesslgauss(x = scale[k], mu = condmean, Q = condprec)
    }
    continue <- isTRUE(all(eigen(fpp, only.values = TRUE)$values < 0))
    if(continue){
      # Mean and scale for the proposal - with a damping factor
      Sig1 <- -as.matrix(solve(fpp))
      mu1 <- c(scale[k] + discount * Sig1 %*% fp);
      scale.p[k] <- TruncatedNormal::mvrandn(l = lbscale, u = ubscale, mu = mu1, Sig = Sig1, n = 1)
      #Taylor approximation from proposal - obtain mean and variance for Laplace approximation
      fp2 <- sapply(k, function(j){mev::gpd.score(par = c(scale.p[j], ifelse(cshape, shape, shape[j])), dat = ldat[[j]])[1]}) +
        gradlgauss(x = scale.p[k], mu = condmean, Q = condprec)
      if(length(k) > 1L){
        fpp2 <- diag(sapply(k, function(j){-mev::gpd.infomat(par = c(scale.p[j], ifelse(cshape, shape,shape[j])), dat = ldat[[j]], method = "obs")[1,1]})) +           hesslgauss(x = scale.p[k], mu = condmean, Q = condprec)
      } else{
        fpp2 <- -mev::gpd.infomat(par = c(scale.p[k], ifelse(cshape, shape, shape[k])), dat = ldat[[k]], method = "obs")[1,1] +
          hesslgauss(x =scale.p[k], mu = condmean, Q = condprec)
      }
      continue <- isTRUE(all(eigen(fpp2, only.values = TRUE)$values < 0))
      if(continue){
        Sig2 <- -as.matrix(solve(fpp2))
        mu2 <- c(scale.p[k] + discount * Sig2 %*% fp2);
       # Jacobian of transformation
        jac <- mgp::.dmvnorm_arma(x = t(as.matrix(as.vector(scale[k]))), mean = as.vector(mu2), sigma = as.matrix(Sig2), logd = TRUE) -
          mgp::.dmvnorm_arma(x = t(as.matrix(as.vector(scale.p[k]))), mean  = as.vector(mu1), sigma = as.matrix(Sig1), logd = TRUE) -
          suppressWarnings(log(TruncatedNormal::mvNcdf(l = lbscale - mu2, u = rep(Inf, length(lbscale)), Sig = as.matrix(Sig2), n = 1000)$prob)) +
          suppressWarnings(log(TruncatedNormal::mvNcdf(l = lbscale - mu1, u = rep(Inf, length(lbscale)), Sig = as.matrix(Sig1), n = 1000)$prob))

        #Evaluate log-density and log-prior
        logdens <- llfun(scale = scale, shape = rep(shape, length.out = length(scale)), ldat = ldat, ind = k)
        logdens.p <- llfun(scale = scale.p, shape = rep(shape, length.out = length(scale)), ldat = ldat, ind = k)
        # Need not be computed for every update of the parameters - so decouple
        if(length(k) > 1L){
          priorA.p <-  dmvnorm.precis(x = log(scale.p[k]), mean = condmean, precis = condprec, logd = TRUE) - sum(log(scale.p[k]))
          priorA <-  dmvnorm.precis(log(scale[k]), mean = condmean, precis = condprec, logd = TRUE) - sum(log(scale[k]))
        } else{
          priorA.p <-  dnorm(log(scale.p[k]), mean = condmean, sd = sqrt(1/condprec), log = TRUE) - log(scale.p[k])
          priorA <-  dnorm(log(scale[k]), mean = condmean, sd = sqrt(1/condprec), log = TRUE) - log(scale[k])
        }
        acpt <- logdens.p - logdens + priorA.p - priorA + jac
        if(log(runif(1)) < acpt){
          scale <- scale.p
        }
      }
    }
  }
  return(scale)
}


#' Update for latent Gaussian model for the scale parameter of a generalized Pareto
#'
#' The scale has a log-Gaussian prior with variance \code{lscale.tausq} and precision (the inverse of the correlation matrix)
#' given by \code{lscale.precis}.
#'
#' @param scale vector of scale parameters for the generalized Pareto distribution
#' @param shape vector of shape parameters for the generalized Pareto distribution
#' @param ldat list with exceedances at each site
#' @param shape.mu mean of the Gaussian process for the shape parameters
#' @param shape.precis precision matrix of the Gaussian process corresponding to the inverse of the correlation matrix
#' @param shape.tausq variance of the Gaussian process for the shape parameters
#' @param lbound lower bound if parameter is truncated below
#' @param ubound upper bound if parameter is truncated above
#' @param mmax vector of maximum of each series in \code{ldat}
#' @param discount numeric giving the discount factor for the Newton  Default to 1.
#' @param maxstep maximum step size for the MCMC (proposal will be at most \code{maxstep} units away from the current value).
#' @return a vector of scale parameters
#' @export
shape.lgm <- function(scale, shape, ldat, shape.mu = NULL, shape.precis = NULL, shape.tausq = NULL, mmax, lbound = -0.5, ubound = 0.5, discount = 1, maxstep = 0.1){
  D <- length(scale)
  stopifnot(length(ldat) == length(scale))
  # likelihood function
  llfun <- function(scale, shape, ind = 1:D, ...){
    scale <- rep(scale, length.out = D)
    shape <- rep(shape, length.out = D)
    sum(sapply(ind, function(i){mev::gpd.ll(par = c(scale[i], shape[i]), dat = ldat[[i]])}))
  }
  if(length(shape) == 1L){
    cshape <- TRUE
  } else{
    stopifnot(length(scale) == length(shape), !is.null(shape.mu), !is.null(shape.precis), !is.null(shape.tausq))
    cshape <- FALSE
  }

  # Componentwise maximum at site to compute bounds for simulating scale
  if(missing(mmax)){
    mmax <- as.vector(unlist(lapply(ldat, max)))
  }

for(k in sample.int(n = length(shape), size = length(shape), replace = FALSE)){
  if(cshape){
    condmean <- 0; condprec <- 25
    lbxi <- pmax(lbound, -min(scale/mmax), shape - maxstep)
    ubxi <- pmin(ubound, shape + maxstep)
  } else{
    condmean <- c(shape.mu[k] - solve(shape.precis[k,k, drop = FALSE]) %*% shape.precis[k,-k, drop = FALSE] %*% (shape[-k] - shape.mu[-k]))
    condprec <- shape.precis[k,k, drop = FALSE]/shape.tausq
    lbxi <- pmax(lbound,-scale[k]/mmax[k])
    ubxi <- pmin(ubound,pmax(shape[k] + maxstep))
  }
  #Copy shape, perform Newton step from current value
  shape.p <- shape
  if(cshape){
    fp <- sum(sapply(1:D, function(j){mev::gpd.score(par = c(scale[j], shape), dat = ldat[[j]])[2]})) +
      -c(condprec %*% (shape[k] - condmean))
    fpp <- sum(sapply(1:D, function(j){-mev::gpd.infomat(par = c(scale[j], shape), dat = ldat[[j]], method = "obs")[2,2]})) - condprec
  } else{
    #Gradient and Hessian at current value
    fp <- sapply(k, function(j){mev::gpd.score(par = c(scale[j], shape[j]), dat = ldat[[j]])[2]}) +
      -c(condprec %*% (shape[k] - condmean))
    if(length(k) > 1L){
      fpp <- diag(sapply(k, function(j){-mev::gpd.infomat(par = c(scale[j], shape[j]), dat = ldat[[j]], method = "obs")[2,2]})) - condprec
    } else{
      fpp <- -mev::gpd.infomat(par = c(scale[k], shape[k]), dat = ldat[[k]], method = "obs")[2,2] - condprec
    }
  }
  continue <- isTRUE(all(eigen(fpp, only.values = TRUE)$values < 0))
  if(continue){
    # Mean and scale for the proposal - with a damping factor to avoid oscillations
    Sig1 <- -as.matrix(solve(fpp))
    mu1 <- shape[k] + discount * Sig1 %*% fp
    shape.p[k] <- TruncatedNormal::mvrandn(l = lbxi, u = ubxi, mu = mu1, Sig = Sig1, n = 1)
    if(cshape){
      fp2 <- sum(sapply(1:D, function(j){mev::gpd.score(par = c(scale[j], shape.p), dat = ldat[[j]])[2]})) +
        -c(condprec %*% (shape.p[k] - condmean))
      fpp2 <- sum(sapply(1:D, function(j){-mev::gpd.infomat(par = c(scale[j], shape.p), dat = ldat[[j]], method = "obs")[2,2]})) - condprec
    } else{
      #Taylor approximation from proposal - obtain mean and variance for Laplace approximation
      fp2 <- sapply(k, function(j){mev::gpd.score(par = c(scale[j], shape.p[j]), dat = ldat[[j]])[2]}) -
        c(condprec %*% (shape.p[k] - condmean))
      if(length(k) > 1L){
        fpp2 <- diag(sapply(k, function(j){-mev::gpd.infomat(par = c(scale[j], shape.p[j]), dat = ldat[[j]], method = "obs")[2,2]})) - condprec
      } else{
        fpp2 <- -mev::gpd.infomat(par = c(scale[k], shape.p[k]), dat = ldat[[k]], method = "obs")[2,2] - condprec
      }
    }
    continue <- isTRUE(all(eigen(fpp2, only.values = TRUE)$values < 0))
    if(continue){
      Sig2 <- -as.matrix(solve(fpp2))
      mu2 <- shape.p[k] + discount * Sig2 %*% fp2
      # Jacobian of transformation

        jac <- mgp::.dmvnorm_arma(x = matrix(shape[k], ncol = length(k)), mean = as.vector(mu2), sigma = as.matrix(Sig2), logd = TRUE) -
          mgp::.dmvnorm_arma(x = matrix(shape.p[k], ncol = length(k)), mean = as.vector(mu1), sigma = as.matrix(Sig1), logd = TRUE) -
          suppressWarnings(log(TruncatedNormal::mvNcdf(l = lbxi - mu2, u = ubxi - mu2, Sig = as.matrix(Sig2), n = 1000)$prob)) +
          suppressWarnings(log(TruncatedNormal::mvNcdf(l = lbxi - mu1, u = ubxi - mu1, Sig = as.matrix(Sig1), n = 1000)$prob))

      if(cshape){
        logprior <- dnorm(shape, 0, sd = sqrt(1/condprec), log = TRUE)
        logprior.p <- dnorm(shape.p, 0, sd = sqrt(1/condprec), log = TRUE)
      } else{
        logprior.p <-  mgp::dmvnorm.precis(shape.p[k], mean = condmean, precis = condprec, logd = TRUE)
        logprior <-  mgp::dmvnorm.precis(shape[k], mean = condmean, precis = condprec, logd = TRUE)
      }
      logdens <- llfun(scale = scale, shape = shape, ind = if(cshape){1:D} else{k})
      logdens.p <- llfun(scale = scale, shape = shape.p, ind = if(cshape){1:D} else{k})
      acpt <- logdens.p - logdens + logprior.p - logprior + jac
      if(log(runif(1)) < acpt){
        shape <- shape.p
      }
    }
  }
}
  return(shape)
}

#' Update function based on Metropolis-Hasting algorithm
#'
#' Proposals are made based on a (conditional) truncated Gaussian distribution at indices \code{ind}
#' given other components.
#' @param cur current value of the vector of parameters of which one component is to be updated.
#' @param lb lower bounds for the parameter of interest. Default to \code{-Inf} if argument is missing.
#' @param ub upper bounds for the parameters of interest. Default to \code{Inf} if argument is missing.
#' @param ind indices of \code{cur} to update.
#' @param prior.fun log prior function, a function of the parameter to update.
#' @param lik.fun log-likelihood function, a function of the parameter to  For dependence parameters, updates of the dependence parameters can be passed as attributes.
#' @param ll value of the log-likelihood at \code{cur}.
#' @param pmu proposal mean; if missing, default to random walk.
#' @param pcov covariance matrix of the proposal for \code{cur}.
#' @param cond logical; should updates be made conditional on other values in \code{cur}. Default to \code{TRUE}.
#' @param ... additional arguments passed to the function, currently ignored.
#' @return a list with components
#' \itemize{
#' \item{\code{ll}}: value of the log-likelihood, with potential additional parameters passed as attributes;
#' \item{\code{cur}}: values of the parameters after update;
#' \item{\code{accept}}: logical indicating the proposal has been accepted (\code{TRUE}) or rejected (\code{FALSE}).
#' }
#' @export
mh.fun <- function(cur, lb, ub, prior.fun, lik.fun, ll, ind, pmu, pcov, cond = TRUE, ...){
  # overwrite missing values with default setting
  if(missing(ind)){
    lc <- length(cur)
    ind <- 1:lc
  } else{
    lc <- length(ind)
  }
  if(missing(lb)){
    lb <- rep(-Inf, length.out = lc)
  }
  if(missing(ub)){
    ub <- rep(Inf, length.out = lc)
  }
  #Sanity checks for length of arguments - only lengths
  pcov <- as.matrix(pcov) # handle scalar case
  stopifnot(length(lb) == lc, lc == length(ub), ncol(pcov) == nrow(pcov), ncol(pcov) == length(cur))
  # Copy value
  prop <- cur
  if(missing(pmu)){
   pmu <- cur
   RW <- TRUE
  } else{
   stopifnot(length(pmu) == length(cur))
  }

  if(!cond){
    sig <- as.matrix(pcov[ind, ind])
  # Sample new proposal from truncated Normal centered at current value
  prop[ind] <- TruncatedNormal::mvrandn(l = lb, mu = pmu[ind], u = ub, Sig = sig, n = 1)
  if(RW){
    # mean is at previous iteration, so density term drops out because of symmetry,
    #but not normalizing constant of truncated dist (since means differ)
    jac <- suppressWarnings(-log(TruncatedNormal::mvNcdf(l = lb - prop[ind], u = ub - prop[ind], Sig = sig, n = 1000)$prob) +
                              log(TruncatedNormal::mvNcdf(l = lb - cur[ind], u = ub - cur[ind], Sig = sig, n = 1000)$prob))
  } else{ #fixed mean, so normalizing constants for truncated components cancels out of the ratio of transition kernel densities
    jac <- mgp::dmvnorm(x = cur[ind], mean = pmu[ind], sigma = sig, logd = TRUE) -
      mgp::dmvnorm(x = prop[ind], mean = pmu[ind], sigma = sig, logd = TRUE)
  }
  } else{# cond == TRUE
   prop[ind] <- rcondmvtnorm(n = 1L, ind = ind, x = cur[-ind], lbound = lb, ubound = ub, mu = pmu, sigma = pcov, model = "norm")
   jac <- dcondmvtnorm(n = 1L, ind = ind, x = cur, lbound = lb, ubound = ub, mu = pmu, sigma = pcov, model = "norm", log = TRUE) -
     dcondmvtnorm(n = 1L, ind = ind, x = prop, lbound = lb, ubound = ub, mu = pmu, sigma = pcov, model = "norm", log = TRUE)

  }
  ll.p <- lik.fun(prop)
  if(missing(ll)){
    ll.c <- lik.fun(cur)
  } else{
    ll.c <- ll
  }
  prior.p <- prior.fun(prop)
  prior.c <- prior.fun(cur)
  acpt <- ll.p - ll.c + prior.p - prior.c + jac
  if(isTRUE(log(runif(1)) < acpt)){
    return(list(ll = ll.p, cur = prop, accept = TRUE))
  } else{
    return(list(ll = ll.c, cur = cur, accept = FALSE))
  }
}


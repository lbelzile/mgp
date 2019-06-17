#' Gaussian sampler using the Cholesky decomposition of the precision matrix
#'
#' @param n sample size
#' @param mu mean vector
#' @param precis precision matrix
#' @return sample of size \code{n} from the multivariate Gaussian distribution
#' @export
rmvnormprec <- function(n, mu = rep(0, ncol(precis)), precis) {
  stopifnot(ncol(precis) == nrow(precis), nrow(precis) == length(mu))
  if (length(mu) > 1) {
    t(mu + backsolve(chol(precis), x = diag(nrow(precis))) %*%
      matrix(rnorm(n * ncol(precis)), ncol = n, nrow = ncol(precis)))
  } else {
    mu + rnorm(n) / sqrt(precis)
  }
}

#' Multivariate normal density function
#'
#' This function returns the log-density for a multivariate Gaussian distribution.
#'
#' @param x matrix or vector of observations
#' @param mean mean vector
#' @param sigma positive definite covariance matrix
#' @param logd logical; whether log-density should be returned (default to \code{FALSE})
#' @return density or log-density of the \code{nrow(x)} sample
#' @export
dmvnorm <- function(x, mean, sigma, logd = FALSE) {
  if (any(missing(x), missing(mean), missing(sigma))) {
    stop("Arguments missing in function call to `dmvnorm`")
  }
  stopifnot(length(mean) == ncol(sigma), nrow(sigma) == ncol(sigma), is.logical(logd))
  if (is.vector(x)) {
    stopifnot(length(x) == length(mean))
    x <- matrix(x, nrow = 1, ncol = length(x))
  } else {
    stopifnot(ncol(x) == length(mean))
  }
  mgp::.dmvnorm_arma(x, as.vector(mean), as.matrix(sigma), logd = as.logical(logd))
}

#' Multivariate normal density function parametrized in terms of a precision matrix
#'
#' This function returns the log-density for a multivariate Gaussian distribution.
#'
#' @inheritParams dmvnorm
#' @param precis precision matrix
#' @return density or log-density of the \code{nrow(x)} sample
#' @export
dmvnorm.precis <- function(x, mean, precis, logd = FALSE) {
  if (is.vector(x)) {
    x <- t(as.matrix(x))
  }
  if (missing(mean)) {
    mean <- rep(0, ncol(x))
  }
  if (missing(precis)) {
    precis <- diag(0, ncol(x))
  }
  rooti <- t(chol(precis))
  quads <- colSums((crossprod(rooti, (t(x) - mean)))^2)
  if (logd) {
    return(-(ncol(x) / 2) * log(2 * pi) + sum(log(diag(rooti))) - .5 * quads)
  } else {
    return(exp(-(ncol(x) / 2) * log(2 * pi) + sum(log(diag(rooti))) - .5 * quads))
  }
}

#' Inverse gamma density, parametrized in terms of shape and rate
#'
#' @param x vector of observations
#' @param shape shape parameter, strictly positive
#' @param scale rate parameter of the gamma distribution, strictly positive
#' @param log logical; should the log-density be returned?
#' @export
dinvgamma <- function(x, shape, scale, log = TRUE) {
  stopifnot(any(c(shape > 0, scale > 0)))
  log.density <- shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - (scale / x)
  if (log) {
    return(log.density)
  } else {
    return(exp(log.density))
  }
}

#' Density of the multivariate Student density
#'
#' The function is adapted from \code{mvtnorm} package and computes the (log)-density
#' of the whole sample.
#'
#' @param x vector or matrix of observations
#' @param mu centrality parameter
#' @param sigma covariance matrix
#' @param df degrees of freedom parameter, default to 1 (Cauchy distribution)
#' @param logd logical; should log-density be returned? default to \code{TRUE}
#' @keywords internal
#' @export
dmvstud <- function(x, mu = rep(0, p), sigma = diag(p), df = 1, logd = TRUE) {
  if (is.vector(x)) {
    x <- matrix(x, ncol = length(x))
  }
  p <- ncol(x)
  stopifnot(length(mu) == ncol(x), df > 0)
  if (is.infinite(df)) {
    return(mgp::dmvnorm(x, mean = mu, sigma = sigma, logd = log))
  }
  if (!missing(sigma)) {
    if (p != ncol(sigma)) {
      stop("x and sigma have non-conforming size")
    }
  }
  dec <- try(chol(sigma))
  if (is.character(sigma)) {
    stop("Could not compute the Cholesky decomposition of the covariance matrix `sigma`")
  }
  R.x_m <- backsolve(dec, t(x) - mu, transpose = TRUE)
  rss <- colSums(R.x_m^2)
  logretval <- lgamma((p + df) / 2) - (lgamma(df / 2) + sum(log(diag(dec))) +
    p / 2 * log(pi * df)) - 0.5 * (df + p) * log1p(rss / df)
  return(ifelse(logd, sum(logretval), exp(sum(logretval))))
}

#' Gaussian sampler using the Cholesky decomposition of the precision matrix
#'
#' @param n sample size
#' @param mu mean vector
#' @param precis precision matrix
#' @return sample of size \code{n} from the multivariate Gaussian distribution
#' @export
rmvnormprec <- function(n, mu = rep(0, ncol(precis)), precis){
  stopifnot(ncol(precis) == nrow(precis), nrow(precis) == length(mu))
  if(length(mu) > 1){
    t(mu + backsolve(chol(precis), x = diag(nrow(precis))) %*%
        matrix(rnorm(n*ncol(precis)), ncol = n, nrow = ncol(precis)))
  } else {
    mu + rnorm(n)/sqrt(precis)
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
dmvnorm <- function(x, mean, sigma, logd = FALSE){
  if(any(missing(x), missing(mean), missing(sigma))){
    stop("Arguments missing in function call to `dmvnorm`")
  }
  stopifnot(length(mean) == ncol(sigma), nrow(sigma) == ncol(sigma), is.logical(logd))
  if(is.vector(x)){
   stopifnot(length(x) == length(mean))
   x <- matrix(x, nrow = 1, ncol = length(x))
  } else{
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
  if(is.vector(x)){
    x <- t(as.matrix(x))
  }
  if(missing(mean)){
    mean <- rep(0, ncol(x))
  }
  if(missing(precis)){
    precis <- diag(0, ncol(x))
  }
  rooti <- t(chol(precis))
  quads <- colSums((crossprod(rooti, (t(x) - mean)))^2)
  if(logd){
    return(-(ncol(x) / 2) * log(2 * pi) + sum(log(diag(rooti))) - .5 * quads)
  } else{
    return(exp(-(ncol(x) / 2) * log(2 * pi) + sum(log(diag(rooti))) - .5 * quads))
  }
}

#' Inverse gamma density, parametrized in terms of shape and rate
#'
#' @param x vector of observations
#' @param shape shape parameter, strictly positive
#' @param rate rate parameter, strictly positive
#' @param log logical; should the log-density be returned?
#' @export
dinvgamma <- function(x, alpha, beta, log = TRUE){
  stopifnot(any(c(alpha > 0, beta > 0)))
  log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - (beta/x)
  if(log){
    return(log.density)
  } else {
    return(exp(log.density))
  }
}

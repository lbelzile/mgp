#' Conditional distribution of Gaussian or Student subvectors
#'
#' The function computes the conditional distribution of the sub-components in \code{ind} given the rest and
#' returns the distribution function of the truncated Gaussian or Student components.
#' The location vector \code{mu} and the scale matrix \code{Sigma} are those of the \eqn{d+p} vector.
#' The routine relies on the CDF approximation based on minimax exponential tilting implemented in the \code{TruncatedNormal} package.
#'
#' @param ind a \code{d} vector of indices with integer entries in \eqn{\{1, \ldots, d+p\}} for which to compute the density
#' @param x a \code{p} vector with the values of process at remaining coordinates
#' @param lbound \code{d} vector of lower bounds for truncation
#' @param ubound \code{d} vector of upper bounds for truncation
#' @param mu \code{d+p} vector of location parameters
#' @param Sigma \code{d+p} by \code{d+p} scale matrix
#' @param df degrees of freedom of the \code{d+p} dimensional Student process
#' @param model string indicating family, either \code{norm} for Gaussian or \code{stud} for Student-t
#' @param n sample size for simulations. Default to 500.
#' @param log logical; should log probability be returned? Default to \code{FALSE}.
#' @return conditional distribution function for the components \code{ind} at \code{x[ind]}
#' @export
pcondmvtnorm <- function(n = 500, ind, x, lbound = rep(-Inf, length(ind)), ubound = rep(Inf, length(ind)),
                     mu, Sigma, df = NULL, model = c("norm", "stud"), log = FALSE){
  if(length(x) + length(ind) != length(mu)){
    stop("Invalid argument")
  }
  stopifnot(length(ind) > 0)
  x <- as.vector(x)
  mu <- as.vector(mu)
  schurcomp <- function(Sigma, ind){
    stopifnot(c(length(ind)>0, ncol(Sigma)-length(ind)>0))
    Sigma[ind, ind, drop = FALSE] - Sigma[ind, -ind, drop = FALSE] %*%
      solve(Sigma[-ind, -ind, drop = FALSE]) %*% Sigma[-ind, ind, drop = FALSE]
  }
  if(all(c(isTRUE(all.equal(lbound, rep(-Inf, length(ind)))),
           isTRUE(all.equal(ubound, rep(Inf, length(ind))))))){
    return(1)
  }
  if(length(x) == 0){ #Unconditional simulation
    proba <- switch(model,
             norm = TruncatedNormal::mvNcdf(n = n, l = lbound - mu, u = ubound - mu, Sig = Sigma)$prob,
             stud = TruncatedNormal::mvTcdf(n = n, l = lbound - mu, u = ubound - mu, Sig = Sigma, df = df)$prob
      )
  } else {
      muC <- c(mu[ind] + Sigma[ind, -ind, drop = FALSE] %*% solve(Sigma[-ind, -ind, drop = FALSE]) %*% (x - mu[-ind]))
      proba <- switch(model,
             norm = TruncatedNormal::mvNcdf(n = n, l = lbound - muC, u = ubound - muC,
                                               Sig = schurcomp(Sigma, ind))$prob,
             TruncatedNormal::mvTcdf(n = n, l = lbound - muC, u = ubound - muC,
                                        Sig = c(df + t(x- mu[-ind]) %*% solve(Sigma[-ind, -ind, drop = FALSE]) %*% (x- mu[-ind]))/
                                          (df + length(x)) * schurcomp(Sigma, ind),
                                        df = df + length(x))$prob)
  }
  if(log){
   log(proba)
  } else{
   proba
  }
}
#' Conditional density of Gaussian or Student subvectors
#'
#' The function computes the conditional density of (truncated) Gaussian or Student components corresponding to indices \code{ind},
#'  given the values at the remaining index.
#' The location vector \code{mu} and the scale matrix \code{Sigma} are those of the \eqn{d+p} vector.
#' The routine relies on the CDF approximation based on minimax exponential tilting implemented in the \code{TruncatedNormal} package.
#'
#' @param ind a \code{d} vector of indices to impute with integer entries in \eqn{\{1, \ldots, d+p\}}
#' @param x a \code{p} vector with the values of process at remaining coordinates
#' @param lbound \code{d} vector of lower bounds
#' @param ubound \code{d} vector of upper bounds for truncated
#' @param mu \code{d+p} vector of location parameters
#' @param Sigma \code{d+p} by \code{d+p} scale matrix
#' @param df degrees of freedom of the \code{d+p} dimensional Student process
#' @param model string indicating family, either \code{norm} for Gaussian or \code{stud} for Student-t
#' @param n sample size for simulations. Default to 500.
#' @param log logical; should log probability be returned? Default to \code{FALSE}.
#' @return the conditional (log) density of the vector \code{x[ind]}
#' @export
dcondmvtnorm <- function(x, ind, lbound = rep(-Inf, length(ind)), ubound = rep(Inf, length(ind)),
                         mu, Sigma, df = NULL, model = c("norm", "stud"), n = 500, log = FALSE){
  if(length(x) != length(mu)){
    stop("Invalid argument")
  }
  #stopifnot(length(ind) > 0)
  x <- as.vector(x)
  mu <- as.vector(mu)
  schurcomp <- function(Sigma, ind){
    stopifnot(c(length(ind)>0, ncol(Sigma)-length(ind)>0))
    Sigma[ind, ind, drop = FALSE] - Sigma[ind, -ind, drop = FALSE] %*%
      solve(Sigma[-ind, -ind, drop = FALSE]) %*% Sigma[-ind, ind, drop = FALSE]
  }
  if(length(x[ind]) == 0){ #Unconditional distribution function
    res <- switch(model,
           norm = mvtnorm::dmvnorm(x = x, mean = mu, sigma = as.matrix(Sigma), log = TRUE),
           stud = mvtnorm::dmvt(x = x - mu, sigma = Sigma, df = df, log = TRUE) #TODO check the centering

    )
  } else{
    muC <- c(mu[ind] + Sigma[ind, -ind, drop = FALSE] %*% solve(Sigma[-ind, -ind, drop = FALSE]) %*% (x[-ind] - mu[-ind]))
    sigmaC <- as.matrix(schurcomp(Sigma, ind))
    if(length(ind)==1L){
      res <- switch(model,
                    norm = mvtnorm::dmvnorm(x = x[-ind], mean = mu[-ind], sigma = solve(Sigma[-ind, -ind, drop = FALSE]), log = TRUE) +
                      TruncatedNormal:::lnNpr(a = (lbound - muC)/sqrt(sigmaC), b = (ubound - muC)/sqrt(sigmaC)),
                    stud = mvtnorm::dmvt(x = x[-ind]- mu[-ind], sigma = Sigma[-ind, -ind, drop = FALSE], df = df, log = TRUE) +
                      log(TruncatedNormal::mvTcdf(n = n, l = lbound - muC, u = ubound - muC,
                                                  Sig = c(df + t(x[-ind]- mu[-ind]) %*% solve(Sigma[-ind, -ind, drop = FALSE]) %*% (x[-ind]- mu[-ind]))/
                                                    (df + length(x[-ind])) * sigmaC,
                                                  df = df + length(x[-ind]))$prob))
    } else {
      res <- switch(model,
                    norm = mvtnorm::dmvnorm(x = x[-ind], mean = mu[-ind], sigma = Sigma[-ind,-ind, drop = FALSE], log = TRUE) +
                      log(TruncatedNormal::mvNcdf(n = n, l = lbound - muC, u = ubound - muC, Sig = sigmaC)$prob),
                    stud = mvtnorm::dmvt(x = x[-ind]- mu[-ind], sigma = Sigma[-ind, -ind, drop = FALSE], df = df, log = TRUE) +
                      log(TruncatedNormal::mvTcdf(n = n, l = lbound - muC, u = ubound - muC,
                                                  Sig = c(df + t(x[-ind]- mu[-ind]) %*% solve(Sigma[-ind, -ind, drop = FALSE]) %*% (x[-ind]- mu[-ind]))/
                                                    (df + length(x[-ind])) * sigmaC,
                                                  df = df + length(x[-ind]))$prob))
  }

  }
  return(ifelse(log, res, exp(res)))
}



#' Conditional samples of Gaussian or Student subvectors
#'
#' The function samples (truncated) Gaussian or Student vectors
#' corresponding to indices \code{ind} given the values at the remaining index.
#' The location vector \code{mu} and the scale matrix \code{Sigma} are those of the \eqn{d+p} vector.
#' The routine relies on the CDF approximation based on minimax exponential tilting implemented in the \code{TruncatedNormal} package.
#'
#' @param ind a \code{d} vector of indices to impute with integer entries in \eqn{\{1, \ldots, d+p\}}
#' @param x a \code{p} vector with the values of process at remaining coordinates
#' @param lbound \code{d} vector of lower bounds
#' @param ubound \code{d} vector of upper bounds for truncated
#' @param mu \code{d+p} vector of location parameters
#' @param Sigma \code{d+p} by \code{d+p} scale matrix
#' @param df degrees of freedom of the \code{d+p} dimensional Student process
#' @param model string indicating family, either \code{norm} for Gaussian or \code{stud} for Student-t
#' @param n sample size for the random vector; default to 1.
#' @return an n by d matrix of conditional simulations
#' @export
rcondmvtnorm <- function(n = 1L, ind, x, lbound = rep(-Inf, length(ind)), ubound = rep(Inf, length(ind)),
                         mu, Sigma, df = NULL, model = c("norm", "stud")){
  model <- match.arg(model)
  if(length(x) + length(ind) != length(mu)){
    stop("Invalid argument")
  }
  stopifnot(length(ind) > 0)
  x <- as.vector(x)
  mu <- as.vector(mu)
  schurcomp <- function(Sigma, ind){
    stopifnot(c(length(ind)>0, ncol(Sigma)-length(ind)>0))
    Sigma[ind, ind, drop=FALSE] - Sigma[ind,-ind, drop=FALSE] %*% solve(Sigma[-ind, -ind, drop=FALSE]) %*% Sigma[-ind,ind, drop=FALSE]
  }
 if(length(x) == 0){ #Unconditional simulation
    switch(model,
           norm = TruncatedNormal::mvrandn(n = n, mu = mu, l = lbound, u = ubound, Sig = Sigma),
           stud = TruncatedNormal::mvrandt(n = n, mu = mu, l = lbound, u = ubound, Sig = Sigma, df = df)
    )
  } else {
    muC <- c(mu[ind] + Sigma[ind, -ind, drop=FALSE] %*% solve(Sigma[-ind, -ind, drop=FALSE]) %*% (x - mu[-ind]))
    switch(model,
           norm = TruncatedNormal::mvrandn(n = n, mu = muC, l = lbound, u = ubound,
                                          Sig = schurcomp(Sigma, ind)),
           stud = TruncatedNormal::mvrandt(n = n, l = lbound, mu = muC, u = ubound,
                                   Sig = c(df + t(x- mu[-ind]) %*% solve(Sigma[-ind, -ind, drop=FALSE]) %*% (x- mu[-ind]))/
                                     (df + length(x)) * schurcomp(Sigma, ind),
                                   df = df + length(x))
    )

  }
}


#' Update range parameter via Metropolis-Hastings
#'
#' This function, adapted from \code{\link[spatial.gev.bma]{gev.update.lambda}},
#' uses the Laplace approximation to build a proposal for the range parameter of an
#' exponential covariance model with precision \eqn{\alpha} and range \eqn{\lambda}.
#'
#' The difference with the \code{spatial.gev.bma} package is that proposals
#' are truncated Gaussian, so may be accepted even if the Newton step is centered at negative values.
#' If the current \eqn{\lambda} value is in a region where the Hessian matrix is not negative, the proposal
#' is automatically rejected. The function \code{step} controls the size of the Newton steps.
#'
#' @param tau current value of centered random effects
#' @param alpha precision of the random effect
#' @param lambda range parameter to be updated
#' @param di matrix of pairwise distances between points
#' @param a shape parameter of the Gamma prior for \eqn{\lambda}
#' @param b rate parameter of the Gamma prior for \eqn{\lambda}
#' @param lb lower bound for admissible value of \eqn{\lambda}. Default to 0.01.
#' @param step size of step in Newton update; default to 0.2.
#' @param maxstep maximum step size for the truncated Normal proposal
#'
#' @return new value of \code{lambda}
#' @export
updaterange <- function (tau, alpha, lambda, di, a, b, lb = 1e-2, step = 0.2, maxstep = Inf){
  stopifnot(maxstep > 0)
  l.curr <- l.prime(tau, alpha, lambda, di, a, b)
  l.double.curr <- l.double.prime(tau, alpha, lambda, di, a, b)
  d.curr <- -l.double.curr
  b.curr <- step * l.curr / d.curr + lambda
if (isTRUE(d.curr > 0)) {
  lambda.new <- TruncatedNormal::mvrandn(pmax(lb, lambda - maxstep), lambda + maxstep, Sig = matrix(1/d.curr), n = 1, mu = b.curr)
  #rnorm(1, b.curr/d.curr, sd = sqrt(1/d.curr))
  l.new <- l.prime(tau, alpha, lambda.new, di, a, b)
  l.double.new <- l.double.prime(tau, alpha, lambda.new, di, a, b)
  d.new <- -l.double.new
  b.new <- step * l.new / d.new + lambda.new
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
      log(suppressWarnings(TruncatedNormal::mvNcdf(l = pmax(lb, lambda - maxstep) - b.curr, u = lambda + maxstep, Sig = 1/d.curr, n = 1)$prob))
    prop.new <- dnorm(lambda, b.new, sd = sqrt(1/d.new), log = TRUE) -
      log(suppressWarnings(TruncatedNormal::mvNcdf(l = pmax(lb, lambda.new - maxstep) - b.new, u = lambda.new + maxstep, Sig = 1/d.new, n = 1)$prob))
    mh <- L.new - L.curr + prior.new - prior.curr + prop.new - prop.curr
    if (isTRUE(log(runif(1)) < mh)) {
      lambda <- lambda.new
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

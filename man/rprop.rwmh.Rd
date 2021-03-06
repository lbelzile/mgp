% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{rprop.rwmh}
\alias{rprop.rwmh}
\title{Proposals for random walk Metropolis-Hastings}
\usage{
rprop.rwmh(cur = NULL, trcur = NULL, cov, lbound = rep(-Inf,
  ncol(cov)), ubound = rep(Inf, ncol(cov)))
}
\arguments{
\item{cur}{vector of current parameter on the original (constrained) scale}

\item{trcur}{vector of current parameter values on the unconstrained scale}

\item{cov}{covariance matrix of the proposal for the \code{trcur} scale}

\item{lbound}{lower bounds for the parameters}

\item{ubound}{upper bounds for the parameter}
}
\value{
a list with components \code{prop}, \code{trprop} and \code{logjac}
}
\description{
This function transforms a vector \code{cur} to an unconstrained scale based on \code{lbound} and \code{ubound},
then samples draws from a multivariate Gaussian vector with covariance matrix \code{cov} centered at the current (unconstrained) value.
The function then returns a list containing the proposal on the unconstrained scale (\code{trprop}),
on the original scale (\code{prop}) along with the log of the jacobian of the transformation (\code{logjac}).
}
\keyword{internal}

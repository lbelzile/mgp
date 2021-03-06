% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{updt.range}
\alias{updt.range}
\title{Update range parameter via Metropolis-Hastings}
\usage{
updt.range(tau, alpha, lambda, di, a, b, lb = 0.01, discount = 0.2,
  maxstep = Inf)
}
\arguments{
\item{tau}{current value of centered random effects}

\item{alpha}{precision of the random effect}

\item{lambda}{range parameter to be updated}

\item{di}{matrix of pairwise distances between points}

\item{a}{shape parameter of the Gamma prior for \eqn{\lambda}}

\item{b}{rate parameter of the Gamma prior for \eqn{\lambda}}

\item{lb}{lower bound for admissible value of \eqn{\lambda}. Default to 0.01.}

\item{discount}{size of discount in Newton update; default to 0.2.}

\item{maxstep}{maximum discount size for the truncated Normal proposal}
}
\value{
new value of \code{lambda}
}
\description{
This function, adapted from \code{\link[spatial.gev.bma]{gev.update.lambda}},
uses a Laplace approximation to build a proposal for the range parameter of an
exponential covariance model with precision \eqn{\alpha} and range \eqn{\lambda}.
}
\details{
The difference with the \code{spatial.gev.bma} package is that proposals
are truncated Gaussian,  so may be accepted even if the mean of the Newton step proposal is negative.
If the current \eqn{\lambda} value is in a region where the Hessian matrix is not negative, the proposal
is automatically rejected. The argument \code{discount} controls the size of the Newton steps.
}
\seealso{
\code{\link[spatial.gev.bma]{gev.update.lambda}}
}
\author{
code from spatial.gev.bma by Alex Lenkoski
}

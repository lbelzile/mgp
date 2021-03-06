% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{lscale.lgm}
\alias{lscale.lgm}
\title{Update for latent Gaussian model for the scale parameter of a generalized Pareto}
\usage{
lscale.lgm(scale, shape, ldat, lscale.mu, lscale.precis, lscale.tausq,
  mmax, discount = 1, maxstep = Inf)
}
\arguments{
\item{scale}{vector of scale parameters for the generalized Pareto distribution}

\item{shape}{vector of shape parameters for the generalized Pareto distribution}

\item{ldat}{list with exceedances at each site}

\item{lscale.mu}{mean of the log-Gaussian process for the scale}

\item{lscale.precis}{precision matrix of the log-Gaussian process corresponding to the inverse of the correlation matrix}

\item{lscale.tausq}{variance of the log-Gaussian process}

\item{mmax}{vector of maximum of each series in \code{ldat}}

\item{discount}{numeric giving the discount factor for the Newton  Default to 1.}

\item{maxstep}{maximum step size for the MCMC (proposal will be at most \code{maxstep} units away from the current value).}
}
\value{
a vector of scale parameters
}
\description{
The scale has a log-Gaussian prior with variance \code{lscale.tausq} and precision (the inverse of the correlation matrix)
given by \code{lscale.precis}.
}

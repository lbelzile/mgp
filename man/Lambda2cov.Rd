% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial.R
\name{Lambda2cov}
\alias{Lambda2cov}
\title{Transform variogram matrix to covariance of conditional random field}
\usage{
Lambda2cov(Lambda, co, subA, subB)
}
\arguments{
\item{Lambda}{Negative definite matrix for the Huesler--Reiss model}

\item{co}{vector of integer with conditioning sites}

\item{subA}{vector of integers with sub-entries (not in \code{co}) for rows}

\item{subB}{vector of integers with sub-entries (not in \code{co}) for columns. If missing, default to \code{subA}.}
}
\description{
The matrix \code{Lambda} is half the semivariogram matrix. The function
returns the conditional covariance with respect to entries in \code{co},
restricted to the \code{subA} rows and the \code{subB} columns.
}
\keyword{internal}

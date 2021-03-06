% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{transfo.pars}
\alias{transfo.pars}
\title{Transform variables to unconstrained space}
\usage{
transfo.pars(cur, lbound, ubound)
}
\arguments{
\item{cur}{vector of current parameter on the original (constrained) scale}

\item{lbound}{lower bounds for the parameters}

\item{ubound}{upper bounds for the parameter}
}
\value{
unconstrained transformed vector
}
\description{
This function uses the transformations detailed in the STAN reference manual (Section 56.4 in version 2.9).
}
\keyword{internal}

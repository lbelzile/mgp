#' Logistic link function
#'
#' Canonical link function for logistic regression with binomial data,
#' \eqn{\log(x)-\log(1-x)}.
#'
#' @param x probability
#' @return transformed probability
#' @export
#' @keywords internal
logit <- function(x) {
  stopifnot(isTRUE(all(x > 0, x < 1)))
  log(x) - log(1 - x)
}
#' Expit inverse link function for logistic regression
#'
#' The expit function is \eqn{\exp(x)/(1+\exp(x))}.
#'
#' @param x parameter
#' @return probability
#' @export
#' @keywords internal
expit <- function(x) {
  1 / (1 + exp(-x))
}

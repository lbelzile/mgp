#' Fit multivariate generalized Pareto models
#'
#' @param xdat \code{n} by \code{p} matrix of observations
#' @param threshold scalar threshold
#' @param mthresh marginal threshold
#' @param thresh functional threshold
#' @param model name of the model, either logistic (\code{"log"}), negative logistic (\code{"neglog"}), Huessler-Reiss (\code{"hr"}) or extremal student (\code{"xstud"})
#' @param method optimization routine, either sequential quadratic programming (\code{"sqp"}), PORT routines (\code{"nlminb"}) or \code{"BFGS"}
#' @param start named list of starting values; default to \code{NULL}
#' @param show logical; if \code{FALSE} (default), the result is returned as an invisible list
#' @param fpar named list of fixed parameters; default to \code{NULL}.
# fit.mgp <- function(xdat,
#                     mthresh,
#                     thresh = 0,
#                     model = c("log","neglog","hr","xstud"),
#                     method = c("sqp","nlminb", "BFGS"),
#                     start = NULL,
#                     show = FALSE,
#                     fpar = NULL,
#                     warnSE = FALSE,
#                     ...){
#   stopifnot("Missing data argument \"xdat\"." = !missing(xdat))
#   # Cast data frame or vector to matrix
#   xdat <- try(as.matrix(xdat))
#   stopifnot("\"xdat\" must be a numeric matrix." = is.matrix(xdat))
#   n <- nrow(xdat)
#   d <- ncol(xdat)
#   stopifnot("\xdat\" must be a matrix with more than one column." = p > 1)
#   # Remove entries with missing values for now
#   xdat <- na.omit(xdat)
#
#
#   # Deal with fixed parameters
#   param_names <- c("loc", "scale", "shape")
#   if(model %in% c("log", "neglog")){
#     param_names <- c(param_names, "dep")
#   # For now, only support exchangeable (equicorrelation)
#   } else if(model %in% c("br","xstud")){
#     param_names <- c(param_names, "cor")
#   }
#   if(model == "xstud"){
#     param_names <- c(param_names, "df")
#   }
#   # Check list elements
#   stopifnot(is.null(fpar) | is.list(fpar))
#   wf <- (param_names %in% names(fpar))
#
#   if(is.list(fpar) && (length(fpar) >= 1L)){ #non-empty list
#     if(is.null(names(fpar))){
#       stop("\"fpar\" must be a named list")
#     }
#     if(!isTRUE(all(names(fpar) %in% param_names))){
#       stop("Unknown fixed parameter: see help for more details")
#     }
#     if("loc" %in% names(fpar)){
#       if(!length(fpar$loc) %in% c(1L, d)){
#       stop("Location vector in \"fpar\" must be of length d or a scalar.")
#       }
#     }
#     if("scale" %in% names(fpar)){
#       if(!length(fpar$loc) %in% c(1L, d)){
#         stop("Scale vector in \"fpar\" must be of length d or a scalar.")
#       }
#     }
#     if("shape" %in% names(fpar)){
#       if(!length(fpar$loc) %in% c(1L, d)){
#         stop("Shape vectorin \"fpar\" must be of length d or a scalar.")
#       }
#     }
#     if("dep" %in% names(fpar)){
#       if(!(model %in% c("log","neglog")) | (length(fpar$dep) > 1L)){
#         stop("Invalid \"dep\" parameter provided in \"fpar\".")
#       }
#       if(fpar$dep < 0 | !is.finite(fpar$dep)){
#         stop("Invalid \"dep\" parameter provided: the parameter must be positive")
#       }
#       if(isTRUE(fpar$dep > 1) & model == "log"){
#         fpar$dep <- 1/fpar$dep
#       }
#     }
#     if("cor" %in% names(fpar)){
#       if(!(model %in% c("hr","xstud")) | (length(fpar$cor) > 1L)){
#         stop("Invalid \"cor\" parameter provided in \"fpar\".")
#       }
#       if(isTRUE(abs(fpar$cor) >= 1) | !is.finite(fpar$cor)){
#         stop("Invalid \"cor\" parameter provided: the value must be between -1 and 1.")
#       }
#     }
#     if("df" %in% names(fpar) & model == "xstud"){
#       if(length(fpar$df) > 1L | !isTRUE(fpar$df > 0) | !isTRUE(is.finite(fpar$df))){
#         stop("Invalid \"df\" parameter provided in \"fpar\".")
#       }
#     }
#   # TODO
#   # Step 1: Compute starting values (using fpar for marginal parameters, if any) for marginals
#   # using the NHPP likelihood, must deal with identifiability constraint (first location parameter?)
#   # Make sure to print this information in the output.
#   # Step 2: Create a wrapper for the log likelihood function with the fixed parameter list
#   # and an additional vector of parameters. Compute the length of the free parameter vector.
#   # Step 3: Create a vector function with the inequality constraint (factoring any fixed parameter!),
#   # including for the marginals.
#   # Step 4: Pass result to optim (for 1D case), alabama or Rsolnp (latter is not yet a dependence)
#   # Step 5: Check convergence and arrange result inside a list
#   # Step 6: Define print methods for the results
#   loglik <- function(par, model)
#   # Compute the constraints for each model
#   if(model == "log"){
#     stopifnot(length(par) == 1L)
#     if(par < 0 | par > 1){
#       return(-1e16)
#     }
#   } else if(model == "br"){
#     stopifnot(length(par) == 1L)
#     if(par < 0 | par > 1){
#       return(-1e16)
#     }
#     Lambda <- matrix(data = par,
#                      nrow = D,
#                      ncol = D)
#     diag(Lambda) <- 0
#   } else if(model == "xstud"){
#     stopifnot(length(par) == 2L)
#     if(any(par < 0) | par[2] > 1){
#       return(-1e16)
#     }
#     df <- par[1]
#     Sigma <- matrix(par[2],
#                     nrow = D,
#                     ncol = D)
#     diag(Sigma) <- 1
#   }
#   mev::clikmgp(
#     dat = apply(data, 1:2, qexp),
#     loc = 0,
#     scale = 1,
#     shape = 0,
#     mthresh = rep(qexp(mthresh), D),
#     thresh = qexp(thresh),
#     par = switch(model,
#                  log = list(alpha = par),
#                  br = list(Lambda = Lambda),
#                  xstud = list(df = df, Sigma = Sigma)),
#     model = model,
#     likt = "mgp",
#     lambdau = 1 - mthresh)
# }

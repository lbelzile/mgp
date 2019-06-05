#' Pairwise extremogram for max-risk functional
#'
#' The function computes the pairwise \eqn{chi} estimates and plots them as a function of the distance between sites.
#' @param dat data matrix
#' @param margp marginal probability above which to threshold observations
#' @param loc matrix of locations (one site per row)
#' @param scale geometric anisotropy scale parameter
#' @param rho geometric anisotropy angle parameter
#' @param plot logical; should a graph of the pairwise estimates against distance? Default to \code{FALSE}
#' @param ... additional arguments passed to plot
#' @return an invisible matrix with pairwise estimates of chi along with distance (unsorted)
#' @export
extremo <- function(dat, margp, loc, scale = 1, rho = 0, plot = FALSE, ...){
  stopifnot(isTRUE(all(margp >= 0, margp < 1, length(margp) == 1, nrow(loc) == ncol(dat))))
  # Local quantile - threshold data
  dat <- as.matrix(dat)
  margthresh <- apply(dat, 2, quantile, margp, na.rm = TRUE)
  dat <- t(t(dat)-margthresh)
  # Keep only instances where there is at least one exceedance
  dimat <- distg(loc = loc, scale = scale, rho = rho)
  excind <- apply(dat, 1, function(x){isTRUE(max(x, na.rm = TRUE) > 0)}) #avoid NA
  dat <- dat[excind,]
  res <- matrix(0, ncol = 4, nrow = choose(ncol(dat),2))
  b <- 0L
  for(i in 1:(ncol(dat)-1)){
    for(j in (i+1):ncol(dat)){
      b <- b + 1L
      subdat <- na.omit(dat[,c(i,j)])
      if(length(subdat) > 0){
        res[b, ] <- c(i,j, nrow(subdat), mean(I(subdat[,2] > 0) * I(subdat[,1] > 0))/(0.5*mean(I(subdat[,1] > 0)) + 0.5*mean(I(subdat[,2] > 0))) )
      } else{
       res[b, ] <- NA
      }
    }
  }
  res <- na.omit(res)
  ellips <- list(...)
  ellips$y <- res[,4]
  ellips$x <- apply(res, 1, function(x){dimat[x[1],x[2]]})
  if(plot){
    if(is.null(ellips$xlab)){ellips$xlab <- "distance"}
    if(is.null(ellips$ylab)){ellips$ylab <- "conditional probability of exceedance"}
    if(is.null(ellips$pch)){ ellips$pch <- 20}
    if(is.null(ellips$yaxs)){ ellips$yaxs <- "i"}
    if(is.null(ellips$xlim)){ ellips$xlim <- c(0, max(ellips$x)*1.02)}
    if(is.null(ellips$xaxs)){ ellips$xaxs <- "i"}
    if(is.null(ellips$bty)){ ellips$bty <- "l"}
    if(is.null(ellips$ylim)){ ellips$ylim <- c(0,1)}
    if(is.null(ellips$pch)){ ellips$pch <- 20}
    if(is.null(ellips$col)){ ellips$col <- grDevices::rgb(0, 0, 0, alpha = 0.25)}
  do.call(what = graphics::plot, args = ellips)
  }
  return(invisible(cbind(dist = ellips$x, prob = ellips$y)))
}

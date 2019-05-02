#' Data augmentation for the Brown-Resnick process
#'
#' This function returns a matrix of data on the data scale (if \code{map == TRUE}) or else unit Pareto observations.
#' Data that fall below a marginal threshold are censored and imputed according to the conditional distribution.
#' For \code{riskr == "sum"}, this require an accept-reject scheme to ensure the points with the simulated components
#' still lies in the risk region defined by the marginal and dependence parameters.
#'
#' @inheritParams clikmgp
#' @param riskr string giving the risk region, either \code{max} or \code{sum}
#' @param map logical; should imputed data be returned on the unit Pareto scale? Default to \code{FALSE}
#' @return a matrix of observations with imputed values in place of censored components.
#' @export
impute <- function(dat, thresh, mthresh, loc = 1, scale = 1, shape = 1, lambdau = 1, riskr = c("max", "sum"), par, map = FALSE, ...){
  ellips <- list(...)
  model <- "br" #only model currently supported
  riskr <- match.arg(riskr)
  # Set dimensions and replicate marginal GP parameters
  N <- nrow(dat); D <- ncol(dat)
  shape <- rep(shape, length.out = D)
  scale <- rep(scale, length.out = D)
  loc  <- rep(loc, length.out = D)
  lambdau <- rep(lambdau, length.out = D)
  if(model == "br"){
  Lambda <- par$Lambda
   if(is.null(Lambda)){
     if(!is.null(ellips$Lambda)){
       Lambda <- ellips$Lambda
       } else{
        stop("User must provide dependence matrix `Lambda`")
      }
    }
  }
  #Copy the matrix and transform data to unit Pareto scale (tdat)
  tdat <- dat
  colmax <- apply(dat, 2, max)
  if(isTRUE(any(colmax > ifelse(shape < 0, loc -scale/shape, Inf)))){
    stop("Some exceedances are above the upper bound of the distribution.")
  }
  #Marginal transformation from GP to unit Pareto
  for(j in 1:D){
    tdat[,j] <-  pmax(0,(1+shape[j]/scale[j]*(dat[,j] - loc[j])))^(1/shape[j])/lambdau[j]
  }
  #Also map the thresholds
  mt <- pmax(0,(1+shape/scale*(mthresh-loc)))^(1/shape)/lambdau
  # Marginal threshold must be such that any component is resolvable
  if(riskr == "sum"){
  stopifnot(isTRUE(all(mthresh <= thresh/D)))
  } else if(riskr == "max"){
    stopifnot(isTRUE(all(mthresh <= thresh)))
  }
  if(is.null(ellips$censored)){
    censored <- t(t(dat) < mthresh)
  } else{
   censored <- ellips$censored
  }
  if(is.null(ellips$numAbovePerRow)){
    numAbovePerRow <- D - rowSums(censored)
  } else{
    numAbovePerRow <- ellips$numAbovePerRow
  }

  #if this is true, only one conditional simulation is needed
  # partSumAbove <- rowSums(dat * I(!censored) ) >= u
  stopifnot(dim(tdat) == dim(censored), length(numAbovePerRow) == N)
  # if(model == "xstud"){
  #   for(i in 1:N){
  #     if(numAbovePerRow[i] < D){
  #       ind = which(censored[i,])
  #       Zin <- Z[i,!censored[i,]]^(1/df)
  #       k <- D - length(ind)
  #       kst <- c(Zin %*% solve(Sigma[-ind, -ind]) %*% Zin) / (k + df)
  #       muC <-  c(Sigma[ind, -ind]  %*% solve(Sigma[-ind, -ind]) %*% Zin) #TODO check this
  #       SigC <- kst * schurcompC(Sigma,  which(!censored[i,]))
  #       inriskregion <- FALSE
  #       while(!inriskregion){
  #         prop <- as.matrix(TruncatedNormal::mvrandt(u = mt[ind]^(1/df), mu = muC, l = rep(-Inf, length(ind)),
  #                                                    Sig =  SigC, df = df + k, n = ifelse(partSumAbove[i], 1, 10)))
  #         for(k in 1:nrow(prop)){
  #           prop <- sign(prop) * abs(prop)^df
  #           Z[i,censored[i,]] <-  prop[k,]
  #           inriskregion <- sum(scale*(pmax(0,Z[,i])^shape-1)/shape + loc) > u
  #           if(inriskregion){ break }
  #         }
  #       }
  #     }
  #   }
  # } else if (model == "br"){ #parametrization with lambda
  nsimu <- switch(riskr, max = 1, sum = 10)
    for(i in 1:N){
      if(numAbovePerRow[i] < D){
        be <- which(censored[i,])
        ab <- which(!censored[i,])
        if(length(ab) < 1){ stop("Non censored component in input")}
        #remove first non-censored, shift indices by 1 - because inference is conditional on one site, so this is a degenerate
        #log-Gaussian process (its value is almost surely 1 at the projection site
        be2 <- be - I(be > ab[1])
        ab2 <- ab[-1] - I(ab[-1] > ab[1])
        SigmaD <- outer(2 * Lambda[ab[1], - ab[1]], 2 * Lambda[ab[1], - ab[1]], "+") - 2 * Lambda[-ab[1], - ab[1]]
        muD <- - 2 * Lambda[ab[1], - ab[1]]
        thD <- log(mt[be]) - log(Z[i, ab[1]])
        inriskregion <- FALSE
        while(!inriskregion){
          if(length(ab) == 1){
            if(length(be) > 1){
              prop <- as.matrix(TruncatedNormal::mvrandn(n = nsimu,
                                                         l = rep(-Inf, length(be2)),
                                                         u = thD, mu = muD, Sig = SigmaD))
            } else { #length be == 1
              prop <- as.matrix(TruncatedNormal::mvrandn(n = 1, u = thD, mu = muD,
                                                         Sig = SigmaD,l =  - log(tdat[i, ab[1]]) +
                        ifelse(partSumAbove[i], -Inf, log((1+shape[be]/scale[be]*(thresh - sum(dat[i,-be]) -loc[be]))^(1/shape[be])/lambdau[be]))))

            }
          } else {# length(ab) > 1
            if(length(be) == 1){
              prop <- t(as.matrix(rcondmvtnorm(n = 1, ind = be2, model = "norm",
                         ubound = thD, mu = muD, Sigma = SigmaD, x = log(Z[i, ab2])-log(tdat[i, ab[1]]),
                         lbound =  - log(tdat[i, ab[1]]) + ifelse(partSumAbove[i], -Inf,
                         log((1+shape[be]/scale[be]*(thresh - sum(dat[i,-be] - loc[be])))^(1/shape[be])/lambdau[be])))))
            } else{
              prop <- t(as.matrix(rcondmvtnorm(n = switch(riskr, max = 1, sum = 5*length(be2)),
                                               ind = be2, x = log(tdat[i, ab2])-log(tdat[i, ab[1]]),
                                               ubound = thD, mu = muD, Sigma = SigmaD, model = "norm")))
            }
          }
          for(k in 1:nrow(prop)){
            tdat[i,censored[i,]] <- exp(log(tdat[i, ab[1]]) + prop[k,])
            if(riskr == "max"){
              inriskregion <- TRUE
            } else if(riskr == "sum"){
            inriskregion <- isTRUE(sum(scale*(tdat[i,]^shape-1)/shape + loc) > u)
            }
            if(inriskregion){break()}
          }
        }
      }
    }
  if(!map){
    #Back transform the observations
    for(j in 1:D){
      dat[censored[,j],j] <-  (tdat[censored[,j],j]^shape[j]-1)/shape[j]*scale[j] + loc[j]
    }
    return(dat)
  } else{
    #Return all observations on unit Pareto scale
    return(tdat)
  }
}

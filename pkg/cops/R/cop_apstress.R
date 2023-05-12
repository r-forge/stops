#' PCOPS version of approximated power stress model.
#'
#' This uses an approximation to power stress that makes use of smacofx as workhorse. Free parameters are kappa, lambda and nu
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of parameters to optimize over. Must be of length three, with the first the kappa argument, the second the lambda argument and the third the nu argument. One cannot supply upsilon and tau as of yet. Defaults to 1 1 1.
#' @param type MDS type. 
#' @param ndim number of dimensions of the target space
#' @param itmaxi number of iterations. default is 1000.
#' @param weightmat (optional) a binary matrix of nonnegative weights.
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param cordweight weight to be used for the cordillera; defaults to 0.5
#' @param q the norm of the cordillera; defaults to 1
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to ndim+1
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the distances (min distance minus max distance). If NULL (default) the cordillera will be normed to each configuration's maximum distance, so an absolute value of goodness-of-clusteredness.
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param normed should the cordillera be normed; defaults to TRUE
#' @param scale should the configuration be scale adjusted
#' @param ... additional arguments to be passed to the fitting procedure
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress-1 of the configuration
#'         \item{stress.m:} default normalized stress (sqrt(stress-1)) 
#'         \item{copstress:} the weighted loss value
#'         \item{OC:} the OPTICS cordillera value
#'         \item{parameters:} the theta parameters used for fitting (kappa, lambda, nu)
#'         \item{fit:} the returned object of the fitting procedure (typically of class smacofB or smacofP)
#'         \item{cordillera:} the cordillera object
#' }
#'
#'@importFrom stats dist as.dist
#'@importFrom smacofx apStressMin 
#'@import cordillera 
#'@import smacof
#' 
#'@keywords multivariate
cop_apstress <- function(dis,theta=c(1,1,1),type="ratio",ndim=2,weightmat=1-diag(nrow(dis)),init=NULL,itmaxi=1000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,normed=TRUE,scale="sd") {
  #if(all.equal(theta,c(1,1,1))) theta <- 1
  #TODO Unfolding
  if(inherits(dis,"dist") || is.data.frame(dis)) dis <- as.matrix(dis)
  if(length(setdiff(unique(unlist(as.vector(weightmat))),c(0,1)))>0) stop("For approximated power stress, only binary weight matrices are allowed.")   
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)<3) theta <- rep(theta,length.out=3)
    kappa <- theta[1]
    lambda <- theta[2]
    nu <- theta[3]
    fit <- smacofx::apStressMin(dis, kappa=kappa, lambda=lambda, nu=nu, type=type, ndim=ndim, weightmat=weightmat, init=init, verbose=isTRUE(verbose==2), itmax=itmaxi,...) #optimize with smacofx::apStressMin
    #fit$stress.1 <- fit$stress
    fit$stress.m <- fit$stress^2
    #fit$deltaorig <-stats::as.dist(dis)
    copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
    out <- list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit, copsobj=copobj) #target functions
   out
}

#' STOPS version of rstress
#'
#' Free parameter is kappa=2r for the fitted distances.
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this must be a scalar of the kappa=2*r transformation for the fitted distances proximities. Defaults to 1. Note that what is returned is r, not kappa.
#' @param type MDS type. Default is "ratio"
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param itmaxi number of iterations.
#' @param ... additional arguments to be passed to the fitting procedure
#' @param structures which structuredness indices to be included in the loss
#' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
#' @param strucpars the parameters for the structuredness indices
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param stoptype How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
#' 
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress-1 value
#'         \item{stress.m:} default normalized stress
#'         \item{stoploss:} the weighted loss value
#'         \item{indices:} the values of the structuredness indices
#'         \item{parameters:} the parameters used for fitting 
#'         \item{fit:} the returned object of the fitting procedure
#'         \item{stopobj:} the stopobj object
#' }
#'
#' @import cordillera
#' @keywords multivariate
#' @importFrom smacofx rStressMin
#' @export
stop_rstress <- function(dis,theta=1,type="ratio",weightmat=NULL,init=NULL,ndim=2,itmaxi=10000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative")) {
  if(inherits(dis,"dist")) dis <- as.matrix(dis)   
  if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  #itmaxi2 <- itmaxi
  theta <- as.numeric(theta)
  #if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(missing(stoptype)) stoptype <- "additive"
  weightmat <- as.matrix(weightmat)
  #if(length(theta)<3) theta <- rep(theta,length.out=3)
  r <- theta[1]/2
  fit <- smacofx::rStressMin(delta=dis,r=r,type=type,weightmat=weightmat,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  ncall <- do.call(substitute,list(fit$call,list(r=r,type=type,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi)))
  fit$call <- ncall                 
  fit$r <- r
  #fit$lambda <- 1
  #fit$nu <- 1
  fit$parameters <- fit$theta  <- fit$pars <- c(r=r)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
  out
}

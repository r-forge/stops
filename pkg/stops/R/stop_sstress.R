#' STOPS version of sstress
#'
#' Free parameter is lambda for the observed proximities. Fitted distances are transformed with power 2, weights have exponent of 1. Note that the lambda here works as a multiplicator of 2 (as sstress has f(delta^2)).
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this must be a scalar of the lambda transformation for the observed proximities. Defaults to 1. Note that the lambda here works as a multiplicator of 2 (as sstress has f(delta^2)). 
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param type MDS type.
#' @param init (optional) initial configuration
#' @param ndim the number of dimensions of the target space
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param itmaxi number of iterations
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
#'         \item{parameters:} the parameters used for fitting (lambda)
#'         \item{fit:} the returned object of the fitting procedure
#'          \item{stopobj:} the stopobj object
#' }
#' @import cordillera
#' @importFrom smacofx powerStressMin
#' @keywords multivariate
#' @export
stop_sstress <- function(dis,theta=1,type=type,weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,itmaxi=100000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)  
  if(missing(stoptype)) stoptype <- "additive"
#  if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  #if(length(theta)<3) theta <- rep(theta, length.out=3)
  lambda <- theta[1]
  flambda <- lambda*2 #sstress is d^2 and delta^2 so f(delta^2)=delta^(2*1); lambda works in factors of 2  
  fit <- smacofx::powerStressMin(delta=dis,type=type,kappa=2,lambda=flambda,nu=1,weightmat=weightmat,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  fit$kappa <- 2
  fit$lambda <- flambda
  #fit$nu <- 1
  fit$parameters <- fit$theta <- fit$pars  <- c(lambda=fit$lambda,kappa=2)#c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu) 
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
  out
}

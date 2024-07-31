#' STOPS versions of smacofSphere models
#'
#' The free parameter is lambda for power transformations the observed proximities. The fitted distances power is internally fixed to 1 and the power for the weights is 1. 
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this must be a scalar of the lambda transformation for the observed proximities. Defaults to 1.
#' @param type MDS type.
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param itmaxi number of iterations
#' @param ... additional arguments to be passed to the fitting procedure
#' @param structures which structuredness indices to be included in the loss
#' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
#' @param strucpars the parameters for the structuredness indices
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param stoptype How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative'
#' @param registry registry object with c-structuredness indices.
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
#'          \item{stopobj:} the stopobj object
#' }
#'
#'@importFrom smacof smacofSphere 
#'@importFrom stats dist as.dist
#'@keywords multivariate
#'@export
stop_smacofSphere <- function(dis,theta=1,type="ratio",ndim=2,weightmat=NULL,init=NULL,itmaxi=1000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative"),registry=struc_reg) {
  if(is.null(init)) init <- "torgerson"
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(missing(stoptype)) stoptype <- "additive"
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
  #kappa first argument, lambda=second
  if(length(theta)>1) stop("There are too many parameters in the theta argument.")
  #if(length(theta)<3) theta <- rep(theta,length.out=3)
  lambda <- theta[1]
  fit <- smacof::smacofSphere(dis^lambda,type=type,ndim=ndim,weightmat=weightmat,init=init,verbose=isTRUE(verbose==2),itmax=itmaxi,...) #optimize with smacof
  #fit$kappa <- 1
  fit$lambda <- lambda
  #fit$nu <- 1
 # fit$stress.1 <- fit$stress
 # fitdis <- as.matrix(fit$confdist)
 # delts <- as.matrix(fit$delta)[-1,-1]
 # fit$obsdiss <- fit$dhat 
 # fit$stress.r <- sum(weightmat*(delts-fitdis)^2)
  fit$stress.m <- fit$stress^2
  fit$parameters <- fit$theta <- fit$pars <- c(lambda=lambda)#c(kappa=fit$kappa,lambda=fit$lambda,rho=fit$nu)
  fit$deltaorig <- stats::as.dist(dis) 
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype,registry=registry)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit,stopobj=stopobj) #target functions
  out
}

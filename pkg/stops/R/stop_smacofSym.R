#' STOPS version of smacofSym models
#'
#' The free parameter is lambda for power transformations the observed proximities. The fitted distances power is internally fixed to 1 and the power for the weights is 1. 
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector; must be a scalar for the lambda (proximity) transformation. Defaults to 1.
#' @param type MDS type. Defaults ot 'ratio'.
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param itmaxi number of iterations
#' @param ... additional arguments to be passed to the fitting
#' @param structures which structuredness indices to be included in the loss
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
#' @param strucpars the parameters for the structuredness indices
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param stoptype How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
#' @param registry registry object with c-structuredness indices.
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress-1 (sqrt(stress.m))
#'         \item{stress.m:} default normalized stress (used for STOPS)
#'         \item{stoploss:} the weighted loss value
#'         \item{indices:} the values of the structuredness indices
#'         \item{parameters:} the parameters used for fitting (lambda) 
#'         \item{fit:} the returned object of the fitting procedure
#'         \item{stopobj:} the stops object
#' }
#' 
#'@keywords multivariate
#'@import smacof
#'@importFrom stats as.dist
#'@export
stop_smacofSym <- function(dis, theta=1, type="ratio", ndim=2,weightmat=1-diag(nrow(dis)),init=NULL,itmaxi=1000,...,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"),stressweight=1,strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative"),registry=struc_reg) {
  theta <- as.numeric(theta)
  if(is.null(init)) init <- "torgerson"  
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  #if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  if(missing(stoptype)) stoptype <- "additive"
  #if(type=="ordinal") warning("Ordinal MDS is invariant to monotonic transformations of the dissimilarities.")
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  #if(length(theta)==1) lambda <- theta
  #if(length(theta)==2) lambda <- theta[2]
  #if(length(theta)==3) lambda <- theta[2]
  lambda <- theta[1]
  fit <- smacof::smacofSym(dis^lambda,type=type,ndim=ndim,weightmat=weightmat,init=init,verbose=isTRUE(verbose==2),itmax=itmaxi,...) #optimize with smacof
  #fit$kappa <- 1
  fit$lambda <- lambda
  #fit$nu <- 1
  #fit$stress.1 <- fit$stress
  #fitdis <- as.matrix(fit$confdist)
  #delts <- as.matrix(fit$dhat) 
  #fit$stress.r <- sum(weightmat*(delts-fitdis)^2)
  fit$stress.m <- fit$stress^2 #fit$stress.r/sum(weightmat*delts^2)
  #fit$pars <- c(lambda=fit$lambda)#c(kappa=fit$kappa,lambda=fit$lambda,rho=fit$nu)
  fit$parameters <- fit$theta <- fit$pars  <- c(lambda=lambda)
  fit$deltaorig <- stats::as.dist(dis)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype,registry=registry)
  out <- list(stress=fit$stress, stress.r=fit$stress.r,stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices,parameters=stopobj$parameters,fit=fit,stopobj=stopobj) #target functions
  out
}


#' STOPS version of approximated power stress models.
#'
#' This uses an approximation to power stress that can make use of smacof as workhorse. Free parameters are kappa, lambda and nu. 
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of parameters to optimize over. Must be of length three, with the first the kappa argument, the second the lambda argument and the third the nu argument. One cannot supply upsilon and tau as of yet. Defaults to 1 1 1.
#' @param type MDS type.
#' @param ndim number of dimensions of the target space
#' @param itmaxi number of iterations. default is 1000.
#' @param weightmat (optional) a binary matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param ... additional arguments to be passed to the fitting procedure
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param structures a character vector listing the structure indices to use. They always are called "cfoo" with foo being the structure.
#' @param strucweight weight to be used for the structures; defaults to 1/number of structures
#' @param strucpars a list of list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appearance in structures vector. See examples.  
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param stoptype which weighting to be used in the multi-objective optimization? Either 'additive' (default) or 'multiplicative'.
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress-1 value (sqrt stress.m)
#'         \item{stress.m:} default normalized stress
#'          \item stoploss: the weighted loss value
#'         \item struc: the structuredness indices
#'         \item parameters: the parameters used for fitting (kappa, lambda, nu)
#'         \item fit: the returned object of the fitting procedure
#'         \item{stopobj:} the stopobj object 
#' }
#'
#'@importFrom stats dist as.dist
#'@importFrom smacofx apStressMin 
#'@import smacof
#' 
#'@keywords multivariate
stop_apstress <- function(dis,theta=c(1,1,1),type="ratio",ndim=2,weightmat= 1-diag(nrow(dis)),init=NULL,itmaxi=1000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative")) {
  if(inherits(dis,"dist") || is.data.frame(dis)) dis <- as.matrix(dis)
  if(missing(stoptype)) stoptype <- "additive"
  if(length(setdiff(unique(unlist(as.vector(weightmat))),c(0,1)))>0) stop("For approximated power stress, only binary weight matrices are allowed.")  
  if(length(setdiff(unique(unlist(as.vector(weightmat))),c(0,1)))>0) stop("For approximated power stress, only binary weight matrices are allowed.")  
  #we allow for theta to be of length three for compatibility in stops; maybe change that in the future 
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)<3) theta <- rep(theta,length.out=3)
  #if(length(theta)==2L) theta <- c(1,theta) 
  kappa <- theta[1]
  lambda <- theta[2]
  nu <- theta[3]
  verbose <- isTRUE(verbose>=2)
  fit <- smacofx::apStressMin(dis, kappa=kappa, lambda=lambda, nu=nu, type=type, ndim=ndim, init=init, verbose=verbose, itmax=itmaxi,...) #optimize with smacof
  ncall <- do.call(substitute,list(fit$call,list(kappa=kappa,lambda=lambda,nu=nu,type=type,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi)))
  fit$call <- ncall                
  #fit$kappa <- 1
  #fit$tau <- tau
  #fit$upsilon <- ups
  #fit$stress.1 <- fit$stress #smacof stress is sqrt(stress.m); for compatibility with powerStressMin we use stress^2 as stress.m 
  #fitdis <- as.matrix(fit$confdist)
  #delts <- as.matrix(fit$delta) 
  #fit$stress.r <- sum(combwght*(delts-fitdis)^2)
  fit$stress.m <- fit$stress^2 
  #fit$pars <- c(tau=fit$tau,upsilon=fit$upsilon) #c(kappa=fit$kappa,tau=fit$tau,upsilon=fit$upsilon)
  #fit$deltaorig <- fit$delta^(1/fit$tau)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
  out 
}

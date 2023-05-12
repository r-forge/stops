#' Another STOPS version of Sammon mapping models (via smacofSym)
#'
#' Uses Smacof, so it can deal with a weight matrix too.  The free parameter is lambda for power transformations of the observed proximities. The fitted distances power is internally fixed to 1 and the power for the weights=delta is -1. 
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this must be a scalar of the lambda transformation for the observed proximities. Defaults to 1.
#' @param type MDS type
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
#' @param stoptype How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative'. 
#' 
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress-1 (sqrt(stress.m))
#'         \item{stress.m:} default normalized stress (used for STOPS)
#'         \item{stoploss:} the weighted loss value
#'         \item{indices:} the values of the structuredness indices
#'         \item{parameters:} the parameters used for fitting (lambda)
#'         \item{fit:} the returned object of the fitting procedure
#'          \item{stopobj:} the stopobj object
#' }
#'
#'
#' @importFrom stats as.dist
#' @importFrom smacof smacofSym
#' @import cordillera
#' 
#'@keywords multivariate
#'@export
stop_sammon2 <- function(dis,theta=1,type="ratio",ndim=2,weightmat=NULL,init=NULL,itmaxi=1000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(is.null(init)) init <- "torgerson" 
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(missing(stoptype)) stoptype <- "additive"
  if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis)) 
  #kappa first argument, lambda=second
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  #if(length(theta)<3) theta <- rep(theta, length.out=3)
  lambda <- theta[1]
  nu <- -1
  elscalw <- dis^(nu*lambda) #the weighting in 
  diag(elscalw) <- 1
  combwght <- stats::as.dist(weightmat*elscalw) #combine the user weights and the elastic scaling weights
  fit <- smacof::smacofSym(dis^lambda,type=type,ndim=ndim,weightmat=combwght,init=init,verbose=isTRUE(verbose==2),itmax=itmaxi,...) #optimize with smacof
  #fit$kappa <- 1
  fit$lambda <- lambda
  #fit$nu <- nu
  #fit$stress.1 <- fit$stress
  #fitdis <- as.matrix(fit$confdist)
  #delts <- as.matrix(fit$delta)
  #fit$obsdiss <- fit$dhat
  #fit$stress.r <- sum(combwght*((delts-fitdis)^2))
  fit$stress.m <- fit$stress^2# fit$stress.r/sum(combwght*delts^2)
  fit$parameters <- fit$theta <- fit$pars <- c(lambda=lambda)#c(kappa=fit$kappa,lambda=fit$lambda,rho=fit$nu)
  fit$deltaorig <- stats::as.dist(dis) 
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit,stopobj=stopobj) #target functions
  out
}


#' STOPS version of strain
#'
#' The free parameter is lambda for power transformations of the observed proximities.
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this must be a scalar of the lambda transformation for the observed proximities.
#' @param type MDS type. Ignored here. 
#' @param ndim number of dimensions of the target space
#' @param itmaxi number of iterations. No effect here.
#' @param add should the dissimilarities be made Euclidean? Defaults to TRUE.
#' @param weightmat (optional) a matrix of nonnegative weights. Not used. 
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param ... additional arguments to be passed to the fitting procedure
#' @param structures which structuredness indices to be included in the loss
#' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
#' @param strucpars the parameters for the structuredness indices
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param stoptype How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative'
#' @param add if TRUE dis is made to Euclidean distances
#' @param registry registry object with c-structuredness indices.
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the badness-of-fit value (this isn't stress here but 1-(sum_ndim(max(eigenvalues,0))/sum_n(max(eigenvalues,0)), 1-GOF[2])
#'         \item{stress.m:} explictly normalized stress (manually calculated)
#'         \item{stoploss:} the weighted loss value
#'         \item{indices:} the values of the structuredness indices
#'         \item{parameters:} the parameters used for fitting (lambda)
#'         \item{fit:} the returned object of the fitting procedure, which is cmdscalex object with some extra slots for the parameters and stresses
#'          \item{stopobj:} the stopobj object
#' }
#'
#' @import cordillera
#' @importFrom stats dist as.dist
#' @importFrom smacofx cmdscale
#' @keywords multivariate
stop_cmdscale <- function(dis,theta=1,type="ratio",weightmat=NULL,ndim=2,init=NULL,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative"),itmaxi=1000,add=TRUE,registry=struc_reg) {
  theta <- as.numeric(theta)
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  type <- match.arg(type,'ratio')
  if(missing(stoptype)) stoptype <- "additive"
  lambda <- theta[1]
  #if(length(theta)==1) lambda <- theta
  #if(length(theta)==2) lambda <- theta[2]
  #if(length(theta)==3) lambda <- theta[2]
  if(missing(add)) add <- TRUE
  fit <- smacofx::cmdscale(dis^lambda,k=ndim,eig=TRUE,add=add,...) 
  fit$lambda <- lambda
  #fit$kappa <- 1
  fit$delta <- stats::as.dist(dis)
  dhat <- stats::as.dist(fit$dhat)
  fitdis <- stats::as.dist(fit$confdist)
  stress.r <- sum((dhat-fitdis)^2)
  stress.n <- stress.r/sum(dhat^2)
  fit$stress <- sqrt(stress.n)
  fit$stress.m <- stress.n
  fit$parameters <- fit$pars <- fit$theta <- c(lambda=lambda)#c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
 #fit$deltaorig <- stats::as.dist(dis)                                      #fit$nu <- 1
 # dis <- stats::as.dist(dis)
 # fitdis <- stats::dist(fit$points)
 # fit$stress.r <- sum((dis^lambda-fitdis)^2)
 # fit$stress.n <- fit$stress.r/sum(dis^(2*lambda))
 # fit$stress <- sqrt(fit$stress.n)
 # fit$stress.m <- fit$stress.n
  #fit$pars <- c(lambda=fit$lambda)#c(kappa=fit$kappa,lambda=fit$lambda,rho=fit$nu)
  #fit$conf <- fit$points
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype,registry=registry)
  list(stress=1-fit$GOF[2],stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj) #target functions
}


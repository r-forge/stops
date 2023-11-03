
#' STOPS version of CLDA with free k.
#'
#' CLDA with free lambda0 and k and 20 epochs. Should we add alpha0?
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of explicit parameters; first is lambda0 for the maximal neighbourhood and second is k for the number of neighbours for the geodesic distance. 
#' @param type MDS type. 
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param ndim number of dimensions of the target space
#' @param itmaxi number of iterations
#' @param ... additional arguments to be passed to the fitting procedure
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param structures a character vector listing the structure indices to use. They always are called "cfoo" with foo being the structure.
#' @param strucweight weight to be used for the structures; defaults to 1/number of structures
#' @param strucpars a list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appearance in structures 
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param stoptype which weighting to be used in the multi-objective optimization? Either 'additive' (default) or 'multiplicative'. 
#'
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress-1 value
#'         \item stress.m: default normalized stress
#'         \item stoploss: the weighted loss value
#'         \item struc: the structuredness indices
#'         \item parameters: the parameters used for fitting (tau)
#'         \item fit: the returned object of the fitting procedure
#'         \item{stopobj:} the stopobj object 
#' }
#' @keywords multivariate
#' @export
stop_cldak <- function(dis,theta=c(3*max(sd(dis)),nrow(dis)/4),type="ratio",weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,itmaxi=10000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(missing(stoptype)) stoptype <- "additive"
  if(length(theta)>4) stop("There are too many parameters in the theta argument.")
  if(length(theta)<4) theta <- rep(theta,length.out=2)
  #if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  wght <- weightmat
  diag(wght) <- 1
  lambda0 <- theta[1]
  k <- theta[2]
  fit <- smacofx::clda(delta=dis,lambda0=lambda0,k=k,Epochs=20,alpha0=0.5,weightmat=wght,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  ncall <- do.call(substitute,list(fit$call,list(lambda0=lambda0,k=k,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi)))
  fit$call <- ncall               
  fit$lambda0 <- theta[1]
  fit$k <- theta[2]
  fit$parameters <- fit$theta <- fit$pars <- c(lambda0=fit$lambda0,k=fit$k)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
  out 
}


#' STOPS version of CLDA with free epsilon.
#'
#' CLDA with free lambda0 and epsilon and 20 epochs. Should we add alpha0?
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of explicit parameters; first is lambda0 for the maximal neighbourhood and second is k for the number of neighbours for the geodesic distance. 
#' @param type MDS type. 
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param ndim number of dimensions of the target space
#' @param itmaxi number of iterations
#' @param ... additional arguments to be passed to the fitting procedure
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param structures a character vector listing the structure indices to use. They always are called "cfoo" with foo being the structure.
#' @param strucweight weight to be used for the structures; defaults to 1/number of structures
#' @param strucpars a list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appearance in structures 
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param stoptype which weighting to be used in the multi-objective optimization? Either 'additive' (default) or 'multiplicative'. 
#'
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress-1 value
#'         \item stress.m: default normalized stress
#'         \item stoploss: the weighted loss value
#'         \item struc: the structuredness indices
#'         \item parameters: the parameters used for fitting (tau)
#'         \item fit: the returned object of the fitting procedure
#'         \item{stopobj:} the stopobj object 
#' }
#' @keywords multivariate
#' @export
stop_cldae <- function(dis,theta=rep(3*max(sd(dis)),2),type="ratio",weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,itmaxi=10000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(missing(stoptype)) stoptype <- "additive"
  if(length(theta)>4) stop("There are too many parameters in the theta argument.")
  if(length(theta)<4) theta <- rep(theta,length.out=2)
  #if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  wght <- weightmat
  diag(wght) <- 1
  lambda0 <- theta[1]
  epsilon <- theta[2]
  fit <- smacofx::clda(delta=dis,lambda0=lambda0,epsilon=epsilon,Epochs=20,alpha0=0.5,weightmat=wght,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  ncall <- do.call(substitute,list(fit$call,list(lambda0=lambda0,epsilon=epsilon,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi)))
  fit$call <- ncall              
  fit$lambda0 <- lambda0
  fit$epsilon <- epsilon
  fit$parameters <- fit$theta <- fit$pars <- c(lambda0=fit$lambda0,epsilon=fit$epsilon)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
  out 
}


#' STOPS version of sparsified multidimensional distance analysis for fixed k and tau
#'
#' smdda with free parameters tau and k.
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of explicit parameters; first is tau for the neighbourhood, second is k. Defaults to 100, 10.
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
stop_smddak <- function(dis,theta=c(100,10),type="ratio",weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,itmaxi=10000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(missing(stoptype)) stoptype <- "additive"
  if(length(theta)>4) stop("There are too many parameters in the theta argument.")
  if(length(theta)<2) theta <- rep(theta,length.out=2)
  tau <- theta[1]
  k <- theta[2]
  #if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  wght <- as.matrix(weightmat)
  diag(wght) <- 1
  fit <- smacofx::smdda(delta=dis,tau=tau,k=k,type=type,weightmat=wght,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  ncall <- do.call(substitute,list(fit$call,list(tau=tau,k=k,type=type,init=init,weightmat=wght,ndim=ndim,verbose=verbose,itmax=itmaxi)))
  fit$call <- ncall                
  fit$tau <- tau
  fit$k <- k
  fit$parameters <- fit$theta <- fit$pars <- c(tau=fit$tau,k=fit$k)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
  out 
}


#' STOPS version of sparsified multidimensional distance analysis for fixed eps and tau
#'
#' smdda with free parameters tau and epsilon.
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of explicit parameters; first is tau for the neighboourhood, second is epsilon for isomapdist. Defaults to 100, 100.
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
stop_smddae <- function(dis,theta=c(100,100),type="ratio",weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,itmaxi=10000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(missing(stoptype)) stoptype <- "additive"
  if(length(theta)>4) stop("There are too many parameters in the theta argument.")
  if(length(theta)<2) theta <- rep(theta,length.out=2)
  tau <- theta[1]
  epsilon <- theta[2]
  #if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  wght <- as.matrix(weightmat)
  diag(wght) <- 1
  fit <- smacofx::smdda(delta=dis,tau=tau,epsilon=epsilon,type=type,weightmat=wght,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  ncall <- do.call(substitute,list(fit$call,list(tau=tau,epsilon=epsilon,type=type,init=init,weightmat=wght,ndim=ndim,verbose=verbose,itmax=itmaxi)))
  fit$call <- ncall                
  fit$tau <- tau
  fit$epsilon <- epsilon
  fit$parameters <- fit$theta <- fit$pars <- c(tau=fit$tau,epsilon=fit$epsilon)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
  out 
}

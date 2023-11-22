#' STOPS version of sammon with powers
#'
#' This is power stress with free kappa and lambda but rho is fixed to -1 and the weights are delta.
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; a vector of length two where the first element is kappa (for the fitted distances), the second lambda (for the observed proximities). If a scalar is given it is recycled for the free parameters.  Defaults to 1 1.
#' @param type MDS type. Defaults to "ratio".
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param ndim number of dimensions of the target space
#' @param itmaxi number of iterations
#' @param ... additional arguments to be passed to the fitting procedure
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param structures which structures to look for
#' @param strucweight weight to be used for the structures; defaults to 0.5
#' @param strucpars a list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appeacrance in structures 
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param stoptype which weighting to be used in the multi-objective optimization? Either 'additive' (default) or 'multiplicative'. 
#'
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item stoploss: the weighted loss value
#'         \item struc: the structuredness indices
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#'         \item{stopobj:} the stopobj object 
#' }
#'
#' @import cordillera
#' @keywords multivariate
#' @export
stop_powersammon <- function(dis,theta=c(1,1),type="ratio",weightmat=NULL,init=NULL,ndim=2,itmaxi=10000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)  
  if(missing(stoptype)) stoptype <- "additive"
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  #if(length(theta)<2) theta <- rep(theta,length.out=2)
  if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  weightmat <- as.matrix(weightmat)
  nu <- -1
  sammwght <-dis^(theta[2])
  diag(sammwght) <- 1
  combwght <- stats::as.dist(sammwght*weightmat) #powerStressMin will turn this into the correct data structure
  kappa <- theta[1]
  lambda<- theta[2]
  fit <- smacofx::powerStressMin(delta=dis,kappa=kappa,lambda=lambda,nu=nu,type=type,weightmat=combwght,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  ncall <- do.call(substitute,list(fit$call,list(kappa=kappa,lambda=lambda,nu=nu,type=type,init=init,weightmat=combwght,ndim=ndim,verbose=verbose,itmax=itmaxi)))
  fit$call <- ncall                
  fit$kappa <- kappa
  fit$lambda <- lambda
  fit$nu <- nu
  fit$parameters <- fit$theta <- fit$pars <- c(kappa=fit$kappa,lambda=fit$lambda)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
  out 
}


#' STOPS version of powerstress
#'
#' Power stress with free kappa and lambda and rho.
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances), the second lambda (for the observed proximities), the third nu (for the weights). If a scalar is given it is recycled.  Defaults to 1 1 1.
#' @param type MDS type. 
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param ndim number of dimensions of the target space
#' @param itmaxi maximum number of iterations
#' @param acc accuracy (defaults to 1e-8)
#' @param ... additional arguments to be passed to the fitting procedure
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param structures a character vector listing the structure indices to use. They always are called "cfoo" with foo being the structure.
#' @param strucweight weight to be used for the structures; defaults to 1/number of structures
#' @param strucpars a list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appearance in structures 
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param stoptype which weighting to be used in the multi-objective optimization? Either 'additive' (default) or 'multiplicative'. 
#' @param registry registry object with c-structuredness indices.
#' 
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress-1 value
#'         \item stress.m: default normalized stress
#'         \item stoploss: the weighted loss value
#'         \item struc: the structuredness indices
#'         \item parameters: the parameters used for fitting (kappa, lambda, nu)
#'         \item fit: the returned object of the fitting procedure
#'         \item{stopobj:} the stopobj object 
#' }
#' @keywords multivariate
#' @export
stop_powerstress <- function(dis,theta=c(1,1,1),type="ratio",weightmat=NULL,init=NULL,ndim=2,itmaxi=10000,acc=1e-8,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative"),registry=struc_reg) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  if(missing(stoptype)) stoptype <- "additive"
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)<3) theta <- rep(theta,length.out=3)
  weightmat <- as.matrix(weightmat)
  #diag(weightmat) <- 1
  kappa <- theta[1]
  lambda <- theta[2]
  nu <- theta[3]
  fit <- smacofx::powerStressMin(delta=dis,kappa=kappa,lambda=lambda,nu=nu,type=type,weightmat=weightmat,init=init,ndim=ndim,verbose=verbose-3,itmax=itmaxi,acc=acc,...)
  ncall <- do.call(substitute,list(fit$call,list(kappa=kappa,lambda=lambda,nu=nu,type=type,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,acc=acc)))
  fit$call <- ncall                 
  fit$kappa <- kappa
  fit$lambda <- lambda
  fit$nu <- nu
  fit$parameters <- fit$theta <- fit$pars <- c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype,registry=registry)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
  out 
}

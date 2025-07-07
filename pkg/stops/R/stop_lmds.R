#' STOPS version of lMDS
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; the first is k (for the neighbourhood), the second tau (for the penalty) . If a scalar is given it is recycled.  Defaults to 2 and 0.5.
#' @param type MDS type. Ignored.
#' @param weightmat (not used) 
#' @param init (optional) initial configuration
#' @param ndim number of dimensions of the target space
#' @param itmaxi maximum number of iterations
#' @param acc accuracy
#' @param ... additional arguments to be passed to the fitting procedure
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param structures which structures to look for
#' @param strucweight weight to be used for the structures; defaults to 0.5
#' @param strucpars a list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appeacrance in structure
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param stoptype which weighting to be used in the multi-objective optimization? Either 'additive' (default) or 'multiplicative'.
#' @param registry registry object with c-structuredness indices.
#'
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress-1
#'         \item stress.m: default normalized stress
#'         \item stoploss: the weighted loss value
#'         \item struc: the structuredness indices
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#'         \item{stopobj:} the stopobj object 
#' }
#' @keywords multivariate
#' @importFrom smacofx lmds 
#' @export
stop_lmds <- function(dis,theta=c(2,0.5),type="ratio",weightmat=NULL,init=NULL,ndim=2,itmaxi=10000,acc=1e-8,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative"),registry=struc_reg) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")|| is.data.frame(dis)) dis <- as.matrix(dis)
  if(missing(stoptype)) stoptype <- "additive"
  #we allow for three parameters in the theta argument
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)<2) theta <- rep(theta,length.out=2)
  #if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  #wght <- weightmat
                                        #diag(wght) <- 1
  k <- theta[1]
  tau <- theta[2]
  verbose <- verbose+2
  fit <- lmds(delta=dis,k=k,tau=tau,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,acc=acc,...)
  ncall <- do.call(substitute,list(fit$call,list(k=k,tau=tau,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,acc=acc)))
  fit$call <- ncall                
  fit$k <- k
  fit$tau <- tau
  fit$parameters <- fit$theta  <- fit$pars  <- c(k=fit$k,tau=fit$tau)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype,registry=registry)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
  out 
}

#' STOPS version of Box Cox Stress
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; the first is mu (for the fitted distances), the second lambda (for the  proximities), the third nu (for the weights). If a scalar is given it is recycled.  Defaults to 1 1 0.
#' @param type MDS type. Is ignored here.  
#' @param weightmat (not used) 
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
#'          \item{stopobj:} the stopobj object 
#' }
#' @keywords multivariate
#' @importFrom smacofx bcmds 
#' @export
stop_bcmds<- function(dis,theta=c(1,1,0),type="ratio",weightmat=NULL,init=NULL,ndim=2,itmaxi=5000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative"),registry=struc_reg) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist") || is.data.frame(dis)) dis <- as.matrix(dis)
  if(missing(stoptype)) stoptype <- "additive"
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)<3) theta <- rep(theta,length.out=3)
  #if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  #wght <- weightmat
                                        #diag(wght) <- 1
  verbose <- verbose+2
  mu <- theta[1]
  lambda <- theta[2]
  rho <- theta[3]
  fit <- smacofx::bcmds(delta=dis,mu=mu,lambda=lambda,rho=rho,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  ncall <- do.call(substitute,list(fit$call,list(mu=mu,lambda=lambda,rho=rho,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi)))
  fit$call <- ncall                
  fit$mu <- mu
  fit$lambda <- lambda
  fit$rho <- rho
  fit$parameters <- fit$theta  <- fit$pars <- c(mu=fit$mu,lambda=fit$lambda,rho=fit$rho)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype,registry=registry)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
  out 
}

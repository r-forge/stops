
#' STOPS version of restricted powerstress

# This is a power stress where kappa and lambda are free to vary but restricted to be equal, so the same exponent will be used for distances and dissimilarities. rho (for the weights) is also free.  
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; the first two arguments are for kappa and lambda and should be equal (for the fitted distances and observed proximities), the third nu (for the weights). Internally the kappa and lambda are equated. If a scalar is given it is recycled (so all elements of theta are equal); if a vector of length 2 is given, it gets expanded to c(theta[1],theta[1],theta[2]). Defaults to 1 1 1.
#' @param type MDS type. Defaults to "ratio".
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param ndim number of dimensions of the target space
#' @param itmaxi number of iterations. default is 10000.
#' @param ... additional arguments to be passed to the fitting procedure powerStressMin
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param structures a character vector listing the structure indices to use. They always are called "cfoo" with foo being the structure.
#' @param strucweight weight to be used for the structures; defaults to 1/number of structures
#' @param strucpars a list of list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appearance in structures vector. See examples.  
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
#'         \item parameters: the parameters used for fitting (kappa=lambda, nu)
#'         \item fit: the returned object of the fitting procedure
#'          \item{stopobj:} the stopobj object 
#' }
#' @importFrom smacofx powerStressMin
#' @keywords multivariate
stop_rpowerstress <- function(dis,theta=c(1,1,1),type="ratio",weightmat=NULL,init=NULL,ndim=2,itmaxi=10000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative"),registry=struc_reg) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(missing(stoptype)) stoptype <- "additive"
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)==1L) theta <- rep(theta,3)
  if(length(theta)==2L) theta <- c(rep(theta[1],2),theta[2])
  if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  wght <- as.matrix(weightmat)
  #diag(wght) <- 1
  expo <- theta[1]
  nu <- theta[3]
  fit <- smacofx::rpStressMin(delta=dis,expo=expo,nu=nu,weightmat=wght,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  ncall <- do.call(substitute,list(fit$call,list(expo=expo,nu=nu,type=type,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi)))
  fit$call <- ncall                 
  fit$kappa <- expo
  fit$lambda <- expo
  fit$nu <- nu
  fit$parameters <- fit$theta <- fit$pars <- c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype,registry=registry)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
  out 
}

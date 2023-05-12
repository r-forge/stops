
#' STOPS version of Sammon mapping
#'
#' Uses smacofx::sammon. The free parameter is lambda for power transformations of the observed proximities. 
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this must be  a scalar of the lambda transformation for the observed proximities. Defaults to 1.
#' @param type MDS type. Ignored here.
#' @param ndim number of dimensions of the target space
#' @param weightmat a matrix of nonnegative weights. Has no effect here.
#' @param init (optional) initial configuration
#' @param itmaxi number of iterations
#' @param ... additional arguments to be passed to the fitting procedure
#' @param structures which structuredness indices to be included in the loss
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
#' @param strucpars the parameters for the structuredness indices
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param stoptype How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
#' 
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress/1 *sqrt stress(
#'         \item{stress.m:} default normalized stress
#'         \item{stoploss:} the weighted loss value
#'         \item{indices:} the values of the structuredness indices
#'         \item{parameters:} the parameters used for fitting 
#'         \item{fit:} the returned object of the fitting procedure  smacofx::sammon 
#'          \item{stopobj:} the stopobj object
#' }
#'
#' @importFrom stats dist as.dist
#' @importFrom smacofx sammon
#' @import cordillera
#' @keywords multivariate
#'
#' 
#' @export
stop_sammon <- function(dis,theta=1,type="ratio",ndim=2,init=NULL,weightmat=NULL,itmaxi=1000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(missing(stoptype)) stoptype <- "additive"
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  type <- match.arg(type,'ratio')
  lambda <- theta[1]
  #if(length(theta)==1L) lambda <- theta
  #if(length(theta)==2L) lambda <- theta[2]
  #if(length(theta)==3L) lambda <- theta[2]
  #nu <- -1
  fit <- smacofx::sammon(dis^lambda,k=ndim,y=init,trace=isTRUE(verbose>1),niter=itmaxi,...)
  fit$lambda <- lambda
  #fit$kappa <- 1
  #fit$nu <- -1
  #N <- length(dis)
  #dis <- stats::as.dist(dis)
  #fitdis <- stats::dist(fit$points)
  #wghts <- 1/dis^lambda 
  #disl <- dis^lambda
  #dhat <-  disl/sqrt(sum(wghts*disl^2))*sqrt(N)
  #lb <- sum(wghts*fitdis*dhat)/sum(wghts*fitdis^2)   #Restrict config so we have a stress in [0,1] just as in smacof. Rest is unchanged. Maybe use this stress for optimization at some point?
  #fitdisnn <- lb*fitdis
  #fit$stress.r <- fit$stress
  #fit$stress <- sqrt(sum(wghts*(dhat-fitdisnn)^2)/N) #sum(wghts*dhat^2))
  #fit$stress.n <- fit$stress.r/sum(dis^lambda)
  #fit$stress.m <- fit$stress #or stress.r
                                        #fit$conf <- fit$points
  fit$delta <- stats::as.dist(dis)
  fit$parameters <- fit$theta <- c(lambda=lambda)
  #fit$pars <- c(lambda=fit$lambda)#c(kappa=fit$kappa,lambda=fit$lambda,rho=fit$nu)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype)
  list(stress=fit$stress, stress.m=fit$stress.m,stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters,  fit=fit,stopobj=stopobj) #target functions
}

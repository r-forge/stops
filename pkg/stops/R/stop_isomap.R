

#' STOPS version of isomap to optimize over integer k.
#'
#' Free parameter is k. 
#' 
#' @details Currently this version is a bit less flexible than the vegan one, as the only allowed parameter for isomap is the theta (k in isomap, no epsilon) and the shortest path is always estimated with argument "shortest". Also note that fragmentedOK is always set to TRUE which means that for theta that is too small only the largest conected group will be analyzed. If that's not wanted just set the theta higher.  
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the number of shortest dissimilarities retained for a point (nearest neighbours), the isomap parameter. Must be a numeric scalar. Defaults to 3.
#' @param type MDS type. Is "ratio". 
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param structures which structuredness indices to be included in the loss
#' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
#' @param strucpars the parameters for the structuredness indices
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param stoptype How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
#' @param itmaxi placeholder for compatibility in stops call; not used
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} Not really stress but 1-GOF[2] where GOF is the second element returned from smacofx::cmdscale (the sum of the first ndim eigenvalues divided by the sum of all absolute eigenvalues).
#'         \item{stress.m:} default normalized stress (sqrt explicitly normalized stress; really the stress this time)
#'         \item{stoploss:} the weighted loss value
#'         \item{indices:} the values of the structuredness indices
#'         \item{parameters:} the parameters used for fitting 
#'         \item{fit:} the returned object of the fitting procedure
#'          \item{stopobj:} the stopobj object
#' }
#'
#' @import cordillera
#' @importFrom stats dist as.dist
#' @importFrom vegan isomap isomapdist
#' @importFrom smacofx cmdscale 
#' @keywords multivariate
#' @export
stop_isomap1 <- function(dis,theta=3,type="ratio",weightmat=NULL,ndim=2,init=NULL,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative"),itmaxi=NULL) {
  theta <- as.numeric(theta)
  type <- "ratio"
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(missing(stoptype)) stoptype <- "additive"
  if(length(theta)<3) lambda <- theta[1] #I just call the k lambda
  #if(length(theta)==2L) lambda <- theta[1]
  #if(length(theta)==3L) lambda <- theta[1]
  disi <- vegan::isomapdist(dis,k=lambda,path="shortest",fragmentedOK=TRUE)
  if(length(disi)==0) stop("The distance matrix is of length 0 for the current k. Consider increasing the lower bound of the search region.")
  fit <- smacofx::cmdscale(disi,k=ndim,eig=TRUE,add=TRUE) 
  fit$k <- lambda
  #fit$kappa <- 1
  #fit$nu <- 1
  dis <- stats::as.dist(disi)
  fitdis <- stats::dist(fit$conf)
  stress.r <- sum((disi-fitdis)^2)
  stress.n <- stress.r/sum(disi^2)
  fit$stress.m <- stress.n
  #fit$stress <- sqrt(fit$stress.n)
  fit$parameters <- fit$theta <- fit$pars <- c(k=fit$k)
  #fit$conf <- fit$points
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype)
  list(stress=fit$stress,stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj) #target functions
}


#' STOPS version of isomap over real epsilon.
#'
#' Free parameter is eps. 
#' 
#' @details Currently this version is a bit less flexible than the vegan one, as the only allowed parameter for isomap is the theta (epsilon in isomap) and the shortest path is always estimated with argument "shortest". Also note that fragmentedOK is always set to TRUE which means that for theta that is too small only the largest conected group will be analyzed. If that's not wanted just set the theta higher.  
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the number of shortest dissimilarities retained for a point (neighbourhood region), the isomap parameter. Defaults to the 0.1 quantile of the empirical distribution of dis.
#' @param type MDS type. Is "ratio". 
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param structures which structuredness indices to be included in the loss
#' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
#' @param strucpars the parameters for the structuredness indices
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param stoptype How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
#' @param itmaxi placeholder for compatibility in stops call; not used
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} Not really stress but 1-GOF[2] where GOF is the second element returned from cmdscale (the sum of the first ndim absolute eigenvalues divided by the sum of all absolute eigenvalues).
#'         \item{stress.m:} default normalized stress (sqrt explicitly normalized stress; really the stress this time)
#'         \item{stoploss:} the weighted loss value
#'         \item{indices:} the values of the structuredness indices
#'         \item{parameters:} the parameters used for fitting 
#'         \item{fit:} the returned object of the fitting procedure
#'          \item{stopobj:} the stopobj object
#' }
#'
#' @import cordillera
#' @importFrom stats dist as.dist quantile
#' @importFrom vegan isomap isomapdist
#' @importFrom smacofx cmdscale
#' @keywords multivariate
#' @export
stop_isomap2 <- function(dis,theta=stats::quantile(dis,0.1),type="ratio",weightmat=NULL,ndim=2,init=NULL,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,stoptype=c("additive","multiplicative"),itmaxi=NULL) {
  theta <- as.numeric(theta)
  type <- "ratio"
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(missing(stoptype)) stoptype <- "additive"
  if(length(theta)<3) lambda <- theta[1]
  #if(length(theta)==2L) lambda <- theta[1]
  #if(length(theta)==3L) lambda <- theta[1]
  disi <- vegan::isomapdist(dis,epsilon=lambda,path="shortest",fragmentedOK=TRUE)
  if(length(disi)==0) stop("The distance matrix is of length 0 for the current epsilon. Consider increasing the lower bound of the search region.")
  fit <- smacofx::cmdscale(disi,k=ndim,eig=TRUE,add=TRUE) 
  fit$eps <- lambda
  #fit$kappa <- 1
  #fit$nu <- 1
  dis <- stats::as.dist(disi)
  fitdis <- stats::dist(fit$conf)
  stress.r <- sum((disi-fitdis)^2)
  stress.n <- stress.r/sum(disi^2)
  fit$stress.m <- stress.n
  #fit$stress <- sqrt(fit$stress.m)
  fit$parameters <- fit$theta <- fit$pars <- c(eps=fit$eps)
  #fit$conf <- fit$points
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),stoptype=stoptype)
  list(stress=fit$stress,stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj) #target functions
}

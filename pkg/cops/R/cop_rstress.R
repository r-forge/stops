#' PCOPS version of rstress
#'
#' Free parameter is kappa for the fitted distances.
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this must be a scalar of the kappa transformation for the fitted distances proximities. Defaults to 1. Note the kappa here differs from Jan's version where the parameter was called r and the relationship is r=kappa/2 or kappa=2r. 
#' @param ndim number of dimensions of the target space
#' @param itmaxi number of iterations. default is 10000.
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param cordweight weight to be used for the cordillera; defaults to 0.5
#' @param q the norm of the cordillera; defaults to 1
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to ndim+1
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the distances (min distance minus max distance). If NULL (default) the cordillera will be normed to each configuration's maximum distance, so an absolute value of goodness-of-clusteredness.
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param normed should the cordillera be normed; defaults to TRUE
#' @param scale should the configuration be scale adjusted
#' @param stresstype which stress to report? Defaults to explicitly normed stress
#' @param ... additional arguments to be passed to the fitting procedure
#' 
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item copstress: the weighted loss value
#'         \item OC: the Optics cordillera value
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#'         \item cordillera: the cordillera object
#' }
#' @keywords multivariate
#' @import cordillera
cop_rstress <- function(dis,theta=c(1,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,itmaxi=10000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale="sd",normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
  if(missing(stresstype)) stresstype <- "default"  
  if(length(theta)>1) stop("There are too many parameters in the theta argument.")
  kappa <- theta
  #if(length(theta)==3L) kappa <- theta[1] 
  fit <- powerStressMin(delta=dis,kappa=kappa,lambda=1,nu=1,weightmat=weightmat,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  if(stresstype=="default") fit$stress.m <- fit$stress.m
  if(stresstype=="stress1") fit$stress.m <- fit$stress.1
  if(stresstype=="rawstress") fit$stress.m <- fit$stress.r
  if(stresstype=="enormstress") fit$stress.m <- fit$stress.en
  if(stresstype=="enormstress1") fit$stress.m <- fit$stress.en1
  if(stresstype=="normstress") fit$stress.m <- fit$stress.n
  fit$kappa <- kappa
  #fit$lambda <- 1
  #fit$nu <- 1
  fit$parameters <- fit$theta <- c(kappa=fit$kappa)#c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
 # fit$deltaorig <- fit$delta^(1/fit$lambda)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit,copsobj=copobj)
  out
}

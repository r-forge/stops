#' PCOPS version of sstress
#'
#'Free parameter is lambda for the observed proximities. Fitted distances are transformed with power 2, weights have exponent of 1. Note that the lambda here works as a multiplicator of 2 (as sstress has f(delta^2)).
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this must be a scalar of the lambda transformation for the observed proximities. Defaults to 1. Note that the lambda here works as a multiplicator of 2 (as sstress has f(delta^2)).
#' @param type MDS type, defaults to "ratio".
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
#' @param ... additional arguments to be passed to the fitting procedure
#' 
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress-1 value
#'         \item stress.m: default normalized stress
#'         \item copstress: the weighted loss value
#'         \item OC: the Optics cordillera value
#'         \item parameters: the parameters used for fitting (lambda)
#'         \item fit: the returned object of the fitting procedure
#'         \item cordillera: the cordillera object
#' }
#' @keywords multivariate
#' @import cordillera
#' @importFrom smacofx powerStressMin
cop_sstress <- function(dis,theta=1,type="ratio",weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,itmaxi=10000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale="sd",normed=TRUE) {
  #if(isTRUE(all.equal(theta,c(1,1,1)))) theta <- 1
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  lambda <- theta[1]
  flambda <- lambda*2 #sstress is d^2 and delta^2 so f(delta^2)=delta^(2*1); lambda works in factors of 2  
  fit <- smacofx::powerStressMin(delta=dis,kappa=2,lambda=flambda,nu=1,type=type,weightmat=weightmat,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  fit$kappa <- 2
  fit$lambda <- lambda
  #fit$nu <- 1
  fit$parameters <- fit$theta <- fit$pars  <- c(lambda=fit$lambda,kappa=2)#c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit,copsobj=copobj)
  out
}

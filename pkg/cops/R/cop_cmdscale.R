#' PCOPS version of strain
#'
#' The free parameter is lambda for power transformations of the observed proximities.
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta  the theta vector of powers; this must be a scalar of the lambda transformation for the observed proximities.
#' @param ndim number of dimensions of the target space
#' @param itmaxi number of iterations. No effect here.
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param cordweight weight to be used for the cordillera; defaults to 0.5
#' @param q the norm of the corrdillera; defaults to 1
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
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item copstress: the weighted loss value
#'         \item OC: the Optics cordillera value
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#'         \item cordillera: the cordillera object
#' }
#' 
#' @importFrom stats dist as.dist
#' @import cordillera
#' @keywords multivariate
cop_cmdscale <- function(dis,theta=1,weightmat=NULL,ndim=2,init=NULL,itmaxi=1000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale="sd",normed=TRUE) {
  if(length(theta)>1) stop("There are too many parameters in the theta argument.")
  if(length(theta)==1L) lambda <- theta
  fit <- smacofx::cmdscale(dis^lambda,k=ndim,eig=TRUE,...) 
  fit$lambda <- lambda
 # fit$kappa <- 1
 # fit$nu <- 1
  dhat <- stats::as.dist(fit$dhat)
  fitdis <- stats::dist(fit$conf)
  stress.r <- sum((dhat-fitdis)^2)
  stress.n <- fit$stress.r/sum(dhat^2)
  fit$stress <- sqrt(fit$stress.n)
  fit$stress.m <- fit$stress.n
  fit$parameters <- fit$pars <- fit$theta <- c(lambda=lambda)#c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  list(stress=fit$GOF,stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit,copsobj=copobj) #target functions
}

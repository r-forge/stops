
#' PCOPS version of restricted powerstress.  
#
#'
#' This is a power stress where kappa and lambda are free to vary but restricted to be equal, so the same exponent will be used for distances and dissimilarities. nu (for the weights) is also free.  
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; the first two arguments are for kappa and lambda and should be equal (for the fitted distances and observed proximities), the third nu (for the weights). Internally the kappa and lambda are equated based on theta[1]. If a scalar is given it is recycled (so all elements of theta are equal); if a vector of length 2 is given, it gets expanded to c(theta[1],theta[1],theta[2]). Defaults to 1 1 1.
#' @param type MDS type. Defaults to "ratio".
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
#' @param scale  should the configuration be scale adjusted
#' @param ... additional arguments to be passed to the fitting procedure
#'
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress1 value (sqrt(stress.m))
#'         \item stress.m: default normalized stress
#'         \item copstress: the weighted loss value
#'         \item OC: the Optics cordillera value
#'         \item parameters: the  explicit  parameters used for fitting (kappa=lambda, nu)
#'         \item fit: the returned object of the fitting procedure
#'         \item cordillera: the cordillera object
#' }
#' @import cordillera
#' @importFrom smacofx powerStressMin
#' @keywords multivariate
cop_rpowerstress <- function(dis,theta=c(1,1,1),type="ratio",weightmat=NULL,init=NULL,ndim=2,itmaxi=10000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale="sd",normed=TRUE) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  #if(length(theta)==3L & theta[1]!=theta[2]) warning("The powers given for kappa and lambda do not agree. The first value in theta will be used for both kappa and lambda.")  
  if(length(theta)==1L) theta <- rep(theta,3)
  if(length(theta)==2L) theta <- c(rep(theta[1],2),theta[2])
  if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  diag(wght) <- 1
  expo <- theta[1]
  nu <- theta[3]
  fit <- smacofx::rpStressMin(delta=dis,expo=expo,nu=nu,weightmat=wght,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  ncall <- do.call(substitute,list(fit$call,list(expo=expo,nu=nu,type=type,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi)))
  fit$call <- ncall                 
  fit$kappa <- expo
  fit$lambda <- expo
  fit$nu <- nu
  fit$parameters <- fit$theta <- fit$pars <- c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
  # fit$deltaorig <- stats::as.dist(dis)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit, copsobj=copobj)
  out 
}

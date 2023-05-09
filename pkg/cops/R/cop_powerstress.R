
## #' COPS version of powerstress
## #'
## #' Power stress with free kappa and lambda and rho (the theta argument).
## #'
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the theta vector of powers; the first is kappa (for the fitted distances), the second lambda (for the observed proximities), the third nu (for the weights). If a scalar is given it is recycled.  Defaults to 1 1 1.
## #' @param ndim number of dimensions of the target space
## #' @param itmaxi number of iterations. default is 10000.
## #' @param weightmat (optional) a matrix of nonnegative weights
## #' @param init (optional) initial configuration
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param cordweight weight to be used for the cordillera; defaults to 0.5
## #' @param q the norm of the cordillera; defaults to 1
## #' @param minpts the minimum points to make up a cluster in OPTICS; defaults to ndim+1
## #' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
## #' @param rang range of the distances (min distance minus max distance). If NULL (default) the cordillera will be normed to each configuration's maximum distance, so an absolute value of goodness-of-clusteredness.
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param normed should the cordillera be normed; defaults to TRUE
## #' @param scale  should the configuration be scale adjusted
## #' @param stresstype which stress to report? Defaults to explicitly normed stress
## #' @param ... additional arguments to be passed to the fitting procedure
## #'
## #' @return A list with the components
## #' \itemize{
## #'         \item stress: the stress
## #'         \item stress.m: default normalized stress
## #'         \item copstress: the weighted loss value
## #'         \item OC: the Optics cordillera value
## #'         \item parameters: the parameters used for fitting (kappa, lambda)
## #'         \item fit: the returned object of the fitting procedure
## #'         \item cordillera: the cordillera object
## #' }
## #' @import cordillera
## #' @keywords multivariate
## cop_powerstress <- function(dis,theta=c(1,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,itmaxi=10000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale="sd",normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
##   if(missing(stresstype)) stresstype <- "default"
##   if(length(theta)>3) stop("There are too many parameters in the theta argument.")
##   if(length(theta)<3) theta <- rep(theta,length.out=3)
##   wght <- weightmat
##   diag(wght) <- 1
##   fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=theta[3],weightmat=wght,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
##   #fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=theta[3],weightmat=wght,ndim=ndim,verbose=verbose,itmax=itmaxi)
##   if(stresstype=="default") fit$stress.m <- fit$stress.m
##   if(stresstype=="stress1") fit$stress.m <- fit$stress.1
##   if(stresstype=="rawstress") fit$stress.m <- fit$stress.r
##   if(stresstype=="normstress") fit$stress.m <- fit$stress.n
##   if(stresstype=="enormstress") fit$stress.m <- fit$stress.en
##   if(stresstype=="enormstress1") fit$stress.m <- fit$stress.en1
##   fit$kappa <- theta[1]
##   fit$lambda <- theta[2]
##   fit$nu <- theta[3]
##   fit$parameters <- fit$theta <- c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
##   copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
##   out <- list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit, copsobj=copobj)
##   out 
## }

#' COPS version of powerstress
#'
#' Power stress with free kappa and lambda and rho (the theta argument) and ratio and interval optimal scaling.
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances), the second lambda (for the observed proximities), the third nu (for the weights). If a scalar is given it is recycled.  Defaults to 1 1 1.
#' @param type MDS type. Default is ratio.
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
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item copstress: the weighted loss value
#'         \item OC: the Optics cordillera value
#'         \item parameters: the parameters used for fitting (kappa, lambda, nu)
#'         \item fit: the returned object of the fitting procedure plus a slot for the original data $deltaorig
#'         \item cordillera: the cordillera object
#' }
#' @import cordillera
#' @importFrom smacofx powerStressMin
#' @keywords multivariate
cop_powerstress <- function(dis,theta=c(1,1,1),type="ratio",weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,itmaxi=10000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale="sd",normed=TRUE) {
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)<3) theta <- rep(theta,length.out=3)
  wght <- weightmat
  diag(wght) <- 1
  fit <- smacofx::powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=theta[3],type=type,weightmat=wght,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  fit$kappa <- theta[1]
  fit$lambda <- theta[2]
  fit$nu <- theta[3]
  fit$parameters <- fit$theta <- c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
  fit$deltaorig <- stats::as.dist(dis)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit, copsobj=copobj)
  out 
}


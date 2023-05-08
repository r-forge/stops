#' PCOPS versions of smacofSphere models
#'
#' The free parameter is lambda for power transformations the observed proximities. The fitted distances power is internally fixed to 1 and the power for the weights is 1. 
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this must be a scalar of the lambda transformation for the observed proximities. Defaults to 1.
#' @param type MDS type
#' @param ndim number of dimensions of the target space
#' @param itmaxi number of iterations. default is 1000.
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
#' @param stresstype which stress to report. Only takes smacofs default stress currrently.
#' @param ... additional arguments to be passed to the fitting procedure
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress
#'         \item{stress.m:} default normalized stress
#'         \item{copstress:} the weighted loss value
#'         \item{OC:} the Optics cordillera value
#'         \item{parameters:} the parameters used for fitting (kappa, lambda)
#'         \item{fit:} the returned object of the fitting procedure
#'         \item{cordillera:} the cordillera object
#' }
#'
#'@import smacof
#'@import cordillera 
#'@importFrom stats dist as.dist
#'@keywords multivariate
cop_smacofSphere <- function(dis,theta=1,type="ratio",ndim=2,weightmat=NULL,init=NULL,itmaxi=1000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,normed=TRUE,scale="sd",stresstype="default") {
                                        #TODO Unfolding
  if(is.null(init)) init <- "torgerson"
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
  #kappa first argument, lambda=second
  if(length(theta)>1) stop("There are too many parameters in the theta argument.")
  lambda <- theta
 # if(length(theta)==2L) lambda <- theta[2]
  #addargs <- list(...)
  #addargs
  fit <- smacof::smacofSphere(dis^lambda,type=type,ndim=ndim,weightmat=weightmat,init=init,verbose=isTRUE(verbose==2),itmax=itmaxi,...) #optimize with smacof
  #fit$kappa <- 1
  fit$lambda <- lambda
  #nu <- 
  fit$stress.1 <- fit$stress
 # fit$stress <- (fit$stress^2)*sum(fit$obsdiss^2) check if this is like below
 # fitdis <- 2*sqrt(sqdist(fit$conf))
  fitdis <- as.matrix(fit$confdist)
  delts <- as.matrix(fit$delta)[-1,-1]
  fit$stress.r <- sum(weightmat*(delts-fitdis)^2)
  fit$stress.m <- fit$stress^2#fit$stress.r/sum(weightmat*delts^2)
  fit$parameters <- fit$theta <- c(lambda=fit$lambda) #c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
  fit$deltaorig <- fit$delta^(1/fit$lambda)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit,copsobj=copobj) #target functions
  out
}

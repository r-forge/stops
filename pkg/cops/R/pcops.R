#' PCOPS versions of smacofSym models
#'
#' The free parameter is lambda for power transformations the observed proximities. The fitted distances power is internally fixed to 1 and the power for the weights is 1. 
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector; must be a scalar for the lambda (proximity) transformation. Defaults to 1.
#' @param ndim number of dimensions of the target space
#' @param itmaxi number of iterations. default is 1000
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
#'         \item{fit:} the returned object of the fitting procedure (which has all smacofB elements and some more)
#'         \item{cordillera:} the cordillera object
#' }
#'
#'@importFrom stats dist as.dist
#'@import cordillera 
#'@import smacof
#' 
#'@keywords multivariate
cop_smacofSym <- function(dis,theta=1,ndim=2,weightmat=NULL,init=NULL,itmaxi=1000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,normed=TRUE,scale="sd",stresstype="default") {
                                        #TODO Unfolding
  if(is.null(init)) init <- "torgerson"
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
  #kappa first argument, lambda=second
  if(length(theta)>1) stop("There are too many parameters in the theta argument.")
  lambda <- theta
 # if(length(theta)==3L) lambda <- theta[2]
 # addargs <- list(...)
 # addargs
  fit <- smacof::smacofSym(dis^lambda,ndim=ndim,weightmat=weightmat,init=init,verbose=isTRUE(verbose==2),itmax=itmaxi,...) #optimize with smacof
  #fit$kappa <- 1
  fit$lambda <- lambda
  #fit$nu <- 1
  fit$stress.1 <- fit$stress
 # fit$stress <- (fit$stress^2)*sum(fit$obsdiss^2) check if this is like below
 # fitdis <- 2*sqrt(sqdist(fit$conf))
  fitdis <- as.matrix(fit$confdist)
  delts <- as.matrix(fit$delta) #That was my choice to not use the normalized deltas but try it on the original; that is scale and unit free as Buja said
 # if(inherits(fit,"smacofSP")) delts <- as.matrix(fit$delta)[-1,-1]
  fit$stress.r <- sum(weightmat*(delts-fitdis)^2)
  fit$stress.m <- fit$stress^2#fit$stress.r/sum(weightmat*delts^2)
  fit$parameters <- fit$theta <- c(lambda=fit$lambda)#c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
  fit$deltaorig <- fit$delta^(1/fit$lambda)  
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit,copsobj=copobj) #target functions
  out
}

#' PCOPS versions of elastic scaling models (via smacofSym)
#'
#' The free parameter is lambda for power transformations the observed proximities. The fitted distances power is internally fixed to 1 and the power for the weights=delta is -2. Allows for a weight matrix because of smacof.
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this must be a scalar of the lambda transformation for the observed proximities. Defaults to 1.
#' @param ndim number of dimensions of the target space
#' @param itmaxi number of iterations. default is 1000.
#' @param weightmat (optional) a matrix of nonnegative weights (NOT the elscal weights)
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
#'
#'@importFrom stats dist as.dist
#'@import smacof
#'@import cordillera 
#' 
#'@keywords multivariate
cop_elastic <- function(dis,theta=1,ndim=2,weightmat=1,init=NULL,itmaxi=1000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,normed=TRUE,scale="sd",stresstype="default") {
                                        #TODO Unfolding
  if(is.null(init)) init <- "torgerson"
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1]) 
  #kappa first argument, lambda=second
  if(length(theta)>1) stop("There are too many parameters in the theta argument.")
  lambda <- theta
  #if(length(theta)==3L) lambda <- theta[2]
  nu <- -2
  #addargs <- list(...)
  #addargs
  elscalw <- dis^(nu*lambda) #the weighting in elastic scaling
  diag(elscalw) <- 1
  combwght <- weightmat*elscalw #combine the user weights and the elastic scaling weights
  fit <- smacof::smacofSym(dis^lambda,ndim=ndim,weightmat=combwght,init=init,verbose=isTRUE(verbose==2),itmax=itmaxi,...) #optimize with smacof
  #fit$kappa <- 1
  fit$lambda <- lambda
  #fit$nu <- nu
  fit$stress.1 <- fit$stress
 # fit$stress <- (fit$stress^2)*sum(fit$obsdiss^2) check if this is like below
 # fitdis <- 2*sqrt(sqdist(fit$conf))
  fitdis <- as.matrix(fit$confdist)
  delts <- as.matrix(fit$delta) 
  fit$stress.r <- sum(combwght*((delts-fitdis)^2))
  fit$stress.m <- fit$stress^2#fit$stress.r/sum(combwght*delts^2)
  fit$parameters <- fit$theta <- c(lambda=fit$lambda)#c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
  fit$deltaorig <- fit$delta^(1/fit$lambda)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit,copsobj=copobj) #target functions
  out
}


#' PCOPS versions of smacofSphere models
#'
#' The free parameter is lambda for power transformations the observed proximities. The fitted distances power is internally fixed to 1 and the power for the weights is 1. 
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this must be a scalar of the lambda transformation for the observed proximities. Defaults to 1.
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
cop_smacofSphere <- function(dis,theta=1,ndim=2,weightmat=NULL,init=NULL,itmaxi=1000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,normed=TRUE,scale="sd",stresstype="default") {
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
  fit <- smacof::smacofSphere(dis^lambda,ndim=ndim,weightmat=weightmat,init=init,verbose=isTRUE(verbose==2),itmax=itmaxi,...) #optimize with smacof
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


#' PCOPS version of Sammon mapping
#'
#' Uses MASS::sammon. The free parameter is lambda for power transformations of the observed proximities. The fitted distances power is internally fixed to 1 and the power for the weights=delta is -1. 
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this must be  a scalar of the lambda transformation for the observed proximities. Defaults to 1.
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
#'
#' 
cop_sammon <- function(dis,theta=1,ndim=2,init=NULL,weightmat=NULL,itmaxi=100,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale="sd",normed=TRUE,stresstype="default") {
  if(length(theta)>1) stop("There are too many parameters in the theta argument.")
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  lambda <- theta
  #if(length(theta)==2L) lambda <- theta[2]
  #if(length(theta)==3L) lambda <- theta[2]
  #nu <- -1
 # if(is.null(init)) init <- cops::cmdscale(dis^lambda,k=ndim)$points
  fit <- cops::sammon(dis^lambda,k=ndim,y=init,trace=isTRUE(verbose>1),niter=itmaxi,...)
  fit$lambda <- lambda
  #fit$kappa <- 1
  #fit$nu <- -1
  dis <- stats::as.dist(dis)
  fitdis <- stats::dist(fit$points)
  fit$stress.r <- sum(((dis^lambda-fitdis)^2)/dis)
  fit$stress.n <- fit$stress.r/sum(dis)
  fit$stress.m <- fit$stress^2#TODO: Passt das?
  fit$conf <- fit$points
  fit$parameters <- fit$theta <- c(lambda=fit$lambda)#c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters,  fit=fit,copsobj=copobj) #target functions
}

#' Another COPS versions of Sammon mapping models (via smacofSym)
#'
#' Uses Smacof, so it can deal with a weight matrix too.  The free parameter is lambda for power transformations of the observed proximities. The fitted distances power is internally fixed to 1 and the power for the weights=delta is -1. 
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta theta the theta vector of powers; this must be a scalar of the lambda transformation for the observed proximities. Defaults to 1. 
#' @param ndim number of dimensions of the target space
#' @param itmaxi number of iterations. default is 1000.
#' @param weightmat (optional) a matrix of nonnegative weights (NOT the sammon weights)
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
#'
#' @importFrom stats dist as.dist
#' @import smacof
#' @import cordillera
#'@keywords multivariate
cop_sammon2 <- function(dis,theta=1,ndim=2,weightmat=NULL,init=NULL,itmaxi=1000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,normed=TRUE,scale="sd",stresstype="default") {
    if(is.null(init)) init <- "torgerson" 
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1]) 
  #kappa first argument, lambda=second
  if(length(theta)>1) stop("There are too many parameters in the theta argument.")
  lambda <- theta
  #if(length(theta)==3L) lambda <- theta[2]
  nu <- -1
  #addargs <- list(...)
  elscalw <- dis^(nu*lambda) #the weighting in elastic scaling
  diag(elscalw) <- 1
  combwght <- weightmat*elscalw #combine the user weights and the elastic scaling weights
  fit <- smacof::smacofSym(dis^lambda,ndim=ndim,weightmat=combwght,init=init,verbose=isTRUE(verbose==2),itmax=itmaxi,...) #optimize with smacof
  #fit$kappa <- 1
  fit$lambda <- lambda
  #fit$nu <- nu
  fit$stress.1 <- fit$stress
 # fit$stress <- (fit$stress^2)*sum(fit$obsdiss^2) check if this is like below
 # fitdis <- 2*sqrt(sqdist(fit$conf))
  fitdis <- as.matrix(fit$confdist)
  delts <- as.matrix(fit$delta) 
  fit$stress.r <- sum(combwght*((delts-fitdis)^2))
  fit$stress.m <- fit$stress^2#fit$stress.r/sum(combwght*delts^2)
  fit$parameters <- fit$theta <- c(lambda=fit$lambda)#c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
  fit$deltaorig <- fit$delta^(1/fit$lambda)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit,copsobj=copobj) #target functions
  out
}



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
#' @param stresstype which stress to report. Only takes cmdscales default stress currrently.
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
cop_cmdscale <- function(dis,theta=c(1,1,1),weightmat=NULL,ndim=2,init=NULL,itmaxi=1000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale="sd",normed=TRUE,stresstype="default") {
  if(length(theta)>1) stop("There are too many parameters in the theta argument.")
  if(length(theta)==1L) lambda <- theta
  #if(length(theta)==2L) lambda <- theta[2]
  #if(length(theta)==3L) lambda <- theta[2]
  fit <- cops::cmdscale(dis^lambda,k=ndim,eig=TRUE,...) 
  fit$lambda <- lambda
 # fit$kappa <- 1
 # fit$nu <- 1
  dis <- stats::as.dist(dis)
  fitdis <- stats::dist(fit$points)
  fit$stress.r <- sum((dis^lambda-fitdis)^2)
  fit$stress.n <- fit$stress.r/sum(dis^(2*lambda))
  fit$stress <- sqrt(fit$stress.n)
  fit$stress.m <- fit$stress.n
  fit$conf <- fit$points
  fit$parameters <- fit$theta <- c(lambda=fit$lambda)#c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  list(stress=fit$GOF,stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit,copsobj=copobj) #target functions
}

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


#' PCOPS version of sstress
#'
#'Free parameter is lambda for the observed proximities. Fitted distances are transformed with power 2, weights have exponent of 1. Note that the lambda here works as a multiplicator of 2 (as sstress has f(delta^2)).
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this must be a scalar of the lambda transformation for the observed proximities. Defaults to 1. Note that the lambda here works as a multiplicator of 2 (as sstress has f(delta^2)). 
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
cop_sstress <- function(dis,theta=c(2,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,itmaxi=10000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale="sd",normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
  if(missing(stresstype)) stresstype <- "default"  
  if(length(theta)>1) stop("There are too many parameters in the theta argument.")
  lambda <- theta
  #if(length(theta)==3L) lambda <- theta[2]
  flambda <- lambda*2 #sstress is d^2 and delta^2 so f(delta^2)=delta^(2*1); lambda works in factors of 2  
  fit <- powerStressMin(delta=dis,kappa=2,lambda=flambda,nu=1,weightmat=weightmat,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  if(stresstype=="default") fit$stress.m <- fit$stress.m
  if(stresstype=="stress1") fit$stress.m <- fit$stress.1
  if(stresstype=="rawstress") fit$stress.m <- fit$stress.r
  if(stresstype=="normstress") fit$stress.m <- fit$stress.n
  if(stresstype=="enormstress") fit$stress.m <- fit$stress.en
  if(stresstype=="enormstress1") fit$stress.m <- fit$stress.en1
  #fit$kappa <- 2
  fit$lambda <- lambda
  #fit$nu <- 1
  fit$parameters <- fit$theta <- c(lambda=fit$lambda)#c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit,copsobj=copobj)
  out
}


#' PCOPS version of powermds
#'
#' This is power stress with free kappa and lambda but nu is fixed to 1, so no weight transformation.
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; a vector of length 2 where the first element is kappa (for the fitted distances), the second lambda (for the observed proximities). If a scalar is given it is recycled.  Defaults to 1.
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
#' @param stresstype which stress to report? Defaults to whatever whim is my default (currently explicitly normed stress)
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
cop_powermds <- function(dis,theta=c(1,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,itmaxi=itmaxi,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale="sd",normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
  if(missing(stresstype)) stresstype <- "default"
  if(length(theta)>2) stop("There are too many parameters in the theta argument.")
  if(length(theta)<2) theta <- rep(theta,length.out=2)
  fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=1,weightmat=weightmat,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  if(stresstype=="default") fit$stress.m <- fit$stress.m
  if(stresstype=="stress1") fit$stress.m <- fit$stress.1
  if(stresstype=="rawstress") fit$stress.m <- fit$stress.r
  if(stresstype=="normstress") fit$stress.m <- fit$stress.n
  if(stresstype=="enormstress") fit$stress.m <- fit$stress.en
  if(stresstype=="enormstress1") fit$stress.m <- fit$stress.en1
  fit$kappa <- theta[1]
  fit$lambda <- theta[2]
  #fit$nu <- 1
   fit$parameters <- fit$theta <- c(kappa=fit$kappa,lambda=fit$lambda)#c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit, copsobj=copobj)
  out 
}

#' PCOPS version of sammon with powers
#'
#' This is power stress with free kappa and lambda but nu is fixed to -1 and the weights are delta.
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; a vector of length two where the first element is kappa (for the fitted distances), the second lambda (for the observed proximities). If a scalar is given it is recycled for the free parameters.  Defaults to 1 1..
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
#' @import cordillera
#' @keywords multivariate
cop_powersammon <- function(dis,theta=c(1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,itmaxi=10000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale="sd",normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
  if(missing(stresstype)) stresstype <- "default"
  if(length(theta)>2) stop("There are too many parameters in the theta argument.")
  if(length(theta)==1L) theta <- rep(theta,2)
  nu <- -1
  sammwght <-dis^(theta[2])
  diag(sammwght) <- 1
  combwght <- sammwght*weightmat
  fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=nu,weightmat=combwght,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  if(stresstype=="default") fit$stress.m <- fit$stress.m
  if(stresstype=="stress1") fit$stress.m <- fit$stress.1
  if(stresstype=="rawstress") fit$stress.m <- fit$stress.r
  if(stresstype=="normstress") fit$stress.m <- fit$stress.n
  if(stresstype=="enormstress") fit$stress.m <- fit$stress.en
  if(stresstype=="enormstress1") fit$stress.m <- fit$stress.en1
  fit$kappa <- theta[1]
  fit$lambda <- theta[2]
  #fit$nu <- nu
  fit$parameters <- fit$theta <- c(kappa=fit$kappa,lambda=fit$lambda)# c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit, copsobj=copobj)
  out 
}

#' PCOPS version of elastic scaling with powers
#'
# This is power stress with free kappa and lambda but nu is fixed to -2 and the weights are delta.
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers;  a vector of length two where the first element is kappa (for the fitted distances), the second lambda (for the observed proximities). If a scalar for the free parameters is given it is recycled.  Defaults to 1 1.
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
#' @import cordillera
#' @keywords multivariate
cop_powerelastic <- function(dis,theta=c(1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,itmaxi=10000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale="sd",normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
  if(missing(stresstype)) stresstype <- "default"
  if(length(theta)>2) stop("There are too many parameters in the theta argument.")
  if(length(theta)==1L) theta <- rep(theta,2)
  nu <- -2
  elawght <- dis^(theta[2])
  diag(elawght) <- 1
  combwght <- elawght*weightmat
  fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=nu,weightmat=combwght,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  if(stresstype=="default") fit$stress.m <- fit$stress.m
  if(stresstype=="stress1") fit$stress.m <- fit$stress.1
  if(stresstype=="rawstress") fit$stress.m <- fit$stress.r
  if(stresstype=="normstress") fit$stress.m <- fit$stress.n
  if(stresstype=="enormstress") fit$stress.m <- fit$stress.en
  if(stresstype=="enormstress1") fit$stress.m <- fit$stress.en1
  fit$kappa <- theta[1]
  fit$lambda <- theta[2]
  #fit$nu <- nu
  fit$parameters <- fit$theta <-  c(kappa=fit$kappa,lambda=fit$lambda) # c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit, copsobj=copobj)
  out 
}

#' COPS version of powerstress
#'
#' Power stress with free kappa and lambda and rho (the theta argument).
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances), the second lambda (for the observed proximities), the third nu (for the weights). If a scalar is given it is recycled.  Defaults to 1 1 1.
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
#' @import cordillera
#' @keywords multivariate
cop_powerstress <- function(dis,theta=c(1,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,itmaxi=10000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale="sd",normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
  if(missing(stresstype)) stresstype <- "default"
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)<3) theta <- rep(theta,length.out=3)
  wght <- weightmat
  diag(wght) <- 1
  fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=theta[3],weightmat=wght,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  #fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=theta[3],weightmat=wght,ndim=ndim,verbose=verbose,itmax=itmaxi)
  if(stresstype=="default") fit$stress.m <- fit$stress.m
  if(stresstype=="stress1") fit$stress.m <- fit$stress.1
  if(stresstype=="rawstress") fit$stress.m <- fit$stress.r
  if(stresstype=="normstress") fit$stress.m <- fit$stress.n
  if(stresstype=="enormstress") fit$stress.m <- fit$stress.en
  if(stresstype=="enormstress1") fit$stress.m <- fit$stress.en1
  fit$kappa <- theta[1]
  fit$lambda <- theta[2]
  fit$nu <- theta[3]
  fit$parameters <- fit$theta <- c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit, copsobj=copobj)
  out 
}


#' PCOPS version of restricted powerstress.  
#
#'
#' This is a power stress where kappa and lambda are free to vary but restricted to be equal, so the same exponent will be used for distances and dissimilarities. nu (for the weights) is also free.  
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; the first two arguments are for kappa and lambda and should be equal (for the fitted distances and observed proximities), the third nu (for the weights). Internally the kappa and lambda are equated. If a scalar is given it is recycled (so all elements of theta are equal); if a vector of length 2 is given, it gets expanded to c(theta[1],theta[1],theta[2]). Defaults to 1 1 1.
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
#' @param stresstype which stress to report? Defaults to explicitly normed stress
#' @param ... additional arguments to be passed to the fitting procedure
#'
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress1 value (sqrt(stress.m))
#'         \item stress.m: default normalized stress
#'         \item copstress: the weighted loss value
#'         \item OC: the Optics cordillera value
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#'         \item cordillera: the cordillera object
#' }
#' @import cordillera
#' @keywords multivariate
cop_rpowerstress <- function(dis,theta=c(1,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,itmaxi=10000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale="sd",normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
  if(missing(stresstype)) stresstype <- "default"
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  #if(length(theta)==3L & theta[1]!=theta[2]) warning("The powers given for kappa and lambda do not agree. The first value in theta will be used for both kappa and lambda.")  
  if(length(theta)==1L) theta <- rep(theta,3)
  if(length(theta)==2L) theta <- c(rep(theta[1],2),theta[2])
  wght <- weightmat
  diag(wght) <- 1
  fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[1],nu=theta[3],weightmat=wght,init=init,ndim=ndim,verbose=verbose,itmax=itmaxi,...)
  #fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=theta[3],weightmat=wght,ndim=ndim,verbose=verbose,itmax=itmaxi)
  if(stresstype=="default") fit$stress.m <- fit$stress.m
  if(stresstype=="stress1") fit$stress.m <- fit$stress.1
  if(stresstype=="rawstress") fit$stress.m <- fit$stress.r
  if(stresstype=="normstress") fit$stress.m <- fit$stress.n
  if(stresstype=="enormstress") fit$stress.m <- fit$stress.en
  if(stresstype=="enormstress1") fit$stress.m <- fit$stress.en1
  fit$kappa <- theta[1]
  fit$lambda <- theta[1]
  fit$nu <- theta[3]
  fit$parameters <- fit$theta <- c(kappa=fit$kappa,lambda=fit$lambda,nu=fit$nu)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit, copsobj=copobj)
  out 
}


#' PCOPS version of approximated power stress model.
#'
#' This uses an approximation to power stress that can make use of smacof as workhorse. Free parameters are tau and upsilon.
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of parameters to optimize over. Must be of length two, with the first the tau argument and the second the upsilon argument. It can also be a scalar of the tau and upsilon transformation for the observed proximities and gets recycled for both ups and tau (so they are equal). Defaults to 1 1.
#' @param ndim number of dimensions of the target space
#' @param itmaxi number of iterations. default is 1000.
#' @param weightmat (optional) a binary matrix of nonnegative weights
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
#'         \item{fit:} the returned object of the fitting procedure (which has all smacofB elements and some more
#'         \item{cordillera:} the cordillera object
#' }
#'
#'@importFrom stats dist as.dist
#'@import cordillera 
#'@import smacof
#' 
#'@keywords multivariate
cop_apstress <- function(dis,theta=c(1,1,1),ndim=2,weightmat=NULL,init=NULL,itmaxi=1000,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,normed=TRUE,scale="sd",stresstype="default") {
  #TODO Unfolding  
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
  if(length(setdiff(unique(unlist(as.vector(weightmat))),c(0,1)))>0) stop("For approximated power stress, only binary weight matrices are allowed.")   
  #kappa first argument, lambda=second
  if(length(theta)>2) stop("There are too many parameters in the theta argument.")
  if(length(theta)==1L) theta <- rep(theta,2)
# tau <- theta
  #if(length(theta)==3L)  {
        tau <- theta[1]
        ups <- theta[2]
   #     }
 # addargs <- list(...)                                        # addargs
  combwght <- weightmat*(dis^ups)
  fit <- smacof::smacofSym(dis^tau,ndim=ndim,weightmat=combwght,init=init,verbose=isTRUE(verbose==2),itmax=itmaxi,...) #optimize with smacof
  #fit$kappa <- 1
  fit$tau <- tau
  fit$upsilon <- ups
  fit$stress.1 <- fit$stress
 # fit$stress <- (fit$stress^2)*sum(fit$obsdiss^2) check if this is like below
 # fitdis <- 2*sqrt(sqdist(fit$conf))
  fitdis <- as.matrix(fit$confdist)
  delts <- as.matrix(fit$delta) #That was my choice to not use the normalized deltas but try it on the original; that is scale and unit free as Buja said
 # if(inherits(fit,"smacofSP")) delts <- as.matrix(fit$delta)[-1,-1]
  fit$stress.r <- sum(combwght*(delts-fitdis)^2)
  fit$stress.m <- fit$stress^2#fit$stress.r/sum(combwght*delts^2)
  fit$parameters <- fit$theta <- c(tau=fit$tau,upsilon=fit$upsilon)#c(kappa=fit$kappa,tau=fit$tau,upsilon=fit$upsilon)
  fit$deltaorig <- fit$delta^(1/fit$tau)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress^2, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit,copsobj=copobj) #target functions
  out
}


#' Calculates copstress for given MDS object 
#'
#' @param obj MDS object (supported are sammon, cmdscale, smacof, rstress, powermds)
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param cordweight weight to be used for the cordillera; defaults to 0.5
#' @param q the norm of the cordillera; defaults to 1
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to 2
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the distances (min distance minus max distance). If NULL (default) the cordillera will be normed to each configuration's maximum distance, so an absolute value of goodness-of-clusteredness. 
#' @param verbose numeric value hat prints information on the fitting process; >2 is very verbose (copstress level), >3 is extremely (up to MDS optimization level)
#' @param normed should the cordillera be normed; defaults to TRUE
#' @param scale should the configuration be scale adjusted.
#' @param init a reference configuration when doing procrustes adjustment
#' @param ... additional arguments to be passed to the cordillera function
#'
#' @return A list with the components
#' \itemize{
#'        \item copstress: the weighted loss value
#'        \item OC: the Optics cordillera value
#'        \item parameters: the parameters used for fitting (kappa, lambda)
#'        \item cordillera: the cordillera object
#' }
#' @keywords multivariate
#' @import cordillera
copstress <- function(obj,stressweight=1,cordweight=5,q=1,minpts=2,epsilon=10,rang=NULL,verbose=0,normed=TRUE,scale=c("std","sd","proc","none"),init,...)
{
        if(missing(scale)) scale <- "sd"
        stressi <- obj$stress.m
        #kappa <- obj$kappa
        #lambda <- obj$lambda
        #nu <- obj$nu
        confs <- obj$conf
        confs <- scale_adjust(confs,init,scale=scale)
        corrd <- cordillera::cordillera(confs,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=FALSE,...)
        struc <- corrd$raw
        maxstruc <- corrd$normi
        if(normed) {
                   struc <- corrd$normed
                   maxstruc <- 1
                   }
        ic <- stressweight*stressi - cordweight*struc
        if(verbose>0) cat("copstress =",ic,"mdsloss =",stressi,"OC =",struc,"theta =",obj$parameters,"\n")
        out <- list(copstress=ic,OC=struc,parameters=obj$parameters,cordillera=corrd)
        out
     }

#' Profile COPS Function (aka COPS Variant 2)
#'
#' Metaparameter selection for MDS models baseed on the Profile COPS approach (COPS Variant 2). It uses copstress for hyperparameter selection. It is a special case of a STOPS model.  
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param loss which loss function to be used for fitting, defaults to strain. Currently allows for the following models:
#' \itemize{
#' \item Power transformations of observed proximities only (theta must be scalar): Strain loss or classical scaling (\code{strain}, workhorse is cmdscale), Kruskall's stress for symmetric matrices (\code{smacofSym} or \code{stress} and \code{smacofSphere} for scaling onto a sphere; workhorse is smacof), Sammon mapping (\code{sammon} or \code{sammon2}; for the earlier the workhorse is sammon from MASS for the latter it is smacof), elastic scaling (\code{elastic}, the workhorse is smacof), Takane et al's s-Stress \code{sstress} (workhorse is powerstressMin)
#' \item Power transformations of fitted distances only (theta must be scalar): De Leeuw's r-stress \code{rstress} (workhorse is powerstressMin)
#' \item Power transformations of fitted distances and observed proximities (theta must be scalar or of length 2): Power MDS (\code{powermds}), Sammon mapping/elastic scaling with powers (\code{powersammon}, \code{powerelastic})
#' \item Power transfomations of fitted distances, observed proximities and weights (theta must be of length 3 at most): powerstress (POST-MDS, \code{powerstress}), restricted powerstress with equal transformations for distances and proximities (\code{rpowerstress}); workhorse is powerstressMin)
#' \item Approximation to power stress (theta must be of length 2): Approximated power stress (\code{apstress}; workhorse is smacof)
#' }
#' @param theta the theta vector of powers; see the corresponding cop_XXX function for which theta are allowed. If a scalar is given as argument, it will be recycled. Defaults to 1.
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights; defaults to 1 for all off diagonals 
#' @param init (optional) initial configuration. If not supplied, the Torgerson scaling result of the dissimilarity matrix dis^theta[2]/enorm(dis^theta[2],weightmat) is used.
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param cordweight weight to be used for the cordillera; if missing gets estimated from the initial configuration so that copstress = 0 for theta=c(1,1) 
#' @param q the norm of the cordillera; defaults to 1
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to ndim+1
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the minimum reachabilities to be considered. If missing it is found from the initial configuration by taking 1.5 times the maximal minimum reachability of the model with theta=c(1,1). If NULL it will be normed to each configuration's minimum and maximum distance, so an absolute value of goodness-of-clusteredness. Note that the latter is not necessarily desirable when comparing configurations for their relative clusteredness. See also \code{\link{cordillera}}.     
#' @param optimmethod What general purpose optimizer to use? Defaults to our adaptive LJ version (ALJ). Also allows particle swarm optimization with s particles ("pso") and simulated annealing ("SANN"), "DIRECT" and "DIRECTL", Hooke-Jeeves ("hjk"), StoGo ("stogo"), and "MADS". We recommend not using SANN and pso with the rstress, sstress and the power stress models. We made good experiences with ALJ, stogo, DIRECT and DIRECTL and also MADS. 
#' @param lower A vector of the lower box contraints of the search region. Its length must match the length of theta.
#' @param upper A vector of the upper box contraints of the search region. Its length must match the length of theta. 
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose. Note that for models with some parameters fixed, the iteration progress of the optimizer shows different values also for the fixed parameters because due to the modular setup we always optimize over a three parameter vector. These values are inconsequential however as internally they will be fixed. 
#' @param normed should the cordillera be normed; defaults to TRUE
#' @param scale should the configuration be scaled and/or centered for calculating the cordillera? "std" standardizes each column of the configurations to mean=0 and sd=1 (typically not a good idea), "sd" scales the configuration by the maximum standard devation of any column (default), "proc" adjusts the fitted configuration to the init configuration (or the Togerson scaling solution if init=NULL). This parameter only has an effect for calculating the cordillera, the fitted and returned configuration is NOT scaled.     
#'@param s number of particles if pso is used
#'@param stresstype what stress to be used for comparisons between solutions. Currently not implemented and pcops uses explicitly normalized stress for copstress (not stress-1). Stress-1 is reported by the print function though.   
#'@param itmaxo iterations of the outer step (optimization over the hyperparmeters; if solver allows it). Defaults to 200.  
#'@param itmaxi iterations of the inner step (optimization of the MDS). Defaults to 10000 (whichis huge).
#'@param acc termination threshold difference of two successive outer minimization steps.
#'@param ... additional arguments to be passed to the optimization procedure
#'
#'@return A list with the components
#'         \itemize{
#'         \item copstress: the weighted loss value
#'         \item OC: the OPTICS cordillera for the scaled configuration (as defined by scale) 
#'         \item optim: the object returned from the optimization procedure
#'         \item stress: the stress (square root of stress.m)
#'         \item stress.m: default normalized stress 
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#'         \item cordillera: the cordillera object
#' }
#' 
#'@examples
#' dis<-as.matrix(smacof::kinshipdelta)
#' set.seed(210485)
#' #configuration is scaled with highest column sd for calculating cordilera 
#' res1<-pcops(dis,loss="strain",lower=0.1,upper=5,minpts=2) 
#' res1
#' summary(res1)
#' plot(res1)
#' 
#' 
#'@importFrom stats dist as.dist optim sd
#'@importFrom pso psoptim
#'@importFrom nloptr direct directL stogo 
#'@importFrom crs snomadr
#'@importFrom dfoptim hjk
#'@import cordillera
#' 
#'@keywords clustering multivariate
#'@export
pcops <- function(dis,loss=c("stress","smacofSym","smacofSphere","strain","sammon","rstress","powermds","sstress","elastic","powersammon","powerelastic","powerstress","sammon2","powerstrain","apstress","rpowerstress"),weightmat=NULL,ndim=2,init=NULL,theta=1,stressweight=1,cordweight,q=2,minpts=ndim+1,epsilon=100,rang,optimmethod=c("ALJ","pso","SANN","DIRECT","DIRECTL","stogo","MADS","hjk"),lower=0.5,upper=5,verbose=0,scale=c("proc", "sd", "none", "std"),normed=TRUE,s=4,stresstype="default",acc=1e-7,itmaxo=200,itmaxi=10000,...)
{
      if(missing(scale)) scale <- "sd"
      if(inherits(dis,"dist")) dis <- as.matrix(dis)
      if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1]) 
      if(missing(loss)) loss <- "strain"
      #if(length(theta)==1L) expo <- theta
      #if(length(theta)>2) expo <- theta[2]
      #if(is.null(init)) init <- cops::torgerson(dis^expo/enorm(dis^expo,weightmat),p=ndim) #like in powerstressMin
      .confin <- init #initialize a configuration
      psfunc <- switch(loss,"strain"=cop_cmdscale,"powerstrain"=cop_cmdscale,"elastic"=cop_elastic,"sstress"=cop_sstress,"stress"=cop_smacofSym,"smacofSym"= cop_smacofSym,"smacofSphere"=cop_smacofSphere,"rstress"=cop_rstress,"powermds"=cop_powermds,"powerstress"=cop_powerstress,"sammon"=cop_sammon,"sammon2"=cop_sammon2,"powersammon"=cop_powersammon,"powerelastic"=cop_powerelastic,"apstress"=cop_apstress,"rpowerstress"=cop_rpowerstress) #choose the stress to minimize
      if(missing(optimmethod)) optimmethod <- "ALJ"
      if(missing(rang)) 
          {
           if(verbose>1) cat ("Fitting configuration for rang. \n")    
           initsol <- do.call(psfunc,list(dis=dis,theta=1,init=.confin,weightmat=weightmat,ndim=ndim,rang=c(0,1),q=q,minpts=minpts,epsilon=epsilon,verbose=verbose-2,scale=scale,normed=normed,stresstype=stresstype,itmaxi=itmaxi))
           init0 <- initsol$fit$conf
           #if(scale=="std") init0 <- scale(init0)
           #if(scale=="sd") init0 <- init0/max(apply(init0,2,sd))
                                        #if(scale=="proc") init0 <- smacof::Procrustes(init,init0)$Yhat
           init0 <- scale_adjust(init0,.confin,scale=scale)
           crp <- cordillera::cordillera(init0,q=q,minpts=minpts,epsilon=epsilon,scale=FALSE)$reachplot
           cin <- max(crp)
           rang <- c(0,1.5*cin) #approximate upper bound by 1.5 times the max distance in the initial config
                 #alternatives: use an adjusted boxplot idea so e.g., rang<-c(quantile(crp,0.25)-exp(-4*robustbase::mc(crp))*1.5*IQR(crp),quantile(crp,0.75)+exp(4*robustbase::mc(crp))*1.5*IQR(crp)
                 #alternatives: use an adjusted boxplot idea so e.g., c(min(crp)-exp(-4*robustbase::mc(crp))*1.5,max(crp)+exp(4*robustbase::mc(crp))*1.5) 
            if(verbose>1) cat("dmax is",max(rang),". rang is",rang,"\n")
           }
      if(is.null(rang) && verbose > 1) cat("rang=NULL which makes the cordillera a goodness-of-clustering relative to the largest distance of each given configuration. \n") 
      if(missing(cordweight))
      {
          if(!exists("initsol")) {
                 if(verbose>1) cat ("Fitting configuration for cordweight. \n")     
                 initsol <- do.call(psfunc,list(dis=dis,theta=1,init=.confin,weightmat=weightmat,ndim=ndim,rang=rang,q=q,minpts=minpts,epsilon=epsilon,verbose=verbose-2,scale=scale,normed=normed,stresstype=stresstype,itmaxi=itmaxi))
          }
            init0 <- initsol$fit$conf
            init0 <- scale_adjust(init0,.confin,scale=scale)
            #if(scale=="std") init0 <- scale(init0)
            #if(scale=="sd") init0 <- init0/max(apply(init0,2,sd))
            #if(scale=="proc") init0 <- smacof::Procrustes(init,init0)$Yhat
            initcorrd <- cordillera::cordillera(init0,q=q,epsilon=epsilon,minpts=minpts,rang=rang,scale=FALSE)$normed 
            if(identical(normed,FALSE)) initcorrd <- cordillera::cordillera(init0,q=q,epsilon=epsilon,minpts=minpts,rang=rang,scale=scale)$raw
            cordweight <- initsol$stress.m/initcorrd  
            if(verbose>1) cat("Weights are stressweight=",stressweight,"cordweight=",cordweight,"\n")
            }
      if(verbose>1) cat("Starting Optimization \n ")
      if(optimmethod=="SANN") {
          opt<- stats::optim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,weightmat=weightmat,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,stresstype=stresstype,itmaxi=itmaxi))$copstress,method="SANN",control=list(maxit=itmaxo,trace=verbose-2,reltol=acc),...)
      }
      if(optimmethod=="pso") {
        addargs <- list(...)
        control <- list(trace=verbose-2,s=s,addargs)
        opt<- pso::psoptim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,weightmat=weightmat,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,stresstype=stresstype,itmaxi=itmaxi))$copstress,lower=lower,upper=upper,control=control)
        thetaopt <- opt$par
       }
      if(optimmethod=="ALJ") {
          opt<- cops::ljoptim(theta, function(theta) do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=theta,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,stresstype=stresstype,itmaxi=itmaxi))$copstress,lower=lower,upper=upper,verbose=verbose-2,itmax=itmaxo,acc=acc,...)
            # opt<- cops::ljoptim(theta, function(theta) do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=theta,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,stresstype=stresstype,itmaxi=itmaxi))$copstress,lower=lower,upper=upper,verbose=verbose-2,itmax=itmaxo,acc=acc)
            thetaopt <- opt$par
      }
      if(optimmethod=="DIRECT") {
          opt<- nloptr::direct(function(theta) do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=theta,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,stresstype=stresstype,itmaxi=itmaxi))$copstress,lower=lower,upper=upper,nl.info=isTRUE(all.equal(verbose-2,0)),control=list(maxeval=itmaxo,xtol_rel=acc),...)
            thetaopt <- opt$par
      }
       if(optimmethod=="stogo") {
           opt<- nloptr::stogo(theta,function(theta) do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=theta,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,stresstype=stresstype,itmaxi=itmaxi))$copstress,lower=lower,upper=upper,nl.info=isTRUE(all.equal(verbose-2,0)),maxeval=itmaxo,xtol_rel=acc,...)
             thetaopt <- opt$par
      }
      if(optimmethod=="DIRECTL") {
          opt<- nloptr::directL(function(theta) do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=theta,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,stresstype=stresstype,itmaxi=itmaxi))$copstress,lower=lower,upper=upper,nl.info=isTRUE(all.equal(verbose-2,0)),control=list(maxeval=itmaxo,xtol_rel=acc),...)
            thetaopt <- opt$par
       }
      if(optimmethod=="MADS") {
      #snomard is super stupid with extra parameters    
      eval.f.pars <- function(x,params)
      {
        psfunc <- params[[1]]
        tmplist <- list(theta=x)
        parlist <- c(tmplist,params[-1])
        tmpo <- do.call(psfunc,parlist)
        return(tmpo$copstress)
       }
       params <- list(psfunc,dis=dis,weightmat=weightmat,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-4,scale=scale,normed=normed,stresstype=stresstype,itmaxi=itmaxi)   
       opt<- crs::snomadr(n=length(theta),x0=theta,eval.f=eval.f.pars,params=params,bbin=0,lb=lower,ub=upper,print.output=isTRUE(all.equal(verbose-2,0)),opts=list("MAX_BB_EVAL"=itmaxo),...)
       thetaopt <- opt$solution
       opt$par <- opt$solution  
       }
      if(optimmethod=="hjk") {
          opt<- dfoptim::hjkb(theta, function(theta) do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=theta,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,stresstype=stresstype,itmaxi=itmaxi))$copstress,lower=lower,upper=upper,control=list(info=isTRUE(all.equal(verbose-2,0)),maxfeval=itmaxo,tol=acc),...)
       thetaopt <- opt$par
       } 
    #refit the optimal version (TODO probably unnecessary if the other functions are properly reimplemented)
    out <- do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=thetaopt,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-2,scale=scale,normed=normed,stresstype=stresstype,itmaxi=itmaxi))
    confopt <- scale_adjust(out$fit$conf,.confin,scale=scale)
    out$OC <- cordillera::cordillera(confopt,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=FALSE)
    #out$copstress <- opt$value 
    out$optim <- opt
    out$stressweight <- stressweight
    out$cordweight <- cordweight
    out$call <- match.call()
    out$optimethod <- optimmethod
    out$losstype <- loss
    out$nobj <- dim(out$fit$conf)[1]
    out$scale <- scale
    if(verbose>1) cat("Found minimum after",opt$counts["function"]," iterations at",round(thetaopt,4),"with copstress=",round(out$copstress,4),"and default scaling loss=",round(out$stress.m,4),"and OC=", round(out$OC$normed,4),". Thanks for your patience. \n")
    class(out) <- c("pcops","stops","cops")
    out
}


#'@export
summary.pcops <- function(object,...)
    {
      sppmat <- NULL
      if(!is.null(object$fit$spp))
      { 
           spp.perc <- object$fit$spp/sum(object$fit$spp) * 100
           sppmat <- cbind(sort(object$fit$spp), sort(spp.perc))
           colnames(sppmat) <- c("SPP", "SPP(%)")
      } 
      res <- list(conf=object$fit$conf,sppmat=sppmat)
      class(res) <- "summary.pcops"
      res
    }

#'@export
print.summary.pcops <- function(x,...)
    {
    cat("\n")
    cat("Configurations:\n")
    print(round(x$conf, 4))
    cat("\n\n")
    if(!is.null(x$sppmat))
     {   
      cat("Stress per point:\n")
      print(round(x$sppmat, 4))
      cat("\n")
     }
    }

#'@export
print.pcops <- function(x,...)
    {
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    cat("Model: P-COPS with", x$loss,"loss function and theta parameter vector =",x$parameters,"\n")
    cat("\n")
    cat("Number of objects:", x$nobj, "\n")
    cat("MDS loss value:", x$stress, "\n")
    cat("OPTICS Cordillera: Raw", x$OC$raw,"Normed", x$OC$normed,"\n")
    cat("Cluster optimized loss (copstress): ", x$copstress, "\n")
    cat("MDS loss weight:",x$stressweight," OPTICS Cordillera weight:",x$cordweight,"\n")
    cat("Number of iterations of",x$optimethod,"optimization:", x$optim$counts["function"], "\n")
    cat("\n")
    }


#'S3 plot method for p-cops objects
#' 
#'@param x an object of class cops
#'@param plot.type String indicating which type of plot to be produced: "confplot", "reachplot", "resplot","transplot", "Shepard", "stressplot" (see details)
#'@param main the main title of the plot
#'@param asp aspect ratio of x/y axis; defaults to NA; setting to 1 will lead to an accurate represenation of the fitted distances. 
#'@param ... Further plot arguments passed: see 'plot.smacof' and 'plot' for detailed information.
#' 
#'Details:
#' \itemize{
#' \item Configuration plot (plot.type = "confplot"): Plots the MDS configurations.
#' \item Reachability plot (plot.type = "confplot"): Plots the OPTICS reachability plot and the OPTICS cordillera 
#' \item Residual plot (plot.type = "resplot"): Plots the dissimilarities against the fitted distances.
#' \item Linearized Shepard diagram (plot.type = "Shepard"): Diagram with the transformed observed dissimilarities against the transformed fitted distance as well as loess smooth and a least squares line.
#' \item Transformation Plot (plot.type = "transplot"): Diagram with the observed dissimilarities (lighter) and the transformed observed dissimilarities (darker) against the fitted distances together with loess smoothing lines 
#' \item Stress decomposition plot (plot.type = "stressplot", only for SMACOF objects in $fit): Plots the stress contribution in of each observation. Note that it rescales the stress-per-point (SPP) from the corresponding smacof function to percentages (sum is 100). The higher the contribution, the worse the fit.
#' \item Bubble plot (plot.type = "bubbleplot", only available for SMACOF objects $fit): Combines the configuration plot with the point stress contribution. The larger the bubbles, the better the fit.
#'}
#'@export
#'@examples
#' dis<-as.matrix(smacof::kinshipdelta)
#' resl<-pcops(dis,loss="strain",lower=0.1,upper=5,minpts=2)
#' plot(resl)
#' plot(resl,plot.type="Shepard")
plot.pcops <- function(x,plot.type=c("confplot"), main, asp=NA,...)
    {
     if(missing(plot.type)) plot.type <- "confplot"  
     if(plot.type=="reachplot") {
        if(missing(main)) main <- paste("Reachability plot")
        else main <- main
        plot(x$OC,main=main,...)
     } else if(inherits(x$fit,"smacofB") && !inherits(x$fit,"smacofP") && plot.type=="transplot"){
     #ok, old code here has side effects: it changes the smacof object in the cops object; not sure we should do that  
       if(missing(main)) main <- paste("Transformation Plot") 
       plot(x$fit,plot.type="transplot",...)
     # invisible(tmp) #I give the changed smacof object back
     }
     else {      
       plot(x$fit,plot.type=plot.type,main=main,asp=asp,...)
   }
 }

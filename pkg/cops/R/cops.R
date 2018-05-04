#' PCOPS versions of smacofSym models
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this is either a scalar of the lambda transformation for the observed proximities, or a vector where the first is the kappa argument for the fitted distances (here internally fixed to 1) and the second the lambda argument and the third the nu argument (here internally fixed to 1). Defaults to 1 1 1 
#' @param ndim number of dimensions of the target space
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
#'         \item{fit:} the returned object of the fitting procedure (which has all smacofB elements and some more
#'         \item{cordillera:} the cordillera object
#' }
#'
#'@importFrom stats dist as.dist
#'@import cordillera 
#'@import smacof
#' 
#'@keywords multivariate
#'@export
cop_smacofSym <- function(dis,theta=c(1,1,1),ndim=2,weightmat=NULL,init=NULL,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,normed=TRUE,scale,stresstype="default") {
  #TODO Unfolding  
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
  #kappa first argument, lambda=second
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  lambda <- theta
  if(length(theta)==3L) lambda <- theta[2]
 # addargs <- list(...)
 # addargs
  fit <- smacof::smacofSym(dis^lambda,ndim=ndim,weightmat=weightmat,init=init,verbose=isTRUE(verbose==2),...) #optimize with smacof
  fit$kappa <- 1
  fit$lambda <- lambda
  fit$nu <- 1
  fit$stress.1 <- fit$stress
 # fit$stress <- (fit$stress^2)*sum(fit$obsdiss^2) check if this is like below
 # fitdis <- 2*sqrt(sqdist(fit$conf))
  fitdis <- as.matrix(fit$confdiss)
  delts <- as.matrix(fit$delta) #That was my choice to not use the normalized deltas but try it on the original; that is scale and unit free as Buja said
 # if(inherits(fit,"smacofSP")) delts <- as.matrix(fit$delta)[-1,-1]
  fit$stress.r <- sum(weightmat*(delts-fitdis)^2)
  fit$stress.m <- fit$stress.r/sum(weightmat*delts^2)
  fit$pars <- c(fit$kappa,fit$lambda,fit$nu)
  fit$deltaorig <- fit$delta^(1/fit$lambda)  
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit,cordillera=copobj) #target functions
  out
}

#' PCOPS versions of elastic scaling models (via smacofSym)
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this is either a scalar of the lambda transformation for the observed proximities, or a vector where the first is the kappa argument for the fitted distances (here internally fixed to 1) and the second the lambda argument and the third the nu argument (here internally fixed to -2). Defaults to 1 1 -2
#' @param ndim number of dimensions of the target space
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
#'@export
cop_elastic <- function(dis,theta=c(1,1,-2),ndim=2,weightmat=1,init=NULL,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,normed=TRUE,scale=3,stresstype="default") {
  #TODO Unfolding  
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1]) 
  #kappa first argument, lambda=second
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  lambda <- theta
  if(length(theta)==3L) lambda <- theta[2]
  nu <- -2
  addargs <- list(...)
  addargs
  elscalw <- dis^(nu*lambda) #the weighting in elastic scaling
  diag(elscalw) <- 1
  combwght <- weightmat*elscalw #combine the user weights and the elastic scaling weights
  fit <- smacof::smacofSym(dis^lambda,ndim=ndim,weightmat=combwght,init=init,verbose=isTRUE(verbose==2),...) #optimize with smacof
  fit$kappa <- 1
  fit$lambda <- lambda
  fit$nu <- nu
  fit$stress.1 <- fit$stress
 # fit$stress <- (fit$stress^2)*sum(fit$obsdiss^2) check if this is like below
 # fitdis <- 2*sqrt(sqdist(fit$conf))
  fitdis <- as.matrix(fit$confdiss)
  delts <- as.matrix(fit$delta) 
  fit$stress.r <- sum(combwght*((delts-fitdis)^2))
  fit$stress.m <- fit$stress.r/sum(combwght*delts^2)
  fit$pars <- c(kappa,lambda,nu)
  fit$deltaorig <- fit$delta^(1/fit$lambda)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit,cordillera=copobj) #target functions
  out
}


#' PCOPS versions of smacofSphere models
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this is either a scalar of the lambda transformation for the observed proximities, or a vector where the first is the kappa argument for the fitted distances (here internally fixed to 1) and the second the lambda argument and teh third the nu argument (here internally fixed to 1). Defaults to 1 1 1
#' @param ndim number of dimensions of the target space
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
#'@export
cop_smacofSphere <- function(dis,theta=c(1,1),ndim=2,weightmat=NULL,init=NULL,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,normed=TRUE,scale=3,stresstype="default") {
  #TODO Unfolding  
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
  #kappa first argument, lambda=second
  if(length(theta)>2) stop("There are too many parameters in the theta argument.")
  lambda <- theta
  if(length(theta)==2L) lambda <- theta[2]
  addargs <- list(...)
  addargs
  fit <- smacof::smacofSphere(dis^lambda,ndim=ndim,weightmat=weightmat,init=init,verbose=isTRUE(verbose==2),...) #optimize with smacof
  fit$kappa <- 1
  fit$lambda <- lambda
  fit$stress.1 <- fit$stress
 # fit$stress <- (fit$stress^2)*sum(fit$obsdiss^2) check if this is like below
 # fitdis <- 2*sqrt(sqdist(fit$conf))
  fitdis <- as.matrix(fit$confdiss)
  delts <- as.matrix(fit$delta)[-1,-1]
  fit$stress.r <- sum(weightmat*(delts-fitdis)^2)
  fit$stress.m <- fit$stress.r/sum(weightmat*delts^2)
  fit$pars <- c(kappa,lambda)
  fit$deltaorig <- fit$delta^(1/fit$lambda)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit,cordillera=copobj) #target functions
  out
}


#' PCOPS version of sammon mapping
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this is either a scalar of the lambda transformation for the observed proximities, or a vector where the first is the kappa argument for the fitted distances (here internally fixed to 1) and the second the lambda argument and the third the nu argument (here internally fixed to -1). Defaults to 1 1 -1 
#' @param ndim number of dimensions of the target space
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
#' @export
cop_sammon <- function(dis,theta=c(1,1,-1),ndim=2,init=NULL,weightmat=NULL,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale=3,normed=TRUE,stresstype="default") {
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(length(theta)==1L) lambda <- theta
  if(length(theta)==2L) lambda <- theta[2]
  if(length(theta)==3L) lambda <- theta[2]
  nu <- -1
 # if(is.null(init)) init <- cops::cmdscale(dis^lambda,k=ndim)$points
  fit <- cops::sammon(dis^lambda,k=ndim,y=init,trace=isTRUE(verbose>1),...)
  fit$lambda <- lambda
  fit$kappa <- 1
  fit$nu <- -1
  dis <- stats::as.dist(dis)
  fitdis <- stats::dist(fit$points)
  fit$stress.r <- sum(((dis^lambda-fitdis)^2)/dis)
  fit$stress.n <- fit$stress.r/sum(dis)
  fit$stress.m <- sqrt(fit$stress)
  fit$conf <- fit$points
#  fit$pars <- c(kappa,lambda,nu)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters,  fit=fit,cordillera=copobj) #target functions
}

#' COPS versions of Sammon mapping models (via smacofSym)
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this is either a scalar of the lambda transformation for the observed proximities, or a vector where the first is the kappa argument for the fitted distances (here internally fixed to 1) and the second the lambda argument, the thrid the nu argiument (fixed to -1). Defaults to 1 1 -1. 
#' @param ndim number of dimensions of the target space
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
#'@export
cop_sammon2 <- function(dis,theta=c(1,1,-1),ndim=2,weightmat=NULL,init=NULL,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,normed=TRUE,scale=3,stresstype="default") {
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1]) 
  #kappa first argument, lambda=second
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  lambda <- theta
  if(length(theta)==3L) lambda <- theta[2]
  nu <- -1
  addargs <- list(...)
  elscalw <- dis^(nu*lambda) #the weighting in elastic scaling
  diag(elscalw) <- 1
  combwght <- weightmat*elscalw #combine the user weights and the elastic scaling weights
  fit <- smacof::smacofSym(dis^lambda,ndim=ndim,weightmat=combwght,init=init,verbose=isTRUE(verbose==2),...) #optimize with smacof
  fit$kappa <- 1
  fit$lambda <- lambda
  fit$nu <- nu
  fit$stress.1 <- fit$stress
 # fit$stress <- (fit$stress^2)*sum(fit$obsdiss^2) check if this is like below
 # fitdis <- 2*sqrt(sqdist(fit$conf))
  fitdis <- as.matrix(fit$confdiss)
  delts <- as.matrix(fit$delta) 
  fit$stress.r <- sum(combwght*((delts-fitdis)^2))
  fit$stress.m <- fit$stress.r/sum(combwght*delts^2)
  fit$pars <- c(kappa,lambda,nu)
  fit$deltaorig <- fit$delta^(1/fit$lambda)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit,cordillera=copobj) #target functions
  out
}



#' PCOPS version of strain
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this is either a scalar of the lambda transformation for the observed proximities, or a vector where the first is the kappa argument for the fitted distances (here internally fixed to 1) and the second and third the lambda and the nu argument (the latter is fixed to 1). Defaults to 1 1 1
#' @param ndim number of dimensions of the target space
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
#' @export
cop_cmdscale <- function(dis,theta=c(1,1,1),weightmat=NULL,ndim=2,init=NULL,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale=3,normed=TRUE,stresstype="default") {
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)==1L) lambda <- theta
  if(length(theta)==2L) lambda <- theta[2]
  if(length(theta)==3L) lambda <- theta[2]
  fit <- cops::cmdscale(dis^lambda,k=ndim,eig=TRUE,...) 
  fit$lambda <- lambda
  fit$kappa <- 1
  fit$nu <- 1
  dis <- stats::as.dist(dis)
  fitdis <- stats::dist(fit$points)
  fit$stress.r <- sum((dis^lambda-fitdis)^2)
  fit$stress.n <- fit$stress.r/sum(dis^(2*lambda))
  fit$stress.m <- sqrt(fit$stress.n)
  fit$conf <- fit$points
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  list(stress=fit$GOF,stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit,cordillera=copobj) #target functions
}

#' PCOPS version of rstress
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this is either a scalar of the kappa transformation for the fitted distances proximities, or a vector where the first is the kappa argument for the fitted distances and the second the lambda argument, the third the nu argument (here internally fixed to 1). Defaults to 1 1 1 
#' @param ndim number of dimensions of the target space
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
#' @export
cop_rstress <- function(dis,theta=c(1,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale=3,normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
  if(missing(stresstype)) stresstype <- "default"  
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  kappa <- theta
  if(length(theta)==3L) kappa <- theta[1] 
  fit <- powerStressMin(delta=dis,kappa=kappa,lambda=1,nu=1,weightmat=weightmat,init=init,ndim=ndim,verbose=verbose,...)
  if(stresstype=="default") fit$stress.m <- fit$stress.m
  if(stresstype=="stress1") fit$stress.m <- fit$stress.1
  if(stresstype=="rawstress") fit$stress.m <- fit$stress.r
  if(stresstype=="enormstress") fit$stress.m <- fit$stress.en
  if(stresstype=="enormstress1") fit$stress.m <- fit$stress.en1
  if(stresstype=="normstress") fit$stress.m <- fit$stress.n
  fit$kappa <- theta[1]
  fit$lambda <- 1
  fit$nu <- 1
 # fit$pars <- c(kappa,lambda)
 # fit$deltaorig <- fit$delta^(1/fit$lambda)
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit,cordillera=copobj)
  out
}


#' PCOPS version of sstress
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this is either a scalar of the lambda transformation for the observed proximities, or a vector where the first is the kappa argument for the fitted distances (here internally fixed to 2) and the second the lambda argument and the third the nu argument (internally fixed to 1). Defaults to 2 1 1
#' @param ndim number of dimensions of the target space
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
#' @export
cop_sstress <- function(dis,theta=c(2,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale=3,normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
  if(missing(stresstype)) stresstype <- "default"  
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  lambda <- theta
  if(length(theta)==3L) lambda <- theta[2]
  flambda <- lambda*2 #sstress is d^2 and delta^2 so f(delta^2)=delta^(2*1); lambda works in factors of 2  
  fit <- powerStressMin(delta=dis,kappa=2,lambda=flambda,nu=1,weightmat=weightmat,init=init,ndim=ndim,verbose=verbose,...)
  if(stresstype=="default") fit$stress.m <- fit$stress.m
  if(stresstype=="stress1") fit$stress.m <- fit$stress.1
  if(stresstype=="rawstress") fit$stress.m <- fit$stress.r
  if(stresstype=="normstress") fit$stress.m <- fit$stress.n
  if(stresstype=="enormstress") fit$stress.m <- fit$stress.en
  if(stresstype=="enormstress1") fit$stress.m <- fit$stress.en1
  fit$kappa <- 2
  fit$lambda <- lambda
  fit$nu <- 1
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit,cordillera=copobj)
  out
}


#' PCOPS version of powermds
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances), the second lambda (for the observed proximities), nu is fixed to 1. If a scalar is given it is recycled.  Defaults to 1 1 1.
#' @param ndim number of dimensions of the target space
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
#' @export
cop_powermds <- function(dis,theta=c(1,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale=3,normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
  if(missing(stresstype)) stresstype <- "default"
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)==1L) theta <- rep(theta,2)
  fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=1,weightmat=weightmat,init=init,ndim=ndim,verbose=verbose,...)
  if(stresstype=="default") fit$stress.m <- fit$stress.m
  if(stresstype=="stress1") fit$stress.m <- fit$stress.1
  if(stresstype=="rawstress") fit$stress.m <- fit$stress.r
  if(stresstype=="normstress") fit$stress.m <- fit$stress.n
  if(stresstype=="enormstress") fit$stress.m <- fit$stress.en
  if(stresstype=="enormstress1") fit$stress.m <- fit$stress.en1
  fit$kappa <- theta[1]
  fit$lambda <- theta[2]
  fit$nu <- 1
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit, cordillera=copobj)
  out 
}

#' PCOPS version of sammon with powers
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances), the second lambda (for the observed proximities), the third nu (fixed to -1). If a scalar is given it is recycled for the free parameters.  Defaults to 1 1 -1.
#' @param ndim number of dimensions of the target space
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
#' @export
cop_powersammon <- function(dis,theta=c(1,1,-1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale=3,normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
  if(missing(stresstype)) stresstype <- "default"
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)==1L) theta <- rep(theta,2)
  nu <- -1
  sammwght <-dis^(theta[2])
  diag(sammwght) <- 1
  combwght <- sammwght*weightmat
  fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=nu,weightmat=combwght,init=init,ndim=ndim,verbose=verbose,...)
  if(stresstype=="default") fit$stress.m <- fit$stress.m
  if(stresstype=="stress1") fit$stress.m <- fit$stress.1
  if(stresstype=="rawstress") fit$stress.m <- fit$stress.r
  if(stresstype=="normstress") fit$stress.m <- fit$stress.n
  if(stresstype=="enormstress") fit$stress.m <- fit$stress.en
  if(stresstype=="enormstress1") fit$stress.m <- fit$stress.en1
  fit$kappa <- theta[1]
  fit$lambda <- theta[2]
  fit$nu <- nu
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit, cordillera=copobj)
  out 
}

#' PCOPS version of elastic scaling with powers
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances), the second lambda (for the observed proximities) and nu as the third (fixed to -2). If a scalar for the free parameters is given it is recycled.  Defaults to 1 1 -2.
#' @param ndim number of dimensions of the target space
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
#' @export
cop_powerelastic <- function(dis,theta=c(1,1,-2),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale=3,normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
  if(missing(stresstype)) stresstype <- "default"
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)==1L) theta <- rep(theta,2)
  nu <- -2
  elawght <- dis^(theta[2])
  diag(elawght) <- 1
  combwght <- elawght*weightmat
  fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=nu,weightmat=combwght,init=init,ndim=ndim,verbose=verbose,...)
  if(stresstype=="default") fit$stress.m <- fit$stress.m
  if(stresstype=="stress1") fit$stress.m <- fit$stress.1
  if(stresstype=="rawstress") fit$stress.m <- fit$stress.r
  if(stresstype=="normstress") fit$stress.m <- fit$stress.n
  if(stresstype=="enormstress") fit$stress.m <- fit$stress.en
  if(stresstype=="enormstress1") fit$stress.m <- fit$stress.en1
  fit$kappa <- theta[1]
  fit$lambda <- theta[2]
  fit$nu <- nu
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit, cordillera=copobj)
  out 
}

#' COPS version of powerstress
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances), the second lambda (for the observed proximities), the third nu (for the weights). If a scalar is given it is recycled.  Defaults to 1 1 1.
#' @param ndim number of dimensions of the target space
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
#' @export
cop_powerstress <- function(dis,theta=c(1,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale=2,normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
  if(missing(stresstype)) stresstype <- "default"
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)==1L) theta <- rep(theta,3)
  wght <- weightmat
  diag(wght) <- 1
  fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=theta[3],weightmat=wght,init=init,ndim=ndim,verbose=verbose,...)
  if(stresstype=="default") fit$stress.m <- fit$stress.m
  if(stresstype=="stress1") fit$stress.m <- fit$stress.1
  if(stresstype=="rawstress") fit$stress.m <- fit$stress.r
  if(stresstype=="normstress") fit$stress.m <- fit$stress.n
  if(stresstype=="enormstress") fit$stress.m <- fit$stress.en
  if(stresstype=="enormstress1") fit$stress.m <- fit$stress.en1
  fit$kappa <- theta[1]
  fit$lambda <- theta[2]
  fit$nu <- theta[3]
  copobj <- copstress(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed,init=init)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, copstress=copobj$copstress, OC=copobj$OC, parameters=copobj$parameters, fit=fit, cordillera=copobj)
  out 
}

#' Calculates copstress for given MDS object 
#'
#' @param obj MDS object (supported are sammon, cmdscale, smacof, rstress, powermds)
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param cordweight weight to be used for the cordillera; defaults to 0.5
#' @param q the norm of the corrdillera; defaults to 1
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
#' @export
copstress <- function(obj,stressweight=1,cordweight=5,q=1,minpts=2,epsilon=10,rang=NULL,verbose=0,normed=TRUE,scale=c("std","sd","proc","none"),init,...)
{
        if(missing(scale)) scale <- "sd"
        stressi <- obj$stress.m
        kappa <- obj$kappa
        lambda <- obj$lambda
        nu <- obj$nu
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
        if(verbose>0) cat("copstress =",ic,"mdsloss =",stressi,"OC =",struc,"kappa =",kappa,"lambda =",lambda,"nu=",nu,"\n")
        out <- list(copstress=ic,OC=struc,parameters=c(kappa=kappa,lambda=lambda,nu=nu),cordillera=corrd)
        out
     }

#' Profile COPS Function (aka COPS Variant 2)
#'
#' Metaparameter selection for MDS models baseed on the Profile COPS approach (COPS Variant 2). It uses copstress for hyperparameter selection. It is a special case of a STOPS model.  
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param loss which loss function to be used for fitting, defaults to strain. Currently allows for the following models:
#' \itemize{
#' \item Power transformations of observed proximities only: Strain loss or classical scaling (\code{strain}, workhorse is cmdscale), Kruskall's stress for symmetric matrices (\code{smacofSym} or \code{stress} and \code{smacofSphere} for scaling onto a sphere; workhorse is smacof), Sammon mapping (\code{sammon} or \code{sammon2}; for the earlier the workhorse is sammon from MASS for the latter it is smacof), elastic scaling (\code{elastic}, the workhorse is smacof), Takane et al's S-Stress \code{sstress} (workhorse is powerStressMin)
#' \item Power transformations of fitted distances only: De Leeuw's r-stress \code{rstress} (workhorse is powerStressMin)
#' \item Power transformations of fitted distances and observed proximities: Powermds \code{powermds}, Sammon mapping and elastic scaling with powers (\code{powersammon}, \code{powerelastic}; workhorse is powerStressMin)
#' }
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances if it exists), the second lambda (for the observed proximities if it exist), the third is nu (for the weights if it exists) . If a scalar is given as argument, it will take the role designated by the loss argument. Defaults to 1 1 1
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights; defaults to 1 for all off diagonals 
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param cordweight weight to be used for the cordillera; if missing gets estimated from the initial configuration so that copstress = 0 for theta=c(1,1) 
#' @param q the norm of the corrdillera; defaults to 1
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to ndim+1
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the minimum reachabilities to be considered. If missing it is found from the initial configuration by taking 1.5 times the maximal minimum reachability of the model with theta=c(1,1). If NULL it will be normed to each configuration's minimum and maximum distance, so an absolute value of goodness-of-clusteredness. Note that the latter is not necessarily desirable when comparing configurations for their relative clusteredness. See also \code{\link{cordillera}}     
#' @param optimmethod What general purpose optimizer to use? Defaults to our adaptive LJ version (ALJ). Also allows particle swarm optimization with s particles ("pso") and simulated annealing ("SANN"). We recommend not using the latter with the rstress, sstress and the power stress models. 
#' @param lower The lower contraints of the search region
#' @param upper The upper contraints of the search region 
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param normed should the cordillera be normed; defaults to TRUE
#' @param scale should the configuration be scaled and/or centered for calculating the cordillera? "std" standardizes each column of the configurations to mean=0 and sd=1, "sd" scales the configuration by the maximum standard devation of any column, "proc" adjusts the fitted configuration to the init configuration (or the Togerson scaling solution if init=NULL). This parameter only has an effect for calculating the cordillera, the fitted and returned configuration is NOT scaled.     
#'@param s number of particles if pso is used
#'@param stresstype what stress to be used for comparisons between solutions 
#'@param ... additional arguments to be passed to the optimization procedure
#'
#'@return A list with the components
#'         \itemize{
#'         \item copstress: the weighted loss value
#'         \item OC: the Optics cordillera
#'         \item optim: the object returned from the optimization procedure
#'         \item stress: the stress
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
#' #only use cordillera for finding the parameters; 
#' res2<-pcops(dis,loss="strain",stressweight=0,lower=0.1,upper=5,minpts=2) 
#' res2
#' summary(res2)
#' plot(res2)
#'
#' #With Procrustes adjustment to res1 
#' res3<-pcops(dis,loss="strain",stressweight=0,lower=0.1,upper=5,minpts=2,
#' init=res1$conf,scale="proc") 
#' res3
#' summary(res3)
#' plot(res3)
#' 
#' par(mfrow=c(1,3))
#' plot(res1,"reachplot")
#' plot(res2,"reachplot")
#' plot(res3,"reachplot") 
#' par(mfrow=c(1,1))
#'
#'@importFrom stats dist as.dist optim sd
#'@importFrom pso psoptim
#'@import cordillera
#' 
#'@keywords clustering multivariate
#'@export
pcops <- function(dis,loss=c("stress","smacofSym","smacofSphere","strain","sammon","rstress","powermds","sstress","elastic","powersammon","powerelastic","powerstress","sammon2","powerstrain"),weightmat=NULL,ndim=2,init=NULL,theta=c(1,1,1),stressweight=1,cordweight,q=2,minpts=ndim+1,epsilon=100,rang,optimmethod=c("ALJ","pso","SANN"),lower=c(1,1,0.5),upper=c(5,5,2),verbose=0,scale=c("proc", "sd", "none", "std"),normed=TRUE,s=4,stresstype="default",...)
{
      if(missing(scale)) scale <- "sd"
      if(inherits(dis,"dist")) dis <- as.matrix(dis)
      if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1]) 
      if(missing(loss)) loss <- "strain"
      if(is.null(init)) init <- torgerson(dis,p=ndim)
      .confin <- init #initialize a configuration
      psfunc <- switch(loss,"strain"=cop_cmdscale,"powerstrain"=cop_cmdscale,"elastic"=cop_elastic,"sstress"=cop_sstress,"stress"=cop_smacofSym,"smacofSym"= cop_smacofSym,"smacofSphere"=cop_smacofSphere,"rstress"=cop_rstress,"powermds"=cop_powermds,"powerstress"=cop_powerstress,"sammon"=cop_sammon,"sammon2"=cop_sammon2,"powersammon"=cop_powersammon,"powerelastic"=cop_powerelastic) #choose the stress to minimize
      if(missing(optimmethod)) optimmethod <- "ALJ"
      if(missing(rang)) 
          {
           if(verbose>1) cat ("Fitting configuration for rang. \n")    
           initsol <- do.call(psfunc,list(dis=dis,theta=c(1,1,1),init=.confin,weightmat=weightmat,ndim=ndim,rang=c(0,1),q=q,minpts=minpts,epsilon=epsilon,verbose=verbose-2,scale=scale,normed=normed,stresstype=stresstype))
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
                 initsol <- do.call(psfunc,list(dis=dis,theta=c(1,1,1),init=.confin,weightmat=weightmat,ndim=ndim,rang=rang,q=q,minpts=minpts,epsilon=epsilon,verbose=verbose-2,scale=scale,normed=normed,stresstype=stresstype))
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
          opt<- stats::optim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,weightmat=weightmat,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,stresstype=stresstype))$copstress,method="SANN",...)
      }
      if(optimmethod=="pso") {
        addargs <- list(...)
        control <- list(trace=verbose-2,s=s,addargs)
        opt<- pso::psoptim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,weightmat=weightmat,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,stresstype=stresstype))$copstress,lower=lower,upper=upper,control=control)
       }
      if(optimmethod=="ALJ") {
      opt<- cops::ljoptim(theta, function(theta) do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=theta,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,stresstype=stresstype))$copstress,lower=lower,upper=upper,verbose=verbose-2,...)
      }
    thetaopt <- opt$par 
    #refit the optimal version (TODO probably unnecessary if the other functions are properly reimplemented)
    out <- do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=thetaopt,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-2,scale=scale,normed=normed,stresstype=stresstype))
    confopt <- scale_adjust(out$fit$conf,.confin,scale=scale)
    out$OC <- cordillera::cordillera(confopt,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=FALSE)
    out$copstress <- opt$value
    out$optim <- opt
    out$stressweight <- stressweight
    out$cordweight <- cordweight
    out$call <- match.call()
    out$optimethod <- optimmethod
    out$losstype <- loss
    out$nobj <- dim(out$fit$conf)[1]
    out$scale <- scale
    if(verbose>1) cat("Found minimum after",opt$counts["function"]," iterations at",round(opt$par,4),"with copstress=",round(out$copstress,4),"and default scaling loss=",round(out$stress.m,4),"and OC=", round(out$OC$normed,4),". Thanks for your patience. \n")
    class(out) <- c("pcops","stops","cops")
    out
}

#' Adjusts a configuration
#' 
#'@param conf a configuration
#'@param ref a reference configuration (only for scale="proc") 
#'@param scale Scale adjustment. "std" standardizes each column of the configurations to mean=0 and sd=1, "sd" scales the configuration by the maximum standard devation of any column, "proc" adjusts the fitted configuration to the reference
#' @return The scale adjusted configuration.
#'
#' @importFrom stats sd
#' @export
scale_adjust <- function(conf,ref,scale=c("sd","std","proc","none"))
{
  #conf target, ref reference
  if(scale=="std") conf <- scale(conf)
  if(scale=="sd") conf <- conf/max(apply(conf,2,stats::sd))
  if(scale=="proc") conf <- smacof::Procrustes(ref,conf)$Yhat
  if(scale=="none") conf <- conf
  conf
}


#'@export
print.cops <- function(x,...)
    {
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    cat("Model: COPS with parameters kappa=",x$par[1],"lambda=",x$par[2],"nu=",x$par[3],"\n")
    cat("\n")
    cat("Number of objects:", x$nobj, "\n")
    cat("Stress of configuration (default normalization):", x$stress, "\n")
    cat("OPTICS Cordillera: Raw", x$OC$raw,"Normed", x$OC$normed,"\n")
    cat("Cluster optimized loss (copstress): ", x$copstress, "\n")
    cat("Stress weight:",x$stressweight," OPTICS Cordillera weight:",x$cordweight,"\n")
    cat("Number of iterations of",x$optimethod,"optimization:", x$niter, "\n")
    cat("\n")
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
    cat("Model: COPS with", x$loss,"loss function and parameters kappa=",x$par[1],"lambda=",x$par[2],"nu=",x$par[3],"\n")
    cat("\n")
    cat("Number of objects:", x$nobj, "\n")
    cat("MDS loss value:", x$stress.m, "\n")
    cat("OPTICS Cordillera: Raw", x$OC$raw,"Normed", x$OC$normed,"\n")
    cat("Cluster optimized loss (copstress): ", x$copstress, "\n")
    cat("MDS loss weight:",x$stressweight," OPTICS Cordillera weight:",x$cordweight,"\n")
    cat("Number of iterations of",x$optimethod,"optimization:", x$optim$counts["function"], "\n")
    cat("\n")
    }

#'@export
#'@importFrom stats coef
coef.cops <- function(object,...)
    {
    return(c(kappa=object$par[1],lambda=object$par[2],nu=object$par[3]))
    }


#'S3 plot method for cops objects
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

#'S3 plot method for cops objects
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
plot.cops <- function(x,plot.type=c("confplot"), main, asp=1,...)
    {
     if(missing(plot.type)) plot.type <- "confplot"  
     if(plot.type=="reachplot") {
        if(missing(main)) main <- paste("Reachability plot")
        else main <- main
        plot(x$OC,main=main,...)
     } else if(inherits(x$fit,"smacofB") && !inherits(x$fit,"smacofP") && plot.type=="transplot"){
       if(missing(main)) main <- paste("Transformation Plot") 
       plot.smacofP(x,plot.type="transplot",asp=asp,...)
     # invisible(tmp) #I give the changed smacof object back
     }
     else {      
       plot.smacofP(x,plot.type=plot.type,main=main,asp=asp,...)
   }
 }

#' An old ratio MDS cops. Is going to be removed, only here for testing.
#' 
#' @param delta numeric matrix or dist object of a matrix of proximities
#' @param kappa power transformation for fitted distances
#' @param lambda power transformation for proximities
#' @param nu power transformation for weights
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances if it exists), the second lambda (for the observed proximities if it exist), the third is nu (for the weights if it exists) . If less than three elements are is given as argument, it will be recycled. Defaults to 1 1 1. Will override any kappa, lmabda, nu parameters if they are given and do not match
#' @param weightmat (optional) a matrix of nonnegative weights; defaults to 1 for all off diagonals
#' @param ndim number of dimensions of the target space
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 0.975
#' @param cordweight weight to be used for the cordillera; defaults to 0.025
#' @param q the norm of the corrdillera; defaults to 1
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to ndim+1
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param dmax The winsorization limit of reachability distances in the OPTICS Cordillera. If supplied, it should be either a numeric value that matches max(rang) or NULL; if NULL it is found as 1.5 times (for kappa >1) or 1 times (for kappa <=1) the maximum reachbility value of the power torgerson model with the same lambda. If dmax and rang are supplied and dmax is not max(rang), a warning is given and rang takes precedence.   
#' @param rang range of the reachabilities to be considered. If missing it is found from the initial configuration by taking 0 as the lower boundary and dmax (see above) as upper boundary. See also \code{\link{cordillera}}     
#' @param optimmethod What optimizer to use? Defaults to NEWUOA, Nelder-Mead is also supported.
#' @param verbose numeric value hat prints information on the fitting process; >2 is very verbose
#' @param normed should the cordillera be normed; defaults to TRUE
#' @param scale Allows to scale the configuration for the OC (the scaled configuration is also returned as $conf). One of none (so no scaling), sd (configuration divided by the maximum standard deviation of the columns), std (standardize all columns !NOTE: This does not preserve the relative distances of the optimal config), proc (procrustes adjustment to the initial fit) and rmsq (configuration divided by the maximum root mean square of the columns). Default is sd.   
#' @param accuracy numerical accuracy, defaults to 1e-8
#' @param itmax maximum number of iterations. Defaults to 100000
#' @param stresstype which stress to use in the copstress. Defaults to stress-1. If anything else is set, explicitly normed stress which is (stress-1)^2. Using stress-1 puts more weight on MDS fit.   
#' @param ... additional arguments to be passed to the optimization procedure
#' 
#' @import cordillera
#' @importFrom stats dist as.dist optim sd
#' @importFrom minqa newuoa
#'
#' @export
copstressMinOLD <- function (delta, kappa=1, lambda=1, nu=1, theta=c(kappa,lambda,nu),weightmat=1-diag(nrow(delta)),  ndim = 2, init=NULL, stressweight=0.975,cordweight=0.025,q=2,minpts=ndim+1,epsilon=10,dmax=NULL,rang,optimmethod=c("Nelder-Mead","Newuoa"),verbose=0,scale=c("sd","rmsq","std","proc","none"),normed=TRUE, accuracy = 1e-7, itmax = 100000, stresstype=c("stress-1","stress"),...)
{
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    kappa <- theta[1]
    lambda <- theta[2]
    nu <- theta[3]
    plot <- FALSE
    if(verbose>0) cat("Minimizing copstress with kappa=",kappa,"lambda=",lambda,"nu=",nu,".\n")
    if(missing(optimmethod)) optimmethod <- "Newuoa"
    if(missing(scale)) scale <- "sd"
    if(missing(stresstype)) stresstype <- "stress-1"
    if(missing(rang))
        #perhaps put this into the optimization function?
    {
        if(is.null(dmax))
           {
           if(verbose>1) cat ("Fitting configuration for rang. \n")    
           #initsol <- cops::powerstressFast(delta,kappa=kappa,lambda=lambda,nu=nu,weightmat=weightmat,ndim=ndim)
           initsol <- torgerson(delta^lambda,p=ndim)
           #init0 <- initsol$conf
           init0 <- initsol
           init0 <- init0/enorm(init0)
           if(scale=="std") init0 <- scale(init0) #standardizes config before cordillera
           if(scale=="none") init0 <- init0
           if(scale=="sd") #scales config to sd=1 for most spread dimension before cordillera
             {
                init0 <- init0/max(apply(init0,2,stats::sd))
             }   
             if(scale=="rmsq") #scales config to rmsq=1 for most spread dimension before cordillera
             {
                 testso <- scale(init0,center=FALSE)
                 init0 <- init0/max(attr(testso,"scaled:scale"))
             }
             if(scale=="proc") #scales config by procrusting to init
             {
                 if(!exists("init")) init <- initsol
                 procr <- smacof::Procrustes(init,init0)
                 init0 <- procr$Yhat
             }
           crp <- cordillera::cordillera(init0,q=q,minpts=minpts,epsilon=epsilon,scale=FALSE)$reachplot
           cin <- max(crp)
           dmax <- ifelse(kappa>1,1.5*cin,1.1*cin)
           #if(verbose > 2) cat("rang, dmax was NULL or missing which makes the cordillera a goodness-of-clustering relative to the largest distance of each given configuration. \n")
           }
        rang <- c(0,dmax) #This sets the range to (0,dmax)
        if(verbose>1) cat("dmax is",dmax,". rang is",rang,".\n")
    }
    if(isTRUE(max(rang)!=dmax)) warning("Note: The supplied dmax and rang do not match. I took supplied rang as rang.\n") 
    r <- kappa/2
    deltaorig <- delta
    delta <- delta^lambda
    weightmato <- weightmat
    weightmat <- weightmat^nu
    weightmat[!is.finite(weightmat)] <- 1 #new
    deltaold <- delta
    delta <- delta / enorm (delta, weightmat) #sum=1
    if(is.null(init))
    {
        if(exists("init0")) init <- init0 else init <- cops::powerStressFast(delta,kappa=kappa,lambda=lambda,nu=nu,ndim=ndim)$conf
    }
    xold <- init
    xold <- xold/enorm(xold) 
    copsf <- function(x,delta,r,ndim,weightmat,stressweight,cordweight,q,minpts,epsilon,rang,scale,normed,init,...)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=ndim)
             delta <- delta/enorm(delta,weightmat)
             x <- x/enorm(x)
             dnew <- sqdist (x)
             rnew <- sum (weightmat * delta * mkPower (dnew, r))
             nnew <- sum (weightmat * mkPower (dnew,  2*r))
             anew <- rnew / nnew
             stressi <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
             if(stresstype=="stress-1") stressi <- sqrt(stressi)
             if(scale=="none") x <- x 
             if(scale=="std") x <- scale(x) #standardizes config before cordillera
             if(scale=="sd") #scales config to sd=1 for most spread dimension before cordillera
             {
                x <- x/max(apply(x,2,stats::sd))
             }   
             if(scale=="rmsq") #scales config to rmsq=1 for most spread dimension before cordillera
             {
                 testso <- scale(x,center=FALSE)
                 x <- x/max(attr(testso,"scaled:scale"))
             }
             if(scale=="proc") #scales config by procrusting to init
             {
                 procr <- smacof::Procrustes(init,x)
                 x <- procr$Yhat
             }
             corrd <- cordillera::cordillera(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=FALSE,...)
             struc <- corrd$raw
             if(normed) {
                        struc <- corrd$normed
                       }
             ic <- stressweight*stressi - cordweight*struc
             if(verbose>1) cat("copstress =",ic,"mdsloss =",stressi,"OC =",struc,"minpts=",minpts,"kappa =",kappa,"lambda =",lambda,"nu=",nu,"\n")
             ic
           }
     if(verbose>1) cat("Starting Minimization with",optimmethod,":\n")
     if(optimmethod=="Newuoa") {
         optimized <- minqa::newuoa(xold,function(par) copsf(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose-2),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- optimized$feval
         ovalue <-optimized$fval
     }
     if(optimmethod=="Nelder-Mead") {
         optimized <- optim(xold,function(par) copsf(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmax,trace=verbose-2),...)
         xnew <- optimized$par
         itel <- optimized$counts[[1]]
         ovalue <-optimized$val 
     }
     xnew <- xnew/enorm(xnew)
     dnew <- sqdist (xnew)
     rnew <- sum (weightmat * delta * mkPower (dnew, r))
     nnew <- sum (weightmat * mkPower (dnew,  2*r))
     anew <- rnew / nnew
     stress <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
     if(stresstype=="stress-1") stress <- sqrt(stress)
     attr(xnew,"dimnames")[[1]] <- rownames(delta)
     attr(xnew,"dimnames")[[2]] <- paste("D",1:ndim,sep="")
     doutm <- (2*sqrt(sqdist(xnew)))^kappa  #fitted powered euclidean distance but times two
     #doutm <- as.matrix(dist(xnew)^kappa)
     deltam <- delta
     deltaorigm <- deltaorig
     deltaoldm <- deltaold
     delta <- stats::as.dist(delta)
     deltaorig <- stats::as.dist(deltaorig)
     deltaold <- stats::as.dist(deltaold)
     doute <- doutm/enorm(doutm)
     doute <- stats::as.dist(doute)
     dout <- stats::as.dist(doutm)
     resmat <- as.matrix(delta - doute)^2
     spp <- colMeans(resmat)
     weightmatm <-weightmat
     weightmat <- stats::as.dist(weightmatm)
     stressen <- sum(weightmat*(doute-delta)^2)#raw stress on the normalized proximities and normalized distances
     if(scale=="std") xnews <- scale(xnew) #standardizes config before cordillera
     if(scale=="sd") #scales config to sd=1 for most spread dimension before cordillera
             {
                xnews <- xnew/max(apply(xnew,2,stats::sd))
             }   
             if(scale=="rmsq") #scales config to rmsq=1 for most spread dimension before cordillera
             {
                 testso <- scale(xnew,center=FALSE)
                 xnews <- xnew/max(attr(testso,"scaled:scale"))
             }
             if(scale=="proc") #scales config by procrusting to init
             {
                 procr <- smacof::Procrustes(init,xnew)
                 xnews <- procr$Yhat
             }
             if(scale=="none") xnews <- xnew #no standardisation
      if(verbose>0) cat("*** stress (both normalized - for COPS/STOPS):",stress,"; stress 1 (both normalized - default reported):",sqrt(stress),"; stress manual (for debug only):",stressen,"; from optimization: ",ovalue,"\n")   
    out <- list(delta=deltaold, obsdiss=delta, confdiss=dout, conf = xnews, confo=xnew, pars=c(kappa,lambda,nu), niter = itel, stress=sqrt(stress), spp=spp, ndim=ndim, model=paste("Copstress",optimmethod), call=match.call(), nobj = dim(xnew)[1], type = "copstress", gamma=NA, stress.m=stress, stress.en=stressen, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat)
    out$par <- theta
    out$loss <- "copstress"
    out$OC <- cordillera::cordillera(out$conf,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=FALSE)
    out$OCv <- ifelse(normed,out$OC$normed,out$OC$raw)
    out$copstress <- ovalue
    out$optim <- optimized
    out$stressweight <- stressweight
    out$cordweight <- cordweight
    out$optimethod <- optimmethod
    out$losstype <- out$loss
    out$nobj <- dim(out$conf)[1]
    class(out) <- c("cops","smacofP","smacofB","smacof")
    out
}



#' Fitting a COPS-C Model (COPS Variant 1).
#'
#' Minimizing Copstress to obtain a clustered MDS configuration with given hyperparameters theta.
#'
#' @param delta numeric matrix or dist object of a matrix of proximities
#' @param kappa power transformation for fitted distances
#' @param lambda power transformation for proximities
#' @param nu power transformation for weights
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances if it exists), the second lambda (for the observed proximities if it exist), the third is nu (for the weights if it exists) . If less than three elements are is given as argument, it will be recycled. Defaults to 1 1 1. Will override any kappa, lmabda, nu parameters if they are given and do not match
#' @param type what type of MDS to fit. Currently one of "ratio" or "interval".
#' @param weightmat (optional) a matrix of nonnegative weights; defaults to 1 for all off diagonals
#' @param ndim number of dimensions of the target space
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 0.975
#' @param cordweight weight to be used for the cordillera; defaults to 0.025
#' @param q the norm of the corrdillera; defaults to 1
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to ndim+1
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param dmax The winsorization limit of reachability distances in the OPTICS Cordillera. If supplied, it should be either a numeric value that matches max(rang) or NULL; if NULL it is found as 1.5 times (for kappa >1) or 1 times (for kappa <=1) the maximum reachbility value of the power torgerson model with the same lambda. If dmax and rang are supplied and dmax is not max(rang), a warning is given and rang takes precedence.   
#' @param rang range of the reachabilities to be considered. If missing it is found from the initial configuration by taking 0 as the lower boundary and dmax (see above) as upper boundary. See also \code{\link{cordillera}}     
#' @param optimmethod What optimizer to use? Choose one of NEWUOA, Nelder-Mead, Hooke-Jeeves (hjk), BB::dfsane, NlcOptim::solnl, Rsolnp::solnp, subplex::subplex, SANN (derivative free) and BFGS (numeric derivative). Usually hkj, BFGS, NEWUOA, subplex and solnl work well (depending on the smoothness of copstress). There are also combinations that are good, like hjk/NEWUOA (twostep1), hjk/BFGS (twostep2), BFGS/hjk (twostep6), BFGS/NEWUOA (twostep5), hjk/solnl (twostep3), hjk/subplex (twostep4). Use twostep1.   
#' @param verbose numeric value hat prints information on the fitting process; >2 is very verbose
#' @param normed should the cordillera be normed; defaults to TRUE
#' @param scale Allows to scale the configuration for the OC (the scaled configuration is also returned as $conf). One of none (so no scaling), sd (configuration divided by the maximum standard deviation of the columns), std (standardize all columns !NOTE: This does not preserve the relative distances of the optimal config), proc (procrustes adjustment to the initial fit) and rmsq (configuration divided by the maximum root mean square of the columns). Default is sd.   
#' @param accuracy numerical accuracy, defaults to 1e-8
#' @param itmax maximum number of iterations. Defaults to 100000
#' @param stresstype which stress to use in the copstress. Defaults to stress-1. If anything else is set, explicitly normed stress which is (stress-1)^2. Using stress-1 puts more weight on MDS fit.   
#' @param ... additional arguments to be passed to the optimization procedure
#'
#'@return A list with the components
#'         \itemize{
#'         \item delta: the original transformed dissimilarities
#'         \item obsdiss: the explicitly normed transformed dissimilarities (which are approximated by the fit)
#'         \item confdiss: the fitted distances
#'         \item conf: the configuration to which the scaling of argument scale was applied
#'         \item confo: the unscaled but explicitly normed configuration returned from the fitting procedure. Scaling applied to confo gives conf.
#'         \item par, pars : the theta vector of powers tranformations (kappa,lambda,nu)
#'         \item niter: number of iterations of the optimizer
#'         \item stress: the square root of explicitly normalized stress (calculated for confo)
#'         \item spp: stress per point
#'         \item ndim: number of dimensions
#'         \item model: Fitted model name with optimizer
#'         \item call: the call
#'         \item nobj: the number of objects
#'         \item type, loss, losstype: stresstype
#'         \item stress.m, stress.en: other ways to calculate the stress
#'         \item deltaorig: the original untransformed dissimilarities  
#'         \item copstress: the copstress loss value
#'         \item resmat: the matrix of residuals
#'         \item weightmat: the matrix of untransformed weights 
#'         \item OC: the (normed) OPTICS Cordillera object (calculated for scaled conf)
#'         \item OCv: the (normed) OPTICS Cordillera value alone (calculated for scaled conf)
#'         \item optim: the object returned from the optimization procedure
#'         \item stressweight, cordweight: the weights of the stress and OC respectively (v_1 and v_2)
#'         \item optimmethod: The optimizer
#' }
#' 
#'@examples
#'dis<-as.matrix(smacof::kinshipdelta)
#'
#'#Copstress with equal weight to stress and cordillera 
#'res1<-copstressMin(dis,stressweight=0.5,cordweight=0.5) 
#'res1
#'summary(res1)
#'plot(res1)  #super clustered
#'
#' @import cordillera
#' @importFrom stats dist as.dist optim sd
#' @importFrom minqa newuoa
#' @importFrom dfoptim hjk
#' @importFrom NlcOptim solnl
#' @importFrom Rsolnp solnp
#' @importFrom subplex subplex
#' @importFrom crs snomadr
#' @importFrom cmaes cma_es
#' 
#' 
#' 
#' @keywords clustering multivariate
#' @export
copstressMin <- function (delta, kappa=1, lambda=1, nu=1, theta=c(kappa,lambda,nu), type=c("ratio","interval"), weightmat=1-diag(nrow(delta)),  ndim = 2, init=NULL, stressweight=0.975,cordweight=0.025,q=1,minpts=ndim+1,epsilon=10,dmax=NULL,rang,optimmethod=c("NelderMead","Newuoa","BFGS","SANN","hjk","solnl","solnp","subplex","snomadr","hjk-Newuoa","hjk-BFGS","BFGS-hjk","Newuoa-hjk","cmaes"),verbose=0,scale=c("sd","rmsq","std","proc","none"),normed=TRUE, accuracy = 1e-7, itmax = 100000, stresstype=c("stress-1","stress"),...)
{
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    ## -- Setup for MDS type
    if(missing(type)) type <- "ratio"
    trans <- type
    if (trans=="ratio"){
    trans <- "none"
    }
    #TODO: if we want other optimal scalings as well
    #else if (trans=="ordinal" & ties=="primary"){
    #trans <- "ordinalp"
  #} else if(trans=="ordinal" & ties=="secondary"){
  #  trans <- "ordinals"
  #} else if(trans=="ordinal" & ties=="tertiary"){
  #  trans <- "ordinalt"
  #} else if(trans=="spline"){
  #  trans <- "mspline"
  #}
    if(type=="interval") theta <- c(1,1,1) #We dont allow powers in intervall MDS as of yet
    kappa <- theta[1]
    lambda <- theta[2]
    nu <- theta[3]
    plot <- FALSE

    n <- dim(delta)[1]
    
    if(verbose>0) cat(paste("Minimizing",type,"copstress with kappa=",kappa,"lambda=",lambda,"nu=",nu,".\n"))
    if(missing(optimmethod)) optimmethod <- "Newuoa"
    if(missing(scale)) scale <- "sd"
    if(missing(stresstype)) stresstype <- "stress-1"

    ##-- Prepare for dissimilarity scaling
    r <- kappa/2
    deltaorig <- delta #delta is diss in smacof, deltaorig is delta in smacof 
    delta <- delta^lambda
    weightmato <- weightmat
    weightmat <- weightmat^nu
    weightmat[!is.finite(weightmat)] <- 1 #new
    deltaold <- delta
    disobj <- smacof::transPrep(as.dist(delta), trans = trans, spline.intKnots = 2, spline.degree = 2)#spline.intKnots = spline.intKnots, spline.degree = spline.degree) #FIXME: only works with dist() style object 
    ## Add an intercept to the spline base transformation
    #if (trans == "mspline") disobj$base <- cbind(rep(1, nrow(disobj$base)), disobj$base)
    
    delta <- delta / enorm (delta, weightmat) #normalize to sum=1
    ## --- starting rang if not given
    if(missing(rang))
        #perhaps put this into the optimization function?
       {
        if(is.null(dmax))
           {
           if(verbose>1) cat ("Fitting configuration for rang. \n")    
           #initsol <- cops::powerstressFast(delta,kappa=kappa,lambda=lambda,nu=nu,weightmat=weightmat,ndim=ndim)
           initsol <- torgerson(delta,p=ndim)
           #init0 <- initsol$conf
           init0 <- initsol
           init0 <- init0/enorm(init0)
           if(scale=="std") init0 <- scale(init0) #standardizes config before cordillera
           if(scale=="none") init0 <- init0
           if(scale=="sd") #scales config to sd=1 for most spread dimension before cordillera
             {
                init0 <- init0/max(apply(init0,2,stats::sd))
             }   
             if(scale=="rmsq") #scales config to rmsq=1 for most spread dimension before cordillera
             {
                 testso <- scale(init0,center=FALSE)
                 init0 <- init0/max(attr(testso,"scaled:scale"))
             }
             if(scale=="proc") #scales config by procrusting to init
             {
                 if(!exists("init")) init <- initsol
                 procr <- smacof::Procrustes(init,init0)
                 init0 <- procr$Yhat
             }
           crp <- cordillera::cordillera(init0,q=q,minpts=minpts,epsilon=epsilon,scale=FALSE)$reachplot
           cin <- max(crp)
           dmax <- ifelse(kappa>1,1.5*cin,1.1*cin)
           #if(verbose > 2) cat("rang, dmax was NULL or missing which makes the cordillera a goodness-of-clustering relative to the largest distance of each given configuration. \n")
           }
        rang <- c(0,dmax) #This sets the range to (0,dmax)
        if(verbose>1) cat("dmax is",dmax,". rang is",rang,".\n")
    }
    if(isTRUE(max(rang)!=dmax)) warning("Note: The supplied dmax and rang do not match. I took supplied rang as rang.\n")

    ## --- starting values
    if(is.null(init))
    {
        if(exists("init0")) init <- init0 else init <- cops::powerStressFast(delta,kappa=kappa,lambda=lambda,nu=nu,ndim=ndim)$conf
    }
    xold <- init
    xold <- xold/enorm(xold)
    copsf <- function(x,delta,disobj,r,n,ndim,weightmat,stressweight,cordweight,q,minpts,epsilon,rang,scale,normed,init,...)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=ndim)
             delta <- delta/enorm(delta,weightmat)             
             x <- x/enorm(x)
             dnew <- sqdist (x)
             e <- as.dist(sqrt(dnew)) #I need the dist(x) here for interval
             #e <- dist(x) #I need the dist(x) here for interval
             dhat <- smacof::transform(e, disobj, w = as.dist(weightmat), normq = 0.5)  ## dhat update FIXME: check if it works okay to have as.dist(weightmat) here
             dhatt <- dhat$res #FIXME: I need the structure here to reconstruct the delta; alternatively turn all into vectors? - check how they do it in smacof
             dhatd <- structure(dhatt, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
             #FIXME: labels
             delta <<- as.matrix(dhatd)
             rnew <- sum (weightmat * delta * mkPower (dnew, r))
             nnew <- sum (weightmat * mkPower (dnew,  2*r))
             anew <- rnew / nnew
             stressi <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
             if(stresstype=="stress-1") stressi <- sqrt(stressi)
             if(scale=="none") x <- x 
             if(scale=="std") x <- scale(x) #standardizes config before cordillera
             if(scale=="sd") #scales config to sd=1 for most spread dimension before cordillera
             {
                x <- x/max(apply(x,2,stats::sd))
             }   
             if(scale=="rmsq") #scales config to rmsq=1 for most spread dimension before cordillera
             {
                 testso <- scale(x,center=FALSE)
                 x <- x/max(attr(testso,"scaled:scale"))
             }
             if(scale=="proc") #scales config by procrusting to init
             {
                 procr <- smacof::Procrustes(init,x)
                 x <- procr$Yhat
             }
             corrd <- cordillera::cordillera(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=FALSE,...)
             struc <- corrd$raw
             if(normed) {
                        struc <- corrd$normed
                       }
             ic <- stressweight*stressi - cordweight*struc
             if(verbose>1) cat("copstress =",ic,"mdsloss =",stressi,"OC =",struc,"minpts=",minpts,"kappa =",kappa,"lambda =",lambda,"nu=",nu,"\n")
             ic
           }
     if(verbose>1) cat("Starting Minimization with",optimmethod,":\n")
    if(optimmethod=="Newuoa") {
         optimized <- minqa::newuoa(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose-2),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- optimized$feval
         ovalue <-optimized$fval
     }
    if(optimmethod=="cmaes") {
         xold <- as.vector(xold)
         optimized <- cmaes::cma_es(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmax),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- optimized$counts
         ovalue <-optimized$value
     }
     if(optimmethod=="NelderMead") {
         optimized <- optim(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmax/16,trace=0,reltol=accuracy),...)
         xnew <- optimized$par
         itel <- optimized$counts[[1]]
         ovalue <-optimized$val 
     }
    if(optimmethod=="BFGS") {
         optimized <- optim(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="BFGS",control=list(maxit=itmax,trace=0,reltol=accuracy),...)
         xnew <- optimized$par
         itel <- optimized$counts[[1]]
         ovalue <-optimized$val 
     }
    if(optimmethod=="SANN") {
         optimized <- optim(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="SANN",control=list(maxit=itmax,trace=0,reltol=accuracy),...)
         xnew <- optimized$par
         itel <- optimized$counts[[1]]
         ovalue <-optimized$val 
     }
     if(optimmethod=="hjk") {
         optimized <- dfoptim::hjk(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfeval=itmax,trace=0),...)
         xnew <- optimized$par
         itel <- optimized$feval
         ovalue <-optimized$value 
     }
    if(optimmethod=="hjk-Newuoa") { #twostep1
         optimized1 <- dfoptim::hjk(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfeval=itmax),...)
         xnew <- optimized1$par
         itel <- optimized1$feval
         optimized <- minqa::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmax-itel,rhoend=accuracy,iprint=verbose-2),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- itel+optimized$feval
         ovalue <-optimized$fval
     }
    if(optimmethod=="hjk-BFGS") { #twostep2
         optimized1 <- dfoptim::hjk(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfeval=itmax,trace=0),...)
         xnew <- optimized1$par
         itel <- optimized1$feval
         optimized <- optim(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="BFGS",control=list(maxit=itmax-itel,trace=0,reltol=accuracy),...)
         xnew <- optimized$par
         itel <- itel+optimized$counts[[1]]
         ovalue <-optimized$val
     }
     if(optimmethod=="BFGS-hjk") { #twostep6
         optimized1 <- optim(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="BFGS",control=list(maxit=itmax,trace=0,reltol=accuracy),...)
         xnew <- optimized1$par
         itel <- optimized1$counts[[1]]
         optimized <- dfoptim::hjk(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfeval=itmax-itel,trace=0),...)
         xnew <- optimized$par
         itel <- itel+optimized$feval
         ovalue <-optimized$val
     }
     if(optimmethod=="BFGS-Newuoa") { #twostep5
         optimized1 <- optim(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="BFGS",control=list(maxfeval=itmax,trace=0),...)
         itel <- optimized1$counts[[1]]
         xnew <- optimized1$par
         optimized <- minqa::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmax-itel,rhoend=accuracy,iprint=verbose-2),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- itel+optimized$feval
         ovalue <-optimized$fval
     }
        if(optimmethod=="hjk-solnl") {#twostep3
         optimized1 <- dfoptim::hjk(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfeval=itmax,trace=0),...)
         xnew <- optimized1$par
         itel <- optimized1$feval
         optimized <- NlcOptim::solnl(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),maxnFun=itmax-itel,tolFun=accuracy,...)
         xnew <- optimized$par
         itel <- itel+optimized$counts[[1]]
         ovalue <-optimized$fn
     }
      if(optimmethod=="hjk-subplex") {#twostep4
         optimized1 <- dfoptim::hjk(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfeval=itmax,trace=0),...)
         xnew <- optimized1$par
         itel <- optimized1$feval
         optimized <- subplex::subplex(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmax-itel+1,reltol=accuracy),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- itel+optimized$count
         ovalue <-optimized$value
     }
  #   if(optimmethod=="dfsane") {
  #       optimized <- BB::dfsane(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmax,trace=0,tol=accuracy),...)
   #      xnew <- optimized$par
    #     itel <- optimized$feval
    #     ovalue <-optimized$residual 
    # }
     if(optimmethod=="solnl") {
         optimized <- NlcOptim::solnl(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),maxnFun=itmax,tolFun=accuracy,...)
         xnew <- optimized$par
         itel <- optimized$counts[[1]]
         ovalue <-optimized$fn 
     }
     ## if(optimmethod=="isres") {
     ##     optimized <- nloptr::isres(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),maxeval=itmax,trace=verbose-2,xtol_rel=accuracy,...)
     ##     xnew <- optimized$par
     ##     itel <- optimized$iter
     ##     ovalue <-optimized$value 
     ## }
    if(optimmethod=="solnp") {
         optimized <- Rsolnp::solnp(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(outer.iter=itmax,trace=0,tol=accuracy),...)
         xnew <- matrix(optimized$pars,ncol=ndim)
         itel <- optimized$nfuneval
         ovalue <-tail(optimized$values,1) 
     }
     if(optimmethod=="subplex") {
         optimized <- subplex::subplex(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmax,reltol=accuracy),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- optimized$count
         ovalue <-optimized$value 
     }
    if(optimmethod=="snomadr") {
         optimized <- crs::snomadr(function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),n=dim(xold)[1],x0=xold,print.output=isTRUE(verbose-2>0),...)
         xnew <- optimized$solution
         itel <- optimized$iterations
         ovalue <-optimized$objective 
     }
     xnew <- xnew/enorm(xnew)
     dnew <- sqdist (xnew)
     rnew <- sum (weightmat * delta * mkPower (dnew, r))
     nnew <- sum (weightmat * mkPower (dnew,  2*r))
     anew <- rnew / nnew
     stress <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
     if(stresstype=="stress-1") stress <- sqrt(stress)
     attr(xnew,"dimnames")[[1]] <- rownames(delta)
     attr(xnew,"dimnames")[[2]] <- paste("D",1:ndim,sep="")
     #doutm <- (2*sqrt(sqdist(xnew)))^kappa  #fitted powered euclidean distance but times two
     doutm <- as.matrix(dist(xnew)^kappa)
     #deltam <- delta
     #deltaorigm <- deltaorig
     #deltaoldm <- deltaold
     delta <- stats::as.dist(delta)
     deltaorig <- stats::as.dist(deltaorig)
     deltaold <- stats::as.dist(deltaold)
     #doute <- doutm/enorm(doutm)
     #doute <- stats::as.dist(doute)
     dout <- doute <- stats::as.dist(doutm)
     resmat <- as.matrix((delta - doute)^2)
     spp <- colMeans(resmat)
     #weightmatm <-weightmat
     weightmat <- stats::as.dist(weightmat)
     stressen <- sum(weightmat*(delta-doute)^2)#raw stress on the normalized proximities and normalized distances
     if(scale=="std") xnews <- scale(xnew) #standardizes config before cordillera
     if(scale=="sd") #scales config to sd=1 for most spread dimension before cordillera
             {
                xnews <- xnew/max(apply(xnew,2,stats::sd))
             }   
             if(scale=="rmsq") #scales config to rmsq=1 for most spread dimension before cordillera
             {
                 testso <- scale(xnew,center=FALSE)
                 xnews <- xnew/max(attr(testso,"scaled:scale"))
             }
             if(scale=="proc") #scales config by procrusting to init
             {
                 procr <- smacof::Procrustes(init,xnew)
                 xnews <- procr$Yhat
             }
             if(scale=="none") xnews <- xnew #no standardisation
      if(verbose>0) cat("*** stress (both normalized - for COPS/STOPS):",stress,"; stress 1 (both normalized - default reported):",sqrt(stress),"; stress manual (for debug only):",stressen,"; from optimization: ",ovalue,"\n")   
    out <- list(delta=deltaold, obsdiss=delta, confdiss=dout, conf = xnews, confo=xnew, pars=c(kappa,lambda,nu), niter = itel, stress=sqrt(stress), spp=spp, ndim=ndim, model=paste(type,"copstress",optimmethod), call=match.call(), nobj = n, type = type, gamma=NA, stress.m=stress, stress.en=stressen, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat)
    out$par <- theta
    out$loss <- "copstress"
    out$OC <- cordillera::cordillera(out$conf,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=FALSE)
    out$OCv <- ifelse(normed,out$OC$normed,out$OC$raw)
    out$copstress <- ovalue
    out$optim <- optimized
    out$stressweight <- stressweight
    out$cordweight <- cordweight
    out$optimethod <- optimmethod
    out$losstype <- out$loss
    #out$nobj <- dim(out$conf)[1]
    class(out) <- c("cops","smacofP","smacofB","smacof")
    out
}


#' High Level COPS Function
#'
#' Minimizing copstress for a clustered MDS configuration. Allows to choose COPS-C (finding a configuration from copstress with cordillera penalty) and profile COPS (finding hyperparameters for MDS models with power transformations). It is wrapper for copstressMin and pcops.
#'
#'@param dis a dissimilarity matrix or a dist object
#'@param variant a character string specifying which variant of COPS to fit. Allowed is any of the following "1","2","Variant1","Variant2","v1","v2","COPS-C","P-COPS","configuration-c","profile","copstress-c","p-copstress". Defaults to "COPS-C".
#'@param ... arguments to be passed to  \code{\link{copstressMin}} (for Variant 1) or \code{\link{pcops}} (for Variant 2).
#'
#'@return For COPS-C Variant 1 see \code{\link{copstressMin}}, for P-COPS Variant 2 see \code{\link{pcops}}
#' 
#'@examples
#'\donttest{
#'dis<-as.matrix(smacof::kinshipdelta)
#'
#'#COPS-C with equal weight to stress and cordillera 
#'res1<-cops(dis,variant="COPS-C",stressweight=0.5,cordweight=0.5,minpts=2) 
#'res1
#'summary(res1)
#'plot(res1)
#'
#'#ratio mds (i.e. COPS-C with stressweight=1 and cordweight=0)
#'res2<-cops(dis,variant="COPS-C",stressweight=1,cordweight=0,minpts=2) 
#'res2
#'summary(res2)
#'plot(res2)
#' 
#'#procrustes adjustment
#'plot(Procrustes(res1$conf,res2$conf))
#'
#'par(mfrow=c(1,2))
#'plot(res1,"reachplot")
#'plot(res2,"reachplot") 
#'par(mfrow=c(1,1))
#'
#' 
#'
#'#s-stress type copstress (i.e. kappa=2, lambda=2)
#'res3<-cops(dis,variant="COPS-C",kappa=2,lambda=2,stressweight=0.5,cordweight=0.5) 
#'res3
#'summary(res3)
#'plot(res3)
#'
#' 
#'#power-stress type profile copstress
#'# search for optimal kappa and lambda between kappa=0.5,lambda=0.5 and kappa=2,lambda=5
#'# nu is fixed on -1
#'ws<-1/dis
#'diag(ws)<-1 
#'res5<-cops(dis,variant="P-COPS",loss="powerstress",theta=c(1.4,3,-1),lower=c(1,0.5,-1),upper=c(3,5,-1),weightmat=ws,stressweight=0.9,cordweight=0.1) 
#'res5
#'summary(res5)
#'plot(res5)
#'}
#'
#' @import cordillera
#' @importFrom stats dist as.dist optim
#' @importFrom minqa newuoa
#' 
#' 
#'@keywords clustering multivariate
#'@export
cops <- function(dis, variant=c("1","2","Variant1","Variant2","v1","v2","COPS-C","P-COPS","configuration-c","profile","copstress-c","p-copstress","COPS-P","copstress-p","cops-c","p-cops"),...)
{
    #variant=c("0","1","2","Variant0","Variant1","Variant2","v0","v1","v2","COPS-0","COPS-C","P-COPS","configuration-0","configuration-c","profile","copstress-0","copstress-c","p-copstress","COPS-P","copstress-p","cops-c","p-cops")
                 if(missing(variant)) variant <- "COPS-C"
                 if(variant%in%c("1","Variant1","v1","configuration-c","COPS-C","copstress-c","cops-c")) out <- copstressMin(dis,...)
#                 if(variant%in%c("0","Variant0","v0","configuration-0","COPS-0","copstress-c")) stop("Not yet implemented") out <- shrinkCopstress0(dis,...)
                 if(variant%in%c("2","Variant2","v2","profile","p-copstress","P-COPS","COPS-P","p-cops","copstress-p")) out <- pcops(dis,...)
                 return(out)
                 }


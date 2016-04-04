#' COPS versions of smacofSym models
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this is either a scalar of the lambda transformation for the observed proximities, or a vector where the first is the kappa argument for the fitted distances (here internally fixed to 1) and the second the lambda argument and the third the nu argument (here internally fixed to 1). Defaults to 1 1 1 
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
#' @param scale should the configuration be scaled to mean=0 and sd=1? Defaults to TRUE
#' @param stresstype which stress to report. Only takes smacofs default stress currrently.
#' @param ... additional arguments to be passed to the fitting procedure
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress
#'         \item{stress.m:} default normalized stress
#'         \item{coploss:} the weighted loss value
#'         \item{OC:} the Optics cordillera value
#'         \item{parameters:} the parameters used for fitting (kappa, lambda)
#'         \item{fit:} the returned object of the fitting procedure (which has all smacofB elements and some more
#'         \item{cordillera:} the cordillera object
#' }
#'
#'@importFrom stats dist as.dist
#'@import smacof
#' 
#'@keywords multivariate
#'@export
cop_smacofSym <- function(dis,theta=c(1,1,1),ndim=2,weightmat=NULL,init=NULL,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,normed=TRUE,scale=TRUE,stresstype="default") {
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
  delts <- as.matrix(fit$delta) #That was my choice to not use the normalized deltas but try it ion the original; that is scale and unit free as Buja said
 # if(inherits(fit,"smacofSP")) delts <- as.matrix(fit$delta)[-1,-1]
  fit$stress.r <- sum(weightmat*(delts-fitdis)^2)
  fit$stress.m <- fit$stress.r/sum(weightmat*delts^2)
  fit$pars <- c(fit$kappa,fit$lambda,fit$nu)
  fit$deltaorig <- fit$delta^(1/fit$lambda)   
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed)
  out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress.m, coploss=copobj$coploss, OC=copobj$OC, parameters=copobj$parameters, fit=fit,cordillera=copobj) #target functions
  out
}

#' COPS versions of elastic scaling models (via smacofSym)
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
#' @param scale should the configuration be scaled to mean=0 and sd=1? Defaults to TRUE
#' @param stresstype which stress to report. Only takes smacofs default stress currrently.
#' @param ... additional arguments to be passed to the fitting procedure
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress
#'         \item{stress.m:} default normalized stress
#'         \item{coploss:} the weighted loss value
#'         \item{OC:} the Optics cordillera value
#'         \item{parameters:} the parameters used for fitting (kappa, lambda)
#'         \item{fit:} the returned object of the fitting procedure
#'         \item{cordillera:} the cordillera object
#' }
#'
#'
#'@importFrom stats dist as.dist
#'@import smacof
#'@keywords multivariate
#'@export
cop_elastic <- function(dis,theta=c(1,1,-2),ndim=2,weightmat=1,init=NULL,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,normed=TRUE,scale=TRUE,stresstype="default") {
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed)
  out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress.m, coploss=copobj$coploss, OC=copobj$OC, parameters=copobj$parameters, fit=fit,cordillera=copobj) #target functions
  out
}


#' COPS versions of smacofSphere models
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
#' @param scale should the configuration be scaled to mean=0 and sd=1? Defaults to TRUE
#' @param stresstype which stress to report. Only takes smacofs default stress currrently.
#' @param ... additional arguments to be passed to the fitting procedure
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress
#'         \item{stress.m:} default normalized stress
#'         \item{coploss:} the weighted loss value
#'         \item{OC:} the Optics cordillera value
#'         \item{parameters:} the parameters used for fitting (kappa, lambda)
#'         \item{fit:} the returned object of the fitting procedure
#'         \item{cordillera:} the cordillera object
#' }
#'
#'@import smacof 
#'@importFrom stats dist as.dist
#'@keywords multivariate
#'@export
cop_smacofSphere <- function(dis,theta=c(1,1),ndim=2,weightmat=NULL,init=NULL,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,normed=TRUE,scale=TRUE,stresstype="default") {
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed)
  out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress.m, coploss=copobj$coploss, OC=copobj$OC, parameters=copobj$parameters, fit=fit,cordillera=copobj) #target functions
  out
}


#' COPS version of sammon mapping
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
#' @param scale should the configuration be scaled to mean=0 and sd=1? Defaults to TRUE
#' @param stresstype which stress to report. Only takes smacofs default stress currrently.
#' @param ... additional arguments to be passed to the fitting procedure
#' 
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item coploss: the weighted loss value
#'         \item OC: the Optics cordillera value
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#'         \item cordillera: the cordillera object
#' }
#'
#' @importFrom stats dist as.dist
#' 
#' @keywords multivariate
#'
#' 
#' @export
cop_sammon <- function(dis,theta=c(1,1,-1),ndim=2,init=NULL,weightmat=NULL,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale=TRUE,normed=TRUE,stresstype="default") {
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(length(theta)==1L) lambda <- theta
  if(length(theta)==2L) lambda <- theta[2]
  if(length(theta)==3L) lambda <- theta[2]
  nu <- -1
 # if(is.null(init)) init <- stops::cmdscale(dis^lambda,k=ndim)$points
  fit <- stops::sammon(dis^lambda,k=ndim,y=init,trace=isTRUE(verbose>1),...)
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed)
  list(stress=fit$stress, stress.m=fit$stress.m, coploss=copobj$coploss, OC=copobj$OC, parameters=copobj$parameters,  fit=fit,cordillera=copobj) #target functions
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
#' @param scale should the configuration be scaled to mean=0 and sd=1? Defaults to TRUE
#' @param stresstype which stress to report. Only takes smacofs default stress currrently.
#' @param ... additional arguments to be passed to the fitting procedure
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress
#'         \item{stress.m:} default normalized stress
#'         \item{coploss:} the weighted loss value
#'         \item{OC:} the Optics cordillera value
#'         \item{parameters:} the parameters used for fitting (kappa, lambda)
#'         \item{fit:} the returned object of the fitting procedure
#'         \item{cordillera:} the cordillera object
#' }
#'
#'
#' @importFrom stats dist as.dist
#' @import smacof
#'@keywords multivariate
#'@export
cop_sammon2 <- function(dis,theta=c(1,1,-1),ndim=2,weightmat=NULL,init=NULL,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,normed=TRUE,scale=TRUE,stresstype="default") {
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed)
  out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress.m, coploss=copobj$coploss, OC=copobj$OC, parameters=copobj$parameters, fit=fit,cordillera=copobj) #target functions
  out
}



#' COPS version of strain
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
#' @param scale should the configuration be scaled to mean=0 and sd=1? Defaults to TRUE
#' @param stresstype which stress to report. Only takes cmdscales default stress currrently.
#' @param ... additional arguments to be passed to the fitting procedure
#'
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item coploss: the weighted loss value
#'         \item OC: the Optics cordillera value
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#'         \item cordillera: the cordillera object
#' }
#' 
#' @importFrom stats dist as.dist
#' @keywords multivariate
#' @export
cop_cmdscale <- function(dis,theta=c(1,1,1),weightmat=NULL,ndim=2,init=NULL,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale=TRUE,normed=TRUE,stresstype="default") {
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)==1L) lambda <- theta
  if(length(theta)==2L) lambda <- theta[2]
  if(length(theta)==3L) lambda <- theta[2]
  fit <- stops::cmdscale(dis^lambda,k=ndim,eig=TRUE,...) 
  fit$lambda <- lambda
  fit$kappa <- 1
  fit$nu <- 1
  dis <- stats::as.dist(dis)
  fitdis <- stats::dist(fit$points)
  fit$stress.r <- sum((dis^lambda-fitdis)^2)
  fit$stress.n <- fit$stress.r/sum(dis^(2*lambda))
  fit$stress.m <- sqrt(fit$stress.n)
  fit$conf <- fit$points
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed)
  list(stress=fit$GOF,stress.m=fit$stress.m, coploss=copobj$coploss, OC=copobj$OC, parameters=copobj$parameters, fit=fit,cordillera=copobj) #target functions
}

#' COPS version of rstress
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
#' @param scale should the configuration be scaled to mean=0 and sd=1? Defaults to TRUE
#' @param stresstype which stress to report? Defaults to explicitly normed stress
#' @param ... additional arguments to be passed to the fitting procedure
#' 
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item coploss: the weighted loss value
#'         \item OC: the Optics cordillera value
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#'         \item cordillera: the cordillera object
#' }
#' @keywords multivariate
#' @export
cop_rstress <- function(dis,theta=c(1,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale=TRUE,normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, coploss=copobj$coploss, OC=copobj$OC, parameters=copobj$parameters, fit=fit,cordillera=copobj)
  out
}


#' COPS version of sstress
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
#' @param scale should the configuration be scaled to mean=0 and sd=1? Defaults to TRUE
#' @param stresstype which stress to report? Defaults to explicitly normed stress
#' @param ... additional arguments to be passed to the fitting procedure
#' 
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item coploss: the weighted loss value
#'         \item OC: the Optics cordillera value
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#'         \item cordillera: the cordillera object
#' }
#' @keywords multivariate
#' @export
cop_sstress <- function(dis,theta=c(2,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale=TRUE,normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, coploss=copobj$coploss, OC=copobj$OC, parameters=copobj$parameters, fit=fit,cordillera=copobj)
  out
}


#' COPS version of powermds
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
#' @param scale should the configuration be scaled to mean=0 and sd=1? Defaults to TRUE
#' @param stresstype which stress to report? Defaults to whatever whim is my default (currently explicitly normed stress)
#' @param ... additional arguments to be passed to the fitting procedure
#'
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item coploss: the weighted loss value
#'         \item OC: the Optics cordillera value
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#'         \item cordillera: the cordillera object
#' }
#' @keywords multivariate
#' @export
cop_powermds <- function(dis,theta=c(1,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale=TRUE,normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, coploss=copobj$coploss, OC=copobj$OC, parameters=copobj$parameters, fit=fit, cordillera=copobj)
  out 
}

#' COPS version of sammon with powers
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
#' @param scale should the configuration be scaled to mean=0 and sd=1? Defaults to TRUE
#' @param stresstype which stress to report? Defaults to explicitly normed stress
#' @param ... additional arguments to be passed to the fitting procedure
#'
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item coploss: the weighted loss value
#'         \item OC: the Optics cordillera value
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#'         \item cordillera: the cordillera object
#' }
#' @keywords multivariate
#' @export
cop_powersammon <- function(dis,theta=c(1,1,-1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale=TRUE,normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, coploss=copobj$coploss, OC=copobj$OC, parameters=copobj$parameters, fit=fit, cordillera=copobj)
  out 
}

#' COPS version of elastic scaling with powers
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
#' @param scale should the configuration be scaled to mean=0 and sd=1? Defaults to TRUE
#' @param stresstype which stress to report? Defaults to explicitly normed stress
#' @param ... additional arguments to be passed to the fitting procedure
#'
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item coploss: the weighted loss value
#'         \item OC: the Optics cordillera value
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#'         \item cordillera: the cordillera object
#' }
#' @keywords multivariate
#' @export
cop_powerelastic <- function(dis,theta=c(1,1,-2),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale=TRUE,normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, coploss=copobj$coploss, OC=copobj$OC, parameters=copobj$parameters, fit=fit, cordillera=copobj)
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
#' @param scale should the configuration be scaled to mean=0 and sd=1? Defaults to TRUE
#' @param stresstype which stress to report? Defaults to explicitly normed stress
#' @param ... additional arguments to be passed to the fitting procedure
#'
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item coploss: the weighted loss value
#'         \item OC: the Optics cordillera value
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#'         \item cordillera: the cordillera object
#' }
#' @keywords multivariate
#' @export
cop_powerstress <- function(dis,theta=c(1,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,cordweight=0.5,q=1,minpts=ndim+1,epsilon=10,rang=NULL,verbose=0,scale=TRUE,normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),scale=scale,normed=normed)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, coploss=copobj$coploss, OC=copobj$OC, parameters=copobj$parameters, fit=fit, cordillera=copobj)
  out 
}

#' Calculates coploss for given MDS object 
#'
#' @param obj MDS object (supported are sammon, cmdscale, smacof, rstress, powermds)
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param cordweight weight to be used for the cordillera; defaults to 0.5
#' @param q the norm of the corrdillera; defaults to 1
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to 2
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the distances (min distance minus max distance). If NULL (default) the cordillera will be normed to each configuration's maximum distance, so an absolute value of goodness-of-clusteredness. 
#' @param verbose numeric value hat prints information on the fitting process; >2 is very verbose (coploss level), >3 is extremely (up to MDS optimization level)
#' @param normed should the cordillera be normed; defaults to TRUE
#' @param scale should the configuration be scaled to mean=0 and sd=1? Defaults to TRUE
#' @param ... additional arguments to be passed to the cordillera function
#'
#' @return A list with the components
#' \itemize{
#'        \item coploss the weighted loss value
#'        \item OC the Optics cordillera value
#'        \item parameters the parameters used for fitting (kappa, lambda)
#'        \item cordillera the cordillera object
#' }
#' @keywords multivariate
#' @export
coploss <- function(obj,stressweight=1,cordweight=0.5,q=1,normed=TRUE,minpts=2,epsilon=10,rang=NULL,verbose=0,scale=TRUE,...)
    {
        stressi <- obj$stress.m
        kappa <- obj$kappa
        lambda <- obj$lambda
        nu <- obj$nu
        confs <- obj$conf 
        corrd <- stops::cordillera(confs,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,...)
        struc <- corrd$raw
        maxstruc <- corrd$normi
        if(normed) {
                   struc <- corrd$normed
                   maxstruc <- 1
                   }
        ic <- stressweight*stressi - cordweight*struc
        if(verbose>0) cat("coploss =",ic,"mdsloss =",stressi,"OC =",struc,"kappa =",kappa,"lambda =",lambda,"nu=",nu,"\n")
        out <- list(coploss=ic,OC=struc,parameters=c(kappa=kappa,lambda=lambda,nu=nu),cordillera=corrd)
        out
     }

#' Profile COPS Function (aka COPS Variant 2)
#'
#' Metaparameter selection for MDS models baseed on the Profile COPS approach (COPS Variant 2). It uses coploss for hyperparameter selection. It is a special case of a STOPS model.  
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
#' @param cordweight weight to be used for the cordillera; if missing gets estimated from the initial configuration so that coploss = 0 for theta=c(1,1) 
#' @param q the norm of the corrdillera; defaults to 1
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to ndim+1
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the minimum reachabilities to be considered. If missing it is found from the initial configuration by taking 1.5 times the maximal minimum reachability of the model with theta=c(1,1). If NULL it will be normed to each configuration's minimum and maximum distance, so an absolute value of goodness-of-clusteredness. Note that the latter is not necessarily desirable when comparing configurations for their relative clusteredness. See also \code{\link{cordillera}}     
#' @param optimmethod What general purpose optimizer to use? Defaults to our adaptive LJ version (ALJ). Also allows particle swarm optimization with s particles ("pso") and simulated annealing ("SANN"). We recommend not using the latter with the rstress, sstress and the power stress models. 
#' @param lower The lower contraints of the search region
#' @param upper The upper contraints of the search region 
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param normed should the cordillera be normed; defaults to TRUE
#' @param scale should the configuration be scaled to mean=0 and sd=1? Defaults to TRUE
#'@param s number of particles if pso is used
#'@param stresstype what stress to be used for comparisons between solutions 
#'@param ... additional arguments to be passed to the optimization procedure
#'
#'@return A list with the components
#'         \itemize{
#'         \item coploss: the weighted loss value
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
#'dis<-as.matrix(smacof::kinshipdelta)
#'set.seed(210485)
#'res1<-pcops(dis,loss="strain",lower=0.1,upper=5,minpts=2) #optimum around lambda=0.16
#'res1
#'summary(res1)
#'plot(res1)
#'#only use cordillera for finding the parameters; 
#'res2<-pcops(dis,loss="strain",stressweight=0,lower=0.1,upper=5,minpts=2) #optimum around lambda=0.156
#'res2
#'summary(res2)
#'plot(res2)
#' 
#'#cordillera value of res1 and res 2 very close but res 2 is a bit more clustered
#'#the reason is the distance between sister and son
#' 
#'#procrustes adjusted
#'resadj<-conf_adjust(res2$fit$points,res1$fit$points)
#'plot(resadj$ref.conf) #res 2
#'plot(resadj$other.conf) #res 1
#' 
#'par(mfrow=c(1,2))
#'plot(res1,"reachplot")
#'plot(res2,"reachplot")
#' par(mfrow=c(1,1))
#'
#'\donttest{
#'# From De Leuuw et al (2016) example 7.2.
#'#They look at different rstress versions and compare how clustered the configuration is
#'#where stress is minimal and that stress is a monotonically increasing function of r;
#' dats <- c(5.63,5.27, 6.72,4.60, 5.64, 5.46,4.80, 6.22, 4.97, 3.20,7.54 ,5.12, 8.13, 7.84 ,7.80, 6.73 ,4.59 ,7.55, 6.73, 7.08, 4.08, 7.18 ,7.22 ,6.90 ,7.28 ,6.96 ,6.34 ,6.88, 6.17, 5.47, 4.67, 6.13, 6.04 ,7.42, 6.36, 7.36)
#'num_cols <- (1 + sqrt(1 + 8*length(dats)))/2 - 1
#'mat <- matrix(0, num_cols, num_cols)
#'mat[row(mat) <= col(mat)] <- dats
#'mat <- t(mat)
#'mat <- rbind(0, mat)
#'mat <- cbind(mat, 0)
#'colnames(mat) <- rownames(mat) <- c(" KVP", "PvdA" , "VVD" , "ARP" , "CHU" , "CPN" , "PSP" ,  "BP", "D66")
#'dobj <- as.dist(mat)
#'dobj
#'#We can do this in one go by setting cordweight to 0 and find that stress is minimal (0.0033) around r~=0.17 (kappa~=0.34)
#'#and that stress appears thus not monotonically increasing in r
#' set.seed(210485)
#' m1 <- pcops(dobj,loss="rstress",lower=c(0.05,1,1),upper=c(5,1,1),verbose=3,cordweight=0,stressweight=1)
#' m1
#'# They observe increasing clustering for larger r which we can again do systematically:
#'# When only clusteredness is of interest, we use cordweight=1 stressweight=0 and try clusters of at least k=2 and k=3 observations
#' set.seed(210485)
#' m2 <- pcops(dobj,loss="rstress",minpts=2,lower=c(0.05,1,1),upper=c(5,1,1),verbose=3,cordweight=1,stressweight=0) 
#' m3 <- pcops(dobj,loss="rstress",minpts=3,lower=c(0.05,1,1),upper=c(5,1,1),verbose=3,cordweight=1,stressweight=0)
#' m2   #r~=1.24
#' m3   #r~=1.39
#'
#'# It is generally better to trade off clusteredness and fit
#' set.seed(210485)
#' m2t <- pcops(dobj,loss="rstress",minpts=2,theta=c(m1$par[1],1,1),lower=c(0.05,1,1),upper=c(5,1,1),verbose=3,cordweight=1/3,stressweight=2/3)
#' m3t <- pcops(dobj,loss="rstress",minpts=3,theta=c(m1$par[1],1,1),lower=c(0.05,1,1),upper=c(5,1,1),verbose=3,cordweight=1/3,stressweight=2/3)
#' m2t #r~=0.08
#' m4t #r~=1.39
#'}
#'
#' @importFrom stats dist as.dist optim
#' @importFrom pso psoptim
#' 
#'@keywords clustering multivariate
#'@export
pcops <- function(dis,loss=c("stress","smacofSym","smacofSphere","strain","sammon","rstress","powermds","sstress","elastic","powersammon","powerelastic","powerstress","sammon2","powerstrain"),weightmat=NULL,ndim=2,init=NULL,theta=c(1,1,1),stressweight=1,cordweight,q=1,minpts=ndim+1,epsilon=10,rang,optimmethod=c("ALJ","pso","SANN"),lower=c(1,1,0.5),upper=c(5,5,2),verbose=0,scale=TRUE,normed=TRUE,s=4,stresstype="default",...)
    {
      if(inherits(dis,"dist")) dis <- as.matrix(dis)
      if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1]) 
      if(missing(loss)) loss <- "strain"
      .confin <- init #initialize a configuration
      psfunc <- switch(loss,"strain"=cop_cmdscale,"powerstrain"=cop_cmdscale,"elastic"=cop_elastic,"sstress"=cop_sstress,"stress"=cop_smacofSym,"smacofSym"= cop_smacofSym,"smacofSphere"=cop_smacofSphere,"rstress"=cop_rstress,"powermds"=cop_powermds,"powerstress"=cop_powerstress,"sammon"=cop_sammon,"sammon2"=cop_sammon2,"powersammon"=cop_powersammon,"powerelastic"=cop_powerelastic) #choose the stress to minimize
      if(missing(optimmethod)) optimmethod <- "ALJ"
      if(missing(rang)) 
          {
           if(verbose>1) cat ("Fitting configuration for rang. \n")    
           initsol <- do.call(psfunc,list(dis=dis,theta=c(1,1,1),init=.confin,weightmat=weightmat,ndim=ndim,rang=c(0,1),q=q,minpts=minpts,epsilon=epsilon,verbose=verbose-2,scale=scale,normed=normed,stresstype=stresstype))
           init0 <- initsol$fit$conf
           if(isTRUE(scale)) init0 <- scale(init0)
           crp <- stops::cordillera(init0,q=q,minpts=minpts,epsilon=epsilon,scale=scale)$reachplot
           cin <- max(crp)
           rang <- c(0,1.5*cin) #approximate upper bound by 1.5 times the max distance in the initial config
                 #alternatives: use an adjusted boxplot idea so e.g., rang<-c(quantile(crp,0.25)-exp(-4*robustbase::mc(crp))*1.5*IQR(crp),quantile(crp,0.75)+exp(4*robustbase::mc(crp))*1.5*IQR(crp)
                 #alternatives: use an adjusted boxplot idea so e.g., c(min(crp)-exp(-4*robustbase::mc(crp))*1.5,max(crp)+exp(4*robustbase::mc(crp))*1.5) 
            if(verbose>1) cat("dmax is",max(rang),". rang is",rang,"\n")
           #.confin <- initsol$fit$conf
           }
      if(is.null(rang) && verbose > 1) cat("rang=NULL which makes the cordillera a goodness-of-clustering relative to the largest distance of each given configuration \n") 
      if(missing(cordweight))
               {
                 if(verbose>1) cat ("Fitting configuration for cordweight. \n")     
                 initsol <- do.call(psfunc,list(dis=dis,theta=c(1,1,1),init=.confin,weightmat=weightmat,ndim=ndim,rang=rang,q=q,minpts=minpts,epsilon=epsilon,verbose=verbose-2,scale=scale,normed=normed,stresstype=stresstype))  
                 initcorrd <- stops::cordillera(initsol$fit$conf,q=q,epsilon=epsilon,minpts=minpts,rang=rang,scale=scale)$normed 
                 if(identical(normed,FALSE)) initcorrd <- stops::cordillera(initsol$fit$conf,q=q,epsilon=epsilon,minpts=minpts,rang=rang,scale=scale)$raw
                cordweight <- initsol$stress.m/initcorrd  
                #cat("stress.m=",initsol$stress.m,"cord=",initcorrd,"cweight=",cordweight,"\n") 
                if(verbose>1) cat("Weights are stressweight=",stressweight,"cordweight=",cordweight,"\n")
             }
      if(verbose>1) cat("Starting Optimization \n ")
      if(optimmethod=="SANN") {
          opt<- stats::optim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,weightmat=weightmat,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,stresstype=stresstype))$coploss,method="SANN",...)
      }
      if(optimmethod=="pso") {
        addargs <- list(...)
        control <- list(trace=verbose-2,s=s,addargs)
        opt<- pso::psoptim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,weightmat=weightmat,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,stresstype=stresstype))$coploss,lower=lower,upper=upper,control=control)
       }
      if(optimmethod=="ALJ") {
      opt<- stops::ljoptim(theta, function(theta) do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=theta,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,stresstype=stresstype))$coploss,lower=lower,upper=upper,verbose=verbose-2,...)
      }
    thetaopt <- opt$par 
    #refit the optimal version (TODO probably unnecessary if the other functions are properly reimplemented)
    out <- do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=thetaopt,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-2,scale=scale,normed=normed,stresstype=stresstype))
    out$OC <- stops::cordillera(out$fit$conf,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale)
    out$coploss <- opt$value
    out$optim <- opt
    out$stressweight <- stressweight
    out$cordweight <- cordweight
    out$call <- match.call()
    out$optimethod <- optimmethod
    out$losstype <- loss
    out$nobj <- dim(out$fit$conf)[1]
    if(verbose>1) cat("Found minimum after",opt$counts["function"]," iterations at",round(opt$par,4),"with coploss=",round(out$coploss,4),"and default scaling loss=",round(out$stress.m,4),"and OC=", round(out$OC$normed,4),". Thanks for your patience. \n")
    class(out) <- c("pcops","cops","stops")
    out
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
    cat("Cluster optimized loss (coploss): ", x$coploss, "\n")
    cat("Stress weight:",x$stressweight," OPTICS Cordillera weight:",x$cordweight,"\n")
    cat("Number of iterations of",x$optimethod,"optimization:", x$niter, "\n")
    cat("\n")
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
    cat("OPTICS cordillera: Raw", x$OC$raw,"Normed", x$OC$normed,"\n")
    cat("Cluster optimized loss (coploss): ", x$coploss, "\n")
    cat("MDS loss weight:",x$stressweight," OPTICS cordillera weight:",x$cordweight,"\n")
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
     #ok, old code here has side effects: it changes the smacof object in the cops object; not sure we should do that  
       if(missing(main)) main <- paste("Transformation Plot") 
       plot.smacofP(x,plot.type="transplot",asp=asp,...)
     # invisible(tmp) #I give the changed smacof object back
     }
     else {      
       plot.smacofP(x,plot.type=plot.type,main=main,asp=asp,...)
   }
 }


#' Fitting a COPS Model. (old version to be discontinued)
#'
#' Minimizing Coploss for a clustered Power Stress MDS configuration with given hyperparameters theta.
#'
#' @param delta numeric matrix or dist object of a matrix of proximities
#' @param kappa power transformation for fitted distances
#' @param lambda power transformation for proximities
#' @param nu power transformation for weights
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances if it exists), the second lambda (for the observed proximities if it exist), the third is nu (for the weights if it exists) . If less than three elements are is given as argument, it will be recycled. Defaults to 1 1 1. Will override any kappa, lmabda, nu parameters if they are given and do not match
#' @param weightmat (optional) a matrix of nonnegative weights; defaults to 1 for all off diagonals
#' @param ndim number of dimensions of the target space
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param cordweight weight to be used for the cordillera; if missing gets estimated from the initial configuration as  
#' @param q the norm of the corrdillera; defaults to 1
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to ndim+1
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the minimum reachabilities to be considered. If missing it is found from the initial configuration by taking 1.5 times the maximal minimum reachability of the initial fit. If NULL it will be normed to each configuration's minimum and maximum distance, so an absolute value of goodness-of-clusteredness. Note that the latter is not necessarily desirable when comparing configurations for their relative clusteredness. See also \code{\link{cordillera}}     
#' @param optimmethod What optimizer to use? Defaults to NEWUOA, Nelder-Mead is also supported.
#' @param verbose numeric value hat prints information on the fitting process; >2 is very verbose
#' @param normed should the cordillera be normed; defaults to TRUE
#' @param scale should the configuration be scaled to mean=0 and sd=1 for calculating the cordillera? Defaults to TRUE
#' @param accuracy numerical accuracy, defulats to 1e-12
#' @param itmax maximum number of iterations. Defaults to 100000
#' @param ... additional arguments to be passed to the optimization procedure
#'
#'@return A list with the components
#'         \itemize{
#'         \item coploss: the weighted loss value
#'         \item OC: the Optics cordillera
#'         \item optim: the object returned from the optimization procedure
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#'         \item cordillera: the cordillera object
#' }
#' 
#'
#'@importFrom stats dist as.dist optim
#'@importFrom minqa newuoa
#' 
#' 
#'@keywords clustering multivariate
#'@export
coplossMinOLD <- function (delta, kappa=1, lambda=1, nu=1, theta=c(kappa,lambda,nu),weightmat=1-diag(nrow(delta)),  ndim = 2, init=NULL, stressweight=1,cordweight,q=1,minpts=ndim+1,epsilon=10,rang=NULL,optimmethod=c("Nelder-Mead","Newuoa"),verbose=0,scale=TRUE,normed=TRUE, accuracy = 1e-7, itmax = 100000,...)
{
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    kappa <- theta[1]
    lambda <- theta[2]
    nu <- theta[3]
    plot <- FALSE
    if(verbose>0) cat("Minimizing coploss with kappa=",kappa,"lambda=",lambda,"nu=",nu,"\n")
    if(missing(optimmethod)) optimmethod <- "Newuoa"
    if(missing(rang))
        #perhaps put this into the optimization function?
          {
           if(verbose>1) cat ("Fitting configuration for rang. \n")    
           initsol <- stops::powerStressFast(delta,kappa=kappa,lambda=lambda,nu=nu,weightmat=weightmat,ndim=ndim)
           init0 <- initsol$conf
           if(isTRUE(scale)) init0 <- scale(init0)
           crp <- stops::cordillera(init0,q=q,minpts=minpts,epsilon=epsilon,scale=scale)$reachplot
           cin <- max(crp)
           rang <- c(0,1.5*cin)  
           if(verbose>1) cat("dmax is",max(rang),". rang is",rang,"\n")
           }
      if(is.null(rang) && verbose > 1) cat("rang=NULL which makes the cordillera a goodness-of-clustering relative to the largest distance of each given configuration \n") 
      if(missing(cordweight))
               {
                 #cordweight how to fix? here we do not fix for lambda=1, kappa=1, nu=1 but for cordweight=0, so it is stress/cord for initial solution
                 if(verbose>1) cat ("Fitting configuration for cordweight. \n")     
                 initsol <- stops::powerStressFast(delta,kappa=kappa,lambda=lambda,nu=nu,weightmat=weightmat,ndim=ndim)
                 initcord <- stops::cordillera(initsol$conf,q=q,epsilon=epsilon,minpts=minpts,rang=rang,scale=scale)
                 initcorrd <- initcord$normed
                 if(identical(normed,FALSE)) initcorrd <- initcord$raw
                 cordweight <- initsol$stress/initcorrd #use stress.m or stress?
                 if(verbose>1) cat("Weights are stressweight=",stressweight,"cordweight=",cordweight,"\n")
             }
    r <- kappa
    p <- ndim
    deltaorig <- delta
    delta <- delta^lambda
    weightmato <- weightmat
    weightmat <- weightmat^nu
    weightmat[!is.finite(weightmat)] <- 1 #new
    deltaold <- delta
    delta <- delta / enorm (delta, weightmat) #sum=1
    xold <- init
    if(is.null(init)) xold <- stops::torgerson (delta, p = p)
    xold <- xold/enorm(xold) 
    copsf <- function(x,delta,p,weightmat,stressweight,cordweight,q,minpts,epsilon,rang,scale,normed=normed,...)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=p)
             delta <- delta/enorm(delta,weightmat)
             x <- x/enorm(x)
             ds <- (2*as.matrix(dist(x)))^kappa
             #ds <- (2*sqrt(sqdist(x)))^kappa
             ds <- ds/enorm(ds)
             #print(ds)
             stressi <- sum(weightmat*(ds-delta)^2)/2
             #stressi <- sum(weightmat*(ds-delta)^2)/sum(weightmat*(ds^2)) # sqrt stress 1 on the normalized transformed proximities and distances; we use this as the value returned by print
             corrd <- stops::cordillera(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale)
 #            corrd <- stops::cordillera(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,plot,scale=scale,...)
             struc <- corrd$raw
             if(normed) {
                        struc <- corrd$normed
                       }
             ic <- stressweight*stressi - cordweight*struc
             if(verbose>2) cat("coploss =",ic,"mdsloss =",stressi,"OC =",struc,"kappa =",kappa,"lambda =",lambda,"nu=",nu,"\n")
             ic
            }
     if(verbose>1) cat("Starting Minimization with",optimmethod,":\"n")
     if(optimmethod=="Newuoa") {
         optimized <- minqa::newuoa(xold,function(par) copsf(par,delta=delta,p=p,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose),...)
         xnew <- matrix(optimized$par,ncol=2)
         itel <- optimized$feval
         ovalue <-optimized$fval
     }
     if(optimmethod=="Nelder-Mead") {
         optimized <- optim(xold,function(par) copsf(par,delta=delta,p=p,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed),control=list(maxit=itmax,trace=verbose),...)
         xnew <- optimized$par
         itel <- optimized$counts[[1]]
         ovalue <-optimized$val 
     }    
     attr(xnew,"dimnames")[[2]] <- paste("D",1:p,sep="")
     attr(xnew,"dimnames")[[1]] <- rownames(delta)
     #doutm <- (2*sqrt(sqdist(xnew)))^kappa  #fitted powered euclidean distance but times two
     doutm <- as.matrix(dist(xnew)^kappa)
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
     #the following stress versions differ by how the distances and proximities are normalized; either both are normalized (stressen,stressen1), only proximities are normalized (stresse, stresse1), nothing is normalized (stressr, stressn, stresss)
     stressr <- sum(weightmat*(dout-deltaold)^2) #raw stress on the observed proximities
     stresse <- sum(weightmat*(dout-delta)^2) #raw stress on the normalized proximities
     stressen <- sum(weightmat*(doute-delta)^2) #raw stress on the normalized proximities and normalized distances 
     stressen1 <- sqrt(sum(weightmat*(doute-delta)^2)/sum(weightmat*(doute^2))) # sqrt stress 1 on the normalized transformed proximities and distances; we use this as the value returned by print
     stress1 <- sqrt(stressr/sum(weightmat*(dout^2)))  #stress 1 on the original proximities 
     stresse1 <- sqrt(stresse/sum(weightmat*(dout^2)))  #stress 1 on the normalized proximities
     stressn <- stressr/(sum(weightmat*deltaold^2)) #normalized to the maximum stress delta^2*lambda as the normalizing constant (was defualt until v. 0.0-16)
     stresss <- sqrt(stressn) #sqrt of stressn
     if(verbose>0) cat("*** stress (both normalized):",stressen,"; stress 1 (both normalized - default reported):",stressen1,"; sqrt raw stress (both normalized):",sqrt(stressen),"; raw stress (original data):",stressr,"; stress 1 (original data):",stress1,"; explicitly normed stress (original data):",stressn,"; sqrt explicitly normed stress (original data - used in STOPS):",stresss,"; raw stress (proximities normalized):",stresse,"; stress 1 (proximities normalized):", stresse1,"; from optimization: ",ovalue,"\n")   
    out <- list(delta=deltaold, obsdiss=delta, confdiss=dout, conf = xnew, pars=c(kappa,lambda,nu), niter = itel, stress=stressen1, spp=spp, ndim=p, model="Coploss NEWUOA", call=match.call(), nobj = dim(xnew)[1], type = "coploss", gamma=NA, stress.m=sqrt(stressen1), stress.r=stressr/2, stress.n=stressn, stress.1=stress1, stress.s=stresss,stress.e=stresse,stress.en=stressen, stress.en1=stressen1,stress.e1=stresse1, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat)
    out$par <- theta
    out$loss <- "coploss"
    out$OC <- stops::cordillera(out$conf,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale)
    out$coploss <- ovalue
    out$optim <- optimized
    out$stressweight <- stressweight
    out$cordweight <- cordweight
    out$call <- match.call()
    out$optimethod <- optimmethod
    out$losstype <- out$loss
    out$nobj <- dim(out$conf)[1]
    class(out) <- c("cops","smacofP","smacofB","smacof")
    out
}


#' Fitting a COPS Model (Variant 1).
#'
#' Minimizing Coploss for a clustered Power Stress MDS configuration with given hyperparameters theta.
#'
#' @param delta numeric matrix or dist object of a matrix of proximities
#' @param kappa power transformation for fitted distances
#' @param lambda power transformation for proximities
#' @param nu power transformation for weights
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances if it exists), the second lambda (for the observed proximities if it exist), the third is nu (for the weights if it exists) . If less than three elements are is given as argument, it will be recycled. Defaults to 1 1 1. Will override any kappa, lmabda, nu parameters if they are given and do not match
#' @param weightmat (optional) a matrix of nonnegative weights; defaults to 1 for all off diagonals
#' @param ndim number of dimensions of the target space
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param cordweight weight to be used for the cordillera; if missing gets estimated from the initial configuration as  
#' @param q the norm of the corrdillera; defaults to 1
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to ndim+1
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the minimum reachabilities to be considered. If missing it is found from the initial configuration by taking 1.5 times the maximal minimum reachability of the initial fit. If NULL it will be normed to each configuration's minimum and maximum distance, so an absolute value of goodness-of-clusteredness. Note that the latter is not necessarily desirable when comparing configurations for their relative clusteredness. See also \code{\link{cordillera}}     
#' @param optimmethod What optimizer to use? Defaults to NEWUOA, Nelder-Mead is also supported.
#' @param verbose numeric value hat prints information on the fitting process; >2 is very verbose
#' @param normed should the cordillera be normed; defaults to TRUE
#' @param scale should the configuration be scaled to mean=0 and sd=1 for calculating the cordillera? Defaults to TRUE
#' @param accuracy numerical accuracy, defaults to 1e-8
#' @param itmax maximum number of iterations. Defaults to 100000
#' @param ... additional arguments to be passed to the optimization procedure
#'
#'@return A list with the components
#'         \itemize{
#'         \item coploss: the weighted loss value
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
#'dis<-as.matrix(smacof::kinshipdelta)
#'
#'#Coploss with equal weight to stress and cordillera 
#'res1<-coplossMin(dis,stressweight=0.5,cordweight=0.5) 
#'res1
#'summary(res1)
#'plot(res1)  #super clustered
#'
#' @importFrom stats dist as.dist optim
#' @importFrom minqa newuoa
#' 
#' 
#'@keywords clustering multivariate
#'@export
coplossMin <- function (delta, kappa=1, lambda=1, nu=1, theta=c(kappa,lambda,nu),weightmat=1-diag(nrow(delta)),  ndim = 2, init=NULL, stressweight=1,cordweight,q=1,minpts=ndim+1,epsilon=10,rang=NULL,optimmethod=c("Nelder-Mead","Newuoa"),verbose=0,scale=TRUE,normed=TRUE, accuracy = 1e-7, itmax = 100000,...)
{
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    kappa <- theta[1]
    lambda <- theta[2]
    nu <- theta[3]
    plot <- FALSE
    if(verbose>0) cat("Minimizing coploss with kappa=",kappa,"lambda=",lambda,"nu=",nu,"\n")
    if(missing(optimmethod)) optimmethod <- "Newuoa"
    if(missing(rang))
        #perhaps put this into the optimization function?
          {
           if(verbose>1) cat ("Fitting configuration for rang. \n")    
           initsol <- stops::powerStressFast(delta,kappa=kappa,lambda=lambda,nu=nu,weightmat=weightmat,ndim=ndim)
           init0 <- initsol$conf
           if(isTRUE(scale)) init0 <- scale(init0)
           crp <- stops::cordillera(init0,q=q,minpts=minpts,epsilon=epsilon,scale=scale)$reachplot
           cin <- max(crp)
           rang <- c(0,1.5*cin)  
           if(verbose>1) cat("dmax is",max(rang),". rang is",rang,"\n")
           }
      if(is.null(rang) && verbose > 1) cat("rang=NULL which makes the cordillera a goodness-of-clustering relative to the largest distance of each given configuration \n") 
      if(missing(cordweight))
               {
                 #cordweight how to fix? here we do not fix for lambda=1, kappa=1, nu=1 but for cordweight=0, so it is stress/cord for initial solution
                 if(verbose>1) cat ("Fitting configuration for cordweight. \n")     
                 initsol <- stops::powerStressFast(delta,kappa=kappa,lambda=lambda,nu=nu,weightmat=weightmat,ndim=ndim)
                 initcord <- stops::cordillera(initsol$conf,q=q,epsilon=epsilon,minpts=minpts,rang=rang,scale=scale)
                 initcorrd <- initcord$normed
                 if(identical(normed,FALSE)) initcorrd <- initcord$raw
                 cordweight <- initsol$stress/initcorrd #use stress.m or stress?
                 if(verbose>1) cat("Weights are stressweight=",stressweight,"cordweight=",cordweight,"\n")
             }
    r <- kappa/2
#    p <- ndim
    deltaorig <- delta
    delta <- delta^lambda
    weightmato <- weightmat
    weightmat <- weightmat^nu
    weightmat[!is.finite(weightmat)] <- 1 #new
    deltaold <- delta
    delta <- delta / enorm (delta, weightmat) #sum=1
    xold <- init
    if(is.null(init)) xold <- stops::powerStressMin(delta,kappa=kappa,lambda=lambda,nu=nu,ndim=ndim)$conf
    xold <- xold/enorm(xold) 
    copsf <- function(x,delta,r,ndim,weightmat,stressweight,cordweight,q,minpts,epsilon,rang,scale,normed,...)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=ndim)
             delta <- delta/enorm(delta,weightmat)
             x <- x/enorm(x)
             #ds <- (2*as.matrix(dist(x)))^kappa
             #ds <- (2*sqrt(sqdist(x)))^kappa
             #ds <- ds/enorm(ds)
             #print(ds)
             dnew <- sqdist (x)
             rnew <- sum (weightmat * delta * mkPower (dnew, r))
             nnew <- sum (weightmat * mkPower (dnew,  2*r))
             anew <- rnew / nnew
             stressi <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
             #stressi <- sum(weightmat*(ds-delta)^2)/2
             #stressi <- sum(weightmat*(ds-delta)^2)/sum(weightmat*(ds^2)) # sqrt stress 1 on the normalized transformed proximities and distances; we use this as the value returned by print
         #    corrd <- stops::cordillera(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale)
             corrd <- stops::cordillera(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,...)
             struc <- corrd$raw
             if(normed) {
                        struc <- corrd$normed
                       }
             ic <- stressweight*stressi - cordweight*struc
             if(verbose>2) cat("coploss =",ic,"mdsloss =",stressi,"OC =",struc,"kappa =",kappa,"lambda =",lambda,"nu=",nu,"\n")
             ic
           }
     if(verbose>1) cat("Starting Minimization with",optimmethod,":\"n")
     if(optimmethod=="Newuoa") {
         optimized <- minqa::newuoa(xold,function(par) copsf(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- optimized$feval
         ovalue <-optimized$fval
     }
     if(optimmethod=="Nelder-Mead") {
         optimized <- optim(xold,function(par) copsf(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed),control=list(maxit=itmax,trace=verbose),...)
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
     stressen <- sum(weightmat*(doute-delta)^2) #raw stress on the normalized proximities and normalized distances 
     if(verbose>1) cat("*** stress (both normalized - for COPS/STOPS):",stress,"; stress 1 (both normalized - default reported):",sqrt(stress),"; stress manual (for debug only):",stressen,"; from optimization: ",ovalue,"\n")   
    out <- list(delta=deltaold, obsdiss=delta, confdiss=dout, conf = xnew, pars=c(kappa,lambda,nu), niter = itel, stress=sqrt(stress), spp=spp, ndim=ndim, model="Coploss NEWUOA", call=match.call(), nobj = dim(xnew)[1], type = "coploss", gamma=NA, stress.m=stress, stress.en=stressen, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat)
    out$par <- theta
    out$loss <- "coploss"
    out$OC <- stops::cordillera(out$conf,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale)
    out$coploss <- ovalue
    out$optim <- optimized
    out$stressweight <- stressweight
    out$cordweight <- cordweight
    out$call <- match.call()
    out$optimethod <- optimmethod
    out$losstype <- out$loss
    out$nobj <- dim(out$conf)[1]
    class(out) <- c("cops","smacofP","smacofB","smacof")
    out
}


#' Fitting a COPS Model by shrinking residuals to Zero (COPS-0).
#'
#' Minimizing coploss by shrinking residulas to zero to achieve a clustered Power Stress MDS configuration with given hyperparameters theta.
#'
#' @param delta numeric matrix or dist object of a matrix of proximities
#' @param kappa power transformation for fitted distances
#' @param lambda power transformation for proximities
#' @param nu power transformation for weights
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances if it exists), the second lambda (for the observed proximities if it exist), the third is nu (for the weights if it exists) . If less than three elements are is given as argument, it will be recycled. Defaults to 1 1 1. Will override any kappa, lmabda, nu parameters if they are given and do not match
#' @param weightmat (optional) a matrix of nonnegative weights; defaults to 1 for all off diagonals
#' @param ndim number of dimensions of the target space
#' @param init (optional) initial configuration
#' @param cordweight weight to be used for the shrinkage; defaults to 1
#' @param q used in cordillera and shrink matrix and controls the effect of using the norm; defaults to 2 (least squares MDS)
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to ndim+1
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the minimum reachabilities to be considered. If missing it is found from the initial configuration by taking 1.5 times the maximal minimum reachability of the initial fit. If NULL it will be normed to each configuration's minimum and maximum distance, so an absolute value of goodness-of-clusteredness. Note that the latter is not necessarily desirable when comparing configurations for their relative clusteredness. See also \code{\link{cordillera}}
#' @param scaleX should X be scaled; defaults to TRUE
#' @param enormX should X be enormed; defaults to FALSE
#' @param scaleB should X be scaled for the shrink matrix; defaults to TRUE.
#' @param scaleC should X be scaled for the OPTICS Cordillera; defaults to TRUE. These parameter lets one tweak the way the shrinkage works and how its quantified; the defaults lead usually to a sensible result. It might be that some scale versions weill be depreciated in future versions. 
#' @param optimmethod What optimizer to use? Defaults to NEWUOA, Nelder-Mead is also supported.
#' @param verbose numeric value hat prints information on the fitting process; >2 is very verbose
#' @param accuracy numerical accuracy, defaults to 1e-8
#' @param itmax maximum number of iterations. Defaults to 100000
#' @param ... additional arguments to be passed to the optimization procedure
#'
#' @return A list with the components
#'         \itemize{
#'         \item coploss: the weighted loss value
#'         \item OC: the Optics cordillera
#'         \item optim: the object returned from the optimization procedure
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#'         \item cordillera: the OPTICS cordillera object
#' }
#' 
#' @examples
#' dis<-as.matrix(smacof::kinshipdelta)
#'
#' #Coploss with shrinkage to 0 
#' res1<-shrinkCoploss(dis,cordweight=1) 
#' res1
#' summary(res1)
#' plot(res1)  #super clustered
#'
#' @importFrom stats dist as.dist optim
#' @importFrom minqa newuoa
#' 
#' 
#' @keywords clustering multivariate
#' @export
shrinkCoploss <- function (delta, kappa=1, lambda=1, nu=1, theta=c(kappa,lambda,nu),weightmat=1-diag(nrow(delta)),  ndim = 2, init=NULL,cordweight=1,q=2,minpts=ndim+1,epsilon=10,rang=NULL,optimmethod=c("Nelder-Mead","Newuoa"),verbose=0,scaleX=TRUE,enormX=FALSE,scaleB=TRUE,scaleC=TRUE,accuracy = 1e-7, itmax = 100000,...)
{
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    kappa <- theta[1]
    lambda <- theta[2]
    nu <- theta[3]
    plot <- FALSE
    if(verbose>0) cat("Minimizing coploss with kappa=",kappa,"lambda=",lambda,"nu=",nu,"\n")
    if(missing(optimmethod)) optimmethod <- "Newuoa"
    if(missing(rang))
        #perhaps put this into the optimization function?
          {
           if(verbose>1) cat ("Fitting configuration for rang. \n")    
           initsol <- stops::powerStressFast(delta,kappa=kappa,lambda=lambda,nu=nu,weightmat=weightmat,ndim=ndim)
           init0 <- initsol$conf
           if(isTRUE(scaleX)) init0 <- scale(init0)
           crp <- stops::cordillera(init0,q=q,minpts=minpts,epsilon=epsilon,scale=scaleC)$reachplot
           cin <- max(crp)
           rang <- c(0,1.5*cin)  
           if(verbose>1) cat("dmax is",max(rang),". rang is",rang,"\n")
           }
      if(is.null(rang) && verbose > 1) cat("rang=NULL which makes the cordillera a goodness-of-clustering relative to the largest distance of each given configuration \n") 
    r <- kappa/2
    deltaorig <- delta
    delta <- delta^lambda
    weightmato <- weightmat
    weightmat <- weightmat^nu
    weightmat[!is.finite(weightmat)] <- 1 #new
    deltaold <- delta
    delta <- delta / enorm (delta, weightmat) #sum=1
    xold <- init
    if(is.null(init)) xold <- stops::powerStressFast(delta,kappa=kappa,lambda=lambda,nu=nu,ndim=ndim)$conf
    if(enormX) xold <- xold/enorm(xold)
    if(scaleX) xold <- scale(xold)
    shrinkcops <- function(x,delta,r,ndim,weightmat,cordweight,q,minpts,epsilon,rang,scaleX,enormX,scaleB,...)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=ndim)
             #try these variants again with kinship and cali:
             #it all about how to normalize so that the shrinkage will play its part
             #take care of r and so if sqrt(sqdist()) use r*2
             #delta enormed, x scaled + enormed; looks good! -> looks best? 
             #delta enormed, x scaled, dnew enormed; looks ok like #2 but a bit better
             #delta enormed, x enormed, dnew normal; looks ok with clusters for kinship but wrong clusters; closest snew and mdsloss
             if(scaleX) x <- scale(x)
             if(enormX) x <- x/enorm(x)
             delta <- delta/enorm(delta,weightmat)
             dnew <- sqdist(x)
             #dnew <- sqrt(sqdist(x))
             #dnew <- dnew/enorm(dnew,weightmat)
             #dnew <- dnew^2
#r <- 2*r
             rnew <- sum (weightmat * delta * mkPower (dnew, r))
             nnew <- sum (weightmat * mkPower (dnew,  2*r))
             anew <- rnew / nnew
             snew <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
             resen <- abs(mkPower(dnew,r)-delta)
             #resen <- abs(dnew-delta)
             #shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,...)
             shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scaleB=scaleB) 
             #shrinkres <- resen-cordweight*resen*shrinkb/(resen+shrinkb)
             shrinkres <- resen*(1-cordweight*(shrinkb/(resen+shrinkb)))
             #shrinkres <- resen
             diag(shrinkres) <- 0
             #TODO check for increasing residual
             ic <- sum(shrinkres^2)
             if(verbose>3) cat("coploss =",ic,"mdslossm =",sum(resen^2),"delta(cop/mds)=",ic-sum(resen^2),"mdslosss =",snew,"delta(mds/sma)=",sum(resen^2)-snew,"\n")
             #delta cops/mds should be positive if cordweight is too high,no?
             ic
           }
    ## shrinkcops <- function(x,delta,r,ndim,weightmat,cordweight,q,minpts,epsilon,rang,...)
    ##        {
    ##          #TODO: Need to decide whether we use dist() and 2*r throughout or sqdist() and r - for compatibility with the other functions especially coplossMin  
    ##          if(!is.matrix(x)) x <- matrix(x,ncol=ndim)
    ##          delta <- delta/enorm(delta,weightmat)
    ##          x <- x/enorm(x)
    ##          dnew <- sqdist(x)   #sqdist so in power only ^r
    ##          #dnew <- as.matrix(dist(x))   #alternative 
    ##          #rnew <- sum (weightmat * delta * mkPower (dnew, r))
    ##          #nnew <- sum (weightmat * mkPower (dnew,  2*r))
    ##          #anew <- rnew / nnew
    ##          resen <- abs(mkPower(dnew,r)-delta)
    ##          #resen <- abs(mkPower(dnew,2*r)-delta) #alternative
    ##          shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,...)
    ##          #shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale) 
    ##          #shrinkres <- resen-((cordweight*resen*shrinkb)/(resen+shrinkb))
    ##          shrinkres <- resen*(1-cordweight*(shrinkb/(resen+shrinkb)))
    ##          diag(shrinkres) <- 0
    ##          #stressi <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
    ##          #corrd <- stops::cordillera(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,...)
    ##          #struc <- corrd$raw
    ##          #if(normed) {
    ##          #           struc <- corrd$normed
    ##          #          }
    ##           #ic <- stressweight*stressi - cordweight*struc
    ##          #
    ##          ic <- sum(shrinkres^2)/2
    ##          if(verbose>2) cat("coploss =",ic,"mdsloss =",sum(resen^2)/2,"kappa =",kappa,"lambda =",lambda,"nu=",nu,"\n")
    ##          ic
    ##        }
     if(verbose>1) cat("Starting Minimization with",optimmethod,":\"n")
     if(optimmethod=="Newuoa") {
         optimized <- minqa::newuoa(xold,function(par) shrinkcops(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat, cordweight=cordweight, q=q,minpts=minpts,epsilon=epsilon,rang=rang,scaleX=scaleX,scaleB=scaleB,enormX=enormX),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- optimized$feval
         ovalue <-optimized$fval
     }
     if(optimmethod=="Nelder-Mead") {
         optimized <- optim(xold,function(par) shrinkcops(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,
                       cordweight=cordweight, q=q,minpts=minpts,epsilon=epsilon,rang=rang,scaleX=scaleX,scaleB=scaleB,enormX=enormX),control=list(maxit=itmax,trace=verbose),...)
         xnew <- optimized$par
         itel <- optimized$counts[[1]]
         ovalue <-optimized$val 
     }
     if(enormX) xnew <- xnew/enorm(xnew)
     if(scaleX) xnew <- scale(xnew)
     #dnew <- as.matrix(dist (xnew)^2) #alternative
     dnew <- sqdist (xnew)
     rnew <- sum (weightmat * delta * mkPower (dnew, r))
     nnew <- sum (weightmat * mkPower (dnew,  2*r))
     anew <- rnew / nnew
     stress <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
     attr(xnew,"dimnames")[[1]] <- rownames(delta)
     attr(xnew,"dimnames")[[2]] <- paste("D",1:ndim,sep="")
     doutm <- (sqrt(dnew))^kappa  #fitted powered euclidean distance
     #doutm <- as.matrix(dist(xnew)^kappa)  #alternative 
     deltam <- delta
     deltaorigm <- deltaorig
     deltaoldm <- deltaold
     resmat <- deltam - doutm
     delta <- stats::as.dist(delta)
     deltaorig <- stats::as.dist(deltaorig)
     deltaold <- stats::as.dist(deltaold)
     doute <- doutm/enorm(doutm)
     doute <- stats::as.dist(doute)
     dout <- stats::as.dist(doutm)
     spp <- colMeans(resmat)
     weightmatm <-weightmat
     weightmat <- stats::as.dist(weightmatm)
     stressen <- sum(weightmatm*resmat^2)/2 #raw stress on the normalized proximities and normalized distances 
     if(verbose>1) cat("*** stress (both normalized - for COPS/STOPS):",stress,"; stress 1 (both normalized - default reported):",sqrt(stress),"; stress manual (for debug only):",stressen,"; from optimization: ",ovalue,"\n")   
    out <- list(delta=deltaold, obsdiss=delta, confdiss=dout, conf = xnew, pars=c(kappa,lambda,nu), niter = itel, stress=sqrt(stress), spp=spp, ndim=ndim, model="Coploss NEWUOA", call=match.call(), nobj = dim(xnew)[1], type = "coploss", gamma=NA, stress.m=stress, stress.en=stressen, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat)
    out$par <- theta
    out$loss <- "coploss"
    out$OC <- stops::cordillera(out$conf,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scaleC)
    out$coploss <- ovalue
    out$optim <- optimized
    out$cordweight <- cordweight
    out$stressweight <- 1
    out$call <- match.call()
    out$optimethod <- optimmethod
    out$losstype <- out$loss
    out$nobj <- dim(out$conf)[1]
    class(out) <- c("cops","smacofP","smacofB","smacof")
    out
}


#' Finding the shrinkage matrix for COPS-0
#'
#' @param x numeric matrix
#' @param q the norm to be used
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to ndim+1
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the minimum reachabilities to be considered. If missing it is found from the initial configuration.
#' @param scaleB shoudl the X matrix be scaled before calculating the shrinkage weights; defaults to TRUE (which is sensible for relative shrinkage)
#' @param ... additional arguments to be passed to the OPTICS algorithm procedure
#' 
#' @importFrom dbscan optics 
#' 
#'@export
shrinkB <- function(x,q=1,minpts=2,epsilon=10,rang=NULL,scaleB=TRUE,...)
     {
       shift1 <- function (v) {
             vlen <- length(v)
             sv <- as.vector(numeric(vlen)) 
             sv[-1] <- v[1:(vlen - 1)]
             sv
        }
        if(scaleB) x <- scale(x)
        N <- dim(x)[1]
        optres<-dbscan::optics(x,minPts=minpts,eps=epsilon,...)
        optord <- optres$order
        optind <- 1:length(optres$order)
        indordered <- optind[optord]
        predec <- shift1(indordered)
        reachdist <- optres$reachdist[optres$order]
        reachdist[!is.finite(reachdist)] <- ifelse(is.null(rang),max(reachdist[is.finite(reachdist)]),max(rang))
        reachdiffs <- c(NA,abs(diff(reachdist)))
        mats <- cbind(indordered,predec,reachdiffs)
        Bmat <- matrix(0,ncol=N,nrow=N)
        for (i in 2:N) {
            indo <- mats[i,]
            Bmat[indo[1],indo[2]] <- Bmat[indo[2],indo[1]] <- indo[3]/(2^(1/q))
        }
        return(Bmat)
     }



#' High Level COPS Function
#'
#' Minimizing coploss for a clustered MDS configuration. Allows to choose COPS-0 (finding a configuration from coploss with residual shrinkage to zero) and COPS-C (finding a configuration from coploss with cordillera penalty) and profile COPS (finding hyperparameters for MDS models with power transformations). It is wrapper for shrinkCoploss, coplossMin and pcops.
#'
#'@param dis a dissimilarity matrix or a dist object
#'@param variant a character string specifying which variant of COPS to fit. Allowed is any of the following "0","1","2","Variant0","Variant1","Variant2","v0","v1","v2","COPS-0","COPS-C","P-COPS","configuration-0","configuration-c","profile","coploss-0","coploss-c","p-coploss". Defaults to "COPS-C".
#'@param ... arguments to be passed to shrinkCoploss (for Variant 0) coplossMin (for Variant 1) or pcops (for Variant 2). See also \code{\link{shrinkCoploss}} \code{\link{coplossMin}} or \code{\link{pcops}}
#'
#'@return For Variant 0 see \code{\link{shrinkCoploss}}, Variant 1 see \code{\link{coplossMin}}, for Variant 2 see \code{\link{pcops}}
#' 
#'@examples
#' \donttest{
#'dis<-as.matrix(smacof::kinshipdelta)
#'
#'#COPS-C with equal weight to stress and cordillera 
#'res1<-cops(dis,variant="COPS-C",stressweight=0.5,cordweight=0.5,minpts=2) 
#'res1
#'summary(res1)
#'plot(res1)
#'
#'#classic mds (i.e. COPS-C with stressweight=1 and cordweight=0 or COPS-0 with cordweight=0)
#'res2<-cops(dis,variant="COPS-C",stressweight=1,cordweight=0,minpts=2) 
#'res2
#'summary(res2)
#'plot(res2)
#' 
#'#procrustes adjusted
#'resadj<-conf_adjust(res2$fit$conf,res1$fit$conf)
#'plot(resadj$ref.conf) #res 2
#'plot(resadj$other.conf) #res 1
#'
#'par(mfrow=c(1,2))
#'plot(res1,"reachplot")
#'plot(res2,"reachplot") 
#'
#'
#'#COPS-0 to improve over an MDS result
#'res0<-powerStressFast(dis)
#'res2a<-cops(dis,variant="COPS-0",cordweight=1,q=2,init=res0$conf,minpts=2) 
#'res2a
#'summary(res2a)
#'plot(res2a)
#'
#'resadj<-conf_adjust(res0$fit$conf,res2a$fit$conf)
#'plot(resadj$ref.conf) #res 0
#'plot(resadj$other.conf) #res 2a
#' 
#'
#'#s-stress type coploss (i.e. kappa=2, lambda=2)
#'res3<-cops(dis,variant="COPS-C",kappa=2,lambda=2,stressweight=0.5,cordweight=0.5) 
#'res3
#'summary(res3)
#'plot(res3)
#'
#'#Sammon stress type coploss
#'ws<-weightmat=1/dis
#'diag(ws)<-1 
#'res4<-cops(dis,variant="COPS-0",nu=-1,weightmat=ws,cordweight=0.5) 
#'res4
#'summary(res4)
#'plot(res4)
#' 
#'#power-stress type profile coploss
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
#' @importFrom stats dist as.dist optim
#' @importFrom minqa newuoa
#' 
#' 
#'@keywords clustering multivariate
#'@export
cops <- function(dis, variant=c("0","1","2","Variant0","Variant1","Variant2","v0","v1","v2","COPS-0","COPS-C","P-COPS","configuration-0","configuration-c","profile","coploss-0","coploss-c","p-coploss","COPS-P","coploss-p"),...)
                 {
                 if(missing(variant)) variant <- "1"
                 if(variant%in%c("1","Variant1","v1","configuration-c","COPS-C","coploss-0")) out <- coplossMin(dis,...)
                 if(variant%in%c("0","Variant0","v0","configuration-0","COPS-0","coploss-c")) out <- shrinkCoploss(dis,...)
                 if(variant%in%c("2","Variant2","v2","profile","p-coploss","P-COPS","COPS-P","coploss-p")) out <- pcops(dis,...)

                 return(out)
                 }

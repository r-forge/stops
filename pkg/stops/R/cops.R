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
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to 2
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the distances (min distance minus max distance). If NULL (default) the cordillera will be normed to each configuration's maximum distance, so an absolute value of goodness-of-clusteredness.
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param plot plot the cordillera
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
#'@keywords multivariate
#'@export
cop_smacofSym <- function(dis,theta=c(1,1,1),ndim=2,weightmat=NULL,init=NULL,...,stressweight=1,cordweight=0.5,q=1,minpts=2,epsilon=10,rang=NULL,verbose=0,plot=FALSE,normed=TRUE,scale=TRUE,stresstype="default") {
  #TODO Unfolding  
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
  #kappa first argument, lambda=second
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  lambda <- theta
  if(length(theta)==3L) lambda <- theta[2]
 # addargs <- list(...)
 # addargs
  fit <- smacofSym(dis^lambda,ndim=ndim,weightmat=weightmat,init=init,verbose=isTRUE(verbose==2),...) #optimize with smacof
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),plot=plot,scale=scale,normed=normed)
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
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to 2
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the distances (min distance minus max distance). If NULL (default) the cordillera will be normed to each configuration's maximum distance, so an absolute value of goodness-of-clusteredness.
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param plot plot the cordillera
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
#'@keywords multivariate
#'@export
cop_elastic <- function(dis,theta=c(1,1,-2),ndim=2,weightmat=1,init=NULL,...,stressweight=1,cordweight=0.5,q=1,minpts=2,epsilon=10,rang=NULL,verbose=0,plot=FALSE,normed=TRUE,scale=TRUE,stresstype="default") {
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
  fit <- smacofSym(dis^lambda,ndim=ndim,weightmat=combwght,init=init,verbose=isTRUE(verbose==2),...) #optimize with smacof
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),plot=plot,scale=scale,normed=normed)
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
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to 2
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the distances (min distance minus max distance). If NULL (default) the cordillera will be normed to each configuration's maximum distance, so an absolute value of goodness-of-clusteredness.
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param plot plot the cordillera
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
#'@keywords multivariate
#'@export
cop_smacofSphere <- function(dis,theta=c(1,1),ndim=2,weightmat=NULL,init=NULL,...,stressweight=1,cordweight=0.5,q=1,minpts=2,epsilon=10,rang=NULL,verbose=0,plot=FALSE,normed=TRUE,scale=TRUE,stresstype="default") {
  #TODO Unfolding  
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
  #kappa first argument, lambda=second
  if(length(theta)>2) stop("There are too many parameters in the theta argument.")
  lambda <- theta
  if(length(theta)==2L) lambda <- theta[2]
  addargs <- list(...)
  addargs
  fit <- smacofSphere(dis^lambda,ndim=ndim,weightmat=weightmat,init=init,verbose=isTRUE(verbose==2),...) #optimize with smacof
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),plot=plot,scale=scale,normed=normed)
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
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to 2
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the distances (min distance minus max distance). If NULL (default) the cordillera will be normed to each configuration's maximum distance, so an absolute value of goodness-of-clusteredness.
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param plot plot the cordillera
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
#' @keywords multivariate
#' @export
cop_sammon <- function(dis,theta=c(1,1,-1),ndim=2,init=NULL,weightmat=NULL,...,stressweight=1,cordweight=0.5,q=1,minpts=2,epsilon=10,rang=NULL,verbose=0,plot=FALSE,scale=TRUE,normed=TRUE,stresstype="default") {
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
  dis <- as.dist(dis)
  fitdis <- dist(fit$points)
  fit$stress.r <- sum(((dis^lambda-fitdis)^2)/dis)
  fit$stress.n <- fit$stress.r/sum(dis)
  fit$stress.m <- sqrt(fit$stress)
  fit$conf <- fit$points
#  fit$pars <- c(kappa,lambda,nu)
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),plot=plot,scale=scale,normed=normed)
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
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to 2
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the distances (min distance minus max distance). If NULL (default) the cordillera will be normed to each configuration's maximum distance, so an absolute value of goodness-of-clusteredness.
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param plot plot the cordillera
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
#'@keywords multivariate
#'@export
cop_sammon2 <- function(dis,theta=c(1,1,-1),ndim=2,weightmat=NULL,init=NULL,...,stressweight=1,cordweight=0.5,q=1,minpts=2,epsilon=10,rang=NULL,verbose=0,plot=FALSE,normed=TRUE,scale=TRUE,stresstype="default") {
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
  fit <- smacofSym(dis^lambda,ndim=ndim,weightmat=combwght,init=init,verbose=isTRUE(verbose==2),...) #optimize with smacof
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),plot=plot,scale=scale,normed=normed)
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
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to 2
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the distances (min distance minus max distance). If NULL (default) the cordillera will be normed to each configuration's maximum distance, so an absolute value of goodness-of-clusteredness.
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param plot plot the cordillera
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
#' @keywords multivariate
#' @export
cop_cmdscale <- function(dis,theta=c(1,1,1),weightmat=NULL,ndim=2,init=NULL,...,stressweight=1,cordweight=0.5,q=1,minpts=2,epsilon=10,rang=NULL,verbose=0,plot=FALSE,scale=TRUE,normed=TRUE,stresstype="default") {
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)==1L) lambda <- theta
  if(length(theta)==2L) lambda <- theta[2]
  if(length(theta)==3L) lambda <- theta[2]
  fit <- stops::cmdscale(dis^lambda,k=ndim,eig=TRUE,...) 
  fit$lambda <- lambda
  fit$kappa <- 1
  fit$nu <- 1
  dis <- as.dist(dis)
  fitdis <- dist(fit$points)
  fit$stress.r <- sum((dis^lambda-fitdis)^2)
  fit$stress.n <- fit$stress.r/sum(dis^(2*lambda))
  fit$stress.m <- sqrt(fit$stress.n)
  fit$conf <- fit$points
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),plot=plot,scale=scale,normed=normed)
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
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to 2
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the distances (min distance minus max distance). If NULL (default) the cordillera will be normed to each configuration's maximum distance, so an absolute value of goodness-of-clusteredness.
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param plot plot the cordillera
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
cop_rstress <- function(dis,theta=c(1,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,cordweight=0.5,q=1,minpts=2,epsilon=10,rang=NULL,verbose=0,plot=FALSE,scale=TRUE,normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),plot=plot,scale=scale,normed=normed)
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
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to 2
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the distances (min distance minus max distance). If NULL (default) the cordillera will be normed to each configuration's maximum distance, so an absolute value of goodness-of-clusteredness.
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param plot plot the cordillera
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
cop_sstress <- function(dis,theta=c(2,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,cordweight=0.5,q=1,minpts=2,epsilon=10,rang=NULL,verbose=0,plot=FALSE,scale=TRUE,normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),plot=plot,scale=scale,normed=normed)
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
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to 2
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the distances (min distance minus max distance). If NULL (default) the cordillera will be normed to each configuration's maximum distance, so an absolute value of goodness-of-clusteredness.
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param plot plot the cordillera
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
cop_powermds <- function(dis,theta=c(1,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,cordweight=0.5,q=1,minpts=2,epsilon=10,rang=NULL,verbose=0,plot=FALSE,scale=TRUE,normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),plot=plot,scale=scale,normed=normed)
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
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to 2
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the distances (min distance minus max distance). If NULL (default) the cordillera will be normed to each configuration's maximum distance, so an absolute value of goodness-of-clusteredness.
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param plot plot the cordillera
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
cop_powersammon <- function(dis,theta=c(1,1,-1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,cordweight=0.5,q=1,minpts=2,epsilon=10,rang=NULL,verbose=0,plot=FALSE,scale=TRUE,normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),plot=plot,scale=scale,normed=normed)
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
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to 2
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the distances (min distance minus max distance). If NULL (default) the cordillera will be normed to each configuration's maximum distance, so an absolute value of goodness-of-clusteredness.
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param plot plot the cordillera
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
cop_powerelastic <- function(dis,theta=c(1,1,-2),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,cordweight=0.5,q=1,minpts=2,epsilon=10,rang=NULL,verbose=0,plot=FALSE,scale=TRUE,normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),plot=plot,scale=scale,normed=normed)
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
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to 2
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the distances (min distance minus max distance). If NULL (default) the cordillera will be normed to each configuration's maximum distance, so an absolute value of goodness-of-clusteredness.
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param plot plot the cordillera
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
cop_powerstress <- function(dis,theta=c(1,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,cordweight=0.5,q=1,minpts=2,epsilon=10,rang=NULL,verbose=0,plot=FALSE,scale=TRUE,normed=TRUE,stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
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
  copobj <- coploss(fit,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=isTRUE(verbose>1),plot=plot,scale=scale,normed=normed)
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
#' @param plot plot the cordillera
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
coploss <- function(obj,stressweight=1,cordweight=0.5,q=1,normed=TRUE,minpts=2,epsilon=10,rang=NULL,verbose=0,plot=FALSE,scale=TRUE,...)
    {
        stressi <- obj$stress.m
        kappa <- obj$kappa
        lambda <- obj$lambda
        nu <- obj$nu
        confs <- obj$conf 
        corrd <- cordillera(confs,q=q,minpts=minpts,epsilon=epsilon,rang=rang,plot=plot,scale=scale,...)
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

#' High Level COPS Function 
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
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to 2
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the minimum reachabilities to be considered. If missing it is found from the initial configuration by taking 1.5 times the maximal minimum reachability of the model with theta=c(1,1). If NULL it will be normed to each configuration's minimum and maximum distance, so an absolute value of goodness-of-clusteredness. Note that the latter is not necessarily desirable when comparing configurations for their relative clusteredness. See also \code{\link{cordillera}}     
#' @param optimmethod What general purpose optimizer to use? Defaults to our adaptive LJ version (ALJ). Also allows particle swarm optimization with s particles ("pso") and simulated annealing ("SANN"). We recommend not using the later with the rstress, sstress and the power stress models. 
#' @param lower The lower contraints of the search region
#' @param upper The upper contraints of the search region 
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param plot plot the cordillera
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
#'\donttest{ 
#'dis<-as.matrix(smacof::kinshipdelta)
#'res1<-cops(dis,loss="strain",lower=0.1,upper=5) #optimum around lambda=0.15
#'res1
#'summary(res1)
#'plot(res1)
#'#only use cordillera for finding the parameters; 
#'res2<-cops(dis,loss="strain",stressweight=0,lower=0.1,upper=5) #optimum around lambda=0.156
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
#'}
#'
#' 
#'@keywords clustering multivariate
#'@export
cops <- function(dis,loss=c("stress","smacofSym","smacofSphere","strain","sammon","rstress","powermds","sstress","elastic","powersammon","powerelastic","powerstress","sammon2","powerstrain"),weightmat=NULL,ndim=2,init=NULL,theta=c(1,1,1),stressweight=1,cordweight,q=1,minpts=2,epsilon=10,rang,optimmethod=c("ALJ","pso","SANN"),lower=c(1,1,0.5),upper=c(5,5,2),verbose=0,plot=FALSE,scale=TRUE,normed=TRUE,s=4,stresstype="default",...)
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
           crp <- cordillera(init0,q=q,minpts=minpts,epsilon=epsilon,scale=scale)$reachplot
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
                 initcorrd <- cordillera(initsol$fit$conf,q=q,epsilon=epsilon,minpts=minpts,rang=rang,scale=scale)$normed 
                 if(identical(normed,FALSE)) initcorrd <- cordillera(initsol$fit$conf,q=q,epsilon=epsilon,minpts=minpts,rang=rang,scale=scale)$raw
                cordweight <- initsol$stress.m/initcorrd  
                #cat("stress.m=",initsol$stress.m,"cord=",initcorrd,"cweight=",cordweight,"\n") 
                if(verbose>1) cat("Weights are stressweight=",stressweight,"cordweight=",cordweight,"\n")
             }
      if(verbose>1) cat("Starting Optimization \n ")
      if(optimmethod=="SANN") {
          opt<- optim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,weightmat=weightmat,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,plot=plot,scale=scale,normed=normed,stresstype=stresstype))$coploss,method="SANN",...)
      }
      if(optimmethod=="pso") {
        addargs <- list(...)
        control <- list(trace=verbose-2,s=s,addargs)
        opt<- pso::psoptim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,weightmat=weightmat,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,plot=plot,scale=scale,normed=normed,stresstype=stresstype))$coploss,lower=lower,upper=upper,control=control)
       }
      if(optimmethod=="ALJ") {
      opt<- ljoptim(theta, function(theta) do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=theta,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,plot=plot,scale=scale,normed=normed,stresstype=stresstype))$coploss,lower=lower,upper=upper,verbose=verbose-2,...)
      }
    thetaopt <- opt$par 
    #refit the optimal version (TODO probably unnecessary if the other functions are properly reimplemented)
    out <- do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=thetaopt,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-2,plot=plot,scale=scale,normed=normed,stresstype=stresstype))
    out$OC <- cordillera(out$fit$conf,q=q,minpts=minpts,epsilon=epsilon,rang=rang,plot=plot,scale=scale)
    out$coploss <- opt$value
    out$optim <- opt
    out$stressweight <- stressweight
    out$cordweight <- cordweight
    out$call <- match.call()
    out$optimethod <- optimmethod
    out$losstype <- loss
    out$nobj <- dim(out$fit$conf)[1]
    if(verbose>1) cat("Found minimum after",opt$counts["function"]," iterations at",round(opt$par,4),"with coploss=",round(out$coploss,4),"and default scaling loss=",round(out$stress.m,4),"and OC=", round(out$OC$normed,4),". Thanks for your patience. \n")
    class(out) <- c("cops","stops")
    out
  }

#'@export
print.cops <- function(x,...)
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
plot.cops <- function(x,plot.type=c("confplot"), main, asp=NA,...)
    {
     if(missing(plot.type)) plot.type <- "confplot"  
     if(plot.type=="reachplot") {
        if(missing(main)) main <- paste("Reachability plot")
        else main <- main
        plot(x$OC,main=main,...)
     } else if(inherits(x$fit,"smacofB") && !inherits(x$fit,"smacofP") && plot.type=="transplot"){
     #ok, old code here has side effects: it changes the smacof object in the cops object; not sure we should do that  
       if(missing(main)) main <- paste("Transformation Plot") 
       plot.smacofP(x$fit,plot.type="transplot",...)
     # invisible(tmp) #I give the changed smacof object back
     }
     else {      
       plot(x$fit,plot.type=plot.type,main=main,asp=asp,...)
   }
 }

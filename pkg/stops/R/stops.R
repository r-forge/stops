
## # Idea for stops function allow an arbitrary number of indices in a weighted multi-objective optimization way; for this use stoplose
## # write stops_foo where foo is the MDS model of interest
## # TODO: also do this with a pareto approach    
## #'  Calculate the weighted multiobjective loss function used in STOPS
## #'
## #' @param obj object returned inside a stop_* function. Uses the stress.m slot for getting the stress. 
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param structures which c-structuredness indices to be included in the loss
## #' @param strucweight the weights of the structuredness indices; defaults to -1/#number of structures
## #' @param strucpars a list of parameters to be passed to the c-structuredness indices in the same order as the values in structures. If the index has no parameters or you want to use the defaults, supply NULL. (alternatively a named list that has the structure name as the element name).
## #' @param type what type of weighted combination should be used? Can be 'additive' or 'multiplicative'.
## #' @param verbose verbose output
## #'
## #' @import cordillera
## #'
## #'
## #' @return a list with calculated stoploss ($stoploss), structuredness indices ($strucinidices) and hyperparameters ($parameters and $theta) 
## #' 
## #' @export
## stoploss<- function(obj,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"),strucweight=rep(-1/length(structures),length(structures)),strucpars,type=c("additive","multiplicative"),verbose=0)
##     {
##         if(missing(strucpars)) strucpars <- vector("list", length(structures))
##         stressi <- obj$stress.m #we use stress.m everytime
##         pars <- obj$pars
##         confs <- obj$conf 
##         if("cclusteredness"%in%structures)
##             {
##               indst <- which(structures=="cclusteredness")  
##               cclusteredness <- do.call(cordillera::cordillera,c(list(confs),strucpars[[indst]]))$normed
##             }
##          if("cregularity"%in%structures)
##             {
##               indst <- which(structures=="cregularity")  
##               cregularity <- do.call(stops::c_regularity,c(list(confs),strucpars[[indst]]))
##             }                          
##         if("clinearity"%in%structures)
##             {
##                indst <- which(structures=="clinearity")
##                clinearity <- do.call(stops::c_linearity,list(confs))
##            }
##         if("cdependence"%in%structures)
##             {
##                indst <- which(structures=="cdependence")
##                cdependence <- do.call(stops::c_dependence,c(list(confs),strucpars[[indst]])) 
##            }
##         if("cmanifoldness"%in%structures)
##             {
##                indst <- which(structures=="cmanifoldness")
##                cmanifoldness <- do.call(stops::c_manifoldness,c(list(confs)))
##            }
##         if("cassociation"%in%structures)
##             {
##                indst <- which(structures=="cassociation")
##                cassociation <- do.call(stops::c_association,c(list(confs),strucpars[[indst]]))
##            }
##         if("cnonmonotonicity"%in%structures)
##             {
##                indst <- which(structures=="cnonmonotonicity")
##                cnonmonotonicity <- do.call(stops::c_nonmonotonicity,c(list(confs),strucpars[[indst]]))
##            }
##         if("cfunctionality"%in%structures)
##             {
##                indst <- which(structures=="cfunctionality")
##                cfunctionality <- do.call(stops::c_functionality,c(list(confs),strucpars[[indst]])) 
##            }
##          if("ccomplexity"%in%structures)
##             {
##                indst <- which(structures=="ccomplexity")
##                ccomplexity <- do.call(stops::c_complexity,c(list(confs),strucpars[[indst]])) 
##            }
##         if("cfaithfulness"%in%structures)
##             {
##                indst <- which(structures=="cfaithfulness")
##                cfaithfulness <- do.call(stops::c_faithfulness,c(list(confs),strucpars[[indst]]))$mda 
##             }
##          if("chierarchy"%in%structures)
##             {
##                indst <- which(structures=="chierarchy")
##                chierarchy <- do.call(stops::c_hierarchy,c(list(confs),strucpars[[indst]]))
##             }
##          if("coutlying"%in%structures)
##             {
##                indst <- which(structures=="coutlying")
##                coutlying <- do.call(stops::c_outlying,c(list(confs),strucpars[[indst]]))
##             }
##          if("cconvexity"%in%structures)
##             {
##                indst <- which(structures=="cconvexity")
##                cconvexity <- do.call(stops::c_convexity,c(list(confs),strucpars[[indst]]))
##             }
##          if("cskinniness"%in%structures)
##             {
##                indst <- which(structures=="cskinniness")
##                cskinniness <- do.call(stops::c_skinniness,c(list(confs),strucpars[[indst]])) 
##             }
##          if("cstringiness"%in%structures)
##             {
##                indst <- which(structures=="cstringiness")
##                cstringiness <- do.call(stops::c_stringiness,c(list(confs),strucpars[[indst]])) 
##             }
##          if("csparsity"%in%structures)
##             {
##                indst <- which(structures=="csparsity")
##                csparsity <- do.call(stops::c_sparsity,c(list(confs),strucpars[[indst]])) 
##             }
##          if("cclumpiness"%in%structures)
##             {
##                indst <- which(structures=="cclumpiness")
##                cclumpiness <- do.call(stops::c_clumpiness,c(list(confs),strucpars[[indst]])) 
##             }
##          if("cstriatedness"%in%structures)
##             {
##                indst <- which(structures=="cstriatedness")
##                cstriatedness <- do.call(stops::c_striatedness,c(list(confs),strucpars[[indst]]))
##             }
##         if("cinequality"%in%structures)
##             {
##                indst <- which(structures=="cinequality")
##                cinequality <- do.call(stops::c_inequality,c(list(confs),strucpars[[indst]]))
##            }
##         ##TODO add more structures
##         struc <- unlist(mget(structures))
##         ic <- stressi*stressweight + sum(struc*strucweight) 
##         if (type =="multiplicative") ic <- exp(stressweight*log(stressi) + sum(strucweight*log(struc))) #is this what we want? stress/structure or do we want stress - prod(structure)
##         if(verbose>0) cat("stoploss =",ic,"mdsloss =",stressi,"structuredness =",struc,"parameters =",pars,"\n")
##         #return the full combi of stress and indices or only the aggregated scalars; for aSTOPS and mSTOPS we want the latter but for a Pareto approach we want the first; get rid of the sums in ic if the first is wanted  
##         out <- list(stoploss=ic,strucindices=struc,parameters=pars,theta=pars)
##         out
##     }

## #' STOPS version of smacofSym models
## #'
## #' The free parameter is lambda for power transformations the observed proximities. The fitted distances power is internally fixed to 1 and the power for the weights is 1. 
## #'
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the theta vector; must be a scalar for the lambda (proximity) transformation. Defaults to 1.
## #' @param ndim number of dimensions of the target space
## #' @param weightmat (optional) a matrix of nonnegative weights
## #' @param init (optional) initial configuration
## #' @param itmax number of iterations
## #' @param ... additional arguments to be passed to the fitting
## #' @param structures which structuredness indices to be included in the loss
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
## #' @param strucpars the parameters for the structuredness indices
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
## #' 
## #' 
## #' @return A list with the components
## #'    \itemize{
## #'         \item{stress:} the stress-1 (sqrt(stress.m))
## #'         \item{stress.m:} default normalized stress (used for STOPS)
## #'         \item{stoploss:} the weighted loss value
## #'         \item{indices:} the values of the structuredness indices
## #'         \item{parameters:} the parameters used for fitting 
## #'         \item{fit:} the returned object of the fitting procedure
## #'         \item{stopobj:} the stops object
## #' }
## #' 
## #'@keywords multivariate
## #'@import smacof
## #'@import cordillera
## #'@export
## stop_smacofSym <- function(dis, theta=1, ndim=2,weightmat=NULL,init=NULL,itmax=1000,...,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"),stressweight=1,strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
##   theta <- as.numeric(theta)
##   if(is.null(init)) init <- "torgerson"  
##   if(inherits(dis,"dist")) dis <- as.matrix(dis)
##   if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
##   if(missing(type)) type <- "additive"
##   if(length(theta)>1) stop("There are too many parameters in the theta argument.")
##   #if(length(theta)==1) lambda <- theta
##   #if(length(theta)==2) lambda <- theta[2]
##   #if(length(theta)==3) lambda <- theta[2]
##   lambda <- theta
##   fit <- smacof::smacofSym(dis^lambda,ndim=ndim,weightmat=weightmat,init=init,verbose=isTRUE(verbose==2),itmax=itmax,...) #optimize with smacof
##   #fit$kappa <- 1
##   fit$lambda <- theta
##   #fit$nu <- 1
##   fit$stress.1 <- fit$stress
##   fitdis <- as.matrix(fit$confdist)
##   delts <- as.matrix(fit$dhat) 
##   fit$stress.r <- sum(weightmat*(delts-fitdis)^2)
##   fit$stress.m <- fit$stress^2 #fit$stress.r/sum(weightmat*delts^2)
##   fit$pars <- c(lambda=fit$lambda)#c(kappa=fit$kappa,lambda=fit$lambda,rho=fit$nu)
##   fit$deltaorig <- fit$delta^(1/fit$lambda)
##   stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
##   out <- list(stress=fit$stress, stress.r=fit$stress.r,stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices,parameters=stopobj$parameters,fit=fit,stopobj=stopobj) #target functions
##   out
## }

## #' STOPS versions of flexsmacof models (models with a parametric f() transformation to be determined)
## #'
## #' Not functional
## #' 
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param transformation function to transform the proximities or distances; need to be parameterized by theta  
## #' @param theta the theta vector of transformations
## #' @param ndim number of dimensions of the target space
## #' @param weightmat (optional) a matrix of nonnegative weights
## #' @param init (optional) initial configuration
## #' @param structures which structuredness indices to be included in the loss
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param strucpars strucpars
## #' @param itmax number of iterations 
## #' @param ... additional arguments to be passed to the fitting function

## #' 
## #' @return A list with the components
## #'    \itemize{
## #'         \item{stress:} the stress
## #'         \item{stress.m:} default normalized stress
## #'         \item{stoploss:} the weighted loss value
## #'         \item{indices:} the values of the structuredness indices
## #'         \item{parameters:} the parameters used for fitting 
## #'         \item{fit:} the returned object of the fitting procedure
## #'         \item{indobj:} the index objects
## #' }
## #' 
## #'@keywords multivariate
## #'@import smacof
## #'@import cordillera 
## #'@export
## stop_flexsmacof <- function(dis,transformation=mkPower2, theta=c(1,1), ndim=2,weightmat=NULL,init=NULL,itmax=1000,...,structures=c("clusteredness","linearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"),stressweight=1,strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0) {
##     if(is.null(init)) init <- "torgerson" 
##   theta <- as.numeric(theta)
##   if(inherits(dis,"dist")) dis <- as.matrix(dis)
##   if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
##   addargs <- list(...)
##   #TODO: Other transformations parametrized by theta; use splines
##   #Transformations must be so that first argument is the dissimilarity matrix and the second the theta parameters 
##   diso <- dis
##   dis <- do.call(transformation,list(diso,theta))
##   diso <- dis
##   dis <- do.call(transformation,list(diso,theta))
##   fit <- smacof::smacofSym(dis,ndim=ndim,weightmat=weightmat,init=init,verbose=isTRUE(verbose==2),itmax=itmax,...) #optimize with smacof
##   fit$stress.1 <- fit$stress
##   fitdis <- as.matrix(fit$confdist)
##   delts <- as.matrix(fit$delta) #That was my choice to not use the normalized deltas but try it ion the original; that is scale and unit free as Buja said
##   fit$stress.r <- sum(weightmat*(delts-fitdis)^2)
##   fit$stress.m <- fit$stress.r/sum(weightmat*delts^2)
##   fit$pars <- theta
##   stopobj <- stoploss(fit,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars)
##   out <- list(stress=fit$stress, stress.r=fit$stress.r/2, stress.m=fit$stress.m/2, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters,fit=fit) #target functions
##   #TODO include the objects of the indices returned as a list? indicesfull=stopobj 
##   out
## }

## #' STOPS versions of elastic scaling models (via smacofSym)
## #'
## #' The free parameter is lambda for power transformations the observed proximities. The fitted distances power is internally fixed to 1 and the power for the weights=delta is -2. Allows for a weight matrix because of smacof.
## #' 
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the theta vector of powers; this must be a scalar of the lambda transformation for the observed proximities. Defaults to 1.
## #' @param ndim number of dimensions of the target space
## #' @param weightmat (optional) a matrix of nonnegative weights (NOT the elscal weights)
## #' @param init (optional) initial configuration
## #' @param itmax number of iterations
## #' @param ... additional arguments to be passed to the fitting procedure
## #' @param structures which structuredness indices to be included in the loss
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
## #' @param strucpars the parameters for the structuredness indices
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
## #' 
## #' 
## #' @return A list with the components
## #'    \itemize{
## #'         \item{stress:} the stress-1 (sqrt(stress.m))
## #'         \item{stress.m:} default normalized stress (used for STOPS)
## #'         \item{stoploss:} the weighted loss value
## #'         \item{indices:} the values of the structuredness indices
## #'         \item{parameters:} the parameters used for fitting 
## #'         \item{fit:} the returned object of the fitting procedure
## #'         \item{stopobj:} the stopobj objects
## #' }
## #'
## #'@importFrom stats dist as.dist
## #'@import smacof 
## #'@keywords multivariate
## #'@export
## stop_elastic <- function(dis,theta=1,ndim=2,weightmat=NULL,init=NULL,itmax=1000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
##                                         #TODO Unfolding
##   if(is.null(init)) init <- "torgerson" 
##   theta <- as.numeric(theta)
##   if(inherits(dis,"dist")) dis <- as.matrix(dis)
##   if(missing(type)) type <- "additive"
##   #kappa first argument, lambda=second
##   if(length(theta)>1) stop("There are too many parameters in the theta argument.")
##   #if(length(theta)<3) theta <- rep(theta,length.out=3)
##   if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
##   lambda <- theta
##   nu <- -2
##   elscalw <- dis^(nu*lambda) #the weighting in elastic scaling
##   diag(elscalw) <- 1
##   combwght <- weightmat*elscalw #combine the user weights and the elastic scaling weights
##   fit <- smacof::smacofSym(dis^lambda,ndim=ndim,weightmat=combwght,init=init,verbose=isTRUE(verbose==2),itmax=itmax,...) #optimize with smacof
##   #fit$kappa <- 1
##   fit$lambda <- lambda
##   #fit$nu <- nu
##   fit$stress.1 <- fit$stress
##   fitdis <- as.matrix(fit$confdist)
##   delts <- as.matrix(fit$delta) 
##   fit$stress.r <- sum(combwght*((delts-fitdis)^2))
##   fit$obsdiss <- fit$dhat
##   fit$stress.m <- fit$stress^2 #fit$stress.r/sum(combwght*delts^2)
##   fit$pars <- c(lambda=fit$lambda)#c(kappa=fit$kappa,lambda=fit$lambda,rho=fit$nu)
##   fit$deltaorig <- fit$delta^(1/fit$lambda)
##   stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
##   out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit,stopobj=stopobj) #target functions
##   out
## }


## #' STOPS versions of smacofSphere models
## #'
## #' The free parameter is lambda for power transformations the observed proximities. The fitted distances power is internally fixed to 1 and the power for the weights is 1. 
## #' 
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the theta vector of powers; this must be a scalar of the lambda transformation for the observed proximities. Defaults to 1.
## #' @param ndim number of dimensions of the target space
## #' @param weightmat (optional) a matrix of nonnegative weights
## #' @param init (optional) initial configuration
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param itmax number of iterations
## #' @param ... additional arguments to be passed to the fitting procedure
## #' @param structures which structuredness indices to be included in the loss
## #' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
## #' @param strucpars the parameters for the structuredness indices
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
## #' 
## #' 
## #' @return A list with the components
## #'    \itemize{
## #'         \item{stress:} the stress
## #'         \item{stress.m:} default normalized stress
## #'         \item{stoploss:} the weighted loss value
## #'         \item{indices:} the values of the structuredness indices
## #'         \item{parameters:} the parameters used for fitting 
## #'         \item{fit:} the returned object of the fitting procedure
## #'          \item{stopobj:} the stopobj object
## #' }
## #'
## #'@import smacof 
## #'@importFrom stats dist as.dist
## #'@keywords multivariate
## #'@export
## stop_smacofSphere <- function(dis,theta=1,ndim=2,weightmat=NULL,init=NULL,itmax=1000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
##                                         #TODO Unfolding
##   if(is.null(init)) init <- "torgerson"
##   theta <- as.numeric(theta)
##   if(inherits(dis,"dist")) dis <- as.matrix(dis)
##   if(missing(type)) type <- "additive"
##   if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
##   #kappa first argument, lambda=second
##   if(length(theta)>1) stop("There are too many parameters in the theta argument.")
##   #if(length(theta)<3) theta <- rep(theta,length.out=3)
##   lambda <- theta
##   fit <- smacof::smacofSphere(dis^lambda,ndim=ndim,weightmat=weightmat,init=init,verbose=isTRUE(verbose==2),itmax=itmax,...) #optimize with smacof
##   #fit$kappa <- 1
##   fit$lambda <- lambda
##   #fit$nu <- 1
##   fit$stress.1 <- fit$stress
##   fitdis <- as.matrix(fit$confdist)
##   delts <- as.matrix(fit$delta)[-1,-1]
##   fit$obsdiss <- fit$dhat 
##   fit$stress.r <- sum(weightmat*(delts-fitdis)^2)
##   fit$stress.m <- fit$stress^2
##   fit$pars <- c(lambda=fit$lambda)#c(kappa=fit$kappa,lambda=fit$lambda,rho=fit$nu)
##   fit$deltaorig <- fit$delta^(1/fit$lambda)
##   stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
##   out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit,stopobj=stopobj) #target functions
##   out
## }


## #' STOPS version of Sammon mapping
## #'
## #' Uses MASS::sammon. The free parameter is lambda for power transformations of the observed proximities. The fitted distances power is internally fixed to 1 and the power for the weights=delta is -1. 
## #' 
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the theta vector of powers; this must be  a scalar of the lambda transformation for the observed proximities. Defaults to 1.
## #' @param ndim number of dimensions of the target space
## #' @param weightmat a matrix of nonnegative weights. Has no effect here.
## #' @param init (optional) initial configuration
## #' @param itmax number of iterations
## #' @param ... additional arguments to be passed to the fitting procedure
## #' @param structures which structuredness indices to be included in the loss
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
## #' @param strucpars the parameters for the structuredness indices
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
## #' 
## #' 
## #' @return A list with the components
## #'    \itemize{
## #'         \item{stress:} the stress
## #'         \item{stress.m:} default normalized stress
## #'         \item{stoploss:} the weighted loss value
## #'         \item{indices:} the values of the structuredness indices
## #'         \item{parameters:} the parameters used for fitting 
## #'         \item{fit:} the returned object of the fitting procedure
## #'          \item{stopobj:} the stopobj object
## #' }
## #'
## #' @importFrom stats dist as.dist
## #' @import cordillera
## #' @keywords multivariate
## #'
## #' 
## #' @export
## stop_sammon <- function(dis,theta=1,ndim=2,init=NULL,weightmat=NULL,itmax=1000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
##   theta <- as.numeric(theta)
##   if(length(theta)>1) stop("There are too many parameters in the theta argument.")
##   if(missing(type)) type <- "additive"
##   if(inherits(dis,"dist")) dis <- as.matrix(dis)
##   if(length(theta)==1L) lambda <- theta
##   #if(length(theta)==2L) lambda <- theta[2]
##   #if(length(theta)==3L) lambda <- theta[2]
##   #nu <- -1
##   fit <- stops::sammon(dis^lambda,k=ndim,y=init,trace=isTRUE(verbose>1),niter=itmax,...)
##   fit$lambda <- lambda
##   #fit$kappa <- 1
##   #fit$nu <- -1
##   N <- length(dis)
##   dis <- stats::as.dist(dis)
##   fitdis <- stats::dist(fit$points)
##   wghts <- 1/dis^lambda 
##   disl <- dis^lambda
##   dhat <-  disl/sqrt(sum(wghts*disl^2))*sqrt(N)
##   lb <- sum(wghts*fitdis*dhat)/sum(wghts*fitdis^2)   #Restrict config so we have a stress in [0,1] just as in smacof. Rest is unchanged. Maybe use this stress for optimization at some point?
##   fitdisnn <- lb*fitdis
##   fit$stress.r <- fit$stress
##   fit$stress <- sqrt(sum(wghts*(dhat-fitdisnn)^2)/N) #sum(wghts*dhat^2))
##   #fit$stress.n <- fit$stress.r/sum(dis^lambda)
##   fit$stress.m <- fit$stress #or stress.r
##   fit$conf <- fit$points
##   fit$pars <- c(lambda=fit$lambda)#c(kappa=fit$kappa,lambda=fit$lambda,rho=fit$nu)
##   stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
##   list(stress=fit$stress, stress.m=fit$stress.m, stress.r=fit$stress.r,stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters,  fit=fit,stopobj=stopobj) #target functions
## }

## #' Another STOPS version of Sammon mapping models (via smacofSym)
## #'
## #' Uses Smacof, so it can deal with a weight matrix too.  The free parameter is lambda for power transformations of the observed proximities. The fitted distances power is internally fixed to 1 and the power for the weights=delta is -1. 
## #'
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the theta vector of powers; this must be a scalar of the lambda transformation for the observed proximities. Defaults to 1.
## #' @param ndim number of dimensions of the target space
## #' @param weightmat (optional) a matrix of nonnegative weights
## #' @param init (optional) initial configuration
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param itmax number of iterations
## #' @param ... additional arguments to be passed to the fitting procedure
## #' @param structures which structuredness indices to be included in the loss
## #' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
## #' @param strucpars the parameters for the structuredness indices
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative'. 
## #' 
## #' 
## #' @return A list with the components
## #'    \itemize{
## #'         \item{stress:} the stress-1 (sqrt(stress.m))
## #'         \item{stress.m:} default normalized stress (used for STOPS)
## #'         \item{stoploss:} the weighted loss value
## #'         \item{indices:} the values of the structuredness indices
## #'         \item{parameters:} the parameters used for fitting 
## #'         \item{fit:} the returned object of the fitting procedure
## #'          \item{stopobj:} the stopobj object
## #' }
## #'
## #'
## #' @importFrom stats dist as.dist
## #' @import smacof
## #' @import cordillera
## #' 
## #'@keywords multivariate
## #'@export
## stop_sammon2 <- function(dis,theta=1,ndim=2,weightmat=NULL,init=NULL,itmax=1000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
##     theta <- as.numeric(theta)
##     if(is.null(init)) init <- "torgerson" 
##   if(inherits(dis,"dist")) dis <- as.matrix(dis)
##   if(missing(type)) type <- "additive"
##   if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1]) 
##   #kappa first argument, lambda=second
##   if(length(theta)>1) stop("There are too many parameters in the theta argument.")
##   #if(length(theta)<3) theta <- rep(theta, length.out=3)
##   lambda <- theta
##   nu <- -1
##   elscalw <- dis^(nu*lambda) #the weighting in 
##   diag(elscalw) <- 1
##   combwght <- weightmat*elscalw #combine the user weights and the elastic scaling weights
##   fit <- smacof::smacofSym(dis^lambda,ndim=ndim,weightmat=combwght,init=init,verbose=isTRUE(verbose==2),itmax=itmax,...) #optimize with smacof
##   #fit$kappa <- 1
##   fit$lambda <- lambda
##   #fit$nu <- nu
##   fit$stress.1 <- fit$stress
##   fitdis <- as.matrix(fit$confdist)
##   delts <- as.matrix(fit$delta)
##   fit$obsdiss <- fit$dhat
##   fit$stress.r <- sum(combwght*((delts-fitdis)^2))
##   fit$stress.m <- fit$stress^2# fit$stress.r/sum(combwght*delts^2)
##   fit$pars <- c(lambda=fit$lambda)#c(kappa=fit$kappa,lambda=fit$lambda,rho=fit$nu)
##   fit$deltaorig <- fit$delta^(1/fit$lambda)
##   stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
##   out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit,stopobj=stopobj) #target functions
##   out
## }



## #' STOPS version of strain
## #'
## #' The free parameter is lambda for power transformations of the observed proximities.
## #'
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the theta vector of powers; this must be a scalar of the lambda transformation for the observed proximities.
## #' @param ndim number of dimensions of the target space
## #' @param weightmat (optional) a matrix of nonnegative weights. Not used. 
## #' @param init (optional) initial configuration
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param ... additional arguments to be passed to the fitting procedure
## #' @param structures which structuredness indices to be included in the loss
## #' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
## #' @param strucpars the parameters for the structuredness indices
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
## #' @param itmax placeholder for compatibility in stops call; not used
## #' 
## #' @return A list with the components
## #'    \itemize{
## #'         \item{stress:} Sqrt of explicitly normalized stress. 
## #'         \item{stress.m:} explictly normalized stress
## #'         \item{stoploss:} the weighted loss value
## #'         \item{indices:} the values of the structuredness indices
## #'         \item{parameters:} the parameters used for fitting 
## #'         \item{fit:} the returned object of the fitting procedure
## #'          \item{stopobj:} the stopobj object
## #' }
## #'
## #' @import cordillera
## #' @importFrom stats dist as.dist
## #' @keywords multivariate
## #' @export
## stop_cmdscale <- function(dis,theta=1,weightmat=NULL,ndim=2,init=NULL,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative"),itmax=NULL) {
##   theta <- as.numeric(theta)
##   if(length(theta)>1) stop("There are too many parameters in the theta argument.")
##   if(missing(type)) type <- "additive"
##   if(length(theta)==1) lambda <- theta
##   #if(length(theta)==2) lambda <- theta[2]
##   #if(length(theta)==3) lambda <- theta[2]
##   fit <- stops::cmdscale(dis^lambda,k=ndim,eig=TRUE,...) 
##   fit$lambda <- lambda
##   #fit$kappa <- 1
##   #fit$nu <- 1
##   dis <- stats::as.dist(dis)
##   fitdis <- stats::dist(fit$points)
##   fit$stress.r <- sum((dis^lambda-fitdis)^2)
##   fit$stress.n <- fit$stress.r/sum(dis^(2*lambda))
##   fit$stress <- sqrt(fit$stress.n)
##   fit$stress.m <- fit$stress.n
##   fit$pars <- c(lambda=fit$lambda)#c(kappa=fit$kappa,lambda=fit$lambda,rho=fit$nu)
##   fit$conf <- fit$points
##   stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
##   list(stress=fit$stress,stress.r=fit$stress.r,stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj) #target functions
## }



## #' STOPS version of isomap to optimize over integer k.
## #'
## #' Free parameter is k. 
## #' 
## #' @details Currently this version is a bit less flexible than the vegan one, as the only allowed parameter for isomap is the theta (k in isomap, no epsilon) and the shortest path is always estimated with argument "shortest". Also note that fragmentedOK is always set to TRUE which means that for theta that is too small only the largest conected group will be analyzed. If that's not wanted just set the theta higher.  
## #' 
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the number of shortest dissimilarities retained for a point (nearest neighbours), the isomap parameter. Must be a numeric scalar. Defaults to 3.
## #' @param ndim number of dimensions of the target space
## #' @param weightmat (optional) a matrix of nonnegative weights
## #' @param init (optional) initial configuration
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param structures which structuredness indices to be included in the loss
## #' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
## #' @param strucpars the parameters for the structuredness indices
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
## #' @param itmax placeholder for compatibility in stops call; not used
## #' 
## #' @return A list with the components
## #'    \itemize{
## #'         \item{stress:} Not really stress but 1-GOF where GOF is the first element returned from cmdscale (the sum of the first ndim absolute eigenvalues divided by the sum of all absolute eigenvalues).
## #'         \item{stress.m:} default normalized stress (sqrt explicitly normalized stress; really the stress this time)
## #'         \item{stoploss:} the weighted loss value
## #'         \item{indices:} the values of the structuredness indices
## #'         \item{parameters:} the parameters used for fitting 
## #'         \item{fit:} the returned object of the fitting procedure
## #'          \item{stopobj:} the stopobj object
## #' }
## #'
## #' @import cordillera
## #' @importFrom stats dist as.dist
## #' @importFrom vegan isomap isomapdist
## #' @keywords multivariate
## #' @export
## stop_isomap1 <- function(dis,theta=3,weightmat=NULL,ndim=2,init=NULL,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative"),itmax=NULL) {
##   theta <- as.numeric(theta)
##   if(length(theta)>1) stop("There are too many parameters in the theta argument.")
##   if(missing(type)) type <- "additive"
##   if(length(theta)==1L) lambda <- theta #I just call the k lambda
##   #if(length(theta)==2L) lambda <- theta[1]
##   #if(length(theta)==3L) lambda <- theta[1]
##   disi <- vegan::isomapdist(dis,k=lambda,path="shortest",fragmentedOK=TRUE)
##   if(length(disi)==0) stop("The distance matrix is of length 0 for the current k. Consider increasing the lower bound of the search region.")
##   fit <- stops::cmdscale(disi,k=ndim,eig=TRUE) 
##   fit$k <- lambda
##   #fit$kappa <- 1
##   #fit$nu <- 1
##   dis <- stats::as.dist(disi)
##   fitdis <- stats::dist(fit$points)
##   fit$stress.r <- sum((disi-fitdis)^2)
##   fit$stress.n <- fit$stress.r/sum(disi^2)
##   fit$stress.m <- fit$stress.n
##   fit$stress <- sqrt(fit$stress.n)
##   fit$pars <- c(k=fit$k)
##   fit$conf <- fit$points
##   stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
##   list(stress=fit$stress,stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj) #target functions
## }


## #' STOPS version of isomap over real epsilon.
## #'
## #' Free parameter is eps. 
## #' 
## #' @details Currently this version is a bit less flexible than the vegan one, as the only allowed parameter for isomap is the theta (epsilon in isomap) and the shortest path is always estimated with argument "shortest". Also note that fragmentedOK is always set to TRUE which means that for theta that is too small only the largest conected group will be analyzed. If that's not wanted just set the theta higher.  
## #' 
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the number of shortest dissimilarities retained for a point (neighbourhood region), the isomap parameter. Defaults to the 0.1 quantile of the empirical distribution of dis.
## #' @param ndim number of dimensions of the target space
## #' @param weightmat (optional) a matrix of nonnegative weights
## #' @param init (optional) initial configuration
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param structures which structuredness indices to be included in the loss
## #' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
## #' @param strucpars the parameters for the structuredness indices
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
## #' @param itmax placeholder for compatibility in stops call; not used
## #' 
## #' @return A list with the components
## #'    \itemize{
## #'         \item{stress:} Not really stress but 1-GOF where GOF is the first element returned from cmdscale (the sum of the first ndim absolute eigenvalues divided by the sum of all absolute eigenvalues).
## #'         \item{stress.m:} default normalized stress (sqrt explicitly normalized stress; really the stress this time)
## #'         \item{stoploss:} the weighted loss value
## #'         \item{indices:} the values of the structuredness indices
## #'         \item{parameters:} the parameters used for fitting 
## #'         \item{fit:} the returned object of the fitting procedure
## #'          \item{stopobj:} the stopobj object
## #' }
## #'
## #' @import cordillera
## #' @importFrom stats dist as.dist quantile
## #' @importFrom vegan isomap isomapdist
## #' @keywords multivariate
## #' @export
## stop_isomap2 <- function(dis,theta=stats::quantile(dis,0.1),weightmat=NULL,ndim=2,init=NULL,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative"),itmax=NULL) {
##   theta <- unique(theta)
##   if(length(theta)>1) stop("There are too many parameters in the theta argument.")
##   if(missing(type)) type <- "additive"
##   if(length(theta)==1L) lambda <- theta
##   #if(length(theta)==2L) lambda <- theta[1]
##   #if(length(theta)==3L) lambda <- theta[1]
##   disi <- vegan::isomapdist(dis,epsilon=lambda,path="shortest",fragmentedOK=TRUE)
##   if(length(disi)==0) stop("The distance matrix is of length 0 for the current epsilon. Consider increasing the lower bound of the search region.")
##   fit <- stops::cmdscale(disi,k=ndim,eig=TRUE) 
##   fit$eps <- lambda
##   #fit$kappa <- 1
##   #fit$nu <- 1
##   dis <- stats::as.dist(disi)
##   fitdis <- stats::dist(fit$points)
##   fit$stress.r <- sum((disi-fitdis)^2)
##   fit$stress.n <- fit$stress.r/sum(disi^2)
##   fit$stress.m <- fit$stress.n
##   fit$stress <- sqrt(fit$stress.m)
##   fit$pars <- c(eps=fit$eps)
##   fit$conf <- fit$points
##   stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
##   list(stress=fit$stress,stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj) #target functions
## }


## #' STOPS version of rstress
## #'
## #' Free parameter is kappa for the fitted distances.
## #'
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the theta vector of powers; this must be a scalar of the kappa transformation for the fitted distances proximities. Defaults to 1. Note the kappa here differs from Jan's version where the parameter was called r and the relationship is r=kappa/2 or kappa=2r.
## #' @param ndim number of dimensions of the target space
## #' @param weightmat (optional) a matrix of nonnegative weights
## #' @param init (optional) initial configuration
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param itmax number of iterations
## #' @param ... additional arguments to be passed to the fitting procedure
## #' @param structures which structuredness indices to be included in the loss
## #' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
## #' @param strucpars the parameters for the structuredness indices
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
## #' 
## #' 
## #' @return A list with the components
## #'    \itemize{
## #'         \item{stress:} the stress
## #'         \item{stress.m:} default normalized stress
## #'         \item{stoploss:} the weighted loss value
## #'         \item{indices:} the values of the structuredness indices
## #'         \item{parameters:} the parameters used for fitting 
## #'         \item{fit:} the returned object of the fitting procedure
## #'          \item{stopobj:} the stopobj object
## #' }
## #'
## #' @import cordillera
## #' @keywords multivariate
## #' @export
## stop_rstress <- function(dis,theta=1,weightmat=NULL,init=NULL,ndim=2,itmax=100000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
##   theta <- as.numeric(theta)
##   if(inherits(dis,"dist")) dis <- as.matrix(dis) 
##   if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
##   if(length(theta)>1) stop("There are too many parameters in the theta argument.")
##   if(missing(type)) type <- "additive"
##   #if(length(theta)<3) theta <- rep(theta,length.out=3)
##   kappa <- theta
##   fit <- powerStressMin(delta=dis,kappa=kappa,lambda=1,nu=1,weightmat=weightmat,init=init,ndim=ndim,verbose=verbose,itmax=itmax,...)
##   fit$kappa <- kappa
##   #fit$lambda <- 1
##   #fit$nu <- 1
##   fit$pars <- c(kappa=fit$kappa)#c(kappa=fit$kappa,lambda=fit$lambda,rho=fit$nu)
##   stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
##   out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
##   out
## }


## #' STOPS version of sstress
## #'
## #' Free parameter is lambda for the observed proximities. Fitted distances are transformed with power 2, weights have exponent of 1. Note that the lambda here works as a multiplicator of 2 (as sstress has f(delta^2)).
## #' 
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the theta vector of powers; this must be a scalar of the lambda transformation for the observed proximities. Defaults to 1. Note that the lambda here works as a multiplicator of 2 (as sstress has f(delta^2)). 
## #' @param weightmat (optional) a matrix of nonnegative weights
## #' @param init (optional) initial configuration
## #' @param ndim the number of dimensions of the target space
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param itmax number of iterations
## #' @param ... additional arguments to be passed to the fitting procedure
## #' @param structures which structuredness indices to be included in the loss
## #' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
## #' @param strucpars the parameters for the structuredness indices
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
## #' 
## #' 
## #' @return A list with the components
## #'    \itemize{
## #'         \item{stress:} the stress
## #'         \item{stress.m:} default normalized stress
## #'         \item{stoploss:} the weighted loss value
## #'         \item{indices:} the values of the structuredness indices
## #'         \item{parameters:} the parameters used for fitting 
## #'         \item{fit:} the returned object of the fitting procedure
## #'          \item{stopobj:} the stopobj object
## #' }
## #' @import cordillera
## #' @keywords multivariate
## #' @export
## stop_sstress <- function(dis,theta=1,weightmat=NULL,init=NULL,ndim=2,itmax=100000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
##   theta <- as.numeric(theta)
##   if(inherits(dis,"dist")) dis <- as.matrix(dis)  
##   if(missing(type)) type <- "additive"
##   if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
##   if(length(theta)>1) stop("There are too many parameters in the theta argument.")
##   #if(length(theta)<3) theta <- rep(theta, length.out=3)
##   lambda <- theta
##   flambda <- lambda*2 #sstress is d^2 and delta^2 so f(delta^2)=delta^(2*1); lambda works in factors of 2  
##   fit <- powerStressMin(delta=dis,kappa=2,lambda=flambda,nu=1,weightmat=weightmat,init=init,ndim=ndim,verbose=verbose,itmax=itmax,...)
##   #fit$kappa <- 2
##   fit$lambda <- flambda
##   #fit$nu <- 1
##   fit$pars <- c(lambda=fit$lambda)# c(kappa=fit$kappa,lambda=fit$lambda,rho=fit$nu)
##   stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
##   out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
##   out
## }

## #' STOPS version of powermds
## #'
## #' This is power stress with free kappa and lambda but rho is fixed to 1, so no weight transformation.
## #'
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the theta vector of powers; a vector of length 2 where the first element is kappa (for the fitted distances), the second lambda (for the observed proximities). If a scalar is given it is recycled.  Defaults to 1.
## #' @param weightmat (optional) a matrix of nonnegative weights
## #' @param init (optional) initial configuration
## #' @param ndim number of dimensions of the target space
## #' @param itmax number of iterations
## #' @param ... additional arguments to be passed to the fitting procedure
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param structures which structures to look for
## #' @param strucweight weight to be used for the structures; defaults to 0.5
## #' @param strucpars a list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appeacrance in structures 
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param type which weighting to be used in the multi-objective optimization? Either 'additive' (default) or 'multiplicative'. 
## #'
## #' @return A list with the components
## #' \itemize{
## #'         \item stress: the stress
## #'         \item stress.m: default normalized stress
## #'         \item stoploss: the weighted loss value
## #'         \item struc: the structuredness indices
## #'         \item parameters: the parameters used for fitting (kappa, lambda)
## #'         \item fit: the returned object of the fitting procedure
## #'         \item{stopobj:} the stopobj object 
## #' }
## #' 
## #' @import cordillera
## #' @keywords multivariate
## #' @export
## stop_powermds <- function(dis,theta=c(1,1),weightmat=NULL,init=NULL,ndim=2,itmax=100000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
##   theta <- as.numeric(theta)
##   if(inherits(dis,"dist")) dis <- as.matrix(dis) 
##   if(missing(type)) type <- "additive"
##   if(length(theta)>2) stop("There are too many parameters in the theta argument.")
##   if(length(theta)<2) theta <- rep(theta,length.out=2)
##   if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
##   fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=1,weightmat=weightmat,init=init,ndim=ndim,verbose=verbose,itmax=itmax,...)
##   fit$kappa <- theta[1]
##   fit$lambda <- theta[2]
##  # fit$nu <- 1
##   fit$pars <- c(kappa=fit$kappa,lambda=fit$lambda)# c(kappa=fit$kappa,lambda=fit$lambda,rho=fit$nu)
##   stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
##   out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
##   out 
## }

## #' STOPS version of sammon with powers
## #'
## #' This is power stress with free kappa and lambda but rho is fixed to -1 and the weights are delta.
## #' 
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the theta vector of powers; a vector of length two where the first element is kappa (for the fitted distances), the second lambda (for the observed proximities). If a scalar is given it is recycled for the free parameters.  Defaults to 1 1.
## #' @param weightmat (optional) a matrix of nonnegative weights
## #' @param init (optional) initial configuration
## #' @param ndim number of dimensions of the target space
## #' @param itmax number of iterations
## #' @param ... additional arguments to be passed to the fitting procedure
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param structures which structures to look for
## #' @param strucweight weight to be used for the structures; defaults to 0.5
## #' @param strucpars a list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appeacrance in structures 
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param type which weighting to be used in the multi-objective optimization? Either 'additive' (default) or 'multiplicative'. 
## #'
## #' @return A list with the components
## #' \itemize{
## #'         \item stress: the stress
## #'         \item stress.m: default normalized stress
## #'         \item stoploss: the weighted loss value
## #'         \item struc: the structuredness indices
## #'         \item parameters: the parameters used for fitting (kappa, lambda)
## #'         \item fit: the returned object of the fitting procedure
## #'         \item{stopobj:} the stopobj object 
## #' }
## #'
## #' @import cordillera
## #' @keywords multivariate
## #' @export
## stop_powersammon <- function(dis,theta=c(1,1),weightmat=NULL,init=NULL,ndim=2,itmax=100000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
##   theta <- as.numeric(theta)
##   if(inherits(dis,"dist")) dis <- as.matrix(dis)  
##   if(missing(type)) type <- "additive"
##   if(length(theta)>2) stop("There are too many parameters in the theta argument.")
##   if(length(theta)<2) theta <- rep(theta,length.out=2)
##   if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
##   nu <- -1
##   sammwght <-dis^(theta[2])
##   diag(sammwght) <- 1
##   combwght <- sammwght*weightmat
##   fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=nu,weightmat=combwght,init=init,ndim=ndim,verbose=verbose,itmax=itmax,...)
##   fit$kappa <- theta[1]
##   fit$lambda <- theta[2]
##   #fit$rho <- nu
##   fit$pars <- c(kappa=fit$kappa,lambda=fit$lambda)#c(kapp=fit$kappa,lambda=fit$lambda,rho=fit$rho)
##   stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
##   out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
##   out 
## }

## #' STOPS version of elastic scaling with powers for proximities and distances
## #'
## #' This is power stress with free kappa and lambda but rho is fixed to -2 and the weights are delta.
## #'
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the theta vector of powers;  a vector of length two where the first element is kappa (for the fitted distances), the second lambda (for the observed proximities). If a scalar for the free parameters is given it is recycled.  Defaults to 1 1.
## #' @param weightmat (optional) a matrix of nonnegative weights
## #' @param init (optional) initial configuration
## #' @param ndim number of dimensions of the target space
## #' @param itmax number of iterations
## #' @param ... additional arguments to be passed to the fitting procedure
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param structures which strcutures to look for
## #' @param strucweight weight to be used for the structures; defaults to 0.5
## #' @param strucpars a list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appeacrance in structures 
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param type which weighting to be used in the multi-objective optimization? Either 'additive' (default) or 'multiplicative'. 
## #'
## #' @return A list with the components
## #' \itemize{
## #'         \item stress: the stress
## #'         \item stress.m: default normalized stress
## #'         \item stoploss: the weighted loss value
## #'         \item struc: the structuredness indices
## #'         \item parameters: the parameters used for fitting (kappa, lambda)
## #'         \item fit: the returned object of the fitting procedure
## #'         \item{stopobj:} the stopobj object 
## #' }
## #' 
## #' @import cordillera
## #' @keywords multivariate
## #' @export
## stop_powerelastic <- function(dis,theta=c(1,1,-2),weightmat=NULL,init=NULL,ndim=2,itmax=100000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
##   theta <- as.numeric(theta)
##   if(inherits(dis,"dist")) dis <- as.matrix(dis)  
##   if(missing(type)) type <- "additive"
##   if(length(theta)>2) stop("There are too many parameters in the theta argument.")
##   if(length(theta)<2) theta <- rep(theta,length.out=2)
##   if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
##   nu <- -2
##   elawght <- dis^(theta[2])
##   diag(elawght) <- 1
##   combwght <- elawght*weightmat
##   fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=nu,weightmat=combwght,init=init,ndim=ndim,verbose=verbose,itmax=itmax,...)
##   fit$kappa <- theta[1]
##   fit$lambda <- theta[2]
##   #fit$nu <- nu
##   fit$pars <- c(kappa=fit$kappa,lambda=fit$lambda)#c(kappa=fit$kappa,lambda=fit$lambda,rho=fit$nu)
##   stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
##   out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
##   out 
## }


## #' STOPS version of powerstress
## #'
## #' Power stress with free kappa and lambda and rho.
## #'
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the theta vector of powers; the first is kappa (for the fitted distances), the second lambda (for the observed proximities), the third nu (for the weights). If a scalar is given it is recycled.  Defaults to 1 1 1.
## #' @param weightmat (optional) a matrix of nonnegative weights
## #' @param init (optional) initial configuration
## #' @param ndim number of dimensions of the target space
## #' @param itmax number of iterations
## #' @param ... additional arguments to be passed to the fitting procedure
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param structures a character vector listing the structure indices to use. They always are called "cfoo" with foo being the structure.
## #' @param strucweight weight to be used for the structures; defaults to 1/number of structures
## #' @param strucpars a list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appeacrance in structures 
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param type which weighting to be used in the multi-objective optimization? Either 'additive' (default) or 'multiplicative'. 
## #'
## #' @return A list with the components
## #' \itemize{
## #'         \item stress: the stress
## #'         \item stress.m: default normalized stress
## #'         \item stoploss: the weighted loss value
## #'         \item struc: the structuredness indices
## #'         \item parameters: the parameters used for fitting (kappa, lambda, nu)
## #'         \item fit: the returned object of the fitting procedure
## #'         \item{stopobj:} the stopobj object 
## #' }
## #' @keywords multivariate
## #' @export
## stop_powerstress <- function(dis,theta=c(1,1,1),weightmat=NULL,init=NULL,ndim=2,itmax=10000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
##   theta <- as.numeric(theta)
##   if(inherits(dis,"dist")) dis <- as.matrix(dis)
##   if(missing(type)) type <- "additive"
##   if(length(theta)>3) stop("There are too many parameters in the theta argument.")
##   if(length(theta)<3) theta <- rep(theta,length.out=3)
##   if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
##   wght <- weightmat
##   diag(wght) <- 1
##   fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=theta[3],weightmat=wght,init=init,ndim=ndim,verbose=verbose,itmax=itmax,...)
##   fit$kappa <- theta[1]
##   fit$lambda <- theta[2]
##   fit$rho <- theta[3]
##   fit$pars <- c(kappa=fit$kappa,lambda=fit$lambda,rho=fit$rho)
##   stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
##   out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
##   out 
## }


## #' STOPS version of Box Cox Stress
## #'
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the theta vector of powers; the first is mu (for the fitted distances), the second lambda (for the  proximities), the third nu (for the weights). If a scalar is given it is recycled.  Defaults to 1 1 0.
## #' @param weightmat (not used) 
## #' @param init (optional) initial configuration
## #' @param ndim number of dimensions of the target space
## #' @param itmax number of iterations
## #' @param ... additional arguments to be passed to the fitting procedure
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param structures which structures to look for
## #' @param strucweight weight to be used for the structures; defaults to 0.5
## #' @param strucpars a list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appeacrance in structures 
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param type which weighting to be used in the multi-objective optimization? Either 'additive' (default) or 'multiplicative'. 
## #'
## #' @return A list with the components
## #' \itemize{
## #'         \item stress: the stress
## #'         \item stress.m: default normalized stress
## #'         \item stoploss: the weighted loss value
## #'         \item struc: the structuredness indices
## #'         \item parameters: the parameters used for fitting (kappa, lambda)
## #'         \item fit: the returned object of the fitting procedure
## #'          \item{stopobj:} the stopobj object 
## #' }
## #' @keywords multivariate
## #' @export
## stop_bcstress <- function(dis,theta=c(1,1,0),weightmat=NULL,init=NULL,ndim=2,itmax=5000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
##   theta <- as.numeric(theta)
##   if(inherits(dis,"dist")) dis <- as.matrix(dis)
##   if(missing(type)) type <- "additive"
##   if(length(theta)>3) stop("There are too many parameters in the theta argument.")
##   if(length(theta)<3) theta <- rep(theta,length.out=3)
##   #if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
##   #wght <- weightmat
##   #diag(wght) <- 1
##   fit <- bcStressMin(delta=dis,mu=theta[1],lambda=theta[2],rho=theta[3],init=init,ndim=ndim,verbose=verbose+2,itmax=itmax,...)
##   fit$mu <- theta[1]
##   fit$lambda <- theta[2]
##   fit$rho <- theta[3]
##   fit$pars <- c(mu=fit$mu,lambda=fit$lambda,rho=fit$rho)
##   stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
##   out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
##   out 
## }

## #' STOPS version of lMDS
## #'
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the theta vector of powers; the first is k (for the neighbourhood), the second tau (for the penalty) . If a scalar is given it is recycled.  Defaults to 2 and 0.5.
## #' @param weightmat (not used) 
## #' @param init (optional) initial configuration
## #' @param ndim number of dimensions of the target space
## #' @param itmax number of iterations 
## #' @param ... additional arguments to be passed to the fitting procedure
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param structures which structures to look for
## #' @param strucweight weight to be used for the structures; defaults to 0.5
## #' @param strucpars a list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appeacrance in structures 
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param type which weighting to be used in the multi-objective optimization? Either 'additive' (default) or 'multiplicative'. 
## #'
## #' @return A list with the components
## #' \itemize{
## #'         \item stress: the stress
## #'         \item stress.m: default normalized stress
## #'         \item stoploss: the weighted loss value
## #'         \item struc: the structuredness indices
## #'         \item parameters: the parameters used for fitting (kappa, lambda)
## #'         \item fit: the returned object of the fitting procedure
## #'         \item{stopobj:} the stopobj object 
## #' }
## #' @keywords multivariate
## #' @export
## stop_lmds <- function(dis,theta=c(2,0.5),weightmat=NULL,init=NULL,ndim=2,itmax=5000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
##   theta <- as.numeric(theta)
##   if(inherits(dis,"dist")) dis <- as.matrix(dis)
##   if(missing(type)) type <- "additive"
##   #we allow for three parameters in the theta argument
##   if(length(theta)>2) stop("There are too many parameters in the theta argument.")
##   if(length(theta)<2) theta <- rep(theta,length.out=2)
##   #if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
##   #wght <- weightmat
##   #diag(wght) <- 1
##   fit <- lmds(delta=dis,k=theta[1],tau=theta[2],init=init,ndim=ndim,verbose=verbose+2,itmax=itmax,...)
##   fit$k <- theta[1]
##   fit$tau <- theta[2]
##   fit$pars  <- c(k=fit$k,tau=fit$tau)
##   stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
##   out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
##   out 
## }


## #' STOPS version of restricted powerstress

## # This is a power stress where kappa and lambda are free to vary but restricted to be equal, so the same exponent will be used for distances and dissimilarities. rho (for the weights) is also free.  
## #'
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the theta vector of powers; the first two arguments are for kappa and lambda and should be equal (for the fitted distances and observed proximities), the third nu (for the weights). Internally the kappa and lambda are equated. If a scalar is given it is recycled (so all elements of theta are equal); if a vector of length 2 is given, it gets expanded to c(theta[1],theta[1],theta[2]). Defaults to 1 1 1.
## #' @param weightmat (optional) a matrix of nonnegative weights
## #' @param init (optional) initial configuration
## #' @param ndim number of dimensions of the target space
## #' @param itmax number of iterations. default is 10000.
## #' @param ... additional arguments to be passed to the fitting procedure powerStressMin
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param structures a character vector listing the structure indices to use. They always are called "cfoo" with foo being the structure.
## #' @param strucweight weight to be used for the structures; defaults to 1/number of structures
## #' @param strucpars a list of list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appearance in structures vector. See examples.  
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param type which weighting to be used in the multi-objective optimization? Either 'additive' (default) or 'multiplicative'.
## #' 
## #' @return A list with the components
## #' \itemize{
## #'         \item stress: the stress
## #'         \item stress.m: default normalized stress
## #'         \item stoploss: the weighted loss value
## #'         \item struc: the structuredness indices
## #'         \item parameters: the parameters used for fitting (kappa=lambda, nu)
## #'         \item fit: the returned object of the fitting procedure
## #'          \item{stopobj:} the stopobj object 
## #' }
## #' @keywords multivariate
## stop_rpowerstress <- function(dis,theta=c(1,1,1),weightmat=NULL,init=NULL,ndim=2,itmax=10000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
##   theta <- as.numeric(theta)
##   if(inherits(dis,"dist")) dis <- as.matrix(dis)
##   if(missing(type)) type <- "additive"
##   if(length(theta)>3) stop("There are too many parameters in the theta argument.")
##   #if(length(theta)==3L & theta[1]!=theta[2]) warning("The powers given for kappa and lambda do not agree. The first value in theta will be used for both kappa and lambda.")  
##   if(length(theta)==1L) theta <- rep(theta,3)
##   if(length(theta)==2L) theta <- c(rep(theta[1],2),theta[2])
##   if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
##   wght <- weightmat
##   diag(wght) <- 1
##   fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[1],nu=theta[3],weightmat=wght,init=init,ndim=ndim,verbose=verbose,itmax=itmax,...)
##   #fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=theta[3],weightmat=wght,ndim=ndim,verbose=verbose,itmax=itmaxi)
##   #if(stresstype=="default") fit$stress.m <- fit$stress.m
##   #if(stresstype=="stress1") fit$stress.m <- fit$stress.1
##   #if(stresstype=="rawstress") fit$stress.m <- fit$stress.r
##   #if(stresstype=="normstress") fit$stress.m <- fit$stress.n
##   #if(stresstype=="enormstress") fit$stress.m <- fit$stress.en
##   #if(stresstype=="enormstress1") fit$stress.m <- fit$stress.en1
##   fit$kappa <- theta[1]
##   fit$lambda <- theta[1]
##   fit$rho <- theta[3]
##   fit$pars <- c(kappa=fit$kappa,lambda=fit$lambda,rho=fit$rho)
##   stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
##   out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
##   out 
## }

## #' STOPS version of approximated power stress models.
## #'
## #' This uses an approximation to power stress that can make use of smacof as workhorse. Free parameters are tau and upsilon.
## #' 
## #' @param dis numeric matrix or dist object of a matrix of proximities
## #' @param theta the theta vector of parameters to optimize over. Must be of length two, with the first the tau argument and the second the upsilon argument. It can also be a scalar of the tau and upsilon transformation for the observed proximities and gets recycled for both ups and tau (so they are equal). Defaults to 1 1.  
## #' @param ndim number of dimensions of the target space
## #' @param itmax number of iterations. default is 1000.
## #' @param weightmat (optional) a binary matrix of nonnegative weights
## #' @param init (optional) initial configuration
## #' @param ... additional arguments to be passed to the fitting procedure
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param structures a character vector listing the structure indices to use. They always are called "cfoo" with foo being the structure.
## #' @param strucweight weight to be used for the structures; defaults to 1/number of structures
## #' @param strucpars a list of list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appearance in structures vector. See examples.  
## #' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
## #' @param type which weighting to be used in the multi-objective optimization? Either 'additive' (default) or 'multiplicative'.
## #' 
## #' @return A list with the components
## #'    \itemize{
## #'         \item{stress:} the stress 1 (sqrt stress.m)
## #'         \item{stress.m:} default normalized stress
## #'          \item stoploss: the weighted loss value
## #'         \item struc: the structuredness indices
## #'         \item parameters: the parameters used for fitting (kappa=1, tau, ups)
## #'         \item fit: the returned object of the fitting procedure
## #'         \item{stopobj:} the stopobj object 
## #' }
## #'
## #'@importFrom stats dist as.dist
## #'@import smacof
## #' 
## #'@keywords multivariate
## stop_apstress <- function(dis,theta=c(1,1),ndim=2,weightmat=NULL,init=NULL,itmax=1000,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
##                                         #TODO Unfolding
##   if(is.null(init)) init <- "torgerson"
##   if(inherits(dis,"dist")) dis <- as.matrix(dis)
##   if(missing(type)) type <- "additive"
##   if(is.null(weightmat)) weightmat <- 1-diag(ncol(dis))
##   if(length(setdiff(unique(unlist(as.vector(weightmat))),c(0,1)))>0) stop("For approximated power stress, only binary weight matrices are allowed.")  
##   #we allow for theta to be of length three for compatibility in stops; maybe change that in the future 
##   if(length(theta)>2) stop("There are too many parameters in the theta argument.")
##   if(length(theta)==1L) theta <- rep(theta,2)
##   #if(length(theta)==2L) theta <- c(1,theta) 
##   tau <- theta[1]
##   ups <- theta[2]
##   combwght <- weightmat*(dis^ups)
##   fit <- smacof::smacofSym(dis^tau,ndim=ndim,weightmat=combwght,init=init,verbose=isTRUE(verbose==2),itmax=itmax,...) #optimize with smacof
##   #fit$kappa <- 1
##   fit$tau <- tau
##   fit$upsilon <- ups
##   fit$stress.1 <- fit$stress #smacof stress is sqrt(stress.m); for compatibility with powerStressMin we use stress^2 as stress.m 
##   fitdis <- as.matrix(fit$confdist)
##   delts <- as.matrix(fit$delta) 
##   fit$stress.r <- sum(combwght*(delts-fitdis)^2)
##   fit$stress.m <- fit$stress^2 #This is for compatibility with powerStressMin
##   fit$pars <- c(tau=fit$tau,upsilon=fit$upsilon) #c(kappa=fit$kappa,tau=fit$tau,upsilon=fit$upsilon)
##   fit$deltaorig <- fit$delta^(1/fit$tau)
##   stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
##   out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
##   out 
## }

## #' MakePower
## #'
## #' @param x matrix
## #' @param theta numeric (power)
## mkPower2<-function(x,theta) {
##     r <- theta
##     if(length(theta) > 1) r <- theta[2] 
##     n<-nrow(x)
##     return(abs((x+diag(n))^r)-diag(n))
## }

#' High Level STOPS Function
#'
#' @description This allows to fit STOPS models as described in Rusch, Mair, Hornik (2023). 
#'
#' @details The combination of c-structurednes indices and stress uses the stress.m values, which are the explictly normalized stresses. Reported however is the stress-1 value which is sqrt(stress.m). 
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param loss which loss function to be used for fitting, defaults to stress. 
#' @param theta hyperparameter vector starting values for the transformation functions. If the length is smaller than the number of hyperparameters for the MDS version the vector gets recycled (see the corresponding stop_XXX function or the vignette for how theta must look like exactly for each loss). If larger than the number of hyperparameters for the MDS method, an error is thrown. If completely missing theta is set to 1 and recycled.      
#' @param structures character vector of which c-structuredness indices should be considered; if missing no structure is considered.
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights; defaults to 1 for all off diagonals 
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param strucweight vector of weights to be used for the c-structuredness indices (in the same order as in structures); defaults to -1/length(structures) for each index 
#' @param strucpars (possibly named with the structure). Metaparameters for the structuredness indices (gamma in the article). It's safest for it be a list of lists with the named arguments for the structuredness indices and the order of the lists must be like the order of structures. So something like this \code{list(list(par1Struc1=par1Struc1,par2Struc1=par2Struc1),list(par1Struc2=par1Struc2,par2Struc2=par2Struc2),...)} where parYStrucX are the named arguments for the metaparameter Y of the structure X the list elements corresponds to. For a structure without parameters, set NULL. Parameters in different list elements parYStrucX can have the same name. For example, say we want to use cclusteredness with metaparameters epsilon=10 and k=4 (and the default for the other parameters), cdependence with no metaparameters and cfaithfulness with metaparameter k=7 one would \code{list(list(epsilon=10,k=4),list(NULL),list(dis=obdiss,k=6))}  for structures vector ("cclusteredness","cdependence","cfaithfulness"). The parameter lists must be in the same ordering as the indices in structures. If missing it is set to NULL and defaults are used. It is also possible to supply a structure's metaparameters as a list of vectors with named elements if the metaparameters are scalars, so like \code{list(c(par1Struc1=parStruc1,par2Struc1=par1Struc1,...),c(par1Struc2=par1Struc2,par2Struc2=par2Struc2,...))}. That can have unintended consequences if the metaparameter is a vector or matrix.  
#' @param optimmethod What solver to use. Currently supported are Bayesian optimization with Gaussian Process priors and Kriging ("Kriging"), Bayesian optimization with treed Gaussian processes with jump to linear models ("tgp"), Adaptive LJ Search ("ALJ"), Particle Swarm optimization ("pso"), simulated annealing ("SANN"), "DIRECT", Stochastic Global Optimization ("stogo"), COBYLA ("cobyla"), Controlled Random Search 2 with local mutation ("crs2lm"), Improved Stochastic Ranking Evolution Strategy ("isres"), Multi-Level Single-Linkage ("mlsl"), Nelder-Mead ("neldermead"), Subplex ("sbplx"), Hooke-Jeeves Pattern Search ("hjk"), CMA-ES ("cmaes"). Defaults to "ALJ" version. tgp, ALJ, Kriging and pso usually work well for relatively low values of itmax. 
#' @param lower The lower contraints of the search region. Needs to be a numeric vector of the same length as the parameter vector theta. 
#' @param upper The upper contraints of the search region. Needs to be a numeric vector of the same length as the parameter vector theta.  
#' @param verbose numeric value hat prints information on the fitting process; >2 is very verbose.
#' @param type which aggregation for the multi objective target function? Either 'additive' (default) or 'multiplicative'
#' @param itmax maximum number of iterations of the outer optimization (for theta) or number of steps of Bayesian optimization; default is 50. We recommend a higher number for ALJ (around 150). Note that due to the inner workings of some solvers, this may or may not correspond to the actual number of function evaluations performed (or PS models fitted). E.g., with tgp the actual number of function evaluation of the PS method is between itmax and 6*itmax as tgp samples 1-6 candidates from the posterior and uses the best candidate. For pso it is the number of particles s times itmax. For cmaes it is usually a bit higher than itmax. This currently may get overruled by a control argument if it is used (and then set to either ewhat is supplie dby control or to the default of the method).    
#' @param itmaxps maximum number of iterations of the inner optimization (to obtain the PS configuration)
#' @param initpoints number of initial points to fit the surrogate model for Bayesian optimization; default is 10.
#' @param model a character specifying the surrogate model to use. For Kriging it specifies the covariance kernel for the GP prior; see \code{\link{covTensorProduct-class}} defaults to "powerexp". For tgp it specifies the non stationary process used see \code{\link{bgp}}, defaults to "btgpllm" 
#' @param control a control argument passed to the outer optimization procedure. Will override any other control arguents passed, especially verbose and itmax. For the efect of control, see the functions pomp::sannbox for SANN and pso::psoptim for pso, cmaes::cma_es for cmaes, dfoptim::hjkb for hjk and the nloptr docs for the algorithms DIRECT, stogo, cobyla, crs2lm, isres, mlsl, neldermead, sbplx.
#' @param ... additional arguments passed to the outer optimization procedures (not fully tested).
#
#'@return A list with the components
#'         \itemize{
#'         \item stoploss: the stoploss value
#'         \item optim: the object returned from the optimization procedure
#'         \item stressweight: the stressweight
#'         \item strucweight: the vector of structure weights
#'         \item call: the call
#'         \item optimmethod: The solver selected
#'         \item losstype: The PS badness-of-fit function
#'         \item nobj: the number of objects in the configuration
#'         \item type: The type of stoploss scalacrisation (additive or multiplicative)
#'         \item fit: The fitted PS object (most importantly $fit$conf the fitted configuration) 
#' }
#' 
#' @examples
#'
#' data(kinshipdelta,package="smacof")
#' strucpar<-list(NULL,NULL) #parameters for indices
#' res1<-stops(kinshipdelta,loss="stress",
#' structures=c("cclumpiness","cassociation"),strucpars=strucpar,
#' lower=0,upper=10,itmax=10)
#' res1
#' 
#' \donttest{
#' data(BankingCrisesDistances)
#' strucpar<-list(c(epsilon=10,minpts=2),NULL) #parameters for indices
#' res1<-stops(BankingCrisesDistances[,1:69],loss="stress",verbose=0,
#' structures=c("cclusteredness","clinearity"),strucpars=strucpar,
#' lower=0,upper=10)
#' res1
#'
#' strucpar<-list(list(alpha=0.6,C=15,var.thr=1e-5,zeta=NULL),
#' list(alpha=0.6,C=15,var.thr=1e-5,zeta=NULL))
#' res1<-stops(BankingCrisesDistances[,1:69],loss="stress",verbose=0,
#' structures=c("cfunctionality","ccomplexity"),strucpars=strucpar,
#' lower=0,upper=10)
#' res1
#' }
#' 
#' @importFrom stats dist as.dist optim
#' @importFrom utils capture.output
#' @importFrom pso psoptim
#' @importFrom DiceOptim EGO.nsteps
#' @importFrom DiceKriging km
#' @importFrom tgp lhs dopt.gp
#' @importFrom pomp sannbox
#' @importFrom nloptr direct stogo cobyla crs2lm isres mlsl neldermead sbplx
#' @importFrom cmaes cma_es
#' @importFrom dfoptim hjkb
#' @import cordillera
#' 
#' @keywords clustering multivariate
#' @export
stops <- function(dis,loss=c("strain","stress","smacofSym","powerstress","powermds","powerelastic","powerstrain","elastic","sammon","sammon2","smacofSphere","powersammon","rstress","sstress","isomap","isomapeps","bcstress","lmds","apstress","rpowerstress"), theta=1, structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"), ndim=2, weightmat=NULL, init=NULL, stressweight=1, strucweight, strucpars, optimmethod=c("SANN","ALJ","pso","Kriging","tgp","DIRECT","stogo","cobyla","crs2lm","isres","mlsl","neldermead","sbplx","hjk","cmaes"), lower, upper, verbose=0, type=c("additive","multiplicative"), initpoints=10, itmax=50,itmaxps=10000, model, control,...)
    {
      #TODO add more transformations for the g() and f() by the transformation argument. We only use power versions right now, flexsmacof will allow for more (splines or a smoother or so)
      if(missing(structures)) {
          structures <- "clinearity"
          strucweight <- 0
      }
      if(missing(strucpars)) strucpars <- vector("list",length(structures))
      if(inherits(dis,"dist")) dis <- as.matrix(dis)
      if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
      if(missing(loss)) loss <- "stress"
      if(missing(type)) type <- "additive"
      #TODO implement a Pareto multiobjective
      .confin <- init #initialize a configuration
      psfunc <- switch(loss, "powerstrain"=stop_cmdscale, "stress"=stop_smacofSym,"smacofSym"=stop_smacofSym,"powerstress"=stop_powerstress,"strain"=stop_cmdscale,"smacofSphere"=stop_smacofSphere,"rstress"=stop_rstress,"sammon"=stop_sammon, "elastic"=stop_elastic, "powermds"=stop_powermds,"powerelastic"=stop_powerelastic,"powersammon"=stop_powersammon,"sammon2"=stop_sammon2,"sstress"=stop_sstress,"isomap"=stop_isomap1,"isomapeps"=stop_isomap2,"bcstress"=stop_bcmds,"bcmds"=stop_bcmds,"lmds"=stop_lmds,"apstress"=stop_apstress,"rpowerstress"=stop_rpowerstress) #choose the stress to minimize
      if(missing(strucweight)) {
         #TODO: automatic handler of setting weights that makes sense
         strucweight <- rep(-1/length(structures),length(structures))
         if(verbose>1) cat("Weights are stressweight=",stressweight,"strucweights=", strucweight,"\n")
      }
      if(missing(optimmethod)) optimmethod <- "ALJ"  
      if(verbose>0) cat("Starting Optimization \n ")
        if(optimmethod=="SANN") {
            ##NOTE: Used optim first, but that has no box constraints. Then tried all kinds but it all sucks.
            ## opt<- stats::optim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmax=itmaxps))$stoploss,method="SANN",...)
            ## GenSA::GenSA
            ## GenSA but can't pass-by-value properly with objects, so needs to be set manuall. Also has other issues with the max.call not always working (not sure why). SUCKS!
           ## if(missing(control)) control <- list(verbose=!isTRUE(verbose==0),smooth=FALSE,max.call=120)
           ##TR: somethings off here, doesn't take the control list from above. I know why:> doesn't allows object passing
           ##control <- list(verbose=TRUE,max.call=10,maxit=10,smooth=FALSE) #takes this one,difference is max.call-itmax vs a real value   
           ## opt <- GenSA::GenSA(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmax=itmaxps))$stoploss,lower=lower,upper=upper,control=list(verbose=TRUE,smooth=FALSE,max.call=9,maxit=9))
            ## optimization::optim_sa
            ## Disadvantage: Can't control the maximum number of iterations
          ##if(missing(control)) control <- list(nlimit=itmax)
          ##control <- list(nlimit=5)
            ##opt <- optimization::optim_sa(start=theta,fun=function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmax=itmaxps))$stoploss,trace=TRUE,lower=lower,upper=upper,control=control)
          ##thetaopt <- opt$par #for optim and GenSA
          ##bestval <- opt$value #for optim and GenSA
          ##bestval <- opt$function_value#for optim_sa
          ##itel <- opt$counts#for optim and GenSA
          ## itel <- dim(opt$trace)[1]#for optim_sa
          ##pomp::sannbox
            if(missing(control)) control <- list(trace=verbose-2,lower=lower,upper=upper,maxit=itmax)
            opt <- pomp::sannbox(par=theta,fn=function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps))$stoploss,control=control)
          thetaopt <- opt$par
          bestval <- opt$value
          itel <- opt$counts[1]
      }
       if(optimmethod=="pso") {
        #addargs <- list(...)
        if(missing(control)) control <- list(trace=verbose-2,s=5,maxit=itmax)
        opt<- pso::psoptim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps))$stoploss,lower=lower,upper=upper,control=control,...)
         thetaopt <- opt$par
         bestval <-  opt$value
         itel <- opt$counts["function"]
       }
      if(optimmethod=="ALJ")  {
        opt <- ljoptim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps))$stoploss,lower=lower,upper=upper,verbose=verbose-2,itmax=itmax,...)
       thetaopt <- opt$par
       bestval <-  opt$value
       itel <- opt$counts["function"] 
    }    
    if(optimmethod=="Kriging")
    {
        if(missing(model)) model <- "powexp"
       # optdim <- 3 #dimensions
       # if(loss%in%c("powerstrain","stress","smacofSym","smacofSphere","strain","sammon","elastic","sammon2","sstress","rstress")) optdim <- 1
       # if(loss%in%c("powermds","powerelastic","powersammon","smacofSphere","strain","sammon","elastic","sammon2")) optdim <- 2
        recto <- cbind(lower,upper)
        Xcand <- tgp::lhs(initpoints*100,recto)
        if(missing(theta)) theta <- Xcand[1,]
        x <- t(theta)
        X <- tgp::dopt.gp(initpoints-1,X=x,Xcand)$XX
        design <- rbind(x,X)
        #design <- data.frame(X) 
        responsec <- apply(design, 1, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps))$stoploss) #support points for fitting kriging model
        if (verbose>1) cat("Kriging Model Fitting","\n")
        surrogatemodel <- DiceKriging::km(~1, design = design, response = responsec,covtype=model,control=list(trace=isTRUE(verbose>3))) #fit the kriging model
        #EGO.nsteps has no verbose argument so I capture.output and return it if desired
        if (verbose>2) cat("EGO (DICE) Optimization","\n")
        logged <- capture.output({
           opt<- DiceOptim::EGO.nsteps(model=surrogatemodel, fun=function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps))$stoploss,lower=lower,upper=upper,nsteps=itmax,...)
       }) #bayesian optimization with gaussian process prior
       if(verbose>2) print(logged)
       thetaopt <- opt$par[which.min(opt$value),] #parameters where best value found (we do not use the last one as that may be worse)
       bestval <- min(opt$value) #best stoploss value
       itel <- itmax
       }
  if(optimmethod=="tgp")
    {
        #if(!isNamespaceLoaded("tgp")) attachNamespace("tgp")
        if(missing(model)) model <- "btgpllm"
        #model <- get(model,envir=getNamespace("tgp"))
        #if(loss%in%c("powerstrain","stress","smacofSym","smacofSphere","strain","sammon","elastic","sammon2","sstress","rstress")) optdim <- 1
        #if(loss%in%c("powermds","powerelastic","powersammon","smacofSphere","strain","sammon","elastic","sammon2")) optdim <- 2
        if (verbose>1) cat("EGO (TGP) Optimization","\n")
        opt <- tgpoptim(theta, fun=function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps))$stoploss,lower=lower,upper=upper,itmax=itmax,initpoints=initpoints,model=model,verbose=verbose-2,...) #bayesian optimization with treed gaussian process prior
       thetaopt <- opt$par #parameters where best value found (we do not use the last one as that may be worse)
       bestval <- opt$value #best stoploss value
       itel <- opt$counts["function"] 
    }
     if(optimmethod=="DIRECT") {
            if (verbose>1) cat("DIRECT Optimization","\n")
            if(missing(control)) control <- list(maxeval=itmax)
          opt<- nloptr::direct(function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps))$stoploss,lower=lower,upper=upper,nl.info=isTRUE(all.equal(verbose-2,0)),control=control,...)
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$iter
     }
      if(optimmethod=="stogo") {
          if (verbose>1) cat("StoGO Optimization","\n")
          #cat(itmax,"\n")
         #   if(missing(control)) control <- list(maxeval=itmax)  
          opt<- nloptr::stogo(theta,function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps))$stoploss,lower=lower,upper=upper,nl.info=isTRUE(verbose>2),maxeval=itmax,...)
          #TODO Issue with maxeval?
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$iter
      }
         if(optimmethod=="cobyla") {
            if (verbose>1) cat("COBYLA Optimization","\n")
            if(missing(control)) control <- list(maxeval=itmax)  
           opt<- nloptr::cobyla(theta,function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps))$stoploss,lower=lower,upper=upper,nl.info=isTRUE(all.equal(verbose-2,0)),control=control,...)
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$iter
         }
        if(optimmethod=="crs2lm") {
            if (verbose>1) cat("crs2lm Optimization","\n")
            #if(missing(control)) control <- list(maxeval=itmax)  
           opt<- nloptr::crs2lm(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps))$stoploss,lower=lower,upper=upper,nl.info=isTRUE(verbose>3),maxeval=itmax,...)
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$iter
        }
        if(optimmethod=="isres") {
            if (verbose>1) cat("isres Optimization","\n")
            #if(missing(control)) control <- list(maxeval=itmax)  
           opt<- nloptr::crs2lm(theta,function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps))$stoploss,lower=lower,upper=upper,nl.info=isTRUE(verbose>3),maxeval=itmax,...)
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$iter
        }
        if(optimmethod=="mlsl") {
            if (verbose>1) cat("MLSL Optimization","\n")
            if(missing(control)) control <- list(maxeval=itmax)  
           opt<- nloptr::mlsl(theta,function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps))$stoploss,lower=lower,upper=upper,nl.info=isTRUE(verbose>3),control=control,...)
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$iter
        }
        if(optimmethod=="neldermead") {
            if (verbose>1) cat("Nelder-Mead Optimization","\n")
            if(missing(control)) control <- list(maxeval=itmax)  
           opt<- nloptr::neldermead(theta,function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps))$stoploss,lower=lower,upper=upper,nl.info=isTRUE(verbose>3),control=control,...)
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$iter
        }
          if(optimmethod=="sbplx") {
            if (verbose>1) cat("Subplex Optimization","\n")
            if(missing(control)) control <- list(maxeval=itmax)  
           opt<- nloptr::sbplx(theta,function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps))$stoploss,lower=lower,upper=upper,nl.info=isTRUE(verbose>3),control=control,...)
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$iter
      }
        if(optimmethod=="hjk") {
            if (verbose>1) cat("Hooke-Jeeves Optimization","\n")
            if(missing(control)) control <- list(info=isTRUE(all.equal(verbose-2,0)),maxfeval=itmax)
          opt<- dfoptim::hjkb(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps))$stoploss,lower=lower,upper=upper,control=control,...)
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$feval
        }
         if(optimmethod=="cmaes") {
            if (verbose>1) cat("CMA-ES Optimization","\n")
            theta <- as.vector(theta)
            if(missing(control)) {
                  #for calculating the correct number of calls to fn
                  N <- length(theta)
                  lambda <- 4 + floor(3 * log(N))
                  control <- list(maxit=ceiling(itmax/lambda))
                }
          opt<- cmaes::cma_es(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps))$stoploss,lower=lower,upper=upper,control=control,...)
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$counts[1]
         }
        #TODO: Streamline the number of function evaluations to be supplied returned for the different solvers. E.g. for cma_es it is itmax*population size. Or for tgp it is also itmax*6 or so. 
    #refit optimal model  
    out <- do.call(psfunc,list(dis=dis,theta=thetaopt,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps))
    out$stoploss <- bestval
    out$theta <- out$parameters
    out$optim <- opt
    out$stressweight <- stressweight
    out$strucweight <- strucweight
    out$call <- match.call()
    out$optimethod <- optimmethod
    out$losstype <- loss
    out$nobj <- dim(out$fit$conf)[1]
    out$type <- type
    if(verbose>1) cat("Found minimum after",itel," iterations at",round(thetaopt,4),"with stoploss=",round(out$stoploss,4),"and default scaling loss=",round(out$stress.m,4),"and c-structuredness indices:",t(data.frame(names(out$strucindices),out$strucindices)),". Thanks for your patience. \n")
    class(out) <- c("stops")
    out
  }

#' S3 print method for stops objects
#'
#'@param x stops object
#'@param ... additional arguments
#'@export
#'@return no return value, just prints
print.stops <- function(x,...)
    {
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    cat("Model:",x$type,"STOPS with", x$loss,"loss function and theta parameter vector",paste("(",paste(attr(x$parameters,"names"),collapse=" "),")",sep="")," = ",x$parameters,"\n")
    cat("\n")
    cat("Number of objects:", x$nobj, "\n")
    cat("MDS loss value:", x$stress.m, "\n")
    cat("C-Structuredness Indices:", t(data.frame(names(x$strucindices),x$strucindices)),"\n")
    cat("Structure optimized loss (stoploss):", x$stoploss, "\n")
    cat("MDS loss weight:",x$stressweight,"c-structuredness weights:",x$strucweight,"\n")
    cat("Number of iterations of",x$optimethod,"optimization:", x$optim$counts["function"], "\n")
    cat("\n")
    }

#' S3 coef method for stops objects
#'
#'@param object object of class stops 
#'@param ... addditional arguments 
#'@export
#'@importFrom stats coef
#'@return a vector of hyperparmeters theta 
coef.stops <- function(object,...)
    {
    return(c(object$par))
    }


#'S3 plot method for stops objects
#' 
#'@param x an object of class stops
#'@param plot.type String indicating which type of plot to be produced: "confplot", "resplot", "Shepard", "stressplot", "bubbleplot" (see details)
#'@param main the main title of the plot
#'@param asp aspect ratio of x/y axis; defaults to NA; setting to 1 will lead to an accurate represenation of the fitted distances. 
#'@param ... Further plot arguments passed: see 'plot.smacof' and 'plot' for detailed information.
#' 
#'Details:
#' \itemize{
#' \item Configuration plot (plot.type = "confplot"): Plots the MDS configurations.
#' \item Residual plot (plot.type = "resplot"): Plots the dissimilarities against the fitted distances.
#' \item Linearized Shepard diagram (plot.type = "Shepard"): Diagram with the transformed observed dissimilarities against the transformed fitted distance as well as loess smooth and a least squares line.
#' \item Stress decomposition plot (plot.type = "stressplot", only for SMACOF objects in $fit): Plots the stress contribution in of each observation. Note that it rescales the stress-per-point (SPP) from the corresponding smacof function to percentages (sum is 100). The higher the contribution, the worse the fit.
#' \item Bubble plot (plot.type = "bubbleplot", only available for SMACOF objects $fit): Combines the configuration plot with the point stress contribution. The larger the bubbles, the better the fit.
#'}
#'
#'@return no return value, just plots
#' 
#'@importFrom graphics plot 
#'@export 
plot.stops <- function(x,plot.type=c("confplot"), main, asp=NA,...)
    {
     if(missing(plot.type)) plot.type <- "confplot"  
      plot(x$fit,plot.type=plot.type,main=main,asp=asp,...)
 }

#' S3 plot3d method for class stops
#'
#' 
#' This methods produces a dynamic 3D configuration plot.
#' @param x object of class stops
#' @param ... Further plot arguments to the method of the class of slot $fit, see \code{\link{plot.smacof}} or \code{\link{plot3d.cmdscaleE}} . Also see 'rgl' in package 'rgl' 
#'
#' @return no return value, just plots
#' 
#' @export
#' @import rgl
plot3d.stops <- function(x,...)
    {
        plot3d(x$fit,...)
    }

#' S3 plot3dstatic method for class stops
#' 
#' This methods produces a static 3D configuration plot.
#' @param x object of class stops
#' @param ... Further plot arguments to the method of the class of slot fit, see \code{\link{plot3dstatic}} or \code{\link{plot3dstatic.cmdscaleE}} . Also see 'scatterplot3d' in package 'scatterplot3d'.
#'
#' @return no return value, just plots
#' 
#' @export
plot3dstatic.stops <- function(x,...)
    {
        plot3dstatic(x$fit,...)
    }

#' S3 residuals method for stops
#'@param object object of class stops
#'@param ... addditional arguments
#'@importFrom stats residuals
#'@export
#'@return a vector of residuals (observed minus fitted distances) 
residuals.stops <- function(object,...)
    {
    stats::residuals(object$fit,...)
    }


#' S3 summary method for stops
#'
#'@param object object of class stops
#'@param ... addditional arguments
#' 
#'@export
#'@return object of class 'summary.stops'
summary.stops <- function(object,...)
    {
      sppmat <- NULL
      if(!is.null(object$fit$spp))
      { 
           spp.perc <- object$fit$spp/sum(object$fit$spp) * 100
           sppmat <- cbind(sort(object$fit$spp), sort(spp.perc))
           colnames(sppmat) <- c("SPP", "SPP(%)")
      } 
      res <- list(conf=object$fit$conf,sppmat=sppmat)
      class(res) <- "summary.stops"
      res
    }

#' S3 print method for summary.stops
#' 
#'@param x object of class summary.stops
#'@param ... additional arguments
#'@export
#'@return no return value, just prints 
print.summary.stops <- function(x,...)
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

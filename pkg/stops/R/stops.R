#What follows are proof of concepts for the STOPS paper 

#Idea for stops function allow an arbitrary number of indices in a weighted multi-objective optimization way; for this use stoplose
# write stops_foo where foo is the MDS model of interest
# also do this with a pareto approach
    
#'  Calculate the weighted multiobjective loss function used in STOPS
#'
#' @param obj MDS object (supported are stop_sammon, stop_cmdscale, stop_smacofSym, stop_rstress, stop_powerstress, stop_smacofSphere) 
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param structures which c-structuredness indices to be included in the loss
#' @param strucweight the weights of the structuredness indices; defaults to -1/#number of structures
#' @param strucpars a list of parameters to be passed to the c-structuredness indices in the same order as the values in structures #(alternatively a named list that has the structure name as the element name)
#' @param type what type of weighted optimization should be used? Can be 'additive' or 'multiplicative'. NOte that for penalizing the mds loss. 
#' @param verbose verbose output
#'
#' @import cordillera
#' 
#' @export
stoploss<- function(obj,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness"),strucweight=rep(-1/length(structures),length(structures)),strucpars,type=c("additive","multiplicative"),verbose=0)
    {
        if(missing(strucpars)) strucpars <- vector("list", length(structures))
        stressi <- obj$stress.m
        pars <- obj$pars
        confs <- obj$conf 
        if("cclusteredness"%in%structures)
            {
              indst <- which(structures=="cclusteredness")  
              cclusteredness <- do.call(cordillera::cordillera,c(list(confs),strucpars[[indst]]))$normed
            }                           
        if("clinearity"%in%structures)
            {
               indst <- which(structures=="clinearity")
               clinearity <- do.call(stops::c_linearity,list(confs))
           }
        if("cdependence"%in%structures)
            {
               indst <- which(structures=="cdependence")
               cdependence <- do.call(stops::c_dependence,c(list(confs),strucpars[[indst]])) 
           }
        if("cmanifoldness"%in%structures)
            {
               indst <- which(structures=="cmanifoldness")
               cmanifoldness <- do.call(stops::c_manifoldness,c(list(confs)))
           }
        if("cassociation"%in%structures)
            {
               indst <- which(structures=="cassociation")
               cassociation <- do.call(stops::c_association,c(list(confs),strucpars[[indst]]))
           }
        if("cnonmonotonicity"%in%structures)
            {
               indst <- which(structures=="cnonmonotonicity")
               cnonmonotonicity <- do.call(stops::c_nonmonotonicity,c(list(confs),strucpars[[indst]]))
           }
        if("cfunctionality"%in%structures)
            {
               indst <- which(structures=="cfunctionality")
               cfunctionality <- do.call(stops::c_functionality,c(list(confs),strucpars[[indst]])) 
           }
         if("ccomplexity"%in%structures)
            {
               indst <- which(structures=="ccomplexity")
               ccomplexity <- do.call(stops::c_complexity,c(list(confs),strucpars[[indst]])) 
           }
        if("cfaithfulness"%in%structures)
            {
               indst <- which(structures=="cfaithfulness")
               cfaithfulness <- do.call(stops::c_faithfulness,c(list(confs),strucpars[[indst]]))$mda 
           }
        ##TODO add more structures
        struc <- unlist(mget(structures))
        ic <- stressi*stressweight + sum(struc*strucweight) 
        if (type =="multiplicative") ic <- exp(stressweight*log(stressi) + sum(strucweight*log(struc))) #is this what we want? stress/structure or do we want stress - prod(structure)
        if(verbose>0) cat("stoploss =",ic,"mdsloss =",stressi,"structuredness =",struc,"parameters =",pars,"\n")
        #return the full combi of stress and indices or only the aggregated scalars; for aSTOPS and mSTOPS we want the latter but for a Pareto approach we want the first; get rid of the sums in ic if the first is wanted  
        out <- list(stoploss=ic,strucindices=struc,parameters=pars)
        out
     }
#' STOPS version of smacofSym models
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector for transformations
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param ... additional arguments to be passed to the fitting
#' @param structures which structuredness indices to be included in the loss
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
#' @param strucpars the parameters for the structuredness indices
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
#' 
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress
#'         \item{stress.m:} default normalized stress
#'         \item{stoploss:} the weighted loss value
#'         \item{indices:} the values of the structuredness indices
#'         \item{parameters:} the parameters used for fitting 
#'         \item{fit:} the returned object of the fitting procedure
#'         \item{indobj:} the index objects
#' }
#' 
#'@keywords multivariate
#'@import smacof
#'@import cordillera
#'@export
stop_smacofSym <- function(dis, theta=c(1,1,1), ndim=2,weightmat=NULL,init=NULL,...,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness"),stressweight=1,strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
  if(missing(type)) type <- "additive"
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)<3) theta <- rep(theta, length.out=3)
  lambda <- theta[2]
  fit <- smacof::smacofSym(dis^lambda,ndim=ndim,weightmat=weightmat,init=init,verbose=isTRUE(verbose==2),...) #optimize with smacof
  fit$kappa <- 1
  fit$lambda <- lambda
  fit$nu <- 1
  fit$stress.1 <- fit$stress
  fitdis <- as.matrix(fit$confdiss)
  delts <- as.matrix(fit$delta) #That was my choice to not use the normalized deltas but try it on the original; that is scale and unit free as Buja said
  fit$stress.r <- sum(weightmat*(delts-fitdis)^2)
  fit$stress.m <- fit$stress.r/sum(weightmat*delts^2)
  fit$pars <- c(fit$kappa,fit$lambda,fit$nu)
  fit$deltaorig <- fit$delta^(1/fit$lambda)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
  out <- list(stress=fit$stress, stress.r=fit$stress.r,stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices,parameters=stopobj$parameters,fit=fit,stopobj=stopobj) #target functions
  out
}

#' STOPS versions of flexsmacof models (models with a parametric f() transformation to be determined)
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param transformation function to transform the proximities or distances; need to be parameterized by theta  
#' @param theta the theta vector of transformations
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param structures which structuredness indices to be included in the loss
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param strucpars strucpars
#' @param ... additional arguments to be passed to the fitting
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress
#'         \item{stress.m:} default normalized stress
#'         \item{stoploss:} the weighted loss value
#'         \item{indices:} the values of the structuredness indices
#'         \item{parameters:} the parameters used for fitting 
#'         \item{fit:} the returned object of the fitting procedure
#'         \item{indobj:} the index objects
#' }
#' 
#'@keywords multivariate
#'@import smacof
#'@import cordillera 
#'@export
stop_flexsmacof <- function(dis,transformation=mkPower2, theta=c(1,1), ndim=2,weightmat=NULL,init=NULL,...,structures=c("clusteredness","linearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness"),stressweight=1,strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
  addargs <- list(...)
  #TODO: Other transformations parametrized by theta; use splines
  #Transformations must be so that first argument is the dissimilarity matrix and the second the theta parameters 
  diso <- dis
  dis <- do.call(transformation,list(diso,theta))
  diso <- dis
  dis <- do.call(transformation,list(diso,theta))
  fit <- smacof::smacofSym(dis,ndim=ndim,weightmat=weightmat,init=init,verbose=isTRUE(verbose==2),...) #optimize with smacof
  fit$stress.1 <- fit$stress
  fitdis <- as.matrix(fit$confdiss)
  delts <- as.matrix(fit$delta) #That was my choice to not use the normalized deltas but try it ion the original; that is scale and unit free as Buja said
  fit$stress.r <- sum(weightmat*(delts-fitdis)^2)
  fit$stress.m <- fit$stress.r/sum(weightmat*delts^2)
  fit$pars <- theta
  stopobj <- stoploss(fit,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars)
  out <- list(stress=fit$stress, stress.r=fit$stress.r/2, stress.m=fit$stress.m/2, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters,fit=fit) #target functions
  #TODO include the objects of the indices returned as a list? indicesfull=stopobj 
  out
}

#' STOPS versions of elastic scaling models (via smacofSym)
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this is either a scalar of the lambda transformation for the observed proximities, or a vector where the first is the kappa argument for the fitted distances (here internally fixed to 1) and the second the lambda argument and the third the nu argument (here internally fixed to -2). Defaults to 1 1 -2
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights (NOT the elscal weights)
#' @param init (optional) initial configuration
#' @param ... additional arguments to be passed to the fitting procedure
#' @param structures which structuredness indices to be included in the loss
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
#' @param strucpars the parameters for the structuredness indices
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
#' 
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress
#'         \item{stress.m:} default normalized stress
#'         \item{stoploss:} the weighted loss value
#'         \item{indices:} the values of the structuredness indices
#'         \item{parameters:} the parameters used for fitting 
#'         \item{fit:} the returned object of the fitting procedure
#'         \item{indobj:} the index objects
#' }
#'
#'@importFrom stats dist as.dist
#'@import smacof 
#'@keywords multivariate
#'@export
stop_elastic <- function(dis,theta=c(1,1,-2),ndim=2,weightmat=NULL,init=NULL,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
  #TODO Unfolding
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(missing(type)) type <- "additive"
  #kappa first argument, lambda=second
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)<3) theta <- rep(theta,length.out=3)
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
  lambda <- theta[2]
  nu <- -2
  elscalw <- dis^(nu*lambda) #the weighting in elastic scaling
  diag(elscalw) <- 1
  combwght <- weightmat*elscalw #combine the user weights and the elastic scaling weights
  fit <- smacof::smacofSym(dis^lambda,ndim=ndim,weightmat=combwght,init=init,verbose=isTRUE(verbose==2),...) #optimize with smacof
  fit$kappa <- 1
  fit$lambda <- lambda
  fit$nu <- nu
  fit$stress.1 <- fit$stress
  fitdis <- as.matrix(fit$confdiss)
  delts <- as.matrix(fit$delta) 
  fit$stress.r <- sum(combwght*((delts-fitdis)^2))
  fit$stress.m <- fit$stress.r/sum(combwght*delts^2)
  fit$pars <- c(fit$kappa,fit$lambda,fit$nu)
  fit$deltaorig <- fit$delta^(1/fit$lambda)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
  out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit,stopobj=stopobj) #target functions
  out
}


#' STOPS versions of smacofSphere models
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this is either a scalar of the lambda transformation for the observed proximities, or a vector where the first is the kappa argument for the fitted distances (here internally fixed to 1) and the second the lambda argument and the third the nu argument (here internally fixed to 1). Defaults to 1 1 1
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param ... additional arguments to be passed to the fitting procedure
#' @param structures which structuredness indices to be included in the loss
#' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
#' @param strucpars the parameters for the structuredness indices
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
#' 
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress
#'         \item{stress.m:} default normalized stress
#'         \item{stoploss:} the weighted loss value
#'         \item{indices:} the values of the structuredness indices
#'         \item{parameters:} the parameters used for fitting 
#'         \item{fit:} the returned object of the fitting procedure
#'         \item{indobj:} the index objects
#' }
#'
#'@import smacof 
#'@importFrom stats dist as.dist
#'@keywords multivariate
#'@export
stop_smacofSphere <- function(dis,theta=c(1,1),ndim=2,weightmat=NULL,init=NULL,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
  #TODO Unfolding
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
   if(missing(type)) type <- "additive"
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
  #kappa first argument, lambda=second
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)<3) theta <- rep(theta,length.out=3)
  lambda <- theta[2]
  fit <- smacof::smacofSphere(dis^lambda,ndim=ndim,weightmat=weightmat,init=init,verbose=isTRUE(verbose==2),...) #optimize with smacof
  fit$kappa <- 1
  fit$lambda <- lambda
  fit$nu <- 1
  fit$stress.1 <- fit$stress
  fitdis <- as.matrix(fit$confdiss)
  delts <- as.matrix(fit$delta)[-1,-1]
  fit$stress.r <- sum(weightmat*(delts-fitdis)^2)
  fit$stress.m <- fit$stress.r/sum(weightmat*delts^2)
  fit$pars <- c(fit$kappa,fit$lambda,fit$nu)
  fit$deltaorig <- fit$delta^(1/fit$lambda)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
  out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit,stopobj=stopobj) #target functions
  out
}


#' STOPS version of sammon mapping
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this is either a scalar of the lambda transformation for the observed proximities, or a vector where the first is the kappa argument for the fitted distances (here internally fixed to 1) and the second the lambda argument and the third the nu argument (here internally fixed to -1). Defaults to 1 1 -1
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param ... additional arguments to be passed to the fitting procedure
#' @param structures which structuredness indices to be included in the loss
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
#' @param strucpars the parameters for the structuredness indices
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
#' 
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress
#'         \item{stress.m:} default normalized stress
#'         \item{stoploss:} the weighted loss value
#'         \item{indices:} the values of the structuredness indices
#'         \item{parameters:} the parameters used for fitting 
#'         \item{fit:} the returned object of the fitting procedure
#'         \item{indobj:} the index objects
#' }
#'
#' @importFrom stats dist as.dist
#' @import cordillera
#' @keywords multivariate
#'
#' 
#' @export
stop_sammon <- function(dis,theta=c(1,1,-1),ndim=2,init=NULL,weightmat=NULL,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(missing(type)) type <- "additive"
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(length(theta)==1L) lambda <- theta
  if(length(theta)==2L) lambda <- theta[2]
  if(length(theta)==3L) lambda <- theta[2]
  nu <- -1
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
  fit$pars <- c(fit$kappa,fit$lambda,fit$nu)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
  list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters,  fit=fit,stopobj=stopobj) #target functions
}

#' STOPS versions of Sammon mapping models (via smacofSym)
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this is either a scalar of the lambda transformation for the observed proximities, or a vector where the first is the kappa argument for the fitted distances (here internally fixed to 1) and the second the lambda argument, the thrid the nu argiument (fixed to -1). Defaults to 1 1 -1.
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param ... additional arguments to be passed to the fitting procedure
#' @param structures which structuredness indices to be included in the loss
#' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
#' @param strucpars the parameters for the structuredness indices
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
#' 
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress
#'         \item{stress.m:} default normalized stress
#'         \item{stoploss:} the weighted loss value
#'         \item{indices:} the values of the structuredness indices
#'         \item{parameters:} the parameters used for fitting 
#'         \item{fit:} the returned object of the fitting procedure
#'         \item{indobj:} the index objects
#' }
#'
#'
#' @importFrom stats dist as.dist
#' @import smacof
#' @import cordillera
#' 
#'@keywords multivariate
#'@export
stop_sammon2 <- function(dis,theta=c(1,1,-1),ndim=2,weightmat=NULL,init=NULL,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(missing(type)) type <- "additive"
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1]) 
  #kappa first argument, lambda=second
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)<3) theta <- rep(theta, length.out=3)
  lambda <- theta[2]
  nu <- -1
  elscalw <- dis^(nu*lambda) #the weighting in elastic scaling
  diag(elscalw) <- 1
  combwght <- weightmat*elscalw #combine the user weights and the elastic scaling weights
  fit <- smacof::smacofSym(dis^lambda,ndim=ndim,weightmat=combwght,init=init,verbose=isTRUE(verbose==2),...) #optimize with smacof
  fit$kappa <- 1
  fit$lambda <- lambda
  fit$nu <- nu
  fit$stress.1 <- fit$stress
  fitdis <- as.matrix(fit$confdiss)
  delts <- as.matrix(fit$delta) 
  fit$stress.r <- sum(combwght*((delts-fitdis)^2))
  fit$stress.m <- fit$stress.r/sum(combwght*delts^2)
  fit$pars <- c(fit$kappa,fit$lambda,fit$nu)
  fit$deltaorig <- fit$delta^(1/fit$lambda)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
  out <- list(stress=fit$stress, stress.r=fit$stress.r, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit,stopobj=stopobj) #target functions
  out
}



#' STOPS version of strain
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this is either a scalar of the lambda transformation for the observed proximities, or a vector where the first is the kappa argument for the fitted distances (here internally fixed to 1) and the second and third the lambda and the nu argument (the latter is fixed to 1). Defaults to 1 1 1
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param ... additional arguments to be passed to the fitting procedure
#' @param structures which structuredness indices to be included in the loss
#' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
#' @param strucpars the parameters for the structuredness indices
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
#' 
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} Not really stress but 1-GOF where GOF is the first element returned from cmdscale (the sum of the first ndim absolute eigenvalues divided by the sum of all absolute eigenvalues)
#'         \item{stress.m:} default normalized stress
#'         \item{stoploss:} the weighted loss value
#'         \item{indices:} the values of the structuredness indices
#'         \item{parameters:} the parameters used for fitting 
#'         \item{fit:} the returned object of the fitting procedure
#'         \item{indobj:} the index objects
#' }
#'
#' @import cordillera
#' @importFrom stats dist as.dist
#' @keywords multivariate
#' @export
stop_cmdscale <- function(dis,theta=c(1,1,1),weightmat=NULL,ndim=2,init=NULL,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(missing(type)) type <- "additive"
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
  fit$pars <- c(fit$kappa,fit$lambda,fit$nu)
  fit$conf <- fit$points
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
  list(stress=1-fit$GOF[1],stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj) #target functions
}



#' STOPS version of isomap.
#'
#' Currently this version is a bit less flexible than the vegan one, as the only allowed parameter for isomap is the theta (k in isomap, no epsilon) and the shortest path is always estimated with argument "shortest". Also note that fragmentedOK is always set to TRUE which means that for theta that is too small only the largest conected group will be analyzed. If that's not wanted just set the theta higher.  
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the number of shortest dissimilarities retained for a point (nearest neighbours), the isomap parameter. Defaults to 3.
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param structures which structuredness indices to be included in the loss
#' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
#' @param strucpars the parameters for the structuredness indices
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
#' 
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} Not really stress but 1-GOF where GOF is the first element returned from cmdscale (the sum of the first ndim absolute eigenvalues divided by the sum of all absolute eigenvalues).
#'         \item{stress.m:} default normalized stress (sqrt explicitly normalized stress; really the stress this time)
#'         \item{stoploss:} the weighted loss value
#'         \item{indices:} the values of the structuredness indices
#'         \item{parameters:} the parameters used for fitting 
#'         \item{fit:} the returned object of the fitting procedure
#'         \item{indobj:} the index objects
#' }
#'
#' @import cordillera
#' @importFrom stats dist as.dist
#' @importFrom vegan isomap isomapdist
#' @keywords multivariate
#' @export
stop_isomap <- function(dis,theta=3,weightmat=NULL,ndim=2,init=NULL,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(missing(type)) type <- "additive"
  if(length(theta)==1L) lambda <- theta
  if(length(theta)==2L) lambda <- theta[1]
  if(length(theta)==3L) lambda <- theta[1]
  disi <- vegan::isomapdist(dis,k=lambda,path="shortest",fragmentedOK=TRUE)
  fit <- stops::cmdscale(disi,k=ndim,eig=TRUE) 
  fit$k <- lambda
  #fit$kappa <- 1
  #fit$nu <- 1
  dis <- stats::as.dist(disi)
  fitdis <- stats::dist(fit$points)
  fit$stress.r <- sum((disi-fitdis)^2)
  fit$stress.n <- fit$stress.r/sum(disi^2)
  fit$stress.m <- sqrt(fit$stress.n)
  fit$pars <- fit$k
  fit$conf <- fit$points
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
  list(stress=1-fit$GOF[1],stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj) #target functions
}


#' STOPS version of rstress
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this is either a scalar of the kappa transformation for the fitted distances proximities, or a vector where the first is the kappa argument for the fitted distances and the second the lambda argument, the third the nu argument (here internally fixed to 1). Defaults to 1 1 1
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param ... additional arguments to be passed to the fitting procedure
#' @param structures which structuredness indices to be included in the loss
#' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
#' @param strucpars the parameters for the structuredness indices
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
#' 
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress
#'         \item{stress.m:} default normalized stress
#'         \item{stoploss:} the weighted loss value
#'         \item{indices:} the values of the structuredness indices
#'         \item{parameters:} the parameters used for fitting 
#'         \item{fit:} the returned object of the fitting procedure
#'         \item{indobj:} the index objects
#' }
#'
#' @import cordillera
#' @keywords multivariate
#' @export
stop_rstress <- function(dis,theta=c(1,1,1),weightmat=NULL,init=NULL,ndim=2,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis) 
  if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(missing(type)) type <- "additive"
  if(length(theta)<3) theta <- rep(theta,length.out=3)
  kappa <- theta[1]
  fit <- powerStressMin(delta=dis,kappa=kappa,lambda=1,nu=1,weightmat=weightmat,init=init,ndim=ndim,verbose=verbose,...)
  fit$kappa <- theta[1]
  fit$lambda <- 1
  fit$nu <- 1
  fit$pars <- c(fit$kappa,fit$lambda,fit$nu)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
  out
}


#' STOPS version of sstress
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; this is either a scalar of the lambda transformation for the observed proximities, or a vector where the first is the kappa argument for the fitted distances (here internally fixed to 2) and the second the lambda argument and the third the nu argument (internally fixed to 1). Defaults to 2 1 1
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param ndim the number of dimensions of the target space
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param ... additional arguments to be passed to the fitting procedure
#' @param structures which structuredness indices to be included in the loss
#' @param strucweight weight to be used for the structuredness indices; ; defaults to 1/#number of structures
#' @param strucpars the parameters for the structuredness indices
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param type How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative' 
#' 
#' 
#' @return A list with the components
#'    \itemize{
#'         \item{stress:} the stress
#'         \item{stress.m:} default normalized stress
#'         \item{stoploss:} the weighted loss value
#'         \item{indices:} the values of the structuredness indices
#'         \item{parameters:} the parameters used for fitting 
#'         \item{fit:} the returned object of the fitting procedure
#'         \item{indobj:} the index objects
#' }
#' @import cordillera
#' @keywords multivariate
#' @export
stop_sstress <- function(dis,theta=c(2,1,1),weightmat=NULL,init=NULL,ndim=2,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)  
  if(missing(type)) type <- "additive"
  if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)<3) theta <- rep(theta, length.out=3)
  lambda <- theta[2]
  flambda <- lambda*2 #sstress is d^2 and delta^2 so f(delta^2)=delta^(2*1); lambda works in factors of 2  
  fit <- powerStressMin(delta=dis,kappa=2,lambda=flambda,nu=1,weightmat=weightmat,init=init,ndim=ndim,verbose=verbose,...)
  fit$kappa <- 2
  fit$lambda <- flambda
  fit$nu <- 1
  fit$pars <- c(fit$kappa,fit$lambda,fit$nu)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
  out
}


#' STOPS version of powermds
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances), the second lambda (for the observed proximities), the third nu (for the weights). If a scalar is given it is recycled.  Defaults to 1 1 1.
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param ndim number of dimensions of the target space
#' @param ... additional arguments to be passed to the fitting procedure
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param structures which structures to look for
#' @param strucweight weight to be used for the structures; defaults to 0.5
#' @param strucpars a list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appeacrance in structures 
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param type which weighting to be used in the multi-objective optimization? Either 'additive' (default) or 'multiplicative'. 
#'
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item stoploss: the weighted loss value
#'         \item struc: the structuredness indices
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#' }
#' 
#' @import cordillera
#' @keywords multivariate
#' @export
stop_powermds <- function(dis,theta=c(1,1,1),weightmat=NULL,init=NULL,ndim=2,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis) 
  if(missing(type)) type <- "additive"
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)<3) theta <- rep(theta,length.out=3)
  if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=1,weightmat=weightmat,init=init,ndim=ndim,verbose=verbose,...)
  fit$kappa <- theta[1]
  fit$lambda <- theta[2]
  fit$nu <- 1
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
  out 
}

#' STOPS version of sammon with powers
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances), the second lambda (for the observed proximities), the third nu (fixed to -1). If a scalar is given it is recycled for the free parameters.  Defaults to 1 1 -1.
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param ndim number of dimensions of the target space
#' @param ... additional arguments to be passed to the fitting procedure
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param structures which strcutures to look for
#' @param strucweight weight to be used for the structures; defaults to 0.5
#' @param strucpars a list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appeacrance in structures 
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param type which weighting to be used in the multi-objective optimization? Either 'additive' (default) or 'multiplicative'. 
#'
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item stoploss: the weighted loss value
#'         \item struc: the structuredness indices
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#' }
#'
#' @import cordillera
#' @keywords multivariate
#' @export
stop_powersammon <- function(dis,theta=c(1,1,-1),weightmat=NULL,init=NULL,ndim=2,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)  
  if(missing(type)) type <- "additive"
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)<3) theta <- rep(theta,length.out=3)
  if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  nu <- -1
  sammwght <-dis^(theta[2])
  diag(sammwght) <- 1
  combwght <- sammwght*weightmat
  fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=nu,weightmat=combwght,init=init,ndim=ndim,verbose=verbose,...)
  fit$kappa <- theta[1]
  fit$lambda <- theta[2]
  fit$nu <- nu
  fit$pars <- c(fit$kappa,fit$lambda,fit$nu)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
  out 
}

#' STOPS version of elastic scaling with powers
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances), the second lambda (for the observed proximities) and nu as the third (fixed to -2). If a scalar for the free parameters is given it is recycled.  Defaults to 1 1 -2.
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param ndim number of dimensions of the target space
#' @param ... additional arguments to be passed to the fitting procedure
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param structures which strcutures to look for
#' @param strucweight weight to be used for the structures; defaults to 0.5
#' @param strucpars a list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appeacrance in structures 
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param type which weighting to be used in the multi-objective optimization? Either 'additive' (default) or 'multiplicative'. 
#'
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item stoploss: the weighted loss value
#'         \item struc: the structuredness indices
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#' }
#' 
#' @import cordillera
#' @keywords multivariate
#' @export
stop_powerelastic <- function(dis,theta=c(1,1,-2),weightmat=NULL,init=NULL,ndim=2,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)  
  if(missing(type)) type <- "additive"
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)<3) theta <- rep(theta,length.out=3)
  if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  nu <- -2
  elawght <- dis^(theta[2])
  diag(elawght) <- 1
  combwght <- elawght*weightmat
  fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=nu,weightmat=combwght,init=init,ndim=ndim,verbose=verbose,...)
  fit$kappa <- theta[1]
  fit$lambda <- theta[2]
  fit$nu <- nu
  fit$pars <- c(fit$kappa,fit$lambda,fit$nu)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
  out 
}


#' STOPS version of powerstress
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances), the second lambda (for the observed proximities), the third nu (for the weights). If a scalar is given it is recycled.  Defaults to 1 1 1.
#' @param weightmat (optional) a matrix of nonnegative weights
#' @param init (optional) initial configuration
#' @param ndim number of dimensions of the target space
#' @param ... additional arguments to be passed to the fitting procedure
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param structures which strcutures to look for
#' @param strucweight weight to be used for the structures; defaults to 0.5
#' @param strucpars a list of parameters for the structuredness indices; each list element corresponds to one index in the order of the appeacrance in structures 
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param type which weighting to be used in the multi-objective optimization? Either 'additive' (default) or 'multiplicative'. 
#'
#' @return A list with the components
#' \itemize{
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item stoploss: the weighted loss value
#'         \item struc: the structuredness indices
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#' }
#' @import cordillera
#' @keywords multivariate
#' @export
stop_powerstress <- function(dis,theta=c(1,1,1),weightmat=NULL,init=NULL,ndim=2,...,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative")) {
  theta <- as.numeric(theta)
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(missing(type)) type <- "additive"
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  if(length(theta)<3) theta <- rep(theta,length.out=3)
  if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
  wght <- weightmat
  diag(wght) <- 1
  fit <- powerStressMin(delta=dis,kappa=theta[1],lambda=theta[2],nu=theta[3],weightmat=wght,init=init,ndim=ndim,verbose=verbose,...)
  fit$kappa <- theta[1]
  fit$lambda <- theta[2]
  fit$nu <- theta[3]
  fit$pars <- c(fit$kappa,fit$lambda,fit$nu)
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucindices=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
  out 
}

#' MakePower
#'
#' @param x matrix
#' @param theta numeric (power)
mkPower2<-function(x,theta) {
    r <- theta
    if(length(theta) > 1) r <- theta[2] 
    n<-nrow(x)
    return(abs((x+diag(n))^r)-diag(n))
}

#' High Level STOPS Function 
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param loss which loss function to be used for fitting, defaults to stress
#' @param theta parameters for the proximiy and distance transformation
#' @param structures what c-structuredness should be considered; if missing no structure is considered.
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights; defaults to 1 for all off diagonals 
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param strucweight weight to be used for the cordillera; defaults to -1/length(structures)
#' @param strucpars list of parameters for the structuredness indices; must be in the same ordering as the indices in structures. If missing it is set to NULL.   
#' @param optimmethod What optimizer to use? Currently supported are Bayesian optimization with Gaussian Process priors and Kriging ("Kriging"), Bayesian optimization with treed Gaussian processes ("tgp"), Adaptive LJ Search ("ALJ"), Particle Swarm optimization ("pso"), simulated annealing ("SANN"). Defaults to ALJ version.
#' @param lower The lower contraints of the search region
#' @param upper The upper contraints of the search region 
#' @param verbose numeric value hat prints information on the fitting process; >2 is very verbose
#' @param type which aggregation for the multi objective target function? Either 'additive' (default) or 'multiplicative'
#' @param s number of particles if pso is used
#' @param itmax maximum number of iterations; number of steps of Bayesian optimization if Kriging or tgp is used; default is 50. Note that with tgp the actual number of evaluation of the MDS method is between itmax and 5*itmax as tgp it samples 1-5 candidates from the posterior and uses the best candidate.
#' @param initpoints number of initial points to fit the surrogate model for bayesian optimization; default is 10
#' @param model a character specifying the surrogate model to use. For Kriging it specifies the covariance kernel for the GP prior; see \code{\link{covTensorProduct-class}} defaults to "powerexp". For tgp it specifies the non stationary process used see \code{\link{bgp}}, defaults to "btgpllm" 
#' @param ... additional arguments to be passed to the optimization procedure
#
#'@return A list with the components
#'         \itemize{
#'         \item stoploss: the weighted loss value
#'         TBD
#' }
#' 
#' @examples
#' \donttest{
#' data(BankingCrisesDistances)
#' strucpar<-list(c(eps=10,minpts=2),NULL)
#' res1<-stops(BankingCrisesDistances[,1:69],loss="stress",verbose=0,
#' structures=c("cclusteredness","clinearity"),
#' strucpars=strucpar)
#' res1
#'
#' strucpar<-list(list(alpha=1,C=15,var.thr=1e-5,eps=NULL),list(alpha=1,C=15,var.thr=1e-5,eps=NULL))
#' res1<-stops(BankingCrisesDistances[,1:69],loss="stress",verbose=0,
#' structures=c("cfunctionality","ccomplexity"),
#' strucpars=strucpar)
#' res1
#' }
#' 
#' @importFrom stats dist as.dist optim
#' @importFrom utils capture.output
#' @importFrom pso psoptim
#' @importFrom DiceOptim EGO.nsteps
#' @importFrom DiceKriging km
#' @importFrom tgp lhs dopt.gp
#' @import cordillera
#' 
#' @keywords clustering multivariate
#' @export
stops <- function(dis,loss=c("strain","stress","smacofSym","powerstress","powermds","powerelastic","powerstrain","elastic","sammon","sammon2","smacofSphere","powersammon","rstress","sstress","isomap"), theta=1, structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness"), ndim=2, weightmat=NULL, init=NULL, stressweight=1, strucweight, strucpars, optimmethod=c("SANN","ALJ","pso","Kriging","tgp"), lower=c(1,1,0.5), upper=c(5,5,2), verbose=0, type=c("additive","multiplicative"),s=5,initpoints=10,itmax=50,model,...)
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
      psfunc <- switch(loss, "powerstrain"=stop_cmdscale, "stress"=stop_smacofSym,"smacofSym"=stop_smacofSym,"powerstress"=stop_powerstress,"strain"=stop_cmdscale,"smacofSphere"=stop_smacofSphere,"rstress"=stop_rstress,"sammon"=stop_sammon, "elastic"=stop_elastic, "powermds"=stop_powermds,"powerelastic"=stop_powerelastic,"powersammon"=stop_powersammon,"sammon2"=stop_sammon2,"sstress"=stop_sstress,"isomap"=stop_isomap) #choose the stress to minimize
      if(missing(strucweight)) {
         #TODO: automatic handler of setting weights that makes sense
         strucweight <- rep(-1/length(structures),length(structures))
         if(verbose>1) cat("Weights are stressweight=",stressweight,"strucweights=", strucweight,"\n")
      }
      if(missing(optimmethod)) optimmethod <- "ALJ"
      if(verbose>0) cat("Starting Optimization \n ")
      if(optimmethod=="SANN") {
       opt<- stats::optim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type))$stoploss,method="SANN",...)
        thetaopt <- opt$par
      }
       if(optimmethod=="pso") {
        addargs <- list(...)
        control <- list(trace=verbose-2,s=s,addargs)
        opt<- pso::psoptim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type))$stoploss,lower=lower,upper=upper,control=control,...)
         thetaopt <- opt$par
         bestval <-  opt$value
         itel <- opt$counts["function"]
       }
      if(optimmethod=="ALJ")  {
        opt <- ljoptim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type))$stoploss,lower=lower,upper=upper,verbose=verbose-2,itmax=itmax,...)
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
        responsec <- apply(design, 1, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type))$stoploss) #support points for fitting kriging model
        if (verbose>1) cat("Kriging Model Fitting","\n")
        surrogatemodel <- DiceKriging::km(~1, design = design, response = responsec,covtype=model,control=list(trace=isTRUE(verbose>3))) #fit the kriging model
        #EGO.nsteps has no verbose argument so I capture.output and return it if desired
        if (verbose>2) cat("EGO (DICE) Optimization","\n")
        logged <- capture.output({
           opt<- DiceOptim::EGO.nsteps(model=surrogatemodel, fun=function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type))$stoploss,lower=lower,upper=upper,nsteps=itmax,...)
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
        opt <- tgpoptim(theta, fun=function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type))$stoploss,lower=lower,upper=upper,itmax=itmax,initpoints=initpoints,model=model,verbose=verbose-2,...) #bayesian optimization with treed gaussian process prior
       thetaopt <- opt$par #parameters where best value found (we do not use the last one as that may be worse)
       bestval <- opt$value #best stoploss value
       itel <- opt$counts["function"] 
       }
    #refit optimal model  
    out <- do.call(psfunc,list(dis=dis,theta=thetaopt,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type))
    out$stoploss <- bestval
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

#'@export
print.stops <- function(x,...)
    {
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    cat("Model:",x$type,"STOPS with", x$loss,"loss function and theta parameters=",x$parameters,"\n")
    cat("\n")
    cat("Number of objects:", x$nobj, "\n")
    cat("MDS loss value:", x$stress.m, "\n")
    cat("C-Structuredness Indices:", t(data.frame(names(x$strucindices),x$strucindices)),"\n")
    cat("Structure optimized loss (stoploss):", x$stoploss, "\n")
    cat("MDS loss weight:",x$stressweight,"c-structuredness weights:",x$strucweight,"\n")
    cat("Number of iterations of",x$optimethod,"optimization:", x$optim$counts["function"], "\n")
    cat("\n")
    }

#'@export
#'@importFrom stats coef
coef.stops <- function(object,...)
    {
    return(c(object$par))
    }


#'S3 plot method for stops objects
#' 
#'@param x an object of class stops
#'@param plot.type String indicating which type of plot to be produced: "confplot", "resplot", "Shepard", "stressplot" (see details)
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
#'@importFrom graphics plot 
#'@export 
plot.stops <- function(x,plot.type=c("confplot"), main, asp=NA,...)
    {
     if(missing(plot.type)) plot.type <- "confplot"  
      plot(x$fit,plot.type=plot.type,main=main,asp=asp,...)
 }

#' 3D plots: plot3d method for class stops
#'
#' 
#' This methods produces a dynamic 3D configuration plot.
#' @param x object of class stops
#' @param ... Further plot arguments to the method of the class of slot fit, see \code{\link{plot.smacof}} or \code{\link{plot3d.cmdscale}} . Also see 'rgl' in package 'rgl' 
#'
#'
#' 
#'@export
#'@import rgl
plot3d.stops <- function(x,...)
    {
        plot3d(x$fit,...)
    }

#' 3D plots: plot3dstatic method for class stops
#' 
#' This methods produces a static 3D configuration plot.
#' @param x object of class stops
#' @param ... Further plot arguments to the method of the class of slot fit, see \code{\link{plot3dstatic}} or \code{\link{plot3dstatic.cmdscale}} . Also see 'scatterplot3d' in package 'scatterplot3d'
#'
#'@export
plot3dstatic.stops <- function(x,...)
    {
        plot3dstatic(x$fit,...)
    }

#'@importFrom stats residuals
#'@export
residuals.stops <- function(object,...)
    {
    stats::residuals(object$fit,...)
    }


#'@export
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

#'@export
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

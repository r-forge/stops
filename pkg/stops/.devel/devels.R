#What follows are proof of concepts for the STOPS paper 

#Idea for stops function allow an arbitrary number of indices in a weighted multi-objective optimization way; for this use stoploss_w
# write stops_foo where foo is the MDS model of interest
# do the same as in cops()
# also do this with a pareto approach

#'c-linearity
#'calculates c-linearity as the multiple correlation
#'
#' @param confs a numeric matrix or data frame
#' @param ... additional arguments to be passed to lm.fit
#'
#' @examples
#' x<-1:10
#' y<-2+3*x+rnorm(10)
#' confs<-cbind(x,y)
#' clinearity(confs)
#' 
clinearity <- function(confs,...)
    {
        y <- confs[,1]
        n <- dim(confs)[1]
        p <- dim(confs)[2]
        x <- confs[,2:p]
        tmp <- lm(y~x)
        out <- sqrt(summary(tmp)$r.squared)
        out
    }
    
#' stoploss_w the weighted multiobjective loss function used in STOPS
#'
#' @param obj MDS object (supported are stop_sammon, stop_cmdscale, stop_smacofSym, stop_rstress, stop_powerstress, stop_smacofSphere) 
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param structures which structuredness indices to be included in the loss
#' @param strucweight the weights of the structuredness indices; defaults to 1/#number of structures
#' @param strucpars the parameters to be passed to the structuredness indices
#' @param ... additional parameters passed to the structuredness indices
stoploss_w<- function(obj,stressweight=1,structures=c("c-clusteredness","c-linearity"),strucweight=rep(1/length(structures),length(structures)),strucpars)
    {
       #TODO: figure out how to best process the parameters for the strcuturedness indices; perhaps a list with each element being the parameters per structure index? What used as defaults? Check Caret  
        stressi <- obj$stress.m
        pars <- object$pars
        confs <- obj$conf 
        if("clusteredness"%in%structures)
            {
              indst <- which(structures=="clusteredness")  
              clusteredness <- do.call(cordillera,c(list(confs),strucpar[[indst]]))$normed
            }                           
        if("linearity"%in%structures)
            {
               indst <- which(structures=="linearity")
               linearity <- do.call(clinearity,c(list(confs),strucpar[[indst]]))
           }
        struc <- get(structures)
        ic <- stressweight*stressi - strucweight*struc 
        if(verbose>0) cat(paste("stoploss =",ic,"mdsloss =",stressi,"structuredness =",struc,"parameters =",pars,"\n"))
        #TODO return the full  iobject of the indicies or only the scalars; I think I want the latter
        out <- list(stoploss=ic,strucindices=struc,parameters=pars)
        out
     }
#' STOPS versions of smacofSym models
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
#' @param q the norm of the corrdillera; defaults to 1
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to 2
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the distances (min distance minus max distance)
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param plot plot the cordillera
#' @param normed should the cordillera be normed; defaults to TRUE
#' @param scale should the configuration be scaled to mean=0 and sd=1? Defaults to TRUE
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
#'@export
stop_smacofSym <- function(dis,transformation=mkPower, theta=c(1,1), ndim=2,weightmat=NULL,init=NULL,...,structures=c("clusteredness","linearity"),stressweight=1,strucweight=rep(1/length(structures),length(structures)),strucpars) {
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
  #TODO: Other transformations; I only use mkPower() now
  #Transformations must be so that first argument is the dissimilarity matrix and the second the theta parameters 
  addargs <- list(...)
  addargs
  diso <- dis
  dis <- do.call(transformation,list(diso,theta))
  fit <- smacofSym(dis,ndim=ndim,weightmat=weightmat,init=init,verbose=isTRUE(verbose==2),...) #optimize with smacof
  fit$stress.1 <- fit$stress
  fitdis <- as.matrix(fit$confdiss)
  delts <- as.matrix(fit$delta) #That was my choice to not use the normalized deltas but try it ion the original; that is scale and unit free as Buja said
  fit$stress.r <- sum(weightmat*(delts-fitdis)^2)
  fit$stress.m <- fit$stress.r/sum(weightmat*delts^2)
  fit$pars <- theta
  stopobj <- stoploss_w(fit,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars)
  out <- list(stress=fit$stress, stress.r=fit$stress.r/2, stress.m=fit$stress.m/2, stoploss=stopobj$stoploss, strucindices=stopoobj$strucindices, parameters=stopobj$parameters,fit=fit) #target functions
  #TODO include the objects of the indices returned as a list? indicesfull=stopobj 
  out
}



#' MakePower
#'
#' @param x matrix
#' @param r numeric (power)
mkPower<-function(x,theta) {
    if(length(theta) > 1) r <- theta[2] 
    n<-nrow(x)
    return(abs((x+diag(n))^r)-diag(n))
}

#' High Level STOPS Function 
#'
#' Currently only COPS, the STOPS model for c-clusteredness is implemented
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param loss which loss function to be used for fitting, defaults to strain
#' @param transformation function to transform the proximities or distances; need to be parameterized by theta  
#' @param theta parameters for the proximiy and distance transformation
#' @param structures what structuredness should be considered 
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights; defaults to 1 for all off diagonals 
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param strucweight weight to be used for the cordillera; defaults to 0.5
#' @param strucpars list of parameters for the structuredness indices; must be in the same ordering as the indices in structures  
#' @param optimmethod What general purpose optimizer to use? defaults to our LJ version
#' @param lower The lower contraints of the search region
#' @param upper The upper contraints of the search region 
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param ... additional arguments to be passed to the optimization procedure
#
#' @return see \code{\link{cops}}
#' 
#' @examples
#' data(BankingCrisesDistances)
#' res1<-stops(BankingCrisesDistances[,1:69],structure=c("clusteredness","linearity"),loss="strain",verbose=0)
#' res1
#'
#' @keywords clustering multivariate
#' @export
stops <- function(dis,loss=c("smacofSym"), transformation=mkPower, theta, structures=c("clusteredness","linearity"), ndim=2, weightmat=1-diag(nrow(dis)), init=NULL, stressweight=1, strucweight, strucpars, optimmethod="LJ", lower, upper, verbose=0, ...)
    {
      if(missing(loss)) loss <- "strain"
      if(missing(optimmethod)) optimmethod <- "LJ"
      .confin <- init #initialize a configuration
      psfunc <- switch(loss, "strain"=stop_cmdscale,"stress"=stop_smacofSym,"smacofSym"= stop_smacofSym,"smacofSphere"=stop_smacofSphere,"rstress"=stop_rstress,"powerstress"=stop_powerstress,"sammon"=stop_sammon) #choose the stress to minimize    
      if(missing(strucweight)) { #automatic handler of setting weights #TODO 
         if(verbose>1) cat("Weights are stressweight=",stressweight,"strucweights=", strucweight,"\n")
      }
      if(verbose>1) cat("Starting Optimization \n ")
      if(optimmethod=="LJ")  opt<- ljoptim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars))$stoploss,lower=lower,upper=upper,verbose=verbose,...)
    thetaopt <- opt$par
    out <- do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars))
    out$stoploss <- opt$value
    out$optim <- opt
    out$stressweight <- stressweight
    out$strucweight <- strucweight
    out$call <- match.call()
    out$optimethod <- optimmethod
    out$losstype <- loss
    out$nobj <- dim(out$fit$conf)[1]
    if(verbose>1) cat("Found minimum after",opt$counts["function"]," iterations at",round(opt$par,4),"with stoploss=",round(out$stoploss,4),"and target loss=",round(out$stress.m,4),"and structuredness indices=", round(out$strucindices,4),". Thanks for your patience. \n")
    class(out) <- c("stops")
    out
  }

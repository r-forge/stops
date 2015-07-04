#What follows are proof of concepts for the STOPS paper 

#Idea for stops function allow an arbitrary number of indices in a weighted multi-objective optimization way; for this use stoplose
# write stops_foo where foo is the MDS model of interest
# also do this with a pareto approach
    
#'  Calculate the weighted multiobjective loss function used in STOPS
#'
#' @param obj MDS object (supported are stop_sammon, stop_cmdscale, stop_smacofSym, stop_rstress, stop_powerstress, stop_smacofSphere) 
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param structures which c-structuredness indices to be included in the loss
#' @param strucweight the weights of the structuredness indices; defaults to 1/#number of structures
#' @param strucpars a list of parameters to be passed to the c-structuredness indices in the same order as the values in structures #(alternatively a named list that has the structure name as the element name)
#' @param type what type of weighted optimization should be used? Can be 'additive' or 'multiplicative'
#' @param verbose verbose output
#' 
stoploss<- function(obj,stressweight=1,structures=c("cclusteredness","clinearity"),strucweight=rep(1/length(structures),length(structures)),strucpars,type=c("additive","multiplicative"),verbose=0)
    {
        #TODO make strucpars defaults
        stressi <- obj$stress.m
        pars <- obj$pars
        confs <- obj$conf 
        if("cclusteredness"%in%structures)
            {
              indst <- which(structures=="cclusteredness")  
              cclusteredness <- do.call(cordillera,c(list(confs),strucpars[[indst]]))$normed
            }                           
        if("clinearity"%in%structures)
            {
               indst <- which(structures=="clinearity")
               clinearity <- do.call(c_linearity,c(list(confs))) #,strucpars[[indst]])) has no strucpars
           }
        ##TODO add more structures
        struc <- unlist(mget(structures))
        ic <- stressi*stressweight - sum(struc*strucweight) 
        if (type =="multiplicative") ic <- exp(stressweight*log(stressi) - sum(strucweight*log(struc))) #is this what we want? stress/structure or do we want stress - prod(structure)
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
#' @param stresstype which stress to report. Only takes smacofs default stress currrently.
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
#'@export
stop_smacofSym <- function(dis, theta=c(1,1,1), ndim=2,weightmat=NULL,init=NULL,...,structures=c("cclusteredness","clinearity"),stressweight=1,strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative"), stresstype="default") {
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
  if(missing(type)) type <- "additive"
  if(length(theta)>3) stop("There are too many parameters in the theta argument.")
  lambda <- theta
  if(length(theta)==3L) lambda <- theta[2]
  fit <- smacofSym(dis^lambda,ndim=ndim,weightmat=weightmat,init=init,verbose=isTRUE(verbose==2),...) #optimize with smacof
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
#'@export
stop_flexsmacof <- function(dis,transformation=mkPower2, theta=c(1,1), ndim=2,weightmat=NULL,init=NULL,...,structures=c("clusteredness","linearity"),stressweight=1,strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0) {
  if(inherits(dis,"dist")) dis <- as.matrix(dis)
  if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
  addargs <- list(...)
  #TODO: Other transformations parametrized by theta; use splines
  #Transformations must be so that first argument is the dissimilarity matrix and the second the theta parameters 
  diso <- dis
  dis <- do.call(transformation,list(diso,theta))
  diso <- dis
  dis <- do.call(transformation,list(diso,theta))
  fit <- smacofSym(dis,ndim=ndim,weightmat=weightmat,init=init,verbose=isTRUE(verbose==2),...) #optimize with smacof
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
#' @param stresstype which stress to report? Defaults to explicitly normalized stress
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
#' @keywords multivariate
#' @export
stop_powerstress <- function(dis,theta=c(1,1,1),weightmat=1-diag(nrow(dis)),init=NULL,ndim=2,...,stressweight=1,structures=c("cclusteredness","clinearity"), strucweight=rep(1/length(structures),length(structures)),strucpars,verbose=0,type=c("additive","multiplicative"),stresstype=c("default","stress1","rawstress","normstress","enormstress","enormstress1")) {
  if(missing(stresstype)) stresstype <- "default"
  if(missing(type)) type <- "additive"
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
  stopobj <- stoploss(fit,stressweight=stressweight,structures=structures,strucweight=strucweight,strucpars=strucpars,verbose=isTRUE(verbose>1),type=type)
  out <- list(stress=fit$stress, stress.m=fit$stress.m, stoploss=stopobj$stoploss, strucs=stopobj$strucindices, parameters=stopobj$parameters, fit=fit, stopobj=stopobj)
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
#' @param transformation function to transform the proximities and/or distances; need to be parameterized by theta; currently not used  
#' @param theta parameters for the proximiy and distance transformation
#' @param structures what c-structuredness should be considered 
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights; defaults to 1 for all off diagonals 
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param strucweight weight to be used for the cordillera; defaults to 0.5
#' @param strucpars list of parameters for the structuredness indices; must be in the same ordering as the indices in structures  
#' @param optimmethod What general purpose optimizer to use? defaults to our ALJ version
#' @param lower The lower contraints of the search region
#' @param upper The upper contraints of the search region 
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param type which aggregation for the multi objective target function? Either 'additive' (default) or 'multiplicative'
#' @param s number of particles if pso is used
#' @param stresstype what stress to be used for comparisons between solutions 
#' @param ... additional arguments to be passed to the optimization procedure
#
#' @return see \code{\link{cops}}
#' 
#' @examples
#' data(BankingCrisesDistances)
#' strucpar<-list(c(eps=10,minpts=2),NULL)
#' res1<-stops(BankingCrisesDistances[,1:69],loss="stress",verbose=0,
#' structures=c("cclusteredness","clinearity"),
#' strucpars=strucpar)
#' res1
#'
#' @keywords clustering multivariate
#' @export
stops <- function(dis,loss=c("stress","smacofSym","powerstress"), transformation=mkPower, theta=1, structures=c("cclusteredness","clinearity"), ndim=2, weightmat=1-diag(nrow(dis)), init=NULL, stressweight=1, strucweight, strucpars, optimmethod=c("SANN","ALJ","pso"), lower=c(1,1,0.5), upper=c(5,5,2), verbose=0, type=c("additive","multiplicative"),s=4,stresstype="default",...)
    {
      #TODO add more transformations for the g() and f() by the transformation argument. We only use power versions right now, flexsmacof will allow for more (splines or a smoother or so)
      if(inherits(dis,"dist")) dis <- as.matrix(dis)
      if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
      if(missing(loss)) loss <- "stress"
      if(missing(type)) type <- "additive"
      #TODO implement a Pareto idea
      .confin <- init #initialize a configuration
      psfunc <- switch(loss, "stress"=stop_smacofSym,"smacofSym"=stop_smacofSym,"powerstress"=stop_powerstress)# "strain"=stop_cmdscale,,"smacofSphere"=stop_smacofSphere,"rstress"=stop_rstress,"sammon"=stop_sammon) #choose the stress to minimize    
      if(missing(strucweight)) {
         #TODO: automatic handler of setting weights that makes sense
         strucweight <- rep(1/length(structures),length(structures))
         if(verbose>1) cat("Weights are stressweight=",stressweight,"strucweights=", strucweight,"\n")
      }
      if(missing(optimmethod)) optimmethod <- "ALJ"
      if(verbose>1) cat("Starting Optimization \n ")
      if(optimmethod=="SANN") {
       opt<- optim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,stresstype=stresstype))$stoploss,method="SANN",...)
      }
       if(optimmethod=="pso") {
        addargs <- list(...)
        control <- list(trace=verbose-2,s=s,addargs)
        opt<- pso::psoptim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,stresstype=stresstype))$stoploss,lower=lower,upper=upper,control=control)
       }
      if(optimmethod=="ALJ")  opt <- ljoptim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,stresstype=stresstype))$stoploss,lower=lower,upper=upper,verbose=verbose-2,...
  )   
    thetaopt <- opt$par
    #refit optimal model  
    out <- do.call(psfunc,list(dis=dis,theta=thetaopt,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,stresstype=stresstype))
    out$stoploss <- opt$value
    out$optim <- opt
    out$stressweight <- stressweight
    out$strucweight <- strucweight
    out$call <- match.call()
    out$optimethod <- optimmethod
    out$losstype <- loss
    out$nobj <- dim(out$fit$conf)[1]
    out$type <- type
    if(verbose>1) cat("Found minimum after",opt$counts["function"]," iterations at",round(opt$par,4),"with stoploss=",round(out$stoploss,4),"and default scaling loss=",round(out$stress.m,4),"and structuredness indices=",round(out$strucindices,4),". Thanks for your patience. \n")
    class(out) <- c("stops")
    out
  }

#'@export
print.stops <- function(x,...)
    {
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    cat("Model:",x$type," STOPS with", x$loss,"loss function and theta parameters=",x$par,"\n")
    cat("\n")
    cat("Number of objects:", x$nobj, "\n")
    cat("MDS loss value:", x$stress.m, "\n")
    cat("C-Structuredness Indices:", t(data.frame(names(x$strucindices),x$strucindices)),"\n")
    cat("Structure optimized loss (stoploss): ", x$stoploss, "\n")
    cat("MDS loss weight: ",x$stressweight,", c-structuredness weights: ",x$strucweight,"\n",sep="")
    cat("Number of iterations of",x$optimethod,"optimization:", x$optim$counts["function"], "\n")
    cat("\n")
    }

#'@export
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
#' @param ... Further plot arguments to the method of the class of slot fit, see \code{\link{plot3d.smacof}} or \code{\link{plot3d.cmdscale}} . Also see 'rgl' in package 'rgl' 
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
#' @param ... Further plot arguments to the method of the class of slot fit, see \code{\link{plot3dstatic.smacof}} or \code{\link{plot3dstatic.cmdscale}} . Also see 'scatterplot3d' in package 'scatterplot3d'
#'
#'@export
#'@import smacof
plot3dstatic.stops <- function(x,...)
    {
        plot3dstatic(x$fit,...)
    }

#'@export
residuals.stops <- function(object,...)
    {
    residuals(object$fit,...)
    }


#'@export
summary.stops <- function(object,...)
    {
      sppmat <- NULL
      if(!is.null(object$fit$spp))
      { 
           spp.perc <- object$fit$spp/sum(object$spp) * 100
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

#What follows are proof of concepts for the STOPS paper 

#Idea for stops function allow an arbitrary number of indices in a weighted multi-objective optimization way; for this use stoplose
# write stops_foo where foo is the MDS model of interest
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
#' @export
c_linearity <- function(confs,...)
    {
        y <- confs[,1]
        n <- dim(confs)[1]
        p <- dim(confs)[2]
        x <- confs[,2:p]
        tmp <- lm(y~x)
        out <- sqrt(summary(tmp)$r.squared)
        out
    }
    
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
#' STOPS versions of smacofSym models
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
stop_flexsmacof <- function(dis,transformation=mkPower, theta=c(1,1), ndim=2,weightmat=NULL,init=NULL,...,structures=c("clusteredness","linearity"),stressweight=1,strucweight=rep(1/length(structures),length(structures)),strucpars) {
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
  stopobj <- stoploss_w(fit,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars)
  out <- list(stress=fit$stress, stress.r=fit$stress.r/2, stress.m=fit$stress.m/2, stoploss=stopobj$stoploss, strucindices=stopoobj$strucindices, parameters=stopobj$parameters,fit=fit) #target functions
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
#' @param r numeric (power)
mkPower<-function(x,theta) {
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
#' res1<-stops(BankingCrisesDistances[,1:69],structures=c("cclusteredness","clinearity"),loss="strain",verbose=0)
#' res1
#'
#' @keywords clustering multivariate
#' @export
stops <- function(dis,loss=c("stress","smacofSym","powerstress"), transformation=mkPower, theta, structures=c("cclusteredness","clinearity"), ndim=2, weightmat=1-diag(nrow(dis)), init=NULL, stressweight=1, strucweight, strucpars, optimmethod=c("SANN","ALJ","pso"), lower=c(1,1,0.5), upper=c(5,5,2), verbose=0, type=c("additive","multiplicative"),s=4,stresstype="default",...)
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
    cat("Model: ,",x$type," STOPS with", x$loss,"loss function and theta parameters=",x$par,"\n")
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


###Estimating configuration from coploss with given theta; this is what Reviewer 1 suggested. He suggested to use Nelder-Mead, but I think using something new like NEWUOA is better.

copslossMin <- function (delta, kappa=1, lambda=1, nu=1, theta=c(kappa,lambda,nu),lambdamax=lambda, weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, stressweight=1,cordweight,q=1,minpts=2,epsilon=10,rang=NULL,optimmethod=c("Nelder-Mead","Newuoa","ALJ"),verbose=0,plot=FALSE,scale=TRUE,normed=TRUE,stresstype="default", eps = 1e-12, itmax = 100000,...)
{
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    kappa <- theta[1]
    lambda <- theta[2]
    nu <- theta[3]
    if(verbose>0) cat("Minimizing coploss with kappa=",kappa,"lambda=",lambda,"nu=",nu,"\n")
    if(missing(optimmethod)) optimmethod <- "Newuoa"
    if(missing(rang))
        #perhaps put this into the optimization function?
          {
           if(verbose>1) cat ("Fitting configuration for rang. \n")    
           initsol <- powerStressFast(delta,kappa=kappa,lambda=lambda,nu=nu,weightmat=weightmat,ndim=ndim)
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
                 initsol <- powerStressFast(delta,kappa=kappa,lambda=lambda,nu=nu,weightmat=weightmat,ndim=ndim)
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
    if(is.null(init)) xold <- torgerson (delta, p = p)
    xold <- xold/enorm(xold) 
    copsf <- function(x,delta,p,weightmat,stressweight,cordweight,q,minpts,epsilon,rang,plot,scale,normed=normed,...)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=p)
             delta <- delta/enorm(delta,weightmat)
             x <- x/enorm(x)
             ds <- (2*as.matrix(dist(x)))^kappa
             #ds <- (2*sqrt(sqdist(x)))^kappa
             ds <- ds/enorm(ds)
             #print(ds)
             stressi <- sum(weightmat*(ds-delta)^2)#/2
             #stressi <- sum(weightmat*(ds-delta)^2)/sum(weightmat*(ds^2)) # sqrt stress 1 on the normalized transformed proximities and distances; we use this as the value returned by print
             corrd <- stops::cordillera(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,plot=plot,scale=scale)
 #            corrd <- stops::cordillera(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,plot=plot,scale=scale,...)
             struc <- corrd$raw
             #maxstruc <- corrd$normi #was tut das hier?
             if(normed) {
                        struc <- corrd$normed
                        #maxstruc <- 1
                       }
             ic <- stressweight*stressi - cordweight*struc
             if(verbose>3) cat("coploss =",ic,"mdsloss =",stressi,"OC =",struc,"kappa =",kappa,"lambda =",lambda,"nu=",nu,"\n")
             #ic
             stressweight*stressi - cordweight*struc
             #check if stress 1 is the same as in smacof and whether config is the same for
             #using dist ; check
             #using matrices ; check 
             #using enorm with dist; check
             #using enorm with matrix
           }
     if(optimmethod=="Newuoa") {
         optimized <- minqa::newuoa(xold,function(par) copsf(par,delta=delta,p=p,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,plot=plot,scale=scale,normed=normed),control=list(maxfun=itmax,rhoend=eps,iprint=verbose))
         xnew <- matrix(optimized$par,ncol=2)
         itel <- optimized$feval
     }
     if(optimmethod=="Nelder-Mead") {
         optimized <- optim(xold,function(par) copsf(par,delta=delta,p=p,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,plot=plot,scale=scale,normed=normed),control=list(maxit=itmax,trace=verbose))
         xnew <- optimized$par
         itel <- optimized$counts[[1]]2
     }
    
     attr(xnew,"dimnames")[[2]] <- paste("D",1:p,sep="")
     #doutm <- (2*sqrt(sqdist(xnew)))^kappa  #fitted powered euclidean distance but times two
     doutm <- as.matrix(dist(x)^kappa)
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
     if(verbose>1) cat("*** stress (both normalized):",stressen,"; stress 1 (both normalized - default reported):",stressen1,"; sqrt raw stress (both normalized):",sqrt(stressen),"; raw stress (original data):",stressr,"; stress 1 (original data):",stress1,"; explicitly normed stress (original data):",stressn,"; sqrt explicitly normed stress (original data - used in STOPS):",stresss,"; raw stress (proximities normalized):",stresse,"; stress 1 (proximities normalized):", stresse1,"; from optimization: ",optimized$fval,"\n")   
    out <- list(delta=deltaold, obsdiss=delta, confdiss=dout, conf = xnew, pars=c(kappa,lambda,nu), niter = itel, stress=defaultstress, spp=spp, ndim=p, model="Coploss NEWUOA", call=match.call(), nobj = dim(xnew)[1], type = "Power Stress", gamma=NA, stress.m=sqrt(stressn), stress.r=stressr/2, stress.n=stressn, stress.1=stress1, stress.s=stresss,stress.e=stresse,stress.en=stressen, stress.en1=stressen1,stress.e1=stresse1, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat)
    class(out) <- c("smacofP","smacofB","smacof")
    out
}


library(smacof)
library(stops)
data(kinshipdelta)
delta <- kinshipdelta
kappa <- 1
ndim <- p <- 2
lambda <- 1
weightmat <- 1-diag(nrow(delta))
init <- NULL
eps <- 1e-10
itmax <- 100000
nu <- 1
verbose <- 2
stressweight <- 0.5
cordweight <- 0.5
q <- 1
minpts <- 2
epsilon <- 10
rang <- c(0,2.5)
plot <- FALSE
scale <- TRUE
normed <- TRUE

stressr 
stresse 
stressen 
stressen1 
stress1 
stresse1 
stressn 
stresss 



delta <- delta/enorm(delta)
delta <- as.dist(delta)
xold <- torgerson(delta)
init <- xold/enorm(xold) 


innerf <- function(x,delta,p=p)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=p)
             #delta <- as.matrix(delta)
             ds <- dist(x)
             #print(ds)
             sum((ds-delta)^2)#/2
           }

innerf(init,delta)
delta
as.matrix(delta)
undebug(innerf)
teso <- optim(init,function(par) innerf(par,delta=delta,p=2),control=list(maxit=100000))
plot(teso$par,type="n")
text(teso$par,labels=1:15)

library(minqa)
teso <- newuoa(init,function(par) innerf(par,delta=delta,p=2),control=list(maxfun=100000)) #am Besten
xnewuoa <- matrix(teso$par,ncol=2)
plot(xnewuoa,type="n")
text(xnewuoa,labels=1:15)

library(optimx)
init1 <- as.numeric(init)
teso <- optimx(init1,function(par) innerf(par,delta=delta,p=2),method=c("newuoa"),itnmax=100000)
plot(teso$par,type="n")
text(teso$par,labels=1:15)

par(mfrow=c(1,2))
plot(res1$points)
plot(scale(res1$points))

res1 <- cmdscale(delta)
res1$points
teso$par
plot(res1$points)
plot(teso$par)
res2 <- smacofSym(delta)
innerf(teso$par,delta)


###Estimating configuration from powerstress with given theta; using something new like NEWUOA is better.


data(kinshipdelta)
delta <- kinshipdelta

resmac <- smacofSym(kinshipdelta)
resmaj <- powerStressMin(kinshipdelta)
resuoa <- powerStressFast(kinshipdelta)

dis<-smacof::kinshipdelta
res1<-powerStressMaj(as.matrix(dis),kappa=2,lambda=1.5)
res2<-powerStressMin(as.matrix(dis),kappa=2,lambda=1.5,eps=1e-8)
res3 <- powerStressMaj(as.matrix(dis),init=res2$conf,kappa=2,lambda=1.5,eps=1e-16)

res1
res2
plot(res3)
res3

plot(res1)
plot(res2)

resmac
resmaj
resuoa

plot(resmac)
dev.new()
plot(resmaj)
dev.new()
plot(resuoa)


par(mfrow=c(1,2))
plot(res1$points)
plot(scale(res1$points))

res1 <- cmdscale(delta)
res1$points
teso$par
plot(res1$points)
plot(teso$par)
res2 <- smacofSym(delta)
innerf(teso$par,delta)


data(kinshipdelta)
delta <- kinshipdelta
system.time(res1 <- powerStressMin(delta)) #0.267
system.time(res2 <- powerStressFast(delta)) #0.268
system.time(res3 <- smacofSym(delta)) 
res1
res2
par(mfrow=c(1,2))
plot(res1)
plot(res2)


data(BankingCrisesDistances)
delta <- BankingCrisesDistances[,1:69]
system.time(res1 <- powerStressMin(delta)) #0.349
system.time(res2 <- powerStressFast(delta)) #0.35
system.time(res3 <- smacofSym(delta)) 
res1
res2
par(mfrow=c(1,2))
plot(res1)
plot(res2)



data(Pendigits500)
delta <- dist(Pendigits500[,1:16])
system.time(resP1 <- powerStressMin(delta)) #9586 ~ 2.5 Stunden #0.249
system.time(resP2 <- powerStressFast(delta)) #1100 ~ 20 min  
system.time(resP3 <- smacofSym(delta)) #9586 ~ 2.5 Stunden
resP1
resP2
resP3
par(mfrow=c(1,2))
plot(resP1)
plot(resP2)
plot(resP3)

data(CAClimateIndicatorsCountyMedian)
delta <- dist(CAClimateIndicatorsCountyMedian[,2:52])
system.time(res1 <- powerStressMin(delta))
system.time(res2 <- powerStressFast(delta))
system.time(res3 <- smacofSym(delta))
res1
res2
res3
par(mfrow=c(1,2))
plot(res1)
plot(res2)
plot(res3)

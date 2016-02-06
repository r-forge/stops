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
system.time(res4 <- coplossMin(delta,cordweight=0,stressweight=1,verbose=2))
system.time(res5 <- coplossMin(delta,cordweight=0.5,stressweight=0.5,verbose=2))
system.time(res6 <- coplossMin(delta,verbose=2))
system.time(res7 <- coplossMin(delta,kappa=2,lambda=2,cordweight=0.5,stressweight=0.5,verbose=2))
system.time(res8 <- coplossMin(delta,kappa=1.2,lambda=4,verbose=2))
res1
res2
res3
res4
res5
res6
res7
res8
par(mfrow=c(2,2))
plot(res1)
plot(res2)
plot(res3)
plot(res4)
plot(res5)
plot(res6)
plot(res7)
plot(res8)

data(BankingCrisesDistances)
delta <- BankingCrisesDistances[,1:69]
system.time(res1 <- powerStressMin(delta)) #stress 0.349 time 88.06s
system.time(res2 <- powerStressFast(delta)) #0.35 14.015s 
system.time(res3 <- smacofSym(delta))  #0.344 0.13s
system.time(cops0505 <- coplossMin(delta,cordweight=0.5,stressweight=0.5,verbose=2,eps=1e-7))
system.time(cops10 <- coplossMin(delta,cordweight=0,stressweight=1,verbose=2,eps=1e-7))
system.time(copsdef <- coplossMin(delta,verbose=2,eps=1e-7))
system.time(copsdefc <- coplossMin(delta,stressweight=max(cops10$stress,cops10$OC$normed)/(cops10$stress+cops10$OC$normed),cordweight=min(cops10$stress,cops10$OC$normed)/(cops10$stress+cops10$OC$normed),verbose=2,eps=1e-7))
system.time(cops0901 <- coplossMin(delta,stressweight=0.9,cordweight=0.1,verbose=2,eps=1e-7))
system.time(cops095005 <- coplossMin(delta,stressweight=0.95,cordweight=0.05,verbose=2,eps=1e-7))

system.time(cops0505.148 <- coplossMin(delta,kappa=1.4,lambda=8,cordweight=0.5,stressweight=0.5,verbose=2,eps=1e-7))
system.time(cops10.148 <- coplossMin(delta,cordweight=0,stressweight=1,verbose=2,eps=1e-7))
system.time(copsdef.148 <- coplossMin(delta,verbose=2,eps=1e-7))
system.time(copsdefc.148 <- coplossMin(delta,kappa=1.4,lambda=8,stressweight=max(cops10.148$stress,cops10.148$OC$normed)/(cops10.148$stress+cops10.148$OC$normed),cordweight=min(cops10.148$stress,cops10.148$OC$normed)/(cops10.148$stress+cops10.148$OC$normed),verbose=2,eps=1e-7))
system.time(cops009005.148 <- coplossMin(delta,kappa=1.4,lambda=8,stressweight=0.95,cordweight=0.05,verbose=2,eps=1e-7))

res1
res2
res3

par(mfrow=c(2,2))
plot(res1)
#plot(res2)
plot(res3)
plot(cops0505)
plot(cops10)
plot(copsdef)
plot(copsdefc)
plot(cops0505.148)
plot(cops10.148)
plot(copsdef.148)
plot(copsdefc.148)
plot(cops0901)


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


#################testing Bayesian Optimziation

library(DiceOptim)
d <- 2
n <- 10
set.seed(0)
design <- optimumLHS(n, d)
design <- data.frame(design)
names(design) <- c("lambda", "kappa")
response.branin <- apply(design, 1, branin)
fitted.model1 <- km(design = design, response = response.branin)

###################################################
x.grid <- y.grid <- seq(0, 1, length = n.grid <- 25)
design.grid <- expand.grid(x.grid, y.grid)
EI.grid <- apply(design.grid, 1, EI, fitted.model1)


###################################################
z.grid <- matrix(EI.grid, n.grid, n.grid)
contour(x.grid, y.grid, z.grid, 20)
points(design[,1], design[,2], pch = 17, col = "blue")


###################################################
design1 <- design
EI.grid1 <- EI.grid

set.seed(0)
design <- randomLHS(n, d)
design <- data.frame(design)
names(design) <- c("x1", "x2")
response.branin <- apply(design, 1, branin)
fitted.model <- km(design = design, response = response.branin)
EI.grid <- apply(design.grid, 1, EI, fitted.model)
design2 <- design
EI.grid2 <- EI.grid

set.seed(100)
design <- randomLHS(n, d)
design <- data.frame(design)
names(design)<- c("x1", "x2")
response.branin <- apply(design, 1, branin)
fitted.model <- km(design = design, response = response.branin)
EI.grid <- apply(design.grid, 1, EI, fitted.model)
design3 <- design
EI.grid3 <- EI.grid
design <- design1
EI.grid <- EI.grid1

###################################################
par(mfrow = c(1,3))
z.grid <- matrix(EI.grid, n.grid, n.grid)
contour(x.grid, y.grid, z.grid, 20)
points(design[ , 1], design[ , 2], pch = 17, col = "blue")
z.grid2 <- matrix(EI.grid2, n.grid, n.grid)
contour(x.grid, y.grid, z.grid2, 20)
points(design2[,1], design2[ , 2], pch = 17, col = "blue")
z.grid3 <- matrix(EI.grid3, n.grid, n.grid)
contour(x.grid, y.grid, z.grid3, 20)
points(design3[ , 1], design3[ , 2], pch = 17, col = "blue")



d <- 2
nsteps <- 10
lower <- rep(0, d)
upper <- rep(1, d)     
oEGO <- EGO.nsteps(model = fitted.model1, fun = branin, nsteps = nsteps, 
  lower, upper, control = list(pop.size = 20, BFGSburnin = 2))

par(mfrow = c(1, 2))
response.grid <- apply(design.grid, 1, branin)
z.grid <- matrix(response.grid, n.grid, n.grid)
contour(x.grid, y.grid, z.grid, 40)
points(design[ , 1], design[ , 2], pch = 17, col = "blue")
points(oEGO$par, pch = 19, col = "red")
text(oEGO$par[ , 1], oEGO$par[ , 2], labels = 1:nsteps, pos = 3)
EI.grid <- apply(design.grid, 1, EI, oEGO$lastmodel)
z.grid <- matrix(EI.grid, n.grid, n.grid)
contour(x.grid, y.grid, z.grid, 20)
points(design[ , 1], design[ , 2], pch = 17, col = "blue")
points(oEGO$par, pch = 19, col = "red")


###################################################
hartman6.log <- function(x) {-log(-hartman6(x))}
data(mydata)
X.total <- matrix(unlist(mydata), 50, 6)
nb <- 50
X <- X.total[1:nb, ]
y <- apply(X, 1, hartman6.log) 


###################################################
hartman6.mm <- km(design = data.frame(X), response = y, 
  control = list(pop.size = 50,  max.generations = 20, wait.generations = 5, 
  BFGSburnin = 5), optim.method = "gen")


###################################################
nsteps <- 50
# don't run
# res.nsteps <- EGO.nsteps(model = hartman6.mm, fun = hartman6.log, 
# nsteps = nsteps, lower = rep(0, 6), upper = rep(1, 6), 
# parinit = rep(0.5, 6), control = list(pop.size = 50, 
# max.generations = 20, wait.generations = 5, 
# BFGSburnin = 5), kmcontrol = NULL)
#
# To be compared with the current minimum of Harman6:
hartman6.min <- -3.32




##########
library(DiceOptim)
library(stops)
d <- 2   #dimensions
n <- 10  #starting points
set.seed(0)
design <- randomLHS(n, d)
#design <- optimumLHS(n, d) #optimal design
design <- data.frame(design)  
names(design) <- c("kappa", "lambda")
data(kinshipdelta)


funo <- function(x) powerStressMin(kinshipdelta,kappa=x[1],lambda=x[2],nu=1,verbose=1)$stress

response.branin <- apply(design, 1, function(x) powerStressMin(kinshipdelta,kappa=x[1],lambda=x[2],nu=1,verbose=2)$stress)
fitted.model1 <- km(design = design, response = response.branin)

###################################################
x.grid <- y.grid <- seq(0.0001, 1, length = n.grid <- 25)
design.grid <- expand.grid(x.grid, y.grid)
EI.grid <- apply(design.grid, 1, EI, fitted.model1)


###################################################
z.grid <- matrix(EI.grid, n.grid, n.grid)
contour(x.grid, y.grid, z.grid, 20)
points(design[,1], design[,2], pch = 17, col = "blue")


###################################################
design1 <- design
EI.grid1 <- EI.grid

set.seed(666)
design <- randomLHS(n, d)
design <- data.frame(design)
names(design) <- c("kappa", "lambda")
response.branin <- apply(design, 1, function(x) powerStressMin(kinshipdelta,kappa=x[1],lambda=x[2],nu=1,verbose=2)$stress)
fitted.model <- km(design = design, response = response.branin)
EI.grid <- apply(design.grid, 1, EI, fitted.model)
design2 <- design
EI.grid2 <- EI.grid

set.seed(100)
design <- randomLHS(n, d)
design <- data.frame(design)
names(design)<- c("kappa", "lambda")
response.branin <- apply(design, 1, function(x) powerStressMin(kinshipdelta,kappa=x[1],lambda=x[2],nu=1,verbose=2)$stress)
fitted.model <- km(design = design, response = response.branin)
EI.grid <- apply(design.grid, 1, EI, fitted.model)
design3 <- design
EI.grid3 <- EI.grid
design <- design1
EI.grid <- EI.grid1

###################################################
par(mfrow = c(1,3))
z.grid <- matrix(EI.grid, n.grid, n.grid)
contour(x.grid, y.grid, z.grid, 20)
points(design[ , 1], design[ , 2], pch = 17, col = "blue")
z.grid2 <- matrix(EI.grid2, n.grid, n.grid)
contour(x.grid, y.grid, z.grid2, 20)
points(design2[,1], design2[ , 2], pch = 17, col = "blue")
z.grid3 <- matrix(EI.grid3, n.grid, n.grid)
contour(x.grid, y.grid, z.grid3, 20)
points(design3[ , 1], design3[ , 2], pch = 17, col = "blue")


library(DiceOptim)
d <- 2   #dimensions
n <- 10
set.seed(0)
design <- randomLHS(n, d) #optimal design
design <- data.frame(design)  
names(design) <- c("kappa", "lambda")
library(stops)
data(kinshipdelta)
funo <- function(x) powerStressMin(kinshipdelta,kappa=x[1],lambda=x[2],nu=1,verbose=1)$stress
funo <- function(x) stop_powermds(dis=kinshipdelta,theta=x,structures="cclusteredness",stressweight=0.8,strucweight=-0.2,strucpars=list(minpts=2,epsilon=10),verbose=1)$stoploss

response.post <- apply(design, 1, funo)
response.postc <- apply(design, 1, funo)

response.ban <- apply(design, 1, fbana)
fitted.model1 <- km(y~1,design = design, response = response.postc, lower=rep(.0001,d), upper=rep(1,d))

###################################################
x.grid <- y.grid <- seq(0.00001, 1, length = n.grid <- 10)
design.grid <- expand.grid(x.grid, y.grid)
EI.grid <- apply(design.grid, 1, EI, fitted.model1)



nsteps <- 50
lower <- rep(0.0001, d)
upper <- rep(1, d)

oEGO <- EGO.nsteps(model = fitted.model1, fun = funo, nsteps = nsteps, lower, upper)

#oEGO <- qEGO.nsteps(model = fitted.model1, fun = funo, nsteps = nsteps, lower, upper)

response.grid <- apply(design.grid, 1, funo)
par(mfrow = c(1, 2))
z.grid <- matrix(response.grid, n.grid, n.grid)
contour(x.grid, y.grid, z.grid,40,main=paste("optimum at",round(min(oEGO$value),6)))
points(design[ , 1], design[ , 2], pch = 17, col = "blue")
points(oEGO$par, pch = 19, col = "red")
text(oEGO$par[, 1], oEGO$par[, 2], labels = 1:nsteps, pos = 3)
points(oEGO$par[which.min(oEGO$value),,drop=FALSE], pch = 19,col="green")
#text(oEGO2$par[ , 1], oEGO2$par[ , 2], labels = 1:nsteps, pos = 3)
EI.grid <- apply(design.grid, 1, EI, oEGO$lastmodel)
z.grid <- matrix(EI.grid, n.grid, n.grid)
contour(x.grid, y.grid, z.grid,40)
points(design[ , 1], design[ , 2], pch = 17, col = "blue")
points(oEGO$par, pch = 19, col = "red")
points(oEGO$par[which.min(oEGO$value),,drop=FALSE], pch = 19, col = "green")

funo <- function(x) stop_powermds(dis=kinshipdelta,theta=x,structures="cclusteredness",stressweight=0.8,strucweight=0.2,strucpars=list(minpts=2,epsilon=10),verbose=1)$stoploss

noise.var <- 0.01
fitted.model2 <- km(y~1, design=design, response=response.postc,
          covtype="gauss", noise.var=rep(noise.var,1,n), 
          lower=rep(.0001,d), upper=rep(1,d), control=list(trace=FALSE))

optim.param <- list()
optim.param$quantile <- 0.95
nsteps <- 100

nEGO <- noisy.optimizer(optim.crit="EQI",model = fitted.model2, optim.param=optim.param, n.ite=nsteps, funnoise = funo, lower=lower, upper=upper,noise.var=noise.var,NoiseReEstimate=FALSE, CovReEstimate=FALSE)

par(mfrow = c(1, 2))
#response.grid <- apply(design.grid, 1, funo)
z.grid <- matrix(response.grid, n.grid, n.grid)
contour(x.grid, y.grid, z.grid,40,main=paste("optimum at",round(nEGO$best.y,6)))
points(design[ , 1], design[ , 2], pch = 17, col = "blue")
points(t(nEGO$history.x), pch = 19, col = "red")
points(nEGO$best.x[1],nEGO$best.x[2], pch = 19, col = "green")
text(t(nEGO$history.x)[ , 1], t(nEGO$history.x)[ , 2], labels = 1:nsteps, pos = 3)
EI.grid <- apply(design.grid, 1, EI, nEGO$model)
z.grid <- matrix(EI.grid, n.grid, n.grid)
contour(x.grid, y.grid, z.grid)
points(design[ , 1], design[ , 2], pch = 17, col = "blue")
points(nEGO$history.x, pch = 19, col = "red")
points(nEGO2$best.x, pch = 19, col = "green")

test <- ljoptim(c(0.5,0.5),funo,lower=lower,upper=upper,itmax=100)


res1<-ljoptim(c(0.8,0.2),fbana,lower=0,upper=1,accd=1e-16,acc=1e-16)

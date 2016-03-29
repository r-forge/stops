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


#################Testing Bayesian Optimziation with DiceOptim

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
n <- 20
lambdamax <- 6
kappamax <- 3
set.seed(0)
design <- maximinLHS(n, d)
design <- optimumLHS(n, d)
design <- design*matrix(c(kappamax,lambdamax),byrow=TRUE,ncol=dim(design)[2],nrow=dim(design)[1]) #to the support of lambda and kappa
design <- data.frame(design)  
names(design) <- c("kappa", "lambda")
library(stops)
data(kinshipdelta)
funo <- function(x) stop_powermds(dis=kinshipdelta,theta=x,structures="cclusteredness",stressweight=0.8,strucweight=-0.2,strucpars=list(minpts=2,epsilon=10),verbose=1)$stoploss

response.postc <- apply(design, 1, funo)

response.postc <- apply(design, 1, do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type))$stoploss)


lower <- rep(0.001, d)
upper <- c(kappamax,lambdamax)

#response.ban <- apply(design, 1, fbana)
fitted.model1 <- km(~1,design = design, response = response.postc, lower=lower, upper=upper) #Uses matern(5/2) as covariance function and constant trend #may be okay  
fitted.model2 <- km(~1,design = design, response = response.postc, covtype="gauss", lower=lower, upper=upper) #Uses gauss as covariance function and constant trend #-> squared exponential probably too smooth
fitted.model3 <- km(~1,design = design, response = response.postc, covtype="matern3_2", lower=lower, upper=upper) #Uses matern(3/2) as covariance function and constant trend #may be better than 5_2
fitted.model4 <- km(~1,design = design, response = response.postc,covtype="exp", lower=lower, upper=upper) #Uses exponential as covariance function and constant trend #likely good as its very rough (is Mattern v=1/2) -> Ornstein Uhlenbeck process but perhaps too slow
fitted.model5 <- km(~1,design = design, response = response.postc, covtype="powexp",control=list(trace=FALSE)) #Uses power exponential as covariance function and constant trend #likely the best as it is estimating the power (for power=1 it is exponential for power=2 it is gauss) and can do anything in between smooth and ragged which we need 
#See http://www.gaussianprocess.org/gpml/chapters/RW4.pdf for covtypes

fitted.model5 <- km(~1,design = design, response = response.postc, covtype="powexp",upper=c(2,2),control=)

fitted.model <- fitted.model5

###################################################
x.grid <- seq(0.001, kappamax, length = n.gridx <- 40)
y.grid <- seq(0.001, lambdamax, length = n.gridy <- 50)
design.grid <- expand.grid(x.grid, y.grid)
EI.grid <- apply(design.grid, 1, EI, fitted.model)

nsteps <- 10
set.seed(210485)
oEGO1 <- EGO.nsteps(model = fitted.model1, fun = funo, nsteps = nsteps, lower, upper)
oEGO2 <- EGO.nsteps(model = fitted.model2, fun = funo, nsteps = nsteps, lower, upper)
oEGO3 <- EGO.nsteps(model = fitted.model3, fun = funo, nsteps = nsteps, lower, upper)
oEGO4 <- EGO.nsteps(model = fitted.model4, fun = funo, nsteps = nsteps, lower, upper)
oEGO5 <- EGO.nsteps(model = fitted.model5, fun = funo, nsteps = nsteps, lower, upper,control=list(trace=FALSE),kmcontrol=list(trace=FALSE))

logs <- capture.output({
    oEGO5 <- EGO.nsteps(model = fitted.model1, fun = funo, nsteps = nsteps,lower,upper=c(2,2));
})
logs
oEGO5

#response.grid <- apply(design.grid, 1, function(x) tryCatch(funo(x),error=function(e) NA))
#save(response.grid,file="respgrid.rda")
load("respgrid.rda")

oEGO <- oEGO3
oEGO <- oEGO5

z.grid <- matrix(response.grid, n.gridx, n.gridy)
contour(x.grid, y.grid, z.grid, 40, main=paste("optimum at",round(min(oEGO$value),6)))
points(design[ , 1], design[ , 2], pch = 17, col = "blue")
points(oEGO$par, pch = 19, col = "red")
#text(oEGO$par[, 1], oEGO$par[, 2], labels = 1:nsteps, pos = 3)
points(oEGO$par[which.min(oEGO$value),,drop=FALSE], pch = 19,col="green")
text(oEGO$par[ , 1], oEGO$par[ , 2], labels = 1:nsteps, pos = 3)
points(test$par[1],test$par[2], pch = 17,col="green")
set.seed(210485)
test <- ljoptim(c(1,1),funo,lower=lower,upper=upper,itmax=100,red=0.99,acc=1e-8,accd=1e-6)

library(rgl)
persp3d(x=x.grid, y=y.grid, z=z.grid,smooth=FALSE,col="lightblue")
#plot3d(x=x.grid, y=y.grid, z=z.grid)

par(mfrow = c(1, 2))
EI.grid <- apply(design.grid, 1, EI, oEGO$lastmodel)
z.grid <- matrix(EI.grid, n.gridx, n.gridy)
contour(x.grid, y.grid, z.grid, 40)
points(design[ , 1], design[ , 2], pch = 17, col = "blue")
points(oEGO$par, pch = 19, col = "red")
points(oEGO$par[which.min(oEGO$value),,drop=FALSE], pch = 19, col = "green")


#test for in STOPS.R
library(DiceOptim)
library(stops)
loss <- "powerstrain" #tested, works
loss <- "powermds" #tested, works
loss <- "powerstress" #tested, works
psfunc <- switch(loss, "powerstrain"=stop_cmdscale, "stress"=stop_smacofSym,"smacofSym"=stop_smacofSym,"powerstress"=stop_powerstress,"strain"=stop_cmdscale,"smacofSphere"=stop_smacofSphere,"rstress"=stop_rstress,"sammon"=stop_sammon, "elastic"=stop_elastic, "powermds"=stop_powermds,"powerelastic"=stop_powerelastic,"powersammon"=stop_powersammon,"sammon2"=stop_sammon2,"sstress"=stop_sstress) #choose the stress to minimize
structures <- c("clinearity","cclusteredness")
strucpars <- vector("list",length(structures))
dis <- kinshipdelta
weightmat <- 1-diag(dim(dis)[1])
type <- "additive"
.confin <- NULL
strucweight <- rep(-1/length(structures),length(structures))
optimmethod <- "BO"
theta <- 1
ndim <- 2
stressweight <- 0.5
verbose <- 4
lower <- c(0.5,0.7,-0.5)
upper <- c(2,4,2)
optdim <- 3 #dimensions
 if(loss%in%c("powerstrain","stress","smacofSym","smacofSphere","strain","sammon","elastic","sammon2","sstress","rstress")) optdim <- 1
        if(loss%in%c("powermds","powerelastic","powersammon","smacofSphere","strain","sammon","elastic","sammon2")) optdim <- 2
        designt <- lhs::optimumLHS(kmpoints, optdim) #optimal latin hypercuve 
        minvals <- matrix(lower,byrow=TRUE,nrow=dim(designt)[1],ncol=dim(designt)[2])
        maxvals <- matrix(upper,byrow=TRUE,nrow=dim(designt)[1],ncol=dim(designt)[2])
        design <- designt*(maxvals-minvals)+minvals #to the support of the minimum/maximum values 
        design <- data.frame(design)  
        responsec <- apply(design, 1, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type))$stoploss)
        surrogatemodel <- DiceKriging::km(~1, design = design, response = responsec,covtype=covtype,control=list(trace=isTRUE(verbose>2)))
        #EGO.nsteps has no verbose argument so I capture.output and return it if desired 
       logged <- capture.output({
           opt<- DiceOptim::EGO.nsteps(model=surrogatemodel, fun=function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type))$stoploss,lower=lower,upper=upper,nsteps=nsteps)
       })
       if(verbose>2) print(logged)
       thetaopt <- opt$par[which.min(opt$value),]
       bestval <- min(opt$value)


d <- 2   #dimensions
n <- 20
lambdamax <- 6
kappamax <- 3
set.seed(0)
design <- maximinLHS(n, d)
design <- optimumLHS(n, d)
design <- design*matrix(c(kappamax,lambdamax),byrow=TRUE,ncol=dim(design)[2],nrow=dim(design)[1]) #to the support of lambda and kappa
design <- data.frame(design)  
names(design) <- c("kappa", "lambda")
library(stops)
data(kinshipdelta)
funo <- function(x) stop_powermds(dis=kinshipdelta,theta=x,structures="cclusteredness",stressweight=0.8,strucweight=-0.2,strucpars=list(minpts=2,epsilon=10),verbose=1)$stoploss

response.postc <- apply(design, 1, funo)

response.postc <- apply(design, 1, do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type))$stoploss)


lower <- rep(0.001, d)
upper <- c(kappamax,lambdamax)

#response.ban <- apply(design, 1, fbana)
fitted.model1 <- km(~1,design = design, response = response.postc, lower=lower, upper=upper) #Uses matern(5/2) as covariance function and constant trend #may be okay  
fitted.model2 <- km(~1,design = design, response = response.postc, covtype="gauss", lower=lower, upper=upper) #Uses gauss as covariance function and constant trend #-> squared exponential probably too smooth
fitted.model3 <- km(~1,design = design, response = response.postc, covtype="matern3_2", lower=lower, upper=upper) #Uses matern(3/2) as covariance function and constant trend #may be better than 5_2
fitted.model4 <- km(~1,design = design, response = response.postc,covtype="exp", lower=lower, upper=upper) #Uses exponential as covariance function and constant trend #likely good as its very rough (is Mattern v=1/2) -> Ornstein Uhlenbeck process but perhaps too slow
fitted.model5 <- km(~1,design = design, response = response.postc, covtype="powexp",control=list(trace=FALSE)) #Uses power exponential as covariance function and constant trend #likely the best as it is estimating the power (for power=1 it is exponential for power=2 it is gauss) and can do anything in between smooth and ragged which we need 
#See http://www.gaussianprocess.org/gpml/chapters/RW4.pdf for covtypes

fitted.model5 <- km(~1,design = design, response = response.postc, covtype="powexp",upper=c(2,2),control=)

fitted.model <- fitted.model5

###################################################
x.grid <- seq(0.001, kappamax, length = n.gridx <- 40)
y.grid <- seq(0.001, lambdamax, length = n.gridy <- 50)
design.grid <- expand.grid(x.grid, y.grid)
EI.grid <- apply(design.grid, 1, EI, fitted.model)

nsteps <- 10
set.seed(210485)
oEGO1 <- EGO.nsteps(model = fitted.model1, fun = funo, nsteps = nsteps, lower, upper)
oEGO2 <- EGO.nsteps(model = fitted.model2, fun = funo, nsteps = nsteps, lower, upper)
oEGO3 <- EGO.nsteps(model = fitted.model3, fun = funo, nsteps = nsteps, lower, upper)
oEGO4 <- EGO.nsteps(model = fitted.model4, fun = funo, nsteps = nsteps, lower, upper)
oEGO5 <- EGO.nsteps(model = fitted.model5, fun = funo, nsteps = nsteps, lower, upper,control=list(trace=FALSE),kmcontrol=list(trace=FALSE))

logs <- capture.output({
    oEGO5 <- EGO.nsteps(model = fitted.model1, fun = funo, nsteps = nsteps,lower,upper=c(2,2));
})
logs
oEGO5

#response.grid <- apply(design.grid, 1, function(x) tryCatch(funo(x),error=function(e) NA))
#save(response.grid,file="respgrid.rda")
load("respgrid.rda")

oEGO <- oEGO3
oEGO <- oEGO5

z.grid <- matrix(response.grid, n.gridx, n.gridy)
contour(x.grid, y.grid, z.grid, 40, main=paste("optimum at",round(min(oEGO$value),6)))
points(design[ , 1], design[ , 2], pch = 17, col = "blue")
points(oEGO$par, pch = 19, col = "red")
#text(oEGO$par[, 1], oEGO$par[, 2], labels = 1:nsteps, pos = 3)
points(oEGO$par[which.min(oEGO$value),,drop=FALSE], pch = 19,col="green")
text(oEGO$par[ , 1], oEGO$par[ , 2], labels = 1:nsteps, pos = 3)
points(test$par[1],test$par[2], pch = 17,col="green")
set.seed(210485)
test <- ljoptim(c(1,1),funo,lower=lower,upper=upper,itmax=100,red=0.99,acc=1e-8,accd=1e-6)

library(rgl)
persp3d(x=x.grid, y=y.grid, z=z.grid,smooth=FALSE,col="lightblue")
#plot3d(x=x.grid, y=y.grid, z=z.grid)

par(mfrow = c(1, 2))
EI.grid <- apply(design.grid, 1, EI, oEGO$lastmodel)
z.grid <- matrix(EI.grid, n.gridx, n.gridy)
contour(x.grid, y.grid, z.grid, 40)
points(design[ , 1], design[ , 2], pch = 17, col = "blue")
points(oEGO$par, pch = 19, col = "red")
points(oEGO$par[which.min(oEGO$value),,drop=FALSE], pch = 19, col = "green")

res1 <- stops(kinshipdelta,loss="powermds",theta=1,structures=c("cclusteredness","clinearity"),optimmethod="BO",verbose=3,lower=c(0.1,0.1),upper=c(2,6),covtype="gauss")

res2 <- stops(kinshipdelta,loss="powermds",theta=1,structures=c("cclusteredness","clinearity"),optimmethod="ALJ",verbose=3,lower=c(0.1,0.1),upper=c(2,6))

res2 <- stops(kinshipdelta,loss="powermds",theta=1,structures=c("cclusteredness","clinearity"),optimmethod="ALJ",verbose=3,lower=c(0.1,0.1),upper=c(2,6))

#################Testing Bayesian Optimziation with tgp
# this allows for nonstaitionay covariance matrices which may be better for our problems with varying degrees of smoothness and jumps
library(tgp)
library(stops)
f <- function(x) stop_powermds(dis=kinshipdelta,theta=x,structures="cclusteredness",stressweight=0.8,strucweight=-0.2,strucpars=list(minpts=2,epsilon=10),verbose=1)$stoploss

d <- 2   #dimensions
lambdamax <- 6
kappamax <- 3
lower <- rep(0.001, d)
upper <- c(kappamax,lambdamax)


n <- 10

set.seed(0)
Xcand <- lhs(500,rect)
X1 <- dopt.gp(n,X=NULL,Xcand)$XX
X <- X1
Z <- apply(X, 1, f)

model <- btgp

nudl <- Z
out <- optim.step.tgp(f, X=X, model=model, Z=nudl, rect=rect, prev=NULL)

out <- progress <- NULL



## plot the progress so far
par(mfrow=c(2,2))
plot(out$obj, layout="surf")
plot(out$obj, layout="as", as="improv")
matplot(progress[,1:nrow(rect)], main="optim results",
xlab="rounds", ylab="x[,1:2]", type="l", lwd=2)
plot(log(progress$improv), type="l", main="max log improv",
xlab="rounds", ylab="max log(improv)")


test1 <- tgpoptim(c(1,1),f,lower=c(0.001,0.001),upper=c(3,6),itmax=10,verbose=2,model=btgpllm)



for(i in 1:10) {
## get recommendations for the next point to sample
out <- optim.step.tgp(f, X=X, Z=Z, rect=rect, model = bgp,prev=out)
## add in the inputs, and newly sampled outputs
X <- rbind(X, out$X)
tmp1 <- apply(out$X,1,f)
Z <- c(Z, tmp1)
## keep track of progress and best optimum
progress <- rbind(progress, out$progress)
print(progress[i,])
}




###################################################
x.grid <- seq(0.001, kappamax, length = n.gridx <- 40)
y.grid <- seq(0.001, lambdamax, length = n.gridy <- 50)

load("respgrid.rda")

z.grid <- matrix(response.grid, n.gridx, n.gridy)
contour(x.grid, y.grid, z.grid, 40, main=paste("optimum at",round(min(progress$z),6)))
points(X1[ , 1], X1[ , 2], pch = 17, col = "blue")
points(progress[,1:nrow(rect)], pch = 19, col = "red")
points(progress[which.min(progress$z),,drop=FALSE], pch = 19,col="green")

text(oEGO$par[ , 1], oEGO$par[ , 2], labels = 1:nsteps, pos = 3)
points(test$par[1],test$par[2], pch = 17,col="green")

set.seed(210485)
test <- ljoptim(c(1,1),funo,lower=lower,upper=upper,itmax=100,red=0.99,acc=1e-8,accd=1e-6)

library(rgl)
persp3d(x=x.grid, y=y.grid, z=z.grid,smooth=FALSE,col="lightblue")

## optimize the simple exponential function
f <- function(x) { exp2d.Z(x)$Z }
## create the initial design with D-optimal candidates
rect <- rbind(c(-2,6), c(-2,6))
Xcand <- lhs(500, rect)
X <- dopt.gp(50, X=NULL, Xcand)$XX
Z <- f(X)
## do 10 rounds of adaptive sampling
out <- progress <- NULL
for(i in 1:10) {
## get recommendations for the next point to sample
out <- optim.step.tgp(f, X=X, Z=Z, rect=rect, prev=out)
## add in the inputs, and newly sampled outputs
X <- rbind(X, out$X)
Z <- c(Z, f(out$X))
## keep track of progress and best optimum
progress <- rbind(progress, out$progress)
print(progress[i,])
}
## plot the progress so far
par(mfrow=c(2,2))
plot(out$obj, layout="surf")
plot(out$obj, layout="as", as="improv")
matplot(progress[,1:nrow(rect)], main="optim results",
xlab="rounds", ylab="x[,1:2]", type="l", lwd=2)
plot(log(progress$improv), type="l", main="max log improv",
xlab="rounds", ylab="max log(improv)")



###
library(stops)


lower <- 0
upper <- 4
set.seed(210485)
resalj <- stops(kinshipdelta,theta=1,loss="stress",structures=c("cclusteredness","cmanifoldness","ccomplexity"),stressweight=1,strucweight=c(-0.7,-0.3,0.1),strucpars=list(list(minpts=3,epsilon=10),list(NULL),list(alpha=1,C=15,var.thr=1e-5,eps=NULL)),verbose=4,lower=lower,upper=upper,optimmethod="ALJ",acc=1e-16,accd=1e-8)

set.seed(210485)
reskrig <- stops(dis,theta=1,loss="sammon",structures=c("cclusteredness","cmanifoldness","ccomplexity"),stressweight=1,strucweight=c(-0.7,-0.3,0.1),strucpars=list(list(minpts=3,epsilon=10),list(NULL),list(alpha=1,C=15,var.thr=1e-5,eps=NULL)),verbose=4,lower=lower,upper=upper,optimmethod="Kriging",model="exp",itmax=40)

set.seed(210485)
restgp <- stops(kinshipdelta,theta=1,loss="stress",structures=c("cclusteredness","cmanifoldness","ccomplexity"),stressweight=1,strucweight=c(-0.7,-0.3,0.1),strucpars=list(list(minpts=3,epsilon=10),list(NULL),list(alpha=1,C=15,var.thr=1e-5,eps=NULL)),verbose=5,lower=lower,upper=upper,optimmethod="tgp",model="btgpllm",itmax=20)

theta <- seq(0,4,by=0.0005)
theta <- c(theta,resalj$par[2],reskrig$par[2],restgp$par[2])
theta <- sort(theta)
pst <- vector("list",length(theta))
for(i in 1:length(theta))
{
pst[[i]] <- stop_smacofSym(dis=kinshipdelta,theta=theta[i],structures=c("cclusteredness","cmanifoldness","ccomplexity"),stressweight=1,strucweight=c(-0.7,-0.3,0.1),strucpars=list(list(minpts=3,epsilon=10),list(NULL),list(alpha=1,C=15,var.thr=1e-5,eps=NULL)),verbose=1)
cat(i,"\n")
}



valos <- lapply(pst,function(x) x$stoploss)
plot(theta,valos,type="l")
abline(v=resalj$par[2],col="red")
abline(v=restgp$par[2],col="green")
abline(v=reskrig$par[2],col="blue")

#########################
library(stops)
data(Pendigits500)
pendss <- Pendigits500
cols <- factor(pendss[,17])
library(colorspace)
levels(cols) <- rainbow_hcl(10,c=70,l=50)



#pendss <- Pendigits500

dis <- dist(pendss[,1:16])
q <- 1
rang <- c(0,0.6)
minpts <- 5
eps <- 10
scale <- TRUE
#dis <- as.matrix(dis)

lower <- 1
upper <- 6
strucweight <- c(-500,-0.5,0.1)
strucpars <- list(list(minpts=minpts,epsilon=eps,rang=rang),list(NULL),list(alpha=1,C=15,var.thr=1e-5,eps=NULL))
structures <- c("cclusteredness","cmanifoldness","ccomplexity")

set.seed(210485)
resalj <- stops(dis,theta=1,loss="sammon",structures=c("cclusteredness","cmanifoldness","ccomplexity"),stressweight=1,strucweight=strucweight,strucpars=strucpars,verbose=4,lower=lower,upper=upper,optimmethod="ALJ",acc=1e-16,accd=1e-6)

set.seed(210485)
reskrig <- stops(dis,theta=1,loss="sammon",structures=c("cclusteredness","cmanifoldness","ccomplexity"),stressweight=1,strucweight=strucweight,strucpars=strucpars,verbose=6,lower=lower,upper=upper,optimmethod="Kriging",model="powexp",itmax=40)

set.seed(210485)
restgp <- stops(dis,theta=1,loss="sammon",structures=c("cclusteredness","cmanifoldness","ccomplexity"),stressweight=1,strucweight=strucweight,strucpars=strucpars,verbose=5,lower=lower,upper=upper,optimmethod="tgp",model="btgpllm",itmax=25)


initsam <- sammon(dis)
samopt <- restgp$fit

names(pendss)[17] <- "digit"
colsopt <- cols1 <- colstraj <- factor(pendss[,17])
library(colorspace)
levels(colsopt) <- rainbow_hcl(10,c=50,l=55)
levels(cols1) <- rainbow_hcl(10,c=50,l=96)
levels(colstraj) <- rainbow_hcl(10,c=50,l=99)
tmp <- conf_adjust(samopt$points,initsam$points)
plot(scale(tmp$other.conf)[,1],scale(tmp$other.conf)[,2],
     col=as.character(colsopt),type="n",xlab="D1",ylab="D2",
     ylim=c(-5,3.7),xlim=c(-3.5,2.5))
points(scale(tmp$other.conf)[,1],scale(tmp$other.conf)[,2],
       col=as.character(cols1),pch=as.character(pendss[,17]))
arrows(scale(tmp$other.conf)[,1],scale(tmp$other.conf)[,2],
       scale(tmp$ref.conf)[,1],scale(tmp$ref.conf)[,2],
       col=as.character(colstraj),length=0.1)
points(scale(tmp$ref.conf)[,1],scale(tmp$ref.conf)[,2],
       col=as.character(colsopt),pch=as.character(pendss[,17]))

library(party)
newdatl1 <- data.frame("D1"=scale(initsam$points)[,1],
                       "D2"=scale(initsam$points)[,2],
                       "label"=factor(pendss[,17]))
m1 <- ctree(label~D1+D2,newdatl1,controls=ctree_control(mincriterion=0.95))
#plot(m1)
newdatlo <- data.frame("D1"=scale(samopt$points)[,1],
                       "D2"=scale(samopt$points)[,2],
                       "label"=factor(pendss[,17]))
m1o <- ctree(label~D1+D2,newdatlo)
#plot(m1o)
mall <- ctree(factor(digit)~.,data=pendss)
library(caret)3
cmato <- confusionMatrix(predict(m1o),pendss[,17])
cmat1 <- confusionMatrix(predict(m1),pendss[,17])
cmatall <- confusionMatrix(predict(mall),pendss[,17])
cmato



###### For values

theta <- seq(5.001,6,by=0.001)
#theta <- c(theta,resalj$par[2],reskrig$par[2],restgp$par[2])
#theta <- sort(theta)
pst <- vector("list",length(theta))
valstop <- valstruc <- valstress <- valstressm <- vector("list",length(theta))
for(i in 1:length(theta))
{
pst <- stop_sammon2(dis,theta=theta[i],structures=structures,init=torgerson(dis),stressweight=1,strucweight=strucweight,strucpars=strucpars,verbose=1)
valstop[[i]] <- pst$stoploss
valstruc[[i]] <- pst$strucindices
valstress[[i]] <- pst$stress
valstressm[[i]] <- pst$stress.m
cat(i,"\n")
}

pstF <- pst


save(valstop,valstruc,valstress,valstressm,file="sammongridresult2_5-6.rda")
load("sammongridresultbits1-3.rda")
valstop1 <- valstop
valstruc1 <- valstruc
valstress1 <- valstress
valstressm1 <- valstressm
load("sammongridresultbits3-5.rda")
valstop2 <- valstop
valstruc2 <- valstruc
valstress2 <- valstress
valstressm2 <- valstressm
load("sammongridresultbits5-6.rda")
valstop3 <- valstop
valstruc3 <- valstruc
valstress3 <- valstress
valstressm3 <- valstressm



load("sammongridresult2_1-2.rda")
valstop1 <- valstop
valstruc1 <- valstruc
valstress1 <- valstress
valstressm1 <- valstressm
load("sammongridresult2_2-3.rda")
valstop2 <- valstop
valstruc2 <- valstruc
valstress2 <- valstress
valstressm2 <- valstressm
load("sammongridresult2_3-5.rda")
valstop3 <- valstop
valstruc3 <- valstruc
valstress3 <- valstress
valstressm3 <- valstressm
load("sammongridresult2_5-6.rda")
valstop4 <- valstop
valstruc4 <- valstruc
valstress4 <- valstress
valstressm4 <- valstressm


save(resalj,reskrig,restgp,valstop,valstruc,valstress,file="~/svn/egmds/paper/sammongridresultbits.rda")

valstop <- c(valstop1,valstop2,valstop3,valstop4)
valstruc <- c(valstruc1,valstruc2,valstruc3,valstruc4)
valclus <- lapply(valstruc, function(x) x["cclusteredness"]) 
valmani <- lapply(valstruc, function(x) x["cmanifoldness"])
valcomp <- lapply(valstruc, function(x) x["ccomplexity"])
valstress <- c(valstress1,valstress2,valstress3,valstress4)
valstressm <- c(valstressm1,valstressm2,valstressm3,valstressm4)



theta <- seq(1,6,by=0.001)
valos <- unlist(valstop1)
plot(theta,valos,type="l")
abline(v=resalj$par[2],col="red")
abline(v=restgp$par[2],col="green")
abline(v=reskrig$par[2],col="blue")


model <- get(btgpllm)
library(tgp)


###################################
#Shrinkage cops
shrinkCoploss <- function (delta, kappa=1, lambda=1, nu=1, theta=c(kappa,lambda,nu),weightmat=1-diag(nrow(delta)),  ndim = 2, init=NULL,cordweight=1,q=2,minpts=ndim+1,epsilon=10,rang=NULL,optimmethod=c("Nelder-Mead","Newuoa"),verbose=0,scale=TRUE,accuracy = 1e-7, itmax = 100000,...)
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
    r <- kappa/2
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
    shrinkcops <- function(x,delta,r,ndim,weightmat,cordweight,q,minpts,epsilon,rang,...)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=ndim)
             delta <- delta/enorm(delta,weightmat)
             x <- x/enorm(x)
             dnew <- sqdist (x)
             #rnew <- sum (weightmat * delta * mkPower (dnew, r))
             #nnew <- sum (weightmat * mkPower (dnew,  2*r))
             #anew <- rnew / nnew
             resen <- abs(mkPower(dnew,2*r)-delta)
             shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,...)
             #shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale) 
             #shrinkres <- resen-cordweight*resen*shrinkb/(resen+shrinkb)
             shrinkres <- resen*(1-cordweight*shrinkb/(resen+shrinkb))
             diag(shrinkres) <- 0
             #stressi <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
             #corrd <- stops::cordillera(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,...)
             #struc <- corrd$raw
             #if(normed) {
             #           struc <- corrd$normed
             #          }
              #ic <- stressweight*stressi - cordweight*struc
             #
             #Tis weird why does it become best if I have very low weight?

             ic <- sum(shrinkres^2)/2
             if(verbose>2) cat("coploss =",ic,"mdsloss =",sum(resen),"kappa =",kappa,"lambda =",lambda,"nu=",nu,"\n")
             ic
           }
     if(verbose>1) cat("Starting Minimization with",optimmethod,":\"n")
     if(optimmethod=="Newuoa") {
         optimized <- minqa::newuoa(xold,function(par) shrinkcops(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,
                       cordweight=cordweight, q=q,minpts=minpts,epsilon=epsilon,rang=rang
                                   ),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- optimized$feval
         ovalue <-optimized$fval
     }
     if(optimmethod=="Nelder-Mead") {
         optimized <- optim(xold,function(par) shrinkcops(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,
                       cordweight=cordweight, q=q,minpts=minpts,epsilon=epsilon,rang=rang
                                   ),control=list(maxit=itmax,trace=verbose),...)
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
     stressen <- sum(weightmatm*resmat)/2 #raw stress on the normalized proximities and normalized distances 
     if(verbose>1) cat("*** stress (both normalized - for COPS/STOPS):",stress,"; stress 1 (both normalized - default reported):",sqrt(stress),"; stress manual (for debug only):",stressen,"; from optimization: ",ovalue,"\n")   
    out <- list(delta=deltaold, obsdiss=delta, confdiss=dout, conf = xnew, pars=c(kappa,lambda,nu), niter = itel, stress=sqrt(stress), spp=spp, ndim=ndim, model="Coploss NEWUOA", call=match.call(), nobj = dim(xnew)[1], type = "coploss", gamma=NA, stress.m=stress, stress.en=stressen, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat)
    out$par <- theta
    out$loss <- "coploss"
    out$OC <- stops::cordillera(out$conf,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale)
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



shrinkB2 <- function(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,...)
     {
       shift1 <- function (v) {
             vlen <- length(v)
             sv <- as.vector(numeric(vlen), mode = storage.mode(v)) 
             sv[-1] <- v[1:(vlen - 1)]
             sv
        }
        #if(scale) x <- scale(x)
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
       #something weird; what happens to the first jump? 
        for (i in 1:N) {
            indo <- mats[i,]
            Bmat[indo[1],indo[2]] <- Bmat[indo[2],indo[1]] <- indo[3]/(2^(1/q))
        }
        return(Bmat)
     }
    


 shrinkcops <- function(x,delta,r,ndim,weightmat,cordweight,q=2,minpts,epsilon,rang,...)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=ndim)
             delta <- delta/enorm(delta,weightmat)
             x <- x/enorm(x)
             dnew <- sqdist(x)
             #rnew <- sum (weightmat * delta * mkPower (dnew, r))
             #nnew <- sum (weightmat * mkPower (dnew,  2*r))
             #anew <- rnew / nnew
             resen <- abs(mkPower(dnew,r)-delta)
             #resen <- abs(dnew-delta)
             #shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,...)
             shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang) 
             #shrinkres <- resen-cordweight*resen*shrinkb/(resen+shrinkb)
             shrinkres <- resen*(1-cordweight*(shrinkb/(resen+shrinkb)))
             #shrinkres <- resen
             diag(shrinkres) <- 0
             #stressi <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
             #
             #Tis weird why does it become best if I have very low weight?
             ic <- sum(shrinkres^2)/2
             if(verbose>1) cat("coploss =",ic,"mdsloss =",sum(resen^2)/2,"\n")
             ic
           }


#checkl with how it is in the cops.R code

cordweight=1
q <- 1
optimized <- minqa::newuoa(xold,function(par) shrinkcops(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat, cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose))
xnew2a <- matrix(optimized$par,ncol=ndim)

optimizedn <- optim(xold,function(par) shrinkcops(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat, cordweight=cordweight, q=q,minpts=minpts,epsilon=epsilon,rang=rang), method="Nelder-Mead")
xnewn <- matrix(optimizedn$par,ncol=ndim)

library(crs)
#q <- 0
optimizeds <- crs::snomadr(eval.f=function(par) shrinkcops(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,cordweight=cordweight,q=1,minpts=minpts,epsilon=epsilon,rang=rang),n=length(xold),x0=xnew1,bbin=rep(0,30),ub=rep(5,30),lb=rep(-5,30,bbout=0))
xnews <- matrix(optimizeds$solution,ncol=ndim)

#For the situation with dist() and *r we get an interesting result with nomads:
# In standard MDS the genders are separated which puts some close relationship stati on opposite sides (grandson/father and daughter mother; not so when shrinking - the task was to order similarities by relationhsip status NOT gender 

xnews <- xnews/enorm(xnews)
plot(xnews)

optimized$fval
optimizeds$objective
optimizedn$value



par(mfrow=c(2,2))
plot(xnew1,asp=1)
text(xnew1,label=rownames(xold),pos=3)
plot(xnew2,asp=1)
text(xnew2,label=rownames(xold),pos=3)
plot(xnews,asp=1)
text(xnews,label=rownames(xold),pos=3)


par(mfrow=c(2,2))
#xold <- xold/enorm(xold)
plot(xold,asp=1)
text(xold,label=rownames(xold),pos=3)
xnew <- xnew/enorm(xnew)
plot(xnew,asp=1)
text(xnew,label=rownames(xold),pos=3)
xnewn <- xnewn/enorm(xnewn)
plot(xnewn,asp=1)
text(xnewn,label=rownames(xold),pos=3)
xnews <- xnews/enorm(xnews)
plot(xnews,asp=1)
text(xnews,label=rownames(xold),pos=3)



cordweight <-0 
scale <- TRUE
x0s <- xnew
scale <- FALSE
x0ns <- xnew

cordweight <-1 
scale <- TRUE
x1s <- xnew
scale <- FALSE
x1ns <- xnew


optimizednos <- optimized
optimizednos
xnos <- matrix(optimizednos$par,ncol=ndim)
xs <- matrix(optimized$par,ncol=ndim)
xs10 <- matrix(optimized$par,ncol=ndim)

library(stops)
scale <- TRUE
delta <- kinshipdelta
weightmat <- 1-diag(15)
xsma <- powerStressMin(delta)
xsma <- powerStressFast(delta)
x <- xsma$conf
ndim <- 2
cordweight <- 1
r <- 0.5
q <- 1
minpts <- 2
epsilon <- 10
rang <- c(0,1.6)
plot(xsma)
xold <- xsma$conf
itmax <- 100000
accuracy <- 1e-12
verbose=2

shrinkcops(x=x,delta=delta,r=r,ndim=2,weightmat=weightmat,cordweight=cordweight,rang=rang,q=q,minpts=minpts,epsilon=epsilon)

xnew <- xnew/enorm(xnew)

xnew2 <- shrinkCoploss(delta,minpts=minpts,q=q,epsilon=epsilon,ndim=ndim,cordweight=0,init=xold,rang=rang)
xnew <- shrinkCoploss(delta,minpts=minpts,q=q,epsilon=epsilon,ndim=ndim,cordweight=cordweight,init=xold)
xsma
xnew
xnew2

plot(xsma)
par(mfrow=c(1,2))
plot(xold,asp=1)
text(xold,label=rownames(xold),pos=3)
plot(xnew2a,asp=1)
text(xnew2a,label=rownames(xold),pos=3)

par(mfrow=c(1,2))
plot(xold)
text(xold,label=rownames(xold),pos=3)
plot(xnew)
text(xnew,label=rownames(xold),pos=3)

par(mfrow=c(2,2))
plot(x0ns)
plot(x1ns)
plot(x0s)
plot(x1s)

text(xnos,label=rownames(xold),pos=3)
plot(xs)
text(xs,label=rownames(xold),pos=3)
plot(xs10)
text(xs10,label=rownames(xold),pos=3)


cordweight=0
m0 <- shrinkCoploss(delta,minpts=minpts,init=xold,epsilon=epsilon,q=q,ndim=ndim,cordweight=cordweight,weightmat=weightmat,rang=rang)
xnew0 <- m0$conf

xnew0 <- xold
cordweight=1
m1 <- shrinkCoploss(delta,minpts=minpts,init=xnew0,epsilon=epsilon,q=q,ndim=ndim,cordweight=cordweight,weightmat=weightmat,rang=rang)
xnew1 <- m1$conf

cordweight=0.5
m05 <- shrinkCoploss(delta,init=xnew0,minpts=minpts,epsilon=epsilon,q=q,ndim=ndim,cordweight=cordweight,weightmat=weightmat,rang=rang)
xnew05 <- m05$conf


cordweight=0.1
m01 <- shrinkCoploss(delta,init=xnew0,minpts=minpts,epsilon=epsilon,q=q,ndim=ndim,cordweight=cordweight,weightmat=weightmat,rang=rang)
xnew01 <- m01$conf


par(mfrow=c(2,2))
plot(xsma)
plot(m0)
plot(m05)
plot(m1)
plot(m01)

par(mfrow=c(1,2))
plot(m05)
plot(m1)



plot(m0$OC)
plot(m05$OC)
plot(m1$OC)





plot(m1)   #it was so that tehre was a bug in shrinkB where I used indo[3]/(2^1/q); this worked much better! Why; check it




####GENERAL CASE
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
#' @param shrinkweight weight to be used for the shrinkage; defaults to 1
#' @param q the norm to be minimized for the matrix estimation; defaults to 2 (least squares MDS)
#' @param p the power of the fitted the minksowski distance; defaults to 2 (Euclidean distance) 
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to ndim+1
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the minimum reachabilities to be considered. If missing it is found from the initial configuration by taking 1.5 times the maximal minimum reachability of the initial fit. If NULL it will be normed to each configuration's minimum and maximum distance, so an absolute value of goodness-of-clusteredness. Note that the latter is not necessarily desirable when comparing configurations for their relative clusteredness. See also \code{\link{cordillera}}     
#' @param optimmethod What optimizer to use? Defaults to NEWUOA, Nelder-Mead is also supported.
#' @param verbose numeric value hat prints information on the fitting process; >2 is very verbose
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
#'         \item cordillera: the OPTICS cordillera object
#' }
#' 
#'@examples
#'dis<-as.matrix(smacof::kinshipdelta)
#'
#'#Coploss with shrinkage to 0 
#'res1<-shrinkCoploss(dis,shrinkweight=1) 
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
shrinkCoploss <- function (delta, kappa=1, lambda=1, nu=1, theta=c(kappa,lambda,nu),weightmat=1-diag(nrow(delta)),  ndim = 2, init=NULL,shrinkweight=1,q=2,p=2,minpts=ndim+1,epsilon=10,rang=NULL,optimmethod=c("Nelder-Mead","Newuoa"),verbose=0,accuracy = 1e-7, itmax = 100000,...)
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
           crp <- stops::cordillera(init0,q=q,minpts=minpts,epsilon=epsilon,scale=TRUE)$reachplot
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
    if(is.null(init)) xold <- stops::powerStressMin(delta,kappa=kappa,lambda=lambda,nu=nu,ndim=ndim)$conf
    xold <- xold/enorm(xold) 
    shrinkcops <- function(x,delta,r,ndim,weightmat,shrinkweight,q,p,minpts,epsilon,rang,...)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=ndim)
             delta <- delta/enorm(delta,weightmat)
             x <- x/enorm(x)
             dnew <- as.matrix(dist(x,method="minkowski",p=p))
             #rnew <- sum (weightmat * delta * mkPower (dnew, r))
             #nnew <- sum (weightmat * mkPower (dnew,  2*r))
             #anew <- rnew / nnew
             resen <- abs(mkPower(dnew,2*r)-delta)
             shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,...)
             shrinkres <- resen-resen*shrinkweight*shrinkb/(resen+shrinkb)
             #shrinkres <- resen*(1-shrinkweight*shrinkb/(resen+shrinkb))
             diag(shrinkres) <- 0
             #stressi <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
             #corrd <- stops::cordillera(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,...)
             #struc <- corrd$raw
             #if(normed) {
             #           struc <- corrd$normed
             #          }
              #ic <- stressweight*stressi - shrinkweight*struc
             #
             ic <- sum((weightmat^(1/q))*shrinkres^q)^(1/q)
             if(verbose>2) cat("coploss =",ic,"mdsloss =",sum((weightmat^(1/q))*resen^q)^(1/q),"kappa =",kappa,"lambda =",lambda,"nu=",nu,"\n")
             ic
           }
     if(verbose>1) cat("Starting Minimization with",optimmethod,":\"n")
     if(optimmethod=="Newuoa") {
         optimized <- minqa::newuoa(xold,function(par) shrinkcops(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,
                                    shrinkweight=shrinkweight, q=q, p=p, minpts=minpts,epsilon=epsilon,rang=rang),
                                    control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- optimized$feval
         ovalue <-optimized$fval
     }
     if(optimmethod=="Nelder-Mead") {
         optimized <- optim(xold,function(par) shrinkcops(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,
                       shrinkweight=shrinkweight, q=q,p=p,minpts=minpts,epsilon=epsilon,rang=rang),control=list(maxit=itmax,trace=verbose),...)
         xnew <- optimized$par
         itel <- optimized$counts[[1]]
         ovalue <-optimized$val 
     }
     xnew <- xnew/enorm(xnew)
     dnew <- dist (xnew,method="minkowski",p=p)
     attr(xnew,"dimnames")[[1]] <- rownames(delta)
     attr(xnew,"dimnames")[[2]] <- paste("D",1:ndim,sep="")
     #doutm <- (2*sqrt(sqdist(xnew)))
     doutm <- as.matrix(dnew)^kappa  #fitted powered euclidean distance but times p
     #doutm <- (p*dnew)^kappa  #fitted powered euclidean distance but times p
     deltam <- delta
     deltaorigm <- deltaorig
     deltaoldm <- deltaold
     delta <- stats::as.dist(delta)
     deltaorig <- stats::as.dist(deltaorig)
     deltaold <- stats::as.dist(deltaold)
     doute <- doutm/enorm(doutm)
     doute <- stats::as.dist(doute)
     dout <- stats::as.dist(doutm)
     resmat <- abs(as.matrix(doute - delta))^q
     spp <- colMeans(resmat)
     weightmatm <-weightmat
     weightmat <- stats::as.dist(weightmatm)
     stress <- stressen <- sum(weightmatm^(1/q)*resmat)^(1/q) #raw stress on the normalized proximities and normalized distances
     shrinkmat <- shrinkB(xnew,q=q,minpts=minpts,epsilon=epsilon,rang=rang)
     if(q==2){
     rnew <- sum (weightmat * delta * mkPower (dout, r))
     nnew <- sum (weightmat * mkPower (dout,  2*r))
     anew <- rnew / nnew
     stress <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
     }
     if(verbose>1) cat("*** stress (both normalized - for COPS/STOPS):",stress,"; stress 1 (both normalized - default reported):",sqrt(stress),"; stress manual (for debug only):",stressen,"; from optimization: ",ovalue,"\n")   
    out <- list(delta=deltaold, obsdiss=delta, confdiss=dout, conf = xnew, pars=c(kappa,lambda,nu), niter = itel, stress=sqrt(stress), spp=spp, ndim=ndim, model="COPS-0 (Shrinkage Coploss)", call=match.call(), nobj = dim(xnew)[1], type = "coploss", gamma=NA, stress.m=stress, stress.en=stressen, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat,shrinkmat=shrinkmat)
    out$par <- theta
    out$loss <- "coploss"
    out$OC <- stops::cordillera(out$conf,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=TRUE)
    out$coploss <- ovalue
    out$optim <- optimized
    out$cordweight <- shrinkweight
    out$stressweight <- 1
    out$call <- match.call()
    out$optimethod <- optimmethod
    out$losstype <- out$loss
    out$nobj <- dim(out$conf)[1]
    class(out) <- c("cops","smacofP","smacofB","smacof")
    out
}


#########################################################################################################
# so the question here is how to normlize the dissimilarities/distances in COPS-0; this clearly has a bearing on the overall result as I no longer normalize the shrinkB matrix; also I observed that after procrustes the shrinked results are like the unshrinked.

#old shrinkCoploss
library(stops)

#scale <- TRUE
delta <- kinshipdelta
weightmat <- 1-diag(15)
#weightmat <- 1-diag(58)
xsma <- powerStressFast(delta)
plot(xsma)
x <- xsma$conf
ndim <- 2
r <- 0.5
q <- 1
minpts <- 2
epsilon <- 10
rang <- c(0,1.6)
xold <- xsma$conf
itmax <- 100000
accuracy <- 1e-12
verbose=2
c0 <- cordillera(xold,q=1,minpts=minpts,epsilon=epsilon,rang=rang)
c0

#old shrinkcoploss
m1 <- shrinkCoploss(delta,q=q,minpts=minpts,epsilon=epsilon,weightmat=weightmat,init=xold,ndim=ndim,cordweight=1,rang=rang,scale=scale)
m1
plot(m1,plot.type="reachplot")
c1 <- cordillera(m1$conf,q=1,minpts=minpts,epsilon=epsilon,rang=rang,scale=TRUE)
plot(c1)

c1 <- cordillera(m1$conf,q=1,minpts=minpts,epsilon=epsilon,rang=c(0,0.05613),scale=FALSE)
plot(c1)
c1

shrinkB <- function(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,...)
     {
       shift1 <- function (v) {
             vlen <- length(v)
             sv <- as.vector(numeric(vlen), mode = storage.mode(v)) 
             sv[-1] <- v[1:(vlen - 1)]
             sv
        }
        #if(scale) x <- scale(x)
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
       #something weird; what happens to the first jump? 
        for (i in 2:N) {
            indo <- mats[i,]
            Bmat[indo[1],indo[2]] <- Bmat[indo[2],indo[1]] <- indo[3]/(2^(1/q))
        }
        return(Bmat)
     }

x <- xold
x <- scale(xold)

 shrinkcops <- function(x,delta,r=0.5,ndim,weightmat,cordweight,q=2,minpts,epsilon,rang,scale=TRUE,...)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=ndim)
             #try these variants again with kinship and cali:
             #it all about how to normalize so that the shrinkage will play its part
             #take care of r and so if sqrt(sqdist()) use r*2
             #see what recovers the xold config with cordweight=0
             #see what works with procrustes
             #delta enormed, dnew enormed; not very clustered; very spread out
             #delta enormed, x scaled; looks good! -> use this? Procrustes gives the same result with cordweight=0
             #delta enormed, x scaled + enormed; looks good! -> looks best? 
             #delta enormed, x scaled, dnew enormed; looks ok like #2 but a bit better
             #delta enormed, x enormed, dnew normal; looks ok with clusters for kinship but wrong clusters; closest snew and mdsloss
             #if(scale) x <- scale(x)
             delta <- delta/enorm(delta,weightmat)
             #x <- x/enorm(x)
             #dnew <- sqdist(x)
             dnew <- sqrt(sqdist(x))
             dnew <- dnew/enorm(dnew,weightmat)
             dnew <- dnew^2
             #r <- 2*r
             rnew <- sum (weightmat * delta * mkPower (dnew, r))
             nnew <- sum (weightmat * mkPower (dnew,  2*r))
             anew <- rnew / nnew
             snew <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
             resen <- abs(mkPower(dnew,r)-delta)
             #resen <- abs(dnew-delta)
             #shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,...)
             shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang) 
             #shrinkres <- resen-cordweight*resen*shrinkb/(resen+shrinkb)
             shrinkres <- resen*(1-cordweight*(shrinkb/(resen+shrinkb)))
             #shrinkres <- resen
             diag(shrinkres) <- 0
             ic <- sum(shrinkres^2)
             if(verbose>1) cat("coploss =",ic,"mdslossm =",sum(resen^2),"delta(cop/mds)=",ic-sum(resen^2),"mdslosss =",snew,"delta(mds/sma)=",sum(resen^2)-snew,"\n")
             #delta cops/mds should be positive if cordweight is too high,no?
             ic
           }


cordweight <- 1
optimized <- minqa::newuoa(x,function(par) shrinkcops(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose))
xnew <- matrix(optimized$par,ncol=ndim)

par(mfrow=c(1,2))
plot(x,asp=1,pch=20)
text(x,label=rownames(xold),pos=3)
plot(xnew,asp=1,pch=20)
text(xnew,label=rownames(xold),pos=3)

cx <- cordillera(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang)
cx
plot(cx)

cn <- cordillera(xnew,q=q,minpts=minpts,epsilon=epsilon,rang=rang)
cn
plot(cn)


xa <- conf_adjust(x,xnew)$ref.conf
xnewa <- conf_adjust(x,xnew)$other.conf

par(mfrow=c(1,2))
plot(xa,asp=1,pch=20)
text(xa,label=rownames(x),pos=3)
plot(xnewa,asp=1)
text(xnewa,label=rownames(xold),pos=3)

xaa <- scale(x)
xnewaa <- scale(xnew)

par(mfrow=c(1,2))
plot(xaa,asp=1,pch=20)
text(xaa,label=rownames(x),pos=3)
plot(xnewaa,asp=1)
text(xnewaa,label=rownames(xold),pos=3)


points(xnew,col="red")
points(xnew2,col="green")
points(xnew3,col="blue")
points(xnew20,col="black")



par(mfrow=c(2,2))
plot(xold,asp=1)
text(xold,label=rownames(xold),pos=3)
plot(xnew,asp=1)
text(xnew,label=rownames(xold),pos=3)
plot(xnew2,asp=1)
text(xnew2,label=rownames(xold),pos=3)
plot(xnew,asp=1)
text(xnew,label=rownames(xold),pos=3)


xnew2a <- conf_adjust(xold,xnew2)$other.conf
#xoldaa <- conf_adjust(xold,xnew2)$ref.conf

xoldaa <- scale(xold)
xnewaa <- scale(xnew)
xnew2aa <- scale(xnew2)


par(mfrow=c(1,2))
plot(xold,asp=1)
text(xold,label=rownames(xold),pos=3)
plot(xolda,asp=1)
text(xolda,label=rownames(xold),pos=3)


plot(xnew2a,asp=1)
text(xnew2a,label=rownames(xold),pos=3)

par(mfrow=c(2,2))
plot(xoldaa,asp=1)
text(xoldaa,label=rownames(xold),pos=3)
points(xnewaa,col="red")
plot(xnewaa,asp=1)
text(xnewaa,label=rownames(xold),pos=3)
plot(xnew2aa,asp=1)
text(xnew2aa,label=rownames(xold),pos=3)


xoldaaa <- conf_adjust(xoldaa,xnewaa)$ref.conf
xnewaaa <- conf_adjust(xoldaa,xnewaa)$other.conf
xnew2aaa <- conf_adjust(xoldaa,xnew2aa)$other.conf
par(mfrow=c(2,2))
plot(xoldaaa,asp=1)
text(xoldaaa,label=rownames(xold),pos=3)
points(xnewaa,col="red")
plot(xnewaaa,asp=1)
text(xnewaaa,label=rownames(xold),pos=3)
plot(xnew2aaa,asp=1)
text(xnew2aaa,label=rownames(xold),pos=3)

#xoldaa <- conf_adjust(xold,xnew2)$ref.conf

#tested it with cordweight=1; there are differences here


cordweight <-0 
scale <- TRUE
x0s <- xnew
scale <- FALSE
x0ns <- xnew

cordweight <-1 
scale <- TRUE
x1s <- xnew
scale <- FALSE
x1ns <- xnew


optimizednos <- optimized
optimizednos
xnos <- matrix(optimizednos$par,ncol=ndim)
xs <- matrix(optimized$par,ncol=ndim)
xs10 <- matrix(optimized$par,ncol=ndim)


xnew <- xnew/enorm(xnew)

xnew2 <- shrinkCoploss(delta,minpts=minpts,q=q,epsilon=epsilon,ndim=ndim,cordweight=0,init=xold,rang=rang)
xnew <- shrinkCoploss(delta,minpts=minpts,q=q,epsilon=epsilon,ndim=ndim,cordweight=cordweight,init=xold)
xsma
xnew
xnew2

plot(xsma)
par(mfrow=c(1,2))
plot(xold,asp=1)
text(xold,label=rownames(xold),pos=3)
plot(xnew2a,asp=1)
text(xnew2a,label=rownames(xold),pos=3)

par(mfrow=c(1,2))
plot(xold)
text(xold,label=rownames(xold),pos=3)
plot(xnew)
text(xnew,label=rownames(xold),pos=3)

par(mfrow=c(2,2))
plot(x0ns)
plot(x1ns)
plot(x0s)
plot(x1s)

text(xnos,label=rownames(xold),pos=3)
plot(xs)
text(xs,label=rownames(xold),pos=3)
plot(xs10)
text(xs10,label=rownames(xold),pos=3)


cordweight=0
m0 <- shrinkCoploss(delta,minpts=minpts,init=xold,epsilon=epsilon,q=q,ndim=ndim,cordweight=cordweight,weightmat=weightmat,rang=rang)
xnew0 <- m0$conf

xnew0 <- xold
cordweight=1
m1 <- shrinkCoploss(delta,minpts=minpts,init=xnew0,epsilon=epsilon,q=q,ndim=ndim,cordweight=cordweight,weightmat=weightmat,rang=rang)
xnew1 <- m1$conf

cordweight=0.5
m05 <- shrinkCoploss(delta,init=xnew0,minpts=minpts,epsilon=epsilon,q=q,ndim=ndim,cordweight=cordweight,weightmat=weightmat,rang=rang)
xnew05 <- m05$conf


cordweight=0.1
m01 <- shrinkCoploss(delta,init=xnew0,minpts=minpts,epsilon=epsilon,q=q,ndim=ndim,cordweight=cordweight,weightmat=weightmat,rang=rang)
xnew01 <- m01$conf


par(mfrow=c(2,2))
plot(xsma)
plot(m0)
plot(m05)
plot(m1)
plot(m01)

par(mfrow=c(1,2))
plot(m05)
plot(m1)



plot(m0$OC)
plot(m05$OC)
plot(m1$OC)





plot(m1)   #it was so that tehre was a bug in shrinkB where I used indo[3]/(2^1/q); this worked much better! Why; check it

#For the situation with dist() and *r we get an interesting result with nomads:
# In standard MDS the genders are separated which puts some close relationship stati on opposite sides (grandson/father and daughter mother; not so when shrinking - the task was to order similarities by relationhsip status NOT gender 


optimizedn <- optim(xold,function(par) shrinkcops(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat, cordweight=cordweight, q=q,minpts=minpts,epsilon=epsilon,rang=rang), method="Nelder-Mead")
xnewn <- matrix(optimizedn$par,ncol=ndim)

library(crs)
#q <- 0
optimizeds <- crs::snomadr(eval.f=function(par) shrinkcops(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,cordweight=cordweight,q=1,minpts=minpts,epsilon=epsilon,rang=rang),n=length(xold),x0=xnew1,bbin=rep(0,30),ub=rep(5,30),lb=rep(-5,30,bbout=0))
xnews <- matrix(optimizeds$solution,ncol=ndim)

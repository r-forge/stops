#' Profile COPS Function (aka COPS Variant 2)
#'
#' Metaparameter selection for MDS models baseed on the Profile COPS approach (COPS Variant 2). It uses copstress for hyperparameter selection of explicit transformations (currently power transformations). It is a special case of a STOPS model and predated it; \code{\link[stops]{stops}} has more functionality and can be seen as the successor. pcops uses explicitly normalized stress for copstress (not stress-1).
#'
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param loss which loss function to be used for fitting, defaults to strain. Currently allows for the following models:
#' \itemize{
#' \item Power transformations of observed proximities only (theta must be scalar): Strain loss or classical scaling (\code{strain}, workhorse is cmdscale), Kruskall's stress for symmetric matrices (\code{smacofSym} or \code{stress} and \code{smacofSphere} for scaling onto a sphere; workhorse is smacof), Sammon mapping (\code{sammon} or \code{sammon2}; for the earlier the workhorse is sammon from MASS for the latter it is smacof), elastic scaling (\code{elastic}, the workhorse is smacof), Takane et al's s-Stress \code{sstress} (workhorse is powerstressMin)
#' \item Power transformations of fitted distances only (theta must be scalar): De Leeuw's r-stress \code{rstress} (workhorse is powerstressMin)
#' \item Power transformations of fitted distances and observed proximities (theta must be scalar or of length 2): Power MDS (\code{powermds}), Sammon mapping/elastic scaling with powers (\code{powersammon}, \code{powerelastic})
#' \item Power transfomations of fitted distances, observed proximities and weights (theta must be of length 3 at most): powerstress (POST-MDS, \code{powerstress}), restricted powerstress with equal transformations for distances and proximities (\code{rpowerstress}); workhorse is powerstressMin)
#' \item Approximation to power stress (theta must be of length 2): Approximated power stress (\code{apstress}; workhorse is smacof)
#' }
#' @param theta the theta vector of free parameters; see above for the number of free parameters for each loss function. If this is a vector with more elements than necessary, it is either cut (so for a vector of length 3 and a function with 2 free parameters, the first two elements of the vector are used) or there will be an error.  If a scalar is given as argument and the number of free parameters is larger than 1, the scalar will be recycled and this will also make the optimizers equate all free parameters. Bottom line is: Make sure to supply a theta of the correct length. Defaults to 1 for all free parameters.
#' @param type MDS type which may be one of "ratio", interval", "ordinal". Defaults to "ratio". Note not all loss arguments support all types; if not there will be an error and infor which types are supported. In that case choose another type.   
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights; defaults to 1 for all off diagonals 
#' @param init (optional) initial configuration. If not supplied, the Torgerson scaling result of the dissimilarity matrix dis^theta[2]/enorm(dis^theta[2],weightmat) is used.
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param cordweight weight to be used for the cordillera; if missing gets estimated from the initial configuration so that copstress = 0 for theta=c(1,1) 
#' @param q the norm of the cordillera; defaults to 1
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to ndim+1
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the minimum reachabilities to be considered. If missing it is found from the initial configuration by taking 1.5 times the maximal minimum reachability of the model with theta=c(1,1). If NULL it will be normed to each configuration's minimum and maximum distance, so an absolute value of goodness-of-clusteredness. Note that the latter is not necessarily desirable when comparing configurations for their relative clusteredness. See also \code{\link{cordillera}}.     
#' @param optimmethod What general purpose optimizer to use? Defaults to our adaptive LJ version (ALJ). Also allows particle swarm optimization with s particles ("pso") and simulated annealing ("SANN"), "DIRECT" and "DIRECTL", Hooke-Jeeves ("hjk"), StoGo ("stogo"), and "MADS". We recommend not using SANN and pso with the rstress, sstress and the power stress models. We made good experiences with ALJ, stogo, DIRECT and DIRECTL and also MADS. 
#' @param lower A vector of the lower box contraints of the search region. Its length must match the length of theta.
#' @param upper A vector of the upper box contraints of the search region. Its length must match the length of theta. 
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose. Note that for models with some parameters fixed, the iteration progress of the optimizer shows different values also for the fixed parameters because due to the modular setup we always optimize over a three parameter vector. These values are inconsequential however as internally they will be fixed. 
#' @param normed should the cordillera be normed; defaults to TRUE
#' @param scale should the configuration be scaled and/or centered for calculating the cordillera? "std" standardizes each column of the configurations to mean=0 and sd=1 (typically not a good idea), "sd" scales the configuration by the maximum standard devation of any column (default), "proc" adjusts the fitted configuration to the init configuration (or the Togerson scaling solution if init=NULL). This parameter only has an effect for calculating the cordillera, the fitted and returned configuration is NOT scaled.     
#'@param s number of particles if pso is used

#'@param itmaxo iterations of the outer step (optimization over the hyperparmeters; if solver allows it). Defaults to 200.  
#'@param itmaxi iterations of the inner step (optimization of the MDS). Defaults to 10000 (whichis huge).
#'@param acc termination threshold difference of two successive outer minimization steps.
#'@param ... additional arguments to be passed to the optimization procedure
#'
#'@return A list with the components
#'         \itemize{
#'         \item copstress: the weighted loss value
#'         \item OC: the OPTICS cordillera for the scaled configuration (as defined by scale) 
#'         \item optim: the object returned from the optimization procedure
#'         \item stress: the stress (square root of stress.m)
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
#' 
#'@importFrom stats dist as.dist optim sd
#'@importFrom pso psoptim
#'@importFrom nloptr direct directL stogo 
#'@importFrom crs snomadr
#'@importFrom dfoptim hjk
#'@import cordillera
#' 
#'@keywords clustering multivariate
#'@export
pcops <- function(dis,loss=c("stress","smacofSym","smacofSphere","strain","sammon","rstress","powermds","sstress","elastic","powersammon","powerelastic","powerstress","sammon2","powerstrain","apstress","rpowerstress"),type="ratio",weightmat=NULL,ndim=2,init=NULL,theta=c(1,1,1),stressweight=1,cordweight,q=2,minpts=ndim+1,epsilon=100,rang,optimmethod=c("ALJ","pso","SANN","DIRECT","DIRECTL","stogo","MADS","hjk"),lower=0.5,upper=5,verbose=0,scale=c("proc", "sd", "none", "std"),normed=TRUE,s=4,acc=1e-7,itmaxo=200,itmaxi=10000,...)
{
      if(missing(scale)) scale <- "sd"
      if(inherits(dis,"dist")) dis <- as.matrix(dis)
      if(is.null(weightmat)) weightmat <- 1-diag(nrow(dis))
      if(missing(loss)) loss <- "strain"
      #if(length(theta)==1L) expo <- theta
      #if(length(theta)>2) expo <- theta[2]
      #if(is.null(init)) init <- cops::torgerson(dis^expo/enorm(dis^expo,weightmat),p=ndim) #like in powerstressMin
      .confin <- init #initialize a configuration
      psfunc <- switch(loss,"strain"=cop_cmdscale,"powerstrain"=cop_cmdscale,"elastic"=cop_elastic,"sstress"=cop_sstress,"stress"=cop_smacofSym,"smacofSym"= cop_smacofSym,"smacofSphere"=cop_smacofSphere,"rstress"=cop_rstress,"powermds"=cop_powermds,"powerstress"=cop_powerstress,"sammon"=cop_sammon,"sammon2"=cop_sammon2,"powersammon"=cop_powersammon,"powerelastic"=cop_powerelastic,"apstress"=cop_apstress,"rpowerstress"=cop_rpowerstress) #choose the stress to minimize
      if(missing(optimmethod)) optimmethod <- "ALJ"
      if(missing(rang)) 
          {
            if(verbose>1) cat ("Fitting configuration for rang. \n")
            initsol <- do.call(psfunc,list(dis=dis,theta=c(1,1,1),init=.confin,weightmat=weightmat,ndim=ndim,rang=c(0,1),q=q,minpts=minpts,epsilon=epsilon,verbose=verbose-2,scale=scale,normed=normed,itmaxi=itmaxi,type=type))
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
                initsol <- do.call(psfunc,list(dis=dis,theta=c(1,1,1),init=.confin,weightmat=weightmat,ndim=ndim,rang=rang,q=q,minpts=minpts,epsilon=epsilon,verbose=verbose-2,scale=scale,normed=normed,itmaxi=itmaxi,type=type))
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
      cat("THETA IS",theta,"\n") 
      if(optimmethod=="SANN") {
          opt<- stats::optim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,weightmat=weightmat,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,itmaxi=itmaxi,type=type))$copstress,method="SANN",control=list(maxit=itmaxo,trace=verbose-2,reltol=acc),...)
      }
      if(optimmethod=="pso") {
        addargs <- list(...)
        control <- list(trace=verbose-2,s=s,addargs)
        opt<- pso::psoptim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,weightmat=weightmat,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,itmaxi=itmaxi,type=type))$copstress,lower=lower,upper=upper,control=control)
        thetaopt <- opt$par
       }
      if(optimmethod=="ALJ") {
          opt<- cops::ljoptim(theta, function(theta) do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=theta,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,itmaxi=itmaxi,type=type))$copstress,lower=lower,upper=upper,verbose=verbose-2,itmax=itmaxo,acc=acc,...)
            # opt<- cops::ljoptim(theta, function(theta) do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=theta,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,stresstype=stresstype,itmaxi=itmaxi))$copstress,lower=lower,upper=upper,verbose=verbose-2,itmax=itmaxo,acc=acc)
            thetaopt <- opt$par
      }
      if(optimmethod=="DIRECT") {
          opt<- nloptr::direct(function(theta) do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=theta,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,itmaxi=itmaxi,type=type))$copstress,lower=lower,upper=upper,nl.info=isTRUE(all.equal(verbose-2,0)),control=list(maxeval=itmaxo,xtol_rel=acc),...)
            thetaopt <- opt$par
      }
       if(optimmethod=="stogo") {
           opt<- nloptr::stogo(theta,function(theta) do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=theta,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,itmaxi=itmaxi,type=type))$copstress,lower=lower,upper=upper,nl.info=isTRUE(all.equal(verbose-2,0)),maxeval=itmaxo,xtol_rel=acc,...)
             thetaopt <- opt$par
      }
      if(optimmethod=="DIRECTL") {
          opt<- nloptr::directL(function(theta) do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=theta,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,itmaxi=itmaxi,type=type))$copstress,lower=lower,upper=upper,nl.info=isTRUE(all.equal(verbose-2,0)),control=list(maxeval=itmaxo,xtol_rel=acc),...)
            thetaopt <- opt$par
       }
      if(optimmethod=="MADS") {
      #snomard is super stupid with extra parameters    
      eval.f.pars <- function(x,params)
      {
        psfunc <- params[[1]]
        tmplist <- list(theta=x)
        parlist <- c(tmplist,params[-1])
        tmpo <- do.call(psfunc,parlist)
        return(tmpo$copstress)
       }
       params <- list(psfunc,dis=dis,weightmat=weightmat,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-4,scale=scale,normed=normed,itmaxi=itmaxi,type=type)   
       opt<- crs::snomadr(n=length(theta),x0=theta,eval.f=eval.f.pars,params=params,bbin=0,lb=lower,ub=upper,print.output=isTRUE(all.equal(verbose-2,0)),opts=list("MAX_BB_EVAL"=itmaxo),...)
       thetaopt <- opt$solution
       opt$par <- opt$solution  
       }
      if(optimmethod=="hjk") {
          opt<- dfoptim::hjkb(theta, function(theta) do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=theta,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-3,scale=scale,normed=normed,itmaxi=itmaxi,type=type))$copstress,lower=lower,upper=upper,control=list(info=isTRUE(all.equal(verbose-2,0)),maxfeval=itmaxo,tol=acc),...)
       thetaopt <- opt$par
       } 
    #refit the optimal version (TODO probably unnecessary if the other functions are properly reimplemented)
    out <- do.call(psfunc,list(dis=dis,weightmat=weightmat,theta=thetaopt,init=.confin,ndim=ndim,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,verbose=verbose-2,scale=scale,normed=normed,itmaxi=itmaxi,type=type))
    confopt <- scale_adjust(out$fit$conf,.confin,scale=scale)
    out$OC <- cordillera::cordillera(confopt,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=FALSE)
    #out$copstress <- opt$value 
    out$optim <- opt
    out$stressweight <- stressweight
    out$cordweight <- cordweight
    out$call <- match.call()
    out$optimethod <- optimmethod
    out$losstype <- loss
    out$type <- type
    out$nobj <- dim(out$fit$conf)[1]
    out$scale <- scale
    if(verbose>1) cat("Found minimum after",opt$counts["function"]," iterations at",round(thetaopt,4),"with copstress=",round(out$copstress,4),"and default scaling loss=",round(out$stress.m,4),"and OC=", round(out$OC$normed,4),". Thanks for your patience. \n")
    class(out) <- c("pcops","stops","cops")
    out
}



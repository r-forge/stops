## #' Fitting a COPS-C Model (COPS Variant 1).
## #'
## #' Minimizing Copstress to obtain a clustered ratio, interval or ordinal PS configuration with given explicit power transformations theta. The function allows mix-and-match of explicit (via theta) and implicit (via type) transformations by setting the kappa, lambda, nu (or theta) and type arguments.
## #'
## #' This is an extremely flexible approach to least squares proximity scaling: It supports ratio power stress; ratio, interval and ordinal r stress and ratio, interval and ordinal MDS with or without a COPS penalty. Famous special cases of these models that can be fitted are multiscale MDS if kappa->0 and delta=log(delta), Alscal MDS (sstress) with lambda=kappa=2, sammon type mapping with weightmat=delta and nu=-1, elastic scaling with weightmat=delta and nu=-2. Due to mix-and-match this function also allows to fit models that have not yet been published, such as for example an "elastic scaling ordinal s-stress with cops penalty".
## #'
## #' If one wants to fit these models without the cops penalty, we recommend to use powerStressMin (for ratio MDS with any power transformation for weights, dissimilarities and distances) or rStressMin (for interval and ordinal MDS with power transformations for distances and weights) as these use majorization.  
## #'
## #' @rdname copstressMin
## #' 
## #' @param delta numeric matrix or dist object of a matrix of proximities
## #' @param kappa power transformation for fitted distances
## #' @param lambda power transformation for proximities (only used if type="ratio")
## #' @param nu power transformation for weights
## #' @param theta the theta vector of powers; the first is kappa (for the fitted distances if it exists), the second lambda (for the observed proximities if it exist and type="ratio"), the third is nu (for the weights if it exists). If less than three elements are is given as argument, it will be recycled. Defaults to 1 1 1. Will override any kappa, lambda, nu parameters if they are given and do not match.
## #' @param type what type of MDS to fit. Currently one of "ratio", "interval" or "ordinal". Default is "ratio".
## #' @param ties the handling of ties for ordinal (nonmetric) MDS. Possible are "primary" (default), "secondary" or "tertiary".
## #' @param weightmat (optional) a matrix of nonnegative weights; defaults to 1 for all off diagonals
## #' @param ndim number of dimensions of the target space
## #' @param init (optional) initial configuration
## #' @param stressweight weight to be used for the fit measure; defaults to 0.975
## #' @param cordweight weight to be used for the cordillera; defaults to 0.025
## #' @param q the norm of the cordillera; defaults to 1
## #' @param minpts the minimum points to make up a cluster in OPTICS, see [dbscan::optics()] where it is called \code{minPts}; defaults to ndim+1.
## #' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked, see [dbscan::optics()]; defaults to 10. Note this means we do not expect any noise objects per default. This number will rarely be exceeded if we standardize the configuration as is the default in cops. However if no standardization is applied or there is a procrustes adjustment to a configuration with variance of 10 or more on any of the axes, it can have the effect of being too small. In that case just set a much higher epsilon.
## #' @param dmax The winsorization limit of reachability distances in the OPTICS Cordillera. If supplied, it should be either a numeric value that matches max(rang) or NULL; if NULL it is found as 1.5 times (for kappa >1) or 1 times (for kappa <=1) the maximum reachbility value of the power torgerson model with the same lambda. If dmax and rang are supplied and dmax is not max(rang), a warning is given and rang takes precedence.   
## #' @param rang range of the reachabilities to be considered. If missing it is found from the initial configuration by taking 0 as the lower boundary and dmax (see above) as upper boundary. See also \code{\link{cordillera}}     
## #' @param optimmethod What optimizer to use? Choose one string of 'Newuoa' (from package minqa), 'NelderMead', 'hjk' (Hooke-Jeeves algorithm from dfoptim), 'solnl' (from nlcOptim), 'solnp' (from Rsolnp), 'subplex' (from subplex), 'SANN' (simulated annealing), 'BFGS', 'snomadr' (from crs), 'genoud' (from rgenoud), 'gensa' (from GenSA), 'cmaes' (from cmaes) and 'direct' (from nloptr). See the according R packages for details on these solvers. There are also combinations that proved to work well good, like 'hjk-Newuoa', 'hjk-BFGS', 'BFGS-hjk', 'Newuoa-hjk', 'direct-Newuoa' and 'direct-BFGS'. Usually everything with hjk, BFGS, Newuoa, subplex and solnl in it work rather well in an acceptable time frame (depending on the smoothness of copstress). Default is 'hjk-Newuoa'.   
## #' @param verbose numeric value hat prints information on the fitting process; >2 is very verbose
## #' @param normed should the cordillera be normed; defaults to TRUE.
## #' @param scale Allows to scale the configuration for the OC (the scaled configuration is also returned as $conf). One of "none" (so no scaling), "sd" (configuration divided by the highest standard deviation of the columns), "std" (standardize all columns !NOTE: This does not preserve the relative distances of the optimal config), "proc" (procrustes adjustment to the initial fit) and "rmsq" (configuration divided by the maximum root mean square of the columns). Default is "sd".   
## #' @param accuracy numerical accuracy, defaults to 1e-7.
## #' @param itmax maximum number of iterations. Defaults to 10000. For the two-step algorithms if itmax is exceeded by the first solver, the second algorithm is run for at least 0.1*itmax (so overall itmax may be exceeded by a factor of 1.1).
## #' @param stresstype which stress to use in the copstress. Defaults to stress-1. If anything else is set, explicitly normed stress which is (stress-1)^2. Using stress-1 puts more weight on MDS fit.   
## #' @param ... additional arguments to be passed to the optimization procedure
## #'
## #' @return A list with the components
## #'         \itemize{
## #'         \item delta: the original transformed dissimilarities
## #'         \item obsdiss: the explicitly normed transformed dissimilarities (which are approximated by the fit)
## #'         \item confdist: the fitted distances
## #'         \item conf: the configuration to which the scaling of argument scale was applied
## #'         \item confo: the unscaled but explicitly normed configuration returned from the fitting procedure. Scaling applied to confo gives conf.
## #'         \item par, pars : the theta vector of powers tranformations (kappa,lambda,nu)
## #'         \item niter: number of iterations of the optimizer. 
## #'         \item stress: the square root of explicitly normalized stress (calculated for confo).
## #'         \item spp: stress per point
## #'         \item ndim: number of dimensions
## #'         \item model: Fitted model name with optimizer
## #'         \item call: the call
## #'         \item nobj: the number of objects
## #'         \item type, loss, losstype: stresstype
## #'         \item stress.m: The stress used for copstress. If stresstype="stress-1" this is like $stress else it is stress^2
## #'         \item stress.en: another ways to calculate the stress
## #'         \item deltaorig: the original untransformed dissimilarities  
## #'         \item copstress: the copstress loss value
## #'         \item resmat: the matrix of residuals
## #'         \item weightmat: the matrix of untransformed weights 
## #'         \item OC: the (normed) OPTICS Cordillera object (calculated for scaled conf)
## #'         \item OCv: the (normed) OPTICS Cordillera value alone (calculated for scaled conf)
## #'         \item optim: the object returned from the optimization procedure
## #'         \item stressweight, cordweight: the weights of the stress and OC respectively (v_1 and v_2)
## #'         \item optimmethod: The solver used 
## #'         \item type: the type of MDS fitted
## #'}
## #'
## #'
## #'
## #' @details
## #' Some optimizers (including the default hjk-Newuoa) will print a warning if itmax is (too) small or if there was no convergence. Consider increasing itmax then.
## #'
## #' For some solvers there also sometimes may be an error starting [smacof::transform()] which comes from the algorithm placing two object at exactly the same place so their fitted distance is 0. This is good from a OPTICS Cordillera point of view (as it is more clustered) which is why some solvers lie to pick that up, but can lead to an issue in the optimal scaling in smacof. This can usually be mitigated when specifying the model by either using less cordweight, less itmax, less accuracy or combining the two offending objects (so include them as a combined row in the distance matrix).
## #'
## #' We might eventually switch to newuoa in nloptr. 
## #' 
## #' @examples
## #' dis<-as.matrix(smacof::kinshipdelta)
## #'
## #' set.seed(1)
## #' ## Copstress with equal weight to stress and cordillera 
## #' res1<-copstressMin(dis,stressweight=0.5,cordweight=0.5,
## #'                   itmax=500) #use higher itmax about 10000 
## #' res1
## #' summary(res1)
## #' plot(res1)  #super clustered
## #'
## #' ##Alias name 
## #' res1<-copsc(dis,stressweight=0.5,
## #'                   cordweight=0.5,itmax=500) 
## #'
## #'
## #' ## Elastic scaling ordinal s-stress with cops penalty
## #' res1<-copsc(dis,type="ordinal",kappa=2,nu=-2,weightmat=dis,
## #'             stressweight=0.5, cordweight=0.5,itmax=500)
## #' 
## #' 
## #' @import cordillera
## #' @importFrom utils tail
## #' @importFrom stats dist as.dist optim sd
## #' @importFrom dfoptim hjk
## #' @importFrom NlcOptim solnl
## #' @importFrom Rsolnp solnp
## #' @importFrom subplex subplex
## #' @importFrom crs snomadr
## #' @importFrom cmaes cma_es
## #' @importFrom rgenoud genoud
## #' @importFrom GenSA GenSA
## #' @importFrom nloptr direct
## #' @importFrom minqa newuoa
## #' 
## #' 
## #' @keywords clustering multivariate
## #' @export
## copstressMinOLD <- function (delta, kappa=1, lambda=1, nu=1, theta=c(kappa,lambda,nu), type=c("ratio","interval","ordinal"), ties="primary", weightmat=1-diag(nrow(delta)),  ndim = 2, init=NULL, stressweight=0.975,cordweight=0.025,q=1,minpts=ndim+1,epsilon=max(10,max(delta)),dmax=NULL,rang,optimmethod=c("NelderMead","Newuoa","BFGS","SANN","hjk","solnl","solnp","subplex","snomadr","hjk-Newuoa","hjk-BFGS","BFGS-hjk","Newuoa-hjk","cmaes","direct","direct-Newuoa","direct-BFGS","genoud","gensa"),verbose=0,scale=c("sd","rmsq","std","proc","none"),normed=TRUE, accuracy = 1e-7, itmax = 10000, stresstype=c("stress-1","stress"),...)
## {
##     if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
##     if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
##     ## -- Setup for MDS type
##     if(missing(type)) type <- "ratio"
##     trans <- type
##     typo <- type
##     if (trans=="ratio"){
##     trans <- "none"
##     }
##     #TODO: if we want other optimal scalings as well
##     else if (trans=="ordinal" & ties=="primary"){
##     trans <- "ordinalp"
##     typo <- "ordinal (primary)"
##    } else if(trans=="ordinal" & ties=="secondary"){
##     trans <- "ordinals"
##     typo <- "ordinal (secondary)"
##   } else if(trans=="ordinal" & ties=="tertiary"){
##     trans <- "ordinalt"
##     typo <- "ordinal (tertiary)"
##   #} else if(trans=="spline"){
##   #  trans <- "mspline"
##    }
##     if(type=="interval" | type =="ordinal") theta <- c(kappa,1,nu) #We dont allow powers for dissimlarities in interval and nonmetric MDS
##     kappa <- theta[1]
##     lambda <- theta[2]
##     nu <- theta[3]
##     #plot <- FALSE

##     n <- dim(delta)[1]
##     labos <- rownames(delta)
    
##     if(verbose>0) cat(paste("Minimizing",type,"copstress with kappa=",kappa,"lambda=",lambda,"nu=",nu,".\n"))
##     if(missing(optimmethod)) optimmethod <- "hjk-Newuoa"
##     if(missing(scale)) scale <- "sd"
##     if(missing(stresstype)) stresstype <- "stress-1"

##     ##-- Prepare for dissimilarity scaling
##     r <- kappa/2
##     deltaorig <- delta #delta is diss in smacof, deltaorig is delta in smacof 
##     delta <- delta^lambda
##     weightmato <- weightmat
##     weightmat <- weightmat^nu
##     weightmat[!is.finite(weightmat)] <- 1 #new
##     deltaold <- delta
##     disobj <- smacof::transPrep(as.dist(delta), trans = trans, spline.intKnots = 2, spline.degree = 2)#spline.intKnots = spline.intKnots, spline.degree = spline.degree) #FIXME: only works with dist() style object 
##     ## Add an intercept to the spline base transformation
##     #if (trans == "mspline") disobj$base <- cbind(rep(1, nrow(disobj$base)), disobj$base)
    
##     delta <- delta / enorm (delta, weightmat) #normalize to sum=1
##     ## --- starting rang if not given
##     if(missing(rang))
##         #perhaps put this into the optimization function?
##        {
##         if(is.null(dmax))
##            {
##            if(verbose>1) cat ("Fitting configuration for rang. \n")    
##            #initsol <- cops::powerstressFast(delta,kappa=kappa,lambda=lambda,nu=nu,weightmat=weightmat,ndim=ndim)
##            initsol <- smacof::torgerson(delta,p=ndim)
##            #init0 <- initsol$conf
##            init0 <- initsol
##            init0 <- init0/enorm(init0)
##            if(scale=="std") init0 <- scale(init0) #standardizes config before cordillera
##            if(scale=="none") init0 <- init0
##            if(scale=="sd") #scales config to sd=1 for most spread dimension before cordillera
##              {
##                 init0 <- init0/max(apply(init0,2,stats::sd))
##              }   
##              if(scale=="rmsq") #scales config to rmsq=1 for most spread dimension before cordillera
##              {
##                  testso <- scale(init0,center=FALSE)
##                  init0 <- init0/max(attr(testso,"scaled:scale"))
##              }
##              if(scale=="proc") #scales config by procrusting to init
##              {
##                  if(!exists("init")) init <- initsol
##                  procr <- smacof::Procrustes(init,init0)
##                  init0 <- procr$Yhat
##              }
##            crp <- cordillera::cordillera(init0,q=q,minpts=minpts,epsilon=epsilon,scale=FALSE)$reachplot
##            cin <- max(crp)
##            dmax <- ifelse(kappa>1,1.5*cin,1.1*cin)
##            #if(verbose > 2) cat("rang, dmax was NULL or missing which makes the cordillera a goodness-of-clustering relative to the largest distance of each given configuration. \n")
##            }
##         rang <- c(0,dmax) #This sets the range to (0,dmax)
##         if(verbose>1) cat("dmax is",dmax,". rang is",rang,".\n")
##     }
##     if(isTRUE(max(rang)!=dmax)) warning("Note: The supplied dmax and rang do not match. I took supplied rang as rang.\n")

##     ## --- starting values
##     if(is.null(init))
##     {
##         if(exists("init0")) init <- init0 else init <- cops::powerStressFast(delta,kappa=kappa,lambda=lambda,nu=nu,ndim=ndim)$conf
##     }
##     xold <- init
##     xold <- xold/enorm(xold)
##     #labs <- row.names(delta)
##     copsf <- function(x,delta,disobj,r,n,ndim,weightmat,stressweight,cordweight,q,minpts,epsilon,rang,scale,normed,init,...)
##     {
##             #init is used here only for Procrustes 
##              if(!is.matrix(x)) x <- matrix(x,ncol=ndim)
##              delta <- delta/enorm(delta,weightmat)             
##              x <- x/enorm(x)
##              dnew <- sqdist (x)
##              e <- as.dist(sqrt(dnew)) #I need the dist(x) here for interval
##              #e <- dist(x) #I need the dist(x) here for interval
##              dhat <- smacof::transform(e, disobj, w = as.dist(weightmat), normq = 0.5)  ## dhat update FIXME: check if it works okay to have as.dist(weightmat) here
##              dhatt <- dhat$res #FIXME: I need the structure here to reconstruct the delta; alternatively turn all into vectors? - check how they do it in smacof
##              dhatd <- structure(dhatt, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
##                                         #FIXME: labels
             
##              delta <<- as.matrix(dhatd)
##              rnew <- sum (weightmat * delta * mkPower (dnew, r))
##              nnew <- sum (weightmat * mkPower (dnew,  2*r))
##              anew <- rnew / nnew
##              stressi <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
##              if(stresstype=="stress-1") stressi <- sqrt(stressi)
##              if(scale=="none") x <- x 
##              if(scale=="std") x <- base::scale(x) #standardizes config before cordillera
##              if(scale=="sd") #scales config to sd=1 for most spread dimension before cordillera
##              {
##                 x <- x/max(apply(x,2,stats::sd))
##              }   
##              if(scale=="rmsq") #scales config to rmsq=1 for most spread dimension before cordillera
##              {
##                  testso <- base::scale(x,center=FALSE)
##                  x <- x/max(attr(testso,"scaled:scale"))
##              }
##              if(scale=="proc") #scales config by procrusting to init
##              {
##                  procr <- smacof::Procrustes(init,x)
##                  x <- procr$Yhat
##              }
##              corrd <- cordillera::cordillera(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=FALSE,...)
##              struc <- corrd$raw
##              if(normed) {
##                         struc <- corrd$normed
##                        }
##              ic <- stressweight*stressi - cordweight*struc
##              if(verbose>1) cat("copstress =",ic,"mdsloss =",stressi,"OC =",struc,"minpts=",minpts,"kappa =",kappa,"lambda =",lambda,"nu=",nu,"\n")
##              ic
##            }
##     if(verbose>1) cat("Starting Minimization with",optimmethod,":\n")
##     if(optimmethod=="Newuoa") {
##          suppressWarnings(optimized <- minqa::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose-2),...))
##          itel <- itel+optimized$feval
##          ovalue <-optimized$fval
##          #optimized <- nloptr::newuoa(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose-2),...)
##          #xnew <- matrix(optimized$par,ncol=ndim)
##          #itel <- optimized$iter
##          #ovalue <-optimized$value
##      }
##        if(optimmethod=="direct") {
##           xold <- as.vector(xold)
##           optimized <- nloptr::direct(function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),lower=rep(5*min(xold),length(xold)),upper=rep(5*max(xold),length(xold)),nl.info=isTRUE(verbose>1),control=list(maxeval=itmax,xtol_rel=accuracy),...)
##  #        optimized <- nloptr::direct(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,#n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=),#lower=rep(5*min(xold),length(xold)),upper=rep(5*max(xold),length(xold)), nl.info=isTRUE(verbose>1),control=list(maxeval=itmax,xtol_rel=accuracy),...)
##          xnew <- matrix(optimized$par,ncol=ndim)
##          itel <- optimized$iter
##          ovalue <-optimized$value
##        }
##         if(optimmethod=="direct-Newuoa") {
##          xold <- as.vector(xold)
##          optimized1 <- nloptr::direct(function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),lower=rep(5*min(xold),length(xold)),upper=rep(5*max(xold),length(xold)),nl.info=isTRUE(verbose>1),control=list(maxeval=itmax,xtol_rel=accuracy),...)
##          xnew <- optimized1$par
##          itel1 <- optimized1$iter
##          itmaxreduced <- itmax-itel1
##          if(itel1>itmax) itmaxreduced <- 0.1*itmax
##          suppressWarnings(optimized <- minqa::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmaxreduced,rhoend=accuracy,iprint=verbose-2),...))
##          itel <- itel1+optimized$feval
##          ovalue <-optimized$fval
##          #optimized <- nloptr::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),nl.info=isTRUE(verbose>2),control=list(maxeval=itmaxreduced,xtol_rel=accuracy),...)
##          #xnew <- matrix(optimized$par,ncol=ndim)
##          #itel <- itel1+optimized$iter
##          #ovalue <-optimized$value
##          }
##          if(optimmethod=="direct-BFGS") {
##          xold <- as.vector(xold)
##          optimized1 <- nloptr::direct(function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),lower=rep(5*min(xold),length(xold)),upper=rep(5*max(xold),length(xold)),nl.info=isTRUE(verbose>1),control=list(maxeval=itmax,xtol_rel=accuracy),...)
##          xnew <- matrix(optimized1$par,ncol=ndim)
##          itel1 <- optimized1$iter
##          itmaxreduced <- itmax-itel1
##          if(itel1>itmax) itmaxreduced <- 0.1*itmax
##          optimized <- optim(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="BFGS",control=list(maxit=itmaxreduced,trace=0,reltol=accuracy),...)
##          xnew <- optimized$par
##          itel <- itel1+optimized$counts[[1]]
##          ovalue <-optimized$val
##      }
##        if(optimmethod=="genoud") {
##          optimized <- rgenoud::genoud(function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),starting.values=as.vector(xold),nvars=length(xold),print.level=verbose,...)
##          xnew <- matrix(optimized$par,ncol=ndim)
##          itel <- optimized$generations
##          ovalue <-optimized$value
##      }
##        if(optimmethod=="gensa") {
##          optimized <- GenSA::GenSA(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),lower=rep(5*min(xold),length(xold)),upper=rep(5*max(xold),length(xold)),control=list(max.call=itmax,verbose=isTRUE(verbose>1),smooth=FALSE),...)
##          xnew <- matrix(optimized$par,ncol=ndim)
##          itel <- optimized$counts
##          ovalue <-optimized$value
##      }
##     if(optimmethod=="cmaes") {
##          xold <- as.vector(xold)
##          optimized <- cmaes::cma_es(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmax),...)
##          xnew <- matrix(optimized$par,ncol=ndim)
##          itel <- optimized$counts[1]
##          ovalue <-optimized$value
##     }
##      if(optimmethod=="cmaes-Newuoa") {
##          xold <- as.vector(xold)
##          optimized1 <- cmaes::cma_es(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmax),...)
##          xnew <- matrix(optimized1$par,ncol=ndim)
##          itel1 <- optimized1$counts[1]
##          itmaxreduced <- itmax-itel1
##          if(itel1>itmax) itmaxreduced <- 0.1*itmax
##          suppressWarnings(optimized <- minqa::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmaxreduced,rhoend=accuracy,iprint=verbose-2),...))
##          xnew <- matrix(optimized$par,ncol=ndim)
##          itel <- itel1+optimized$feval
##          ovalue <-optimized$fval
##          #optimized <- nloptr::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),nl.info=isTRUE(verbose>2),control=list(maxeval=itmaxreduced,xtol_rel=accuracy),...)
##          #xnew <- matrix(optimized$par,ncol=ndim)
##          #itel <- itel1+optimized$iter
##          #ovalue <-optimized$value
##      }
##      if(optimmethod=="Newuoa-cmaes") {
##           suppressWarnings(optimized1 <- minqa::newuoa(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose-2),...))
##          #optimized1 <- nloptr::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),nl.info=isTRUE(verbose>2),control=list(maxeval=itmax,xtol_rel=accuracy),...)
##          xnew <- optimized1$par
##          itel1 <- optimized1$feval
##          itmaxreduced <- itmax-itel1
##          if(itel1>itmax) itmaxreduced <- 0.1*itmax
##          optimized <- cmaes::cma_es(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmaxreduced),...)
##          xnew <- matrix(optimized$par,ncol=ndim)
##          itel <- optimized$counts[1]+itel1
##          ovalue <-optimized$value
##      }
##      if(optimmethod=="NelderMead") {
##          optimized <- optim(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmax,trace=0,reltol=accuracy),...)
##          xnew <- optimized$par
##          itel <- optimized$counts[[1]]
##          ovalue <-optimized$val 
##      }
##     if(optimmethod=="BFGS") {
##          optimized <- optim(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="BFGS",control=list(maxit=itmax,trace=0,reltol=accuracy),...)
##          xnew <- optimized$par
##          itel <- optimized$counts[[1]]
##          ovalue <-optimized$val 
##      }
##     if(optimmethod=="SANN") {
##          optimized <- optim(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="SANN",control=list(maxit=itmax,trace=0,reltol=accuracy),...)
##          xnew <- optimized$par
##          itel <- optimized$counts[[1]]
##          ovalue <-optimized$val 
##      }
##      if(optimmethod=="hjk") {
##          optimized <- dfoptim::hjk(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfeval=itmax),...)
##          xnew <- optimized$par
##          itel <- optimized$feval
##          ovalue <-optimized$value 
##      }
##     if(optimmethod=="hjk-Newuoa") { #twostep1
##          #cat("Before HJK Iterations:", itmax,"\n")
##          optimized1 <- dfoptim::hjk(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfeval=itmax),...)
##          xnew <- optimized1$par
##          itel1 <- optimized1$feval
##          #cat("OpVal after HJK:",optimized1$value,"\n")
##          #cat("After HJK Iterations:", itel1,"\n")
##          itmaxreduced <- itmax-itel1
##          if(itel1>itmax) itmaxreduced <- 0.1*itmax
##          #cat("Before Newuoa Iterations:", itmaxreduced,"\n")
##          suppressWarnings(optimized <- minqa::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmaxreduced,rhoend=accuracy,iprint=verbose-2),...))
##          xnew <- matrix(optimized$par,ncol=ndim)
##          itel <- itel1+optimized$feval
##          ovalue <-optimized$fval
##          #optimized <- nloptr::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),nl.info=isTRUE(verbose>2),control=list(maxeval=itmaxreduced,xtol_rel=accuracy),...)
##          ##BUG: Was init= before in all Newuoa (no argument). Changed in 1.3-6
##          #xnew <- matrix(optimized$par,ncol=ndim) #
##          ##test function for return value fr <- function(x) { sum((x - matrix(c(8,5,3,4,6,7),ncol=2))^2) } if returned is matrix(c(8,5,3,4,6,7),ncol=2) it is cool
##          #cat("Newuoa Iterations:", optimized$iter,"\n")
##          #cat("Newuoa Iterations:", optimized$feval,"\n")
##          #itel <- itel1+optimized$iter
##          #cat("After Newuoa Iterations:", itel,"\n")
##          #ovalue <-optimized$value
##          #cat("OpVal after all:",ovalue,"\n")
##      }
##     if(optimmethod=="hjk-BFGS") { #twostep2
##          optimized1 <- dfoptim::hjk(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfeval=itmax,trace=0),...)
##          xnew <- optimized1$par
##          itel1 <- optimized1$feval
##          itmaxreduced <- itmax-itel1
##          if(itel1>itmax) itmaxreduced <- 0.1*itmax
##          optimized <- optim(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="BFGS",control=list(maxit=itmaxreduced,trace=0,reltol=accuracy),...)
##          xnew <- optimized$par
##          itel <- itel1+optimized$counts[[1]]
##          ovalue <-optimized$val
##      }
##      if(optimmethod=="BFGS-hjk") { #twostep6
##          optimized1 <- optim(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="BFGS",control=list(maxit=itmax,trace=0),...)
##          xnew <- optimized1$par
##          itel1 <- optimized1$counts[[1]]
##          itmaxreduced <- itmax-itel1
##          if(itel1>itmax) itmaxreduced <- 0.1*itmax
##          optimized <- dfoptim::hjk(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfeval=itmaxreduced,trace=0,tol=accuracy),...)
##          xnew <- optimized$par
##          itel <- itel1+optimized$feval
##          ovalue <-optimized$val
##      }
##      if(optimmethod=="BFGS-Newuoa") { #twostep5
##          optimized1 <- optim(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="BFGS",control=list(maxfeval=itmax,trace=0),...)
##          itel1 <- optimized1$counts[[1]]
##          xnew <- optimized1$par
##          itmaxreduced <- itmax-itel1
##          if(itel1>itmax) itmaxreduced <- 0.1*itmax
##          suppressWarnings(optimized <- minqa::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmaxreduced,rhoend=accuracy,iprint=verbose-2),...))
##          xnew <- matrix(optimized$par,ncol=ndim)
##          itel <- itel1+optimized$feval
##          ovalue <-optimized$fval
##          #optimized <- nloptr::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),nl.info=isTRUE(verbose>2),control=list(maxeval=itmaxreduced,xtol_rel=accuracy),...)
##          #xnew <- matrix(optimized$par,ncol=ndim) #
##          #itel <- itel1+optimized$iter
##          #ovalue <-optimized$value
##      }
##         if(optimmethod=="hjk-solnl") {#twostep3
##          optimized1 <- dfoptim::hjk(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfeval=itmax,trace=0),...)
##          xnew <- optimized1$par
##          itel1 <- optimized1$feval
##          itmaxreduced <- itmax-itel1
##          if(itel1>itmax) itmaxreduced <- 0.1*itmax
##          optimized <- NlcOptim::solnl(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),maxnFun=itmaxreduced,tolFun=accuracy,...)
##          xnew <- optimized$par
##          itel <- itel1+optimized$counts[[1]]
##          ovalue <-optimized$fn
##      }
##       if(optimmethod=="hjk-subplex") {#twostep4
##          optimized1 <- dfoptim::hjk(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfeval=itmax,trace=0),...)
##          xnew <- optimized1$par
##          itel1 <- optimized1$feval
##          itmaxreduced <- itmax-itel1
##          if(itel1>itmax) itmaxreduced <- 0.1*itmax
##          optimized <- subplex::subplex(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmaxreduced,reltol=accuracy),...)
##          xnew <- matrix(optimized$par,ncol=ndim)
##          itel <- itel1+optimized$count
##          ovalue <-optimized$value
##      }
##   #   if(optimmethod=="dfsane") {
##   #       optimized <- BB::dfsane(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmax,trace=0,tol=accuracy),...)
##    #      xnew <- optimized$par
##     #     itel <- optimized$feval
##     #     ovalue <-optimized$residual 
##     # }
##      if(optimmethod=="solnl") {
##          optimized <- NlcOptim::solnl(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),maxnFun=itmax,tolFun=accuracy,...)
##          xnew <- optimized$par
##          itel <- optimized$counts[[1]]
##          ovalue <-optimized$fn 
##      }
##      ## if(optimmethod=="isres") {
##      ##     optimized <- nloptr::isres(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),maxeval=itmax,trace=verbose-2,xtol_rel=accuracy,...)
##      ##     xnew <- optimized$par
##      ##     itel <- optimized$iter
##      ##     ovalue <-optimized$value 
##      ## }
##     if(optimmethod=="solnp") {
##          optimized <- Rsolnp::solnp(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(outer.iter=itmax,trace=0,tol=accuracy),...)
##          xnew <- matrix(optimized$pars,ncol=ndim)
##          itel <- optimized$nfuneval
##          ovalue <-utils::tail(optimized$values,1) 
##      }
##      if(optimmethod=="subplex") {
##          optimized <- subplex::subplex(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmax,reltol=accuracy),...)
##          xnew <- matrix(optimized$par,ncol=ndim)
##          itel <- optimized$count
##          ovalue <-optimized$value 
##      }
##     if(optimmethod=="snomadr") {
##         copsf2 <- function(x,params)
##         {
##         copsf(x,delta=params[[1]],disobj=params[[2]],r=params[[3]],n=params[[4]],ndim=params[[5]],weightmat=params[[6]],stressweight=params[[7]],cordweight=params[[8]],q=params[[9]],minpts=params[[10]],epsilon=params[[11]],rang=params[[12]],scale=params[[13]],normed=params[[14]],init=params[[15]])
##         }      
##        params <- list(delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init)
##         optimized <- crs::snomadr(copsf2,params=params,,n=dim(xold)[1],x0=xold,print.output=isTRUE(verbose-2>0),...)
##          xnew <- optimized$solution
##          itel <- optimized$iterations
##          ovalue <-optimized$objective 
##     }
##  #     if(optimmethod=="snomadr-Newuoa") {
##  #        optimized1 <- crs::snomadr(function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),n=dim(xold)[1],x0=xold,print.output=isTRUE(verbose-2>0),...)
## #         xnew <- optimized1$solution
## #         itel <- optimized1$iterations
## #         optimized <- minqa::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmax-itel,rhoend=accuracy,iprint=verbose-2),...)
## #         xnew <- matrix(optimized$par,ncol=ndim)
## #         itel <- itel+optimized$feval
## #         ovalue <-optimized$fval
## #      }
## #        if(optimmethod=="snomadr-BFGS") {
## #         optimized1 <- crs::snomadr(function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=w#eightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),n=dim(xold)[1],x0=xold,print.output=isTRUE(verbose-2>0),...)
## #         xnew <- optimized1$solution
## #         itel <- optimized1$iterations
## #        optimized <- optim(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="BFGS",control=list(maxit=itmax-itel,trace=0,reltol=accuracy),...)
## #         xnew <- optimized$par
## #         itel <- itel+optimized$counts[[1]]
## #         ovalue <-optimized$val
##                                         #     }
##      rownames(delta) <- labos
##      xnew <- xnew/enorm(xnew)
##      dnew <- sqdist (xnew)
##      rnew <- sum (weightmat * delta * mkPower (dnew, r))
##      nnew <- sum (weightmat * mkPower (dnew,  2*r))
##      anew <- rnew / nnew
##      stress.m <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
##      stress <- sqrt(stress.m)
##      if(stresstype=="stress-1") {
##          stress.m <- sqrt(stress.m)
##          }
##      attr(xnew,"dimnames")[[1]] <- rownames(delta)
##      attr(xnew,"dimnames")[[2]] <- paste("D",1:ndim,sep="")
##      #doutm <- (2*sqrt(sqdist(xnew)))^kappa  #fitted powered euclidean distance but times two
##      doutm <- as.matrix(dist(xnew)^kappa)
##      #deltam <- delta
##      #deltaorigm <- deltaorig
##      #deltaoldm <- deltaold
##      delta <- stats::as.dist(delta)
##      deltaorig <- stats::as.dist(deltaorig)
##      deltaold <- stats::as.dist(deltaold)
##      #doute <- doutm/enorm(doutm)
##      #doute <- stats::as.dist(doute)
##      dout <- doute <- stats::as.dist(doutm)
##      resmat <- as.matrix((delta - doute)^2)
##      spp <- colMeans(resmat)
##      #weightmatm <-weightmat
##      weightmat <- stats::as.dist(weightmat)
##      stressen <- sum(weightmat*(delta-doute)^2)#raw stress on the normalized proximities and normalized distances
##      if(scale=="std") xnews <- base::scale(xnew) #standardizes config before cordillera
##      if(scale=="sd") #scales config to sd=1 for most spread dimension before cordillera
##              {
##                 xnews <- xnew/max(apply(xnew,2,stats::sd))
##              }   
##              if(scale=="rmsq") #scales config to rmsq=1 for most spread dimension before cordillera
##              {
##                  testso <- base::scale(xnew,center=FALSE)
##                  xnews <- xnew/max(attr(testso,"scaled:scale"))
##              }
##              if(scale=="proc") #scales config by procrusting to init
##              {
##                  procr <- smacof::Procrustes(init,xnew)
##                  xnews <- procr$Yhat
##              }
##              if(scale=="none") xnews <- xnew #no standardisation
##       if(verbose>0) cat("*** stress (both normalized - for COPS/STOPS):",stress.m,"; stress 1 (both normalized - default reported):",stress,"; stress manual (for debug only):",stressen,"; from optimization: ",ovalue,"\n")   
##     out <- list(delta=deltaold, obsdiss=delta, confdist=dout, conf = xnews, confo=xnew, pars=c(kappa,lambda,nu), niter = itel, stress=stress, spp=spp, ndim=ndim, model=paste(typo,"copstress",optimmethod), call=match.call(), nobj = n, type = type, ties=ties, gamma=NA, stress.m=stress.m, stress.en=stressen, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat)
##     out$parameters <- out$theta <- theta
##     out$loss <- "copstress"
##     out$OC <- cordillera::cordillera(out$conf,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=FALSE)
##     out$OCv <- ifelse(normed,out$OC$normed,out$OC$raw)
##     out$copstress <- ovalue
##     out$optim <- optimized
##     out$stressweight <- stressweight
##     out$cordweight <- cordweight
##     out$optimethod <- optimmethod
##     out$losstype <- out$loss
##     out$typo <- typo
##     #out$nobj <- dim(out$conf)[1]
##     class(out) <- c("copsc","cops","smacofP","smacofB","smacof")
##     out
## }

#' Fitting a COPS-C Model (COPS Variant 1).
#'
#' Minimizing Copstress to obtain a clustered ratio, interval or ordinal PS configuration with given explicit power transformations theta. The function allows mix-and-match of explicit (via theta) and implicit (via type) transformations by setting the kappa, lambda, nu (or theta) and type arguments.
#'
#' This is an extremely flexible approach to least squares proximity scaling: It supports ratio power stress; ratio, interval and ordinal r stress and ratio, interval and ordinal MDS with or without a COPS penalty. Famous special cases of these models that can be fitted are multiscale MDS if kappa->0 and delta=log(delta), Alscal MDS (sstress) with lambda=kappa=2, sammon type mapping with weightmat=delta and nu=-1, elastic scaling with weightmat=delta and nu=-2. Due to mix-and-match this function also allows to fit models that have not yet been published, such as for example an "elastic scaling ordinal s-stress with cops penalty".
#'
#' If one wants to fit these models without the cops penalty, we recommend to use powerStressMin (for ratio MDS with any power transformation for weights, dissimilarities and distances) or rStressMin (for interval and ordinal MDS with power transformations for distances and weights) as these use majorization.  
#'
#' @rdname copstressMin
#' 
#' @param delta numeric matrix or dist object of a matrix of proximities
#' @param kappa power transformation for fitted distances
#' @param lambda power transformation for proximities (only used if type="ratio" or "interval")
#' @param nu power transformation for weights
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances if it exists), the second lambda (for the observed proximities if it exist and type="ratio" or "interval"), the third is nu (for the weights if it exists). If less than three elements are is given as argument, it will be recycled. Defaults to 1 1 1. Will override any kappa, lambda, nu parameters if they are given and do not match.
#' @param type what type of MDS to fit. Currently one of "ratio", "interval" or "ordinal". Default is "ratio".
#' @param ties the handling of ties for ordinal (nonmetric) MDS. Possible are "primary" (default), "secondary" or "tertiary".
#' @param weightmat (optional) a matrix of nonnegative weights; defaults to 1 for all off diagonals
#' @param ndim number of dimensions of the target space
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 0.975
#' @param cordweight weight to be used for the cordillera; defaults to 0.025
#' @param q the norm of the cordillera; defaults to 1
#' @param minpts the minimum points to make up a cluster in OPTICS, see [dbscan::optics()] where it is called \code{minPts}; defaults to ndim+1.
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked, see [dbscan::optics()]; defaults to 10. Note this means we do not expect any noise objects per default. This number will rarely be exceeded if we standardize the configuration as is the default in cops. However if no standardization is applied or there is a procrustes adjustment to a configuration with variance of 10 or more on any of the axes, it can have the effect of being too small. In that case just set a much higher epsilon.
#' @param dmax The winsorization limit of reachability distances in the OPTICS Cordillera. If supplied, it should be either a numeric value that matches max(rang) or NULL; if NULL it is found as 1.5 times (for kappa >1) or 1 times (for kappa <=1) the maximum reachbility value of the power torgerson model with the same lambda. If dmax and rang are supplied and dmax is not max(rang), a warning is given and rang takes precedence.   
#' @param rang range of the reachabilities to be considered. If missing it is found from the initial configuration by taking 0 as the lower boundary and dmax (see above) as upper boundary. See also \code{\link{cordillera}}     
#' @param optimmethod What optimizer to use? Choose one string of 'Newuoa' (from package minqa), 'NelderMead', 'hjk' (Hooke-Jeeves algorithm from dfoptim), 'solnl' (from nlcOptim), 'solnp' (from Rsolnp), 'subplex' (from subplex), 'SANN' (simulated annealing), 'BFGS', 'snomadr' (from crs), 'genoud' (from rgenoud), 'gensa' (from GenSA), 'cmaes' (from cmaes) and 'direct' (from nloptr). See the according R packages for details on these solvers. There are also combinations that proved to work well good, like 'hjk-Newuoa', 'hjk-BFGS', 'BFGS-hjk', 'Newuoa-hjk', 'direct-Newuoa' and 'direct-BFGS'. Usually everything with hjk, BFGS, Newuoa, subplex and solnl in it work rather well in an acceptable time frame (depending on the smoothness of copstress). Default is 'hjk-Newuoa'.   
#' @param verbose numeric value hat prints information on the fitting process; >2 is very verbose
#' @param normed should the Cordillera be normed; defaults to TRUE.
#' @param scale Scale the configuration for calculation of the OC (the scaled configuration is returned as $confs, the unscaled as $conf, so manual calculation of the OC must be done with $confs) . One of "none" (so no scaling), "sd" (configuration divided by the highest standard deviation of the columns), "proc" (procrustes adjustment to the initial fit) and "rmsq" (configuration divided by the maximum root mean square of the columns). Default is "sd".   
#' @param accuracy numerical accuracy, defaults to 1e-7.
#' @param itmax maximum number of iterations. Defaults to 10000. For the two-step algorithms if itmax is exceeded by the first solver, the second algorithm is run for at least 0.1*itmax (so overall itmax may be exceeded by a factor of 1.1).
#' @param stresstype which stress to use in the copstress. Defaults to stress-1. If anything else is set, explicitly normed stress which is (stress-1)^2 is used. Using stress-1 puts more weight on MDS fit.
#' @param principal If ‘TRUE’, principal axis transformation is applied to the final configuration.
#' @param ... additional arguments to be passed to the optimization procedure
#'
#' @return A copsc object (inheriting from smacofP). A list with the components
#'         \itemize{
#'         \item delta: the original untransformed dissimilarities
#'         \item tdelta: the explicitly transformed dissimilarities 
#'         \item dhat: the explicitly transformed dissimilarities (dhats), optimally scaled and normalized (which are approximated by the fit)
#'         \item confdist: Configuration distances, the fitted distances
#'         \item conf: the configuration (normed)
#'         \item sconf: the scaled configuration as specified in scale. Scaling applied to conf gives sconf.
#'         \item parameters, par, pars : the theta vector of powers tranformations (kappa, lambda, nu)
#'         \item niter: number of iterations of the optimizer. 
#'         \item stress: the square root of explicitly normalized stress (calculated for confo).
#'         \item spp: stress per point
#'         \item ndim: number of dimensions
#'         \item model: Fitted model name
#'         \item call: the call
#'         \item nobj: the number of objects
#'         \item type, loss, losstype: stresstype
#'         \item stress.m: The stress used for copstress. If stresstype="stress-1" this is like $stress else it is stress^2
#'         \item copstress: the copstress loss value
#'         \item resmat: the matrix of residuals
#'         \item weightmat: the matrix of untransformed weights
#'         \item tweightmat: the transformed weighting matrix (here weightmat^nu) 
#'         \item OC: the (normed) OPTICS Cordillera object (calculated for scaled conf)
#'         \item OCv: the (normed) OPTICS Cordillera value alone (calculated for scaled conf)
#'         \item optim: the object returned from the optimization procedure
#'         \item stressweight, cordweight: the weights of the stress and OC respectively (v_1 and v_2)
#'         \item optimmethod: The solver used 
#'         \item type: the type of MDS fitted
#'}
#'
#'
#'
#' @details
#' Some optimizers (including the default hjk-Newuoa) will print a warning if itmax is (too) small or if there was no convergence. Consider increasing itmax then.
#'
#' For some solvers there also sometimes may be an error starting [smacof::transform()] which comes from the algorithm placing two object at exactly the same place so their fitted distance is 0. This is good from a OPTICS Cordillera point of view (as it is more clustered) which is why some solvers lie to pick that up, but can lead to an issue in the optimal scaling in smacof. This can usually be mitigated when specifying the model by either using less cordweight, less itmax, less accuracy or combining the two offending objects (so include them as a combined row in the distance matrix).
#'
#' We might eventually switch to newuoa in nloptr. 
#' 
#' @examples
#' dis<-as.matrix(smacof::kinshipdelta)
#'
#' set.seed(1)
#' ## Copstress with equal weight to stress and cordillera 
#' res1<-copstressMin(dis,stressweight=0.5,cordweight=0.5,
#'                   itmax=500) #use higher itmax about 10000 
#' res1
#' summary(res1)
#' plot(res1)  #super clustered
#'
#' ##Alias name 
#' res1<-copsc(dis,stressweight=0.5,
#'                   cordweight=0.5,itmax=500) 
#'
#'
#' ## Elastic scaling ordinal s-stress with cops penalty
#' res1<-copsc(dis,type="ordinal",kappa=2,nu=-2,weightmat=dis,
#'             stressweight=0.5, cordweight=0.5,itmax=500)
#' 
#' 
#' @import cordillera
#' @importFrom utils tail
#' @importFrom stats dist as.dist optim sd
#' @importFrom dfoptim hjk
#' @importFrom NlcOptim solnl
#' @importFrom Rsolnp solnp
#' @importFrom subplex subplex
#' @importFrom crs snomadr
#' @importFrom cmaes cma_es
#' @importFrom rgenoud genoud
#' @importFrom GenSA GenSA
#' @importFrom nloptr direct
#' @importFrom minqa newuoa
#' 
#' 
#' @keywords clustering multivariate
#' @export
copstressMin <- function (delta, kappa=1, lambda=1, nu=1, theta=c(kappa,lambda,nu), type=c("ratio","interval","ordinal"), ties="primary", weightmat=1-diag(nrow(delta)),  ndim = 2, init=NULL, stressweight=0.975,cordweight=0.025,q=1,minpts=ndim+1,epsilon=max(10,max(delta)),dmax=NULL,rang,optimmethod=c("NelderMead","Newuoa","BFGS","SANN","hjk","solnl","solnp","subplex","snomadr","hjk-Newuoa","hjk-BFGS","BFGS-hjk","Newuoa-hjk","cmaes","direct","direct-Newuoa","direct-BFGS","genoud","gensa"),verbose=0,scale=c("sd","rmsq","proc","none"),normed=TRUE, accuracy = 1e-7, itmax = 10000, stresstype=c("stress-1","stress"),principal=FALSE,...)
{
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    ## -- Setup for MDS type
    if(missing(type)) type <- "ratio"
    type <- match.arg(type, c("ratio", "interval", "ordinal",several.ok = FALSE)) 
    trans <- type
    typo <- type
    if (trans=="ratio"){
    trans <- "none"
    }
    #TODO: if we want other optimal scalings as well
    else if (trans=="ordinal" & ties=="primary"){
    trans <- "ordinalp"
    typo <- "ordinal (primary)"
   } else if(trans=="ordinal" & ties=="secondary"){
    trans <- "ordinals"
    typo <- "ordinal (secondary)"
  } else if(trans=="ordinal" & ties=="tertiary"){
    trans <- "ordinalt"
    typo <- "ordinal (tertiary)"
  #} else if(trans=="spline"){
  #  trans <- "mspline"
   }
    if(type =="ordinal") theta <- c(kappa,1,nu) #We dont allow powers for dissimilarities in nonmetric MDS
    kappa <- theta[1]
    lambda <- theta[2]
    nu <- theta[3]
    #plot <- FALSE

    n <- nrow(delta)
    if (ndim > (n - 1)) stop("Maximum number of dimensions is n-1!")
    if(is.null(rownames(delta))) rownames(delta) <- 1:n 
    labos <- rownames(delta)
    
    if(verbose>0) cat(paste("Minimizing",type,"copstress with kappa=",kappa,"lambda=",lambda,"nu=",nu,".\n"))
    if(missing(optimmethod)) optimmethod <- "hjk-Newuoa"
    if(missing(scale)) scale <- "sd"
    if(missing(stresstype)) stresstype <- "stress-1"

    ##-- Prepare for dissimilarity scaling
    r <- kappa/2
    deltaorig <- delta 
    delta <- delta^lambda
    weightmato <- weightmat
    weightmat <- weightmat^nu
    weightmat[!is.finite(weightmat)] <- 1 #new
    deltaold <- delta
    disobj <- smacof::transPrep(as.dist(delta), trans = trans, spline.intKnots = 2, spline.degree = 2)#spline.intKnots = spline.intKnots, spline.degree = spline.degree) #FIXME: only works with dist() style object 
    ## Add an intercept to the spline base transformation
    #if (trans == "mspline") disobj$base <- cbind(rep(1, nrow(disobj$base)), disobj$base)
    delta <- delta / enorm (delta, weightmat) #normalize to sum=1
    ## --- starting rang if not given
    if(missing(rang))
        #perhaps put this into the optimization function?
       {
        if(is.null(dmax))
        {
          if(is.null(init))
           { 
           if(verbose>1) cat ("Fitting configuration for finding rang argument. \n")    
           initsol <- smacof::torgerson(delta,p=ndim)  
           init0 <- initsol
           } else init0 <- init

           init0 <- init0/enorm(init0)
          # if(scale=="std") init0 <- scale(init0) #standardizes config before cordillera
           if(scale=="none") init0 <- init0
           if(scale=="sd") #scales config to sd=1 for most spread dimension before cordillera
             {
                init0 <- init0/max(apply(init0,2,stats::sd))
             }   
             if(scale=="rmsq") #scales config to rmsq=1 for most spread dimension before cordillera
             {
                 testso <- scale(init0,center=FALSE)
                 init0 <- init0/max(attr(testso,"scaled:scale"))
             }
             if(scale=="proc") #scales config by procrusting to init
             {
                 if(is.null(init)) init <- initsol
                 procr <- smacof::Procrustes(init,init0)
                 init0 <- procr$Yhat
             }
          crp <- cordillera::cordillera(init0,q=q,minpts=minpts,epsilon=epsilon,scale=FALSE)$reachplot
          cin <- max(crp)
          dmax <- ifelse(kappa>1,1.5*cin,1.1*cin)
          }
        rang <- c(0,dmax) #This sets the range to (0,dmax)
        if(verbose>1) cat("rang is",rang,".\n")
    }
    if(isTRUE(max(rang)!=dmax)) warning("Note: The supplied dmax and rang do not match. I took supplied rang as rang.\n")

    ## --- starting values
    if(is.null(init))
    {
        if(exists("init0")) init <- init0 else init <- smacof::torgerson(delta,p=ndim)
        #was: cops::powerStressFast(delta,kappa=kappa,lambda=lambda,nu=nu,ndim=ndim)$conf
    }
    xstart <- xold <- init
    xold <- xold/enorm(xold)
    dhat2 <- NA
    #labs <- row.names(delta)
    copsf <- function(x,delta,disobj,r,n,ndim,weightmat,stressweight,cordweight,q,minpts,epsilon,rang,scale,normed,init,...)
    {
            #init is used here only for Procrustes 
             if(!is.matrix(x)) x <- matrix(x,ncol=ndim)
             delta <- delta/enorm(delta,weightmat)             
             x <- x/enorm(x)
             dnew <- sqdist (x)
             e <- as.dist(sqrt(dnew)) #I need the dist(x) here for interval
             #e <- dist(x) #I need the dist(x) here for interval
             dhat2 <<- smacof::transform(e, disobj, w = as.dist(weightmat), normq = 0.5)  ##I use <<- to change this also in the parent environment because the copsf function should only return a scalar for the optimizers but I the last dhats2 and also delta later on
             dhatt <- dhat2$res #FIXME: I need the structure here to reconstruct the delta; alternatively turn all into vectors? - check how they do it in smacof
             dhatd <- structure(dhatt, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)   
             delta <<- as.matrix(dhatd) ##I use <<- to save this also in the parent environment
             rnew <- sum (weightmat * delta * mkPower (dnew, r))
             nnew <- sum (weightmat * mkPower (dnew,  2*r))
             anew <- rnew / nnew
             stressi <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
             if(stresstype=="stress-1") stressi <- sqrt(stressi)
             if(scale=="none") x <- x 
             #if(scale=="std") x <- base::scale(x) #standardizes config before cordillera
             if(scale=="sd") #scales config to sd=1 for most spread dimension before cordillera
             {
                x <- x/max(apply(x,2,stats::sd))
             }   
             if(scale=="rmsq") #scales config to rmsq=1 for most spread dimension before cordillera
             {
                 testso <- base::scale(x,center=FALSE)
                 x <- x/max(attr(testso,"scaled:scale"))
             }
             if(scale=="proc") #scales config by procrusting to init
             {
                 procr <- smacof::Procrustes(init,x)
                 x <- procr$Yhat
             }
             corrd <- cordillera::cordillera(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=FALSE,...)
             struc <- corrd$raw
             if(normed) {
                        struc <- corrd$normed
                       }
             ic <- stressweight*stressi - cordweight*struc
             if(verbose>1) cat("copstress =",ic,"mdsloss =",stressi,"OC =",struc,"minpts=",minpts,"kappa =",kappa,"lambda =",lambda,"nu=",nu,"\n")
             ic
           }
    if(verbose>1) cat("Starting Minimization with",optimmethod,":\n")
    if(optimmethod=="Newuoa") {
         suppressWarnings(optimized <- minqa::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose-2),...))
         itel <- itel+optimized$feval
         ovalue <-optimized$fval
         #optimized <- nloptr::newuoa(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose-2),...)
         #xnew <- matrix(optimized$par,ncol=ndim)
         #itel <- optimized$iter
         #ovalue <-optimized$value
     }
       if(optimmethod=="direct") {
          xold <- as.vector(xold)
          optimized <- nloptr::direct(function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),lower=rep(5*min(xold),length(xold)),upper=rep(5*max(xold),length(xold)),nl.info=isTRUE(verbose>1),control=list(maxeval=itmax,xtol_rel=accuracy),...)
 #        optimized <- nloptr::direct(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,#n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=),#lower=rep(5*min(xold),length(xold)),upper=rep(5*max(xold),length(xold)), nl.info=isTRUE(verbose>1),control=list(maxeval=itmax,xtol_rel=accuracy),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- optimized$iter
         ovalue <-optimized$value
       }
        if(optimmethod=="direct-Newuoa") {
         xold <- as.vector(xold)
         optimized1 <- nloptr::direct(function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),lower=rep(5*min(xold),length(xold)),upper=rep(5*max(xold),length(xold)),nl.info=isTRUE(verbose>1),control=list(maxeval=itmax,xtol_rel=accuracy),...)
         xnew <- optimized1$par
         itel1 <- optimized1$iter
         itmaxreduced <- itmax-itel1
         if(itel1>itmax) itmaxreduced <- 0.1*itmax
         suppressWarnings(optimized <- minqa::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmaxreduced,rhoend=accuracy,iprint=verbose-2),...))
         itel <- itel1+optimized$feval
         ovalue <-optimized$fval
         #optimized <- nloptr::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),nl.info=isTRUE(verbose>2),control=list(maxeval=itmaxreduced,xtol_rel=accuracy),...)
         #xnew <- matrix(optimized$par,ncol=ndim)
         #itel <- itel1+optimized$iter
         #ovalue <-optimized$value
         }
         if(optimmethod=="direct-BFGS") {
         xold <- as.vector(xold)
         optimized1 <- nloptr::direct(function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),lower=rep(5*min(xold),length(xold)),upper=rep(5*max(xold),length(xold)),nl.info=isTRUE(verbose>1),control=list(maxeval=itmax,xtol_rel=accuracy),...)
         xnew <- matrix(optimized1$par,ncol=ndim)
         itel1 <- optimized1$iter
         itmaxreduced <- itmax-itel1
         if(itel1>itmax) itmaxreduced <- 0.1*itmax
         optimized <- optim(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="BFGS",control=list(maxit=itmaxreduced,trace=0,reltol=accuracy),...)
         xnew <- optimized$par
         itel <- itel1+optimized$counts[[1]]
         ovalue <-optimized$val
     }
       if(optimmethod=="genoud") {
         optimized <- rgenoud::genoud(function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),starting.values=as.vector(xold),nvars=length(xold),print.level=verbose,...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- optimized$generations
         ovalue <-optimized$value
     }
       if(optimmethod=="gensa") {
         optimized <- GenSA::GenSA(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),lower=rep(5*min(xold),length(xold)),upper=rep(5*max(xold),length(xold)),control=list(max.call=itmax,verbose=isTRUE(verbose>1),smooth=FALSE),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- optimized$counts
         ovalue <-optimized$value
     }
    if(optimmethod=="cmaes") {
         xold <- as.vector(xold)
         optimized <- cmaes::cma_es(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmax),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- optimized$counts[1]
         ovalue <-optimized$value
    }
     if(optimmethod=="cmaes-Newuoa") {
         xold <- as.vector(xold)
         optimized1 <- cmaes::cma_es(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmax),...)
         xnew <- matrix(optimized1$par,ncol=ndim)
         itel1 <- optimized1$counts[1]
         itmaxreduced <- itmax-itel1
         if(itel1>itmax) itmaxreduced <- 0.1*itmax
         suppressWarnings(optimized <- minqa::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmaxreduced,rhoend=accuracy,iprint=verbose-2),...))
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- itel1+optimized$feval
         ovalue <-optimized$fval
         #optimized <- nloptr::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),nl.info=isTRUE(verbose>2),control=list(maxeval=itmaxreduced,xtol_rel=accuracy),...)
         #xnew <- matrix(optimized$par,ncol=ndim)
         #itel <- itel1+optimized$iter
         #ovalue <-optimized$value
     }
     if(optimmethod=="Newuoa-cmaes") {
          suppressWarnings(optimized1 <- minqa::newuoa(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose-2),...))
         #optimized1 <- nloptr::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),nl.info=isTRUE(verbose>2),control=list(maxeval=itmax,xtol_rel=accuracy),...)
         xnew <- optimized1$par
         itel1 <- optimized1$feval
         itmaxreduced <- itmax-itel1
         if(itel1>itmax) itmaxreduced <- 0.1*itmax
         optimized <- cmaes::cma_es(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmaxreduced),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- optimized$counts[1]+itel1
         ovalue <-optimized$value
     }
     if(optimmethod=="NelderMead") {
         optimized <- optim(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmax,trace=0,reltol=accuracy),...)
         xnew <- optimized$par
         itel <- optimized$counts[[1]]
         ovalue <-optimized$val 
     }
    if(optimmethod=="BFGS") {
         optimized <- optim(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="BFGS",control=list(maxit=itmax,trace=0,reltol=accuracy),...)
         xnew <- optimized$par
         itel <- optimized$counts[[1]]
         ovalue <-optimized$val 
     }
    if(optimmethod=="SANN") {
         optimized <- optim(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="SANN",control=list(maxit=itmax,trace=0,reltol=accuracy),...)
         xnew <- optimized$par
         itel <- optimized$counts[[1]]
         ovalue <-optimized$val 
     }
     if(optimmethod=="hjk") {
         optimized <- dfoptim::hjk(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfeval=itmax),...)
         xnew <- optimized$par
         itel <- optimized$feval
         ovalue <-optimized$value 
     }
    if(optimmethod=="hjk-Newuoa") { #twostep1
         #cat("Before HJK Iterations:", itmax,"\n")
         optimized1 <- dfoptim::hjk(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfeval=itmax),...)
         xnew <- optimized1$par
         itel1 <- optimized1$feval
         #cat("OpVal after HJK:",optimized1$value,"\n")
         #cat("After HJK Iterations:", itel1,"\n")
         itmaxreduced <- itmax-itel1
         if(itel1>itmax) itmaxreduced <- 0.1*itmax
         #cat("Before Newuoa Iterations:", itmaxreduced,"\n")
         suppressWarnings(optimized <- minqa::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmaxreduced,rhoend=accuracy,iprint=verbose-2),...))
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- itel1+optimized$feval
         ovalue <-optimized$fval
         #optimized <- nloptr::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),nl.info=isTRUE(verbose>2),control=list(maxeval=itmaxreduced,xtol_rel=accuracy),...)
         ##BUG: Was init= before in all Newuoa (no argument). Changed in 1.3-6
         #xnew <- matrix(optimized$par,ncol=ndim) #
         ##test function for return value fr <- function(x) { sum((x - matrix(c(8,5,3,4,6,7),ncol=2))^2) } if returned is matrix(c(8,5,3,4,6,7),ncol=2) it is cool
         #cat("Newuoa Iterations:", optimized$iter,"\n")
         #cat("Newuoa Iterations:", optimized$feval,"\n")
         #itel <- itel1+optimized$iter
         #cat("After Newuoa Iterations:", itel,"\n")
         #ovalue <-optimized$value
         #cat("OpVal after all:",ovalue,"\n")
     }
    if(optimmethod=="hjk-BFGS") { #twostep2
         optimized1 <- dfoptim::hjk(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfeval=itmax,trace=0),...)
         xnew <- optimized1$par
         itel1 <- optimized1$feval
         itmaxreduced <- itmax-itel1
         if(itel1>itmax) itmaxreduced <- 0.1*itmax
         optimized <- optim(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="BFGS",control=list(maxit=itmaxreduced,trace=0,reltol=accuracy),...)
         xnew <- optimized$par
         itel <- itel1+optimized$counts[[1]]
         ovalue <-optimized$val
     }
     if(optimmethod=="BFGS-hjk") { #twostep6
         optimized1 <- optim(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="BFGS",control=list(maxit=itmax,trace=0),...)
         xnew <- optimized1$par
         itel1 <- optimized1$counts[[1]]
         itmaxreduced <- itmax-itel1
         if(itel1>itmax) itmaxreduced <- 0.1*itmax
         optimized <- dfoptim::hjk(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfeval=itmaxreduced,trace=0,tol=accuracy),...)
         xnew <- optimized$par
         itel <- itel1+optimized$feval
         ovalue <-optimized$val
     }
     if(optimmethod=="BFGS-Newuoa") { #twostep5
         optimized1 <- optim(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="BFGS",control=list(maxfeval=itmax,trace=0),...)
         itel1 <- optimized1$counts[[1]]
         xnew <- optimized1$par
         itmaxreduced <- itmax-itel1
         if(itel1>itmax) itmaxreduced <- 0.1*itmax
         suppressWarnings(optimized <- minqa::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmaxreduced,rhoend=accuracy,iprint=verbose-2),...))
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- itel1+optimized$feval
         ovalue <-optimized$fval
         #optimized <- nloptr::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),nl.info=isTRUE(verbose>2),control=list(maxeval=itmaxreduced,xtol_rel=accuracy),...)
         #xnew <- matrix(optimized$par,ncol=ndim) #
         #itel <- itel1+optimized$iter
         #ovalue <-optimized$value
     }
        if(optimmethod=="hjk-solnl") {#twostep3
         optimized1 <- dfoptim::hjk(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfeval=itmax,trace=0),...)
         xnew <- optimized1$par
         itel1 <- optimized1$feval
         itmaxreduced <- itmax-itel1
         if(itel1>itmax) itmaxreduced <- 0.1*itmax
         optimized <- NlcOptim::solnl(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),maxnFun=itmaxreduced,tolFun=accuracy,...)
         xnew <- optimized$par
         itel <- itel1+optimized$counts[[1]]
         ovalue <-optimized$fn
     }
      if(optimmethod=="hjk-subplex") {#twostep4
         optimized1 <- dfoptim::hjk(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfeval=itmax,trace=0),...)
         xnew <- optimized1$par
         itel1 <- optimized1$feval
         itmaxreduced <- itmax-itel1
         if(itel1>itmax) itmaxreduced <- 0.1*itmax
         optimized <- subplex::subplex(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmaxreduced,reltol=accuracy),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- itel1+optimized$count
         ovalue <-optimized$value
     }
  #   if(optimmethod=="dfsane") {
  #       optimized <- BB::dfsane(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmax,trace=0,tol=accuracy),...)
   #      xnew <- optimized$par
    #     itel <- optimized$feval
    #     ovalue <-optimized$residual 
    # }
     if(optimmethod=="solnl") {
         optimized <- NlcOptim::solnl(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),maxnFun=itmax,tolFun=accuracy,...)
         xnew <- optimized$par
         itel <- optimized$counts[[1]]
         ovalue <-optimized$fn 
     }
     ## if(optimmethod=="isres") {
     ##     optimized <- nloptr::isres(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),maxeval=itmax,trace=verbose-2,xtol_rel=accuracy,...)
     ##     xnew <- optimized$par
     ##     itel <- optimized$iter
     ##     ovalue <-optimized$value 
     ## }
    if(optimmethod=="solnp") {
         optimized <- Rsolnp::solnp(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(outer.iter=itmax,trace=0,tol=accuracy),...)
         xnew <- matrix(optimized$pars,ncol=ndim)
         itel <- optimized$nfuneval
         ovalue <-utils::tail(optimized$values,1) 
     }
     if(optimmethod=="subplex") {
         optimized <- subplex::subplex(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxit=itmax,reltol=accuracy),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- optimized$count
         ovalue <-optimized$value 
     }
    if(optimmethod=="snomadr") {
        copsf2 <- function(x,params)
        {
        copsf(x,delta=params[[1]],disobj=params[[2]],r=params[[3]],n=params[[4]],ndim=params[[5]],weightmat=params[[6]],stressweight=params[[7]],cordweight=params[[8]],q=params[[9]],minpts=params[[10]],epsilon=params[[11]],rang=params[[12]],scale=params[[13]],normed=params[[14]],init=params[[15]])
        }      
       params <- list(delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init)
        optimized <- crs::snomadr(copsf2,params=params,,n=dim(xold)[1],x0=xold,print.output=isTRUE(verbose-2>0),...)
         xnew <- optimized$solution
         itel <- optimized$iterations
         ovalue <-optimized$objective 
    }
 #     if(optimmethod=="snomadr-Newuoa") {
 #        optimized1 <- crs::snomadr(function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),n=dim(xold)[1],x0=xold,print.output=isTRUE(verbose-2>0),...)
#         xnew <- optimized1$solution
#         itel <- optimized1$iterations
#         optimized <- minqa::newuoa(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),control=list(maxfun=itmax-itel,rhoend=accuracy,iprint=verbose-2),...)
#         xnew <- matrix(optimized$par,ncol=ndim)
#         itel <- itel+optimized$feval
#         ovalue <-optimized$fval
#      }
#        if(optimmethod=="snomadr-BFGS") {
#         optimized1 <- crs::snomadr(function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=w#eightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),n=dim(xold)[1],x0=xold,print.output=isTRUE(verbose-2>0),...)
#         xnew <- optimized1$solution
#         itel <- optimized1$iterations
#        optimized <- optim(xnew,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,normed=normed,init=init),method="BFGS",control=list(maxit=itmax-itel,trace=0,reltol=accuracy),...)
#         xnew <- optimized$par
#         itel <- itel+optimized$counts[[1]]
#         ovalue <-optimized$val
                                        #     }
     xnew <- xnew/enorm(xnew)
     rownames(delta) <- labos
     dnew <- sqdist (xnew)
     rnew <- sum (weightmat * delta * mkPower (dnew, r))
     nnew <- sum (weightmat * mkPower (dnew,  2*r))
     anew <- rnew / nnew
     stress.m  <- snew <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
     stress <- sqrt(stress.m)
     if(stresstype=="stress-1") {
         stress.m <- sqrt(stress.m)
     }
     attr(xnew,"dimnames")[[1]] <- rownames(delta)
     attr(xnew,"dimnames")[[2]] <- paste("D",1:ndim,sep="")
     #doutm <- (2*sqrt(sqdist(xnew)))^kappa  #fitted powered euclidean distance but times two
     #doutm <- as.matrix(dist(xnew)^kappa)
     doutm <- mkPower(sqdist(xnew),r)
     #deltam <- delta
     #deltaorigm <- deltaorig
     #deltaoldm <- deltaold
     delta <- stats::as.dist(delta)
     deltaorig <- stats::as.dist(deltaorig)
     deltaold <- stats::as.dist(deltaold)
     #doute <- doutm/enorm(doutm)
     #doute <- stats::as.dist(doute)
     dout <- stats::as.dist(doutm)
     #resmat <- as.matrix((delta - doute)^2)
     #spp <- colMeans(resmat)
     #weightmatm <-weightmat
     weightmat <- stats::as.dist(weightmat)
     spoint <- spp(delta, dout, weightmat)
     resmat<-spoint$resmat
     rss <- sum(spoint$resmat[lower.tri(spoint$resmat)])
     spp <- spoint$spp
     #stressen <- sum(weightmat*(delta-doute)^2)#raw stress on the normalized proximities and normalized distances
     #if(scale=="std") xnews <- base::scale(xnew) #standardizes config before cordillera
     if(scale=="sd") #scales config to sd=1 for most spread dimension before cordillera
             {
                xnews <- xnew/max(apply(xnew,2,stats::sd))
             }   
             if(scale=="rmsq") #scales config to rmsq=1 for most spread dimension before cordillera
             {
                 testso <- base::scale(xnew,center=FALSE)
                 xnews <- xnew/max(attr(testso,"scaled:scale"))
             }
             if(scale=="proc") #scales config by procrusting to init
             {
                 procr <- smacof::Procrustes(init,xnew)
                 xnews <- procr$Yhat
             }
    if(scale=="none") xnews <- xnew #no standardisation
    if (principal) {
        xnew_svd <- svd(xnew)
        xnew <- xnew %*% xnew_svd$v
        xnews_svd <- svd(xnews)
        xnews <- xnews %*% xnews_svd$v
    }
    if(verbose>0) cat("*** Stress:",stress.m,"; Stress-1:",stress,"; from optimization: ",ovalue,"\n")
    #like smacofP
    out <- list(delta=deltaorig, tdelta=deltaold, dhat=delta, confdist=dout, iord=dhat2$iord.prim, conf = xnew, stress=stress, spp=spp, ndim=ndim, weightmat=weightmato, resmat=resmat, rss=rss, init=xstart, model="COPS-C", niter = itel, nobj = n, type = type, call=match.call(), stress.m=snew, alpha = anew, sigma = snew, pars=c(kappa=kappa,lambda=lambda,nu=nu),parameters=c(kappa=kappa,lambda=lambda,nu=nu),theta=c(kappa=kappa,lambda=lambda,nu=nu),tweightmat=weightmat)
    #extra slots for class copsc
    out$gamma <- NA 
    out$ties <- ties 
    out$sconf <- xnews 
    #out$loss <- "copstress"
    out$OC <- cordillera::cordillera(out$conf,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=FALSE)
    out$OCv <- ifelse(normed,out$OC$normed,out$OC$raw)
    out$copstress <- ovalue
    out$optim <- optimized
    out$stressweight <- stressweight
    out$cordweight <- cordweight
    out$optimethod <- optimmethod
    #out$typo <- typo
    #out$nobj <- dim(out$conf)[1]
    class(out) <- c("copsc","cops","smacofP","smacofB","smacof")
    out
}

#' @rdname copstressMin
#' @export 
copsc <- copstressMin

#' @rdname copstressMin
#' @export 
copStressMin <- copstressMin


#'@export
print.copsc <- function(x,...)
    {
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    cat("Model:",x$type,"COPS-C with parameter vector =",x$parameters,"\n")
    cat("\n")
    cat("Number of objects:", x$nobj, "\n")
    cat("Stress-1 value of configuration:", round(x$stress,5), "\n")
    cat("OPTICS Cordillera: Raw", round(x$OC$raw,5),"Normed", round(x$OC$normed,5),"\n")
    cat("Cluster optimized loss (copstress): ", round(x$copstress,5), "\n")
    cat("Stress weight:", x$stressweight," OPTICS Cordillera weight:",x$cordweight,"\n")
    cat("Number of iterations of",x$optimethod,"optimization:", x$niter, "\n")
    cat("\n")
    }


######################### FOR PARALLELIZATION EVENTUALLY #################################
#' Power Stress SMACOF
#'
#' An implementation to minimize power stress by minimization-majorization. Usually more accurate but slower than powerStressFast. Uses for loops.
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param kappa power of the transformation of the fitted distances; defaults to 1
#' @param lambda the power of the transformation of the proximities; defaults to 1
#' @param nu the power of the transformation for weightmat; defaults to 1 
#' @param weightmat a matrix of finite weights
#' @param init starting configuration
#' @param ndim dimension of the configuration; defaults to 2
#' @param acc numeric accuracy of the iteration
#' @param itmax maximum number of iterations
#' @param verbose should iteration output be printed; if > 1 then yes
#'
#' @return a smacofP object (inheriting form smacofB, see \code{\link{smacofSym}}). It is a list with the components
#' \itemize{
#' \item delta: Observed dissimilarities, not normalized
#' \item obsdiss: Observed dissimilarities, normalized 
#' \item confdiss: Configuration dissimilarities, NOT normalized 
#' \item conf: Matrix of fitted configuration, NOT normalized
#' \item stress: Default stress  (stress 1; sqrt of explicitly normalized stress)
#' \item spp: Stress per point (based on stress.en) 
#' \item ndim: Number of dimensions
#' \item model: Name of smacof model
#' \item niter: Number of iterations
#' \item nobj: Number of objects
#' \item type: Type of MDS model
#' }
#' and some additional components
#' \itemize{
#' \item stress.m: default stress for the COPS and STOP defaults to the explicitly normalized stress on the normalized, transformed dissimilarities
#' \item stress.en: a manually calculated stress on the normalized, transformed dissimilarities and normalized transformed distances which is not correct
#' \item deltaorig: observed, untransformed dissimilarities
#' \item weightmat: weighting matrix 
#'}
#'
#' @importFrom stats dist as.dist
#' 
#' @seealso \code{\link{smacofSym}}
#' 
#' 
#' @export
powerStressMin2 <- function (delta, kappa=1, lambda=1, nu=1, weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc= 1e-10, itmax = 100000, verbose = FALSE) {
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    if(verbose>0) cat("Minimizing powerStress with kappa=",kappa,"lambda=",lambda,"nu=",nu,"\n")
    r <- kappa/2
    p <- ndim
    deltaorig <- delta
    delta <- delta^lambda
    #weightmato <- weightmat
    weightmat <- weightmat^nu
    weightmat[!is.finite(weightmat)] <- 1
    deltaold <- delta
    delta <- delta / enorm (delta, weightmat)
    itel <- 1
    xold <- init
    if(is.null(init)) xold <- cops::torgerson (delta, p = p)
    xold <- xold / enorm (xold)
    n <- nrow (xold)
    nn <- diag (n)
    dold <- sqdist (xold)
    rold <- sum (weightmat * delta * mkPower (dold, r))
    nold <- sum (weightmat * mkPower (dold, 2 * r))
    aold <- rold / nold
    sold <- 1 - 2 * aold * rold + (aold ^ 2) * nold
    insideloop <- function(dold,itmax,weightmat,delta,r,aold,nn,xold,sold,verbose,itel,acc)
            {           
    for(i in seq(1,itmax,by=1))
      {
      p1 <- mkPower (dold, r - 1)
      p2 <- mkPower (dold, (2 * r) - 1)
      by <- mkBmat (weightmat * delta * p1)
      cy <- mkBmat (weightmat * p2)
      ga <- 2 * sum (weightmat * p2)
      be <- (2 * r - 1) * (2 ^ r) * sum (weightmat * delta)
      de <- (4 * r - 1) * (4 ^ r) * sum (weightmat)
      if (r >= 0.5) {
        my <- by - aold * (cy - de * nn)
      }
      if (r < 0.5) {
        my <- (by - be * nn) - aold * (cy - ga * nn)
      }
      xnew <- my %*% xold
      xnew <- xnew / enorm (xnew)
      dnew <- sqdist (xnew)
      rnew <- sum (weightmat * delta * mkPower (dnew, r))
      nnew <- sum (weightmat * mkPower (dnew, 2 * r))
      anew <- rnew / nnew
      snew <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
      if(is.na(snew)) #if there are issues with the values
          {
              snew <- sold
              dnew <- dold
              anew <- aold
              xnew <- xold
          }   
      if (verbose>2) {
        cat (
          formatC (itel, width = 4, format = "d"),
          formatC (
            sold,
            digits = 10,
            width = 13,
            format = "f"
          ),
          formatC (
            snew,
            digits = 10,
            width = 13,
            format = "f"
          ),
          "\n"
        )
      }
      if ((sold - snew) < acc)
        break ()
      itel <- i
      xold <- xnew
      dold <- dnew
      sold <- snew
      aold <- anew
     }
     list(itel=itel,xnew=xnew,snew=snew,anew=anew,dnew=dnew)
     }
     islr<-insideloop(dold,itmax,weightmat,delta,r,aold,nn,xold,sold,verbose,itel,acc)
     xnew <- islr$xnew
     dnew <- islr$dnew
     anew <- islr$anew
     snew <- islr$snew
     itel <- islr$itel
     attr(xnew,"dimnames")[[2]] <- paste("D",1:p,sep="")
     xnew <- xnew/enorm(xnew)
     doutm <- mkPower(sqdist(xnew),r)
     deltam <- delta
     delta <- stats::as.dist(delta)
     deltaorig <- stats::as.dist(deltaorig)
     deltaold <- stats::as.dist(deltaold)
     doute <- doutm/enorm(doutm) #this is an issue here!
     doute <- stats::as.dist(doute)
     dout <- stats::as.dist(doutm)
     weightmatm <-weightmat
     resmat <- weightmatm*as.matrix((delta - doute)^2) #BUG
     spp <- colMeans(resmat) #BUG
     weightmat <- stats::as.dist(weightmatm)
     stressen <- sum(weightmat*(doute-delta)^2)
     if(verbose>1) cat("*** stress (both normalized):",snew, "; stress 1 (both normalized - default reported):",sqrt(snew),"; manual stress (only for debug):",stressen, "\n")  
    out <- list(delta=deltaold, obsdiss=delta, confdiss=dout, conf = xnew, pars=c(kappa,lambda,nu), niter = itel, spp=spp, ndim=p, model="Power Stress SMACOF", call=match.call(), nobj = dim(xnew)[1], type = "Power Stress", stress=sqrt(snew), stress.m=snew,stress.en=stressen, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat, alpha = anew, sigma = snew)
    class(out) <- c("smacofP","smacofB","smacof")
    out
  }


########################### COPS-0 ################################################
                                        #TODO: Some exampes code to be transferred to cops high level


#' Fitting a COPS-0 Model by shrinking residuals to Zero (COPS-0).
#'
#' Minimizing copstress by shrinking residulas to zero to achieve a clustered Power Stress MDS configuration with given hyperparameters theta.
#'
#' @param delta numeric matrix or dist object of a matrix of proximities
#' @param kappa power transformation for fitted distances
#' @param lambda power transformation for proximities
#' @param nu power transformation for weights
#' @param theta the theta vector of powers; the first is kappa (for the fitted distances if it exists), the second lambda (for the observed proximities if it exist), the third is nu (for the weights if it exists) . If less than three elements are is given as argument, it will be recycled. Defaults to 1 1 1. Will override any kappa, lmabda, nu parameters if they are given and do not match
#' @param weightmat (optional) a matrix of nonnegative weights; defaults to 1 for all off diagonals
#' @param ndim number of dimensions of the target space
#' @param init (optional) initial configuration
#' @param cordweight weight to be used for the shrinkage; defaults to 1
#' @param q used in cordillera and shrink matrix and controls the effect of using the norm; defaults to 2 (least squares MDS)
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to ndim+1
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the minimum reachabilities to be considered. If missing it is found from the initial configuration by taking 1.5 times the maximal minimum reachability of the initial fit. If NULL it will be normed to each configuration's minimum and maximum distance, so an absolute value of goodness-of-clusteredness. Note that the latter is not necessarily desirable when comparing configurations for their relative clusteredness. See also \code{\link{cordillera}}
#' @param scaleX should X be scaled; defaults to TRUE
#' @param enormX should X be enormed; defaults to FALSE
#' @param scaleB should X be scaled for the shrink matrix; defaults to TRUE.
#' @param scaleC should X be scaled for the OPTICS Cordillera; defaults to TRUE. These parameter lets one tweak the way the shrinkage works and how its quantified; the defaults lead usually to a sensible result. It might be that some scale versions weill be depreciated in future versions. 
#' @param optimmethod What optimizer to use? Defaults to NEWUOA, Nelder-Mead is also supported.
#' @param verbose numeric value hat prints information on the fitting process; >2 is very verbose
#' @param accuracy numerical accuracy, defaults to 1e-8
#' @param itmax maximum number of iterations. Defaults to 100000
#' @param normed should cordillera and shrink matrix be normed; defaults to TRUE
#' @param ... additional arguments to be passed to the optimization procedure
#'
#' @return A list with the components
#'         \itemize{
#'         \item copstress: the weighted loss value
#'         \item OC: the Optics cordillera
#'         \item optim: the object returned from the optimization procedure
#'         \item stress: the stress
#'         \item stress.m: default normalized stress
#'         \item parameters: the parameters used for fitting (kappa, lambda)
#'         \item fit: the returned object of the fitting procedure
#'         \item cordillera: the OPTICS cordillera object
#' }
#' 
#' @examples
#' dis<-as.matrix(smacof::kinshipdelta)
#'
#' #Copstress with shrinkage to 0 
#' res1<-shrinkCopstress0(dis,cordweight=1,minpts=2) 
#' res1
#' summary(res1)
#' plot(res1)  #super clustered
#'
#' @import cordillera
#' @importFrom stats dist as.dist optim
#' @importFrom minqa newuoa
#'
#' #'#COPS-0 to improve over an MDS result
#'res0<-powerStressFast(dis)
#'res2a<-cops(dis,variant="COPS-0",cordweight=1,q=2,init=res0$conf,minpts=2) 
#'res2a
#'summary(res2a)
#'plot(res2a)

#'
#' #'#Sammon stress type copstress
#'ws<-1/dis
#'diag(ws)<-1 
#'res4<-cops(dis,variant="COPS-0",nu=-1,weightmat=ws,cordweight=0.5) 
#'res4
#'summary(res4)
#'plot(res4)
#' @keywords clustering multivariate
#' @export
shrinkCopstress0 <- function (delta, kappa=1, lambda=1, nu=1, theta=c(kappa,lambda,nu),weightmat=1-diag(nrow(delta)),  ndim = 2, init=NULL,cordweight=1,q=2,minpts=ndim+1,epsilon=10,rang=NULL,optimmethod=c("Nelder-Mead","Newuoa"),verbose=0,scaleX=TRUE,enormX=FALSE,scaleB=TRUE,scaleC=TRUE,accuracy = 1e-7, itmax = 100000,normed=2,...)
{
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    kappa <- theta[1]
    lambda <- theta[2]
    nu <- theta[3]
    plot <- FALSE
    if(verbose>0) cat("Minimizing copstress with kappa=",kappa,"lambda=",lambda,"nu=",nu,"\n")
    if(missing(optimmethod)) optimmethod <- "Newuoa"
    if(missing(rang))
        #perhaps put this into the optimization function?
          {
           if(verbose>1) cat ("Fitting configuration for rang. \n")    
           initsol <- cops::powerStressFast(delta,kappa=kappa,lambda=lambda,nu=nu,weightmat=weightmat,ndim=ndim)
           init0 <- initsol$conf
           if(isTRUE(scaleX)) init0 <- scale(init0)
           crp <- cordillera::cordillera(init0,q=q,minpts=minpts,epsilon=epsilon,scale=scaleC)$reachplot
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
    if(is.null(init)) xold <- cops::powerStressFast(delta,kappa=kappa,lambda=lambda,nu=nu,ndim=ndim)$conf
    if(enormX) xold <- xold/enorm(xold)
    if(scaleX) xold <- scale(xold)
    shrinkcops <- function(x,delta,r,ndim,weightmat,cordweight,q,minpts,epsilon,rang,scaleX,enormX,scaleB,normed,...)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=ndim)
             #try these variants again with kinship and cali:
             #it all about how to normalize so that the shrinkage will play its part
             #take care of r and so if sqrt(sqdist()) use r*2
             #delta enormed, x scaled + enormed; looks good! -> looks best? 
             #delta enormed, x scaled, dnew enormed; looks ok like #2 but a bit better
             #delta enormed, x enormed, dnew normal; looks ok with clusters for kinship but wrong clusters; closest snew and mdsloss
             if(scaleX) x <- scale(x)
             if(enormX) x <- x/enorm(x)
             delta <- delta/enorm(delta,weightmat)
             dnew <- sqdist(x)
             #dnew <- sqrt(sqdist(x))
             #dnew <- dnew/enorm(dnew,weightmat)
             #dnew <- dnew^2
#r <- 2*r
             rnew <- sum (weightmat * delta * mkPower (dnew, r))
             nnew <- sum (weightmat * mkPower (dnew,  2*r))
             anew <- rnew / nnew
             snew <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
             resen <- abs(mkPower(dnew,r)-delta)
             #resen <- abs(dnew-delta)
             #shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,...)
             shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scaleB=scaleB,normed=normed,...) 
             #shrinkres <- resen-cordweight*resen*shrinkb/(resen+shrinkb)
             shrinkres <- resen*(1-cordweight*(shrinkb/(resen+shrinkb)))
             #shrinkres <- resen
             diag(shrinkres) <- 0
             #TODO check for increasing residual
             ic <- sum(shrinkres^2)
             if(verbose>3) cat("copstress =",ic,"mdslossm =",sum(resen^2),"delta(cop/mds)=",ic-sum(resen^2),"mdslosss =",snew,"delta(mds/sma)=",sum(resen^2)-snew,"\n")
             #delta cops/mds should be positive if cordweight is too high,no?
             ic
           }
     if(verbose>1) cat("Starting Minimization with",optimmethod,":\"n")
     if(optimmethod=="Newuoa") {
         optimized <- minqa::newuoa(xold,function(par) shrinkcops(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat, cordweight=cordweight, q=q,minpts=minpts,epsilon=epsilon,rang=rang,scaleX=scaleX,scaleB=scaleB,enormX=enormX,normed=normed),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- optimized$feval
         ovalue <-optimized$fval
     }
     if(optimmethod=="Nelder-Mead") {
         optimized <- optim(xold,function(par) shrinkcops(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,
                       cordweight=cordweight, q=q,minpts=minpts,epsilon=epsilon,rang=rang,scaleX=scaleX,scaleB=scaleB,enormX=enormX,normed=normed),control=list(maxit=itmax,trace=verbose),...)
         xnew <- optimized$par
         itel <- optimized$counts[[1]]
         ovalue <-optimized$val 
     }
     if(enormX) xnew <- xnew/enorm(xnew)
     if(scaleX) xnew <- scale(xnew)
     #dnew <- as.matrix(dist (xnew)^2) #alternative
     dnew <- sqdist (xnew)
     rnew <- sum (weightmat * delta * mkPower (dnew, r))
     nnew <- sum (weightmat * mkPower (dnew,  2*r))
     anew <- rnew / nnew
     stress <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
     attr(xnew,"dimnames")[[1]] <- rownames(delta)
     attr(xnew,"dimnames")[[2]] <- paste("D",1:ndim,sep="")
     doutm <- (sqrt(dnew))^kappa  #fitted powered euclidean distance
     #doutm <- as.matrix(dist(xnew)^kappa)  #alternative 
     deltam <- delta
     deltaorigm <- deltaorig
     deltaoldm <- deltaold
     resmat <- deltam - doutm
     delta <- stats::as.dist(delta)
     deltaorig <- stats::as.dist(deltaorig)
     deltaold <- stats::as.dist(deltaold)
     doute <- doutm/enorm(doutm)
     doute <- stats::as.dist(doute)
     dout <- stats::as.dist(doutm)
     spp <- colMeans(resmat)
     weightmatm <-weightmat
     weightmat <- stats::as.dist(weightmatm)
     stressen <- sum(weightmatm*resmat^2)/2 #raw stress on the normalized proximities and normalized distances 
     if(verbose>1) cat("*** stress (both normalized - for COPS/STOPS):",stress,"; stress 1 (both normalized - default reported):",sqrt(stress),"; stress manual (for debug only):",stressen,"; from optimization: ",ovalue,"\n")   
    out <- list(delta=deltaold, obsdiss=delta, confdiss=dout, conf = xnew, pars=c(kappa,lambda,nu), niter = itel, stress=sqrt(stress), spp=spp, ndim=ndim, model="Copstress NEWUOA", call=match.call(), nobj = dim(xnew)[1], type = "copstress", gamma=NA, stress.m=stress, stress.en=stressen, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat)
    out$par <- theta
    out$loss <- "copstress"
    out$OC <- cordillera::cordillera(out$conf,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scaleC)
    out$copstress <- ovalue
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


#' Finding the shrinkage matrix for COPS-0
#'
#' @param x numeric matrix
#' @param q the norm to be used
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to ndim+1
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the minimum reachabilities to be considered. If missing it is found from the initial configuration.
#' @param scaleB should the X matrix be scaled before calculating the shrinkage weights; defaults to TRUE (which is sensible for relative shrinkage)
#' @param normed should the reachability differences be normed? (1=to sup Gamma 2 = to dmax)?
#' @param ... additional arguments to be passed to the OPTICS algorithm procedure
#' 
#' @importFrom dbscan optics
#' 
#'@export
shrinkB <- function(x,q=1,minpts=2,epsilon=10,rang=NULL,scaleB=TRUE,normed=2,...)
     {
       shift1 <- function (v) {
             vlen <- length(v)
             sv <- as.vector(numeric(vlen)) 
             sv[-1] <- v[1:(vlen - 1)]
             sv
        }
        if(scaleB) x <- scale(x)
        N <- dim(x)[1]
        optres<-dbscan::optics(x,minPts=minpts,eps=epsilon,...)
        optord <- optres$order
        optind <- 1:length(optres$order)
        indordered <- optind[optord]
        predec <- shift1(indordered)
        reachdist <- optres$reachdist[optres$order]
        reachdist[!is.finite(reachdist)] <- ifelse(is.null(rang),max(reachdist[is.finite(reachdist)]),min(max(rang),max(reachdist[is.finite(reachdist)])))
        reachdiffs <- c(NA,abs(diff(reachdist)))
        mats <- cbind(indordered,predec,reachdiffs)
        Bmat <- matrix(0,ncol=N,nrow=N)
        for (i in 2:N) {
            indo <- mats[i,]
            Bmat[indo[1],indo[2]] <- Bmat[indo[2],indo[1]] <- indo[3]/(2^(1/q))
        }
        if(normed==1) Bmat <- Bmat/(((rang[2]-rang[1])^q)*(ceiling((N-1)/minpts)+floor((N-1)/minpts)))
        if(normed==2) Bmat <- (Bmat*2^(1/q))/(rang[2]-rang[1])^q 
        return(Bmat)
     }




###################################
#Shrinkage cops older version
shrinkCoploss0 <- function (delta, kappa=1, lambda=1, nu=1, theta=c(kappa,lambda,nu),weightmat=1-diag(nrow(delta)),  ndim = 2, init=NULL,cordweight=1,q=2,minpts=ndim+1,epsilon=10,rang=NULL,optimmethod=c("Nelder-Mead","Newuoa"),verbose=0,scale=TRUE,accuracy = 1e-7, itmax = 100000,...)
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
    shrinkcops0 <- function(x,delta,r,ndim,weightmat,cordweight,q,minpts,epsilon,rang,...)
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
         optimized <- minqa::newuoa(xold,function(par) shrinkcops0(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,
                       cordweight=cordweight, q=q,minpts=minpts,epsilon=epsilon,rang=rang
                                   ),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- optimized$feval
         ovalue <-optimized$fval
     }
     if(optimmethod=="Nelder-Mead") {
         optimized <- optim(xold,function(par) shrinkcops0(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,
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
    

####Another COPS model
 
#' Fitting a COPS Model by penalizing residuals (COPS-C1).
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
#' @param cordweight weight to be used for the shrinkage; defaults to 1
#' @param q used in cordillera and shrink matrix and controls the effect of using the norm; defaults to 2 (least squares MDS)
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to ndim+1
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the minimum reachabilities to be considered. If missing it is found from the initial configuration by taking 1.5 times the maximal minimum reachability of the initial fit. If NULL it will be normed to each configuration's minimum and maximum distance, so an absolute value of goodness-of-clusteredness. Note that the latter is not necessarily desirable when comparing configurations for their relative clusteredness. See also \code{\link{cordillera}}
#' @param scaleX should X be scaled; defaults to TRUE
#' @param enormX should X be enormed; defaults to FALSE
#' @param scaleB should X be scaled for the shrink matrix; defaults to TRUE.
#' @param scaleC should X be scaled for the OPTICS Cordillera; defaults to TRUE. These parameter lets one tweak the way the shrinkage works and how its quantified; the defaults lead usually to a sensible result. It might be that some scale versions weill be depreciated in future versions. 
#' @param optimmethod What optimizer to use? Defaults to NEWUOA, Nelder-Mead is also supported.
#' @param verbose numeric value hat prints information on the fitting process; >2 is very verbose
#' @param accuracy numerical accuracy, defaults to 1e-8
#' @param itmax maximum number of iterations. Defaults to 100000
#' @param normed should cordillera and shrink matrix be normed; defaults to TRUE
#' @param ... additional arguments to be passed to the optimization procedure
#'
#' @return A list with the components
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
#' @examples
#' dis<-as.matrix(smacof::kinshipdelta)
#'
#' #Coploss with shrinkage to 0 
#' res1<-shrinkCoploss(dis,cordweight=1,minpts=2) 
#' res1
#' summary(res1)
#' plot(res1)  #super clustered
#'
#' @importFrom stats dist as.dist optim
#' @importFrom minqa newuoa
#' 
#' 
#' @keywords clustering multivariate
#' @export
shrinkCoploss1 <- function (delta, kappa=1, lambda=1, nu=1, theta=c(kappa,lambda,nu),weightmat=1-diag(nrow(delta)),  ndim = 2, init=NULL,cordweight=1,q=2,minpts=ndim+1,epsilon=10,rang=NULL,optimmethod=c("Nelder-Mead","Newuoa"),verbose=0,scaleX=TRUE,enormX=FALSE,scaleB=TRUE,scaleC=TRUE,accuracy = 1e-7, itmax = 100000,normed=2,...)
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
           if(isTRUE(scaleX)) init0 <- scale(init0)
           crp <- stops::cordillera(init0,q=q,minpts=minpts,epsilon=epsilon,scale=scaleC)$reachplot
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
    if(is.null(init)) xold <- stops::powerStressFast(delta,kappa=kappa,lambda=lambda,nu=nu,ndim=ndim)$conf
    if(enormX) xold <- xold/enorm(xold)
    if(scaleX) xold <- scale(xold)
    shrinkcops1 <- function(x,delta,r,ndim,weightmat,cordweight,q,minpts,epsilon,rang,scaleX,enormX,scaleB,normed,...)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=ndim)
             if(scaleX) x <- scale(x)
             if(enormX) x <- x/enorm(x)
             N <- dim(x)[1]
             delta <- delta/enorm(delta,weightmat)
             dnew <- sqdist(x)
             rnew <- sum (weightmat * delta * mkPower (dnew, r))
             nnew <- sum (weightmat * mkPower (dnew,  2*r))
             anew <- rnew / nnew
             snew <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
             resen <- abs(mkPower(dnew,r)-delta)
             #corrd <- stops::cordillera(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,...)
             shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scaleB=scaleB,normed=normed,...)
             corrd <- sum(shrinkb)
             #shrinks <- matrix(corrd,nrow=dim(shrinkb)[1],ncol=dim(shrinkb)[2])
             #diag(shrinks) <- 0
             shrinkres <- resen-cordweight*sqrt(corrd/(N^2-N))
             diag(shrinkres) <- 0
             #TODO check for increasing residual
             ic <- sum(shrinkres^2)
             if(verbose>3) cat("coploss =",ic,"mdslossm =",sum(resen^2),"delta(cop/mds)=",ic-sum(resen^2),"mdslosss =",snew,"delta(mds/sma)=",sum(resen^2)-snew,"\n")
             ic
           }
     if(verbose>1) cat("Starting Minimization with",optimmethod,":\"n")
     if(optimmethod=="Newuoa") {
         optimized <- minqa::newuoa(xold,function(par) shrinkcops1(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat, cordweight=cordweight, q=q,minpts=minpts,epsilon=epsilon,rang=rang,scaleX=scaleX,scaleB=scaleB,enormX=enormX,normed=normed),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose),...)
         xnew <- matrix(optimized$par,ncol=ndim)
         itel <- optimized$feval
         ovalue <-optimized$fval
     }
     if(optimmethod=="Nelder-Mead") {
         optimized <- optim(xold,function(par) shrinkcops1(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,
                       cordweight=cordweight, q=q,minpts=minpts,epsilon=epsilon,rang=rang,scaleX=scaleX,scaleB=scaleB,enormX=enormX,normed=normed),control=list(maxit=itmax,trace=verbose),...)
         xnew <- optimized$par
         itel <- optimized$counts[[1]]
         ovalue <-optimized$val 
     }
     if(enormX) xnew <- xnew/enorm(xnew)
     if(scaleX) xnew <- scale(xnew)
     #dnew <- as.matrix(dist (xnew)^2) #alternative
     dnew <- sqdist (xnew)
     rnew <- sum (weightmat * delta * mkPower (dnew, r))
     nnew <- sum (weightmat * mkPower (dnew,  2*r))
     anew <- rnew / nnew
     stress <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
     attr(xnew,"dimnames")[[1]] <- rownames(delta)
     attr(xnew,"dimnames")[[2]] <- paste("D",1:ndim,sep="")
     doutm <- (sqrt(dnew))^kappa  #fitted powered euclidean distance
     #doutm <- as.matrix(dist(xnew)^kappa)  #alternative 
     deltam <- delta
     deltaorigm <- deltaorig
     deltaoldm <- deltaold
     resmat <- deltam - doutm
     delta <- stats::as.dist(delta)
     deltaorig <- stats::as.dist(deltaorig)
     deltaold <- stats::as.dist(deltaold)
     doute <- doutm/enorm(doutm)
     doute <- stats::as.dist(doute)
     dout <- stats::as.dist(doutm)
     spp <- colMeans(resmat)
     weightmatm <-weightmat
     weightmat <- stats::as.dist(weightmatm)
     stressen <- sum(weightmatm*resmat^2)/2 #raw stress on the normalized proximities and normalized distances 
     if(verbose>1) cat("*** stress (both normalized - for COPS/STOPS):",stress,"; stress 1 (both normalized - default reported):",sqrt(stress),"; stress manual (for debug only):",stressen,"; from optimization: ",ovalue,"\n")   
    out <- list(delta=deltaold, obsdiss=delta, confdiss=dout, conf = xnew, pars=c(kappa,lambda,nu), niter = itel, stress=sqrt(stress), spp=spp, ndim=ndim, model="Coploss NEWUOA", call=match.call(), nobj = dim(xnew)[1], type = "coploss", gamma=NA, stress.m=stress, stress.en=stressen, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat)
    out$par <- theta
    out$loss <- "coploss"
    out$OC <- stops::cordillera(out$conf,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scaleC)
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


shrinkcops1 <- function(x,delta,r,ndim,weightmat,cordweight,q,minpts,epsilon,rang,scaleX=TRUE,enormX=TRUE,scaleB=TRUE,normed=TRUE,...)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=ndim)
             #if(scaleX) x <- scale(x)
             #if(enormX) x <- x/enorm(x)
             N <- dim(x)[1]
             delta <- delta/enorm(delta,weightmat)
             dnew <- sqdist(x)
             rnew <- sum (weightmat * delta * mkPower (dnew, r))
             nnew <- sum (weightmat * mkPower (dnew,  2*r))
             anew <- rnew / nnew
             snew <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
             resen <- abs(mkPower(dnew,r)-delta)
             #corrd <- stops::cordillera(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,...)
             #corrd <- stops::cordillera(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale)
             shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scaleB=scaleB,normed=normed,...)
            # shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scaleB=scaleB,normed=normed)
             corrd <- sum(shrinkb)
             #shrinks <- matrix(corrd,nrow=dim(shrinkb)[1],ncol=dim(shrinkb)[2])
             #diag(shrinks) <- 0
             shrinkres <- resen-cordweight*sqrt(corrd/(N^2-N))
             diag(shrinkres) <- 0
             #TODO check for increasing residual
             ic <- sum(shrinkres^2)
             if(verbose>1) cat("coploss =",ic,"mdslossm =",sum(resen^2)/2,"delta(cop/mds)=",ic-sum(resen^2),"mdslosss =",snew,"delta(mds/sma)=",sum(resen^2)-snew,"\n")
             ic
           }

#checkl with how it is in the cops.R code

cordweight=1
q <- 1
optimized <- minqa::newuoa(xold,function(par) shrinkcops1(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat, cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose))
xnew1<- matrix(optimized$par,ncol=ndim)

optimizedn <- optim(xold,function(par) shrinkcops1(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat, cordweight=cordweight, q=q,minpts=minpts,epsilon=epsilon,rang=rang), method="Nelder-Mead")
xnewn <- matrix(optimizedn$par,ncol=ndim)

library(crs)
#q <- 0
optimizeds <- crs::snomadr(eval.f=function(par) shrinkcops1(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,cordweight=cordweight,q=1,minpts=minpts,epsilon=epsilon,rang=rang),n=length(xold),x0=xnew1,bbin=rep(0,30),ub=rep(5,30),lb=rep(-5,30,bbout=0))
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

shrinkcops0(x=x,delta=delta,r=r,ndim=2,weightmat=weightmat,cordweight=cordweight,rang=rang,q=q,minpts=minpts,epsilon=epsilon)

shrinkcops1(x=x,delta=delta,r=r,ndim=2,weightmat=weightmat,cordweight=cordweight,rang=rang,q=q,minpts=minpts,epsilon=epsilon)

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




####Str Embed
library(stops)
scale <- TRUE
delta <- kinshipdelta
weightmat <- 1-diag(15)
xsma <- powerStressMin(delta)
#xsma <- powerStressFast(delta)
x <- xsma$conf
ndim <- 2
cordweight <- 1
r <- 0.5
q <- 1
minpts <- 2
epsilon <- 10
rang <- c(0,1.6)
#plot(xsma)
xold <- xsma$conf
itmax <- 100000
accuracy <- 1e-12
verbose=2
nameso <- rownames(xold)

shrinkemb <- function(x,delta,r,ndim,type=c("additive","sqadditive","multiplicative"),weightmat,cordweight,q=2,minpts,epsilon,rang,...)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=ndim)
             x <- scale(x)
             delta <- delta/enorm(delta,weightmat)
             x <- x/enorm(x)
             dnew <- sqdist(x)
             rnew <- sum (weightmat * delta * mkPower (dnew, r))
             nnew <- sum (weightmat * mkPower (dnew,  2*r))
             anew <- rnew / nnew
             resen <- abs(mkPower(dnew,r)-delta)
             shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,...)
             #shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang) 
             #shrinkres <- resen-cordweight*resen*shrinkb/(resen+shrinkb)
             #cat(c(min(resen[lower.tri(resen)]),max(shrinkb[lower.tri(shrinkb)])),"\n")
             if(type=="additive") shrinkres <- resen-cordweight*sqrt(shrinkb)
             if(type=="sqadditive") shrinkres <- resen^2-cordweight*shrinkb
             if(type=="multiplicative") shrinkres <- resen*(1-cordweight*(shrinkb/(abs(resen)+shrinkb)))
             #shrinkres <- resen
             diag(shrinkres) <- 0
             #so there is the problem of shrinkres being below 0; how to deal with that? Either check for that and if thats the case, then the unshrunk sum of residuals is used; or fix the negative residual to zero.
             if(any(shrinkres<0)) {
                                   tshrinkres <- shrinkres[lower.tri(shrinkres)]
                                   t1 <- tshrinkres[which(tshrinkres<0)]
                                   t2 <- which(shrinkres<0,arr.ind=TRUE)
                                   if(verbose>4) cat("Residuals<0: Value", t1,"Indices",t2,"\n")
                                   #return(sum(resen^2))
                                   shrinkres[t2] <- 0
                                   if(verbose>10) cat("Residuals still problematic?",any(shrinkres<0),"\n")
             }
             
             snew <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
             if(type=="sqadditive") shrinkres <- sqrt(shrinkres)
             #Tis weird why does it become best if I have very low weight?
             ic <- sum(shrinkres^2)
             if(verbose>2) cat("coploss =",ic,"mdslossm =",sum(resen^2),"delta(cop/mds)=",ic-sum(resen^2),"mdslosss =",snew,"delta(mds/sma)=",sum(resen^2)-snew,"\n")
             #cop must be smaller than mds
             ic
           }


cordweight <- 0
q <- 1
optimized <- minqa::newuoa(xold,function(par) shrinkemb(par,delta=delta,r=r,ndim=ndim,type="additive",weightmat=weightmat, cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose))
xold<- matrix(optimized$par,ncol=ndim)
cordweight <- 0.1
optimized <- minqa::newuoa(xold,function(par) shrinkemb(par,delta=delta,r=r,ndim=ndim,type="additive",weightmat=weightmat, cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose))
xnewa<- matrix(optimized$par,ncol=ndim)
par(mfrow=c(1,2))
plot(xold,asp=1)
text(xold,label=nameso,pos=3)
plot(xnewa,asp=1)
text(xnewa,label=nameso,pos=3)


cordweight <- 0.05
optimized <- minqa::newuoa(xold,function(par) shrinkemb(par,delta=delta,r=r,ndim=ndim,type="sqadditive",weightmat=weightmat, cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose))
xnewas<- matrix(optimized$par,ncol=ndim)
par(mfrow=c(1,2))
plot(xold,asp=1)
text(xold,label=nameso,pos=3)
plot(xnewas,asp=1)
text(xnewas,label=nameso,pos=3)

cordweight <- 1
optimized <- minqa::newuoa(xold,function(par) shrinkemb(par,delta=delta,r=r,ndim=ndim,type="multiplicative",weightmat=weightmat, cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose))
xnewm<- matrix(optimized$par,ncol=ndim)
par(mfrow=c(1,2))
plot(xold,asp=1)
text(xold,label=nameso,pos=3)
plot(xnewm,asp=1)
text(xnewm,label=nameso,pos=3)


par(mfrow=c(2,2))
plot(xold,asp=1,main="orig")
text(xold,label=nameso,pos=3)
plot(xnewa,asp=1,main="add")
text(xnewa,label=nameso,pos=3)
plot(xnewas,asp=1,main="addsq")
text(xnewas,label=nameso,pos=3)
plot(xnewm,asp=1,main="mult")
text(xnewm,label=nameso,pos=3)





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



optimizedn <- optim(xold,function(par) shrinkemb(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat, cordweight=cordweight, q=q,minpts=minpts,epsilon=epsilon,rang=rang), method="Nelder-Mead")
xnewn <- matrix(optimizedn$par,ncol=ndim)



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

#Cali
soviagg <- read.csv("~/svn/egmds/paper/CaliMedAgg.csv")       
data(CAClimateIndicatorsCountyMedian)
sovisel <- CAClimateIndicatorsCountyMedian[,-1]
sovisel <- soviagg[,c(23:34,40:77)]       
sovisel <- apply(sovisel,2,function(x) (x-min(x))/(max(x)-min(x)))
dis <- dist(sovisel)
delta <- as.matrix(dis)
rang <- c(0,1.2) #fixed range 
minpts <- 3
epsilon <- 10
q <- 1
weightmat <- 1-diag(58)

#kinship
delta <- kinshipdelta
weightmat <- 1-diag(15)
q <- 1
minpts <- 2
epsilon <- 10
rang <- c(0,1.6)

#Banking
data(BankingCrisesDistances)
delta <- BankingCrisesDistances[,-70]
weightmat <- 1-diag(69)
q <- 1
minpts <- 5
epsilon <- 10
rang <- c(0,0.7)
minpts <- 2


#for both
ndim <- 2
r <- 0.5
xsma <- powerStressFast(delta)
x <- xsma$conf
itmax <- 100000
accuracy <- 1e-12
verbose=2
xold <- x
normed <- 2

shrinkB <- function(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=TRUE,normed=TRUE,...)
     {
       shift1 <- function (v) {
             vlen <- length(v)
             sv <- as.vector(numeric(vlen), mode = storage.mode(v)) 
             sv[-1] <- v[1:(vlen - 1)]
             sv
        }
        if(scale) x <- scale(x)
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

x <- scale(xold)
shrinkcops <- function(x,delta,r=0.5,ndim,weightmat,cordweight,q=2,minpts,epsilon,rang,scaleX=TRUE,enormX=FALSE,scaleB=TRUE,normed=TRUE,...)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=ndim)
             #try these variants again with kinship and cali:
             #it all about how to normalize so that the shrinkage will play its part
             #take care of r and so if sqrt(sqdist()) use r*2
             #delta enormed, x scaled + enormed; looks good! -> looks best? 
             #delta enormed, x scaled, dnew enormed; looks ok like #2 but a bit better
             #delta enormed, x enormed, dnew normal; looks ok with clusters for kinship but wrong clusters; closest snew and mdsloss
             if(scaleX) x <- scale(x)
             delta <- delta/enorm(delta,weightmat)
             if(enormX) x <- x/enorm(x)
             dnew <- sqdist(x)
             #dnew <- sqrt(sqdist(x))
             #dnew <- dnew/enorm(dnew,weightmat)
             #dnew <- dnew^2
#r <- 2*r
             rnew <- sum (weightmat * delta * mkPower (dnew, r))
             nnew <- sum (weightmat * mkPower (dnew,  2*r))
             anew <- rnew / nnew
             snew <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
             resen <- abs(mkPower(dnew,r)-delta)
             #resen <- abs(dnew-delta)
             #shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=scale,...)
             shrinkb <- shrinkB(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scaleB=scaleB,normed=normed) 
             #shrinkres <- resen-cordweight*resen*shrinkb/(resen+shrinkb)
             shrinkres <- resen*(1-cordweight*(shrinkb/(resen+shrinkb)))
             #shrinkres <- resen
             diag(shrinkres) <- 0
             #TODO check for increasing residual
             ic <- sum(shrinkres^2)
             if(verbose>1) cat("coploss =",ic,"mdslossm =",sum(resen^2),"delta(cop/mds)=",ic-sum(resen^2),"mdslosss =",snew,"delta(mds/sma)=",sum(resen^2)-snew,"\n")
             #delta cops/mds should be positive if cordweight is too high,no?
             ic
           }

cordweight <- 1
optimized <- minqa::newuoa(x,function(par) shrinkcops(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,normed=normed),control=list(maxfun=itmax,rhoend=accuracy,iprint=verbose))
xnew <- matrix(optimized$par,ncol=ndim)
verbose <- 2
xnew2 <- xnew



library(crs)
#q <- 0
optimizeds <- crs::snomadr(eval.f=function(par) shrinkcops(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat,cordweight=cordweight,q=1,minpts=minpts,epsilon=epsilon,rang=rang),n=length(xold),x0=xold, bbin=rep(0,116),ub=rep(3,116),lb=rep(-3,116),bbout=0)
xnew <- matrix(optimizeds$solution,ncol=ndim)


xnew <- xnew/enorm(xnew)

par(mfrow=c(1,2))
plot(x,asp=1,pch=20)
text(scale(x),label=rownames(x),pos=3)
plot(xnew,asp=1,pch=20)
text(scale(xnew),label=rownames(x),pos=3)

plot(scale(x),asp=1,pch=20)
#text(x,label=rownames(x),pos=3)
plot(scale(xnew),asp=1,pch=20)

cx <- cordillera(x,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=TRUE)
cx
plot(cx)

cn <- cordillera(xnew,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=TRUE)
cn
plot(cn)


xa <- conf_adjust(x,xnew)$ref.conf
xnewa <- conf_adjust(x,xnew)$other.conf

par(mfrow=c(1,2))
plot(xa,asp=1,pch=20)
text(xa,label=rownames(x),pos=3)
plot(xnewa,asp=1,pch=20)
text(xnewa,label=rownames(x),pos=3)

cxa <- cordillera(xa,q=q,minpts=minpts,epsilon=epsilon,rang=rang)
cxa
plot(cxa)

cna <- cordillera(xnewa,q=q,minpts=minpts,epsilon=epsilon,rang=rang)
cna
plot(cna)

xaa <- scale(xa)
xnewaa <- scale(xnewa)

par(mfrow=c(1,2))
plot(xaa,asp=1,pch=20)
points(xnewaa,col="green")

points(xnew,col="green")

text(xaa,label=rownames(x),pos=3)
plot(xnewaa,asp=1)
text(xnewaa,label=rownames(xold),pos=3)


points(xnewa,col="red")

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


cordweight=1
m01 <- shrinkCoploss(delta,init=xold,minpts=minpts,epsilon=epsilon,q=q,ndim=ndim,cordweight=cordweight,weightmat=weightmat,rang=rang,verbose=2)
m01
plot(m01)

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


rs <- delta/enorm(delta)

r2 <- as.matrix(dist(xnew))
res <- abs(rs-r2)
diag(res) <- 1


shrinks <- shrinkB(xnew,q=q,minpts=minpts,epsilon=epsilon,rang=rang)

ft <- (res*shrinks)/(res+shrinks)

res-res*cw*ft>0

res>res*cw*ft

res/ft>res*cw

1/ft>cw


#' Double centering of a matrix
#'
#' @param x numeric matrix
#' @return the double centered matrix
doubleCenter <- function(x) {
        n <- dim(x)[1]
        m <- dim(x)[2]
        s <- sum(x)/(n*m)
        xr <- rowSums(x)/m
        xc <- colSums(x)/n
        return((x-outer(xr,xc,"+"))+s)
    }

# #' Torgerson scaling
# #'
# #' @param delta symmetric, numeric matrix of distances
# #' @param p target space dimensions
# #' @return a n x p matrix (the configuration)
# #' @export
# #' @examples
# #' dis<-as.matrix(smacof::kinshipdelta)
# #' res<-torgerson(dis)
#torgerson <- function(delta, p = 2) {
#    z <- eigen(-doubleCenter((as.matrix (delta) ^ 2)/2))
#    v <- pmax(z$values,0)
#    return(z$vectors[,1:p]%*%diag(sqrt(v[1:p])))
#}

#' Explicit Normalization
#' Normalizes distances
#' @param x numeric matrix 
#' @param w weight
#' @return a constant 
enorm <- function (x, w=1) {
    return (sqrt (sum (w * (x ^ 2))))
}

#' Squared distances
#'
#' @param x numeric matrix
#' @return squared distance matrix
sqdist <- function (x) {
    s <- tcrossprod (x)
    v <- diag (s)
    return (outer (v, v, "+") - 2 * s)
}

#' Squared p-distances
#'
#' @param x numeric matrix
#' @param p p>0 the Minkoswki distance
#' @return squared Minkowski distance matrix
pdist <- function (x,p) {
    s <- tcrossprod (x)
    v <- diag (s)
    return (outer (v, v, "+") - 2 * s)
}

#' Auxfunction1
#' 
#' only used internally
#' @param x matrix
mkBmat <- function (x) {
    d <- rowSums (x)
    x <- -x
    diag (x) <- d
    return (x)
}


#' Take matrix to a power 
#'
#' @param x matrix
#' @param r numeric (power)
#' @return a matrix
mkPower<-function(x,r) {
    n<-nrow(x)
    tmp <- abs((x+diag(n))^r)-diag(n)
    return(tmp)
}


#' Secular Equation 
#'
#' @param a matrix
#' @param b matrix
#'
#' @importFrom stats uniroot
secularEq<-function(a,b) {
    n<-dim(a)[1]
    eig<-eigen(a)
    eva<-eig$values
    eve<-eig$vectors
    beta<-drop(crossprod(eve, b))
    f<-function(mu) {
        return(sum((beta/(eva+mu))^2)-1)
    }
    lmn<-eva [n]
    uup<-sqrt(sum(b^2))-lmn
    ulw<-abs(beta [n])-lmn
    rot<-stats::uniroot(f,lower= ulw,upper= uup)$root
    cve<-beta/(eva+rot)
    return(drop(eve%*%cve))
}    


#' Power stress minimization by NEWUOA (nloptr)
#'
#' An implementation to minimize power stress by a derivative-free trust region optimization algorithm (NEWUOA). Much faster than majorizing as used in powerStressMin but perhaps less accurate. 
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param kappa power of the transformation of the fitted distances; defaults to 1
#' @param lambda the power of the transformation of the proximities; defaults to 1
#' @param nu the power of the transformation for weightmat; defaults to 1 
#' @param weightmat a matrix of finite weights
#' @param init starting configuration
#' @param ndim dimension of the configuration; defaults to 2
#' @param acc  The smallest value of the trust region radius that is allowed. If not defined, then 1e-6 will be used.
#' @param itmax maximum number of iterations. Default is 10000.
#' @param verbose should iteration output be printed; if > 1 then yes
#'
#' @return a smacofP object (inheriting form smacofB, see \code{\link{smacofSym}}). It is a list with the components
#' \itemize{
#' \item delta: Observed dissimilarities, not normalized
#' \item obsdiss: Observed dissimilarities, normalized 
#' \item confdist: Configuration dissimilarities, NOT normalized 
#' \item conf: Matrix of fitted configuration, NOT normalized
#' \item stress: Default stress (stress 1, square root of the explicitly normalized stress on the normalized, transformed dissimilarities)  
#' \item spp: Stress per point (based on stress.en) 
#' \item ndim: Number of dimensions
#' \item model: Name of smacof model
#' \item niter: Number of iterations
#' \item nobj: Number of objects
#' \item type: Type of MDS model
#' }
#' and some additional components
#' \itemize{
#' \item gamma: Empty
#' \item stress.m: default stress for the COPS and STOP. Defaults to the explicitly normalized stress on the normalized, transformed dissimilarities
#' \item stress.en: explicitly stress on the normalized, transformed dissimilarities and normalized transformed distances
#' \item deltaorig: observed, untransformed dissimilarities
#' \item weightmat: weighting matrix 
#'}
#'
#' @importFrom stats dist as.dist
#' @importFrom minqa newuoa
#' 
#' @seealso \code{\link{smacofSym}}
#' 
#' @examples
#' dis<-smacof::kinshipdelta
#' res<-powerStressFast(as.matrix(dis),kappa=2,lambda=1.5)
#' res
#' summary(res)
#' plot(res)
#' 
#' @export
powerStressFast <- function (delta, kappa=1, lambda=1, nu=1, weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc = 1e-6, itmax = 10000, verbose = FALSE)
{
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    if(verbose>0) cat("Minimizing powerstress by NEWUOA with kappa=",kappa,"lambda=",lambda,"nu=",nu,"\n")
    type <- "ratio"
    r <- kappa/2
    p <- ndim
    deltaorig <- delta
    delta <- delta^lambda
    weightmato <- weightmat
    weightmat <- weightmat^nu
    weightmat[!is.finite(weightmat)] <- 1 #new
    deltaold <- delta
    delta <- delta / enorm (delta, weightmat) #sum=1
    xold <- init
    if(is.null(init)) xold <- smacof::torgerson (delta, p = ndim)
    xold <- xold/enorm(xold) 
    stressf <- function(x,delta,r,ndim,weightmat)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=ndim)
             delta <- delta/enorm(delta,weightmat)
             x <- x/enorm(x)
             #adapted from powerStressMin 
             dnew <- sqdist (x)
             rnew <- sum (weightmat * delta * mkPower (dnew, r))
             nnew <- sum (weightmat * mkPower (dnew,  2*r))
             anew <- rnew / nnew
             snew <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
             snew
           }
     suppressWarnings(optimized <- minqa::newuoa(xold,function(par) stressf(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat),control=list(maxfun=itmax,rhoend=acc,iprint=verbose-2)))
      xnew <- matrix(optimized$par,ncol=ndim)
      xnew <- xnew/enorm(xnew)
             #adapted from powerStressMin 
      dnew <- sqdist (xnew)
      rnew <- sum (weightmat * delta * mkPower (dnew, r))
      nnew <- sum (weightmat * mkPower (dnew,  2*r))
      anew <- rnew / nnew
      stress <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
     #xnew <- optimized$par
      attr(xnew,"dimnames")[[1]] <- rownames(delta)
      itel <- optimized$feval
      attr(xnew,"dimnames")[[2]] <- paste("D",1:ndim,sep="")
      doutm <- (2*sqrt(sqdist(xnew)))^kappa  #fitted powered euclidean distance but times two
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
      stressen <- sum(weightmat*(doute-delta)^2) #raw stress on the normalized proximities and normalized distances 
      if(verbose>1) cat("*** stress (both normalized - for COPS/STOPS):",stress,"; stress 1 (both normalized - default reported):",sqrt(stress),"; stress manual (for debug only):",stressen,"; from optimization: ",optimized$fval,"\n")   
    out <- list(delta=deltaold, obsdiss=delta, confdist=dout, conf = xnew, pars=c(kappa,lambda,nu), niter = itel, stress=stress, spp=spp, ndim=p, model="Power-Stress NEWUOA", call=match.call(), nobj = dim(xnew)[1], type = type, gamma = NA, stress.m=sqrt(stress), stress.en=stressen, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat)
    class(out) <- c("smacofP","smacofB","smacof")
    out
}

#' Power Stress SMACOF
#'
#' An implementation to minimize power stress by minimization-majorization. Usually more accurate but slower than powerStressFast. Uses a repeat loop.
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param kappa power of the transformation of the fitted distances; defaults to 1
#' @param lambda the power of the transformation of the proximities; defaults to 1
#' @param nu the power of the transformation for weightmat; defaults to 1 
#' @param weightmat a matrix of finite weights
#' @param init starting configuration
#' @param ndim dimension of the configuration; defaults to 2
#' @param acc numeric accuracy of the iteration. Default is 1e-6.
#' @param itmax maximum number of iterations. Default is 10000.
#' @param verbose should iteration output be printed; if > 1 then yes
#'
#' @return a smacofP object (inheriting form smacofB, see \code{\link{smacofSym}}). It is a list with the components
#' \itemize{
#' \item delta: Observed dissimilarities, not normalized
#' \item obsdiss: Observed dissimilarities, normalized 
#' \item confdist: Configuration dissimilarities, NOT normalized 
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
#' @section Note:
#' The functionality related to power stress and the smacofP class is also available in the stops package (\code{\link[stops]{powerStressMin}}). Expect masking when both are loaded.   
#'
#' @importFrom stats dist as.dist
#' 
#' @seealso \code{\link{smacofSym}}
#' 
#' @examples
#' dis<-smacof::kinshipdelta
#' res<-powerStressMin(as.matrix(dis),kappa=2,lambda=1.5,itmax=1000)
#' res
#' summary(res)
#' plot(res)
#' 
#' @export
powerStressMin <- function (delta, kappa=1, lambda=1, nu=1, weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE) {
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    if(verbose>0) cat("Minimizing powerStress with kappa=",kappa,"lambda=",lambda,"nu=",nu,"\n")
    type <- "ratio"
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
    if(is.null(init)) xold <- smacof::torgerson (delta, p = p)
    xold <- xold / enorm (xold)
    n <- nrow (xold)
    nn <- diag (n)
    dold <- sqdist (xold)
    rold <- sum (weightmat * delta * mkPower (dold, r))
    nold <- sum (weightmat * mkPower (dold, 2 * r))
    aold <- rold / nold
    sold <- 1 - 2 * aold * rold + (aold ^ 2) * nold
    repeat {
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
#      if(is.na(snew)) #to avoid numerical issues if there are zeros somewhere 
#         {
            #  break ()
            #  xnew <- xold
            #  dnew <- dold
            #  sold <- snew
            #  aold <- anew
         # }
      if ((itel == itmax) || ((sold - snew) < acc))
        break ()
      itel <- itel + 1
      xold <- xnew
      dold <- dnew
      sold <- snew
      aold <- anew
     }
     attr(xnew,"dimnames")[[2]] <- paste("D",1:p,sep="")
     xnew <- xnew/enorm(xnew)
     doutm <- mkPower(sqdist(xnew),r)
     deltam <- delta
     delta <- stats::as.dist(delta)
     deltaorig <- stats::as.dist(deltaorig)
     deltaold <- stats::as.dist(deltaold)
     #doute <- doutm/enorm(doutm) #this is an issue here!
     #doute <- stats::as.dist(doute)
     dout <- stats::as.dist(doutm)
     weightmatm <-weightmat
     #resmat <- weightmatm*as.matrix((delta - doute)^2) #BUG
     resmat <- weightmatm*as.matrix((deltam - doutm)^2) #BUG
     spp <- colMeans(resmat) #BUG
     weightmat <- stats::as.dist(weightmatm)
     #stressen <- sum(weightmat*(doute-delta)^2)
     if(verbose>1) cat("*** Stress:",snew, "; Stress 1 (default reported):",sqrt(snew), "\n")  
    out <- list(delta=deltaold, obsdiss=delta, confdist=dout, conf = xnew, parameters=c(kappa,lambda,nu), pars=c(kappa,lambda,nu), theta=c(kappa,lambda,nu), niter = itel, spp=spp, ndim=p, model="Power-Stress SMACOF", call=match.call(), nobj = dim(xnew)[1], type = type, stress=sqrt(snew), stress.m=snew, stress.en=stressen, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat, alpha = anew, sigma = snew)
    class(out) <- c("smacofP","smacofB","smacof")
    out
}

#' @rdname powerStressMin
#' @export
powerstressMin <- powerStressMin

#' @rdname powerStressMin
#' @export
postmds <- powerStressMin

#' @rdname powerStressMin
#' @export
pstressMin <- powerStressMin

#' @rdname powerStressMin
#' @export
pStressMin <- powerStressMin

#' @rdname powerStressMin
#' @export
pstressmds <- powerStressMin


#' Power Stress SMACOF
#'
#' An implementation to minimize power stress by minimization-majorization. Usually more accurate but slower than powerStressFast. Uses a repeat loop.
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param kappa power of the transformation of the fitted distances; defaults to 1
#' @param lambda the power of the transformation of the proximities; defaults to 1
#' @param nu the power of the transformation for weightmat; defaults to 1 
#' @param weightmat a matrix of finite weights
#' @param init starting configuration
#' @param ndim dimension of the configuration; defaults to 2
#' @param acc numeric accuracy of the iteration. Default is 1e-6.
#' @param itmax maximum number of iterations. Default is 10000.
#' @param verbose should iteration output be printed; if > 1 then yes
#'
#' @return a smacofP object (inheriting form smacofB, see \code{\link{smacofSym}}). It is a list with the components
#' \itemize{
#' \item delta: Observed input dissimilarities, not normalized
#' \item tdelta: Explicitly transformed input dissimilarities, not normalized
#' \item dhats: Explicitly transformed dissimilarities, optimally scaled and normalized 
#' \item confdist: Configuration dissimilarities
#' \item conf: Matrix of fitted configuration
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
#' \item weightmat: weighting matrix 
#'}
#'
#' @section Note:
#' The functionality related to power stress and the smacofP class is also available in the stops package (\code{\link[stops]{powerStressMin}}). Expect masking when both are loaded.   
#'
#' @importFrom stats dist as.dist
#' 
#' @seealso \code{\link{smacofSym}}
#' 
#' @examples
#' dis<-smacof::kinshipdelta
#' res<-powerStressMin(as.matrix(dis),kappa=2,lambda=1.5,itmax=1000)
#' res
#' summary(res)
#' plot(res)
#' 
#' @export
powerStressMin2 <- function (delta, kappa=1, lambda=1, nu=1, weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE) {
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    if(verbose>0) cat("Minimizing powerStress with kappa=",kappa,"lambda=",lambda,"nu=",nu,"\n")
    type <- "ratio"
    r <- kappa/2
    p <- ndim
    deltaorig <- delta  #untransformed non-normalized input dissimilarities
    delta <- delta^lambda 
    #weightmato <- weightmat
    weightmat <- weightmat^nu
    weightmat[!is.finite(weightmat)] <- 1
    deltaold <- delta #explicitly transformed non-normalized input dissimilarities
    delta <- delta / enorm (delta, weightmat)
    itel <- 1
    xold <- init
    if(is.null(init)) xold <- smacof::torgerson (delta, p = p)
    xold <- xold / enorm (xold)
    n <- nrow (xold)
    nn <- diag (n)
    dold <- sqdist (xold)
    rold <- sum (weightmat * delta * mkPower (dold, r))
    nold <- sum (weightmat * mkPower (dold, 2 * r))
    aold <- rold / nold
    sold <- 1 - 2 * aold * rold + (aold ^ 2) * nold
    repeat {
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
#      if(is.na(snew)) #to avoid numerical issues if there are zeros somewhere 
#         {
            #  break ()
            #  xnew <- xold
            #  dnew <- dold
            #  sold <- snew
            #  aold <- anew
         # }
      if ((itel == itmax) || ((sold - snew) < acc))
        break ()
      itel <- itel + 1
      xold <- xnew
      dold <- dnew
      sold <- snew
      aold <- anew
     }
     attr(xnew,"dimnames")[[2]] <- paste("D",1:p,sep="")
     xnew <- xnew/enorm(xnew)
     doutm <- mkPower(sqdist(xnew),r)
     deltam <- delta
     delta <- stats::as.dist(delta)
     deltaorig <- stats::as.dist(deltaorig)
     deltaold <- stats::as.dist(deltaold)
     #doute <- doutm/enorm(doutm) #this is an issue here!
     #doute <- stats::as.dist(doute)
     dout <- stats::as.dist(doutm)
     weightmatm <-weightmat
     #resmat <- weightmatm*as.matrix((delta - doute)^2) #BUG
     resmat <- weightmatm*as.matrix((deltam - doutm)^2) #BUG
     spp <- colMeans(resmat) #BUG
     weightmat <- stats::as.dist(weightmatm)
     #stressen <- sum(weightmat*(doute-delta)^2)
     if(verbose>1) cat("*** Stress:",snew, "; Stress 1 (default reported):",sqrt(snew), "\n")  
    out <- list(delta=deltaorig, tdelta=deltaold, dhats=delta, confdist=dout, conf = xnew, parameters=c(kappa=kappa,lambda=lambda,nu=nu), pars=c(kappa=kappa,lambda=lambda,nu=nu), theta=c(kappa=kappa,lambda=lambda,nu=nu), niter = itel, spp=spp, ndim=p, model="Power-Stress SMACOF", call=match.call(), nobj = dim(xnew)[1], type = type, stress=sqrt(snew), stress.m=snew, stress.en=stressen,resmat=resmat,weightmat=weightmat, alpha = anew, sigma = snew)
    class(out) <- c("smacofP","smacofB","smacof")
    out
}

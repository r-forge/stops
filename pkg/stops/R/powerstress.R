#' Power stress
#'
#' A rudimentary implemenation to minimize power stress
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param kappa power of the transformation of the fitted distances; defaults to 1
#' @param lambda the power of the transformation of the proximities; defaults to 1
#' @param weightmat a matrix of finite weights
#' @param init starting configuration
#' @param ndim dimension of the configuration; defaults to 2
#' @param eps numeric accuracy of the iteration
#' @param itmax maximum number of iterations
#' @param verbose should iteration output be printed; if > 1 then yes
#' @param stresstype the stress type to be reported by default; defaults to stress-1 for smacof compatibility
#' @param lambdamax the maximum power of the transformation of the proximities if there are more than one lambda - an upper bound idea; defaults to lambda
#' 
#' @return a smacofB object (see \code{\link{smacofSym}}). It is a list with the components
#' \itemize{
#' \item delta: Observed dissimilarities, not normalized
#' \item obsdiss: Observed dissimilarities, normalized
#' \item confdiss: Configuration dissimilarities
#' \item conf: Matrix of fitted configurations
#' \item stress: Default stresstype 
#' \item spp: Stress per point
#' \item ndim: Number of dimensions
#' \item model: Name of smacof model
#' \item niter: Number of iterations
#' \item nobj: Number of objects
#' \item type: Type of MDS model
#' }
#' and some additional components
#' \itemize{
#' \item gamma: The majorizing function at convergence
#' \item stress.r: raw stress on the observed, transformed dissimilarities
#' \item stress.m: explicitly normalized stress on the observed, transformed dissimilarities
#' \item stress.n: explicitly normalized stress on the observed, transformed dissimilarities
#' \item stress.1: implicitly normalized stress on the observed, transformed dissimilarities
#' \item stress.b: explicitly and implicitly normalized stress on the observed, transformed dissimilarities
#' \item stress.e: raw stress on the normalized, transformed dissimilarities
#' \item stress.e1: implicitly normalized stress on the normalized, transformed dissimilarities
#' \item stress.be: explicitly and implicitly normalized stress on the normalized, transformed dissimilarities
#' \item stress.co: correlation of dissimilarities and fitted distances
#' \item deltaorig: observed, untransformed dissimilarities
#'}
#'
#' @seealso \code{\link{smacofSym}}
#' 
#' @examples
#' library(smacof)
#' data(kinshipdata)
#' res<-powerStressMin(as.matrix(kinshipdelta),kappa=2,lambda=1.5)
#' res
#' summary(res)
#' plot(res)
#' 
#' @export
powerStressMin <- function (delta, kappa=1, lambda=1, lambdamax=lambda, weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, eps = 1e-10, itmax = 100000, verbose = FALSE, stresstype=stresse1) { 
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    if(verbose>0) cat("Minimizing powerStress with kappa=",kappa,"lambda=",lambda,"\n")
    r <- kappa/2
    p <- ndim
    deltaorig <- delta
    delta <- delta^lambda
    deltaold <- delta
    delta <- delta / enorm (delta, weightmat)
    itel <- 1
    xold <- init
    if(is.null(init)) xold <- torgerson (delta, p = p)
    xnorm <- enorm(xold)
    xold <- xold/enorm(xold) 
    n <- nrow (xold)
    k <- sum(weightmat) * ((4*r)-1)*(2^(2*r))
    l <- sum(weightmat) * ((2*r)-1)*(2^r)
    dold <- sqdist (xold)
    rold <- sum(weightmat*delta * mkPower(dold,r))
    nold <- sqrt(sum(weightmat*mkPower(dold, 2 * r)))
    lold <- rold / nold
    repeat {
       by <- mkBmat(weightmat*delta * mkPower(dold, r - 1))
       cy <- mkBmat(weightmat*mkPower(dold, (2*r)-1))
       if (r>=0.5) {
           my <- by - (lold/nold)*(cy-(k*diag(n)))
           xnew <- my %*% xold
           xnorm <- enorm(xnew)
           xnew <- xnew / enorm(xnew)
       }
       if (r<0.5) {
           gy <- as.vector((by-(l*diag(n))) %*% xold)
           ey <-  kronecker(diag(p), (lold/nold)*cy)
           xnew <- matrix(secularEq(ey,gy),n,p)
       }
       dnew <- sqdist (xnew)
       rnew <- sum (weightmat*delta * mkPower(dnew,r))
       nnew <- sqrt(sum (weightmat*mkPower(dnew, 2 * r)))
       lnew <- rnew / nnew
       if (verbose>2) {
          cat (formatC (itel, width = 4, format = "d"),
          formatC (lold, digits = 10, width = 13, format = "f"),
          formatC (lnew, digits = 10, width = 13, format = "f"), "\n")
       }
       if ((itel == itmax) || ((lnew - lold) < eps)) break ()
       itel <- itel + 1
       xold <- xnew
       dold <- dnew
       lold <- lnew
     }
     attr(xnew,"dimnames")[[2]] <- paste("D",1:p,sep="")
     doutm <- (2*sqrt(sqdist(xnew)))^kappa  #fitted powered euclidean distance but times two as is convention
     deltam <- delta
     deltaorigm <- deltaorig
     deltaoldm <- deltaold
     delta <- as.dist(delta)
     deltaorig <- as.dist(deltaorig)
     deltaold <- as.dist(deltaold)
     dout <- as.dist(doutm)
     resmat <- as.matrix(delta - dout)^2
     spp <- colMeans(resmat)
     weightmatm <-weightmat
     weightmat <- as.dist(weightmatm)
     stresso <- sum(weightmat*(dout-deltaorig^lambdamax)^2) #orig stress mit max lambda
     stressr <- sum(weightmat*(dout-deltaold)^2) #raw stress
     stresse <- sum(weightmat*(dout-delta)^2) #enormed raw stress
     stress1 <- sqrt(stressr/sum(weightmat*(dout^2)))  #implicitly normed stress for original data 
     stresse1 <- sqrt(stresse/sum(weightmat*(dout^2)))  #implicitly normed stress for enormed data
     stressn <- stressr/(sum(weightmat*deltaold^2)) #normalized to the maximum stress delta^2*lambda as the normalizing constant
     stressno <- stresso/(sum(weightmat*deltaorig^(2*lambdamax)))  #normalized to the maximum lambda; perhaps for optimization
     stressb <-stressr/(sum(weightmat*deltaold^2)*sum(weightmat*dout^2))  #normalized to both d and delta for the real observations 
     stressbe <-stressr/(sum(weightmat*delta^2)*sum(weightmat*dout^2))#normalized to both d and delta for delta
     stresscor <- cor(as.vector(weightmat*dout),as.vector(weightmat*deltaold)) #correlation of fitted and observed; is this good?
     stresscore <- cor(as.vector(weightmat*dout),as.vector(weightmat*delta)) #correlation of fitted and observed
     if(verbose>1) cat("***raw stress:",stressr,"; stress1:",stress1,"; enormed stress:",stressn,"; bstress:",stressb,"; estress:",stresse,"; estress1:",stresse1,"; bestress:",stressbe,"stresscor:",stresscor,"stresscore:",stresscore,"\n")
     out <- list(delta=deltaold, obsdiss=delta, confdiss=dout, conf = xnew, pars=c(kappa,lambda), niter = itel, stress=stresstype, spp=spp, ndim=p, model="Power Stress SMACOF", call=match.call(), nobj = dim(xnew)[1], type = "Power Stress", gamma = c(lold,lnew), stress.m=stressn, stress.r=stressr, stress.n=stressn, stress.1=stress1, stress.b=stressb, stress.e=stresse,stress.e1=stresse1,stress.be=stressbe,stress.co=stresscor, deltaorig=as.dist(deltaorig),resmat=resmat)
    class(out) <- c("smacofB","smacof")
    out
 }
#' Torgerson scaling
#'
#' @param delta symmetric, numeric matrix of distances
#' @param p target space dimensions
#'
#' @export
torgerson <- function(delta, p = 2) {
    doubleCenter <- function(x) {
        n <- dim(x)[1]
        m <- dim(x)[2]
        s <- sum(x)/(n*m)
        xr <- rowSums(x)/m
        xc <- colSums(x)/n
        return((x-outer(xr,xc,"+"))+s)
    }
    z <- eigen(-doubleCenter((as.matrix (delta) ^ 2)/2))
    v <- pmax(z$values,0)
    return(z$vectors[,1:p]%*%diag(sqrt(v[1:p])))
}

#' Explicit Norm
#'
#' @param x numeric matrix 
#' @param w weight
enorm <- function (x, w=1) {
    return (sqrt (sum (w * (x ^ 2))))
}

#' Inverse explicit norm
#'
#' @param x numeric matrix 
#' @param w weight
ienorm <- function(x,w=1){
    return(sum(w*sqrt(x))^2)
}

#' Squared distances
#'
#' @param x numeric matrix
#' @export
sqdist <- function (x) {
    s <- tcrossprod (x)
    v <- diag (s)
    return (outer (v, v, "+") - 2 * s)
}

#' Auxfunction1
#'
#' @param x matrix
mkBmat <- function (x) {
    d <- rowSums (x)
    x <- -x
    diag (x) <- d
    return (x)
}

#' MakePower
#'
#' @param x matrix
#' @param r numeric (power)
#'
#' @export
mkPower<-function(x,r) {
    n<-nrow(x)
    return(abs((x+diag(n))^r)-diag(n))
}

#' Secular Equation 
#'
#' @param a matrix
#' @param b matrix
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
    rot<-uniroot(f,lower= ulw,upper= uup)$root
    cve<-beta/(eva+rot)
    return(drop(eve%*%cve))
}    

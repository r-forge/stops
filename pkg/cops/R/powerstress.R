#' Power Stress SMACOF
#'
#' An old implemenation to minimize power stress by nested majorization. So slow you should get a coffee. Possibly to be discontinued soon. 
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param kappa power of the transformation of the fitted distances; defaults to 1
#' @param lambda the power of the transformation of the proximities; defaults to 1
#' @param nu the power of the transformation for weightmat; defaults to 1 
#' @param weightmat a matrix of finite weights
#' @param init starting configuration
#' @param ndim dimension of the configuration; defaults to 2
#' @param eps numeric accuracy of the iteration
#' @param itmax maximum number of iterations
#' @param verbose should iteration output be printed; if > 1 then yes
#' @param defaultstress the stress type to be reported by default; defaults to stress-1 for smacof compatibility
#' @param lambdamax the maximum power of the transformation of the proximities if there are more than one lambda - an upper bound idea; defaults to lambda
#'
#' @return a smacofP object (inheriting form smacofB, see \code{\link{smacofSym}}). It is a list with the components
#' \itemize{
#' \item delta: Observed dissimilarities, not normalized
#' \item obsdiss: Observed dissimilarities, normalized 
#' \item confdiss: Configuration dissimilarities, NOT normalized 
#' \item conf: Matrix of fitted configuration, NOT normalized
#' \item stress: Default stress  
#' \item spp: Stress per point (based on stress.en) 
#' \item ndim: Number of dimensions
#' \item model: Name of smacof model
#' \item niter: Number of iterations
#' \item nobj: Number of objects
#' \item type: Type of MDS model
#' }
#' and some additional components
#' \itemize{
#' \item gamma: The majorizing function at convergence
#' \item stress.m: default stress for the COPS and STOP defaults to square root of the explicitly normalized stress on the normalized, transformed dissimilarities
#' \item stress.r: raw stress on the observed, transformed dissimilarities
#' \item stress.s: square root of the explicitly normalized stress on the observed, transformed dissimilarities
#' \item stress.n: explicitly normalized stress on the observed, transformed dissimilarities
#' \item stress.1: implicitly normalized stress on the observed, transformed dissimilarities
#' \item stress.e: raw stress on the normalized, transformed dissimilarities
#' \item stress.en: stress on the normalized, transformed dissimilarities and normalized transformed distances
#' \item stress.en1: stress 1 on the normalized, transformed dissimilarities and normalized transformed distances; this is reported by the print ,ethod
#' \item stress.e1: implicitly normalized stress on the normalized, transformed dissimilarities
#' \item deltaorig: observed, untransformed dissimilarities
#' \item weightmat: weighting matrix 
#'}
#'
#' @importFrom stats dist as.dist
#' 
#' @seealso \code{\link{smacofSym}}
#' 
#' @examples
#' dis<-smacof::kinshipdelta
#' res<-powerStressMin(as.matrix(dis),kappa=2,lambda=1.5)
#' res
#' summary(res)
#' plot(res)
#' 
#' @export
powerStressMinOLD <- function (delta, kappa=1, lambda=1, nu=1,lambdamax=lambda, weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, eps = 1e-10, itmax = 100000, verbose = FALSE, defaultstress=stressen1) {
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    if(verbose>0) cat("Minimizing powerStress with kappa=",kappa,"lambda=",lambda,"nu=",nu,"\n")
    r <- kappa/2
    p <- ndim
    deltaorig <- delta
    delta <- delta^lambda
    weightmato <- weightmat
    weightmat <- weightmat^nu
    weightmat[!is.finite(weightmat)] <- 1 #new
    deltaold <- delta
    delta <- delta / enorm (delta, weightmat)
    itel <- 1
    xold <- init
    if(is.null(init)) xold <- cops::torgerson (delta, p = p)
    #xnorm <- enorm(xold)
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
           #xnorm <- enorm(xnew)
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
     #xnew <- xnew/enorm(xnew)
     doutm <- (2*sqrt(sqdist(xnew)))^kappa  #fitted powered euclidean distance but times two
     deltam <- delta
     deltaorigm <- deltaorig
     deltaoldm <- deltaold
     delta <- stats::as.dist(delta)
     deltaorig <- stats::as.dist(deltaorig)
     deltaold <- stats::as.dist(deltaold)
     doute <- doutm/enorm(doutm)
     #doute <- doutm
     doute <- stats::as.dist(doute)
     dout <- stats::as.dist(doutm)
     weightmatm <-weightmat
     resmat <- weightmatm*as.matrix((delta - doute)^2)
     spp <- colMeans(resmat)
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
     if(verbose>1) cat("*** stress (both normalized):",stressen,"; stress 1 (both normalized - default reported):",stressen1,"; sqrt raw stress (both normalized):",sqrt(stressen),"; raw stress (original data):",stressr,"; stress 1 (original data):",stress1,"; explicitly normed stress (original data):",stressn,"; sqrt explicitly normed stress (original data - used in STOPS):",stresss,"; raw stress (proximities normalized):",stresse,"; stress 1 (proximities normalized):", stresse1,"\n")   
    out <- list(delta=deltaold, obsdiss=delta, confdiss=dout, conf = xnew, pars=c(kappa,lambda,nu), niter = itel, stress=defaultstress, spp=spp, ndim=p, model="Power Stress SMACOF", call=match.call(), nobj = dim(xnew)[1], type = "Power Stress", gamma = c(lold,lnew), stress.m=sqrt(stressn), stress.r=stressr/2, stress.n=stressn, stress.1=stress1, stress.s=stresss,stress.e=stresse,stress.en=stressen, stress.en1=stressen1,stress.e1=stresse1, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat)
    class(out) <- c("smacofP","smacofB","smacof")
    out
 }


#' double centering 
#'
#' @param x numeric matrix
doubleCenter <- function(x) {
        n <- dim(x)[1]
        m <- dim(x)[2]
        s <- sum(x)/(n*m)
        xr <- rowSums(x)/m
        xc <- colSums(x)/n
        return((x-outer(xr,xc,"+"))+s)
    }

#' Torgerson scaling
#'
#' @param delta symmetric, numeric matrix of distances
#' @param p target space dimensions
#' @export
torgerson <- function(delta, p = 2) {
    z <- eigen(-doubleCenter((as.matrix (delta) ^ 2)/2))
    v <- pmax(z$values,0)
    return(z$vectors[,1:p]%*%diag(sqrt(v[1:p])))
}

#' Explicit Norm
#'
#' @param x numeric matrix 
#' @param w weight
#' @export
enorm <- function (x, w=1) {
    return (sqrt (sum (w * (x ^ 2))))
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

#' Squared distances
#'
#' @param x numeric matrix
#' @param p p>0 the minkoswki distance 
#' @export
pdist <- function (x,p) {
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
mkPowerALTERN<-function(x,r) {
    n<-nrow(x)
    tmp <- abs((x+diag(n))^r)-diag(n)
    tmp[!is.finite(tmp)] <- 1e+100
    return(tmp)
}

#' MakePower Old
#'
#' @param x matrix
#' @param r numeric (power)
#'
#' @export
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

#'S3 plot method for smacofP objects
#' 
#'@param x an object of class smacofP 
#'@param plot.type String indicating which type of plot to be produced: "confplot", "resplot", "Shepard", "stressplot","transplot", "bubbleplot" (see details)
#'@param plot.dim  dimensions to be plotted in confplot; defaults to c(1, 2)
#'@param main plot title
#'@param xlab label of x axis
#'@param ylab label of y axis
#'@param xlim scale of x axis
#'@param ylim scale of y axis
#'@param col vector of colors for the points
#'@param bubscale Scaling factor (size) for the bubble plot
#'@param label.conf List with arguments for plotting the labels of the configurations in a configuration plot (logical value whether to plot labels or not, label position, label color)
#'@param identify If 'TRUE', the 'identify()' function is called internally that allows to add configuration labels by mouse click
#'@param type What type of plot should be drawn (see also 'plot')
#'@param legend Flag whether legends should be drawn for plots that have legends
#'@param legpos Position of legend in plots with legends 
#'@param pch  Plot symbol
#'@param asp  Aspect ratio; defaults to 1 so distances between x and y are represented accurately; can lead to slighlty weird looking plots if the variance on one axis is much smaller than on the other axis; use NA if the standard type of R plot is wanted where the ylim and xlim arguments define the aspect ratio - but then the distances seen are no longer accurate
#'@param loess should loess fit be added to Shepard plot 
#'@param ... Further plot arguments passed: see 'plot.smacof' and 'plot' for detailed information.
#' 
#'Details:
#' \itemize{
#' \item  Configuration plot (plot.type = "confplot"): Plots the MDS configurations.
#'  \item Residual plot (plot.type = "resplot"): Plots the dissimilarities against the fitted distances.
#'  \item Linearized Shepard diagram (plot.type = "Shepard"): Diagram with the transformed observed dissimilarities against the transformed fitted distance as well as loess curve and a least squares line.
#'  \item Transformation Plot (plot.type = "transplot"): Diagram with the observed dissimilarities (lighter) and the transformed observed dissimilarities (darker) against the fitted distances together with the nonlinear regression curve 
#'  \item Stress decomposition plot (plot.type = "stressplot"): Plots the stress contribution in of each observation. Note that it rescales the stress-per-point (SPP) from the corresponding smacof function to percentages (sum is 100). The higher the contribution, the worse the fit.
#'  \item Bubble plot (plot.type = "bubbleplot"): Combines the configuration plot with the point stress contribution. The larger the bubbles, the better the fit.
#' }
#'
#' @importFrom graphics plot text identify legend
#' @importFrom stats loess lm predict 
#' 
#' @export 
plot.smacofP <- function (x, plot.type = "confplot", plot.dim = c(1, 2), bubscale = 5, col, label.conf = list(label = TRUE, pos = 3, col = 1, cex = 0.8), identify = FALSE, type = "p", pch = 20, asp = 1, main, xlab, ylab, xlim, ylim, legend = TRUE , legpos, loess=TRUE, ...)
{
    x1 <- plot.dim[1]
    y1 <- plot.dim[2]
    if (type == "n") 
        label.conf$pos <- NULL
    if (plot.type == "confplot") {
        if(missing(col)) col <- 1
        if (missing(main)) 
            main <- paste("Configuration Plot")
        else main <- main
        if (missing(xlab)) 
            xlab <- paste("Configurations D", x1, sep = "")
        else xlab <- xlab
        if (missing(ylab)) 
            ylab <- paste("Configurations D", y1, sep = "")
        else ylab <- ylab
        if (missing(xlim)) xlim <- range(x$conf[, x1])
        if (missing(ylim)) ylim <- range(x$conf[, y1]) 
        graphics::plot(x$conf[, x1], x$conf[, y1], main = main, type = type, 
            xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, 
            pch = pch, asp = asp, col = col, ...)
        if (label.conf[[1]]) 
            graphics::text(x$conf[, x1], x$conf[, y1], labels = rownames(x$conf), 
                cex = label.conf$cex, pos = label.conf$pos, col = label.conf$col)
        if (identify) {
            graphics::identify(x$conf[, x1], x$conf[, y1], labels = rownames(x$conf), 
                cex = label.conf$cex, pos = label.conf$cex, col = label.conf$col)
        }
    }
    if (plot.type == "Shepard") {
        if(missing(col)) col <- c("grey60","grey50","black")
        if (missing(main)) 
            main <- paste("Linearized Shepard Diagram")
        else main <- main
        if (missing(xlab)) 
            xlab <- "Transformed Dissimilarities"
        else xlab <- xlab
        if (missing(ylab)) 
            ylab <- "Transformed Configuration Distances"
        else ylab <- ylab
        if (missing(xlim)) 
            xlim <- range(as.vector(x$delta))
        if (missing(ylim))
            ylim <- range(as.vector(x$confdiss))
        #delta=dhats
        #proximities=obsdiss
        #distances=confdiss
        graphics::plot(as.vector(x$delta), as.vector(x$confdiss), main = main, type = "p", pch=20, cex = 0.75, xlab = xlab, ylab = ylab, col = col[1], xlim = xlim, ylim = ylim, ...)
        #graphics::plot(as.vector(x$delta), as.vector(x$confdiss), main = main, type = "p", cex = 0.75, xlab = xlab, ylab = ylab, col = col[1], xlim = xlim, ylim = ylim)
        #graphics::points(as.vector(x$delta), ),col=col[2],pch=19)
        #graphics::plot(as.vector(x$delta), as.vector(x$obsdiss),col=col[2],pch=20)
        if(loess) {
                   pt <- predict(stats::loess(x$confdiss~-1+x$delta))
                   graphics::lines(x$delta[order(x$delta)],pt[order(x$delta)],col=col[2],type="b",pch=20,cex=0.25)
        }
        ptl <- predict(stats::lm(x$confdiss~-1+x$delta))
        graphics::lines(x$delta[order(x$delta)],ptl[order(x$delta)],col=col[3],type="b",pch=20,cex=0.25)
       # graphics::abline(stats::lm(x$confdiss~-1+x$delta),type="b") #no intercept for fitting
    }
    if (plot.type == "transplot") {
             if(missing(col)) col <- c("grey40","grey70","grey30")#,"grey50")
             kappa <- x$pars[1]
             deltao <- as.vector(x$deltaorig)
             deltat <- as.vector(x$delta)
             dreal <- as.vector(x$confdiss)^(1/kappa)
             if (missing(main)) main <- paste("Transformation Plot")
             else main <- main
             if (missing(ylab)) ylab <- "Dissimilarities"
             else xlab <- xlab
             if (missing(xlab))  xlab <- "Untransformed Configuration Distances"
             else ylab <- ylab
             if (missing(ylim))  ylim <- c(min(deltat,deltao),max(deltat,deltao))
             if (missing(xlim))  xlim <- c(min(dreal^kappa,dreal),max(dreal^kappa,dreal))
            graphics::plot(dreal, deltao, main = main, type = "p", cex = 0.75, xlab = xlab, ylab = ylab, col = col[2], xlim = xlim, ylim = ylim,pch=20)
            #graphics::plot(deltat,dreal, main = main, type = "p", cex = 0.75, xlab = ylab, ylab = xlab, col = col[2], xlim = ylim, ylim = xlim,pch=20)
            graphics::points(dreal, deltat, type = "p", cex = 0.75, col = col[1],pch=20)
            pt <- predict(stats::lm(deltat~-1+I(dreal^kappa))) #with intercept forcing thorugh 0
            #pt2 <- predict(stats::lm(deltat~I(dreal^kappa))) #with intercept not forcing thorugh 0 
            #po <- predict(stats::lm(deltao~-1+I(dreal^kappa))) #with intercept
            #lines(deltat[order(deltat)],pt[order(deltat)],col=col[1],type="b",pch=20,cex=0.5)
            #lines(deltao[order(deltao)],po[order(deltao)],col=col[2],type="b",pch=20,cex=0.5)
            #graphics::lines(dreal[order(dreal)],po[order(dreal)],col=col[4])
            graphics::lines(dreal[order(dreal)],pt[order(dreal)],col=col[3],type="b",pch=19,cex=0.1)
            #graphics::lines(dreal[order(dreal)],po[order(dreal)],col=col[4],type="b",pch=19,cex=0.25) 
            if(legend) {
                if(missing(legpos)) legpos <- "topleft" 
                graphics::legend(legpos,legend=c("Transformed","Untransformed"),col=col[1:2],pch=1)
            }
         }
     if (plot.type == "resplot") {
        if(missing(col)) col <- "darkgrey" 
        if (missing(main)) 
            main <- paste("Residual plot")
        else main <- main
        if (missing(xlab)) 
            xlab <- "Normalized Dissimilarities"
        else xlab <- xlab
        if (missing(ylab)) 
            ylab <- "Configuration Distances"
        else ylab <- ylab
        if (missing(xlim)) 
            xlim <- range(as.vector(x$obsdiss))
        if (missing(ylim)) 
            ylim <- range(as.vector(x$confdiss))
        graphics::plot(as.vector(x$obsdiss), as.vector(x$confdiss), main = main, 
            type = "p", col = col, xlab = xlab, ylab = ylab, 
            xlim = xlim, ylim = ylim, ...)
        abline(lm(x$confdiss~x$obsdiss))
    }
    if (plot.type == "stressplot") {
        if(missing(col)) col <- "lightgray"
        if (missing(main)) 
            main <- paste("Stress Decomposition Chart")
        else main <- main
        if (missing(xlab)) 
            xlab <- "Objects"
        else xlab <- xlab
        if (missing(ylab)) 
            ylab <- "Stress Proportion (%)"
        else ylab <- ylab
        spp.perc <- sort((x$spp/sum(x$spp) * 100), decreasing = TRUE)
        xaxlab <- names(spp.perc)
        if (missing(xlim)) 
            xlim1 <- c(1, length(spp.perc))
        else xlim1 <- xlim
        if (missing(ylim)) 
            ylim1 <- range(spp.perc)
        else ylim1 <- ylim
        plot(1:length(spp.perc), spp.perc, xaxt = "n", type = "p", 
            xlab = xlab, ylab = ylab, main = main, xlim = xlim1, 
            ylim = ylim1, ...)
        text(1:length(spp.perc), spp.perc, labels = xaxlab, pos = 3, cex = 0.8)
        for (i in 1:length(spp.perc)) lines(c(i, i), c(spp.perc[i],0), col=col, lty = 2)
    }
    if (plot.type == "bubbleplot") {
        if(missing(col)) col <- 1
        if (missing(main)) 
            main <- paste("Bubble Plot")
        else main <- main
        if (missing(xlab)) 
            xlab <- paste("Configurations D", x1, sep = "")
        else xlab <- xlab
        if (missing(ylab)) 
            ylab <- paste("Configurations D", y1, sep = "")
        else ylab <- ylab
        if (missing(xlim)) 
            xlim <- range(x$conf[, x1]) * 1.1
        if (missing(ylim)) 
            ylim <- range(x$conf[, y1]) * 1.1
        spp.perc <- x$spp/sum(x$spp) * 100
        bubsize <- (max(spp.perc) - spp.perc + 1)/length(spp.perc) * bubscale
        plot(x$conf, cex = bubsize, main = main, xlab = xlab, 
            ylab = ylab, xlim = xlim, ylim = ylim, ...)
        xylabels <- x$conf
        ysigns <- sign(x$conf[, y1])
        xylabels[, 2] <- (abs(x$conf[, y1]) - (x$conf[, y1] * (bubsize/50))) * ysigns
        text(xylabels, rownames(x$conf), pos = 1, cex = 0.7)
    }
 }


#'@export
summary.smacofP <- function(object,...)
    {
      spp.perc <- object$spp/sum(object$spp) * 100
      sppmat <- cbind(sort(object$spp), sort(spp.perc))
      colnames(sppmat) <- c("SPP", "SPP(%)") 
      res <- list(conf=object$conf,sppmat=sppmat)
      class(res) <- "summary.smacofP"
      res
    }

#'@export
print.summary.smacofP <- function(x,...)
    {
    cat("\n")
    cat("Configurations:\n")
    print(round(x$conf, 4))
    cat("\n\n")
    cat("Stress per point:\n")
    print(round(x$sppmat, 4))
    cat("\n")
    }


#' Power stress minimization by NEWUOA (old version)
#'
#' An implementation to minimize power stress by a derivative-free trust region optimization algorithm (NEWUOA). Much faster than majorizing but perhaps slighly less accurate. 
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param kappa power of the transformation of the fitted distances; defaults to 1
#' @param lambda the power of the transformation of the proximities; defaults to 1
#' @param nu the power of the transformation for weightmat; defaults to 1 
#' @param weightmat a matrix of finite weights
#' @param init starting configuration
#' @param ndim dimension of the configuration; defaults to 2
#' @param eps  The smallest value of the trust region radius that is allowed. If not defined, then 1e-10 will be used.
#' @param itmax maximum number of iterations
#' @param verbose should iteration output be printed; if > 1 then yes
#' @param defaultstress the stress type to be reported by default; defaults to stress-1 for smacof compatibility
#' @param lambdamax the maximum power of the transformation of the proximities if there are more than one lambda - an upper bound idea; defaults to lambda
#'
#' @return a smacofP object (inheriting form smacofB, see \code{\link{smacofSym}}). It is a list with the components
#' \itemize{
#' \item delta: Observed dissimilarities, not normalized
#' \item obsdiss: Observed dissimilarities, normalized 
#' \item confdiss: Configuration dissimilarities, NOT normalized 
#' \item conf: Matrix of fitted configuration, NOT normalized
#' \item stress: Default stress  
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
#' \item stress.m: default stress for the COPS and STOP defaults to square root of the explicitly normalized stress on the normalized, transformed dissimilarities
#' \item stress.r: raw stress on the observed, transformed dissimilarities
#' \item stress.s: square root of the explicitly normalized stress on the observed, transformed dissimilarities
#' \item stress.n: explicitly normalized stress on the observed, transformed dissimilarities
#' \item stress.1: implicitly normalized stress on the observed, transformed dissimilarities
#' \item stress.e: raw stress on the normalized, transformed dissimilarities
#' \item stress.en: stress on the normalized, transformed dissimilarities and normalized transformed distances
#' \item stress.en1: stress 1 on the normalized, transformed dissimilarities and normalized transformed distances; this is reported by the print ,ethod
#' \item stress.e1: implicitly normalized stress on the normalized, transformed dissimilarities
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
powerStressFastOLD <- function (delta, kappa=1, lambda=1, nu=1,lambdamax=lambda, weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, eps = 1e-12, itmax = 100000, verbose = FALSE, defaultstress=stressen1)
{
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    if(verbose>0) cat("Minimizing powerstress by NEWUOA with kappa=",kappa,"lambda=",lambda,"nu=",nu,"\n")
    p <- ndim
    deltaorig <- delta
    delta <- delta^lambda
    weightmato <- weightmat
    weightmat <- weightmat^nu
    weightmat[!is.finite(weightmat)] <- 1 #new
    deltaold <- delta
    delta <- delta / enorm (delta, weightmat) #sum=1
    xold <- init
    if(is.null(init)) xold <- cops::torgerson (delta, p = p)
    xold <- xold/enorm(xold) 
    stressf <- function(x,delta,kappa,p,weightmat)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=p)
             delta <- delta/enorm(delta,weightmat)
             x <- x/enorm(x)
             ds <- as.matrix((dist(x))^kappa)
             ds <- ds/enorm(ds,weightmat)
             #print(ds)
             sum(weightmat*(ds-delta)^2)#/(sum(weightmat*(ds^2)))
             #sum((ds-delta)^2)/2
             #check if stress 1 is the same as in smacof and whether config is the same for
             #using dist 
             #using matrices
             #using enorm with dist
             #using enorm with matrix
           }
     optimized <- minqa::newuoa(xold,function(par) stressf(par,delta=delta,kappa=kappa,p=p,weightmat=weightmat),control=list(maxfun=itmax,rhoend=eps,iprint=verbose))
     #optimized <- optim(xold,function(par) stressf(par,delta=delta,p=p,weightmat=weightmat),control=list(maxit=itmax))
     xnew <- matrix(optimized$par,ncol=2)
     #xnew <- optimized$par
     attr(xnew,"dimnames")[[1]] <- rownames(delta)
     itel <- optimized$feval
     attr(xnew,"dimnames")[[2]] <- paste("D",1:p,sep="")
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
     #the following stress versions differ by how the distances and proximities are normalized; either both are normalized (stressen,stressen1), only proximities are normalized (stresse, stresse1), nothing is normalized (stressr, stressn, stresss)
     stressr <- sum(weightmat*(dout-deltaold)^2) #raw stress on the observed proximities
     stresse <- sum(weightmat*(dout-delta)^2) #raw stress on the normalized proximities
     stressen <- sum(weightmat*(doute-delta)^2)*2 #raw stress on the normalized proximities and normalized distances 
     stressen1 <- sqrt(sum(weightmat*(doute-delta)^2)/sum(weightmat*(doute^2))) # sqrt stress 1 on the normalized transformed proximities and distances; we use this as the value returned by print
     stress1 <- sqrt(stressr/sum(weightmat*(dout^2)))  #stress 1 on the original proximities 
     stresse1 <- sqrt(stresse/sum(weightmat*(dout^2)))  #stress 1 on the normalized proximities
     stressn <- stressr/(sum(weightmat*deltaold^2)) #normalized to the maximum stress delta^2*lambda as the normalizing constant (was defualt until v. 0.0-16)
     stresss <- sqrt(stressn) #sqrt of stressn
     if(verbose>1) cat("*** stress (both normalized):",stressen,"; stress 1 (both normalized - default reported):",stressen1,"; sqrt raw stress (both normalized):",sqrt(stressen),"; raw stress (original data):",stressr,"; stress 1 (original data):",stress1,"; explicitly normed stress (original data):",stressn,"; sqrt explicitly normed stress (original data - used in STOPS):",stresss,"; raw stress (proximities normalized):",stresse,"; stress 1 (proximities normalized):", stresse1,"; from optimization: ",optimized$fval,"\n")   
    out <- list(delta=deltaold, obsdiss=delta, confdiss=dout, conf = xnew, pars=c(kappa,lambda,nu), niter = itel, stress=defaultstress, spp=spp, ndim=p, model="Power Stress NEWUOA", call=match.call(), nobj = dim(xnew)[1], type = "Power Stress", gamma = NA, stress.m=stressen, stress.r=stressr/2, stress.n=stressn, stress.1=stress1, stress.s=stresss,stress.e=stresse,stress.en=stressen, stress.en1=stressen1,stress.e1=stresse1, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat)
    class(out) <- c("smacofP","smacofB","smacof")
    out
}

#' Power stress minimization by NEWUOA
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
#' @param acc  The smallest value of the trust region radius that is allowed. If not defined, then 1e-10 will be used.
#' @param itmax maximum number of iterations
#' @param verbose should iteration output be printed; if > 1 then yes
#'
#' @return a smacofP object (inheriting form smacofB, see \code{\link{smacofSym}}). It is a list with the components
#' \itemize{
#' \item delta: Observed dissimilarities, not normalized
#' \item obsdiss: Observed dissimilarities, normalized 
#' \item confdiss: Configuration dissimilarities, NOT normalized 
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
powerStressFast <- function (delta, kappa=1, lambda=1, nu=1, weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc = 1e-12, itmax = 100000, verbose = FALSE)
{
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    if(verbose>0) cat("Minimizing powerstress by NEWUOA with kappa=",kappa,"lambda=",lambda,"nu=",nu,"\n")
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
    if(is.null(init)) xold <- cops::torgerson (delta, p = ndim)
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
     optimized <- minqa::newuoa(xold,function(par) stressf(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat),control=list(maxfun=itmax,rhoend=acc,iprint=verbose))
     #optimized <- optim(xold,function(par) stressf(par,delta=delta,p=p,weightmat=weightmat),control=list(maxit=itmax))
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
    out <- list(delta=deltaold, obsdiss=delta, confdiss=dout, conf = xnew, pars=c(kappa,lambda,nu), niter = itel, stress=sqrt(stress), spp=spp, ndim=p, model="Power Stress NEWUOA", call=match.call(), nobj = dim(xnew)[1], type = "Power Stress", gamma = NA, stress.m=stress, stress.en=stressen, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat)
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
#' @examples
#' dis<-smacof::kinshipdelta
#' res<-powerStressMin(as.matrix(dis),kappa=2,lambda=1.5,verbose=3)
#' res
#' summary(res)
#' plot(res)
#' 
#' @export
powerStressMin <- function (delta, kappa=1, lambda=1, nu=1, weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc= 1e-10, itmax = 100000, verbose = FALSE) {
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
#' @examples
#' dis<-smacof::kinshipdelta
#' system.time(powerStressMin2(as.matrix(dis),kappa=2,lambda=1.5))
#' res2
#' summary(res2)
#' plot(res)
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

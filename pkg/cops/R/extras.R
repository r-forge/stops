#' Wrapper to \code{cmdscale} for S3 class
#'
#' @details overloads base::cmdscale and adds class attributes for which there are methods. The functionality is duplicated in the stops package.    
#'  
#' @param d a distance structure such as that returned by 'dist' or a full symmetric matrix containing the dissimilarities
#' @param k the maximum dimension of the space which the data are to be represented in
#' @param eig indicates whether eigenvalues should be returned.
#' @param ... additional parameters passed to cmdscale. See \code{\link{cmdscale}} 
#'
#' @return Object of class "cmdscaleE' and 'cmdscale' extending \code{\link{cmdscale}}. This wrapper only adds an extra slot to the list with the call, adds column labels to the $points and assigns S3 class 'cmdscaleE' and 'cmdscale'.
#'
#' @importFrom stats cmdscale as.dist dist
#'

#' 
#' @export
#' @examples
#' dis<-as.matrix(smacof::kinshipdelta)
#' res<-cmdscale(dis)
cmdscale <- function(d,k=2,eig=TRUE,...)
    {
     out <- stats::cmdscale(d,k=k,eig=eig,...)
     colnames(out$points) <- paste("D",1:k,sep="")
     out$call <- match.call()
     out$delta <- as.dist(d) #FIX
     out$confdist <- dist(out$points)
     class(out) <- c("cmdscaleE","cmdscale")
     out
 }

#'Wrapper to \code{sammon} for S3 class
#'
#' @details overloads MASS::sammon and adds class attributes for which there are methods. The functionality is duplicated in the stops package.    
#' 
#' @param d a distance structure such as that returned by 'dist' or a full symmetric matrix.  Data are assumed to be dissimilarities or relative distances, but must be positive except for self-distance.  This can contain missing values.
#' @param y An initial configuration. If NULL, 'cmdscale' is used to provide the classical solution.  (If there are missing values in 'd', an initial configuration must be provided.)  This must not have duplicates.
#' @param k The dimension of the configuration
#' @param ... Additional parameters passed to \code{sammon}, see \code{\link{sammon}}  
#'
#' @return Object of class sammonE' inheriting from \code{\link{sammon}}. This wrapper only adds an extra slot to the list with the call, adds column labels to the $points and assigns S3 classes 'sammonE', 'sammon' and 'cmdscale'. It also adds a slot obsdiss with normalized dissimilarities.
#'
#' @importFrom MASS sammon
#' @importFrom stats as.dist dist 
#' 
#' @export
#' @examples
#' dis<-as.matrix(smacof::kinshipdelta)
#' res<-sammon(dis)
sammon <- function(d,y=NULL,k=2,...)
    {
     if(is.null(y)) y <- cops::cmdscale(d,k,eig=TRUE)$points
     out <- MASS::sammon(d,y=as.matrix(y),k=k,...)
     colnames(out$points) <- paste("D",1:k,sep="") 
     out$call <- match.call()
     out$delta <- as.dist(d)
     out$obsdiss <- as.dist(d/enorm(d))
     out$confdist <- dist(out$points)
     out$stress.m <- sqrt(out$stress)
     class(out) <- c("sammonE","sammon","cmdscale")
     out
 }

#'@export
print.sammon <- function(x,...)
    {
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    cat("Model: Sammon Scaling \n")
    cat("Number of objects:", dim(x$points)[1], "\n")
    cat("Stress:", x$stress, "\n")
    cat("\n")
    }

#'@export
print.cmdscale <- function(x,...)
    {
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    cat("Model: Torgerson-Gower Scaling \n")
    cat("Number of objects:", dim(x$points)[1], "\n")
    cat("GOF:", x$GOF, "\n")
    cat("\n")
    }

#'@export
summary.cmdscale <- function(object,...)
    {
    cat("\n")
    cat("Configurations:\n")
    print(round(object$points, 4))
    cat("\n\n")
    }

#'@export
summary.sammon <- function(object,...)
    {
    cat("\n")
    cat("Configurations:\n")
    print(round(object$points, 4))
    cat("\n\n")
    }


#'@importFrom graphics plot abline lines text identify legend points
#'@importFrom stats predict loess lm
#'@export
plot.cmdscale <- function(x, plot.type = c("confplot"), plot.dim = c(1, 2), col, label.conf = list(label = TRUE, pos = 3, col = 1, cex = 0.8), identify = FALSE, type = "p", pch = 20, asp = 1, main, xlab, ylab, xlim, ylim, legpos,...)
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
        if (missing(xlim)) 
            xlim <- range(x$points[, x1])
        if (missing(ylim)) 
            ylim <- range(x$points[, y1])
        graphics::plot(x$points[, x1], x$points[, y1], main = main, type = type, 
            xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, 
            pch = pch, asp = asp, col = col, ...)
        if (label.conf[[1]]) 
            graphics::text(x$points[, x1], x$points[, y1], labels = rownames(x$points), 
                cex = label.conf$cex, pos = label.conf$pos, col = label.conf$col)
        if (identify) {
            graphics::identify(x$points[, x1], x$points[, y1], labels = rownames(x$points), 
                cex = label.conf$cex, pos = label.conf$cex, col = label.conf$col)
        }
    }
    if (plot.type %in% c("Shepard","resplot")) {
        delt <- as.vector(x$delta) #new in 1-2.0
        confd <- as.vector(x$confdist) #new in 1-2.0
        if(missing(col)) col <- c("grey60","grey65") 
        if (missing(main)) 
            main <- ifelse(plot.type=="Shepard",paste("Linearized Shepard Diagram"),paste("Residual plot"))
        else main <- main
        if (missing(xlab)) 
            xlab <- "Transformed Dissimilarities"
        else xlab <- xlab
        if (missing(ylab)) 
            ylab <- "Transformed Configuration Distances"
        else ylab <- ylab
        if (missing(xlim)) 
            xlim <- range(delt)
        if (missing(ylim)) 
            ylim <- range(confd)
        graphics::plot(delt, confd, main = main, 
            type = "p", pch = ifelse(plot.type=="Shepard",20,1), cex = ifelse(plot.type=="Shepard",0.75,1), xlab = xlab, ylab = ylab, 
            col = col[1], xlim = xlim, ylim = ylim, ...)
        if(plot.type=="Shepard") {
           pt <- predict(stats::loess(confd~delt))
           graphics::lines(delt[order(delt)],pt[order(delt)],col=col[2],type="b",pch=20,cex=0.5)
        }
     graphics::abline(stats::lm(confd~delt))
    }
    if (plot.type == "transplot") {
             if(missing(col)) col <- c("grey40","grey70","grey30","grey60")
             if(is.null(x$lambda)) x$lambda <- 1
             deltao <- as.vector(x$delta^(1/x$lambda))
             deltat <- as.vector(x$delta)
             dreal <- as.vector(x$confdist)
             if (missing(main)) main <- paste("Transformation Plot")
             else main <- main
             if (missing(ylab)) ylab <- "Dissimilarities"
             else ylab <- ylab
             if (missing(xlab)) xlab <- "Configuration Distances"
             else xlab <- xlab
             if (missing(ylim)) ylim <- c(min(deltat,deltao),max(deltat,deltao))
             if (missing(xlim)) xlim <- range(as.vector(dreal))
            graphics::plot(dreal, deltao, main = main, type = "p", cex = 0.75, xlab = xlab, ylab = ylab, col = col[2], xlim = xlim, ylim = ylim, ...)
            graphics::points(dreal, deltat, type = "p", cex = 0.75, col = col[1])
            pt <- predict(stats::lm(deltat~dreal))
            po <- predict(stats::lm(deltao~dreal))
            graphics::lines(dreal[order(dreal)],pt[order(dreal)],col=col[3])
            graphics::lines(dreal[order(dreal)],po[order(dreal)],col=col[4])
            if(missing(legpos)) legpos <- "topleft" 
            graphics::legend(legpos,legend=c("Transformed","Untransformed"),col=col[1:2],pch=1)
         }
    if (plot.type == "stressplot") {
        warning("Plot not supported for this class. Please use a stress-based MDS.")
    }
    if (plot.type == "bubbleplot") {
        warning("Plot not supported for this class. Please use a stress-based MDS.")
    }
    }


#' A static 3d plot S3 generic 
#'
#' A static 3d plot 
#'
#' @title plot3dstatic: static 3D plots 
#' @param x object
#' @param plot.dim dimensions to plot
#' @param main main title
#' @param xlab label for x axis 
#' @param ylab label for y axis
#' @param zlab label for z axis
#' @param col color
#' @param ... other arguments
#' 
#' @export
plot3dstatic <- function(x, plot.dim = c(1,2,3), main, xlab, ylab, zlab, col, ...) UseMethod("plot3dstatic")

#' 3D plots: plot3dstatic method for class cmdscale
#'
#' 
#' This methods produces a static 3D configuration plot.
#' @param x object of class cmdscale
#' @param plot.dim vector of length 3 with dimensions to be plotted
#' @param main plot title
#' @param xlab label of x axis
#' @param ylab label of y axis
#' @param zlab label of z axis
#' @param col color of the text labels
#' @param ... Further plot arguments passed: see 'scatterplot3d' in package 'scatterplot3d' for detailed information.
#'
#'@export
#'@import scatterplot3d
plot3dstatic.cmdscaleE <- function (x, plot.dim = c(1, 2, 3), main, xlab, ylab, zlab, col,...) 
{
    ndim <- dim(x$points)[2]
    if (ndim < 3) 
        stop("No 3D plots can be drawn for ndim < 3 !")
    if (length(plot.dim) != 3) 
        stop("plot.dim must be of length 3!")
    pd1 <- plot.dim[1]
    pd2 <- plot.dim[2]
    pd3 <- plot.dim[3]
    if (pd3 > ndim) 
        stop("Only", ndim, "dimensions were extracted!")
    if (missing(xlab)) 
        xlab <- paste("Dimension", pd1)
    if (missing(ylab)) 
        ylab <- paste("Dimension", pd2)
    if (missing(zlab)) 
        zlab <- paste("Dimension", pd3)
    x1 <- x$points[, pd1]
    y1 <- x$points[, pd2]
    z1 <- x$points[, pd3]
    if (missing(main)) 
        main1 <- "Configuration Plot"
    else main1 <- main
    if (missing(col)) 
        col1 <- "blue"
    else col1 <- col
    pr <- scatterplot3d::scatterplot3d(x1, y1, z1, type = "n", main = main1, xlab = xlab, ylab = ylab, zlab = zlab, ...)
    text(pr$xyz.convert(x1, y1, z1), labels = rownames(x$points), col = col1)
}

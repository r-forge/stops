#'Wrapper to \code{cmdscale} for S3 class
#'
#' @param d a distance structure such as that returned by 'dist' or a full symmetric matrix containing the dissimilarities
#' @param k the maximum dimension of the space which the data are to be represented in
#' @param eig indicates whether eigenvalues should be returned.
#' @param ... additional parameters passed to cmdscale. See \code{\link{cmdscale}} 
#'
#' @return See \code{\link{cmdscale}}. This wrapper only adds an extra slot to the list with the call, adds column labels to the $points and assigns S3 class 'cmdscale'   
#' @export
cmdscale <- function(d,k=2,eig=TRUE,...)
    {
     out <- stats::cmdscale(d,k=k,eig=eig,...)
     colnames(out$points) <- paste("D",1:k,sep="")
     out$call <- match.call()
     class(out) <- "cmdscale"
     out
 }

#'Wrapper to \code{sammon} for S3 class
#'
#' @param d a distance structure such as that returned by 'dist' or a full symmetric matrix.  Data are assumed to be dissimilarities or relative distances, but must be positive except for self-distance.  This can contain missing values.
#' @param y An initial configuration. If NULL, 'cmdscale' is used to provide the classical solution.  (If there are missing values in 'd', an initial configuration must be provided.)  This must not have duplicates.
#' @param k The dimension of the configuration
#' @param ... Additional parameters passed to \code{sammon}, see \code{\link{sammon}}  
#'
#' @return See \code{\link{sammon}}. This wrapper only adds an extra slot to the list with the call, adds column labels to the $points and assigns S3 classes 'sammon', 'cmdscale'   
#'@export
sammon <- function(d,y=NULL,k=2,...)
    {
     if(is.null(y)) y <- stops::cmdscale(d,k,eig=TRUE)$points
     out <- MASS::sammon(d,y=as.matrix(y),k=k,...)
     colnames(out$points) <- paste("D",1:k,sep="") 
     out$call <- match.call()
     class(out) <- c("sammon","cmdscale")
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

#'@export
plot.cmdscale <- function(x, plot.type = "confplot", plot.dim = c(1, 2), col = 1, label.conf = list(label = TRUE, pos = 3, col = 1, cex = 0.8), identify = FALSE, type = "p", pch = 20, asp = 1, main, xlab, ylab, xlim, ylim, ...)
    {
    x1 <- plot.dim[1]
    y1 <- plot.dim[2]
    if (type == "n") 
        label.conf$pos <- NULL
    if (plot.type == "confplot") {
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
        plot(x$points[, x1], x$points[, y1], main = main, type = type, 
            xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, 
            pch = pch, asp = asp, col = col, ...)
        if (label.conf[[1]]) 
            text(x$points[, x1], x$points[, y1], labels = rownames(x$points), 
                cex = label.conf$cex, pos = label.conf$pos, col = label.conf$col)
        if (identify) {
            identify(x$points[, x1], x$points[, y1], labels = rownames(x$points), 
                cex = label.conf$cex, pos = label.conf$cex, col = label.conf$col)
        }
    }
    if (plot.type == "Shepard") {
        plot(-10:10,-10:10,type="n",axes=FALSE, xlab="",ylab="")
        replicate(10,text(runif(1,-10,10),runif(1,-10,10),"NOT SUPPORTED. USE SMACOF!",cex=runif(1,max=3)))
    }
    if (plot.type == "resplot") {
        plot(-10:10,-10:10,type="n",axes=FALSE, xlab="",ylab="")
        replicate(10,text(runif(1,-10,10),runif(1,-10,10),"NOT SUPPORTED. USE SMACOF!",cex=runif(1,max=3)))
    }
    if (plot.type == "stressplot") {
        plot(-10:10,-10:10,type="n",axes=FALSE, xlab="",ylab="")
        replicate(10,text(runif(1,-10,10),runif(1,-10,10),"NOT SUPPORTED. USE SMACOF!",cex=runif(1,max=3)))
    }
    if (plot.type == "bubbleplot") {
        plot(-10:10,-10:10,type="n",axes=FALSE, xlab="",ylab="")
        replicate(10,text(runif(1,-10,10),runif(1,-10,10),"NOT SUPPORTED. USE SMACOF!",cex=runif(1,max=3)))
    }
}

#' 3D plots: plot3d method for class cmdscale
#'
#' This methods produces a dynamic 3D configuration plot.
#'
#' 
#' @param x object of class cmdscale
#' @param plot.dim vector of length 3 with dimensions to be plotted
#' @param xlab label of x axis
#' @param ylab label of y axis
#' @param zlab label of z axis
#' @param col color of the text labels
#' @param main plot title
#' @param bgpng Background image from rgl library; 'NULL' for white background
#' @param ax.grid If 'TRUE', axes grid is plotted.
#' @param sphere.rgl If 'TRUE', rgl sphere (background) is plotted.
#' @param ... Further plot arguments passed: see 'plot3d' in package 'rgl' for detailed information.
#' 
#'@export
#'@import rgl
plot3d.cmdscale <- function (x, plot.dim = c(1, 2, 3), xlab, ylab, zlab, col, main, bgpng = NULL, ax.grid = TRUE, sphere.rgl = FALSE,...) 
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
    x1 <- x$points[, pd1]
    y1 <- x$points[, pd2]
    z1 <- x$points[, pd3]
    if (missing(xlab)) 
        xlab <- paste("Dimension", pd1)
    if (missing(ylab)) 
        ylab <- paste("Dimension", pd2)
    if (missing(zlab)) 
        zlab <- paste("Dimension", pd3)
    if (is.null(bgpng)) {
        texture1 <- NULL
    }
    else {
        texture1 <- system.file(paste("textures/", bgpng, sep = ""), 
            package = "rgl")
    }
    if (missing(main)) 
        main1 <- "Configuration Plot"
    else main1 <- main
    if (missing(col)) 
        col1 <- "blue"
    else col1 <- col
    open3d()
    bg3d(sphere = sphere.rgl, texture = texture1, back = "filled", 
        color = "white")
    text3d(x1, y1, z1, texts = rownames(x$points), col = col1, 
        alpha = 1, ...)
    axes3d(c("x", "y", "z"), labels = TRUE, color = "black", 
        alpha = 1)
    title3d(xlab = xlab, ylab = ylab, zlab = zlab, main = main1, 
        color = "black", alpha = 1)
}

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
#'@import smacof
plot3dstatic.cmdscale <- function (x, plot.dim = c(1, 2, 3), main, xlab, ylab, zlab, col,...) 
{
    ndim <- dim(x$points)[2]
    options(locatorBell = FALSE)
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
#'procruster: a procrustes function 
#'
#'@param x mumeric matrix
#'
#'@export
procruster <- function (x) 
{
    sx <- svd(x)
    return(tcrossprod(sx$u, sx$v))
}

#'conf_adjust: a function to procrustes adjust two matrices
#'
#'@param conf1 reference configuration, a numeric matrix
#'@param conf2 another configuration, a numeric matrix 
#'@param verbose should adjustment be output; default to FALSE
#'@param eps numerical accuracy
#'@param itmax maximum number of iterations
#' 
#'@export
conf_adjust<- function(conf1,conf2,verbose = FALSE,eps = 1e-12, itmax = 100)
 {
x0 <- conf1
n <- nrow(x0)
ndim <- 2
metric <- TRUE
xx <- conf2
kk <- diag(ndim)
cc <- matrix(0, n, ndim)
bb <- matrix(0, n, ndim)
yy <- xx
oloss <- Inf
itel <- 1
repeat {
y0 <- matrix(0, n, ndim)
y0 <- y0 + xx %*% kk
y0 <- ((n - 1) * y0)/(n * (n - 2))
zz <- matrix(0, n, ndim)
zz <- zz + xx %*% kk
xz <- crossprod(xx, zz)
kk <- procruster(xz)
nloss <- 0
 for (i in 1:n) {
   yy <- xx %*% kk    
   yy[i,] <- n * y0[i, ]/(n - 1)
   yy <- yy - outer(rep(1, n), y0[i, ]/(n - 1))
   nloss <- nloss + sum((y0 - yy)^2)
  }
        if (verbose) {
            cat("Iteration: ", formatC(itel, digits = 3, width = 3), 
                "Old Loss: ", formatC(oloss, digits = 10, width = 15, 
                  format = "f"), "New Loss: ", formatC(nloss, 
                  digits = 10, width = 15, format = "f"), "\n")
        }
        if (((oloss - nloss) < eps) || (itel == itmax)) {
            (break)()
        }
        itel <- itel + 1
        oloss <- nloss
}
x0 <- x0 %*% procruster(crossprod(x0,y0))
result <- list(ref.conf = x0, other.conf = yy, comparison.conf = y0)
result
}

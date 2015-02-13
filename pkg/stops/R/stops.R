#' High Level STOPS Function 
#'
#' Currently only COPS, the STOPS model for c-clusteredness is implemented
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param structure what structuredness should be considered 
#' @param ... additional arguments to be passed to the work horses
#
#' @return see \code{\link{cops}}
#' 
#' @examples
#' data(BankingCrisesDistances)
#' res1<-stops(BankingCrisesDistances[,1:69],structure="clusteredness",loss="strain",verbose=0)
#' res1
#'
#' @keywords clustering multivariate
#' @export
stops <- function(dis,structure="clusteredness",...)
    {
     if(structure=="clusteredness")
         {
             out <- cops(dis,...)
             out$call <- match.call()
             return(out)
         } else
             {
                 cat(" Only c-clusteredness STOPS models (COPS) have yet been implemented \n")
             return(NULL)   
             }
 }

    
#'@export
plot.stops <- function(x,plot.type="confplot",...)
    {
     plot(x$fit,plot.type=plot.type,...)
    }

#' 3D plots: plot3d method for class stops
#'
#' 
#' This methods produces a dynamic 3D configuration plot.
#' @param x object of class stops
#' @param ... Further plot arguments to the method of the class of slot fit, see \code{\link{plot3d.smacof}} or \code{\link{plot3d.cmdscale}} . Also see 'rgl' in package 'rgl' 
#'
#'
#' 
#'@export
#'@import rgl
plot3d.stops <- function(x,...)
    {
        plot3d(x$fit,...)
    }

#' 3D plots: plot3dstatic method for class stops
#' 
#' This methods produces a static 3D configuration plot.
#' @param x object of class stops
#' @param ... Further plot arguments to the method of the class of slot fit, see \code{\link{plot3dstatic.smacof}} or \code{\link{plot3dstatic.cmdscale}} . Also see 'scatterplot3d' in package 'scatterplot3d'
#'
#'@export
#'@import smacof
plot3dstatic.stops <- function(x,...)
    {
        plot3dstatic(x$fit,...)
    }

#'@export
residuals.stops <- function(object,...)
    {
    residuals(object$fit,...)
    }

#'@export
print.stops <- function(x,...)
    {
    print(x,...)
    }

#'@export
coef.stops <- function(object,...)
    {
     coef(object,...)
    }

#'@export
summary.stops <- function(object,...)
    {
    summary(object,...)
    }

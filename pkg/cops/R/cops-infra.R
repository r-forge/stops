
#' High Level COPS Function
#'
#' About the function cops: The high level function allows for minimizing copstress for a clustered MDS configuration. Allows to choose COPS-C (finding a configuration from copstress with cordillera penalty) and profile COPS (finding hyperparameters for MDS models with power transformations). It is a wrapper for copstressMin and pcops.
#'
#'@param dis a dissimilarity matrix or a dist object
#'@param variant a character string specifying which variant of COPS to fit. Allowed is any of the following "1","2","Variant1","Variant2","v1","v2","COPS-C","P-COPS","configuration-c","profile","copstress-c","p-copstress". Defaults to "COPS-C".
#'@param ... arguments to be passed to  \code{\link{copstressMin}} (for Variant 1) or \code{\link{pcops}} (for Variant 2).
#'
#'@return For COPS-C Variant 1 see \code{\link{copstressMin}}, for P-COPS Variant 2 see \code{\link{pcops}}
#' 
#'@examples
#'dis<-as.matrix(smacof::kinshipdelta)
#'
#'#COPS-C with equal weight to stress and cordillera 
#'res1<-cops(dis,variant="COPS-C",stressweight=0.75,cordweight=0.25,
#'           minpts=2,itmax=1000) #use higher itmax in real
#'res1
#'summary(res1)
#'plot(res1)
#'plot(res1,"reachplot")
#' 
#'
#'\donttest{
#'#s-stress type copstress (i.e. kappa=2, lambda=2)
#'res3<-cops(dis,variant="COPS-C",kappa=2,lambda=2,stressweight=0.5,cordweight=0.5) 
#'res3
#'summary(res3)
#'plot(res3)
#'
#' 
#'# power-stress type profile copstress
#'# search for optimal kappa and lambda between
#'# kappa=0.5,lambda=0.5 and kappa=2,lambda=5
#'# nu is fixed on -1
#'ws<-1/dis
#'diag(ws)<-1 
#'res5<-cops(dis,variant="P-COPS",loss="powerstress",
#'           theta=c(1.4,3,-1), lower=c(1,0.5,-1),upper=c(3,5,-1),
#'           weightmat=ws, stressweight=0.9,cordweight=0.1) 
#'res5
#'summary(res5)
#'plot(res5)
#'}
#'
#' 
#' 
#'@keywords clustering multivariate
#'@export
cops <- function(dis, variant=c("1","2","Variant1","Variant2","v1","v2",
        "COPS-C","P-COPS","configuration-c","profile",
        "copstress-c","p-copstress","COPS-P",
        "copstress-p","cops-c","p-cops","copsc",
        "pcops"),...)
{
    #variant=c("0","1","2","Variant0","Variant1","Variant2","v0","v1","v2","COPS-0","COPS-C","P-COPS","configuration-0","configuration-c","profile","copstress-0","copstress-c","p-copstress","COPS-P","copstress-p","cops-c","p-cops")
                 if(missing(variant)) variant <- "COPS-C"
                 if(variant%in%c("1","Variant1","v1","configuration-c","COPS-C","copstress-c","cops-c","copsc")) out <- copstressMin(dis,...)
#                 if(variant%in%c("0","Variant0","v0","configuration-0","COPS-0","copstress-c")) stop("Not yet implemented") out <- shrinkCopstress0(dis,...)
                 if(variant%in%c("2","Variant2","v2","profile","p-copstress","P-COPS","COPS-P","p-cops","copstress-p","pcops")) out <- pcops(dis,...)
                 return(out)
                 }



#'S3 plot method for cops objects
#' 
#'@param x an object of class cops
#'@param plot.type String indicating which type of plot to be produced: "confplot", "reachplot", "resplot","transplot", "Shepard", "stressplot" (see details)
#'@param main the main title of the plot
#'@param asp aspect ratio of x/y axis; defaults to NA; setting to 1 will lead to an accurate represenation of the fitted distances. 
#'@param ... Further plot arguments passed: see 'plot.smacof' and 'plot' for detailed information.
#' 
#'Details:
#' \itemize{
#' \item Configuration plot (plot.type = "confplot"): Plots the MDS configurations.
#' \item Reachability plot (plot.type = "confplot"): Plots the OPTICS reachability plot and the OPTICS cordillera 
#' \item Residual plot (plot.type = "resplot"): Plots the dissimilarities against the fitted distances.
#' \item Linearized Shepard diagram (plot.type = "Shepard"): Diagram with the transformed observed dissimilarities against the transformed fitted distance as well as loess smooth and a least squares line.
#' \item Transformation Plot (plot.type = "transplot"): Diagram with the observed dissimilarities (lighter) and the transformed observed dissimilarities (darker) against the fitted distances together with loess smoothing lines 
#' \item Stress decomposition plot (plot.type = "stressplot", only for SMACOF objects in $fit): Plots the stress contribution in of each observation. Note that it rescales the stress-per-point (SPP) from the corresponding smacof function to percentages (sum is 100). The higher the contribution, the worse the fit.
#' \item Bubble plot (plot.type = "bubbleplot", only available for SMACOF objects $fit): Combines the configuration plot with the point stress contribution. The larger the bubbles, the better the fit.
#'} 
#'@export
#'@examples
#'dis<-as.matrix(smacof::kinshipdelta)
#'set.seed(1)
#'resl<-copstressMin(dis,itmax=500)
#'plot(resl)
plot.cops <- function(x,plot.type=c("confplot"), main, asp=1,...)
    {
     if(missing(plot.type)) plot.type <- "confplot"  
     if(plot.type=="reachplot") {
        if(missing(main)) main <- paste("Reachability Plot")
        else main <- main
        plot(x$OC,main=main,...)
     } else if(inherits(x$fit,"smacofB") && !inherits(x$fit,"smacofP") && plot.type=="transplot"){
       if(missing(main)) main <- paste("Transformation Plot") 
       plot.smacofP(x,plot.type="transplot",asp=asp,...)
     # invisible(tmp) #I give the changed smacof object back
     }
     else {      
       plot.smacofP(x,plot.type=plot.type,main=main,asp=asp,...)
   }
 }


#'@importFrom stats coef
#'@export 
coef.cops <- function(object,...)
{
  return(object$parameters)
 }


#'@export
print.cops <- function(x,...)
    {
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    cat("Model: COPS with parameter vector=",x$parameters,"\n")
    cat("\n")
    cat("Number of objects:", x$nobj, "\n")
    cat("Stress of configuration (default normalization):", x$stress, "\n")
    cat("OPTICS Cordillera: Raw", x$OC$raw,"Normed", x$OC$normed,"\n")
    cat("Cluster optimized loss (copstress): ", x$copstress, "\n")
    cat("Stress weight:",x$stressweight," OPTICS Cordillera weight:",x$cordweight,"\n")
    cat("Number of iterations of",x$optimethod,"optimization:", x$niter, "\n")
    cat("\n")
    }

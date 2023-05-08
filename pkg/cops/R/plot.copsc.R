#'S3 plot method for cops objects
#' 
#'@param x an object of class cops
#'@param plot.type String indicating which type of plot to be produced: "confplot", "reachplot", "resplot","transplot", "Shepard", "stressplot", "bubblepot", "histogram"  (see details)
#'@param main the main title of the plot
#'@param asp aspect ratio of x/y axis; defaults to 1; setting to 1 will lead to an accurate representation of the fitted distances. 
#'@param ... Further plot arguments passed: see 'plot.smacofP' and 'plot' for detailed information.
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
#'
#' #@importFrom smacofx plot
#' 
#'@export
#'@examples
#'dis<-as.matrix(smacof::kinshipdelta)
#'set.seed(1)
#'resl<-copstressMin(dis,itmax=500)
#'plot(resl)
plot.copsc <- function(x,plot.type, main, asp=1,...)
    {
     if(missing(plot.type)) plot.type <- "confplot"  
     if(plot.type=="reachplot") {
        if(missing(main)) main <- paste("Reachability Plot")
        else main <- main
        plot(x$OC,main=main,...)
     }
     else {
         class(x) <- c("smacofP",class(x))
         plot(x,plot.type=plot.type,main=main,asp=asp,...)
         class(x) <- class(x)[-1]
   }
 }

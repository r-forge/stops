

#'S3 plot method for p-cops objects
#' 
#' @param x an object of class cops
#' @param plot.type String indicating which type of plot to be produced: can be one of "confplot", "reachplot", "resplot", "transplot", "Shepard", "stressplot", "bubbleplot", "histogram" *see plot.smacofP, plot.smacofB and plot.copsc). Note that not all plots might be available for all losses.  
#' @param main the main title of the plot
#' @param asp aspect ratio of x/y axis; defaults to 1; setting to 1 will lead to an accurate represenation of the fitted distances. 
#' @param ... Further plot arguments passed: see 'plot.smacofP', 'plot.smacofB' and 'plot' for detailed information.
#' 
#' @details
#' See plot.smacofP
#'
#' #@importFrom smacofx plot
#' 
#' @export
#' @examples
#' dis<-as.matrix(smacof::kinshipdelta)
#' resl<-pcops(dis,loss="stress",lower=0.1,upper=5,minpts=2)
#' plot(resl,plot.type="confplot")
#' plot(resl,plot.type="reachplot")
#' plot(resl,plot.type="Shepard")
#' plot(resl,plot.type="transplot")
#' plot(resl,plot.type="stressplot")
#' plot(resl,plot.type="bubbleplot")
#' plot(resl,plot.type="histogram")
plot.pcops <- function(x, plot.type, main, asp=1,...)
    {
     if(missing(plot.type)) plot.type <- "confplot"  
     if(plot.type=="reachplot") {
        if(missing(main)) main <- paste("Reachability plot")
        else main <- main
        plot(x$OC,main=main,...)
     } else if(inherits(x$fit,"smacofB") && !inherits(x$fit,"smacofP") && plot.type=="transplot"){
         if(missing(main)) main <- paste("Transformation Plot")
         class(x$fit) <- c("cmdscalex",class(x))
         plot(x$fit,plot.type="transplot",main=main,...)
         class(x$fit) <- class(x$fit)[-1]
     }
     else {      
       plot(x$fit,plot.type=plot.type,main=main,asp=asp,...)
   }
 }

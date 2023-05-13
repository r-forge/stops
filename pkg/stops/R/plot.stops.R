#'S3 plot method for stops objects
#' 
#'@param x an object of class stops
#'@param plot.type String indicating which type of plot to be produced: "confplot", "resplot", "Shepard", "stressplot", "bubbleplot" (see details)
#'@param main the main title of the plot
#'@param asp aspect ratio of x/y axis; defaults to 1; setting to 1 will lead to an accurate represenation of the fitted distances. 
#'@param ... Further plot arguments passed: see 'plot.smacof' and 'plot' for detailed information.
#' 
#'Details:
#' See plot.smacofP
#'
#'@return no return value, just plots
#' 
#'@importFrom graphics plot 
#'@export 
plot.stops <- function(x,plot.type="confplot", main, asp=1,...)
    {
      if(missing(plot.type)) plot.type <- "confplot"
      if(inherits(x$fit,"smacofB") && !inherits(x$fit,"smacofP") && plot.type=="transplot"){
         #this is unelegant, but not sure how to call a method if it is not of the correct class   
         if(missing(main)) main <- paste("Transformation Plot")
         class(x$fit) <- c("cmdscalex",class(x))
         plot(x$fit,plot.type="transplot",main=main,...)
         class(x$fit) <- class(x$fit)[-1]
     } else {
         plot(x$fit,plot.type=plot.type,main=main,asp=asp,...)
    }
 }

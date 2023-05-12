#' S3 summary method for stops
#'
#'@param object object of class stops
#'@param ... addditional arguments
#' 
#'@export
#'@return object of class 'summary.stops'
summary.stops <- function(object,...)
    {
      sppmat <- NULL
      if(!is.null(object$fit$spp))
      { 
           spp.perc <- object$fit$spp/sum(object$fit$spp) * 100
           sppmat <- cbind(sort(object$fit$spp), sort(spp.perc))
           colnames(sppmat) <- c("SPP", "SPP(%)")
      } 
      res <- list(conf=object$fit$conf,sppmat=sppmat)
      class(res) <- "summary.stops"
      res
    }

#' S3 print method for summary.stops
#' 
#'@param x object of class summary.stops
#'@param ... additional arguments
#'@export
#'@return no return value, just prints 
print.summary.stops <- function(x,...)
    {
    cat("\n")
    cat("Configurations:\n")
    print(round(x$conf, 4))
    cat("\n\n")
    if(!is.null(x$sppmat))
     {   
      cat("Stress per point:\n")
      print(round(x$sppmat[,2], 4))
      cat("\n")
     }
    }

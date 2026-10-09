
#'@export
summary.pcops <- function(object,...)
    {
      sppmat <- NULL
      if(!is.null(object$fit$spp))
      { 
           spp.perc <- object$fit$spp/sum(object$fit$spp) * 100
           sppmat <- cbind(sort(object$fit$spp), sort(spp.perc))
           colnames(sppmat) <- c("SPP", "SPP(%)")
      } 
      res <- list(conf=object$fit$conf,sppmat=sppmat)
      class(res) <- "summary.pcops"
      res
    }

#'@export
#'@rdname print.cops 
print.summary.pcops <- function(x, digits = max(4L, getOption("digits") - 4L), ...)
    {
    cat("\n")
    cat("Configurations:\n")
    print(round(x$conf, digits = digits))
    cat("\n\n")
    if(!is.null(x$sppmat))
     {   
      cat("Stress per point (%):\n")
      print(round(x$sppmat[,2], digits = digits))
      cat("\n")
     }
    invisible(x)
    }

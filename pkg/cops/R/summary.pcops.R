
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
print.summary.pcops <- function(x,...)
    {
    cat("\n")
    cat("Configurations:\n")
    print(round(x$conf, 4))
    cat("\n\n")
    if(!is.null(x$sppmat))
     {   
      cat("Stress per point (%):\n")
      print(round(x$sppmat[,2], 4))
      cat("\n")
     }
    }

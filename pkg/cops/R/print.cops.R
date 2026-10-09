#' Print methods for 'cops' objects
#'
#' @export
#' @param x an object of class \code{"cops"} or \code{"pcops"} or their summaries.
#' @param digits integer; number of significant digits used for printing.
#' @param ... further arguments passed to or from other methods.
print.cops <- function(x, digits = max(4L, getOption("digits") - 4L), ...)
    {
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    cat("Model: COPS with parameter vector=",round(x$parameters, digits = digits),"\n")
    cat("\n")
    cat("Number of objects:", x$nobj, "\n")
    cat("Stress of configuration (default normalization):", round(x$stress, digits = digits), "\n")
    cat("OPTICS Cordillera: Raw", round(x$OC$raw, digits = digits),"Normed", round(x$OC$normed, digits = digits),"\n")
    cat("Cluster optimized loss (copstress): ", round(x$copstress, digits = digits), "\n")
    cat("Stress weight:",x$stressweight," OPTICS Cordillera weight:",x$cordweight,"\n")
    cat("Number of iterations of",x$optimethod,"optimization:", x$niter, "\n")
    cat("\n")
    invisible(x)
    }

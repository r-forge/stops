
#' @export
#' @rdname print.cops
print.copsc <- function(x, digits = max(4L, getOption("digits") - 4L), ...)
    {
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    cat("Model:",x$type,"COPS-C with parameter vector =",round(x$parameters,digits=digits),"\n")
    cat("\n")
    cat("Number of objects:", x$nobj, "\n")
    cat("Stress-1 value of configuration:", round(x$stress,digits=digits), "\n")
    cat("OPTICS Cordillera: Raw", round(x$OC$raw,digits=digits),"Normed", round(x$OC$normed,digits=digits),"\n")
    cat("Cluster optimized loss (copstress): ", round(x$copstress,digits=digits), "\n")
    cat("Stress weight:", x$stressweight," OPTICS Cordillera weight:",x$cordweight,"\n")
    cat("Number of iterations of",x$optimethod,"optimization:", x$niter, "\n")
    cat("\n")
    invisible(x)
    }

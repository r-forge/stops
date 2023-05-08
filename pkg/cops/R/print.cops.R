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

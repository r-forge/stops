
#'@export
print.copsc <- function(x,...)
    {
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    cat("Model:",x$type,"COPS-C with parameter vector =",x$parameters,"\n")
    cat("\n")
    cat("Number of objects:", x$nobj, "\n")
    cat("Stress-1 value of configuration:", round(x$stress,5), "\n")
    cat("OPTICS Cordillera: Raw", round(x$OC$raw,5),"Normed", round(x$OC$normed,5),"\n")
    cat("Cluster optimized loss (copstress): ", round(x$copstress,5), "\n")
    cat("Stress weight:", x$stressweight," OPTICS Cordillera weight:",x$cordweight,"\n")
    cat("Number of iterations of",x$optimethod,"optimization:", x$niter, "\n")
    cat("\n")
    }

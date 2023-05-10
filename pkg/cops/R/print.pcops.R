

#'@export
print.pcops <- function(x,...)
    {
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    cat("Model: P-COPS with",x$type, x$loss,"loss function and theta parameter vector =",x$parameters,"\n")
    cat("\n")
    cat("Number of objects:", x$nobj, "\n")
    cat("MDS loss value:", x$stress, "\n")
    cat("OPTICS Cordillera: Raw", x$OC$raw,"Normed", x$OC$normed,"\n")
    cat("Cluster optimized loss (copstress): ", x$copstress, "\n")
    cat("MDS loss weight:",x$stressweight," OPTICS Cordillera weight:",x$cordweight,"\n")
    cat("Number of iterations of",x$optimethod,"optimization:", x$optim$counts["function"], "\n")
    cat("\n")
    }

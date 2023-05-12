#' S3 print method for stops objects
#'
#'@param x stops object
#'@param ... additional arguments
#'@export
#'@return no return value, just prints
print.stops <- function(x,...)
    {
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    cat("Model:",x$stoptype ,"STOPS with", x$type, x$loss,"loss function and theta parameter vector",paste("(",paste(attr(x$parameters,"names"),collapse=" "),")",sep="")," = ",x$parameters,"\n")
    cat("\n")
    cat("Number of objects:", x$nobj, "\n")
    cat("MDS loss value:", x$stress.m, "\n")
    cat("C-Structuredness Indices:", t(data.frame(names(x$strucindices),x$strucindices)),"\n")
    cat("Structure optimized loss (stoploss):", x$stoploss, "\n")
    cat("MDS loss weight:",x$stressweight,"c-structuredness weights:",x$strucweight,"\n")
    cat("Number of iterations of",x$optimethod,"optimization:", x$optim$counts["function"], "\n")
    cat("\n")
    }

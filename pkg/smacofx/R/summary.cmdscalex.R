#'@export
summary.cmdscalex <- function(object,...)
    {
        out <- object$points
        class(out) <- "summary.cmdscalex"
        out
    }


#'@export
print.summary.cmdscalex <- function(x,...)
{
    if(missing(digits)) digits <- 4
    cat("\n")
    cat("Configurations:\n")
    print(round(x, digits = digits))
    cat("\n")
    invisible(x)
    }

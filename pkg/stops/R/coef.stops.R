

#' S3 coef method for stops objects
#'
#'@param object object of class stops 
#'@param ... addditional arguments 
#'@export
#'@importFrom stats coef
#'@return a vector of hyperparmeters theta 
coef.stops <- function(object,...)
    {
    return(c(object$parameters))
    }

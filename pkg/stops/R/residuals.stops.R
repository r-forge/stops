#' S3 residuals method for stops
#'@param object object of class stops
#'@param ... addditional arguments
#'@importFrom stats residuals
#'@export
#'@return a vector of residuals (observed minus fitted distances) 
residuals.stops <- function(object,...)
    {
    stats::residuals(object$fit,...)
    }

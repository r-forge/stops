#'@importFrom stats coef
#'@export 
coef.cops <- function(object,...)
{
  return(object$parameters)
 }


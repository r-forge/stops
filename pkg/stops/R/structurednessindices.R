#'c-linearity
#'calculates c-linearity as the multiple correlation
#'
#' @param confs a numeric matrix or data frame
#' @param ... additional arguments to be passed to lm.fit
#'
#'
#' @importFrom stats lm summary.lm
#' 
#' @examples
#' x<-1:10
#' y<-2+3*x+rnorm(10)
#' confs<-cbind(x,y)
#' c_linearity(confs)
#' @export
c_linearity <- function(confs,...)
    {
        y <- confs[,1]
        n <- dim(confs)[1]
        p <- dim(confs)[2]
        x <- confs[,2:p]
        tmp <- stats::lm(y~x,...)
        out <- sqrt(summary(tmp)$r.squared)
        out
    }
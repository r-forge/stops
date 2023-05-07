#' Auxfunction1
#' 
#' only used internally
#' @param x matrix
mkBmat <- function (x) {
    d <- rowSums (x)
    x <- -x
    diag (x) <- d
    return (x)
}

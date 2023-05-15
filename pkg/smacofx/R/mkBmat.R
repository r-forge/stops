#' Auxfunction1
#' 
#' only used internally
#' @param x matrix
mkBmat <- function (x) {
    d <- rowSums (x)
    x <- -x
    diag (x) <- d
    if(any(is.na(x))) stop("An error from deep within the bowels: NAs in the MY matrix. This is likely because of too many 0 dissimilarities after transformation and it might help to remove offending observations.")
    if(any(!is.finite(x))) x[!is.finite(x)] <- sign(x[!is.finite(x)])*(max(x[is.finite(x)]))#New as this was a bug
    return (x)
}

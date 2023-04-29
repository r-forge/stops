#' Squared distances
#'
#' @param x numeric matrix
#' @return squared distance matrix
sqdist <- function (x) {
    s <- tcrossprod (x)
    v <- diag (s)
    return (outer (v, v, "+") - 2 * s)
}

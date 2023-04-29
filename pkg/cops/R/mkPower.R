
#' Take matrix to a power 
#'
#' @param x matrix
#' @param r numeric (power)
#' @return a matrix
mkPower<-function(x,r) {
    n<-nrow(x)
    tmp <- abs((x+diag(n))^r)-diag(n)
    return(tmp)
}

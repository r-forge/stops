#' Helper function to conduct jackknife MDS
#'
#' Function deletes every object row and columns once and fits the MDS in object and returns the configuration. The deleted row is set to 0 in the configuration. Is meant for smacofx functions, but should also work for every smacof models. 
#' 
#' @param object  Object of class smacofP if used as method or another object inheriting from smacofB. Note we assume the MDS model was fitted on a symmetric matrix/data frame or dist object 
#' @param delta the data (symmetric matrix, data frame or dist object)
#' @param ndim target dimension of the mds
#' @param type type of MDS
#'
#' @return An array of size n with n coonfigurations
smacofxDeleteOne <- function (object, delta, ndim, type) {
    ocall<-as.call(object$call) 
    n <- nrow (delta)
    x <- array (0, c (n, ndim, n))
    for (i in 1:n) {
        #xi <- smacofSym(delta[-i, -i], ndim = ndim, type = type)$conf
        ocall$delta<-delta[-i, -i]
        xi <- eval(ocall)
        xconf <- xi$conf
        x[((1 : n)[-i]), (1 : ndim), i] <- xconf
        x[i, (1 : ndim), i] <- 0
        }
    return (x)
}

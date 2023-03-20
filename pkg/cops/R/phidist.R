#' Calculating the pairwise phi distance matrix between n vectors
#'
#' The pairwise phi-distance of two vectors x and y is sqrt(sum(((x[i]-y[i])^2)/((x[i]+y[i])*(sum(x)+sum(y))))). The function calculates this for all pairs of rows of a matrix or data frame X.   
#'
#' @param X an n times p numeric matrix or data frame 
#'
#' @return a symmetric n times n matrix of pairwise phi distance (between rows of X) with 0 in the main diagonal. Is an object of class distance and matrix. 
#' 
#' @importFrom analogue distance
#' @export
phidistance <- function(X)
{
      summ<- apply(X,1,sum)
      NN <- outer(summ,summ,"+")
      disttemp <- analogue::distance(X,method="SQeuclidean")
      out <- disttemp/NN
      out <- sqrt(out)
      attr(out,"method") <- 'phi'
      out
}

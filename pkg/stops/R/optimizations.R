#' (Adaptive) Version of Luus-Jaakola Optimization
#'
#' Adaptive means that the search space reduction factors in the number of iterations; makes convergence faster at about 100 iterations 
#' 
#' @param x starting values
#' @param fun function to optimize
#' @param red value of the reduction of the search region
#' @param ... additional arguments to be passed to the function to be optimized
#' @param lower The lower contraints of the search region
#' @param upper The upper contraints of the search region 
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param acc if the numerical accuracy of two successive target function values is below this, stop the optimization; defaults to 1e-6
#' @param accd if the width of the search space is below this, stop the optimization; defaults to 1e-4  
#' @param itmax maximum number of iterations
#' @param adaptive should the adaptive version be used? defaults to TRUE.
#' 
#' @return A list with the components (\code{\link{optim}})
#' \itemize{
#'      \item par The position of the optimimum in the search space (parameters that minimize the function; argmin fun)
#'      \item value The value of the objective function at the optimum (min fun)
#'      \item counts The number of iterations performed at convergence with entries fnction for the number of iterations and gradient which is always NA at the moment
#'      \item convergence 0 successful completion by the accd or acc criterion, 1 indicate iteration limit was reached, 99 is a problem 
#'      \item message is NULL (only for compatibility or future use)
#' }
#'
#' @importFrom stats runif
#' 
#' @examples
#' fbana <- function(x) {
#' x1 <- x[1]
#' x2 <- x[2]
#' 100 * (x2 - x1 * x1)^2 + (1 - x1)^2
#' }
#' res1<-ljoptim(c(-1.2,1),fbana,lower=-5,upper=5,accd=1e-16,acc=1e-16)
#' res1
#'
#' set.seed(210485)
#' fwild <- function (x) 10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80
#' plot(fwild, -50, 50, n = 1000, main = "ljoptim() minimising 'wild function'")
#' res2<-ljoptim(50, fwild,lower=-50,upper=50,adaptive=FALSE,accd=1e-16,acc=1e-16)
#' points(res2$par,res2$value,col="red",pch=19)
#' res2
#' 
#' @export
ljoptim <- function(x,fun,...,red=ifelse(adaptive,0.99,0.95),lower,upper,acc=1e-6,accd=1e-4,itmax=1000,verbose=0,adaptive=TRUE) {
    #addargs <- ...
    x[is.na(x)] <- stats::runif(sum(is.na(x)),min=min(lower),max=max(upper))
   # cat("------ x:",x,"\n")
       #http://en.wikipedia.org/wiki/Luus%E2%80%93Jaakola
    n <- length(x)
    d1 <- upper-lower
    x <- ifelse(x<lower,lower+stats::runif(length(x))*d1,ifelse(x>upper,upper-stats::runif(length(x))*d1,x)) 
    if(verbose>2) cat("------ x:",x,"\n")
    itel <- 1
    fold <- do.call(fun,list(x,...))
   # cat("fold",fold,"\n")
    convo <- 99L
    d <- d1
    repeat {
       addi <- sapply(d,function(f) runif(1,min=-f,max=f)) #use triangular or trapezoidal distribution?
       #addi <- rtriangle(n,a=-d,b=d) #use triangular or trapezoidal distribution? 
       #a <- rcauchy(n,scale=d/8)
       y <- x+addi
     #  cat("------ y:",y,"\n")
       y <- ifelse(y<lower,lower+stats::runif(length(y))*d1,ifelse(y>upper,upper-stats::runif(length(y))*d1,y))
       if(verbose>2) cat("------ y:",y,"\n")
       fnew <- do.call(fun,list(y,...))
      # cat("fnew",fnew,"\n")
       if (verbose>0) {
          cat (
          formatC (itel, digits =0, width = 8, format = "d"),
          formatC (x, digits = 8, width = 10, format = "f"),
          formatC (fold, digits = 10, width = 13, format = "f"),    
          formatC (y, digits = 8, width = 10, format = "f"),
          formatC (fnew, digits = 10, width = 13, format = "f"),
          formatC (d, digits = 10, width = 13, format = "f"), "\n")
      }  
     if ((itel == itmax) || (abs(fold - fnew) < acc) || (max(d) < accd)) {
         convo <- 0L
         if(itel==itmax) convo <- 1L
         break ()
     }
     itel <- itel + 1
     if(fnew<fold) {
               x <- y
               fold <- fnew
             #  d <- max(upper-lower)
           }
           else {
                 redi <- red
                 if(adaptive) {
                        allits <- min(floor((log(accd)-log(max(upper-lower)))/log(red)),itmax)
                        redi <- red*((allits+1-itel)/allits)  #adaptive contraction
                    }
                 d <- redi*d
                }
   }
  return(list(par=x,value=fold,counts=c(`function`=itel,gradient=NA),convergence=convo,message=NULL))
}

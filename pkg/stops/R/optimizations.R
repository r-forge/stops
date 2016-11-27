#' (Adaptive) Version of Luus-Jaakola Optimization
#'
#' Adaptive means that the search space reduction factors in the number of iterations; makes convergence faster at about 100 iterations 
#' 
#' @param x optional starting values
#' @param fun function to minimize
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


#' Bayesian Optimization by a (treed) Bayesian Gaussian Process Prior (with jumps to linear models) surrogate model
#' Essentially a wrapper for the functionality in tgp that has the same slots as optim with defaults for STOPS models.
#'
#' 
#' @param x optional starting values
#' @param fun function to minimize
#' @param ... additional arguments to be passed to the function to be optimized
#' @param initpoints the number of points to sample initially to fit the surrogate model
#' @param lower The lower contraints of the search region
#' @param upper The upper contraints of the search region 
#' @param verbose numeric value hat prints information on the fitting process; >2 is extremely verbose
#' @param acc if the numerical accuracy of two successive target function values is below this, stop the optimization; defaults to 1e-8
#' @param itmax maximum number of iterations
#' @param model which surrogate model class to use (currently uses defaults only, will extend this to tweak the model)
#' 
#' @return A list with the components (for compatiility with \code{\link{optim}})
#' \itemize{
#'      \item par The position of the optimimum in the search space (parameters that minimize the function; argmin fun)
#'      \item value The value of the objective function at the optimum (min fun)
#'      \item counts The number of iterations performed at convergence with entries fnction for the number of iterations and gradient which is always NA at the moment
#'      \item convergence 0 successful completion by the accd or acc criterion, 1 indicate iteration limit was reached, 99 is a problem 
#'      \item message is NULL (only for compatibility or future use)
#'      \item history the improvement history
#'      \item tgpout the output of the tgp model    
#' }
#'
#' @importFrom tgp lhs dopt.gp optim.step.tgp
#' 
#' @examples
#' \donttest{
#' fbana <- function(x) {
#' x1 <- x[1]
#' x2 <- x[2]
#' 100 * (x2 - x1 * x1)^2 + (1 - x1)^2
#' }
#' #Of course not clever as the function is simple enough to evaluate to use other optim 
#' res1<-tgpoptim(c(-1.2,1),fbana,lower=c(-5,-5),upper=c(5,5),acc=1e-16,itmax=100)
#' res1
#'
#' set.seed(210485)
#' fwild <- function (x) 10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80
#' plot(fwild, -50, 50, n = 1000, main = "Bayesian GP Optimization minimizing 'wild function'")
#' res2<-tgpoptim(50, fwild,lower=-50,upper=50,adaptive=FALSE,acc=1e-16,itmax=100,model=bgp)
#' points(res2$par,res2$value,col="red",pch=19)
#' res2
#' }
#' @export
tgpoptim <- function(x,fun,...,initpoints=10,lower,upper,acc=1e-8,itmax=10,verbose=0,model="bgp") {
    #TODO: add control for the tgp models...
        #if(!isNamespaceLoaded("tgp")) attachNamespace("tgp")
        rect <- cbind(lower,upper)
        Xcand <- tgp::lhs(initpoints*100,rect)
        if(missing(x)) x <- Xcand[1,]
        x <- t(x)
        X <- tgp::dopt.gp(initpoints-1,X=x,Xcand)$XX
        X <- rbind(x,X)
        if(verbose>0) cat("Evaluating function at initial points","\n")
        Z <- apply(X, 1, fun)
        out <- progress <-  NULL
        convo <- 99L
        itel <- 1
        fold <- min(Z)
        model <- get(model,envir=getNamespace("tgp"))
        if(verbose>0) cat("Starting Bayesian Optimization","\n")
        repeat {
          ## get recommendations for the next point to sample
          out <- tgp::optim.step.tgp(fun, X=X, Z=Z, rect=rect, model=model, prev=out,verb=verbose-3) 
          ## add in the inputs, and newly sampled outputs
          X <- rbind(X, out$X)
          ztmp <- apply(out$X,1,fun)
          Z <- c(Z, ztmp)
          ## keep track of progress and optimum
          fnew <- out$progress$z 
          progress <- rbind(progress, out$progress)
          if(verbose>1) print(progress[itel,])
          if ((itel == itmax) || (abs(fold - fnew) < acc)) {
            convo <- 0L
            if(itel==itmax) convo <- 1L
            break ()
          }
          itel <- itel + 1
          fold <- fnew
        }
        rets <- progress[which.min(progress$z),1:length(x)]
        return(list(par=as.numeric(rets),value=min(progress$z),counts=c(`function`=itel,gradient=NA),convergence=convo,message=NULL,history=progress,tgpout=out))
}


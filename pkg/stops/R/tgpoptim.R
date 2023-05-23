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
#'      \item par The position of the optimum in the search space (parameters that minimize the function; argmin fun). 
#'      \item value The value of the objective function at the optimum (min fun). Note we do not use the last value in the candidate list but the best candidate (which can but need not coincide). 
#'      \item svalue The value of the surrogate objective function at the optimal parameters
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
#' res1<-tgpoptim(c(-1.2,1),fbana,lower=c(-5,-5),upper=c(5,5),acc=1e-16,itmax=20)
#' res1
#'
#' fwild <- function (x) 10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80
#' plot(fwild, -50, 50, n = 1000, main = "Bayesian GP Optimization minimizing 'wild function'")
#' set.seed(210485)
#' res2<-tgpoptim(50, fwild,lower=-50,upper=50,acc=1e-16,itmax=20,model="btgpllm")
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
        Z <- apply(X, 1, fun) #Capital Z ist real function value, small z is surrogate function value
        out <- progress <-  NULL
        convo <- 99L
        itel <- 1
        fold <- sfold <- min(Z)
        model <- get(model,envir=getNamespace("tgp"))
        if(verbose>0) cat("Starting Bayesian Optimization","\n")
        repeat {
          ## get recommendations for the next point to sample
          out <- tgp::optim.step.tgp(fun, X=X, Z=Z, rect=rect, model=model, prev=out,verb=verbose-3) 
          ## add in the inputs, and newly sampled outputs
          X <- rbind(X, out$X)
          Ztmp <- apply(out$X,1,fun)
          Z <- c(Z, Ztmp)
          ## keep track of progress and optimum
          sfnew <- out$progress$z #candidate point for optimum from surrogate function 
          fnew <- apply(out$progress[1,seq.int(1,dim(X)[2]),drop=FALSE],1,fun) #real function value at surrogate optimum
          out$progress$Z <- fnew 
          progress <- rbind(progress, out$progress) #progress of optimum from surrogate #TODO: extend this by real value?
#          progress <- rbind(progress,out$progress)
          if(verbose>1) print(progress[itel,])
          if ((itel == itmax) || (abs(sfold - sfnew) < acc)) { #stop if itmax is reached or if two successive surrogate function values are essentially equal 
            convo <- 0L
            if(itel==itmax) convo <- 1L
            break ()
          }
          itel <- itel + 1
          sfold <- sfnew
        }
        rets <- progress[which.min(progress$Z),1:length(x)] #we use the best candidate of all, needs not be the last nor necessarily the same as the best surrogate value 
        #rets <- progress[itel,1:length(x)] #last surrogate function candidate
        valsorig <- fnew #real function value at best candidate
       # valsorig <- NULL
        return(list(par=as.numeric(rets),value=min(progress$Z),svalue=progress$z[which.min(progress$Z)],counts=c(`function`=itel,gradient=NA),convergence=convo,message=NULL,history=progress,tgpout=out))
}


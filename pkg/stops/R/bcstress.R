#' An MDS version for minimizing BoxCox Stress (Chen & Buja 2013)
#'
#' Based on code by Lisha Chen.
#' 
#' @param delta dissimilarity or distance matrix
#' @param init initial configuration. If NULL a classical scaling solution is used. 
#' @param ndim the dimension of the configuration
#' @param lambda lambda parameter 
#' @param mu mu parameter
#' @param nu the nu parameter
#' @param itmax number of optimizing iterations, defaults to 10000.
#' @param verbose prints progress if > 3. 
#'
#' @examples
#' dis<-smacof::kinshipdelta
#' res< bcStressMin(as.matrix(dis),mu=2,lambda=1.5,nu=0)
#' res
#' summary(res)
#' plot(res)
#' 
#' @export
bcStressMin <- function(delta,init=NULL,verbose=0,ndim=2,lambda=1,mu=1,nu=0,itmax=10000)
 #Not the local MDS version, there is no neighbour concept here.  
{
  Do <- delta 
  d <- ndim
  X1 <- init
  niter <- itmax
  lambdaorig <- lambda
  lambda <- 1/lambda
  n <- nrow(Do)
  Dnu <- Do^nu
  Dnulam <- Do^(nu+1/lambda)
  diag(Dnu) <- 0
  diag(Dnulam) <- 0
  t <- 0 #not needed
  Grad <- matrix (0, nrow=n, ncol= d)
   if(is.null(X1))
    {
      cmd <- cmds(Do)  
      X1 <- cmd$vec[,1:d]%*%diag(cmd$val[1:d])+enorm(Do)/n/n*0.01*matrix(rnorm(n*d),nrow=n,ncol=d)
    }
  D1 <- as.matrix(dist(X1))
  
  X1 <- X1*enorm(Do)/enorm(D1)

  s1 <- Inf #stressinit
  s0 <- 2 #stress check for update
  stepsize <-0.1
  i <- 0

while ( stepsize > 1E-5 && i < niter)
  {
    if (s1 >= s0 && i>1)
    {
       #stepsize if stress (s1) >= old s1 (=s0) oder >2   
       stepsize<- 0.5*stepsize
       X1 <- X0 - stepsize*normgrad
     }
    else 
      {
      #stepsize if stress (s1) < old s1       
      stepsize <- 1.05*stepsize
      X0 <- X1
      D1mu2 <- D1^(mu-2)
      diag(D1mu2) <- 0
      D1mulam2 <- D1^(mu+1/lambda-2)
      diag(D1mulam2) <- 0
      M <- Dnu*D1mulam2-D1mu2*Dnulam    
      E <- matrix(rep(1,n*d),n,d)
      Grad <- X0*(M%*%E)-M%*%X0
      normgrad <- (enorm(X0)/enorm(Grad))*Grad         
      X1 <- X0 - stepsize*normgrad
     }
    i <- i+1
    s0 <- s1 
    D1 <- as.matrix(dist(X1))
    D1mulam <- D1^(mu+1/lambda)
    Domulam <- Do^(mu+1/lambda) #new
    diag(D1mulam) <- 0
    D1mu <- D1^mu
    Domu <- Do^mu #new
    diag(D1mu) <- 0
    diag(Domu) <- 0 #new
    diag(D1) <- 0
    if(mu+1/lambda==0)
      {
       diag(D1)<-1
       diag(Do)<-1 #new
       #stress   
       s1 <- sum(Dnu*log(D1))-sum((D1mu-1)*Dnulam)/mu 
       #NEW: added normalization. Looks legit. 
       normo <- sum(Dnu*log(Do))-sum((Domu-1)*Dnulam)/mu 
       s1n <- 1-s1/normo
      }

    if(mu==0)
      {
      diag(D1)<-1
      diag(Do)<-1 #new
      #stress
      s1 <- sum(Dnu*(D1mulam-1))/(mu+1/lambda) -sum(log(D1)*Dnulam)
      #NEW: norming
      normo <- sum(Dnu*(Domulam-1))/(mu+1/lambda) -sum(log(Do)*Dnulam)
      s1n <- 1-s1/normo
      }
    if(mu!=0&(mu+1/lambda)!=0)
    {
        s1 <- sum(Dnu*(D1mulam-1))/(mu+1/lambda)-sum((D1mu-1)*Dnulam)/mu
        #NEW: norming:
       normo <- sum(Dnu*(Domulam-1))/(mu+1/lambda)-sum((Domu-1)*Dnulam)/mu #-t*sum((D1mu-1)*(1-Inb1))/mu
       s1n <- 1-s1/normo
    }
    ## Printing and Plotting
    if(verbose > 3 && (i+1)%%100/verbose==0)
      {
        print (paste("niter=",i+1," stress=",round(s1,5)," stressn=",round(sqrt(s1n),5), sep=""))
      }

  }
  result <- list()
  result$conf <- X1 #new
  result$confdist <- D1
  result$delta <- Do
  result$mu <- mu
  result$lambda <- lambda
  result$nu <- nu
  result$theta <- c(mu,lambda,nu)
  result$stress.r <- s1
  result$stress.m <- s1n
  result$stress <- sqrt(s1n)
  return(result)
}

#' normalization function
#norm <- function(x) sqrt(sum(x^2))

#' Classical Scaling
#' @param Do dissimilarity matrix
cmds <- function(Do)
  {
    n <- nrow(Do)
    J <- diag(rep(1,n)) - 1/n * rep(1,n) %*% t(rep(1,n))
    B <- - 1/2 * J %*% (Do^2) %*% J
    pc <- eigen(B)
    return(pc)
  }

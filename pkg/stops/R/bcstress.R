#' An MDS version for minimizing BoxCox Stress (Chen & Buja 2013)
#'
#' 
#' @param delta dissimilarity or distance matrix
#' @param init initial configuration. If NULL a classical scaling solution is used. 
#' @param ndim the dimension of the configuration
#' @param mu mu parameter. Should be 0 or larger for everything working ok. If mu<0 it works but the model is strange and use normalized stress at your own risk.
#' @param lambda lambda parameter. Should be equal or larger than -mu for everything working ok. If lambda < -mu it works but the model is strange and use normalized stress at your own risk.
#' @param nu the nu parameter.
#' @param itmax number of optimizing iterations, defaults to 10000.
#' @param verbose prints progress if > 3.
#' @param addD0 a small number that's added for D(X)=0 for numerical evaluation of worst fit (numerical reasons, see details). If addD0=0 the normalized stress for mu!=0 and mu+lambda!=0 is correct, but will give useless normalized stress for mu=0 or mu+lambda!=0. 
#'
#' @details For numerical reasons with certain parameter combinations, the normalized stress uses a configuration as worst result where every d(X) is 0+addD0. The same number is not added to the delta so there is a small inaccuracy of the normalized stress (but negligible of min(delta)>>addD0). Also, for mu<0 or mu+lambda<0 the normalization can not generally be trusted (in the worst case of D(X)=0 one would have an 0^(-a)).    
#' 
#' 
#' @author Lisha Chen & Thomas Rusch
#' 
#' @examples
#' dis<-smacof::kinshipdelta
#' res<-bcStressMin(as.matrix(dis),mu=2,lambda=1.5,nu=0)
#' res
#' summary(res)
#' plot(res)
#' 
#' @export
bcStressMin <- function(delta,init=NULL,verbose=0,ndim=2,mu=1,lambda=1,nu=0,itmax=10000,addD0=1e-4)
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
      X1 <- cmd$vec[,1:d]%*%diag(cmd$val[1:d])+enorm(Do)/n/n*0.01*matrix(stats::rnorm(n*d),nrow=n,ncol=d)
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
    diag(D1mulam) <- 0
    D1mu <- D1^mu
    diag(D1mu) <- 0
    if(mu+1/lambda==0)
      {
       diag(D1)<-1
       s1 <- sum(Dnu*log(D1))-sum((D1mu-1)*Dnulam)/mu 
      }

    if(mu==0)
      {
      diag(D1)<-1
      s1 <- sum(Dnu*(D1mulam-1))/(mu+1/lambda) -sum(log(D1)*Dnulam)
      }
    if(mu!=0&(mu+1/lambda)!=0)
    {
        s1 <- sum(Dnu*(D1mulam-1))/(mu+1/lambda)-sum((D1mu-1)*Dnulam)/mu     
    }
    ## Printing and Plotting
    if(verbose > 3 & (i+1)%%100/verbose==0)
      {
        print (paste("niter=",i+1," stress=",round(s1,5), sep=""))
      }

  }
  #New For normalization of stress; if it doesn't work normalize the X1 to be smaller than Do
  X1 <- X1*sum(Do*D1)/sum(D1^2)
  D1 <- as.matrix(dist(X1)) 
  D0 <- D1*0+addD0 #for numerical reasons for mu=0 and mu+lambda=0 all get an extra p for "plus"
  Dop <- Do# we could also add +1e-4 here for reasons of comparability with the D0
  diag(Dop) <- 0
  diag(D0) <- 0
  Dopmulam <- Dop^(mu+1/lambda) #new
  D0mulam <- D0^(mu+1/lambda) #new
  diag(Dopmulam) <- 0
  diag(D0mulam) <- 0
  Dopmu <- Dop^mu #new
  D0mu <- D0^mu #new
  diag(Dopmu) <- 0 #new
  diag(D0mu) <- 0 #new
  Dpnu <-  Dop^nu
  Dpnulam <- Dop^(nu+1/lambda)
  diag(Dpnu) <- 0
  diag(Dpnulam) <- 0
  if(mu+1/lambda==0) {
       diag(D0) <- 1
       diag(Do)<-1 #new
       norm0 <- sum(Dnu*log(D0))-sum((D0mu-1)*Dnulam)/mu
       normo <- sum(Dnu*log(Do))-sum((Dopmu-1)*Dnulam)/mu
       s1n <- (s1-normo)/(norm0-normo)       
       }
  if(mu==0) {
       diag(D0) <- 1
       diag(Do)<-1 #new
       norm0 <- sum(Dnu*(D0mulam-1))/(mu+1/lambda) - sum(log(D0)*Dnulam)
       normo <- sum(Dnu*(Dopmulam-1))/(mu+1/lambda) - sum(log(Do)*Dnulam)
      s1n <- (s1-normo)/(norm0-normo)
      }
  if(mu!=0&(mu+1/lambda)!=0)
  {
      #D0 <- D1*0
      #D0mu <- D0^mu # infinity
      #diag(D0mu) <- 0 
      #D0mulam <- D0^(mu+1/lambda) 
      #diag(D0mulam) <- 0
      #Domulam <- Do^(mu+1/lambda) 
      #diag(Domulam) <- 0
      #Domu <- Do^mu #new
      #diag(Domu) <- 0 #new
      norm0 <- sum(Dpnu*(D0mulam-1))/(mu+1/lambda)-sum((D0mu-1)*Dpnulam)/mu
      normo <- sum(Dpnu*(Dopmulam-1))/(mu+1/lambda)-sum((Dopmu-1)*Dpnulam)/mu
      #norm0 <- sum(Dnu*(D0mulam-1))/(mu+1/lambda)-sum((D0mu-1)*Dnulam)/mu
      #normo <- sum(Dnu*(Domulam-1))/(mu+1/lambda)-sum((Domu-1)*Dnulam)/mu
      s1n <- (s1-normo)/(norm0-normo)
      # s1n <- 1-s1/normo #normalized stress
  }
  result <- list()
  result$conf <- X1 #new
  result$confdist <- stats::as.dist(D1)
  result$delta <- stats::as.dist(Do)
  result$obsdiss <- stats::as.dist(Do)
  result$mu <- mu
  result$lambda <- lambda
  result$nu <- nu
  result$pars <- c(mu,lambda,nu)
  result$model<- "Box-Cox Stress MDS"
  result$call <- match.call()
  result$ndim <- ndim
  result$nobj <- n
  result$niter <- i
  result$theta <- c(mu,lambda,nu)
  result$stress.r <- s1
  result$stress.m <- s1n
  result$stress <- sqrt(s1n)
  result$type <- "Box-Cox Stress"
  class(result) <- c("smacofP","smacofB","smacof")
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



## function (delta, init = NULL, verbose = 0, ndim = 2, lambda = 1, 
##     mu = 1, nu = 0, itmax = 10000) 
## {
##     Do <- delta
##     d <- ndim
##     X1 <- init
##     niter <- itmax
##     lambdaorig <- lambda
##     lambda <- 1/lambda
##     n <- nrow(Do)
##     Dnu <- Do^nu
##     Dnulam <- Do^(nu + 1/lambda)
##     diag(Dnu) <- 0
##     diag(Dnulam) <- 0
##     t <- 0
##     Grad <- matrix(0, nrow = n, ncol = d)
##     if (is.null(X1)) {
##         cmd <- cmds(Do)
##         X1 <- cmd$vec[, 1:d] %*% diag(cmd$val[1:d]) + enorm(Do)/n/n * 
##             0.01 * matrix(stats::rnorm(n * d), nrow = n, ncol = d)
##     }
##     D1 <- as.matrix(dist(X1))
##     X1 <- X1 * enorm(Do)/enorm(D1)
##     s1 <- Inf
##     s0 <- 2
##     stepsize <- 0.1
##     i <- 0
##     while (stepsize > 1e-05 && i < niter) {
##         if (s1 >= s0 && i > 1) {
##             stepsize <- 0.5 * stepsize
##             X1 <- X0 - stepsize * normgrad
##         }
##         else {
##             stepsize <- 1.05 * stepsize
##             X0 <- X1
##             D1mu2 <- D1^(mu - 2)
##             diag(D1mu2) <- 0
##             D1mulam2 <- D1^(mu + 1/lambda - 2)
##             diag(D1mulam2) <- 0
##             M <- Dnu * D1mulam2 - D1mu2 * Dnulam
##             E <- matrix(rep(1, n * d), n, d)
##             Grad <- X0 * (M %*% E) - M %*% X0
##             normgrad <- (enorm(X0)/enorm(Grad)) * Grad
##             X1 <- X0 - stepsize * normgrad
##         }
##         i <- i + 1
##         s0 <- s1
##         D1 <- as.matrix(dist(X1))
##         D1mulam <- D1^(mu + 1/lambda)
##         diag(D1mulam) <- 0
##         D1mu <- D1^mu
##         diag(D1mu) <- 0
##         if (mu + 1/lambda == 0) {
##             diag(D1) <- 1
##             s1 <- sum(Dnu * log(D1)) - sum((D1mu - 1) * Dnulam)/mu
##         }
##         if (mu == 0) {
##             diag(D1) <- 1
##             s1 <- sum(Dnu * (D1mulam - 1))/(mu + 1/lambda) - 
##                 sum(log(D1) * Dnulam)
##         }
##         if (mu != 0 & (mu + 1/lambda) != 0) {
##             s1 <- sum(Dnu * (D1mulam - 1))/(mu + 1/lambda) - 
##                 sum((D1mu - 1) * Dnulam)/mu
##         }
##         if (verbose > 3 & (i + 1)%%100/verbose == 0) {
##             print(paste("niter=", i + 1, " stress=", round(s1, 
##                 5), sep = ""))
##         }
##     }
##    # X1 <- X1 * sum(Do * D1)/sum(D1^2)
##     D1 <- as.matrix(dist(X1))
##     D0 <- D1 * 0 + 1e-04
##     Dop <- Do + 1e-04
##     Dopmulam <- Dop^(mu + 1/lambda)
##     Dopmu <- Dop^mu
##     diag(Dopmu) <- 0
##     diag(Dop) <- 1
##     diag(D0) <- 1
##     D0mu <- D0^mu
##     Dopmu <- Dop^mu
##     diag(Dopmu) <- 0
##     Dpnu <- Dop^nu
##     Dpnulam <- Dop^(nu + 1/lambda)
##     if (mu + 1/lambda == 0) {
##         norm0 <- sum(Dpnu * log(D0)) - sum((D0mu - 1) * Dpnulam)/mu
##         normo <- sum(Dpnu * log(Dop)) - sum((Dopmu - 1) * Dpnulam)/mu
##         s1n <- (s1 - normo)/(norm0 - normo)
##     }
##     if (mu == 0) {
##         norm0 <- sum(Dpnu * (D0mulam - 1))/(mu + 1/lambda) - 
##             sum(log(D0) * Dpnulam)
##         normo <- sum(Dpnu * (Dopmulam - 1))/(mu + 1/lambda) - 
##             sum(log(Dop) * Dpnulam)
##         s1n <- (s1 - normo)/(norm0 - normo)
##     }
##     if (mu != 0 & (mu + 1/lambda) != 0) {
##         norm0 <- sum(Dpnu * (D0mulam - 1))/(mu + 1/lambda) - 
##             sum((D0mu - 1) * Dpnulam)/mu
##         normo <- sum(Dpnu * (Dopmulam - 1))/(mu + 1/lambda) - 
##             sum((Dopmu - 1) * Dpnulam)/mu
##         s1n <- (s1 - normo)/(norm0 - normo)
##     }
##     result <- list()
##     result$conf <- X1
##     result$confdist <- stats::as.dist(D1)
##     result$delta <- stats::as.dist(Do)
##     result$obsdiss <- stats::as.dist(Do)
##     result$mu <- mu
##     result$lambda <- lambda
##     result$nu <- nu
##     result$pars <- c(mu, lambda, nu)
##     result$model <- "Box-Cox Stress MDS"
##     result$call <- match.call()
##     result$ndim <- ndim
##     result$nobj <- n
##     result$niter <- i
##     result$theta <- c(mu, lambda, nu)
##     result$stress.r <- s1
##     result$stress.m <- s1n
##     result$stress <- sqrt(s1n)
##     result$type <- "Box-Cox Stress"
##     class(result) <- c("smacofP", "smacofB", "smacof")
##     return(result)
## }

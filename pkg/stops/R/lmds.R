#' An MDS version for local MDS (Chen & Buja 2006)
#'
#' Based on code by Lisha Chen.
#' 
#' @param delta dissimilarity or distance matrix
#' @param init initial configuration. If NULL a classical scaling solution is used. 
#' @param ndim the dimension of the configuration
#' @param k the k neighbourhood parameter parameter 
#' @param tau the penalty parameter 
#' @param itmax number of optimizing iterations, defaults to 10000.
#' @param verbose prints progress if > 3. 
#'
#' @examples
#' dis<-smacof::kinshipdelta
#' res< lmds(as.matrix(dis),k=2,tau=0.1)
#' res
#' summary(res)
#' plot(res)
#' 
#' @export
lmds <- function(delta,init=NULL,ndim=3,k=2,tau=1,
                   itmax=1000,verbose=0)
{
    Do <- delta
    X1 <- init
    lambda <- 1
    mu <- 1
    nu <- 0
    niter <- itmax
    d <- ndim
    n <- nrow(Do)

    #New: make the neighbourhood graph adjacency matrix Inb
    Daux <- apply(Do,2,sort)[k+1,]
    Inb <- ifelse(Do>Daux, 0, 1)
    
   #New I don't think we need this 
   #k.v <-  apply(Inb, 1, sum)
   #k <- (sum(k.v)-n)/n
   #Inb.sum <- matrix(rep(k.v, n),ncol=n)
   #Mka <- 0
   Inb1 <- pmax(Inb,t(Inb)) # expanding the neighbors for symmetry.

   Dnu <- ifelse(Inb1==1, Do^nu, 0)
   Dnulam <- ifelse(Inb1==1, Do^(nu+1/lambda), 0)
   diag(Dnu) <- 0
   diag(Dnulam) <- 0

  cc <- (sum(Inb1)-n)/n/n*median(Dnulam[Dnulam!=0])
  t <- tau*cc
  
  Grad <- matrix (0, nrow=n, ncol= d)
  #if(is.null(X1) & random.start==1) X1 <- matrix(rnorm(n*d),nrow=n,ncol=d)
  if(is.null(X1))
    {
      cmd <- cmds(Do)
      X1 <- cmd$vec[,1:d]%*%diag(cmd$val[1:d])+
        norm(Do)/n/n*0.01*matrix(rnorm(n*d),nrow=n,ncol=d)
    }
  D1 <- as.matrix(dist(X1))
  X1 <- X1*enorm(Do)/enorm(D1)
  s1 <- Inf
  s0 <- 2
  stepsize <-0.1
  i <- 0

while ( stepsize > 1E-5 && i < niter)
  {
    if (s1 >= s0 && i>1)
     {
       stepsize<- 0.5*stepsize
       X1 <- X0 - stepsize*normgrad
     }
    else 
    {
      stepsize <- 1.05*stepsize
      X0 <- X1
      D1mu2 <- D1^(mu-2)
      diag(D1mu2) <- 0
      D1mulam2 <- D1^(mu+1/lambda-2)
      diag(D1mulam2) <- 0
      M <- Dnu*D1mulam2-D1mu2*(Dnulam+t*(!Inb1))
      E <- matrix(rep(1,n*d),n,d)
      Grad <- X0*(M%*%E)-M%*%X0
      normgrad <- (norm(X0)/norm(Grad))*Grad
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
        s1 <- sum(Dnu*log(D1))-sum((D1mu-1)*Dnulam)/mu -t*sum((D1mu-1)*(1-Inb1))/mu
        normo <-sum(Dnu*log(Do))-sum((Domu-1)*Dnulam)/mu -t*sum((Domu-1)*(1-Inb1))/mu
        s1n <- 1-s1/normo
      }

    if(mu==0)
      {
          diag(D1)<-1
          diag(Do)<-1 #new
          s1 <- sum(Dnu*(D1mulam-1))/(mu+1/lambda) -sum(log(D1)*Dnulam)-t*sum(log(D1)*(1-Inb1))
          normo <- sum(Dnu*(Domulam-1))/(mu+1/lambda) -sum(log(Do)*Dnulam)-t*sum(log(Do)*(1-Inb1))
          s1n <- 1-s1/normo
      }

    if(mu!=0&(mu+1/lambda)!=0)
        {
            s1 <- sum(Dnu*(D1mulam-1))/(mu+1/lambda)-sum((D1mu-1)*Dnulam)/mu-t*sum((D1mu-1)*(1-Inb1))/mu
            normo <- sum(Dnu*(Domulam-1))/(mu+1/lambda)-sum((Domu-1)*Dnulam)/mu-t*sum((Domu-1)*(1-Inb1))/mu
            s1n <- 1-s1/normo
        }
    ## Printing and Plotting
     if(verbose > 3 & (i+1)%%100/verbose==0)
      {
        print (paste("niter=",i+1," stress=",round(s1,5)," stressn=",round(sqrt(s1n),5), sep=""))
      }

  }
  result <- list()
  result$conf <- X1 #new
  result$confdist <- D1
  result$delta <- Do
  result$k <- k
  result$tau <- tau
  result$theta <- c(k,tau)
  result$stress.r <- s1
  result$stress.m <- s1n
  result$stress <- sqrt(s1n)
  return(result)
}



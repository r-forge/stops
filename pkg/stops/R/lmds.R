#' An MDS version for local MDS (Chen & Buja 2006)
#'
#' @importFrom stats median
#' @importFrom stats rnorm
#' 
#' @param delta dissimilarity or distance matrix
#' @param init initial configuration. If NULL a classical scaling solution is used. 
#' @param ndim the dimension of the configuration
#' @param k the k neighbourhood parameter parameter 
#' @param tau the penalty parameter 
#' @param itmax number of optimizing iterations, defaults to 10000.
#' @param verbose prints progress if > 3. 
#'
#' @author Lisha Chen & Thomas Rusch
#' 
#' @examples
#' dis<-smacof::kinshipdelta
#' res<- lmds(as.matrix(dis),k=2,tau=0.1)
#' res
#' summary(res)
#' plot(res)
#' 
#' @export
lmds <- function(delta,init=NULL,ndim=3,k=2,tau=1,
                   itmax=5000,verbose=0)
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

  cc <- (sum(Inb1)-n)/n/n*stats::median(Dnulam[Dnulam!=0])
  t <- tau*cc
  
  Grad <- matrix (0, nrow=n, ncol= d)
  if(is.null(X1))
    {
      cmd <- cmds(Do)
      X1 <- cmd$vec[,1:d]%*%diag(cmd$val[1:d])+
        norm(Do)/n/n*0.01*matrix(stats::rnorm(n*d),nrow=n,ncol=d)
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
    diag(D1mulam) <- 0
    D1mu <- D1^mu
    diag(D1mu) <-0
#   diag(D1) <- 1
    #D0mu <- D0^mu #new
    #diag(D0) <- 1
    #diag(D0mu) <- 1
    #TR: Can't happen here as mu=lambda=1. Only needed if BC stress and LMDS are combined.
    #if(mu+1/lambda==0)
    #{
    #    diag(D1)<-1
    #    diag(Do)<-1 #new
    #    D0 <- D0+1e-6 #new: we need this as a stand in for log(0) so we make log(1e-6)
    #    D0mu <- D0^mu #new
    #    diag(D0) <- 1
    #    s1 <-   sum(Dnu*log(D1))-sum((D1mu-1)*Dnulam)/mu  -t*sum((D1mu-1)*(1-Inb1))/mu
    #    normop <-sum(Dnu*log(Do))-sum((Domu-1)*Dnulam)/mu -t*sum((Domu-1)*(1-Inb1))/mu
    #    normo0 <-sum(Dnu*log(D0))-sum((D0mu-1)*Dnulam)/mu -t*sum((D0mu-1)*(1-Inb1))/mu
    #    #s1n <- 1-s1/normo
    #    s1n <- (s1-normop)/(normo0-normop)
    #  }
#    if(mu==0)
#      {
#          diag(D1)<-1
#          diag(Do)<-1 #new
#          
#          D0 <- D0+1e-6 #new: we need this as a stand in for log(0) so we make log(1e-6)
#          diag(D0) <- 1
#          
#          s1 <-    sum(Dnu*(D1mulam-1))/(mu+1/lambda) -sum(log(D1)*Dnulam)-t*sum(log(D1)*(1-Inb1))
#          normop <- sum(Dnu*(Domulam-1))/(mu+1/lambda) -sum(log(Do)*Dnulam)-t*sum(log(Do)*(1-Inb1))
#          normo0 <- 0#sum(Dnu*(D0mulam-1))/(mu+1/lambda) -sum(log(D0)*Dnulam)-t*sum(log(D0)*(1-Inb1))
#          #s1n <- 1-s1/normo
#          s1n <- (s1-normop)/(normo0-normop)
#      }
    #if(mu!=0&(mu+1/lambda)!=0)
    #    {
    s1 <-    sum(Dnu*(D1mulam-1))/(mu+1/lambda)-sum((D1mu-1)*Dnulam)/mu-t*sum((D1mu-1)*(1-Inb1))/mu
     #   }
    ## Printing and Plotting
     if(verbose > 3 & (i+1)%%100/verbose==0)
      {
        print (paste("niter=",i+1," stress=",round(s1,5), sep=""))
      }

  }
  #For normalization
  #prelims  
  X1a <- X1*enorm(Do)/enorm(D1)
  D1a <- as.matrix(dist(X1a))
  D0 <- D1a*0
    
  D1mulama <- D1a^(mu+1/lambda)
  Domulam <- Do^(mu+1/lambda) #new
  D0mulam <- D0^(mu+1/lambda) 

  diag(D1mulama) <- 0
  diag(Domulam) <- 0
  diag(D0mulam) <- 0
    
  D1mua <- D1a^mu
  Domu <- Do^mu #new
  D0mu <- D0^mu #new

  diag(D1mua)<-0
  diag(Domu) <- 0 #new
  diag(D0mu) <- 0
    
  s1e <-    sum(Dnu*(D1mulama-1))/(mu+1/lambda)-sum((D1mua-1)*Dnulam)/mu-t*sum((D1mua-1)*(1-Inb1))/mu   #stress with normed X (X1a)  
  normop <- sum(Dnu*(Domulam-1))/(mu+1/lambda)-sum((Domu-1)*Dnulam)/mu-t*sum((Domu-1)*(1-Inb1))/mu   #best case 
  normo0 <- sum(Dnu*(D0mulam-1))/(mu+1/lambda)-sum((D0mu-1)*Dnulam)/mu-t*sum((D0mu-1)*(1-Inb1))/mu #worst case
  s1n <- (s1e-normop)/(normo0-normop) #normalized stress
 
  result <- list()
  result$conf <- X1 #new
  #result$confn <- X1a #new  
  result$confdist <- stats::as.dist(D1)
  result$delta <- stats::as.dist(Do)
  result$obsdiss <- stats::as.dist(Do)
  #result$Domulam <- Domulam  #should be the same across tau
  #result$D0mulam <- D0mulam #shoudl be the same across tau
  #result$Domu <- Domu#shoudl be the same across tau
  #result$D0mu <- D0mu#shoudl be the same across tau
  #result$D1mulam <- D1mulama #should not be the same across tau 
  #result$D1mu <- D1mu  #should not be the same across tau
  result$k <- k
  result$tau <- tau
  result$pars <- c(k,tau)
  result$stress.r <- s1
  #result$stress.e <- s1e  #for debug
  #result$normop <- normop  #for debug
  #result$normo0 <- normo0  #for debug  
  result$stress.m <- s1n
  result$call <- match.call()
  result$ndim <- ndim
  result$nobj <- n
  result$niter <- i
  result$stress <- sqrt(s1n)
  result$model<- "Local MDS"
  class(result) <- c("smacofP","smacofB","smacof")
  return(result)
}



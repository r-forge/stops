#'c-linearity
#'calculates c-linearity as the maximum multiple correlation
#'
#' @param confs a numeric matrix or data frame
#'
#' @importFrom stats lm summary.lm
#' 
#' @examples
#' x<-1:10
#' y<-2+3*x+rnorm(10)
#' z<- sin(y-x)
#' confs<-cbind(z,y,x)
#' c_linearity(confs)
#' @export
c_linearity <- function(confs)
    {
        confs <- scale(confs)
        p <- dim(confs)[2]
        tmp <- numeric()
        for(i in 1:p)
            {
             y <- confs[,i]
#        n <- dim(confs)[1]
             x <- confs[,-i]
             tmp[i] <- summary(stats::lm(y~x))$r.squared
            # rsq[i] <- summary(tmp[i])$r.squared
           }
        out <- sqrt(max(tmp))
        out
    }


#'c-dependence
#'calculates c-dependence as the distance correlation 
#'
#' @param confs a numeric matrix or data frame with two columns
#' @param index exponent on Euclidean distance, in (0,2]
#'
#'
#' @importFrom energy dcor
#' 
#' @examples
#' x<-1:10
#' y<-2+3*x+rnorm(10)
#' confs<-cbind(x,y)
#' c_dependence(confs,1.5)
#' @export
c_dependence <- function(confs,index=1)
    {
        if(dim(confs)[2]>2) stop("Distance correlation is not defined for more than two column vectors.")
        y <- confs[,1]
        x <- confs[,2]
        out <- energy::dcor(x,y,index)
        out
    }


#'c-manifoldness
#'calculates c-manifoldness as the maximal correlation coefficient
#'
#' @param confs a numeric matrix or data frame with two columns
#'
#' @importFrom acepack ace
#' @importFrom stats cor
#' 
#' @examples
#' x<--100:100
#' y<-sqrt(100^2-x^2)
#' confs<-cbind(x,y)
#' c_manifoldness(confs)
#' @export
c_manifoldness <- function(confs)
    {
        if(dim(confs)[2]>2) stop("Maximal correlation is not available for more than two column vectors. You can use cassociation.")
        #it is not symmetric so we use dim1 to explain dim2
        x <- confs[,1]
        y <- confs[,2]
        tmp <- acepack::ace(x,y)
        out <- stats::cor(tmp$tx,tmp$ty)
        out
    }


#'wrapper for getting the mine coefficients
#'
#' @param confs a numeric matrix or data frame with two columns
#' @param master the master column 
#' @param alpha an optional number of cells allowed in the X-by-Y search-grid. Default value is 1
#' @param C an optional number determining the starting point of the X-by-Y search-grid. When trying to partition the x-axis into X columns, the algorithm will start with at most C X clumps. Default value is 15. 
#' @param var.thr minimum value allowed for the variance of the input variables, since mine can not be computed in case of variance close to 0. Default value is 1e-5.
#' @param eps integer in [0,1] (?).  If NULL (default) it is set to 1-MIC. It can be set to zero for noiseless functions, but the default choice is the most appropriate parametrization for general cases (as stated in Reshef et al. SOM). It provides robustness.
#' 
#' @importFrom minerva mine
c_mine <- function(confs,master=NULL,alpha=1,C=15,var.thr=1e-5,eps=NULL)
    {
        #if(dim(confs)[2]>2) stops("MINE is not defined for more than two column vectors.")
        out <- minerva::mine(x=confs,master=master,alpha=alpha,C=C,var.thr=var.thr,eps=eps)
        out
    }

#' c-association
#' calculates the c-association based on the maximal information coefficient 
#' We define c-association as the maximum association between any two dimensions - or use the determinant of the MINE matrix - TBD
#' 
#' @param confs a numeric matrix or data frame with two columns
#' @param alpha an optional number of cells allowed in the X-by-Y search-grid. Default value is 1
#' @param C an optional number determining the starting point of the X-by-Y search-grid. When trying to partition the x-axis into X columns, the algorithm will start with at most C X clumps. Default value is 15. 
#' @param var.thr minimum value allowed for the variance of the input variables, since mine can not be computed in case of variance close to 0. Default value is 1e-5.
#' @param eps integer in [0,1] (?).  If NULL (default) it is set to 1-MIC. It can be set to zero for noiseless functions, but the default choice is the most appropriate parametrization for general cases (as stated in Reshef et al). It provides robustness.
#' 
#' @importFrom minerva mine
#' 
#' @examples
#' x<-seq(-3,3,length.out=200)
#' y<-sqrt(3^2-x^2)
#' z<- sin(y-x)
#' confs<-cbind(x,y,z)
#' c_association(confs)
#' @export
c_association <- function(confs,alpha=1,C=15,var.thr=1e-5,eps=NULL)
    {
        tmp <- c_mine(confs=confs,master=NULL,alpha=alpha,C=C,var.thr=var.thr,eps=eps)$MIC
        tmp <- tmp[lower.tri(tmp)] #to get rid of the main diagonal
        out <- max(tmp) #the question is how to aggregate this for more than two dimensions, I now use the maximum so the maximum association of any two dimensions is looked at - but perhaps a harmonic mean or even the arithmetic one might be better 
        out
    }

#' c-nonmonotonicity
#' calculates the c-nonmonotonicity based on the maximum asymmetric score 
#' We define c-nonmonotonicity as the maximum nonmonotonicity between any two dimensions - or use the determinant of the MINE matrix? - TBD
#' 
#' @param confs a numeric matrix or data frame with two columns
#' @param alpha an optional number of cells allowed in the X-by-Y search-grid. Default value is 1
#' @param C an optional number determining the starting point of the X-by-Y search-grid. When trying to partition the x-axis into X columns, the algorithm will start with at most C X clumps. Default value is 15. 
#' @param var.thr minimum value allowed for the variance of the input variables, since mine can not be computed in case of variance close to 0. Default value is 1e-5.
#' @param eps integer in [0,1] (?).  If NULL (default) it is set to 1-MIC. It can be set to zero for noiseless functions, but the default choice is the most appropriate parametrization for general cases (as stated in Reshef et al. SOM). It provides robustness.
#' 
#' @importFrom minerva mine
#' 
#' @examples
#' x<-seq(-3,3,length.out=200)
#' y<-sqrt(3^2-x^2)
#' z<- sin(y-x)
#' confs<-cbind(x,y,z)
#' c_nonmonotonicity(confs)
#' @export
c_nonmonotonicity <- function(confs,alpha=1,C=15,var.thr=1e-5,eps=NULL)
    {
        tmp <- c_mine(confs=confs,master=NULL,alpha=alpha,C=C,var.thr=var.thr,eps=eps)$MAS
        tmp <- tmp[lower.tri(tmp)] #to get rid of the main diagonal
        out <- max(tmp) #the question is how to aggregate this for more than two dimensions, I now use the maximum so the maximum association of any two dimensions is looked at - but perhaps a harmonic mean or even the arithmetic one might be better 
        out
    }

#' c-functionality
#' calculates the c-functionality based on the maximum edge value 
#' We define c-functionality as the mean functionality between any two dimensions - or use the determinant of the MINE matrix? - TBD
#' 
#' @param confs a numeric matrix or data frame with two columns
#' @param alpha an optional number of cells allowed in the X-by-Y search-grid. Default value is 1
#' @param C an optional number determining the starting point of the X-by-Y search-grid. When trying to partition the x-axis into X columns, the algorithm will start with at most C X clumps. Default value is 15. 
#' @param var.thr minimum value allowed for the variance of the input variables, since mine can not be computed in case of variance close to 0. Default value is 1e-5.
#' @param eps integer in [0,1] (?).  If NULL (default) it is set to 1-MIC. It can be set to zero for noiseless functions, but the default choice is the most appropriate parametrization for general cases (as stated in Reshef et al. SOM). It provides robustness.
#' 
#' @importFrom minerva mine
#' 
#' @examples
#' x<-seq(-3,3,length.out=200)
#' y<-sqrt(3^2-x^2)
#' z<- sin(y-x)
#' confs<-cbind(x,y,z)
#' c_functionality(confs)
#' @export
c_functionality <- function(confs,alpha=1,C=15,var.thr=1e-5,eps=NULL)
    {
        tmp <- c_mine(confs=confs,master=NULL,alpha=alpha,C=C,var.thr=var.thr,eps=eps)$MEV
        tmp <- tmp[lower.tri(tmp)] #to get rid of the main diagonal
        out <- mean(tmp) #the question is how to aggregate this for more than two dimensions, I now use the mean
        out
    }

#' c-complexity
#' calculates the c-complexity based on the minimum cell number
#' We define c-complexity as the maximum minimum cell number between any two dimensions - or use the determinant of the MINE matrix? - TBD
#'
#' @param confs a numeric matrix or data frame with two columns
#' @param alpha an optional number of cells allowed in the X-by-Y search-grid. Default value is 1
#' @param C an optional number determining the starting point of the X-by-Y search-grid. When trying to partition the x-axis into X columns, the algorithm will start with at most C X clumps. Default value is 15. 
#' @param var.thr minimum value allowed for the variance of the input variables, since mine can not be computed in case of variance close to 0. Default value is 1e-5.
#' @param eps integer in [0,1] (?).  If NULL (default) it is set to 1-MIC. It can be set to zero for noiseless functions, but the default choice is the most appropriate parametrization for general cases (as stated in Reshef et al. SOM). It provides robustness.
#' 
#' @importFrom minerva mine
#' 
#' @examples
#' x<-seq(-3,3,length.out=200)
#' y<-sqrt(3^2-x^2)
#' z<- sin(y-x)
#' confs<-cbind(x,y,z)
#' c_complexity(confs)
#' @export
c_complexity <- function(confs,alpha=1,C=15,var.thr=1e-5,eps=NULL)
    {
        tmp <- c_mine(confs=confs,master=NULL,alpha=alpha,C=C,var.thr=var.thr,eps=eps)$MCN
        tmp <- tmp[lower.tri(tmp)]
        out <- max(tmp)
        out
    }

#' c-faithfulness 
#' calculates the c-faithfulness based on the index by Chen and Buja 2013 (M_adj) with equal input neigbourhoods 
#'
#' @param confs a numeric matrix or a dist object
#' @param obsdiss a symmetric numeric matrix or a dist object
#' @param k the number of nearest neighbours to be looked at
#' @param ... additional arguments passed to dist()  
#' 
#' @examples
#' delts<-smacof::kinshipdelta
#' dis<-smacofSym(delts)$confdist
#' c_faithfulness(dis,delts,k=3)
#' @export
c_faithfulness<- function(confs,obsdiss,k=3,...)
{
    if(inherits(obsdiss,"dist")) obsdiss <- as.matrix(obsdiss)
    tdiss <- apply(obsdiss,2,sort)[k+1,] 
    nnmat <- ifelse(obsdiss>tdiss, 0, 1)
    n <- nrow(nnmat)
    kv <- apply(nnmat,1,sum)
    kvmat <- matrix(rep(kv,n),ncol=n,byrow=TRUE)
    confdist <- as.matrix(dist(confs,...))
    confrk <- apply(confdist, 2, rank)
    nnconf <- ifelse(confrk>kvmat, 0,1)
    nk <- (apply(nnconf*nnmat, 2, sum)-1)/(kv-1)
    mkadj <- mean(nk)-(sum(nnmat)-n)/n/n
    res <- list(mda=mkadj,nk=nk)
    return(res)
  }

#'calculate k nearest neighbours from a distance matrix
#' @param dis distance matrix
#' @param k number of nearest neighbours (Note that with a tie, the function returns the alphanumerically first one!)
#' @export
knn_dist <- function(dis,k)
    {
      dis <- as.matrix(dis)  
      #if(!isSymmetric(dis)) stop("Distance matrix is not symmetric.")  
      n <- nrow(dis)
      knnindex <- matrix(0, nrow = n, ncol = k+1)
      knndist  <- knnindex
      for(i in 1:n){
          knnindex[i,] = order(dis[i,])[1:(k+1)]
          knndist[i,] = dis[i,knnindex[i,]]
      }
      list(index=knnindex[,-1],distance=knndist[,-1])
   }



#' c-clusteredness 
#' calculates c-clusteredness as the OPTICS cordillera. The higher the more clustered. 
#' 
#' @param confs a numeric matrix or a dist object
#' @param q The norm used for the Cordillera. Defaults to 2. 
#' @param minpts The minimum number of points that must make up a cluster in OPTICS (corresponds to k in the paper). It is passed to \code{\link{optics}} where it is called minPts. Defaults to 2.
#' @param epsilon The epsilon parameter for OPTICS (called epsilon_max in the paper). Defaults to 2 times the maximum distance between any two points.
#' @param distmeth The distance to be computed if X is not a symmetric matrix or a dist object (otherwise ignored). Defaults to Euclidean distance. 
#' @param dmax The winsorization value for the highest allowed reachability. If used for comparisons this should be supplied. If no value is supplied, it is NULL (default), then dmax is taken from the data as minimum of epsilon or the largest reachability.
#' @param digits The precision to round the raw Cordillera and the norm factor. Defaults to 10.
#' @param scale Should X be scaled if it is an asymmetric matrix or data frame? Can take values TRUE or FALSE or a numeric value. If TRUE or 1, standardisation is to mean=0 and sd=1. If 2, no centering is applied and scaling of each column is done with the root mean square of each column. If 3, no centering is applied and scaling of all columns is done as X/max(standard deviation(allcolumns)). If 4, no centering is applied and scaling of all columns is done as X/max(rmsq(allcolumns)). If FALSE, 0 or any other numeric value, no standardisation is applied. Defaults to 4. 
#' @param ... Additional arguments to be passed to \code{\link{optics}}
#' 
#' @examples
#' delts<-smacof::kinshipdelta
#' dis<-smacofSym(delts)$confdist
#' c_clusteredness(dis,delts,k=3)
#' @export
c_clusteredness<- function(confs,minpts=2,q=2,epsilon=2*max(dist(confs)),distmeth="euclidean",dmax=NULL,digits=10,scale=4,...)
{
    out <- cordillera::cordillera(confs,minpts=minpts,q=q,epsilon=epsilon,distmeth=distmeth,dmax=dmax,digits=digits,scale=4,...)$normed
    return(out)
  }


#' c-regularity 
#' calculates c-regularity as 1 - OPTICS cordillera for k=2. The higher the more regular. 
#' 
#' @param confs a numeric matrix or a dist object
#' @param q The norm used for the Cordillera. Defaults to 2. 
#' @param epsilon The epsilon parameter for OPTICS (called epsilon_max in the paper). Defaults to 2 times the maximum distance between any two points.
#' @param distmeth The distance to be computed if X is not a symmetric matrix or a dist object (otherwise ignored). Defaults to Euclidean distance. 
#' @param dmax The winsorization value for the highest allowed reachability. If used for comparisons this should be supplied. If no value is supplied, it is NULL (default), then dmax is taken from the data as minimum of epsilon or the largest reachability.
#' @param digits The precision to round the raw Cordillera and the norm factor. Defaults to 10.
#' @param scale Should X be scaled if it is an asymmetric matrix or data frame? Can take values TRUE or FALSE or a numeric value. If TRUE or 1, standardisation is to mean=0 and sd=1. If 2, no centering is applied and scaling of each column is done with the root mean square of each column. If 3, no centering is applied and scaling of all columns is done as X/max(standard deviation(allcolumns)). If 4, no centering is applied and scaling of all columns is done as X/max(rmsq(allcolumns)). If FALSE, 0 or any other numeric value, no standardisation is applied. Defaults to 4. 
#' @param ... Additional arguments to be passed to \code{\link{optics}}
#' 
#' @examples
#' hpts <- sp:::genHexGrid(dx=0.9, ll=c(-2, -2), ur=c(2, 2))
#' plot(hpts[,1],hpts[,2],pch=19,asp=1)
#' c_regularity(hpts)
#' @export
c_regularity<- function(confs,q=2,epsilon=2*max(dist(confs)),distmeth="euclidean",dmax=NULL,digits=10,scale=4,...)
{
    out <- 1-cordillera::cordillera(confs,minpts=2,q=q,epsilon=epsilon,distmeth=distmeth,dmax=dmax,digits=digits,scale=4,...)$normed
    return(out)
  }

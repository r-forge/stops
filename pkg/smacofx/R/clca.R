#' Curvilinear Component Ananlysis.
#'
#'
#' A wrapper to \code{\link[ProjectionBasedClustering]{CCA}} to return a 'smacofP' object. Note this functionality is rather rudimentary.   
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param Epochs Scalar; gives the number of passes through the data.
#' @param lambda0 the boundary/neighbourhood parameter(s) (called lambda_y in the original paper). It is supposed to be a numeric scalar. It defaults to the 90\% quantile of delta. 
#' @param alpha0 (scalar) initial step size, 0.5 by default
#' @param weightmat not used
#' @param init starting configuration, not used
#' @param ndim dimension of the configuration; defaults to 2
#' @param acc numeric accuracy of the iteration; not used
#' @param itmax maximum number of iterations. Not used.
#' @param verbose should iteration output be printed; not used
#' @param method Distance calculation; currently not used.
#' @param principal If 'TRUE', principal axis transformation is applied to the final configuration

#'
#' @return a smacofP object. It is a list with the components
#' \itemize{
#' \item delta: Observed, untransformed dissimilarities
#' \item tdelta: Observed explicitly transformed dissimilarities, normalized
#' \item dhat: Explicitly transformed dissimilarities (dhats), optimally scaled and normalized 
#' \item confdist: Configuration dissimilarities
#' \item conf: Matrix of fitted configuration
#' \item stress: Default stress (stress-1; sqrt of explicitly normalized stress)
#' \item spp: Stress per point 
#' \item ndim: Number of dimensions
#' \item model: Name of model
#' \item niter: Number of iterations (training length)
#' \item nobj: Number of objects
#' \item type: Type of MDS model. Only ratio here.
#' \item weightmat: weighting matrix as supplied
#' \item stress.m: Default stress (stress-1^2)
#' \item tweightmat: transformed weighting matrix; it is weightmat but containing all the 0s for the distances set to 0. 
#'}
#'
#'
#' @details
#' The solution is found by "quasi-majorization", which means that the majorization is only real majorization once the weightmat no longer changes. This typically happens after a few iterations. Due to that it can be that in the beginning the stress may not decrease monotonically and that there's a chance it might never. 
#' 
#' If tau is too small it may happen that all distances for one i to all j are zero and then there will be an error, so make sure to set a larger tau.
#'
#' In the standard functions 'spmds' and 'smds' we keep tau fixed throughout. This means that if tau is large enough, then the result is the same as the corresponding MDS. In the orginal publication the idea was that of a self-organizing map which decreased tau over epochs (i.e., passes through the data). This can be achieved with our function 'so_spmds' 'so_smds' which creates a vector of decreasing tau values, calls the function 'spmds' with the first tau, then supplies the optimal configuration obtained as the init for the next call with the next tau and so on. 
#'
#' 
#' @importFrom stats dist as.dist quantile
#' @importFrom ProjectionBasedClustering CCA
#' 
#' @examples
#' dis<-smacof::morse
#' res<-clca(dis,lambda0=0.4)
#' res
#' summary(res)
#' par(mfrow=c(1,2))
#' plot(res)
#' par(mfrow=c(1,1))
#'
#' ##which d_{ij}(X) exceeded tau at convergence (i.e., have been set to 0)?
#' res$tweighmat
#' res2$tweightmat
#'
#' \dontrun{
#' ## Self-organizing map style (as in the original publication)
#' #run the som-style (p)smds 
#' sommod1<-so_spmds(dis,tau=0.2,kappa=0.5,lambda=2,epochs=20,verbose=1)
#' sommod2<-so_smds(dis,tau=0.2,epochs=20,verbose=1)
#' sommod1
#' sommod2
#' }
#' 
clca <- function (delta, Epochs=20,  alpha0 = 0.5, lambda0, ndim = 2, weightmat=1-diag(nrow(delta)), init= NULL, acc=1e-06, itmax=10000, verbose=0, method = "euclidean", principal=FALSE)
{
    cc <- match.call()
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    type <- "ratio"
    labos <- rownames(delta)
    if(missing(lambda0)) lambda0 <- 3*max(stats::sd(as.matrix(delta)))
    potency_curve <-  function(v0, vn, l) return(v0 * (vn/v0)^((0:(l - 1))/(l - 1)))
    n <- nrow (delta)
    p <- ndim
    lambo <- potency_curve(lambda0, 0.01, Epochs*n)
    if(verbose>0) cat("Minimizing",type,"CLCA Stress in", Epochs,"epochs with lambda0=",lambda0,"and stepsize alpha0", alpha0, "\n")
    deltaorig <- delta
    #delta <- delta/enorm(delta,weightmat)
    deltaold <- delta
    tmp <- ProjectionBasedClustering::CCA(DataOrDistances=delta,Epochs=Epochs,alpha0=alpha0,lambda0=lambda0,OutputDimension=ndim,PlotIt=FALSE,method=method)
    out <- list()
    xnew <- tmp$ProjectedPoints #are already normalized somehow
    #xnew <- xnew/enorm(xnew) 
    ## relabeling as they were removed in the optimal scaling
    attr(xnew,"dimnames")[[1]] <- labos
    attr(xnew,"dimnames")[[2]] <- paste("D",1:ndim,sep="")
    deltaorig <- stats::as.dist(deltaorig)
    deltaold <- stats::as.dist(deltaold)
    delta <- stats::as.dist(delta)
    dout <- stats::dist(xnew)
    weightmato <- stats::as.dist(weightmat)
    weightmat <- stats::as.dist(weightmat)
    weightmat[dout>tail(lambo,1)] <- 0
    spoint <- spp(delta, dout, weightmat)
    resmat<-spoint$resmat
    rss <- sum(spoint$resmat[lower.tri(spoint$resmat)])
    spp <- spoint$spp
    #spp <- colMeans(resmat)
     if (principal) {
        xnew_svd <- svd(xnew)
        xnew <- xnew %*% xnew_svd$v
     }
    snew <- tmp$Error
    itel  <- Epochs*n
    if(verbose>1) cat("*** Stress:",snew, "; Stress 1 (default reported):",sqrt(snew), "\n")  
    out <- list(delta=deltaorig, dhat=delta, confdist=dout, iord=sort(delta), conf = xnew, stress=sqrt(snew), spp=spp,  ndim=p, weightmat=weightmato, resmat=resmat, rss=rss, init=init, model="CLCA", niter = itel, nobj = dim(xnew)[1], type = type, call=cc, stress.m=snew, alpha = NA, sigma = snew, tdelta=deltaold, parameters=c(lambda0=lambda0,alpha0=alpha0),pars=c(lambda0=lambda0,alpha0=alpha0),theta=c(lambda0=lambda0,alpha0=alpha0),tweightmat=weightmat,lambo=lambo)
    class(out) <- c("smacofP","smacofB","smacof")
    out
   }

#' Curvilinear Distance Analysis with or without power transformations either as self-organizing or not
#'
#' An implementation of curvilinear distance analysis (CLDA) by majorization with ratio, interval and ordinal optimal scaling for dissimilarities and optional power transformations. In 'pclda' the logic is that we first transform to geodesic distance, then apply the explicit power transformation and then the implicit optimal scaling.
#'
#' There is a wrapper 'clda' where the exponents are 1, which is standard CLDA but extend to allow optimal scaling. Different from the original article the neighborhood parameter tau is kept fixed in 'pclda' and 'clda'. The functions 'so_pclda' and 'so_clda' implement the self-organising principle of the original article, where the CLCA is repeatedly fitted for a decreasing sequence of taus.
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param lambda exponent of the power transformation of the dissimilarities; defaults to 1, which is also the setup of 'clda'
#' @param kappa exponent of the power transformation of the fitted distances; defaults to 1, which is also the setup of 'clda'.
#' @param nu exponent of the power of the weighting matrix; defaults to 1 which is also the setup for 'clca'. 
#' @param tau the boundary/neighbourhood parameter(s) (called lambda in the original paper). For 'pclda' and 'clda' it is supposed to be a numeric scalar (if a sequence is supplied the maximum is taken as tau) and all the transformed fitted distances exceeding tau are set to 0 via the weightmat (assignment can change between iterations).  It defaults to the 90\% quantile of the enormed (power transformed) geodesic distances of delta. For 'so_pclca' tau is supposed to be either a user supplied decreasing sequence of taus or if a scalar the maximum tau from which a decreasing sequence of taus is generated automatically as 'seq(from=tau,to=tau/epochs,length.out=epochs)' and then used in sequence.
#' @param epsilon  Shortest dissimilarity retained.
#' @param k Number of shortest dissimilarities retained for a point. If both 'epsilon' and 'k' are given, 'epsilon' will be used.
#' @param path Method used in 'stepacross' to estimate the shortest path, with alternatives '"shortest"' and '"extended"'.
#' @param fragmentedOK  What to do if dissimilarity matrix is fragmented. If 'TRUE', analyse the largest connected group, otherwise stop with error.
#' @param type what type of MDS to fit. Currently one of "ratio", "interval" or "ordinal". Default is "ratio".
#' @param ties the handling of ties for ordinal (nonmetric) MDS. Possible are "primary" (default), "secondary" or "tertiary".
#' @param weightmat a matrix of finite weights. 
#' @param init starting configuration
#' @param ndim dimension of the configuration; defaults to 2
#' @param acc numeric accuracy of the iteration. Default is 1e-6.
#' @param itmax maximum number of iterations. Default is 10000.
#' @param verbose should iteration output be printed; if > 1 then yes
#' @param principal If 'TRUE', principal axis transformation is applied to the final configuration
#' @param epochs for 'so_pclca' and tau being scalar, it gives the number of passes through the data. The sequence of taus created is 'seq(tau,tau/epochs,length.out=epochs)'. If tau is of length >1, this argument is ignored.
#'
#' @return a smacofP object (inheriting from smacofB, see \code{\link{smacofSym}}). It is a list with the components
#' \itemize{
#' \item delta: Observed, untransformed dissimilarities
#' \item tdelta: Observed explicitly transformed dissimilarities, normalized
#' \item dhat: Explicitly transformed dissimilarities (dhats), optimally scaled and normalized 
#' \item confdist: Configuration dissimilarities
#' \item conf: Matrix of fitted configuration
#' \item stress: Default stress  (stress 1; sqrt of explicitly normalized stress)
#' \item spp: Stress per point 
#' \item ndim: Number of dimensions
#' \item model: Name of smacof model
#' \item niter: Number of iterations
#' \item nobj: Number of objects
#' \item type: Type of MDS model
#' \item weightmat: weighting matrix as supplied
#' \item stress.m: Default stress (stress-1^2)
#' \item tweightmat: transformed weighting matrix; it is weightmat but containing all the 0s for the distances set to 0. 
#'}
#'
#'
#' @details
#' The geodesic distances are calculated via 'vegan::isomapdist', see \code{\link[vegan]{isomapdist}} for a documentation of what these distances do. The functions of '(p)clda' are just a wrapper for '(p)clca' applied to the geodesic distances obtained via isomapdist. 
#' 
#' If tau is too small it may happen that all distances for one i to all j are zero and then there will be an error, so make sure to set a larger tau.
#'
#' In the standard functions 'pclda' and 'clda' we keep tau fixed throughout. This means that if tau is large enough, then the result is the same as the corresponding MDS. In the orginal publication the idea was that of a self-organizing map which decreased tau over epochs (i.e., passes through the data). This can be achieved with our function 'so_pclda' 'so_clda' which creates a vector of decreasing tau values, calls the function (p)clda with the first tau, then supplies the optimal configuration obtained as the init for the next call with the next tau and so on. 
#' 
#' If tau is too low, there will be an error. 
#'
#' 
#' @importFrom stats dist as.dist quantile
#' @importFrom smacof transform transPrep
#' @importFrom vegan isomapdist 
#' 
#' @examples
#' dis<-smacof::morse
#' res<-pclda(dis,kappa=2,lambda=2,tau=0.4,k=5,itmax=1000)
#' res
#' #already many parameters 
#' coef(res)
#'
#' res2<-clda(dis,type="interval",tau=0.4,epsilon=1,itmax=1000)
#' res2
#' summary(res)
#' par(mfrow=c(1,2))
#' plot(res)
#' plot(res2)
#' par(mfrow=c(1,1))
#'
#' ##which d_{ij}(X) exceeded tau at convergence (i.e., have been set to 0)?
#' res$tweighmat
#' res2$tweightmat
#'
#' \dontrun{
#' ## Self-organizing map style (as in the original publication)
#' #run the som-style (p)clda 
#' sommod1<-so_pclda(dis,tau=0.2,k=5,kappa=0.5,lambda=2,epochs=20,verbose=1)
#' sommod2<-so_clda(dis,tau=0.2,epsilon=1,epochs=20,verbose=1)
#' sommod1
#' sommod2
#' }
#' 
#' @export
pclda <- function (delta, lambda=1, kappa=1, nu=1, tau, type=c("ratio","interval","ordinal"), ties="primary", epsilon, k, path="shortest", fragmentedOK=FALSE, weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE, principal=FALSE) {
                                        #Isomap distances
    cc <- match.call()
    #Conditions for when epsilon or k are given
    if (!missing(epsilon) && !missing(k)) message("Both epsilon and k given, using epsilon.") 
    if (missing(epsilon) && missing(k)) epsilon <- stats::quantile(delta,0.5) 
    if(missing(epsilon) && !missing(k)) delta <- vegan::isomapdist(delta,k=k,path=path,fragmentedOK=fragmentedOK) else delta <- vegan::isomapdist(delta,epsilon=epsilon,path=path,fragmentedOK=fragmentedOK)

    #default tau
    if(missing(tau)) tau <- stats::quantile((delta^lambda/enorm(delta^lambda,weightmat)),0.9)
    #run clca
    out <- pclca(delta=delta, lambda=lambda, kappa=kappa, nu=nu, tau=tau, type=type, ties=ties, weightmat=weightmat, init=init, ndim=ndim, acc=acc, itmax=itmax, verbose=verbose-2, principal=principal)
    #postprocess
    out$model= "power CLDA"
    out$call <- cc
    isocrit <- attr(delta,"criterion")
    isocritval <- attr(delta,"critval")
    out$parameters  <- c(kappa=kappa,lambda=lambda,nu=nu,tau=tau,isocritval)
    names(out$parameters)[5] <- isocrit
    out$theta <- out$pars <- out$parameters
    class(out) <- c("smacofP","smacofB","smacof")
    out
  }


#' @rdname pclda
#' @export
clda <- function(delta, tau=stats::quantile(delta,0.9), type=c("ratio","interval","ordinal"), ties="primary", epsilon, k, path="shortest", fragmentedOK=FALSE, weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE, principal=FALSE) {
    cc <- match.call()
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("delta is not symmetric.\n")
    out <- pclda(delta=delta, lambda=1, kappa=1, nu=1, tau=tau, type=type, ties=ties, epsilon=epsilon, k=k, path=path, fragmentedOK=fragmentedOK, weightmat=weightmat, init=init, ndim=ndim, acc=acc, itmax=itmax, verbose=verbose, principal=principal)
    out$model <- "CLDA"
    out$call <- cc
    paro <- out$parameters[-(1:3)]
    out$parameters <- out$theta <- out$pars <- paro
    out
}

#' @rdname pclda
#' @export
so_pclda <- function(delta, kappa=1, lambda=1, nu=1, tau=max(delta), epochs=10, type=c("ratio","interval","ordinal"), ties="primary", epsilon, k, path="shortest", fragmentedOK=FALSE, weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE, principal=FALSE) {
    cc <- match.call()
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("delta is not symmetric.\n")
    if(length(tau)<2)
       {
         taumax <- tau
         taumin <- tau/epochs
         taus <- seq(taumax,taumin,length.out=epochs)
       } else taus <- tau
    if(any(diff(taus)>0)) taus <- sort(taus,decreasing=TRUE)
    finconf <- init
    for(i in 1:length(taus))
    {
      if(verbose>0) cat(paste0("Epoch ",i,": tau=",taus[i],"\n"))  
      tmp<-pclda(delta=delta, lambda=lambda, kappa=kappa, nu=nu, tau=taus[i], type=type, ties=ties, epsilon=epsilon, k=k, path=path, fragmentedOK=fragmentedOK, weightmat=weightmat, init=finconf, ndim=ndim, verbose=verbose-1, acc=acc, itmax=itmax, principal=principal)
      finconf<-tmp$conf
      finmod<-tmp
    }
    finmod$call  <- cc
    finmod$model  <- "SO-pCLDA"
    return(finmod)
}

#' @rdname pclda
#' @export
so_clda <- function(delta, tau=max(delta), epochs=10, type=c("ratio","interval","ordinal"), ties="primary", epsilon, k, path="shortest", fragmentedOK=FALSE, weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE, principal=FALSE) {
    cc <- match.call()
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("delta is not symmetric.\n")
    if(length(tau)<2)
       {
         taumax <- tau
         taumin <- tau/epochs
         taus <- seq(taumax,taumin,length.out=epochs)
       } else taus <- tau
    if(any(diff(taus)>0)) taus <- sort(taus,decreasing=TRUE)
    finconf <- init
    for(i in 1:length(taus))
    {
      if(verbose>0) cat(paste0("Epoch ",i,": tau=",taus[i],"\n"))  
      tmp<-clda(delta=delta, tau=taus[i], type=type, ties=ties, epsilon=epsilon, k=k, path=path, fragmentedOK=fragmentedOK, weightmat=weightmat, init=finconf, ndim=ndim, verbose=verbose-1,  acc=acc, itmax=itmax, principal=principal)
      finconf<-tmp$conf
      finmod<-tmp
    }
    finmod$call  <- cc
    finmod$model  <- "SO-pCLDA"
    return(finmod)
    }

#' Curvilinear Component Analysis with powers
#'
#' An implementation of curvilinear component analysis by majorization with ratio, interval and ordinal optimal scaling for dissimilarities and power transformations for fitted distances. Uses a repeat loop.
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param kappa power of the transformation of the fitted distances; defaults to 1 for standard stress
#' @param tau the boundary parameter. Transformed fitted distances exceeding the parameter are set to 0 via the weightmat. 
#' @param type what type of MDS to fit. Currently one of "ratio", "interval" or "ordinal". Default is "ratio".
#' @param ties the handling of ties for ordinal (nonmetric) MDS. Possible are "primary" (default), "secondary" or "tertiary".
#' @param weightmat a matrix of finite weights. 
#' @param init starting configuration
#' @param ndim dimension of the configuration; defaults to 2
#' @param acc numeric accuracy of the iteration. Default is 1e-6.
#' @param itmax maximum number of iterations. Default is 10000.
#' @param verbose should iteration output be printed; if > 1 then yes
#' @param principal If ‘TRUE’, principal axis transformation is applied to the final configuration
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
#' If tau is too small it may happen that all distances for one i to all j are zero and then there will be an error, so make sure to set a larger tau.    
#'
#' @importFrom stats dist as.dist
#' @importFrom smacof transform transPrep spp
#' @seealso \code{\link{smacofSym}}
#' 
#' @examples
#' dis<-smacof::kinshipdelta
#' res<-pclca(as.matrix(dis),type="interval",kappa=2,tau=0.4,itmax=1000)
#' res
#' summary(res)
#' plot(res)
#'
#' ##which d_{ij}(X) exceeded tau at convergence (i.e., have been set to 0)?
#' res$tweighmat
#' 
#' @export
pclca <- function (delta, kappa=1, tau=0.1, type=c("ratio","interval","ordinal"), ties="primary", weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE, principal=FALSE) {
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("delta is not symmetric.\n")
    if(inherits(weightmat,"dist") || is.data.frame(weightmat)) weightmat <- as.matrix(weightmat)
    if(!isSymmetric(weightmat)) stop("weightmat is not symmetric.\n")
    r <- kappa/2
    ## -- Setup for MDS type
    if(missing(type)) type <- "ratio"
    type <- match.arg(type, c("ratio", "interval", "ordinal",several.ok = FALSE)) 
    #    "mspline"), several.ok = FALSE)
    trans <- type
    typo <- type
    if (trans=="ratio"){
    trans <- "none"
    }
    else if (trans=="ordinal" & ties=="primary"){
    trans <- "ordinalp"
    typo <- "ordinal (primary)"
   } else if(trans=="ordinal" & ties=="secondary"){
    trans <- "ordinals"
    typo <- "ordinal (secondary)"
  } else if(trans=="ordinal" & ties=="tertiary"){
    trans <- "ordinalt"
    typo <- "ordinal (tertiary)"
  #} else if(trans=="spline"){
  #  trans <- "mspline"
  }
    if(verbose>0) cat(paste("Fitting",type,"CLCA with kappa=",kappa,"and tau=",tau,"\n"))
    n <- nrow (delta)
    normi <- 0.5
    ##normi <- n #if normi=n we can use the iord structure in plot.smacofP
    ## but the problem is we don't get the correct stress then anymore.
    p <- ndim
    if (p > (n - 1)) stop("Maximum number of dimensions is n-1!")
    if(is.null(rownames(delta))) rownames(delta) <- 1:n 
    labos <- rownames(delta) #labels
    deltaorig <- delta
    #delta <- delta^lambda
    weightmato <- weightmat
    weightmat[!is.finite(weightmat)] <- 0
    delta <- delta / enorm (delta, weightmat)
    disobj <- smacof::transPrep(as.dist(delta), trans = trans, spline.intKnots = 2, spline.degree = 2)#spline.intKnots = spline.intKnots, spline.degree = spline.degree) #FIXME: only works with dist() style object 
    ## Add an intercept to the spline base transformation
                                        #if (trans == "mspline") disobj$base <- cbind(rep(1, nrow(disobj$base)), disobj$base)
    #delta <- delta / enorm (delta, weightmat)
    deltaold <- delta
    itel <- 1
    ##Starting Configs
    xold  <- init
    if(is.null(init)) xold <- smacof::torgerson (delta, p = p)
    xstart <- xold
    xold <- xold / enorm (xold) 
    nn <- diag (n)
    dold <- sqdist (xold)
    weightmat[sqrt(dold)>tau] <- 0
   ##first optimal scaling
    eold <- as.dist(sqrt(dold))
    dhat <- smacof::transform(eold, disobj, w = as.dist(weightmat), normq = normi)
    dhatt <- dhat$res #I need the structure here to reconstruct the delta; alternatively turn all into vectors? - checked how they do it in smacof
    dhatd <- structure(dhatt, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
    #FIXME: labels
    delta <- as.matrix(dhatd)
    rold <- sum (weightmat * delta * mkPower (dold, r))
    nold <- sum (weightmat * mkPower (dold, 2 * r))
    aold <- rold / nold
    sold <- 1 - 2 * aold * rold + (aold ^ 2) * nold
    ## Optimizing
    repeat {
      p1 <- mkPower (dold, r - 1)
      p2 <- mkPower (dold, (2 * r) - 1)
      by <- mkBmat (weightmat * delta * p1)
      cy <- mkBmat (weightmat * p2)
      ga <- 2 * sum (weightmat * p2)
      be <- (2 * r - 1) * (2 ^ r) * sum (weightmat * delta)
      de <- (4 * r - 1) * (4 ^ r) * sum (weightmat)
      if (r >= 0.5) {
        my <- by - aold * (cy - de * nn)
      }
      if (r < 0.5) {
        my <- (by - be * nn) - aold * (cy - ga * nn)
      }
      xnew <- my %*% xold
      xnew <- xnew / enorm (xnew)
      dnew <- sqdist (xnew)
      weightmat <- weightmato #new
      weightmat[!is.finite(weightmat)] <- 0
      weightmat[sqrt(dnew)>tau] <- 0
      ##optimal scaling
      e <- as.dist(sqrt(dnew)) #I need the dist(x) here for interval
      dhat2 <- smacof::transform(e, disobj, w = as.dist(weightmat), normq = normi)  ## dhat update
      dhatt <- dhat2$res 
      dhatd <- structure(dhatt, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
      delta <- as.matrix(dhatd)
      #delta <- as.matrix(dhatt) #In cops this is <<- because we need to change it outside of copsf() but here we don't need that 
      rnew <- sum (weightmat * delta * mkPower (dnew, r))
      nnew <- sum (weightmat * mkPower (dnew, 2 * r))
      anew <- rnew / nnew
      snew <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
      if(is.na(snew)) #if there are issues with the values
          {
              snew <- sold
              dnew <- dold
              anew <- aold
              xnew <- xold
          }   
      if (verbose>2) {
        cat (
          formatC (itel, width = 4, format = "d"),
          formatC (
            sold,
            digits = 10,
            width = 13,
            format = "f"
          ),
          formatC (
            snew,
            digits = 10,
            width = 13,
            format = "f"
          ),
          "\n"
        )
      }
      if ((itel == itmax) || ((sold - snew) < acc))
        break ()
      itel <- itel + 1
      xold <- xnew
      dold <- dnew
      sold <- snew
      aold <- anew
    }
    xnew <- xnew/enorm(xnew)
    ## relabeling as they were removed in the optimal scaling
    rownames(delta) <- labos
    attr(xnew,"dimnames")[[1]] <- rownames(delta)
    attr(xnew,"dimnames")[[2]] <- paste("D",1:p,sep="")
    doutm <- mkPower(sqdist(xnew),r)
    deltam <- delta
    #delta <- structure(delta, Size = n, call = quote(as.dist.default(m=b)),
    #                   class = "dist", Diag = FALSE, Upper = FALSE)
    delta <- stats::as.dist(delta)
    deltaorig <- stats::as.dist(deltaorig)
    deltaold <- stats::as.dist(deltaold)
    #doute <- doutm/enorm(doutm) #this is an issue here!
    #doute <- stats::as.dist(doute)
    dout <- stats::as.dist(doutm)
    weightmatm <-weightmat
    #resmat <- weightmatm*as.matrix((delta - doute)^2) #Old version 
    #resmat <- weightmatm*as.matrix((deltam - doutm)^2)
    weightmat <- stats::as.dist(weightmatm)
    #spp <- colMeans(resmat)
    spoint <- spp(delta, dout, weightmat)
    resmat<-spoint$resmat
    rss <- sum(spoint$resmat[lower.tri(spoint$resmat)])
    spp <- spoint$spp
    #spp <- colMeans(resmat)
    if (principal) {
        xnew_svd <- svd(xnew)
        xnew <- xnew %*% xnew_svd$v
    }
    #stressen <- sum(weightmat*(doute-delta)^2)
    if(verbose>1) cat("*** Stress:",snew, "; Stress 1 (default reported):",sqrt(snew),"\n")
    #delta is input delta, tdelta is input delta with explicit transformation and normalized, dhat is dhats 
    out <- list(delta=deltaorig, dhat=delta, confdist=dout, iord=dhat2$iord.prim, conf = xnew, stress=sqrt(snew), spp=spp,  ndim=p, weightmat=weightmato, resmat=resmat, rss=rss, init=xstart, model="power CLCA", niter = itel, nobj = dim(xnew)[1], type = type, call=match.call(), stress.m=snew, alpha = anew, sigma = snew, tdelta=deltaold, parameters=c(kappa=kappa,tau=tau), pars=c(kappa=kappa,tau=tau), theta=c(kappa=kappa, tau=tau),tweightmat=weightmat)
    class(out) <- c("smacofP","smacofB","smacof")
    out
  }


#' @rdname pclca
#' @export
clca <- function(delta, tau=0.1, type=c("ratio","interval","ordinal"), ties="primary", weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE, principal=FALSE) {
    out <- pclca(delta=delta, kappa=1, tau=tau, type=type, ties=ties, weightmat=weightmat, init=init, ndim=ndim, acc=acc, itmax=itmax, verbose=verbose, principal=principal)
    out$model <- "CLCA"
    }

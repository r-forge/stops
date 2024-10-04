#' SMACOF with optimal power of dissimilarities
#'
#' An implementation to minimize stress for dissimilarities to a power by majorization with ratio optimal scaling. The optimal power transformation of the dissimilarities is found by an inner optimization step.
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param type what type of MDS to fit. Currently "ratio".
#' @param ties the handling of ties for ordinal (nonmetric) MDS. Possible are "primary" (default), "secondary" or "tertiary".
#' @param weightmat a matrix of finite weights. 
#' @param init starting configuration
#' @param ndim dimension of the configuration; defaults to 2
#' @param acc numeric accuracy of the iteration. Default is 1e-6.
#' @param itmax maximum number of iterations. Default is 10000.
#' @param verbose should iteration output be printed; if > 1 then yes
#' @param principal If 'TRUE', principal axis transformation is applied to the final configuration
#' @param interval the box constraints c(upper, lower), within which to look for the optimal power transformation lambda. Defaults to c(0,4).
#'
#' @return a 'smacofP' object (inheriting from 'smacofB', see \code{\link[smacof]{smacofSym}}). It is a list with the components
#' \itemize{
#' \item delta: Observed, untransformed dissimilarities
#' \item tdelta: Observed explicitly transformed dissimilarities, normalized
#' \item dhat: Explicitly transformed dissimilarities (dhats), optimally scaled and normalized 
#' \item confdist: Transformed fitted configuration distances
#' \item iord: optimal scaling ordering
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
#' \item tweightmat: transformed weighting matrix (here NULL)
#'}
#'
#'
#' @importFrom stats dist as.dist
#' 
#' @seealso \code{\link[smacof]{smacofSym}}
#' 
#' @examples
#' dis<-smacof::kinshipdelta
#' res<-smacofPO(as.matrix(dis),itmax=1000)
#' res
#' summary(res)
#' plot(res)
#' 

smacofPOE <- function(delta, type="ratio", weightmat=1-diag(nrow(delta)), init=NULL, ndim=2, itmax = 1000, acc = 1e-10, verbose = FALSE, principal=FALSE, interval = c(0, 4)) {
        #TODO: add ndim>2
        #TODO: optimal scaling transformations; do we need to do this only to get a good iord?
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("delta is not symmetric.\n")
    if(inherits(weightmat,"dist") || is.data.frame(weightmat)) weightmat <- as.matrix(weightmat)
    if(!isSymmetric(weightmat)) stop("weightmat is not symmetric.\n")
    ## -- Setup for MDS type
    if(missing(type)) type <- "ratio"
    type <- match.arg(type, c("ratio"))
    trans <- "none"
    nobj <- nrow(delta)
    if(is.null(rownames(delta))) rownames(delta) <- 1:n 
    labos <- rownames(delta) #labels
    deltaorig <- delta
    weightmat[!is.finite(weightmat)] <- 0
    ##TODO: check out when to use delta^r and when to use ep  
    ep <- delta ^ r
    ep <- ep / enorm (ep, weightmat)
    #disobj <- smacof::transPrep(as.dist(delta), trans = trans, spline.intKnots = 2, spline.degree = 2)
    deltaold <- ep
    xe <- init
    if (is.null(xe)) {
      dd <- ep ^ 2 #or delta instead of ep?
      rd <- rowSums(dd) / nobj
      sd <- mean(delta)
      ce <- -.5 * (dd - outer(rd, rd) + sd)
      ee <- eigen(ce)
      xe <- ee$vectors[, 1:2] %*% diag(sqrt(ee$values[1:2]))
    }
    de <- as.matrix(dist(xe))
    if (interval[1] == interval[2]) {
      r <- interval[1]
      fixed <- TRUE
    } else {
      r <- (interval[1] + interval[2]) / 2
      fixed <- FALSE
    }
    g <- function(r, delta, de, weightmat) {
        ##adapted this to enorm the delta
        ep <- delta ^ r
        ep <- ep / enorm(ep, weightmat)
        return(sum(weightmat*(ep - de) ^ 2))
    }
    #ep <- delta ^ r
    sold <- sum(weightmat*(ep - de) ^ 2)
    itel <- 1
    repeat {
      b <- -ep * ifelse(de == 0, 0, 1 / de)
      diag(b) <- -rowSums(b)
      xe <- (b %*% xe) / nobj
      de <- as.matrix(dist(xe))
      ##not sure if we need this optimal scaling as we just do ratio: iord should then simply be 1:nobj 
#     dhat2 <- smacof::transform(dist(xe), disobj, w = as.dist(weightmat), normq = 0.5)  ## dhat update
      if (!fixed) {
        r <- optimize(g, interval = interval, delta = delta, de = de, weightmat = weightmat)$minimum
      }
      ep <- delta ^ r
      ep <- ep / enorm(ep, weightmat)
#      dhatt <- dhat2$res 
#      dhatd <- structure(dhatt, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
#      ep <- as.matrix(dhatd)
      smid <- sum(weightmat*(ep - de) ^ 2)   
      snew <- sum(weightmat*(ep - de) ^ 2)
      if (verbose) {
        cat(
          "itel ",
          formatC(itel, format = "d"),
          "sold ",
          formatC(sold, digits = 6, format = "f"),
          "smid ",
          formatC(smid, digits = 6, format = "f"),
          "snew ",
          formatC(snew, digits = 6, format = "f"),
          "pow  ",
          formatC(r, digits = 6, format = "f"),
          "\n"
        )
      }
      if (((sold - snew) < acc) || (itel == itmax)) {
        break
      }
      itel <- itel + 1
      sold <- snew
    }
     if (principal) {
        xe_svd <- svd(xe)
        xe <- xe %*% xe_svd$v
    }
    out <- list(delta=deltaorig,
      tdelta = ep,
      dhat=tdelta,
      confdist = de,
      iord=seq(1,nobj,by=1),
      conf = xe,
      stress = sqrt(snew),
      spp=NULL,
      ndim=ndim,
      model= c("SMACOF Power MDS"),
      niter = itel,
      nobj=nobj,
      type=type,
      weightmat=weightmat,
      stress.m=snew,
      tweightmat=NULL,
      parameters=c(lambda=r),
      pars=c(lambda=r),
      theta=c(lambda=r),
      call=match.call())
      class(out)<-c("smacofP","smacofB")
      return(out)
}

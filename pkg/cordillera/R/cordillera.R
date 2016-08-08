#' Calculates the Optics cordillera
#'
#' Calculates the OPTICS cordillera as (Rusch et al., 2015). Needs ELKI >0.6.0 - Only tested with the ubuntu binaries. This is the old version that relied on optics; since there is now an R package for optics, the code has been refactored.
#'
#' @param confs numeric matrix or data frame. This should probably be scaled to have mean=0 and variance=1.
#' @param q  the norm of the cordillera. Defaults to 1.
#' @param minpts the minpts argument to elki. Defaults to 2.
#' @param epsilon The epsilon parameter for OPTICS. Defaults to 2 times the range of x.
#' @param rang  A range of values that makes up dmax (dmax=max-min). If used for comparions this should be supplied. If no value is supplied, it is NULL (default). Then dmax is taken from the data as the difference between the largest reachability and the smallest reachability. In the latter case the cordillera is a "goodness-of-clusteredness" given the data, in the former case it also takes into account that the gaps between certain points may be larger which stands for more structure.
#' @param ylim The borders for the cordillera plot
#' @param plot plot the reachability and the raw cordillera
#' @param digits round the raw corrdilrra and the norm factor to these digits. Defaults to 10.
#' @param path the path for storing the temporary files I/O files for optics. Defaults to the current working directory.
#' @param scale Should the confs be scaled to mean 0 and sd 1? Defaults to TRUE
#' @param ... Additional arguments to be passed to optics
#' 
#' @return A list with the elements
#'      \itemize{
#'        \item $raw... The raw cordillera
#'        \item $norm... The normalization constant
#'        \item $normfac... The normalization factor (the number of times that dmax is taken)
#'        \item $dmax... The maximum distance used for maximum structure
#'        \item $normed... The normed cordillera (raw/norm)
#'        \item $optics... The optics object
#' }
#' @section Warning: It may happen that the (normed) cordillera cannot be calculated properly (e.g. division by zero, inifnite raw cordillera, q value to high etc.). A warning will be printed and the normed cordillera is either 0, 1 (if infinity is involved) or NA. In that case one needs to check one or more of the following reachability values returned from optics, minpts, eps, the raw cordillera, dmax or the normalization factor.
#'
#' @importFrom graphics barplot lines
#' @keywords clustering multivariate
#' 
#' @export
e_cordillera <- function(confs,q=1,minpts=2,epsilon,rang=NULL,digits=10,path=getwd(),plot=FALSE,ylim,scale=TRUE,...)
    {
        if(scale) confs <- scale(confs)
        if(missing(epsilon)) epsilon <- 2*diff(range(confs))
        optres <- e_optics(confs,minpts=minpts,epsilon=epsilon,path=path,...)
        res <- optres[["clusterobjectorder"]] 
        reachind <- pmatch("reachability",res[1,]) #check where the reachabilty values are to be found 
        tmp <- sapply(strsplit(res[,reachind],split='=',fixed=TRUE), function(x) x[2]) #extract the numeric values for the reachability
        #cat(tmp,"\n")
        #TODO newer elki returns unicode \342 instead of infinity so we need to work around that, works now but we suppress coercion warning
        tmp <- suppressWarnings(as.numeric(tmp))
        tmp[is.na(tmp)] <- Inf 
        if(is.null(rang)) rang <- c(0,min(epsilon,max(tmp[is.finite(tmp)]))) #This sets the range to min to max if max < eps or eps if eps < max; so this allows to do robustness stuff.
        # is eps if no range is given and eps < max reachability or if eps < range and < max reachability
        # is max range if range is given and < max reachability, eps  
        # is max reachability if max reachabity is < eps and < max range or range is not given
        maxi <- min(max(tmp[is.finite(tmp)]),max(rang)) #robustness
        #maxi <- max(rang)
       # tmp <- ifelse(is.infinite(tmp),maxi,tmp)
        tmp[tmp>maxi] <- maxi #robustness
        if(isTRUE(plot)) {
            if(missing(ylim)) ylim <- c(0,max(rang))
            idind <- pmatch("ID",res[1,]) #check where the ids are
            indtmp <- as.numeric(sapply(strsplit(res[,idind],split='=',fixed=TRUE),function(x) x[2]))
            bp <- graphics::barplot(tmp,names.arg=indtmp,col="lightgrey",border="white",ylim=ylim)
            newpoints <- rep(c(max(tmp),min(tmp)),length.out=length(tmp))
            graphics::lines(x=bp,y=tmp,col="black",lwd=1)
        }
       reachdiff <- diff(tmp) #the distance in reachability from one point to the next, basically the enveloping of the reachability plot (no need for Pythagoras as we have constant difference between each succeissve point) -> the longer the better 
#       reachdiff <- reachdiff/(max(tmp)-min(tmp))
       n <- dim(confs)[1]
       avgsidist <- round(sum(abs(reachdiff)^q,na.rm=TRUE),digits) #raw cordillera; round to three digits
       #if(missing(range)) <- range(confs)
       mdif <- abs(rang[2]-rang[1]) #dmax 
       normfac <- 2*ceiling((n-1)/minpts) #the norm factor
       normi <- round(mdif^q*normfac,digits=digits) #the normalization contant for even division
       if(!isTRUE(all.equal((n-1)%%minpts,0))) normi <- round(normi - mdif^q,digits=digits) #the correction to the upper bound if n-1/p is not fully met
       struc <- (avgsidist/normi)^(1/q) #the normed cordillera
       if(!is.finite(struc) || (normi < 1e-8)) warning(paste("I encountered a numeric problem. Most likely there was a division by a value near zero. Check whether there is something odd with these values (e.g., they are below 1e-8): Raw cordillera",avgsidist,", normalization constant",normi,"."))
       normstruc <- min(abs(struc),1) #if someone supplied a range that is to small we output 1. Also a warning? 
       if(is.na(normstruc)) normstruc <- avgsidist
       out <- list("reachplot"=tmp,"raw"=avgsidist,"norm"=normi,"normfac"=normfac,"dmax"=mdif,"normed"=normstruc,"optics"=optres)
       class(out) <- "cordillera"
       out
    }

#' Print method for OPTICs cordilleras 
#'
#' Prints the raw and normalized cordillera
#' 
#' @param x an object of class optics
#' @param ... additional arguments passed to print
#' 
#' @export
print.cordillera <- function(x,...)
    {
      out <- c(raw=x$raw,normed=x$normed)
      print(out)
    }

#' @export
summary.cordillera <- function(object,...)
    {
      out <- list(raw=object$raw,normed=object$normed,normalization=object$norm,minpts=object$optics$minPts,epsilon=object$optics$eps)
      class(out) <- "summary.cordillera"
      out
    }

#' @export
print.summary.cordillera <- function(x,...)
    {
        cat("\n")
        cat("  OPTICS cordillera values with minpts=",x$minpts,"and epsilon=",x$epsilon,"\n")
        cat("\n")
        cat("   Raw OC:",x$raw,"\n","  Normalization:",x$normalization,"\n","  Normed OC:",x$normed,"\n")
        cat("\n")
    }

#' Plot method for OPTICS cordilleras. Deprecated. 
#'
#' Plots the reachability plot and adds the cordillera to it (as a line). In this plot the cordillera is proportional to the real value. 
#' 
#' @param x an object of class cordillera
#' @param colbp color of the barplot.
#' @param coll color of the cordillera line
#' @param liwd width of the cordillera line
#' @param ylim ylim for the barplots
#' @param legend draw legend
#' @param ... additional arguments passed to barplot or lines
#'
#' @importFrom graphics barplot plot legend par
#' 
#' @export
oldcordilleraplot <- function(x,colbp="lightgrey",coll="black",liwd=1.5,legend=FALSE,ylim,...)
  {
            tmp <- x$reachplot
            rang <- c(0,min(x$optics$eps,max(tmp[is.finite(tmp)])))
            if(missing(ylim)) ylim <- c(0,max(rang))
            res <- x$optics[["clusterobjectorder"]] 
            idind <- pmatch("ID",res[1,]) #check where the ids are
            indtmp <- as.numeric(sapply(strsplit(res[,idind],split='=',fixed=TRUE),function(x) x[2]))
            bp <- graphics::barplot(tmp,names.arg=indtmp,col=colbp,border=graphics::par("bg"),ylim=ylim,...)
            newpoints <- rep(c(max(tmp),min(tmp)),length.out=length(tmp))
            graphics::lines(x=bp,y=tmp,col=coll,lwd=liwd,...)
            if(legend) graphics::legend("topright",legend="OC",col=coll,lty=1)
     }



#' Plot method for OPTICS cordilleras 
#'
#' Plots the reachability plot and adds the cordillera to it (as a line). In this plot the cordillera is proportional to the real value. 
#' 
#' @param x an object of class cordillera
#' @param colbp color of the barplot.
#' @param coll color of the cordillera line
#' @param liwd width of the cordillera line
#' @param legend draw legend
#' @param ylim ylim for the barplots
#' @param ... additional arguments passed to barplot or lines
#'
#' @importFrom graphics barplot plot legend par
#' 
#' @export
plot.cordillera <- function(x,colbp="lightgrey",coll="black",liwd=1.5,legend=FALSE,ylim,...)
  {
            tmp <- x$reachplot
            rang <- c(0,min(x$optics$eps,max(tmp[is.finite(tmp)])))
            if(missing(ylim)) ylim <- c(0,max(rang))
            indtmp <- x$optics$order
            bp <- graphics::barplot(tmp,names.arg=indtmp,col=colbp,border=graphics::par("bg"),ylim=ylim,...)
            newpoints <- rep(c(max(tmp),min(tmp)),length.out=length(tmp))
            graphics::lines(x=bp,y=tmp,col=coll,lwd=liwd,...)
            if(legend) graphics::legend("topright",legend="OC",col=coll,lty=1)
     }


#' Calculates the Optics Cordillera 
#'
#' Calculates the OPTICS cordillera as in Rusch et al., 2015. Based on optics in dbscan package.
#'
#' @param confs numeric matrix or data frame. This should probably be scaled to have mean=0 and variance=1.
#' @param q  the norm of the cordillera. Defaults to 1.
#' @param minpts the minpts argument to elki. Defaults to 2.
#' @param epsilon The epsilon parameter for OPTICS. Defaults to 2 times the range of x.
#' @param rang  A range of values that makes up dmax (dmax=max-min). If used for comparions this should be supplied. If no value is supplied, it is NULL (default). Then dmax is taken from the data as the difference between the largest reachability and the smallest reachability. In the latter case the cordillera is a "goodness-of-clusteredness" given the data, in the former case it also takes into account that the gaps between certain points may be larger which stands for more structure.
#' @param digits round the raw corrdilrra and the norm factor to these digits. Defaults to 10.
#' @param scale Should the confs be scaled to mean 0 and sd 1? Defaults to TRUE
#' @param ... Additional arguments to be passed to optics
#' 
#' @return A list with the elements
#'      \itemize{
#'        \item $raw... The raw cordillera
#'        \item $norm... The normalization constant
#'        \item $normfac... The normalization factor (the number of times that dmax is taken)
#'        \item $dmax... The maximum distance used for maximum structure
#'        \item $normed... The normed cordillera (raw/norm)
#'        \item $optics... The optics object
#' }
#' @section Warning: It may happen that the (normed) cordillera cannot be calculated properly (e.g. division by zero, inifnite raw cordillera, q value to high etc.). A warning will be printed and the normed cordillera is either 0, 1 (if infinity is involved) or NA. In that case one needs to check one or more of the following reachability values returned from optics, minpts, eps, the raw cordillera, dmax or the normalization factor.
#'
#' @importFrom dbscan optics
#' @keywords clustering multivariate
#' 
#' @examples
#' data(iris)
#' res<-princomp(iris[,1:4])
#' #2 dim goodness-of-clusteredness with clusters of at least 2 points
#' cres2<-cordillera(res$scores[,1:2])
#' cres2
#' summary(cres2)
#' plot(cres2)
#' 
#' #4 dim goodness-of-clusteredness with clusters of at least 3 points for PCA
#' cres4<-cordillera(res$scores[,1:4],minpts=3,epsilon=13) 
#' #4 dim goodness-of-clusteredness with clusters of at least 3 points for original data
#' cres<-cordillera(iris[,1:4],minpts=3,epsilon=13,rang=c(0,cres4$dmax))
#' #There is a bit more clusteredness for the PCA result
#' summary(cres4)
#' summary(cres)
#' plot(cres4)
#' plot(cres)
#' @export
cordillera <- function(confs,q=1,minpts=2,epsilon,rang=NULL,digits=10,scale=TRUE,...)
    {
        if(scale) confs <- scale(confs)
        if(missing(epsilon)) epsilon <- 2*diff(range(confs))
        optres <- dbscan::optics(confs,minPts=minpts,eps=epsilon,...)
        #optresdbsc <- dbscan::optics(confs,minPts=minpts,eps=epsilon)
        #optresstops <- stops::optics(confs,minpts=minpts,epsilon=epsilon)
        res <- optres$order
        tmp <- optres$reachdist[res]
        tmp[is.na(tmp)] <- Inf 
        if(is.null(rang)) rang <- c(0,min(epsilon,max(tmp[is.finite(tmp)]))) #This sets the range to min to max if max < eps or eps if eps < max; so this allows to do robustness stuff.
        # is eps if no range is given and eps < max reachability or if eps < range and < max reachability
        # is max range if range is given and < max reachability, eps  
        # is max reachability if max reachabity is < eps and < max range or range is not given
       maxi <- min(max(tmp[is.finite(tmp)]),max(rang)) #robustness
       tmp[tmp>maxi] <- maxi #robustness
       reachdiff <- diff(tmp) #the distance in reachability from one point to the next, basically the enveloping of the reachability plot (no need for Pythagoras as we have constant difference between each succeissve point) -> the longer the better 
       n <- dim(confs)[1]
       avgsidist <- round(sum(abs(reachdiff)^q,na.rm=TRUE),digits) #raw cordillera; round to three digits
       mdif <- abs(rang[2]-rang[1]) #dmax 
       normfac <- 2*ceiling((n-1)/minpts) #the norm factor
       normi <- round(mdif^q*normfac,digits=digits) #the normalization contant for even division
       if(!isTRUE(all.equal((n-1)%%minpts,0))) normi <- round(normi - mdif^q,digits=digits) #the correction to the upper bound if n-1/p is not fully met
       struc <- (avgsidist/normi)^(1/q) #the normed cordillera
       if(!is.finite(struc) || (normi < 1e-8)) warning(paste("I encountered a numeric problem. Most likely there was a division by a value near zero. Check whether there is something odd with these values (e.g., they are below 1e-8): Raw cordillera",avgsidist,", normalization constant",normi,"."))
       normstruc <- min(abs(struc),1) #if someone supplied a range that is to small we output 1. Also a warning? 
       if(is.na(normstruc)) normstruc <- avgsidist
       out <- list("reachplot"=tmp,"raw"=avgsidist,"norm"=normi,"normfac"=normfac,"dmax"=mdif,"normed"=normstruc,"optics"=optres)
       class(out) <- "cordillera"
       out
    }

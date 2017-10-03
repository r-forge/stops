## #' OPTICS function
## #'
## #' A rudimentary I/O interface to OPTICS in ELKI reroots intput from R to elki via the command line, runs elki.optics on an Input file (configs.txt), sets up a temporary directory opticsout, reads in the contents of the elki output file clusterobjectorder.txt and deletes I/O files and directories. Returns the contents or clusterobjectorder.txt as a data frame.
## #' needs ELKI > 0.6.0 - Only tested with the ubuntu trusty binaries very rudimentary as of yet. API of ELKI is not fixed (will be with release 1.0.0). Is no longer used for the rest as it has been refactored to use dbscan::optics. That's why its now called e_optics.

## #' @param x numeric matrix or data frame
## #' @param epsilon The epsilon parameter for OPTICS. Defaults to 2 times the range of x.
## #' @param minpts the minpts argument to elki. Defaults to 2.
## #' @param path the path for storing the temporary files I/O files. Defaults to the current working directory.
## #' @param keep should the optics results form elki be stored in path. If TRUE, it will make a directory ./opticsout and a file config.txt
## #'
## #' @return a list with the contents of the elki output file as a data frame as element
## #' \itemize{
## #'  \item $clusterobjectorder, which is the clusterobjectorder file from ELKI. The first column is the OPTICS ordering (so the ID of the points in successive order), followed by the data in x, the next column lists the reachability followed by the predecessor of each point x[i,] in the ordering in the last column  
## #'  \item $eps
## #'  \item $minpts
## #' }
## #'
## #' @importFrom utils write.table read.table
## #' @keywords clustering multivariate
## #' @examples
## #' \donttest{
## #' data(iris)
## #' res<-e_optics(iris[,1:4],minpts=2,epsilon=100)
## #' print(res)
## #' summary(res)
## #' plot(res,withlabels=TRUE)
## #' }
## #'
## #' 
## #' @export
## e_optics <- function(x,minpts=2,epsilon,path=getwd(),keep=FALSE)
##   {
##  # TODO: add distance function handler
##  # add support for row names 
##   if(missing(epsilon)) epsilon <- 2*diff(range(x))
##   utils::write.table(x,file=paste(path,"/configs.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep=" ")
##   system(paste("elki-cli -dbc.in ",path,"/configs.txt -algorithm clustering.OPTICS -optics.minpts ",minpts," -optics.epsilon ",epsilon," -resulthandler ResultWriter -out ",path,"/opticsout -out.silentoverwrite",sep=""))
##   cluobjord <- utils::read.table(paste(path,"/opticsout/clusterobjectorder.txt",sep=""),header=FALSE,fill=TRUE,stringsAsFactors=FALSE)
##   if(!keep) {
##       unlink(paste(path,"/opticsout",sep=""),recursive=TRUE)
##       file.remove(paste(path,"/configs.txt",sep=""))
##   }
##   out <- list("clusterobjectorder"=cluobjord,"epsilon"=epsilon, "minpts"=minpts)
##   class(out) <- "opticse"
##   out
## }
## #' Print method for OPTICS results
## #'
## #' Prints the object ordering and the reachabilities as directly read in from ELKI.
## #' 
## #' @param x an object of class optics
## #' @param ... additional arguments passed to print
## #'
## #' @return a data frame with the observation order and the reachabilities (invisible)
## #' 
## #' @export
## print.opticse <- function(x,...)
##     {
##        outo <- x$clusterobjectorder[,c(1,dim(x$clusterobjectorder)[2]-1)]
##        attr(outo,"names") <- c("observation","reachability")
##        print(outo)
##        invisible(outo)
##     }
## #' Summary method for OPTICs results
## #'
## #' Displays summaries of the reachability plot. Currently its the five points summary of the reachabilities and a stem and leaf display. The latter should not be confused with the reachability plot. If you need the altter, use plot().
## #' @param object an object of class optics
## #' @param ... additional arguments passed to summary.numeric
## #'
## #' @return an object of class summary.optics wit the reachabilities, the summary and minpts and epsilon parameters 
## #' 
## #' @export
## summary.opticse <- function(object,...)
##     {
##         res <- object[["clusterobjectorder"]] 
##         reachind <- pmatch("reachability",res[1,]) #check where the reachabilty values are to be found 
##         tmp <- sapply(strsplit(res[,reachind],split='=',fixed=TRUE), function(x) x[2]) #extract the numeric values
##         tmp <- suppressWarnings(as.numeric(tmp))
##         minpts <- object$minpts
##         eps <- object$eps
##         summs <- summary(tmp,...)
##         out <- list(reachabilities=tmp,summ=summs,minpts=minpts,epsilon=eps)
##         class(out) <- "summary.opticse"
##         out
##     }
## #' Print method for OPTICS summary 
## #' Displays summaries of the reachability plot. Currently its the five points summary of the reachabilities and a stem and leaf display. The latter should not be confused with the reachability plot. If you need the latter, use plot()
## #'
## #' 
## #' @param x an object of class summary.optics
## #' @param fiven should the 5 point summary be printed. Default is TRUE.
## #' @param stemd should the stema dn leaf plot be printed. Default is TRUE.
## #' @param ... additional arguments passed to \code{\link{stem}}
## #'
## #' @importFrom graphics stem
## #' 
## #'@export
## print.summary.opticse <- function(x,fiven=TRUE,stemd=TRUE,...)
##     {
##         cat("\n")
##         cat(" An OPTICS results with minpts=",x$minpts,"and epsilon=",x$epsilon,"\n",fill=TRUE)
##         if(fiven)
##             {
##              cat(" Five Point Summary of the Minimum Reachabilities:",fill=TRUE)
##              print(x$summ)
##              cat("\n")
##          }
##         if(stemd)
##             {
##              cat(" Stem and Leaf Display of the Minimum Reachabilities:",fill=TRUE)
##              graphics::stem(x$reachabilities,...)
##              cat("\n")
##           }
##     }

## #' Plot method for OPTICS results
## #'
## #' Displays the reachability plot. Points with undefined/infinite minimum reachabilities are colored lighter by default.
## #' @param x an object of class optics
## #' @param withlabels flags whether point labels should be drawn. Defaults to TRUE.
## #' @param col the color of the bars of finite/defined reachabilities. Defaults to grey. 
## #' @param colna the color of the bars for the points with infinite undefined reachabilities. Defaults to a lighter grey. 
## #' @param border the color of the bar borders. Defaults to par("bg")
## #' @param names.arg ... The arguments to be passed as names
## #' @param ... additional arguments passed to barplot
## #'
## #' @importFrom graphics barplot par
## #' 
## #' @export
## plot.opticse <- function(x,withlabels=FALSE,col="grey55",colna="grey80",border=graphics::par("bg"),names.arg,...)
##   {
##   if(withlabels & missing(names.arg)) names.arg <- x$clusterobjectorder$V1    
##   res <- x[["clusterobjectorder"]] 
##   reachind <- pmatch("reachability",res[1,]) #check where the reachabilty values are to be found 
##   tmp <- sapply(strsplit(res[,reachind],split='=',fixed=TRUE), function(x) x[2]) #extract the numeric values
##   reachabilities<- suppressWarnings(as.numeric(tmp))
##   naind <- is.na(reachabilities)
##   reachabilities[naind] <- min(x$eps,max(reachabilities[!naind]))
##   if(withlabels) {
##       graphics::barplot(reachabilities,names.arg=names.arg,border=border,col=c(rep(colna,length(reachabilities[naind])),rep(col,length(reachabilities[!naind]))),...)
##   } else
##       {
##       graphics::barplot(reachabilities,border=border,col=c(rep(colna,length(reachabilities[naind])),rep(col,length(reachabilities[!naind]))),...)
##   }
## }

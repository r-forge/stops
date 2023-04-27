#'S3 plot method for smacofP objects
#' 
#'@param x an object of class smacofP 
#'@param plot.type String indicating which type of plot to be produced: "confplot", "resplot", "Shepard", "stressplot","transplot", "bubbleplot" (see details)
#'@param plot.dim  dimensions to be plotted in confplot; defaults to c(1, 2)
#'@param main plot title
#'@param xlab label of x axis
#'@param ylab label of y axis
#'@param xlim scale of x axis
#'@param ylim scale of y axis
#'@param col vector of colors for the points
#'@param bubscale Scaling factor (size) for the bubble plot
#'@param label.conf List with arguments for plotting the labels of the configurations in a configuration plot (logical value whether to plot labels or not, label position, label color)
#'@param identify If 'TRUE', the 'identify()' function is called internally that allows to add configuration labels by mouse click
#'@param type What type of plot should be drawn (see also 'plot')
#'@param legend Flag whether legends should be drawn for plots that have legends
#'@param legpos Position of legend in plots with legends 
#'@param pch  Plot symbol
#'@param asp  Aspect ratio; defaults to 1 so distances between x and y are represented accurately; can lead to slighlty weird looking plots if the variance on one axis is much smaller than on the other axis; use NA if the standard type of R plot is wanted where the ylim and xlim arguments define the aspect ratio - but then the distances seen are no longer accurate
#'@param loess should loess fit be added to Shepard plot 
#'@param ... Further plot arguments passed: see 'plot.smacof' and 'plot' for detailed information.
#' 
#'@details
#' \itemize{
#' \item  Configuration plot (plot.type = "confplot"): Plots the MDS configurations.
#'  \item Residual plot (plot.type = "resplot"): Plots the dissimilarities against the fitted distances.
#'  \item Linearized Shepard diagram (plot.type = "Shepard"): Diagram with the transformed observed dissimilarities against the transformed fitted distance as well as loess curve and a least squares line. The fitted lines do not have an intercept.
#'  \item Transformation Plot (plot.type = "transplot"): Diagram with the observed dissimilarities (lighter) and the transformed observed dissimilarities (darker) against the fitted distances together with the nonlinear regression curve (corresponding to the power transformations and with no intercept).
#'  \item Stress decomposition plot (plot.type = "stressplot"): Plots the stress contribution in of each observation. Note that it rescales the stress-per-point (SPP) from the corresponding smacof function to percentages (sum is 100). The higher the contribution, the worse the fit.
#'  \item Bubble plot (plot.type = "bubbleplot"): Combines the configuration plot with the point stress contribution. The larger the bubbles, the better the fit.
#' }
#'
#' @importFrom graphics plot text identify legend
#' @importFrom stats loess lm predict 
#'
#' @return no return value; just plot for class 'smacofP' (see details)
#' 
#' @export
#' @examples
#' dis<-as.matrix(smacof::kinshipdelta)
#' res<-powerStressMin(dis)
#' plot(res)
#' plot(res,"reachplot")
#' plot(res,"Shepard")
#' plot(res,"resplot")
#' plot(res,"transplot")
#' plot(res,"stressplot")
#' plot(res,"bubbleplot")
plot.smacofP <- function (x, plot.type = "confplot", plot.dim = c(1, 2), bubscale = 5, col, label.conf = list(label = TRUE, pos = 3, col = 1, cex = 0.8), identify = FALSE, type = "p", pch = 20, asp = 1, main, xlab, ylab, xlim, ylim, legend = TRUE , legpos, loess=TRUE, ...)
{
    x1 <- plot.dim[1]
    y1 <- plot.dim[2]
    if (type == "n") 
        label.conf$pos <- NULL
    if (plot.type == "confplot") {
        if(missing(col)) col <- 1
        if (missing(main)) 
            main <- paste("Configuration Plot")
        else main <- main
        if (missing(xlab)) 
            xlab <- paste("Configurations D", x1, sep = "")
        else xlab <- xlab
        if (missing(ylab)) 
            ylab <- paste("Configurations D", y1, sep = "")
        else ylab <- ylab
        if (missing(xlim)) xlim <- range(x$conf[, x1])
        if (missing(ylim)) ylim <- range(x$conf[, y1]) 
        graphics::plot(x$conf[, x1], x$conf[, y1], main = main, type = type, 
            xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, 
            pch = pch, asp = asp, col = col, ...)
        if (label.conf[[1]]) 
            graphics::text(x$conf[, x1], x$conf[, y1], labels = rownames(x$conf), 
                cex = label.conf$cex, pos = label.conf$pos, col = label.conf$col)
        if (identify) {
            graphics::identify(x$conf[, x1], x$conf[, y1], labels = rownames(x$conf), 
                cex = label.conf$cex, pos = label.conf$cex, col = label.conf$col)
        }
    }
    if (plot.type == "Shepard") {
        delts <- as.vector(x$delta)
        confd <- as.vector(x$confdist)
        if(missing(col)) col <- c("grey60","grey50","black")
        if (missing(main)) 
            main <- paste("Linearized Shepard Diagram")
        else main <- main
        if (missing(xlab)) 
            xlab <- "Transformed Dissimilarities"
        else xlab <- xlab
        if (missing(ylab)) 
            ylab <- "Transformed Configuration Distances"
        else ylab <- ylab
        if (missing(xlim)) 
            xlim <- range(as.vector(x$delta))
        if (missing(ylim))
            ylim <- range(as.vector(x$confdist))
        #delta=dhats
        #proximities=obsdiss
        #distances=confdist
        graphics::plot(delts, confd, main = main, type = "p", pch=20, cex = 0.75, xlab = xlab, ylab = ylab, col = col[1], xlim = xlim, ylim = ylim, ...)
        #graphics::plot(as.vector(x$delta), as.vector(x$confdist), main = main, type = "p", cex = 0.75, xlab = xlab, ylab = ylab, col = col[1], xlim = xlim, ylim = ylim)
        #graphics::points(as.vector(x$delta), ),col=col[2],pch=19)
        #graphics::plot(as.vector(x$delta), as.vector(x$obsdiss),col=col[2],pch=20)
        if(loess) {
                   pt <- predict(stats::loess(confd~-1+delts))
                   graphics::lines(delts[order(delts)],pt[order(delts)],col=col[2],type="b",pch=20,cex=0.25)
        }
        ptl <- predict(stats::lm(confd~-1+delts))
        graphics::lines(delts[order(delts)],ptl[order(delts)],col=col[3],type="b",pch=20,cex=0.25)
       # graphics::abline(stats::lm(x$confdist~-1+x$delta),type="b") #no intercept for fitting
    }
    if (plot.type == "transplot") {
             if(missing(col)) col <- c("grey40","grey70","grey30")#,"grey50")
             kappa <- x$pars[1]
             deltao <- as.vector(x$deltaorig)
             deltat <- as.vector(x$delta)
             dreal <- as.vector(x$confdist)^(1/kappa)
             if (missing(main)) main <- paste("Transformation Plot")
             else main <- main
             if (missing(ylab)) ylab <- "Dissimilarities"
             else xlab <- xlab
             if (missing(xlab))  xlab <- "Untransformed Configuration Distances"
             else ylab <- ylab
             if (missing(ylim))  ylim <- c(min(deltat,deltao),max(deltat,deltao))
             if (missing(xlim))  xlim <- c(min(dreal^kappa,dreal),max(dreal^kappa,dreal))
            graphics::plot(dreal, deltao, main = main, type = "p", cex = 0.75, xlab = xlab, ylab = ylab, col = col[2], xlim = xlim, ylim = ylim,pch=20)
            #graphics::plot(deltat,dreal, main = main, type = "p", cex = 0.75, xlab = ylab, ylab = xlab, col = col[2], xlim = ylim, ylim = xlim,pch=20)
            graphics::points(dreal, deltat, type = "p", cex = 0.75, col = col[1],pch=20)
            pt <- predict(stats::lm(deltat~-1+I(dreal^kappa))) #with intercept forcing thorugh 0
            #pt2 <- predict(stats::lm(deltat~I(dreal^kappa))) #with intercept not forcing thorugh 0 
            #po <- predict(stats::lm(deltao~-1+I(dreal^kappa))) #with intercept
            #lines(deltat[order(deltat)],pt[order(deltat)],col=col[1],type="b",pch=20,cex=0.5)
            #lines(deltao[order(deltao)],po[order(deltao)],col=col[2],type="b",pch=20,cex=0.5)
            #graphics::lines(dreal[order(dreal)],po[order(dreal)],col=col[4])
            graphics::lines(dreal[order(dreal)],pt[order(dreal)],col=col[3],type="b",pch=19,cex=0.1)
            #graphics::lines(dreal[order(dreal)],po[order(dreal)],col=col[4],type="b",pch=19,cex=0.25) 
            if(legend) {
                if(missing(legpos)) legpos <- "topleft" 
                graphics::legend(legpos,legend=c("Transformed","Untransformed"),col=col[1:2],pch=1)
            }
         }
if (plot.type == "resplot") {
        obsd <- as.vector(x$obsdiss)
        confd <- as.vector(x$confdist)
        if(missing(col)) col <- "darkgrey" 
        if (missing(main)) 
            main <- paste("Residual plot")
        else main <- main
        if (missing(xlab)) 
            xlab <- "Normalized Dissimilarities"
        else xlab <- xlab
        if (missing(ylab)) 
            ylab <- "Configuration Distances"
        else ylab <- ylab
        if (missing(xlim)) 
            xlim <- range(obsd)
        if (missing(ylim)) 
            ylim <- range(confd)
        graphics::plot(obsd, confd, main = main, 
            type = "p", col = col, xlab = xlab, ylab = ylab, 
            xlim = xlim, ylim = ylim, ...)
        abline(lm(confd~obsd))
    }
    if (plot.type == "stressplot") {
        if(missing(col)) col <- "lightgray"
        if (missing(main)) 
            main <- paste("Stress Decomposition Chart")
        else main <- main
        if (missing(xlab)) 
            xlab <- "Objects"
        else xlab <- xlab
        if (missing(ylab)) 
            ylab <- "Stress Proportion (%)"
        else ylab <- ylab
        spp.perc <- sort((x$spp/sum(x$spp) * 100), decreasing = TRUE)
        xaxlab <- names(spp.perc)
        if (missing(xlim)) 
            xlim1 <- c(1, length(spp.perc))
        else xlim1 <- xlim
        if (missing(ylim)) 
            ylim1 <- range(spp.perc)
        else ylim1 <- ylim
        plot(1:length(spp.perc), spp.perc, xaxt = "n", type = "p", 
            xlab = xlab, ylab = ylab, main = main, xlim = xlim1, 
            ylim = ylim1, ...)
        text(1:length(spp.perc), spp.perc, labels = xaxlab, pos = 3, cex = 0.8)
        for (i in 1:length(spp.perc)) lines(c(i, i), c(spp.perc[i],0), col=col, lty = 2)
    }
    if (plot.type == "bubbleplot") {
        if(missing(col)) col <- 1
        if (missing(main)) 
            main <- paste("Bubble Plot")
        else main <- main
        if (missing(xlab)) 
            xlab <- paste("Configurations D", x1, sep = "")
        else xlab <- xlab
        if (missing(ylab)) 
            ylab <- paste("Configurations D", y1, sep = "")
        else ylab <- ylab
        if (missing(xlim)) 
            xlim <- range(x$conf[, x1]) * 1.1
        if (missing(ylim)) 
            ylim <- range(x$conf[, y1]) * 1.1
        spp.perc <- x$spp/sum(x$spp) * 100
        bubsize <- (max(spp.perc) - spp.perc + 1)/length(spp.perc) * bubscale
        plot(x$conf, cex = bubsize, main = main, xlab = xlab, 
            ylab = ylab, xlim = xlim, ylim = ylim, ...)
        xylabels <- x$conf
        ysigns <- sign(x$conf[, y1])
        xylabels[, 2] <- (abs(x$conf[, y1]) - (x$conf[, y1] * (bubsize/50))) * ysigns
        text(xylabels, rownames(x$conf), pos = 1, cex = 0.7)
    }
 }

#'@export
summary.smacofP <- function(object,...)
    {
      spp.perc <- object$spp/sum(object$spp) * 100
      sppmat <- cbind(sort(object$spp), sort(spp.perc))
      colnames(sppmat) <- c("SPP", "SPP(%)") 
      res <- list(conf=object$conf,sppmat=sppmat)
      class(res) <- "summary.smacofP"
      res
    }

#'@export
print.summary.smacofP <- function(x,...)
    {
    cat("\n")
    cat("Configurations:\n")
    print(round(x$conf, 4))
    cat("\n\n")
    cat("Stress per point:\n")
    print(round(x$sppmat, 4))
    cat("\n")
    }

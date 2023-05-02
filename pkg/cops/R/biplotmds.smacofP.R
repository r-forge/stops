#' S3 method for smacofP objects
#'
#' @param object An object of class smacofP
#' @param extvar Data frame with external variables.
#' @param scale if TRUE  external variables are standardized internally.
#'
#' @details
#' If a model for individual differences is provided, the external
#' variables are regressed on the group stimulus space
#' configurations. In the biplot only the relative length of the
#' vectors and their direction matters. Using the vecscale argument the
#' user can control for the relative length of the vectors. If
#' ‘vecscale = NULL’, the ‘vecscale()’ function from the ‘candisc’
#' package is used which tries to automatically calculate the scale
#' factor so that the vectors approximately fill the same space as
#' the configuration. In this method vecscale should usually be smaller than in smacof by a factor of 0.1.
#'
#' @return
#' Returns an object belonging to classes ‘"mlm"’ and ‘"mdsbi"’. See ‘lm’ for details.
#' R2vec: Vector containing the R2 values. See also \code{\link[smacof]{biplotmds}} for the plot method. 
#'
#' @examples
#' ## see smacof::biplotmds for more
#' res <- powerStressMin(morse,kappa=0.5,lambda=2)
#' fitbi <- biplotmds(res, morsescales[,2:3])
#' plot(fitbi, main = "MDS Biplot", vecscale = 0.03)
#' 
biplotmds.smacofP <- function(object,extvar,scale=TRUE)
{
    #klassi <- class(object)
    class(object) <- class(object)[3]
    out <- smacof::biplotmds(object,extvar=extvar,scale=scale)
    out$conf <- object$conf
    return(out)
}

## #' S3 plot method for class mdsbiRS which simply rescales the vecscal as smacofP confs are on a different scale than the smacof ones
## #' @importFrom candisc vecscale
## plot.mdsbiRS <- function(x,vecscale = NULL, plot.dim = c(1,2), col = 1, 
##                           label.conf = list(label = TRUE, pos = 3, col = 1, cex = 0.8), 
##                           vec.conf = list(col = 1, cex = 0.8, length = 0.1), 
##                           identify = FALSE, type = "p", pch = 20, 
##                           asp = 1, main, xlab, ylab, xlim, ylim, ...)
## {
##     x1 <- plot.dim[1]
##     y1 <- plot.dim[2]
##     if (is.null(label.conf$label)) 
##         label.conf$label <- TRUE
##     if (is.null(label.conf$pos)) 
##         label.conf$pos <- 3
##     if (is.null(label.conf$col)) 
##         label.conf$col <- 1
##     if (is.null(label.conf$cex)) 
##         label.conf$cex <- 0.8
##     if (identify) 
##         label.conf$label <- FALSE
##     if (is.null(vec.conf$length)) 
##         vec.conf$length <- 0.1
##     if (is.null(vec.conf$col)) 
##         vec.conf$col <- 1
##     if (is.null(vec.conf$cex)) 
##         vec.conf$cex <- 0.8
##     if (missing(main)) 
##         main <- paste("Configuration Plot")
##     else main <- main
##     if (missing(xlab)) 
##         xlab <- paste("Dimension", x1, sep = " ")
##     else xlab <- xlab
##     if (missing(ylab)) 
##         ylab <- paste("Dimension", y1, sep = " ")
##     else ylab <- ylab
##     #class(x) <- class(x)[-1]
##     if (is.null(vecscale) && ncol(x$coef) == 1) vecscale <- 1
##     if (is.null(vecscale)) vecscale <- candisc::vecscale(x$coef)
##     #x$coef <- vecscale * x$coef
##     vecscale <- vecscale/2
##     x$coef <- vecscale * x$coef
##     temp <- rbind(x$conf, t(x$coef))
##     if (missing(xlim)) 
##         xlim <- range(temp[, x1]) * 1.1
##     if (missing(ylim)) 
##         ylim <- range(temp[, y1]) * 1.1
##     plot(x$conf[, x1], x$conf[, y1], main = main, type = type, 
##         xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, pch = pch, 
##         asp = asp, col = col, ...)
##     if (label.conf$label) {
##         if (label.conf$pos[1] == 5) {
##             thigmophobe.labels(x$conf[, x1], x$conf[, y1], labels = rownames(x$conf), 
##                 cex = label.conf$cex, text.pos = NULL, col = label.conf$col)
##         }
##         else {
##             text(x$conf[, x1], x$conf[, y1], labels = rownames(x$conf), 
##                 cex = label.conf$cex, pos = label.conf$pos, col = label.conf$col)
##         }
##     }
##     if (identify) {
##         identify(x$conf[, x1], x$conf[, y1], labels = rownames(x$conf), 
##             cex = label.conf$cex, pos = label.conf$pos, col = label.conf$col)
##     }
##     abline(h = 0, v = 0, lty = 2, col = "darkgray")
##     n <- ncol(x$coef)
##     xycoor <- t(x$coef[c(x1, y1), ])
##     posvec <- apply(xycoor, 1, sign)[2, ] + 2
##     for (i in 1:n) {
##         arrows(0, 0, x$coef[x1, i], x$coef[y1, i], length = vec.conf$length, 
##             col = vec.conf$col)
##         text(x$coef[x1, i], x$coef[y1, i], labels = colnames(x$coef)[i], 
##             col = vec.conf$col, pos = posvec[i], cex = vec.conf$cex)
##     }
## }



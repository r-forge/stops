#' cordillera: The OPTICS Cordillera
#' 
#' A package for calculating the OPTICS Cordillera. The package contains various functions, methods and classes for calculating and plotting the OPTICS Cordillera and an interface to ELKI's OPTICS.  
#'
#' The stops package provides these main functions:
#' \itemize{
#' \item cordillera() ... OPTICS Cordillera using dbscan OPTICS implementation
#' }
#'
#' Methods: 
#' For most of the objects returned by the high-level functions S3 classes and methods for standard generics were implemented, including print, summary, plot.
#'
#' References:
#' \itemize{
#' \item Rusch, T., Hornik, K., & Mair, P. (2018) Assessing and quantifying clusteredness: The OPTICS Cordillera, Journal of Computational and Graphical Statistics. 27 (1), 220-233. \doi{10.1080/10618600.2017.1349664}
#' }
#' 
#'Authors: Thomas Rusch 
#'
#'Maintainer: Thomas Rusch
#'
#'
#' @examples
#' data(CAClimateIndicatorsCountyMedian)
#' 
#' res<-princomp(CAClimateIndicatorsCountyMedian[,3:52])
#' res
#' summary(res)
#' 
#' library(scatterplot3d)
#' scatterplot3d(res$scores[,1:3])
#' 
#' irisrep3d<-res$scores[,1:3]
#' irisrep2d<-res$scores[,1:2]
#'
#' 
#' #OPTICS in dbscan version
#' library(dbscan)
#' ores<-optics(irisrep2d,minPts=15,eps=100)
#' plot(ores)
#
#' #OPTICS cordillera for the 2D representation
#' cres2d<-cordillera(irisrep2d,minpts=15)
#' cres2d
#' summary(cres2d)
#' plot(cres2d)
#'
#' #OPTICS cordillera for the 3D representation
#' cres3d<-cordillera(irisrep3d,minpts=15)
#' cres3d
#' summary(cres3d)
#' plot(cres3d)
#'
#' @keywords internal
#' @aliases cordillera-package
"_PACKAGE"
NULL

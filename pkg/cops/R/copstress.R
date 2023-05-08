#' Calculates copstress for given MDS object 
#'
#' @param obj MDS object (supported are sammon, cmdscale, smacof, rstress, powermds)
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param cordweight weight to be used for the cordillera; defaults to 0.5
#' @param q the norm of the cordillera; defaults to 1
#' @param minpts the minimum points to make up a cluster in OPTICS; defaults to 2
#' @param epsilon the epsilon parameter of OPTICS, the neighbourhood that is checked; defaults to 10
#' @param rang range of the distances (min distance minus max distance). If NULL (default) the cordillera will be normed to each configuration's maximum distance, so an absolute value of goodness-of-clusteredness. 
#' @param verbose numeric value hat prints information on the fitting process; >2 is very verbose (copstress level), >3 is extremely (up to MDS optimization level)
#' @param normed should the cordillera be normed; defaults to TRUE
#' @param scale should the configuration be scale adjusted.
#' @param init a reference configuration when doing procrustes adjustment
#' @param ... additional arguments to be passed to the cordillera function
#'
#' @return A list with the components
#' \itemize{
#'        \item copstress: the weighted loss value
#'        \item OC: the Optics cordillera value
#'        \item parameters: the parameters used for fitting (kappa, lambda)
#'        \item cordillera: the cordillera object
#' }
#' @keywords multivariate
#' @import cordillera
copstress <- function(obj,stressweight=1,cordweight=5,q=1,minpts=2,epsilon=10,rang=NULL,verbose=0,normed=TRUE,scale=c("std","sd","proc","none"),init,...)
{
        if(missing(scale)) scale <- "sd"
        stressi <- obj$stress.m
        #kappa <- obj$kappa
        #lambda <- obj$lambda
        #nu <- obj$nu
        confs <- obj$conf
        confs <- scale_adjust(confs,init,scale=scale)
        corrd <- cordillera::cordillera(confs,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scale=FALSE,...)
        struc <- corrd$raw
        maxstruc <- corrd$normi
        if(normed) {
                   struc <- corrd$normed
                   maxstruc <- 1
                   }
        ic <- stressweight*stressi - cordweight*struc
        if(verbose>0) cat("copstress =",ic,"mdsloss =",stressi,"OC =",struc,"theta =",obj$parameters,"\n")
        out <- list(copstress=ic,OC=struc,parameters=obj$parameters,cordillera=corrd)
        out
     }

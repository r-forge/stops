## #' Idea for stops function allow an arbitrary number of indices in a weighted multi-objective optimization way; for this use stoploss
## #' write stops_foo where foo is the MDS model of interest
## #' TODO: also do this with a pareto approach    
## #'  Calculate the weighted multiobjective loss function used in STOPS
## #'
## #' @param obj object returned inside a stop_* function. Uses the stress.m slot for getting the stress. 
## #' @param stressweight weight to be used for the fit measure; defaults to 1
## #' @param structures which c-structuredness indices to be included in the loss
## #' @param strucweight the weights of the structuredness indices; defaults to -1/#number of structures
## #' @param strucpars a list of parameters to be passed to the c-structuredness indices in the same order as the values in structures. If the index has no parameters or you want to use the defaults, supply NULL. (alternatively a named list that has the structure name as the element name).
## #' @param stoptype what type of weighted combination should be used? Can be 'additive' or 'multiplicative'.
## #' @param verbose verbose output
## #'
## #' @import cordillera
## #'
## #'
## #' @return a list with calculated stoploss ($stoploss), structuredness indices ($strucinidices) and hyperparameters ($parameters and $theta) 
## #' 
## #' @export
## stoplossOLD<- function(obj,stressweight=1,structures=c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality"),strucweight=rep(-1/length(structures),length(structures)),strucpars,stoptype=c("additive","multiplicative"),verbose=0)
##     {
##         if(missing(strucpars)) strucpars <- vector("list", length(structures))
##         stressi <- obj$stress.m #we use stress.m everytime
##         pars <- obj$pars
##         confs <- obj$conf 
##         if("cclusteredness"%in%structures)
##             {
##               indst <- which(structures=="cclusteredness")  
##               cclusteredness <- do.call(cordillera::cordillera,c(list(confs),strucpars[[indst]]))$normed
##             }
##          if("cregularity"%in%structures)
##             {
##               indst <- which(structures=="cregularity")  
##               cregularity <- do.call(stops::c_regularity,c(list(confs),strucpars[[indst]]))
##             }                          
##         if("clinearity"%in%structures)
##             {
##                indst <- which(structures=="clinearity")
##                clinearity <- do.call(stops::c_linearity,list(confs))
##            }
##         if("cdependence"%in%structures)
##             {
##                indst <- which(structures=="cdependence")
##                cdependence <- do.call(stops::c_dependence,c(list(confs),strucpars[[indst]])) 
##            }
##         if("cmanifoldness"%in%structures)
##             {
##                indst <- which(structures=="cmanifoldness")
##                cmanifoldness <- do.call(stops::c_manifoldness,c(list(confs)))
##            }
##         if("cassociation"%in%structures)
##             {
##                indst <- which(structures=="cassociation")
##                cassociation <- do.call(stops::c_association,c(list(confs),strucpars[[indst]]))
##            }
##         if("cnonmonotonicity"%in%structures)
##             {
##                indst <- which(structures=="cnonmonotonicity")
##                cnonmonotonicity <- do.call(stops::c_nonmonotonicity,c(list(confs),strucpars[[indst]]))
##            }
##         if("cfunctionality"%in%structures)
##             {
##                indst <- which(structures=="cfunctionality")
##                cfunctionality <- do.call(stops::c_functionality,c(list(confs),strucpars[[indst]])) 
##            }
##          if("ccomplexity"%in%structures)
##             {
##                indst <- which(structures=="ccomplexity")
##                ccomplexity <- do.call(stops::c_complexity,c(list(confs),strucpars[[indst]])) 
##            }
##         if("cfaithfulness"%in%structures)
##             {
##                indst <- which(structures=="cfaithfulness")
##                cfaithfulness <- do.call(stops::c_faithfulness,c(list(confs),strucpars[[indst]]))$mda 
##             }
##          if("chierarchy"%in%structures)
##             {
##                indst <- which(structures=="chierarchy")
##                chierarchy <- do.call(stops::c_hierarchy,c(list(confs),strucpars[[indst]]))
##             }
##          if("coutlying"%in%structures)
##             {
##                indst <- which(structures=="coutlying")
##                coutlying <- do.call(stops::c_outlying,c(list(confs),strucpars[[indst]]))
##             }
##          if("cconvexity"%in%structures)
##             {
##                indst <- which(structures=="cconvexity")
##                cconvexity <- do.call(stops::c_convexity,c(list(confs),strucpars[[indst]]))
##             }
##          if("cskinniness"%in%structures)
##             {
##                indst <- which(structures=="cskinniness")
##                cskinniness <- do.call(stops::c_skinniness,c(list(confs),strucpars[[indst]])) 
##             }
##          if("cstringiness"%in%structures)
##             {
##                indst <- which(structures=="cstringiness")
##                cstringiness <- do.call(stops::c_stringiness,c(list(confs),strucpars[[indst]])) 
##             }
##          if("csparsity"%in%structures)
##             {
##                indst <- which(structures=="csparsity")
##                csparsity <- do.call(stops::c_sparsity,c(list(confs),strucpars[[indst]])) 
##             }
##          if("cclumpiness"%in%structures)
##             {
##                indst <- which(structures=="cclumpiness")
##                cclumpiness <- do.call(stops::c_clumpiness,c(list(confs),strucpars[[indst]])) 
##             }
##          if("cstriatedness"%in%structures)
##             {
##                indst <- which(structures=="cstriatedness")
##                cstriatedness <- do.call(stops::c_striatedness,c(list(confs),strucpars[[indst]]))
##             }
##         if("cinequality"%in%structures)
##             {
##                indst <- which(structures=="cinequality")
##                cinequality <- do.call(stops::c_inequality,c(list(confs),strucpars[[indst]]))
##             }
##          if("cshepardness"%in%structures)
##             {
##               #indst <- which(structures=="cshepardness")
##               cshepardness <- do.call(stops::c_shepardness,list(obj))
##            }
##         ##TODO add more structures
##         struc <- unlist(mget(structures))
##         ic <- stressi*stressweight + sum(struc*strucweight) 
##         if (stoptype =="multiplicative") ic <- stressi^stressweight*prod(struc^strucweight) 
##         if(verbose>0) cat("stoploss =",ic,"mdsloss =",stressi,"structuredness =",struc,"parameters =",pars,"\n")
##         #TODO: return the individual values for a Pareto approach
##         out <- list(stoploss=ic,strucindices=struc,parameters=pars,theta=pars)
##         out
##     }


# Idea for stops function allow an arbitrary number of indices in a weighted multi-objective optimization way; for this use stoploss
# write stops_foo where foo is the MDS model of interest
# TODO: also do this with a pareto approach    
#' Calculate the weighted multiobjective loss function used in STOPS
#'
#' @param obj object returned inside a stop_* function. Uses the stress.m slot for getting the stress. 
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param structures which c-structuredness indices to be included in the loss
#' @param strucweight the weights of the structuredness indices; defaults to -1/#number of structures
#' @param strucpars a list of parameters to be passed to the c-structuredness indices in the same order as the values in structures. If the index has no parameters or you want to use the defaults, supply NULL. (alternatively a named list that has the structure name as the element name).
#' @param stoptype what type of weighted combination should be used? Can be 'additive' or 'multiplicative'.
#' @param verbose verbose output
#' @param registry an object of class registry. This can be used to add additional c-structuredness indices. Defaults ot the registry created via .onLoad in zzz.R 
#'
#' @import cordillera
#'
#'
#' @return a list with calculated stoploss ($stoploss), structuredness indices ($strucinidices) and hyperparameters ($parameters and $theta) 
#' 
#' @export
stoploss<- function(obj,stressweight=1,structures,strucweight=rep(-1/length(structures),length(structures)),strucpars,stoptype=c("additive","multiplicative"),verbose=0,registry=struc_reg)
    {
        if(missing(strucpars)) strucpars <- vector("list", length(structures))
        stressi <- obj$stress.m #we use stress.m everytime
        pars <- obj$pars
        confs <- obj$conf
        struc <- rep(NA,length(structures))
        for(i in 1:length(structures))
        {
          struc[i] <- do.call(registry$get_entry(structures[i])$index,c(list(confs),strucpars[[i]]))
          names(struc)[i]  <- registry$get_entry(structures[i])$name  
        }  
        #struc <- unlist(mget(structures))
        ic <- stressi*stressweight + sum(struc*strucweight) 
        if (stoptype =="multiplicative") ic <- stressi^stressweight*prod(struc^strucweight) 
        if(verbose>0) cat("stoploss =",ic,"mdsloss =",stressi,"structuredness =",struc,"parameters =",pars,"\n")
        #TODO: return the individual values for a Pareto approach
        out <- list(stoploss=ic,strucindices=struc,parameters=pars,theta=pars)
        out
    }

#' stops: structure optimized proximity scaling
#' 
#' A package for "structure optimized proximity scaling" (STOPS), a collection of methods that fit nonlinear distance transformations in multidimensional scaling (MDS) and trade-off the fit with structure considerations to find optimal parameters. This includes the COPS method. The package contains various functions, wrappers, methods and classes for fitting, plotting and displaying different MDS models (classical scaling, Sammon mapping, elastic scaling, symmetric SMACOF, spherical SMACOF, sstress, rstress, powerStress, sammon mapping with powers and elastic scaling with powers), metaheuristics of multi-objective optimization, functions for various structuredness indices (including the OPTICS cordillera) and an interface to ELKI's OPTICS.  
#'
#' The stops package provides five categories of important functions:
#'
#' Models & Algorithms:
#' \itemize{
#' \item cops() and stops() ... which fit COPS and STOPS models as described in Rusch et al. (2015) and Rusch et al. (2016). By setting cordweight or strucweight to zero they can also be used to fit metric MDS for many different models, see below.  
#' \item powerStressMin()... a workhorse for fitting s-stress, r-stress (de Leeuw, 2014), powerStress (e.g., Rusch et al., 2015a), Sammon mapping with power transformations (powersammon) and elastic scaling with power transformation (powerelastic). They can most conveniently be accessed via the cops or stops functions and setting stressweight=1 and cordweight or strucweight=0 or by the dedicated functions starting with cops_XXX where XXX is the method and setting stressweight=1 and cordweight=0. It uses the nested majorization algorithm for r-stress of De Leeuw(2014).
#' \item e_optics() ... An interface to ELKI's implementation of the OPTICS; DEPRECATED
#' }
#'
#' Structuredness Indices:
#' \itemize{
#' \item cordillera() ... OPTICS cordillera
#' \item c_linearity() ... Multiple correlation of a configuration
#' }
#' 
#' Optimization functions:
#' \itemize{
#' \item ljoptim() ... An (adaptive) version of the Luus-Jakola random search
#'}
#' Wrappers and convenience functions:
#' \itemize{
#' \item conf_adjust(): procrustes adjustment of configurations 
#' \item cmdscale(), sammon(): wrappers that return S3 objects
#' \item coploss() ... a function to calculate coploss (Rusch et al., 2015a)
#' \item cop_smacofSym(), cop_sammon(), cop_cmdscale(), cop_rstress(), cop_powerstress(),cop_smacofSphere(), cop_sammon2(), cop_elastic(), cop_sstress(), cop_powerelastic(),cop_powersammon(): cop versions of these MDS models.
#' \item stop_smacofSym(), stop_powerstress(), stop_flexsmacof(): stop versions of these MDS models.
#' \item stoploss() ... a function to calculate stoploss (Rusch et al., 2015b)
#'}
#'
#' Methods: 
#' For most of the objects returned by the high-level functions S3 classes and methods for standard generics were implemented, including print, summary, plot, plot3d, plot3dstatic.   
#'
#' References:
#' \itemize{
#' \item Rusch, T., Mair, P. \& Hornik, K. (2015) COPS: Cluster optimized proximity scaling, Report 2015/1, Discussion Paper Series, Center for Empirical Research Methods, WU Vienna University of Economics and Business.
#' \item Rusch, T., Mair, P. \& Hornik, K. (2016a) Assessing and quantifying clusteredness: The OPTICS Cordillera, Report 2016/1, Discussion Paper Series, Center for Empirical Research Methods, WU Vienna University of Economics and Business.
#' \item Rusch, T., Mair, P. \& Hornik, K. (2016b) Structure based hyperparameter selection for nonlinear dimension reduction: The Structure Optimized Proximity Scaling (STOPS) framework, Report 2016/2, Discussion Paper Series, Center for Empirical Research Methods, WU Vienna University of Economics and Business. \emph{forthcoming} 
#' }
#' 
#'Authors: Thomas Rusch, Jan de Leeuw, Patrick Mair
#'
#'Maintainer: Thomas Rusch
#'
#'
#' @examples
#' data(BankingCrisesDistances)
#' 
#' #shorthand function for COPS variant 2 (hyperparameter)
#' res<-pcops(BankingCrisesDistances[,1:69])
#' res
#' summary(res)
#' plot(res)
#' plot(res,"reachplot")
#' plot(res,"transplot")
#' plot(res,"Shepard")
#'
#'\donttest{
#' #OPTICS
#' ores<-e_optics(res$conf,minpts=2,epsilon=100)
#' ores
#' summary(ores)
#' plot(ores)
#' }
#
#' #OPTICS cordillera
#' cres<-cordillera(res$fit$conf)
#' cres
#' summary(cres)
#' plot(cres)
#'
#'\donttest{
#' #STOPS
#' strucpars<-list(c(epsilon=10,minpts=2),NULL)
#' res<-stops(BankingCrisesDistances[,1:69],
#' structures=c("cclusteredness","clinearity"),strucpars=strucpars)
#' res
#' summary(res)
#' plot(res)
#' plot(res,"Shepard")
#' }
#' @docType package
#' @name stops
NULL

#' cops: cluster optimized proximity scaling
#' 
#' Cluster optimized proximity scaling (COPS) refers to multidimensional scaling methods that aim at pronouncing the clustered appearance of the configuration. They achieve this by transforming proximities/distances with power functions and augment the fitting criterion with a clusteredness index, the OPTICS Cordillera. There are two variants: One for finding the configuration directly for given parameters (COPS-C), and one for using the augmented fitting criterion to find optimal parameters (P-COPS). The package contains various functions, wrappers, methods and classes for fitting, plotting and displaying different MDS models in a COPS framework like Torgerson scaling, SMACOF, Sammon mapping, elastic scaling, symmetric SMACOF, spherical SMACOF, sstress, rstress, powermds, power elastic scaling, power sammon mapping, powerstress. All of these models can also solely be fit as MDS with power transformations. The package further contains functions for optimization (Adaptive LJ Algorithmus, particle swarm optimization).
#'
#' The cops package provides five categories of important functions:
#'
#' Models & Algorithms:
#' \itemize{
#' \item cops() ... which fit COPS models as described in Rusch et al. (2015) and Rusch et al. (2016). By setting cordweight to zero they can also be used to fit metric MDS for many different models, see below.  
#' \item powerStressMin()... a workhorse for fitting s-stress, r-stress (de Leeuw, 2014), powerStress (e.g., Rusch et al., 2015a), Sammon mapping with power transformations (powersammon) and elastic scaling with power transformation (powerelastic). They can most conveniently be accessed via the cops functions and setting stressweight=1 and cordweight or by the dedicated functions starting with cops_XXX where XXX is the method and setting stressweight=1 and cordweight=0. It uses the nested majorization algorithm for r-stress of De Leeuw(2014).
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
#' \item cop_smacofSym(), cop_sammon(), cop_cmdscale(), cop_rstress(), cop_powerstress(), cop_smacofSphere(), cop_sammon2(), cop_elastic(), cop_sstress(), cop_powerelastic(), cop_powersammon(): cop versions of these MDS models.
#'}
#'
#' Methods: 
#' For most of the objects returned by the high-level functions S3 classes and methods for standard generics were implemented, including print, summary, plot, plot3d, plot3dstatic.   
#'
#' References:
#' \itemize{
#' \item Rusch, T., Mair, P. \& Hornik, K. (2015) COPS: Cluster optimized proximity scaling, Report 2015/1, Discussion Paper Series, Center for Empirical Research Methods, WU Vienna University of Economics and Business.
#' \item Rusch, T., Mair, P. \& Hornik, K. (2016a) Assessing and quantifying clusteredness: The OPTICS Cordillera, Report 2016/1, Discussion Paper Series, Center for Empirical Research Methods, WU Vienna University of Economics and Business.
#' }
#' 
#'Authors: Thomas Rusch, Jan de Leeuw, Patrick Mair
#'
#'Maintainer: Thomas Rusch
#'
#'
#' @examples
#' library(cordillera)
#' data(BankingCrisesDistances)
#'
#' #shorthand function for COPS-C (finding configuration with coploss)
#' res<-coplossMin(BankingCrisesDistances[,1:69],stressweight=0.98,cordweight=0.02)
#' res
#' summary(res)
#' plot(res)
#' plot(res,"reachplot")
#' plot(res,"transplot")
#' plot(res,"Shepard")

#' #shorthand function for P-COPS (hyperparameter search for powerstress)
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
#' ores<-cordillera::e_optics(res$conf,minpts=2,epsilon=100)
#' ores
#' summary(ores)
#' plot(ores)
#' }
#
#' #OPTICS cordillera
#' cres<-cordillera::cordillera(res$fit$conf)
#' cres
#' summary(cres)
#' plot(cres)
#'
#' @docType package
#' @name cops
NULL

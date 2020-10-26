#' cops: cluster optimized proximity scaling
#' 
#' Cluster optimized proximity scaling (COPS) refers to multidimensional scaling methods that aim at pronouncing the clustered appearance of the configuration. They achieve this by transforming proximities/distances with power functions and augment the fitting criterion with a clusteredness index, the OPTICS Cordillera (Rusch, Hornik & Mair 2018). There are two variants: One for finding the configuration directly for given parameters (COPS-C), and one for using the augmented fitting criterion to find optimal parameters for the power transformations (P-COPS). The package contains various functions, wrappers, methods and classes for fitting, plotting and displaying different MDS models in a COPS framework like Torgerson scaling, SMACOF, Sammon mapping, elastic scaling, symmetric SMACOF, spherical SMACOF, sstress, rstress, powermds, power elastic scaling, power sammon mapping, powerstress. All of these models can also solely be fit as MDS with power transformations. The package further contains functions for optimization (Adaptive LJ Algorithmus).
#'
#' The cops package provides five categories of important functions:
#'
#' Models & Algorithms:
#' \itemize{
#' \item cops() ... high level interface to fit COPS models as described in Rusch et al. (2020). By setting cordweight to zero they can also be used to fit metric MDS for many different models, see below.
#' \item copstressMin()... The workhorse for fitting a COPS-C model. Can also be called directly.
#' \item pcops()... The workhorse for fitting a P-COPS model. Can also be called directly.
#' \item powerStressMin()... a workhorse for fitting s-stress, r-stress (de Leeuw, 2014), p-stress (e.g., Rusch et al., 2015), Sammon mapping with power transformations (powersammon) and elastic scaling with power transformation (powerelastic). They can conveniently also be fitted via the cops functions and setting stressweight=1 and cordweight or by the dedicated functions starting with cops_XXX where XXX is the method and setting stressweight=1 and cordweight=0. It uses the nested majorization algorithm for r-stress of De Leeuw (2014).
#' }
#'
#' Optimization functions:
#' \itemize{
#' \item ljoptim() ... An (adaptive) version of the Luus-Jakola random search
#'}
#' 
#' Wrappers and convenience functions:
#' \itemize{
#' \item conf_adjust(): procrustes adjustment of configurations 
#' \item cmdscale(), sammon(): wrappers for that return S3 objects to be used with cops
#' \item copstress() ... a function to calculate copstress (Rusch et al., 2020)
#' \item cop_smacofSym(), cop_sammon(), cop_cmdscale(), cop_rstress(), cop_powerstress(), cop_smacofSphere(), cop_sammon2(), cop_elastic(), cop_sstress(), cop_powerelastic(), cop_powersammon(): cop versions of these MDS models.
#'}
#'
#' Methods: 
#' For most of the objects returned by the high-level functions S3 classes and methods for standard generics were implemented, including print, summary, plot, plot3d, plot3dstatic.   
#'
#' References:
#' \itemize{
#' \item Rusch, T., Mair, P. & Hornik, K. (2015) COPS: Cluster optimized proximity scaling, Report 2015/1, Discussion Paper Series, Center for Empirical Research Methods, WU Vienna University of Economics and Business.
#' \item Rusch, T., Mair, P. & Hornik, K. (2018) Assessing and quantifying clusteredness: The OPTICS Cordillera. Journal of Computational and Graphical Statistics, 27 (1), 220-233. \url{http://dx.doi.org/10.1080/10618600.2017.1349664}
#' \item Rusch, T., Mair, P. & Hornik, K. (2020) Cluster optimized proximity scaling. Journal of Computational and Graphical Statistics, forthcoming.
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
#' \donttest{
#' #shorthand function for COPS-C (finding configuration with copstress)
#' res<-cops(BankingCrisesDistances[,1:69],variant="COPS-C",
#'           stressweight=0.98,cordweight=0.02)
#' res
#' summary(res)
#' plot(res)
#' plot(res,"reachplot")
#' plot(res,"transplot")
#' plot(res,"Shepard")

#' #shorthand function for P-COPS (hyperparameter search for powerstress)
#' res<-cops(BankingCrisesDistances[,1:69],variant="P-COPS")
#' res
#' summary(res)
#' plot(res)
#' plot(res,"reachplot")
#' plot(res,"transplot")
#' plot(res,"Shepard")
#'}
#'
#' @docType package
#' @name cops
NULL

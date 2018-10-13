#' stops: structure optimized proximity scaling
#' 
#' A package for "structure optimized proximity scaling" (STOPS), a collection of methods that fit nonlinear distance transformations in multidimensional scaling (MDS) and trade-off the fit with structure considerations to find optimal parameters or optimal configurations. This includes the three variants of cluster optimized proximity scaling (COPS). The package contains various functions, wrappers, methods and classes for fitting, plotting and displaying different MDS models in a STOPS framework like Torgerson scaling, SMACOF, Sammon mapping, elastic scaling, symmetric SMACOF, spherical SMACOF, sstress, rstress, powermds, power elastic scaling, power sammon mapping, powerstress, isomap. All of these models can also be fit as MDS variants (i.e., no structuredness). The package further contains functions for optimization (Adaptive LJ and for Bayesian optimization with treed Gaussian process with jump to linear models) and functions for various structuredness indices
#'
#' The stops package provides five categories of important functions:
#'
#' Models & Algorithms:
#' \itemize{
#' \item stops() ... which fits STOPS models as described in Rusch et al. (2018). By setting cordweight or strucweight to zero they can also be used to fit metric MDS for many different models, see below.  
#' \item powerStressMin()... a workhorse for fitting s-stress, r-stress (de Leeuw, 2014), Sammon mapping with power transformations (powersammon) and elastic scaling with power transformation (powerelastic). They can most conveniently be accessed via the cops or stops functions and setting stressweight=1 and cordweight or strucweight=0 or by the dedicated functions starting with cops_XXX where XXX is the method and setting stressweight=1 and cordweight=0. It uses the nested majorization algorithm for r-stress of De Leeuw(2014).
#' }
#'
#' Structuredness Indices:
#' \itemize{
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
#' \item stop_smacofSym(), stop_sammon(), stop_cmdscale(), stop_rstress(), stop_powerstress(),stop_smacofSphere(), stop_sammon2(), stop_elastic(), stop_sstress(), stop_powerelastic(), stop_powersammon(), stop_isomap(): stop versions of these MDS models.
#' \item stoploss() ... a function to calculate stoploss (Rusch et al., 2018)
#'}
#'
#' Methods: 
#' For most of the objects returned by the high-level functions S3 classes and methods for standard generics were implemented, including print, summary, plot, plot3d, plot3dstatic.   
#'
#' References:
#' \itemize{
#' \item Rusch, T., Mair, P. \& Hornik, K. (2018) Structure based hyperparameter selection for nonlinear dimension reduction: The Structure Optimized Proximity Scaling (STOPS) framework, Report 2018/2, Discussion Paper Series, Center for Empirical Research Methods, WU Vienna University of Economics and Business. \emph{forthcoming} 
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
#'\donttest{
#' #STOPS
#' strucpars<-list(c(epsilon=10,minpts=2),NULL)
#' res<-stops(BankingCrisesDistances[,1:69],loss="strain",
#' structures=c("cclusteredness","clinearity"),strucpars=strucpars,optimmethod="Kriging")
#' res
#' summary(res)
#' plot(res)
#' plot(res,"Shepard")
#' }
#' @docType package
#' @name stops
NULL

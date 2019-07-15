#' stops: structure optimized proximity scaling
#' 
#' A package for "structure optimized proximity scaling" (STOPS), a collection of methods that fit nonlinear distance transformations in multidimensional scaling (MDS) and trade-off the fit with structure considerations to find optimal parameters or optimal configurations. This includes the three variants of cluster optimized proximity scaling (COPS). The package contains various functions, wrappers, methods and classes for fitting, plotting and displaying different MDS models in a STOPS framework like Torgerson scaling, SMACOF, Sammon mapping, elastic scaling, symmetric SMACOF, spherical SMACOF, sstress, rstress, powermds, power elastic scaling, power sammon mapping, powerstress, isomap. All of these models can also be fit as MDS variants (i.e., no structuredness). The package further contains functions for optimization (Adaptive LJ and for Bayesian optimization with treed Gaussian process with jump to linear models) and functions for various structuredness indices
#'
#' The stops package provides five categories of important functions:
#'
#' Models & Algorithms:
#' \itemize{
#' \item stops() ... which fits STOPS models as described in Rusch et al. (2019). By setting cordweight or strucweight to zero they can also be used to fit metric MDS for many different models, see below.  
#' \item powerStressMin()... a workhorse for fitting s-stress, r-stress (De Leeuw, 2014), Sammon mapping with power transformations (powersammon) and elastic scaling with power transformation (powerelastic). They can most conveniently be accessed via the stops functions and setting stressweight=1 and cordweight or strucweight=0 or by the dedicated functions starting with stops_foo where foo is the method and setting stressweight=1 and strucweight=0. It uses the nested majorization algorithm for r-stress of De Leeuw(2014).
#' \item bcStressMin()... a workhorse for fitting Box-Cox stress (Chen & Buja, 2013).
#' \item lmds()... a workhorse for the local MDS of Chen & Buja (2008).
#' }
#'
#' 
#' Structuredness Indices:
#' Various c-structuredness as c_foo(), where foo is the name of the structuredness.  See Rusch et al. (2019). 
#' 
#' Optimization functions:
#' \itemize{
#' \item ljoptim() ... An (adaptive) version of the Luus-Jakola random search
#'}
#' Wrappers and convenience functions:
#' \itemize{
#' \item conf_adjust(): procrustes adjustment of configurations 
#' \item cmdscale(), sammon(): wrappers that return S3 objects
#' \item stop_smacofSym(), stop_sammon(), stop_cmdscale(), stop_rstress(), stop_powerstress(),stop_smacofSphere(), stop_sammon2(), stop_elastic(), stop_sstress(), stop_powerelastic(), stop_powersammon(),  stop_powermds(), stop_isomap(), stop_isomapeps(), stop_bcstress(), stops_lmds(): stop versions of these MDS models.
#' \item stoploss() ... a function to calculate stoploss (Rusch et al., 2019)
#'}
#'
#' Methods: 
#' For most of the objects returned by the high-level functions S3 classes and methods for standard generics were implemented, including print, summary, plot, plot3d, plot3dstatic.   
#'
#' References:
#' \itemize{
#' \item Rusch, T., Mair, P. \& Hornik, K. (2019) Structure based hyperparameter selection for nonlinear dimension reduction: The Structure Optimized Proximity Scaling (STOPS) framework, Report 2019/1, Discussion Paper Series, Center for Empirical Research Methods, WU Vienna University of Economics and Business. \emph{forthcoming} 
#' }
#' 
#'Authors: Thomas Rusch, Lisha Chen, Jan de Leeuw, Patrick Mair, Kurt Hornik
#'
#'Maintainer: Thomas Rusch
#'
#'
#' @examples
#' data(kinshipdelta)
#'
#'\donttest{
#' 
#' strucpars<-list(list(c(epsilon=10,minpts=2)),NULL)
#' dissm<-as.matrix(kinshipdelta)
#' #STOPS with strain
#' resstrain<-stops(dissm,loss="strain",
#' structures=c("cclusteredness","cdependence"),strucpars=strucpars,optimmethod="ALJ",lower=0,upper=10)
#' resstrain
#' summary(resstrain)
#' plot(resstrain)
#' plot(resstrain,"Shepard")
#'
#' #STOPS with stress
#' resstress<-stops(dissm,loss="stress",
#' structures=c("cclusteredness","cdependence"),strucpars=strucpars,optimmethod="ALJ",lower=0,upper=10)
#' resstress
#' summary(resstress)
#' plot(resstress)
#' plot(res,"Shepard")
#'
#' #STOPS with powerstress
#' respstress<-stops(dissm,,loss="powerstress",
#' structures=c("cclusteredness","cdependence"),strucpars=strucpars,weightmat=dissm,optimmethod="ALJ",lower=c(0,0,1),upper=c(10,10,10))
#' respstress
#' summary(respstress)
#' plot(respstress)
#' plot(respstress,"Shepard")
#'
#' #STOPS with bcstress
#' resbcstress<-stops(dissm,loss="bcstress",
#' structures=c("cclusteredness","cdependence"),strucpars=strucpars,optimmethod="ALJ", lower=c(0,1,0),upper=c(10,10,10))
#' resbcstress
#' summary(resbcstress)
#' plot(resbcstress)
#' plot(resbcstress,"Shepard")
#'
#' #STOPS with lmds
#' reslmds<-stops(dissm,loss="lmds",
#' structures=c("cclusteredness","clinearity"),strucpars=strucpars,optimmethod="ALJ",lower=c(2,0),upper=c(10,2))
#' reslmds
#' summary(reslmds)
#' plot(reslmds)
#' plot(reslmds,"Shepard")
#'
#' #STOPS with Isomap (the epsilon version)
#' resiso<-stops(dissm,loss="isomapeps",
#' structures=c("cclusteredness","clinearity"),strucpars=strucpars,optimmethod="ALJ",lower=47,upper=120)
#' resiso
#' summary(resiso)
#' plot(resiso)
#' plot(resiso,"Shepard")
#' }
#' @docType package
#' @name stops
NULL

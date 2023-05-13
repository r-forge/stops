#' stops: structure optimized proximity scaling
#' 
#' A package for "structure optimized proximity scaling" (STOPS), a collection of methods that fit nonlinear distance transformations in multidimensional scaling (MDS) and trade-off the fit with structure considerations to find optimal parameters or optimal configurations. The package contains various functions, wrappers, methods and classes for fitting, plotting and displaying different MDS models in a STOPS framework like Torgerson scaling, SMACOF, Sammon mapping, elastic scaling, symmetric SMACOF, spherical SMACOF, sstress, rstress, powermds, power elastic scaling, power sammon mapping, power stress, Isomap, approximate power stress, restricted power stress. All of these models can also be fit as MDS variants (i.e., no structuredness). The package further contains functions for optimization (Adaptive Luus-Jaakola and for Bayesian optimization with treed Gaussian process with jump to linear models) and functions for various structuredness indices
#'
#' The stops package provides five categories of important functions:
#'
#' Models & Algorithms:
#' \itemize{
#' \item stops() ... which fits STOPS models as described in Rusch et al. (2023). By setting cordweight or strucweight to zero they can also be used to fit metric MDS for many different models, see below.  
#' \item powerStressMin()... a workhorse for fitting many stresses, including s-stress, r-stress (De Leeuw, 2014), Sammon mapping with power transformations (powersammon), elastic scaling with power transformation (powerelastic), power stress. They can most conveniently be accessed via the stops functions and setting stressweight=1 and cordweight or strucweight=0 or by the dedicated functions starting with stop_foo where foo is the method and setting stressweight=1 and strucweight=0. It uses the nested majorization algorithm for r-stress of De Leeuw(2014).
#' \item bcStressMin()... a workhorse for fitting Box-Cox stress (Chen & Buja, 2013).
#' \item lmds()... a workhorse for the local MDS of Chen & Buja (2008).
#' }
#'
#' 
#' Structuredness Indices:
#' Various c-structuredness as c_foo(), where foo is the name of the structuredness.  See Rusch et al. (2023). 
#' 
#' Optimization functions:
#' \itemize{
#' \item ljoptim() ... An (adaptive) version of the Luus-Jakola random search
#'}
#' Wrappers and convenience functions:
#' \itemize{
#' \item conf_adjust(): procrustes adjustment of configurations 
#' \item cmdscale(), sammon(): wrappers that return S3 objects
#' \item stop_smacofSym(), stop_sammon(), stop_cmdscale(), stop_rstress(), stop_powerstress(),stop_smacofSphere(), stop_sammon2(), stop_elastic(), stop_sstress(), stop_powerelastic(), stop_powersammon(),  stop_powermds(), stop_isomap(), stop_isomapeps(), stop_bcstress(), stop_lmds(), stop_apstress(),stops_rpowerstress(): stop versions of these MDS models.
#' \item stoploss() ... a function to calculate stoploss (Rusch et al., 2023)
#'}
#'
#' Methods: 
#' For most of the objects returned by the high-level functions S3 classes and methods for standard generics were implemented, including print, summary, plot, plot3d, plot3dstatic.   
#'
#' References:
#' \itemize{
#' \item Rusch, T., Mair, P., & Hornik, K. (2023). Structure-based hyperparameter selection with Bayesian optimization in multidimensional scaling. Statistics & Computing, 33, [28]. https://doi.org/10.1007/s11222-022-10197-w
#' }
#' 
#'Authors: Thomas Rusch, Lisha Chen, Jan de Leeuw, Patrick Mair, Kurt Hornik
#'
#'Maintainer: Thomas Rusch
#'
#'
#' @examples
#' #data(Swissroll)
#' #dissm<-as.matrix(dist(Swissroll[,1:3]))
#' #cols<-Swissroll[,4]
#'
#' dissm<-as.matrix(smacof::morse)
#'
#' #Setting up structurenedness parameters
#' strucpars<-list(list(epsilon=10,minpts=2,scale=3),list(NULL))
#'
#' #STOPS with strain
#' resstrain<-stops(dissm,loss="strain",theta=1,structures=c("cclusteredness","cdependence"),
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5,itmax=10)
#' resstrain
#' summary(resstrain)
#' plot(resstrain)
#' #Fun fact: With strain clinearity must be 0 as the
#' #two principal axes are orthogonal
#' #and this can't be changed by taking
#' #the dissimilarities to a power
#' 
#'
#' \donttest{
#' #STOPS with stress or smacofSym
#' resstress<-stops(dissm,loss="stress",theta=1,
#' structures=c("cclusteredness","cdependence"),
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5)
#' resstress
#' summary(resstress)
#' plot(resstress)
#' plot(resstress,"Shepard")
#'
#' #STOPS with smacofSphere
#' ressph<-stops(dissm,loss="smacofSphere",theta=1,
#' structures=c("cclusteredness","cdependence"),
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5)
#' ressph
#' summary(ressph)
#' plot(ressph)
#'
#' #STOPS with sammon
#' ressam<-stops(dissm,loss="sammon",theta=1,
#' structures=c("cclusteredness","cdependence"),
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5)
#' ressam
#' summary(ressam)
#' plot(ressam)
#' plot(ressam,"transplot")
#'
#' #STOPS with sammon2
#' ressam<-stops(dissm,loss="sammon2",theta=1,
#' structures=c("cclusteredness","cdependence"),
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5)
#' ressam
#' summary(ressam)
#' plot(ressam)
#' plot(ressam,"Shepard")
#'
#' #STOPS with elastic 
#' ressam<-stops(dissm,loss="elastic",theta=1,
#' structures=c("cclusteredness","cdependence"),
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5)
#' ressam
#' summary(ressam)
#' plot(ressam)
#' plot(ressam,"transplot")
#'
#' #STOPS with sstress
#' resss<-stops(dissm,loss="sstress",theta=1,
#' structures=c("cclusteredness","cdependence"),
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5)
#' resss
#' summary(resss)
#' plot(resss)
#' plot(resss,"Shepard")
#' plot(resss,"transplot")
#'
#' #STOPS with powerstress
#' respstress<-stops(dissm,loss="powerstress",
#' structures=c("cclusteredness","cdependence"),theta=c(1,1,1),
#' strucpars=strucpars,weightmat=dissm,
#' itmaxps=1000,optimmethod="ALJ",lower=c(0.5,0.5,1),upper=c(5,5,5))
#' respstress
#' summary(respstress)
#' plot(respstress)
#'
#' #STOPS with restricted powerstress
#' respstressr<-stops(dissm,loss="powerstress",theta=c(1,1),
#' structures=c("cclusteredness","cdependence"),
#' strucpars=strucpars,weightmat=dissm,
#' itmaxps=1000,optimmethod="ALJ",lower=c(0.5,0.5),upper=c(5,5))
#' respstressr
#' summary(respstressr)
#' plot(respstressr)
#'
#' #STOPS with powermds
#' respmds<-stops(dissm,loss="powermds",
#' structures=c("cclusteredness","cdependence"),theta=c(1,1),
#' strucpars=strucpars,weightmat=dissm,
#' itmaxps=1000,optimmethod="ALJ",lower=c(0.5,0.5),upper=c(5,5))
#' respmds
#' summary(respmds)
#' plot(respmds)
#'
#' #STOPS with powersammon
#' respmds<-stops(dissm,loss="powersammon",theta=c(1,1),
#' structures=c("cclusteredness","cdependence"),
#' strucpars=strucpars,weightmat=dissm,
#' itmaxps=1000,optimmethod="ALJ",lower=c(0.5,0.5,1),upper=c(5,5,5))
#' respmds
#' summary(respmds)
#' plot(respmds)
#'
#' #STOPS with powerelastic
#' respmds<-stops(dissm,loss="powerelastic",theta=c(1,1),
#' structures=c("cclusteredness","cdependence"),
#' strucpars=strucpars,weightmat=dissm,
#' itmaxps=1000,optimmethod="ALJ",lower=c(0.5,0.5,1),upper=c(5,5,5))
#' respmds
#' summary(respmds)
#' plot(respmds)
#'
#' #STOPS with ordinal rstress 
#' resr<-stops(dissm,loss="rstress",type="ordinal",theta=1,
#' structures=c("cclusteredness","cdependence"),
#' strucpars=strucpars,weightmat=dissm,
#' itmaxps=1000,optimmethod="ALJ",lower=0.5,upper=5)
#' resr
#' summary(resr)
#' plot(resr)
#' 
#' #STOPS with approximated powerstress
#' respstressa<-stops(dissm,loss="powerstress",
#' structures=c("cclusteredness","cdependence"),
#' strucpars=strucpars,weightmat=dissm,theta=c(1,1,1),
#' itmaxps=1000,optimmethod="ALJ",lower=c(0.5,0.5,1),upper=c(5,5,5))
#' respstressa
#' summary(respstressa)
#' plot(respstressa,"transplot")
#' 
#' #STOPS with bcmds
#' resbcstress<-stops(dissm,loss="bcmds",
#' structures=c("cclusteredness","cdependence"),
#' strucpars=strucpars,optimmethod="ALJ",lower=c(0.5,1,0.5),upper=c(5,5,5))
#' resbcstress
#' summary(resbcstress)
#' plot(resbcstress)
#'
#' #STOPS with lmds
#' reslmds<-stops(dissm,loss="lmds",theta=c(1,1),
#' structures=c("cclusteredness","clinearity"),
#' strucpars=strucpars,optimmethod="ALJ",lower=c(2,0),upper=c(5,1))
#' reslmds
#' summary(reslmds)
#' plot(reslmds)
#'
#' #STOPS with Isomap (the k version)
#' resiso2<-stops(dissm,loss="isomap",theta=5,
#' structures=c("cclusteredness","clinearity"),
#' strucpars=strucpars,optimmethod="ALJ",lower=3,upper=10)
#' resiso2
#' summary(resiso2)
#' plot(resiso2)
#'
#' #STOPS with Isomap (the eps version)
#' resiso<-stops(dissm,loss="isomapeps",
#' structures=c("cclusteredness","clinearity"),theta=40,
#' strucpars=strucpars,optimmethod="ALJ",lower=50,upper=120)
#' resiso
#' summary(resiso)
#' plot(resiso)
#' 
#' }
#' @docType package
#' @name stops
NULL


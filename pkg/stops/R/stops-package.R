#' stops: Structure Optimized Proximity Scaling
#' 
#' A package for "structure optimized proximity scaling" (STOPS), a collection of methods that fit nonlinear distance transformations in multidimensional scaling (MDS) and trade-off the fit with structure considerations to find optimal parameters or optimal configurations. The package contains various functions, wrappers, methods and classes for fitting, plotting and displaying different MDS models in a STOPS framework like Torgerson scaling, SMACOF, Sammon mapping, elastic scaling, symmetric SMACOF, spherical SMACOF, sstress, rstress, powermds, power elastic scaling, power sammon mapping, power stress, Isomap, approximate power stress, restricted power stress. All of these models can also be fit as MDS variants (i.e., no structuredness). The package further contains functions for optimization (Adaptive Luus-Jaakola and for Bayesian optimization with treed Gaussian process with jump to linear models) and functions for various structuredness indices
#'
#' The stops package provides five categories of important functions:
#'
#' Models & Algorithms:
#' \itemize{
#' \item stops() ... which fits STOPS models as described in Rusch et al. (2023). By setting cordweight or strucweight to zero they can also be used to fit metric MDS for many different models, see below.  
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
#'Authors: Thomas Rusch, Patrick Mair, Kurt Hornik
#'
#'Maintainer: Thomas Rusch
#'
#'
#' @examples
#' #data(Swissroll)
#' #dissm<-as.matrix(dist(Swissroll[,1:3]))
#' #cols<-Swissroll[,4]
#' #structures<-c("cregularity","cdependence")
#' #strucweight<-c(-0.5,0.5)
#' #strucpars<-list(list(epsilon=10,minpts=2,scale=3),list(NULL))
#' 
#' dissm<-as.matrix(smacof::morse)
#'
#' #Setting up structurenedness parameters
#' strucpars<-list(list(epsilon=10,scale=3),list(NULL))
#' structures<-c("cclusteredness","cdependence")
#' 
#' #STOPS with strain
#' resstrain<-stops(dissm,loss="strain",theta=1,structures=structures,
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
#' structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5)
#' resstress
#' summary(resstress)
#' plot(resstress)
#' plot(resstress,"Shepard")
#'
#' #STOPS with smacofSphere
#' ressph<-stops(dissm,loss="smacofSphere",theta=1,
#' structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5)
#' ressph
#' summary(ressph)
#' plot(ressph)
#'
#' #STOPS with sammon
#' ressam<-stops(dissm,loss="sammon",theta=1,
#' structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5)
#' ressam
#' summary(ressam)
#' plot(ressam)
#' plot(ressam,"transplot")
#'
#' #STOPS with sammon2
#' ressam<-stops(dissm,loss="sammon2",theta=1,
#' structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5)
#' ressam
#' summary(ressam)
#' plot(ressam)
#' plot(ressam,"Shepard")
#'
#' #STOPS with elastic 
#' ressam<-stops(dissm,loss="elastic",theta=1,
#' structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5)
#' ressam
#' summary(ressam)
#' plot(ressam)
#' plot(ressam,"transplot")
#'
#' #STOPS with sstress
#' resss<-stops(dissm,loss="sstress",theta=1,
#' structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5)
#' resss
#' summary(resss)
#' plot(resss)
#' plot(resss,"Shepard")
#' plot(resss,"transplot")
#'
#' #STOPS with powerstress
#' respstress<-stops(dissm,loss="powerstress",
#' structures=structures,theta=c(1,1,1),
#' strucpars=strucpars,weightmat=dissm,
#' itmaxps=1000,optimmethod="ALJ",lower=c(0.5,0.5,1),upper=c(5,5,5))
#' respstress
#' summary(respstress)
#' plot(respstress)
#'
#' #STOPS with restricted powerstress
#' respstressr<-stops(dissm,loss="powerstress",theta=c(1,1),
#' structures=structures,
#' strucpars=strucpars,weightmat=dissm,
#' itmaxps=1000,optimmethod="ALJ",lower=c(0.5,0.5),upper=c(5,5))
#' respstressr
#' summary(respstressr)
#' plot(respstressr)
#'
#' #STOPS with powermds
#' respmds<-stops(dissm,loss="powermds",
#' structures=structures,theta=c(1,1),
#' strucpars=strucpars,weightmat=dissm,
#' itmaxps=1000,optimmethod="ALJ",lower=c(0.5,0.5),upper=c(5,5))
#' respmds
#' summary(respmds)
#' plot(respmds)
#'
#' #STOPS with powersammon
#' respmds<-stops(dissm,loss="powersammon",theta=c(1,1),
#' structures=structures,
#' strucpars=strucpars,weightmat=dissm,
#' itmaxps=1000,optimmethod="ALJ",lower=c(0.5,0.5,1),upper=c(5,5,5))
#' respmds
#' summary(respmds)
#' plot(respmds)
#'
#' #STOPS with powerelastic
#' respmds<-stops(dissm,loss="powerelastic",theta=c(1,1),
#' structures=structures,
#' strucpars=strucpars,weightmat=dissm,
#' itmaxps=1000,optimmethod="ALJ",lower=c(0.5,0.5,1),upper=c(5,5,5))
#' respmds
#' summary(respmds)
#' plot(respmds)
#'
#' #STOPS with ordinal rstress 
#' resr<-stops(dissm,loss="rstress",type="ordinal",theta=1,
#' structures=structures,
#' strucpars=strucpars,weightmat=dissm,
#' itmaxps=1000,optimmethod="ALJ",lower=0.5,upper=5)
#' resr
#' summary(resr)
#' plot(resr)
#' 
#' #STOPS with approximated powerstress
#' respstressa<-stops(dissm,loss="powerstress",
#' structures=structures,
#' strucpars=strucpars,weightmat=dissm,theta=c(1,1,1),
#' itmaxps=1000,optimmethod="ALJ",lower=c(0.5,0.5,1),upper=c(5,5,5))
#' respstressa
#' summary(respstressa)
#' plot(respstressa,"transplot")
#' 
#' #STOPS with bcmds
#' resbcstress<-stops(dissm,loss="bcmds",
#' structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=c(0.5,1,0.5),upper=c(5,5,5))
#' resbcstress
#' summary(resbcstress)
#' plot(resbcstress)
#'
#' #STOPS with lmds
#' reslmds<-stops(dissm,loss="lmds",theta=c(1,1),
#' structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=c(2,0),upper=c(5,1))
#' reslmds
#' summary(reslmds)
#' plot(reslmds)
#'
#' #STOPS with Isomap (the k version)
#' resiso2<-stops(dissm,loss="isomap",theta=5,
#' structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=3,upper=10)
#' resiso2
#' summary(resiso2)
#' plot(resiso2)
#'
#' #STOPS with Isomap (the eps version)
#' resiso<-stops(dissm,loss="isomapeps",
#' structures=structures,
#' theta=40,
#' strucpars=strucpars,optimmethod="ALJ",lower=50,upper=120)
#' resiso
#' summary(resiso)
#' plot(resiso)
#'
#' #STOPS with CLCA 
#' resclca<-stops(dissm,loss="clca",theta=0.3,
#' structures=structures, strucpars=strucpars,
#' strucweight=strucweight,lower=0.1,upper=5,
#' optimmethod="pso",itmax=20)
#' resclca
#' summary(resclca)
#' plot(resclca)
#'
#' #STOPS with pCLCA 
#' respclca<-stops(dissm,loss="pclca",theta=c(1,1,1,1),
#' structures=structures,strucpars=strucpars,
#' strucweight=strucweight,lower=c(0.1,0.1,0.1,0.1),upper=c(5,5,5,5),
#' optimmethod="ALJ",itmax=7)
#' respclca
#' summary(respclca)
#' plot(respclca)
#'
#' #STOPS with CLDA in k 
#' rescldak<-stops(dissm,loss="clda_k",theta=c(1,5),
#' structures=structures,strucpars=strucpars,
#' strucweight=strucweight,lower=c(0.2,4),upper=c(4,20),
#' optimmethod="pso",itmax=20)
#' rescldak
#' summary(rescldak)
#' plot(rescldak)
#'
#' #STOPS with CLDA in eps
#' rescldae<-stops(dissm,loss="clda_eps",theta=c(1,1),
#' structures=structures,strucpars=strucpars,
#' strucweight=strucweight,lower=c(0.2,0.5),upper=c(4,6),
#' optimmethod="SANN",itmax=20)
#' rescldae
#' summary(rescldae)
#' plot(rescldae)
#'
#' #STOPS with pCLDA with k (five parameters already..)
#' respcldak<-stops(dissm,loss="pclda_k",theta=c(1,1,1,1,5),
#' structures=structures,strucpars=strucpars,
#' strucweight=strucweight,lower=c(0.1,0.1,0.1,0.1,4),upper=c(5,5,5,5,20),
#' optimmethod="tgp",itmax=10)
#' respcldak
#' summary(respcldak)
#' plot(respcldak)
#'
#' #STOPS with pCLDA with eps (five parameter already..)
#' respcldae<-stops(dissm,loss="pclda_eps",theta=c(1,1,1,1,0.2),
#' structures=structures,strucpars=strucpars,
#' strucweight=strucweight,lower=c(0.1,0.1,0.1,0.1,0.1),upper=c(5,5,5,5,0.5),
#' optimmethod="tgp",itmax=20)
#' respcldae
#' summary(respcldae)
#' plot(respcldae)
#'
#' 
#' 
#' }
#' @docType package
#' @name stops
NULL


#' stops: Structure Optimized Proximity Scaling
#' 
#' A package for "structure optimized proximity scaling" (STOPS), a collection of methods that fit nonlinear distance transformations in multidimensional scaling (MDS) and trade-off the fit with structure considerations to find optimal parameters or optimal configurations. The package contains various functions, wrappers, methods and classes for fitting, plotting and displaying different MDS models in a STOPS framework like Torgerson scaling, SMACOF, Sammon mapping, elastic scaling, symmetric SMACOF, spherical SMACOF, sstress, rstress, powermds, power elastic scaling, power sammon mapping, power stress, Isomap, approximate power stress, restricted power stress. All of these models can also be fit as MDS variants (i.e., no structuredness). The package further contains functions for optimization (Adaptive Luus-Jaakola and for Bayesian optimization with treed Gaussian process with jump to linear models) and functions for various structuredness indices
#'
#' The stops package provides five categories of important functions:
#'
#' Models:
#' \itemize{
#' \item stops ... which fits STOPS models as described in Rusch et al. (2023) via argument loss with ratio, interval or ordinal optimal scaling (not all optimal scalings can be used with all loss arguments). By setting cordweight or strucweight to zero they can also be used to fit MDS for many different models:
#' \itemize{
#' \item loss="stress": One parameter theta, power transformations of dissimilarities. Via stop_smacofSym.
#' \item loss="sammon": One parameter theta, power transformations of dissimilarities. Via stop_sammon.
#' \item loss="strain": One parameter theta, power transformations of dissimilarities. Via stop_cmdscale.
#' \item loss="rstress": One parameter theta, power transformations of fitted distances. Via stop_rstress.
#' \item loss="smacofSphere": One parameter theta, power transformations of dissimilarities. Via stop_smacofSphere.
#' \item loss="sammon2": One parameter theta, power transformations of dissimilarities. Via stop_sammon2.
#' \item loss="elastic": One parameter theta, power transformations of dissimilarities. Via stop_elastic.
#' \item loss="sstress": One parameter theta, power transformations of dissimilarities. Via stop_sstress.
#' \item loss="isomap_k": One parameter theta, k neighborhood for geodesic distances. Via stop_isomap1.
#' \item loss="isomap_eps": One parameter theta, epsilon neighborhood for geodesic distances. Via stop_isomap2.
#' \item loss="smds": One parameter theta, neighborhood parameter tau. Via stop_smds.
#' \item loss="clca": One parameter theta, neighborhood parameters lambda0. Via stop_clca.
#' \item loss="powerelastic": Two parameter theta, power transformations of dissimilarities and fitted distances. Via stop_powerelastic.
#' \item loss="powersammon": Two parameter theta, power transformations of dissimilarities and fitted distances. Via stop_powersammon.
#' \item loss="powermds": Two parameter theta, power transformations of dissimilarities and fitted distances. Via stop_powermds.
#' \item loss="rpstress": Two parameter theta, power transformations of dissimilarities/fitted distances and weights. Via stop_rpowerstress.
#' \item loss="smdda_k": Two parameter theta, neighborhood parameters k (geodesic distance) and tau. Via stop_smmdak.
#' \item loss="smdda_eps": Two parameter theta, neighborhood parameters eps (geodesic distance) and tau. stop_smddae.
#' \item loss="lmds": Two parameter theta, neighbourhood parameters tau and k. Via stop_lmds.
#' \item loss="clda_eps": Two parameter theta, neighborhood parameters lambda0 and epsilon (geodesic distance). Via stop_cldae.
#' \item loss="clda_k": Two parameter theta, neighborhood parameters lambda0 and k (geodesic distance). Via stop_cldak.
#' \item loss="powerstress": Three parameter theta, power transformations of dissimilarities, fitted distances and weights. Via stop_powerstress.
#' \item loss="bcmds": Three parameter theta, Box-Cox transformations of dissimilarities, fitted distances and weights. Via stop_bcmds. 
#' \item loss="apstress": Three parameter theta, power transformations of dissimilarities, fitted distances and weights. Via stop_apstress.
#' \item loss="spmds": Four parameter theta, power transformations of dissimilarities, fitted distances and weights and neighborhood parametr tau. Via stop_spmds.
#' \item loss="spmdda_k": Five parameter theta, power transformations for dissimilarities, fitted distances and weights and neighborhood parameters k (geodesic distance) and tau. Via stop_spmddak.
#' \item loss="spmdda_eps": Five parameter theta, power transformations for dissimilarities, fitted distances and weights and neighborhood parameters eps (geodesic distance) and tau. Via stop_spmddae.
#' }
#' }
#'
#' Structuredness Indices:
#' Various c-structuredness as c_foo(), where foo is the name of the structuredness.  See Rusch et al. (2023). 
#' 
#' Optimization functions:
#' \itemize{
#' \item ljoptim() ... An (adaptive) version of the Luus-Jakola random search
#'}
#'
#' Methods: 
#' For most of the objects returned by the high-level functions S3 classes and methods for standard generics were implemented, including print, summary, plot, coef (extracting the hyperparameetrs).   
#'
#' References:
#' \itemize{
#' \item Rusch, T., Mair, P., & Hornik, K. (2023). Structure-based hyperparameter selection with Bayesian optimization in multidimensional scaling. Statistics & Computing, 33, [28]. https://doi.org/10.1007/s11222-022-10197-w
#' }
#' 
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
#' resstressm<-stops(dissm,loss="stress",theta=1,structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5,itmax=10,stoptype="multiplicative")
#' resstressm
#' plot(resstressm)
#'
#' \donttest{
#' #STOPS with stress or smacofSym
#' im<-2 #this is the itmax argument used for testing; use higher itmax in practice 
#' resstress<-stops(dissm,loss="stress",theta=1,
#' structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5,itmax=im)
#' resstress
#' summary(resstress)
#' plot(resstress)
#' plot(resstress,"Shepard")
#'
#' #STOPS with smacofSphere
#' ressph<-stops(dissm,loss="smacofSphere",theta=1,
#' structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5,itmax=im)
#' ressph
#' summary(ressph)
#' plot(ressph)
#'
#' #STOPS with sammon
#' ressam<-stops(dissm,loss="sammon",theta=1,
#' structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5,itmax=im)
#' ressam
#' summary(ressam)
#' plot(ressam)
#' plot(ressam,"transplot")
#'
#' #STOPS with sammon2
#' ressam<-stops(dissm,loss="sammon2",theta=1,
#' structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5,itmax=im)
#' ressam
#' summary(ressam)
#' plot(ressam)
#' plot(ressam,"Shepard")
#'
#' #STOPS with elastic 
#' ressam<-stops(dissm,loss="elastic",theta=1,
#' structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5,itmax=im)
#' ressam
#' summary(ressam)
#' plot(ressam)
#' plot(ressam,"transplot")
#'
#' #STOPS with sstress
#' resss<-stops(dissm,loss="sstress",theta=1,
#' structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=0.5,upper=5,itmax=im)
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
#' itmaxps=1000,optimmethod="ALJ",lower=c(0.5,0.5,1),upper=c(5,5,5),itmax=im)
#' respstress
#' summary(respstress)
#' plot(respstress)
#'
#' #STOPS with restricted powerstress
#' respstressr<-stops(dissm,loss="powerstress",theta=c(1,1),
#' structures=structures,
#' strucpars=strucpars,weightmat=dissm,
#' itmaxps=1000,optimmethod="ALJ",lower=c(0.5,0.5),upper=c(5,5),itmax=im)
#' respstressr
#' summary(respstressr)
#' plot(respstressr)
#'
#' #STOPS with powermds
#' respmds<-stops(dissm,loss="powermds",
#' structures=structures,theta=c(1,1),
#' strucpars=strucpars,weightmat=dissm,
#' itmaxps=1000,optimmethod="ALJ",lower=c(0.5,0.5),upper=c(5,5),itmax=im)
#' respmds
#' summary(respmds)
#' plot(respmds)
#'
#' #STOPS with powersammon
#' respmds<-stops(dissm,loss="powersammon",theta=c(1,1),
#' structures=structures,
#' strucpars=strucpars,weightmat=dissm,
#' itmaxps=1000,optimmethod="ALJ",lower=c(0.5,0.5),upper=c(5,5),itmax=im)
#' respmds
#' summary(respmds)
#' plot(respmds)
#'
#' #STOPS with powerelastic
#' respmds<-stops(dissm,loss="powerelastic",theta=c(1,1),
#' structures=structures,
#' strucpars=strucpars,weightmat=dissm,
#' itmaxps=1000,optimmethod="ALJ",lower=c(0.5,0.5,1),upper=c(5,5,5),itmax=im)
#' respmds
#' summary(respmds)
#' plot(respmds)
#'
#' #STOPS with ordinal rstress 
#' resr<-stops(dissm,loss="rstress",type="ordinal",theta=1,
#' structures=structures,
#' strucpars=strucpars,weightmat=dissm,
#' itmaxps=1000,optimmethod="ALJ",lower=0.5,upper=5,itmax=im)
#' resr
#' summary(resr)
#' plot(resr)
#' 
#' #STOPS with approximated powerstress
#' respstressa<-stops(dissm,loss="powerstress",
#' structures=structures,
#' strucpars=strucpars,weightmat=dissm,theta=c(1,1,1),
#' itmaxps=1000,optimmethod="ALJ",lower=c(0.5,0.5,1),upper=c(5,5,5),itmax=im)
#' respstressa
#' summary(respstressa)
#' plot(respstressa,"transplot")
#' 
#' #STOPS with bcmds
#' resbcstress<-stops(dissm,loss="bcmds",
#' structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=c(0.5,1,0.5),upper=c(5,5,5),itmax=im)
#' resbcstress
#' summary(resbcstress)
#' plot(resbcstress)
#'
#' #STOPS with lmds
#' reslmds<-stops(dissm,loss="lmds",theta=c(1,1),
#' structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=c(2,0),upper=c(5,1),itmax=im)
#' reslmds
#' summary(reslmds)
#' plot(reslmds)
#'
#' #STOPS with Isomap (the k version)
#' resiso2<-stops(dissm,loss="isomap",theta=5,
#' structures=structures,
#' strucpars=strucpars,optimmethod="ALJ",lower=3,upper=10,itmax=im)
#' resiso2
#' summary(resiso2)
#' plot(resiso2)
#'
#' #STOPS with Isomap (the eps version)
#' resiso<-stops(dissm,loss="isomapeps",
#' structures=structures,
#' theta=40,
#' strucpars=strucpars,optimmethod="ALJ",lower=50,upper=120,itmax=im)
#' resiso
#' summary(resiso)
#' plot(resiso)
#'
#' strucweight<-c(-0.5,-0.5)
#' 
#' #STOPS with smds
#' resclca<-stops(dissm,loss="smds",theta=0.3,
#' structures=structures, strucpars=strucpars,
#' strucweight=strucweight,lower=0.1,upper=5,
#' optimmethod="pso",itmax=im*4)
#' resclca
#' summary(resclca)
#' plot(resclca)
#'
#' #STOPS with spmds
#' respclca<-stops(dissm,loss="spmds",theta=c(1,1,1,1),
#' structures=structures,strucpars=strucpars,
#' strucweight=strucweight,lower=c(0.1,0.1,0.1,0.1),upper=c(5,5,5,5),
#' optimmethod="ALJ",itmax=im)
#' respclca
#' coef(respclca)
#' summary(respclca)
#' plot(respclca)
#'
#' #STOPS with smdda and k 
#' rescldak<-stops(dissm,loss="smdda_k",theta=c(1,5),
#' structures=structures,strucpars=strucpars,
#' strucweight=strucweight,lower=c(0.2,4),upper=c(4,20),
#' optimmethod="pso",itmax=im*4)
#' rescldak
#' summary(rescldak)
#' plot(rescldak)
#'
#' #STOPS with smdda in eps
#' set.seed(123)
#' rescldae<-stops(dissm,loss="smdda_eps",theta=c(1,2),
#' structures=structures,strucpars=strucpars,
#' strucweight=strucweight,lower=c(0.2,1),upper=c(4,10),
#' optimmethod="SANN",itmax=10*im,stoptype="multiplicative")
#' rescldae
#' summary(rescldae)
#' plot(rescldae)
#'
#' #STOPS with spmdda with k (five parameters already..)
#' respcldak<-stops(dissm,loss="spmdda_k",theta=c(1,1,1,1,5),
#' structures=structures,strucpars=strucpars,
#' strucweight=strucweight,lower=c(0.8,0.8,0.8,0.8,4),upper=c(5,5,5,5,20),
#' optimmethod="tgp",itmax=im)
#' respcldak
#' summary(respcldak)
#' plot(respcldak)
#'
#' #STOPS with spmdda with eps (five parameter already..)
#' set.seed(123)
#' respcldae<-stops(dissm,loss="spmdda_eps",theta=c(1,1,1,1,2),
#' structures=structures,strucpars=strucpars,
#' strucweight=strucweight,lower=c(0.8,0.8,0.8,0.8,1),upper=c(5,5,5,5,10),
#' optimmethod="tgp",itmax=im)
#' respcldae
#' summary(respcldae)
#' plot(respcldae)
#'
#' #STOPS with clca 
#' set.seed(123)
#' resclca<-stops(dissm,loss="clca",theta=c(1),
#' structures=structures,strucpars=strucpars,
#' strucweight=strucweight,lower=0.1,upper=5,itmax=im)
#' resclca
#' summary(resclca)
#' plot(resclca)
#'
#'
#' }
#' @keywords internal
#' @aliases stops-package
"_PACKAGE" 
NULL


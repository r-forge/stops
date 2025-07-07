
#' High Level STOPS Function
#'
#' @description This allows to fit STOPS models as described in Rusch, Mair, Hornik (2023). 
#'
#' @details The combination of c-structurednes indices and stress uses the stress.m values, which are the explictly normalized stresses. Reported however is the stress-1 value which is sqrt(stress.m). 
#' 
#' @param dis numeric matrix or dist object of a matrix of proximities
#' @param loss which loss function to be used for fitting, defaults to stress. 
#' @param theta hyperparameter vector starting values for the transformation functions. If the length is smaller than the number of hyperparameters for the MDS version the vector gets recycled (see the corresponding stop_XXX function or the vignette for how theta must look like exactly for each loss). If larger than the number of hyperparameters for the MDS method, an error is thrown. If completely missing theta is set to 1 and recycled.
#' @param type type of MDS optimal scaling (implicit transformation). One of "ratio", "interval", "mspline" or "ordinal". Default is "ratio". Not every type can be used with every loss, only ratio works with all.
#' @param structures character vector of which c-structuredness indices should be considered; if missing no structure is considered.
#' @param ndim number of dimensions of the target space
#' @param weightmat (optional) a matrix of nonnegative weights; defaults to 1 for all off diagonals 
#' @param init (optional) initial configuration
#' @param stressweight weight to be used for the fit measure; defaults to 1
#' @param strucweight vector of weights to be used for the c-structuredness indices (in the same order as in structures); defaults to -1/length(structures) for each index 
#' @param strucpars (possibly named with the structure). Metaparameters for the structuredness indices (gamma in the article). It's safest for it be a list of lists with the named arguments for the structuredness indices and the order of the lists must be like the order of structures. So something like this \code{list(list(par1Struc1=par1Struc1,par2Struc1=par2Struc1),list(par1Struc2=par1Struc2,par2Struc2=par2Struc2),...)} where parYStrucX are the named arguments for the metaparameter Y of the structure X the list elements corresponds to. For a structure without parameters, set NULL. Parameters in different list elements parYStrucX can have the same name. For example, say we want to use cclusteredness with metaparameters epsilon=10 and k=4 (and the default for the other parameters), cdependence with no metaparameters and cfaithfulness with metaparameter k=7 one would \code{list(list(epsilon=10,k=4),list(NULL),list(dis=obdiss,k=6))}  for structures vector ("cclusteredness","cdependence","cfaithfulness"). The parameter lists must be in the same ordering as the indices in structures. If missing it is set to NULL and defaults are used. It is also possible to supply a structure's metaparameters as a list of vectors with named elements if the metaparameters are scalars, so like \code{list(c(par1Struc1=parStruc1,par2Struc1=par1Struc1,...),c(par1Struc2=par1Struc2,par2Struc2=par2Struc2,...))}. That can have unintended consequences if the metaparameter is a vector or matrix.  
#' @param optimmethod What solver to use. Currently supported are Bayesian optimization with Gaussian Process priors and Kriging ("Kriging", for which the archived package 'DiceOptim' must be installed), Bayesian optimization with treed Gaussian processes with jump to linear models ("tgp", see \code{\link[tgp]{dopt.gp}}), Adaptive LJ Search ("ALJ"), Particle Swarm optimization ("pso", see \code{\link[pso]{psoptim}}), simulated annealing ("SANN", \code{\link[stats]{optim}}), "direct (\code{\link[nloptr]{direct}})", Stochastic Global Optimization ("stogo", \code{\link[nloptr]{stogo}}), COBYLA ("cobyla", \code{\link[nloptr]{cobyla}}), Controlled Random Search 2 with local mutation ("crs2lm", \code{\link[nloptr]{crs2lm}}), Improved Stochastic Ranking Evolution Strategy ("isres", \code{\link[nloptr]{isres}}), Multi-Level Single-Linkage ("mlsl", \code{\link[nloptr]{mlsl}}), Nelder-Mead ("neldermead", \code{\link[nloptr]{neldermead}}), Subplex ("sbplx", \code{\link[nloptr]{sbplx}}), Hooke-Jeeves Pattern Search ("hjk", \code{\link[dfoptim]{hjk}}), CMA-ES ("cmaes", \code{\link[cmaes]{cma_es}}). Defaults to "ALJ" version. "tgp", "ALJ", "Kriging" and "pso" usually work well for relatively low values of 'itmax'. 
#' @param lower The lower contraints of the search region. Needs to be a numeric vector of the same length as the parameter vector theta. 
#' @param upper The upper contraints of the search region. Needs to be a numeric vector of the same length as the parameter vector theta.  
#' @param verbose numeric value hat prints information on the fitting process; >2 is very verbose.
#' @param stoptype which aggregation for the multi objective target function? Either 'additive' (default) or 'multiplicative'
#' @param itmax maximum number of iterations of the outer optimization (for theta) or number of steps of Bayesian optimization; default is 50. We recommend a higher number for ALJ (around 150). Note that due to the inner workings of some solvers, this may or may not correspond to the actual number of function evaluations performed (or PS models fitted). E.g., with tgp the actual number of function evaluation of the PS method is between itmax and 6*itmax as tgp samples 1-6 candidates from the posterior and uses the best candidate. For pso it is the number of particles s times itmax. For cmaes it is usually a bit higher than itmax. This currently may get overruled by a control argument if it is used (and then set to either ewhat is supplie dby control or to the default of the method).    
#' @param itmaxps maximum number of iterations of the inner optimization (to obtain the PS configuration)
#' @param accps accuracy of the inner optimization (to obtain the PS configuration)
#' @param initpoints number of initial points to fit the surrogate model for Bayesian optimization; default is 10.
#' @param model a character specifying the surrogate model to use. For "Kriging" it specifies the covariance kernel for the GP prior; see \code{\link[DiceKriging]{covTensorProduct-class}} defaults to "powerexp". For "tgp" it specifies the non stationary process used see \code{\link[tgp]{bgp}}, defaults to "btgpllm" 
#' @param control a control argument passed to the outer optimization procedure. Will override any other control arguents passed, especially verbose and itmax. For the effect of control, see the functions pomp::sannbox for SANN and pso::psoptim for pso, cmaes::cma_es for cmaes, dfoptim::hjkb for hjk and the nloptr docs for the algorithms direct, stogo, cobyla, crs2lm, isres, mlsl, neldermead, sbplx.
#' @param registry an object of class registry containing the c-structuredness indices. Defaults to the what is created .onLoad. 
#' @param ... additional arguments passed to the outer optimization procedures (not fully tested).
#' 
#
#'@return A list with the components
#'         \itemize{
#'         \item stoploss: the stoploss value
#'         \item optim: the object returned from the optimization procedure
#'         \item stressweight: the stressweight
#'         \item strucweight: the vector of structure weights
#'         \item call: the call
#'         \item optimmethod: The solver selected
#'         \item loss: The PS badness-of-fit function
#'         \item nobj: the number of objects in the configuration
#'         \item type: The type of stoploss scalacrisation (additive or multiplicative)
#'         \item fit: The fitted PS object (most importantly $fit$conf the fitted configuration)
#'          \item stoptype: Type of stoploss combinatio
#'    
#' }
#' 
#' @examples
#' data(kinshipdelta,package="smacof")
#' strucpar<-list(NULL,NULL) #parameters for indices
#' res1<-stops(kinshipdelta,loss="stress",
#' structures=c("cclumpiness","cassociation"),strucpars=strucpar,
#' lower=0,upper=10,itmax=10)
#' res1
#' 
#' \donttest{
#' #use higher itmax in general, we use 5 just to shorten the tests
#' data(BankingCrisesDistances)
#' strucpar<-list(c(epsilon=10,minpts=2),NULL) #parameters for indices
#' res1<-stops(BankingCrisesDistances[,1:69],loss="stress",verbose=0,
#' structures=c("cclusteredness","clinearity"),strucpars=strucpar,
#' lower=0,upper=10,itmax=5)
#' res1
#'
#' strucpar<-list(list(alpha=0.6,C=15,var.thr=1e-5,zeta=NULL),
#' list(alpha=0.6,C=15,var.thr=1e-5,zeta=NULL))
#' res1<-stops(BankingCrisesDistances[,1:69],loss="stress",verbose=0,
#' structures=c("cfunctionality","ccomplexity"),strucpars=strucpar,
#' lower=0,upper=10,itmax=5)
#' res1
#' }
#' 
#' @importFrom stats dist as.dist optim
#' @importFrom utils capture.output
#' @importFrom pso psoptim
# @importFrom DiceOptim EGO.nsteps
# @importFrom DiceKriging km
#' @importFrom tgp lhs dopt.gp
#' @importFrom pomp sannbox
#' @importFrom nloptr direct stogo cobyla crs2lm isres mlsl neldermead sbplx
#' @importFrom cmaes cma_es
#' @importFrom dfoptim hjkb
#' @import cordillera
#' 
#' @keywords clustering multivariate
#' @export
stops <- function(dis,loss="stress", theta=1, type="ratio",structures, ndim=2, weightmat=NULL, init=NULL, stressweight=1, strucweight, strucpars, optimmethod=c("SANN","ALJ","pso","Kriging","tgp","direct","stogo","cobyla","crs2lm","isres","mlsl","neldermead","sbplx","hjk","cmaes"), lower, upper, verbose=0, stoptype=c("additive","multiplicative"), initpoints=10, itmax=50,itmaxps=10000, accps=1e-8, model, control,registry=struc_reg,...)
    {
      if(missing(structures)) {
          structures <- "clinearity"
          strucweight <- 0
      }
      ## allowed losses
      loss <- match.arg(loss,c("strain","stress","smacofSym","powerstress","powermds","powerelastic","powerstrain","elastic","sammon","sammon2","smacofSphere","powersammon","rstress","sstress","isomap","isomapeps","isomap_eps","isomap_k","bcstress","bcmds","lmds","apstress","rpowerstress","smds","spmds","smdda_k","smdda_eps","spmdda_k","spmdda_eps","clca","clda_k","clda_eps"),several.ok=FALSE)
      ## allowed structures
      ##structures <- match.arg(structures,c("cclusteredness","clinearity","cdependence","cmanifoldness","cassociation","cnonmonotonicity","cfunctionality","ccomplexity","cfaithfulness","cregularity","chierarchy","cconvexity","cstriatedness","coutlying","cskinniness","csparsity","cstringiness","cclumpiness","cinequality","cshepardness"),several.ok=TRUE)
      ##Write the above to match against the registry names
     
      if(missing(strucpars)) strucpars <- vector("list",length(structures))
      if(inherits(dis,"dist")) dis <- as.matrix(dis)
      if(is.null(weightmat)) weightmat <- 1-diag(dim(dis)[1])
      #if(missing(loss)) loss <- "stress"
      if(missing(stoptype)) stoptype <- "additive"
      #TODO implement a Pareto multiobjective
      .confin <- init #initialize a configuration

      psfunc <- switch(loss, "powerstrain"=stop_cmdscale, "stress"=stop_smacofSym,"smacofSym"=stop_smacofSym,"powerstress"=stop_powerstress,"strain"=stop_cmdscale,"smacofSphere"=stop_smacofSphere,"rstress"=stop_rstress,"sammon"=stop_sammon, "elastic"=stop_elastic, "powermds"=stop_powermds,"powerelastic"=stop_powerelastic,"powersammon"=stop_powersammon,"sammon2"=stop_sammon2,"sstress"=stop_sstress,"isomap"=stop_isomap1,"isomap_k"=stop_isomap1,"isomapeps"=stop_isomap2,"isomap_eps"=stop_isomap2,"bcstress"=stop_bcmds,"bcmds"=stop_bcmds,"lmds"=stop_lmds,"apstress"=stop_apstress,"rpowerstress"=stop_rpowerstress,"smds"=stop_smds,"spmds"=stop_spmds,"smdda_k"=stop_smddak,"smdda_eps"=stop_smddae,"spmdda_k"=stop_spmddak,"spmdda_eps"=stop_spmddae,"clca"=stop_clca,"clda_k"=stop_cldak,"clda_eps"=stop_cldae) #choose the stress to minimize
      if(missing(strucweight)) {
         #TODO: automatic handler of setting weights that makes sense
         strucweight <- rep(-1/length(structures),length(structures))
         if(verbose>1) cat("Weights are stressweight=",stressweight,"strucweights=", strucweight,"\n")
      }
      if(missing(optimmethod)) optimmethod <- "ALJ"  
      if(verbose>0) cat("Starting Optimization \n ")
        if(optimmethod=="SANN") {
            if(missing(control)) control <- list(trace=verbose-2,lower=lower,upper=upper,maxit=itmax)
            opt <- pomp::sannbox(par=theta,fn=function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps,stoptype=stoptype,registry=registry,acc=accps))$stoploss,control=control)
          thetaopt <- opt$par
          bestval <- opt$value
          itel <- opt$counts[1]
      }
       if(optimmethod=="pso") {
        #addargs <- list(...)
        if(missing(control)) control <- list(trace=verbose-2,s=5,maxit=itmax)
        opt<- pso::psoptim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps,stoptype=stoptype,registry=registry,acc=accps))$stoploss,lower=lower,upper=upper,control=control,...)
         thetaopt <- opt$par
         bestval <-  opt$value
         itel <- opt$counts["function"]
       }
      if(optimmethod=="ALJ")  {
        opt <- ljoptim(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps,stoptype=stoptype,registry=registry,acc=accps))$stoploss,lower=lower,upper=upper,verbose=verbose-2,itmax=itmax,...)
       thetaopt <- opt$par
       bestval <-  opt$value
       itel <- opt$counts["function"] 
    }    
    if(optimmethod=="Kriging")
    {
        if (!requireNamespace("DiceOptim", quietly = TRUE)) stop("Package \"DiceOptim\" must be installed to use this optimazation method. I suggest you use 'optimmethod=tgp'.", call. = FALSE)
        if(missing(model)) model <- "powexp"
       # optdim <- 3 #dimensions
       # if(loss%in%c("powerstrain","stress","smacofSym","smacofSphere","strain","sammon","elastic","sammon2","sstress","rstress")) optdim <- 1
       # if(loss%in%c("powermds","powerelastic","powersammon","smacofSphere","strain","sammon","elastic","sammon2")) optdim <- 2
        recto <- cbind(lower,upper)
        Xcand <- tgp::lhs(initpoints*100,recto)
        if(missing(theta)) theta <- Xcand[1,]
        x <- t(theta)
        X <- tgp::dopt.gp(initpoints-1,X=x,Xcand)$XX
        design <- rbind(x,X)
        #design <- data.frame(X) 
        responsec <- apply(design, 1, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps,stoptype=stoptype,registry=registry,acc=accps))$stoploss) #support points for fitting kriging model
        if (verbose>1) cat("Kriging Model Fitting","\n")
        surrogatemodel <- DiceKriging::km(~1, design = design, response = responsec,covtype=model,control=list(trace=isTRUE(verbose>3))) #fit the kriging model
        #EGO.nsteps has no verbose argument so I capture.output and return it if desired
        if (verbose>2) cat("EGO (DICE) Optimization","\n")
        logged <- capture.output({
           opt<- DiceOptim::EGO.nsteps(model=surrogatemodel, fun=function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps,stoptype=stoptype,registry=registry,acc=accps))$stoploss,lower=lower,upper=upper,nsteps=itmax,...)
       }) #bayesian optimization with gaussian process prior
       if(verbose>2) print(logged)
       thetaopt <- opt$par[which.min(opt$value),] #parameters where best value found (we do not use the last one as that may be worse)
       bestval <- min(opt$value) #best stoploss value
       itel <- itmax
       }
  if(optimmethod=="tgp")
    {
        #if(!isNamespaceLoaded("tgp")) attachNamespace("tgp")
        if(missing(model)) model <- "btgpllm"
        #model <- get(model,envir=getNamespace("tgp"))
        #if(loss%in%c("powerstrain","stress","smacofSym","smacofSphere","strain","sammon","elastic","sammon2","sstress","rstress")) optdim <- 1
        #if(loss%in%c("powermds","powerelastic","powersammon","smacofSphere","strain","sammon","elastic","sammon2")) optdim <- 2
        if (verbose>1) cat("EGO (TGP) Optimization","\n")
        opt <- tgpoptim(theta, fun=function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps,stoptype=stoptype,registry=registry,acc=accps))$stoploss,lower=lower,upper=upper,itmax=itmax,initpoints=initpoints,model=model,verbose=verbose-2,...) #bayesian optimization with treed gaussian process prior
       thetaopt <- opt$par #parameters where best value found (we do not use the last one as that may be worse)
       bestval <- opt$value #best stoploss value
       itel <- opt$counts["function"] 
    }
     if(optimmethod=="direct") {
            if (verbose>1) cat("DIRECT Optimization","\n")
            if(missing(control)) control <- list(maxeval=itmax)
          opt<- nloptr::direct(function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps,stoptype=stoptype,registry=registry,acc=accps))$stoploss,lower=lower,upper=upper,nl.info=isTRUE(all.equal(verbose-2,0)),control=control,...)
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$iter
     }
      if(optimmethod=="stogo") {
          if (verbose>1) cat("StoGO Optimization","\n")
          #cat(itmax,"\n")
         #   if(missing(control)) control <- list(maxeval=itmax)  
          opt<- nloptr::stogo(theta,function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps,stoptype=stoptype,registry=registry,acc=accps))$stoploss,lower=lower,upper=upper,nl.info=isTRUE(verbose>2),maxeval=itmax,...)
          #TODO Issue with maxeval?
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$iter
      }
         if(optimmethod=="cobyla") {
            if (verbose>1) cat("COBYLA Optimization","\n")
            if(missing(control)) control <- list(maxeval=itmax)  
           opt<- nloptr::cobyla(theta,function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps,stoptype=stoptype,registry=registry,acc=accps))$stoploss,lower=lower,upper=upper,nl.info=isTRUE(all.equal(verbose-2,0)),control=control,...)
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$iter
         }
        if(optimmethod=="crs2lm") {
            if (verbose>1) cat("crs2lm Optimization","\n")
            #if(missing(control)) control <- list(maxeval=itmax)  
           opt<- nloptr::crs2lm(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps,stoptype=stoptype,registry=registry,acc=accps))$stoploss,lower=lower,upper=upper,nl.info=isTRUE(verbose>3),maxeval=itmax,...)
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$iter
        }
        if(optimmethod=="isres") {
            if (verbose>1) cat("isres Optimization","\n")
            #if(missing(control)) control <- list(maxeval=itmax)  
           opt<- nloptr::crs2lm(theta,function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps,stoptype=stoptype,registry=registry,acc=accps))$stoploss,lower=lower,upper=upper,nl.info=isTRUE(verbose>3),maxeval=itmax,...)
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$iter
        }
        if(optimmethod=="mlsl") {
            if (verbose>1) cat("MLSL Optimization","\n")
            if(missing(control)) control <- list(maxeval=itmax)  
           opt<- nloptr::mlsl(theta,function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps,stoptype=stoptype,registry=registry,acc=accps))$stoploss,lower=lower,upper=upper,nl.info=isTRUE(verbose>3),control=control,...)
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$iter
        }
        if(optimmethod=="neldermead") {
            if (verbose>1) cat("Nelder-Mead Optimization","\n")
            if(missing(control)) control <- list(maxeval=itmax)  
           opt<- nloptr::neldermead(theta,function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps,stoptype=stoptype,registry=registry,acc=accps))$stoploss,lower=lower,upper=upper,nl.info=isTRUE(verbose>3),control=control,...)
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$iter
        }
          if(optimmethod=="sbplx") {
            if (verbose>1) cat("Subplex Optimization","\n")
            if(missing(control)) control <- list(maxeval=itmax)  
           opt<- nloptr::sbplx(theta,function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps,stoptype=stoptype,registry=registry))$stoploss,lower=lower,upper=upper,nl.info=isTRUE(verbose>3),control=control,...)
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$iter
      }
        if(optimmethod=="hjk") {
            if (verbose>1) cat("Hooke-Jeeves Optimization","\n")
            if(missing(control)) control <- list(info=isTRUE(all.equal(verbose-2,0)),maxfeval=itmax)
          opt<- dfoptim::hjkb(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps,stoptype=stoptype,registry=registry,acc=accps))$stoploss,lower=lower,upper=upper,control=control,...)
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$feval
        }
         if(optimmethod=="cmaes") {
            if (verbose>1) cat("CMA-ES Optimization","\n")
            theta <- as.vector(theta)
            if(missing(control)) {
                  #for calculating the correct number of calls to fn
                  N <- length(theta)
                  lambda <- 4 + floor(3 * log(N))
                  control <- list(maxit=ceiling(itmax/lambda))
                }
          opt<- cmaes::cma_es(theta, function(theta) do.call(psfunc,list(dis=dis,theta=theta,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps,stoptype=stoptype,registry=registry,acc=accps))$stoploss,lower=lower,upper=upper,control=control,...)
            thetaopt <- opt$par
            bestval <- opt$value
            itel <- opt$counts[1]
         }
        #TODO: Streamline the number of function evaluations to be supplied returned for the different solvers. E.g. for cma_es it is itmax*population size. Or for tgp it is also itmax*6 or so. 
    #refit optimal model  
    out <- do.call(psfunc,list(dis=dis,theta=thetaopt,ndim=ndim,weightmat=weightmat,init=.confin,structures=structures,stressweight=stressweight,strucweight=strucweight,strucpars=strucpars,verbose=verbose-3,type=type,itmaxi=itmaxps,stoptype=stoptype,registry=registry,acc=accps))
    out$stoploss <- bestval
    out$theta <- out$parameters
    out$optim <- opt
    out$stressweight <- stressweight
    out$strucweight <- strucweight
    out$call <- match.call()
    out$optimethod <- optimmethod
    out$loss <- loss
    out$nobj <- dim(out$fit$conf)[1]
    out$type <- type
    out$stoptype  <- stoptype    
    if(verbose>1) cat("Found minimum after",itel," iterations at",round(thetaopt,4),"with stoploss=",round(out$stoploss,4),"and default scaling loss=",round(out$stress.m,4),"and c-structuredness indices:",t(data.frame(names(out$strucindices),out$strucindices)),". Thanks for your patience. \n")
    class(out) <- c("stops")
    out
  }


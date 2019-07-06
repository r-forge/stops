
library(stops)
data(kinshipdelta)


diss <- as.matrix(kinshipdelta)
set.seed(210485)

#############################
# Test all losses dispatch
structures <- c("clinearity")
optstrain <- stops(diss,loss="strain",verbose=3,structures=structures,upper=5,lower=0) #works
optstress <- stops(diss,loss="stress",verbose=3,structures=structures,upper=5,lower=1) #not works
optsmacofSym <- stops(diss,loss="smacofSym",verbose=3,structures=structures,upper=5,lower=0) #not works
optsmacofSph <- stops(diss,loss="smacofSphere",verbose=3,structures=structures,upper=5,lower=0) #not works
optsammon <- stops(diss,loss="sammon",verbose=3,structures=structures,upper=5,lower=0) #works
optsammon2 <- stops(diss,loss="sammon2",verbose=3,structures=structures,upper=5,lower=0) #not works
optelastic <- stops(diss,loss="elastic",verbose=3,structures=structures,upper=5,lower=0) #not works
optsstress <- stops(diss,loss="sstress",verbose=3,structures=structures,upper=5,lower=0) #works
optrstress <- stops(diss,loss="rstress",verbose=3,structures=structures,upper=5,lower=0) #works
optpowermds <- stops(diss,loss="powermds",verbose=3,structures=structures,upper=c(5,5),lower=c(0,0)) #works
optpowerstress <- stops(diss,loss="powerstress",verbose=3,structures=structures,upper=c(5,5,5),lower=c(0,0,0)) #works
optpowersammon <- stops(diss,loss="powersammon",verbose=3,structures=structures,upper=c(5,5,5),lower=c(0,0,0)) #works
optpowerelastic <- stops(diss,loss="powerelastic",verbose=3,structures=structures,upper=c(5,5,5),lower=c(0,0,0)) #works
optpowerstrain <- stops(diss,loss="powerstrain",verbose=3,structures=structures,upper=5,lower=0) #works
optisomapk <- stops(diss,loss="isomap",verbose=3,structures=structures,upper=10,lower=2) #works
optisomapeps <- stops(diss,loss="isomapeps",verbose=3,structures=structures,lower=50,upper=100) #works
optbc <- stops(diss,loss="bcstress",verbose=3,structures=structures,lower=c(0,0,0),upper=c(5,5,5)) #works
optlmds <- stops(diss,loss="lmds",theta=1,verbose=3,structures=structures,lower=c(1,0),upper=c(20,2)) #works

# print method
optstrain
optstress
optsmacofSym
optsammon
optsammon2
optelastic
optsstress
optrstress
optpowermds
optpowerstress
optpowersammon
optpowerelastic
optpowerstrain
optelastic
optisomapk
optisomapeps
optbc
optlmds

# plot method
plot(optstrain)
plot(optstress)
plot(optsmacofSym)
plot(optsammon)
plot(optsammon2)
plot(optelastic)
plot(optsstress)
plot(optrstress)
plot(optpowermds)
plot(optpowerstress)
plot(optpowersammon)
plot(optpowerelastic)
plot(optpowerstrain)
plot(optelastic)
plot(optisomapk)
plot(optisomapeps)
plot(optbc)
plot(optlmds)


###########
#test stoploss combination
verbose <- 3
stressweight <- 1
structures <- c("clinearity","cdependence")
type <- c("additive")
s1a <- stoploss(diss,stressweight=stressweight,structures=structures,strucweight=strucweight,type="additive",verbose=verbose)
s1a
#test stoploss with clinearity multiplicative
s1m <- stoploss(diss,stressweight=stressweight,structures=structures,strucweight=strucweight,type="multiplicative",verbose=verbose)
s1m


#########
# Test optimizers
sannres <- stoploss(diss,stressweight=stressweight,structures=structures,strucweight=strucweight,type="additive",verbose=verbose,optimmethod="SANN")
aljres <- stoploss(diss,stressweight=stressweight,structures=structures,strucweight=strucweight,type="additive",verbose=verbose,optimmethod="ALJ")
psores <- stoploss(diss,stressweight=stressweight,structures=structures,strucweight=strucweight,type="additive",verbose=verbose,optimmethod="pso")
krigres <- stoploss(diss,stressweight=stressweight,structures=structures,strucweight=strucweight,type="additive",verbose=verbose,optimmethod="Kriging")
tgpres <- stoploss(diss,stressweight=stressweight,structures=structures,strucweight=strucweight,type="additive",verbose=verbose,optimmethod="tgp")

#########
# Test different weights
structures <- c("clinearity","cdependence")
halfres <- stoploss(diss,stressweight=1,structures=structures,strucweight=c(-0.5,-0.5),verbose=verbose,optimmethod="ALJ")
fullres <- stoploss(diss,stressweight=0.5,structures=structures,strucweight=c(-1,-10),verbose=verbose,optimmethod="ALJ")
all.equal(halfres$conf,fullres$conf)


######
# Test  different theta lengths
structures <- c("clinearity","cdependence")
#theta scalar
thetas <- stoploss(diss,stressweight=1,loss="powerstress",theta=1,structures=structures,strucweight=c(-0.5,-0.5),verbose=verbose,optimmethod="ALJ",upper=5,lower=0)
#theta vector of 2
thetav <- stoploss(diss,stressweight=1,loss="powerstress",theta=c(2,2),structures=structures,strucweight=c(-0.5,-0.5),verbose=verbose,optimmethod="ALJ",upper=5,lower=0)
#theta vector of 3
thetav <- stoploss(diss,stressweight=1,loss="powerstress",theta=c(2,2,4),structures=structures,strucweight=c(-0.5,-0.5),verbose=verbose,optimmethod="ALJ",upper=5,lower=0)





#test stop_smacofSym
sres <- stop_smacofSym(dis,structures=structures,strucpars=strucpars)
sres <- stop_smacofSym(dis,structures=structures,strucpars=strucpars,type="multiplicative")
sres <- stop_smacofSym(dis,theta=c(1,2,1),structures=structures,strucpars=strucpars,type="multiplicative")
sres <- stop_smacofSym(dis,theta=c(2,2,2),structures=structures,strucpars=strucpars,type="multiplicative")
sres <- stop_smacofSym(dis,theta=2,structures=structures,strucpars=strucpars,type="multiplicative")

#test stop_powerstress
pres <- stop_powerstress(dis,structures=structures,strucpars=strucpars)
pres <- stop_powerstress(dis,structures=structures,strucpars=strucpars,type="multiplicative")
pres <- stop_powerstress(dis,theta=c(1,2,1),structures=structures,strucpars=strucpars,type="multiplicative")
pres <- stop_powerstress(dis,theta=c(2,2,2),structures=structures,strucpars=strucpars,type="additive")
pres <- stop_powerstress(dis,theta=2,structures=structures,strucpars=strucpars,type="multiplicative")
pres <- stop_powerstress(dis,theta=c(2,2,-2),weightmat=dis,structures=structures,strucpars=strucpars,type="additive")


library(stops)
data(kinshipdelta)
resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity"),loss="stress",verbose=0,strucpars=list(c(eps=1000,minpts=2),NULL),type="additive",strucweight=c(-0.5,-0.5))

resa1<-stops(kinshipdelta,structures=c("cclusteredness","clinearity"),loss="stress",verbose=0,strucpars=list(c(eps=1000,minpts=2),NULL),type="additive")

resa2<-stops(kinshipdelta,structures=c("cclusteredness","clinearity"),loss="stress",verbose=0,strucpars=list(c(eps=1000,minpts=2),NULL),type="additive",strucweight=c(0.5,0.5))

resm<-stops(kinshipdelta,structures=c("cclusteredness","clinearity"),loss="stress",verbose=0,strucpars=list(c(eps=1000,minpts=2),NULL),type="multiplicative")

resm<-stops(kinshipdelta,structures=c("cclusteredness","clinearity"),loss="stress",verbose=0,strucpars=list(c(eps=1000,minpts=2),NULL),type="multiplicative",strucweight=c(-0.5,-0.5))

resm1<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cmanifoldness"),loss="stress",verbose=0,strucpars=list(c(eps=1000,minpts=2),NULL),type="multiplicative",strucweight=c(0.5,0.5))


#all stops losses
library(stops)
data(kinshipdelta)

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cfunctionality"),loss="stress",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,-0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cdependence"),loss="powerstress",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cassociation","cmanifoldness","clinearity"),loss="powermds",verbose=3,strucpars=list(NULL,NULL,NULL),type="additive",strucweight=c(-0.5,-0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cfunctionality"),loss="strain",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cfunctionality"),loss="powerstrain",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cfunctionality"),loss="powersammon",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cfunctionality"),loss="rstress",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cfunctionality"),loss="sstress",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cfunctionality"),loss="powerelastic",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cfunctionality"),loss="sammon",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","cfunctionality"),loss="smacofSphere",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,structures=c("cclusteredness","clinearity","ccomplexity"),loss="sammon2",verbose=3,strucpars=list(c(eps=1000,minpts=2),NULL,NULL),type="additive",strucweight=c(-0.5,0.5,0.5))
resa

resa<-stops(kinshipdelta,theta=c(1,1,1),structures=c("cfunctionality","clinearity"),loss="elastic",verbose=3,strucpars=c(list(NULL),list(NULL)),type="additive",strucweight=c(-0.5,0.5),lower=c(1,0.5,-2),upper=c(3,4,-2))
resa




#test for kapp=lambda=nu=1
dis <- smacof::kinshipdelta
kappa <- 1
lambda <- 1
nu <- 1
eps <- 10
minp <- 2
rango <- c(0,1.3) 
fit <- powerStressMin(dis,kappa=kappa,lambda=lambda,nu=nu)
c1<-cordillera(fit$conf,epsilon=eps,minpts=minp,rang=rango)

expect_equal(coploss(fit,epsilon=eps,minpts=minp,rang=rango,stressweight=0,cordweight=1),stoploss(fit,strucpars=list(list(epsilon=eps,minpts=minp,rang=rango)),structures="cclusteredness",stressweight=0,strucweight=-1))


#Bug form using the wrong list - huge cordillera - becuae this way q is set to 10 as we only say "eps" - epsilon is no problem?
resa<-stops(kinshipdelta,structures=c("cclusteredness"),loss="powermds",verbose=3,strucpars=list(eps=10,minpts=2,rang=c(0,1.3)),type="additive",strucweight=c(-1),stressweight=0)
resa
resa<-stops(kinshipdelta,structures=c("cclusteredness"),loss="powermds",verbose=3,strucpars=list(epsilon=10,minpts=2,rang=c(0,1.3)),type="additive",strucweight=c(-1),stressweight=0)
resa


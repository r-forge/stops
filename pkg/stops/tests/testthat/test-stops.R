
library(stops)
data(BankingCrisesDistances)
set.seed(210485)
optstrain <- cops(BankingCrisesDistances[,1:69],loss="strain",verbose=1)
optstress <- cops(BankingCrisesDistances[,1:69],loss="stress",verbose=3)
optsammon <- cops(BankingCrisesDistances[,1:69],loss="sammon",verbose=3)
optelastic <- cops(BankingCrisesDistances[,1:69],loss="elastic",verbose=3)
optsstress <- cops(BankingCrisesDistances[,1:69],loss="sstress",verbose=3)
optrstress <- cops(BankingCrisesDistances[,1:69],loss="rstress",verbose=3)
optpowerstress <- cops(BankingCrisesDistances[,1:69],loss="powerstress",verbose=3)
optpowersammon <- cops(BankingCrisesDistances[,1:69],loss="powersammon",verbose=3)
optpowerelastic <- cops(BankingCrisesDistances[,1:69],loss="powerelastic",verbose=3)

optstrain
optstress
optsammon
optelastic
optsstress
optrstress
optpowerstress
optpowersammon
optpowerelastic



plot(optstrain)
plot(optstress)
plot(optsammon)
plot(optelastic)
plot(optsstress)p
lot(optrstress)
plot(optpowerstress)
plot(optpowersammon)
plot(optpowerelastic)


data(kinshipdelta)
dis <- as.matrix(kinshipdelta)
res <- smacofSym(kinshipdelta)
conf <- res$conf

#test c_linearity
c_linearity(conf)

#test stoploss with clinearity additive
stressweight <- 1
structures <- c("clinearity")
strucweight <- rep(1/length(structures),length(structures))/2
type <- c("additive")
verbose <- 1
res$pars <- c(1,1,1)
res$stress.m <- res$stress
s1a <- stoploss(res,stressweight=stressweight,structures=structures,strucweight=strucweight,type=type,verbose=verbose)
s1a
#test stoploss with clinearity multiplicative
type <- c("multiplicative")
s1m <- stoploss(res,stressweight=stressweight,structures=structures,strucweight=strucweight,type=type,verbose=verbose)
s1m

#test with cclusteredness and clinearity additive
structures <- c("cclusteredness","clinearity")
strucpars <- list(c(eps=100,minpts=2),c(NULL))
strucweight <- rep(1/length(structures),length(structures))
type <- c("additive")
s2a <- stoploss(res,stressweight=stressweight,structures=structures,strucpars=strucpars,strucweight=strucweight,type=type,verbose=verbose)
s2a

#multiplicative
type <- c("multiplicative")
s2b <- stoploss(res,stressweight=stressweight,structures=structures,strucpars=strucpars,strucweight=strucweight,type=type,verbose=verbose)
s2b

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

#test stops
stopres <- stops(dis,theta=c(1,1,1),structures=structures,strucpars=strucpars,verbose=4)
stopres <- stops(dis,theta=c(1,1,1),lower=c(1,0.2,1),upper=c(1,10,1),structures=structures,strucpars=strucpars,verbose=4)
#more weight for clinearity
strucweights <- c(0.5,10)
stopres <- stops(dis,theta=c(1,1,1),lower=c(1,1,1),upper=c(1,10,1),structures=structures,strucpars=strucpars,strucweight=strucweights,verbose=4)


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





#there seems to be a bug in stops or upstream as
#- a<-stops(dis) cordillera(a$fit$conf) do not agree
#- b<-cops(dis) a nd b do not agree

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


library(stops)
data(Pendigits500)
strucpars1 <- list(list(k=10))
diso <- dist(Pendigits500[,1:16])
m1 <- stops(diso,structures=c("cfaithfulness"),loss="sammon",strucpars=strucpars1,verbose=3,stressweight=0,upper=c(5,10,1))
strucpars2 <- list(list(minpts=10,epsilon=1000))
m2 <- stops(diso,structures=c("cclusteredness"),loss="sammon",strucpars=strucpars2,verbose=3,stressweight=0,upper=c(5,10,1))
par(mfrow=c(1,2))
plot(m1,col=Pendigits500[,17],label.conf=list(label=FALSE))
reso <- sammon(diso^8)
plot(m2,col=Pendigits500[,17],label.conf=list(label=FALSE))
library(party)
library(caret)
ct1 <- ctree(factor(Pendigits500[,17])~m1$fit$conf[,1]+m1$fit$conf[,2])
confusionMatrix(predict(ct1),factor(Pendigits500[,17]))
ct2 <- ctree(factor(Pendigits500[,17])~m2$fit$conf[,1]+m2$fit$conf[,2])
confusionMatrix(predict(ct2),factor(Pendigits500[,17]))

####################################################
######### stops functions ##########################
####################################################
source("../code/optics.R")
source("../code/stops.R")

#fSmac
library(smacof)
data(kinshipdelta)
dis <- kinshipdelta
#defaults
res1 <- cop_smacof(dis)
#verbose
res1 <- cop_smacof(dis,verbose=1)
res1 <- cop_smacof(dis,verbose=2)
#plot and single theta argument, ndim argument
res2 <- cop_smacof(dis,3,ndim=3,a=0.5,q=1,minpts=2,epsilon=10,plot=TRUE)
#vector theta argument
res3 <- cop_smacof(dis,c(2,3),ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE) #the kappa value gets ignored as it should be
res3
#cordillera meta parameters: a, q, minpts, epsilon and rang 
res2 <- cop_smacof(dis,3,ndim=2,a=0.1,q=2,minpts=3,epsilon=1,verbose=0,plot=TRUE)
res2
res2a <- cop_smacof(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE)
res2a
res2b <- cop_smacof(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=c(0,1),plot=TRUE)
res2b
plot(res2a$fit)
plot(res2$fit)
#weightmat and init
w <- 1-diag(nrow(dis))
w[c(1,2),c(1,2)] <- 0
res <- cop_smacof(dis,3)
res2 <- cop_smacof(dis,3,weightmat=w)
res
res2
res2a <- cop_smacof(dis,3,init=res2$fit$conf,weightmat=w)
res2a

#in optimization
rang <- c(0,1.45)
smacop <- ljoptim(1,function(lambda) cop_smacof(dis,lambda,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)$copstress,lower=c(0.5,0.5),upper=c(5,5))
smacop
cop_smacof(dis,smacop$par,ndim=2,a=0.5,q=1,minpts=2,eps=10,verbose=0,rang=rang,plot=TRUE)
cop_smacof(dis,4.9,ndim=2,a=0.5,q=1,minpts=2,eps=10,verbose=0,rang=rang,plot=TRUE)

################## sammon
#initsam <- sammon(dis)
#colsopt <- cols1 <- colstraj <- factor(pendss[,17])
#library(colorspace)
#levels(colsopt) <- rainbow_hcl(10,c=70,l=40)
#levels(cols1) <- rainbow_hcl(10,c=70,l=80)
#levels(colstraj) <- rainbow_hcl(10,c=70,l=90)
#plot(initsam$points[,1],initsam$points[,2],col=as.character(colsopt),type="n",xlab="D1",ylab="D2")
#points(initsam$points[,1],initsam$points[,2],col=as.character(cols1),pch=as.character(pendss[,17]))
#par(mfrow=c(1,1))
#c1 <- cordillera(scale(initsam$points),minpts=5,eps=100,verbose=TRUE)

#fSamm
library(smacof)
data(kinshipdelta)
dis <- kinshipdelta
#defaults
res1 <- cop_sammon(dis)
res1
#verbose
res1 <- cop_sammon(dis,verbose=1)
res1 <- cop_sammon(dis,verbose=2)
#plot and single theta argument, ndim argument
res2 <- cop_sammon(dis,3,ndim=3,a=0.5,q=1,minpts=2,epsilon=10,plot=TRUE)
#vector theta argument
res3 <- cop_sammon(dis,c(2,3),ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE) #the kappa value gets ignored as it should be
res3
#cordillera meta parameters: a, q, minpts, epsilon and rang 
res2 <- cop_sammon(dis,3,ndim=2,a=0.1,q=2,minpts=3,epsilon=1,verbose=0,plot=TRUE)
res2
res2a <- cop_sammon(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE)
res2a
res2b <- cop_sammon(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=c(0,1),plot=TRUE)
res2b
plot(res2a$fit$conf)
plot(res2$fit$conf)
#weightmat and init
res2 <- cop_sammon(dis,3,verbose=1)
res2
res2a <- cop_sammon(dis,3,init=res2$fit$conf,verbose=1)
res2a

#in optimization
rang <- c(0,1.45)
samcop <- ljoptim(1,function(lambda) cop_sammon(dis,lambda,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)$copstress,lower=c(0.5,0.5),upper=c(5,5))
samcop
cop_sammon(dis,smacop$par,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)
cop_sammon(dis,1,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)


#rang <- c(0,1.5)
#lambda_sam <- optim(1,function(lambda) cop_sammon  (dis,lambda,a=0.5,q=1,minpts=5,eps=10,verbose=1,plot=TRUE,rang=rang)$copStress,method="Brent",lower=0.5,upper=4,control=list(reltol=0.001))
#lopt <- lambda_sam$par
#samopt <- sammon(dis^lopt)

#fCmdscale
library(smacof)
data(kinshipdelta)
dis <- kinshipdelta
#defaults
res1 <- cop_cmdscale(dis)
res1
#verbose
res1 <- cop_cmdscale(dis,verbose=1)
res1 <- cop_cmdscale(dis,verbose=2)
#plot and single theta argument, ndim argument
res2 <- cop_cmdscale(dis,3,ndim=3,a=0.5,q=1,minpts=2,epsilon=10,plot=TRUE)
#vector theta argument
res3 <- cop_cmdscale(dis,c(2,3),ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE) #the kappa value gets ignored as it should be
res3
#cordillera meta parameters: a, q, minpts, epsilon and rang 
res2 <- cop_cmdscale(dis,3,ndim=2,a=0.1,q=2,minpts=3,epsilon=1,verbose=0,plot=TRUE)
res2
res2a <- cop_cmdscale(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE)
res2a
res2b <- cop_cmdscale(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=c(0,1),plot=TRUE)
res2b
plot(res2a$fit$conf)
plot(res2$fit$conf)
#weightmat and init
res2 <- cop_cmdscale(dis,3,verbose=1)
res2
res2a <- cop_cmdscale(dis,3,init=res2$fit$conf,verbose=1)
res2a

#in optimization
rang <- c(0,1.45)
cmdscalecop <- ljoptim(1,function(lambda) cop_cmdscale(dis,lambda,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)$copstress,lower=c(0.5,0.5),upper=c(5,5))
cmdscalecop
res1 <- cop_cmdscale(dis,cmdscalecop$par,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)
res2 <- cop_cmdscale(dis,1,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)
plot(res1$fit$conf)
dev.new()
plot(res2$fit$conf)

#frStress
source("../code/optics.R")
source("../code/stops.R")
source("../code/rstress.R")
library(smacof)
data(kinshipdelta)
dis <- kinshipdelta
#defaults
res1 <- cop_rstress(dis)
#verbose
res1 <- cop_rstress(dis,verbose=1)
res1 <- cop_rstress(dis,verbose=2)
#plot and single theta argument, ndim argument
res2 <- cop_rstress(dis,3,ndim=3,a=0.5,q=1,minpts=2,epsilon=10,plot=TRUE)
#vector theta argument
res3 <- cop_rstress(dis,c(2,3),ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE) #the lambda value gets ignored as it should be
res3
#cordillera meta parameters: a, q, minpts, epsilon and rang 
res2 <- cop_rstress(dis,3,ndim=2,a=0.1,q=2,minpts=3,epsilon=1,verbose=0,plot=TRUE)
res2
res2a <- cop_rstress(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE)
res2a
res2b <- cop_rstress(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=c(0,1),plot=TRUE)
res2b
plot(res2x$fit$conf)
dev.new()
plot(res2$fit$conf)
res2 <- cop_rstress(dis,3,ndim=2,a=0.1,q=1,minpts=4,epsilon=10,verbose=0,plot=TRUE)
res2
res2x <- fpowerStress(dis,3,ndim=2,a=0.1,q=1,minpts=4,epsilon=10,verbose=0,plot=TRUE)
#weightmat and init
w <- 1-diag(nrow(dis))
w[c(1,2),c(1,2)] <- 0
res <- cop_rstress(dis,3)
res2 <- cop_rstress(dis,3,weightmat=w)
res
res2
res2a <- cop_rstress(dis,3,init=res2$fit$conf,weightmat=w) 
res2a

#in optimization
rang <- c(0,1.45)
powercop <- ljoptim(c(1,1),function(theta) cop_rstress(dis,theta,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)$copstress,lower=c(1,0.5),upper=c(5,5))
powercop
resopt <- cop_rstress(dis,powercop$par,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)
res1 <- cop_rstress(dis,c(1,1),ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)

res1
resopt
plot(resopt$fit$conf)
plot(res1$fit$conf)


#cop_powerstress
source("../code/optics.R")
source("../code/stops.R")
source("../code/rstress.R")
library(smacof)
data(kinshipdelta)
dis <- kinshipdelta
#defaults
res1 <- cop_powerstress(dis)
#verbose
res1 <- cop_powerstress(dis,verbose=1)
res1 <- cop_powerstress(dis,verbose=2)
#plot and single theta argument, ndim argument
res2 <- cop_powerstress(dis,3,ndim=3,a=0.5,q=1,minpts=2,epsilon=10,plot=TRUE)
#vector theta argument
res3 <- cop_powerstress(dis,c(2,3),ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE) #the kappa value gets ignored as it should be
res3
#cordillera meta parameters: a, q, minpts, epsilon and rang 
res2 <- cop_powerstress(dis,3,ndim=2,a=0.1,q=2,minpts=3,epsilon=1,verbose=0,plot=TRUE)
res2
res2a <- cop_powerstress(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE)
res2a
res2b <- cop_powerstress(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=c(0,1),plot=TRUE)
res2b
plot(res2a$fit$conf)
plot(res2$fit$conf)
#weightmat and init
w <- 1-diag(nrow(dis))
w[c(1,2),c(1,2)] <- 0
res <- cop_powerstress(dis,3)
res2 <- cop_powerstress(dis,3,weightmat=w)
res
res2
res2a <- cop_powerstress(dis,3,init=res2$fit$conf,weightmat=w) 
res2a

#in optimization
rang <- c(0,1.45)
powercop <- ljoptim(c(1,1),function(theta) cop_powerstress(dis,theta,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)$copstress,lower=c(1,0.5),upper=c(5,5))
powercop
resopt <- cop_powerstress(dis,powercop$par,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)
res1 <- cop_powerstress(dis,c(1,1),ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)

res1
resopt
plot(resopt$fit$conf)
plot(res1$fit$conf)


#coploss
res1 <- cop_powerstress(dis,c(1,1),ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=FALSE)
fit <- res1$fit
coploss(fit,a=a,q=q,normed=FALSE,minpts=2,epsilon=10,rang=rang,verbose=1,plot=FALSE,scale=TRUE)
coploss(fit,a=a,q=q,normed=TRUE,minpts=2,epsilon=10,rang=rang,verbose=1,plot=FALSE,scale=TRUE)
coploss(fit,a=a,q=q,normed=TRUE,minpts=3,epsilon=10,rang=rang,verbose=1,plot=FALSE,scale=TRUE)
coploss(fit,a=a,q=q,normed=TRUE,minpts=2,epsilon=0.5,rang=rang,verbose=1,plot=TRUE,scale=TRUE)
coploss(fit,a=a,q=q,normed=TRUE,minpts=2,epsilon=10,rang=c(0,0.5),verbose=1,plot=FALSE,scale=TRUE)
coploss(fit,a=a,q=q,normed=TRUE,minpts=2,epsilon=10,rang=rang,verbose=1,plot=FALSE,scale=FALSE)

#cops
dis <- kinshipdelta
test1 <- cops(dis,verbose=2)
plot(test1$fit$conf[,1],test1$fit$conf[,2])
test2 <- cops(dis,loss="strain",verbose=2)
plot(test2$fit$conf[,1],test2$fit$conf[,2])
test3 <- cops(dis,loss="stress",verbose=2)
plot(test3)
test3 <- cops(dis,loss="stress",smacoffunc="smacofSphere",verbose=2)
plot(test3)
test4 <- cops(dis,loss="rstress",verbose=2)
plot(test4$fit$conf[,1],test4$fit$conf[,2])
test5 <- cops(dis,loss="sammon",verbose=2)
plot(test5$fit$conf[,1],test5$fit$conf[,2])

test6 <- cops(dis,loss="strain",optimmethod="pso")
test7 <- cops(dis,loss="strain",ndim=3)
test9 <- cops(dis,loss="strain",ndim=2,q=2,minpts=3)
test10 <- cops(dis,loss="strain",ndim=2,q=2,epsilon=2)
test11 <- cops(dis,loss="strain",ndim=2,q=2,epsilon=2,verbose=2)
test12 <- cops(dis,loss="strain",ndim=2,q=2,epsilon=2,plot=TRUE)
test13 <- cops(dis,loss="strain",ndim=2,q=2,epsilon=2,scale=FALSE)



library(MASS)
library(stops)
data(BankingCrisesDistances)
set.seed(210485)
opto <- cops(BankingCrisesDistances[,1:69],theta=1,loss="sammon",verbose=2,acc=1e-16,accd=1e-12,cordweight=0.5,lower=0.5,upper=5)

opto <- cops(BankingCrisesDistances[,1:69],theta=1,loss="sammon",verbose=2,acc=1e-12,accd=1e-12,lower=0.5,upper=5,scale=FALSE)

init <- cmdscale(BankingCrisesDistances[,1:69],k=2)
init2 <- cmdscale(BankingCrisesDistances[,1:69]^opto$par[2],k=2)
res1 <- cop_sammon(BankingCrisesDistances[,1:69],cordweight=0.5,init=init,rang=c(0,opto$cordillera$cordillera$dmax),verbose=1)
resopt <- cop_sammon(BankingCrisesDistances[,1:69],theta=opto$par[2],init=init2,cordweight=0.00141365,rang=c(0,opto$cordillera$cordillera$dmax),verbose=3)
resoptx <- cop_sammon(BankingCrisesDistances[,1:69]^opto$par[2],theta=1,init=init2,cordweight=0.5,rang=c(0,opto$cordillera$cordillera$dmax),verbose=3)
#da passt noch was nicht! resopt ist tatsaechlich besser aber warum ist es anders als das andere?
resopto <- cop_sammon(BankingCrisesDistances[,1:69],theta=opto$par[2],cordweight=0.5,rang=c(0,opto$cordillera$cordillera$dmax),verbose=3)
resopto <- MASS::sammon(BankingCrisesDistances[,1:69]^opto$par[2])
reso1 <- MASS::sammon(BankingCrisesDistances[,1:69])


opto$par[2]
cropt <- cordillera(resopt$fit$points)
cr1 <- cordillera(res1$fit$points)
cr12 <- cordillera(res1$fit$points,rang=c(0,opto$cordillera$cordillera$dmax))
cropt2 <- cordillera(resopt$fit$points,rang=c(0,opto$cordillera$cordillera$dmax))

cropt <- cordillera(resopt$fit$points,rang=c(0,opto$cordillera$cordillera$dmax),plot=TRUE)
cr1<- cordillera(res1$fit$points,rang=c(0,opto$cordillera$cordillera$dmax),plot=TRUE)

plot(resopt$fit$points)
plot(res1$fit$points)





library(smacof)
library(stops)
data(kinshipdelta)
resp <- powerStressMin(kinshipdelta)
plot(resp$confdiss,resp$obsdiss)

plot(unclass(resp$confdiss),unclass(resp$obsdiss),ylim=c(0,max(resp$confdiss,resp$obsidss)),xlim=c(0,max(resp$confdiss,resp$obsidss)))

ekman <-
structure(c(0.86, 0.42, 0.42, 0.18, 0.06, 0.07, 0.04, 0.02, 0.07,
0.09, 0.12, 0.13, 0.16, 0.5, 0.44, 0.22, 0.09, 0.07, 0.07, 0.02,
0.04, 0.07, 0.11, 0.13, 0.14, 0.81, 0.47, 0.17, 0.1, 0.08, 0.02,
0.01, 0.02, 0.01, 0.05, 0.03, 0.54, 0.25, 0.1, 0.09, 0.02, 0.01,
0, 0.01, 0.02, 0.04, 0.61, 0.31, 0.26, 0.07, 0.02, 0.02, 0.01,
0.02, 0, 0.62, 0.45, 0.14, 0.08, 0.02, 0.02, 0.02, 0.01, 0.73,
0.22, 0.14, 0.05, 0.02, 0.02, 0, 0.33, 0.19, 0.04, 0.03, 0.02,
0.02, 0.58, 0.37, 0.27, 0.2, 0.23, 0.74, 0.5, 0.41, 0.28, 0.76,
0.62, 0.55, 0.85, 0.68, 0.76), Size = 14L,
call = quote(as.dist.default(m = b)), class = "dist",
Diag = FALSE, Upper = FALSE, Labels = c(434,
445, 465, 472, 490, 504, 537, 555, 584, 600, 610, 628, 651, 674
))
ekman <- as.matrix (1-ekman)
wave <- row.names(ekman)




#these are r stresses 
e05<-powerStressMin(ekman,kappa=2*0.1,lambda=1,verbose=2,eps=1e-8)
e05
shep(e05)
plot(e05$conf,type="n")
text(e05$conf,labels=row.names(ekman))
e10<-powerStressMin(ekman,lambda=1,kappa=2*.10,verbose=2,eps=1e-10)
e10
shep(e10)
e25<-powerStressMin(ekman,kappa=2*0.25,lambda=1,verbose=2)
e25

e25<-powerStressMin(ekman,kappa=2*0.5,lambda=1,verbose=1)
e25

e25s<-smacofSym(ekman)
e25s

shep <- function(obj,orig=ekman)
    {
        dout <- (2*sqrt(sqdist(obj$conf)))^obj$pars[1]
        plot(orig^obj$pars[2],dout,main=paste("stressn:",round(obj$stress.m,4),"stresse",round(obj$stress.e,4)))
    }


shep2 <- function(obj)
    {
        dout <- obj$confdiss
        plot(as.matrix(obj$delta),as.matrix(dout),main=paste("stressn:",round(obj$stress.m,4),"stresse",round(obj$stress,4)))
        sum((obj$delta-dout)^2)/sum(obj$delta^2)
    }

shep2(e25)

k1 <- smacofSym(kinshipdelta)
k2 <- powerStressMin(kinshipdelta,verbose=2)
delt
k2 <- powerStressMin(delt,verbose=2)


k1
k2
plot(k1)
dev.new()
plot(k2)
summary(k1)
summary(k2)



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


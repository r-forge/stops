##### Cordillera: last passed 2.12.2014
source("optics.R")
library(smacof)
data(kinshipdelta)
dis <- kinshipdelta
lamb <- 1
eps <- 20 #find defaults
minpts <- 2 #find defaults - (1 makes no snese, 2 is hierarcical sinlge linkage) rule of thumb dimension + 1, so in our case almost always 3; the more n the larger too 
q <- 1
verbose <- TRUE
rescoy <- smacofSym(dis^lamb, itmax = 5000)
confs <- rescoy$conf #make matrix from fitted distances; the diagonal is zero
c1 <- cordillera(confs)
c1
c1 <- cordillera(confs,epsilon=eps)
c1
c1 <- cordillera(confs,q=q,plot=TRUE)
c1
c1 <- cordillera(confs,minpts=3)
c1
c1 <- cordillera(confs,rang=c(0,10),minpts=3,plot=TRUE)
c1
c1 <- cordillera(confs,rang=c(0,10),minpts=3,ylim=c(0,1.5),plot=FALSE)
c1
c1 <- cordillera(confs,digits=1,minpts=3,plot=FALSE)
c1
c1 <- cordillera(confs,digits=1,minpts=3,path="~")
c1
#figure 3: But no rang argument so they are normed to the observed maximum, not the maximum of the last configuration.
library(MASS)
par(mfrow=c(3,2))
q <- 1
eps <- 2
x <- c(-0.25,-0.125,0.125,0.25,0,0,0,0)
y <- c(0,0,0,0,-0.25,-0.125,0.125,0.25)
x <- c(-0.375,-0.125,0.125,0.375,-0.375,-0.125,0.125,0.375)
x <- x*2
y <- c(-0.125,-0.125,-0.125,-0.125,0.125,0.125,0.125,0.125)
y <- y*2
c1 <- cordillera(cbind(x,y),eps=eps,q=q)
eqscplot(cbind(x,y),ylim=c(-1,1),xlim=c(-1,1),pch=19,main=paste("normed structure=",round(c1$normed,2)),xlab="D1",ylab="D2")
obj <- 8
 mat <- matrix(rep(1,obj^2),ncol=obj)
 diag(mat) <- 0
 mds <- smacofSym(mat)
# plot(mds)
 mdsmat <-as.matrix(mds$conf)
 confdiss <- as.matrix(mds$confdiss)
c1 <- cordillera(mdsmat,eps=eps,q=q)
eqscplot(mdsmat,ylim=c(-1,1),xlim=c(-1,1),pch=19,main=paste("normed structure=",round(c1$normed,2)),xlab="D1",ylab="D2")
#more strcuture
mdsmat <- cbind(rep(c(0,0,-0.5,0.5),2),rep(c(0.5,-0.5,0,0),2))
set.seed(123)
mdsmat <- jitter(mdsmat,factor=3.5) 
c1 <- cordillera(mdsmat,eps=eps,q=q)
eqscplot(mdsmat,ylim=c(-1,1),xlim=c(-1,1),pch=19,main=paste("normed structure=",round(c1$normed,2)),xlab="D1",ylab="D2")
mdsmat <- cbind(rep(c(0,0,-0.5,0.5),2),rep(c(0.5,-0.5,0,0),2))
set.seed(123)
mdsmat <- jitter(mdsmat,factor=1.7) 
c1 <- cordillera(mdsmat,eps=eps,q=q)
eqscplot(mdsmat,ylim=c(-1,1),xlim=c(-1,1),pch=19,main=paste("normed structure=",round(c1$normed,2)),xlab="D1",ylab="D2")
#more structure
mdsmat <- cbind(rep(c(0,0,-0.5,0.5),2),rep(c(0.5,-0.5,0,0),2))
set.seed(210485)
mdsmat <- jitter(mdsmat,factor=1) 
c1 <- cordillera(mdsmat,eps=eps,q=q)
eqscplot(mdsmat,ylim=c(-1,1),xlim=c(-1,1),pch=19,xlab="D1",ylab="D2",main=paste("normed structure=",round(c1$normed,2)))
ind <- matrix(c(1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1),ncol=8)
fake <- matrix(1,nrow=8,ncol=8)
dissi <- fake*(1-ind)
exdiss <- as.dist(dissi)
exmds <- smacofSym(exdiss)
mdsmat <- exmds$conf
c1 <- cordillera(mdsmat,q=q,eps=eps)
eqscplot(mdsmat,ylim=c(-1,1),xlim=c(-1,1),pch=19,cex=2,xlab="D1",ylab="D2",main=paste("normed structure=",round(c1$normed,2)))
par(mfrow=c(1,1))

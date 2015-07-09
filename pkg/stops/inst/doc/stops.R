## ------------------------------------------------------------------------
library(stops)

## ------------------------------------------------------------------------
dis<-as.matrix(smacof::kinshipdelta)

## ------------------------------------------------------------------------
res1<-powerStressMin(dis,kappa=1,lambda=1)
res1

## ------------------------------------------------------------------------
res2<-powerStressMin(dis,kappa=1,lambda=1,nu=-1,weightmat=dis)
res2

## ------------------------------------------------------------------------
res2a<-sammon(dis)
res2a

## ------------------------------------------------------------------------
res3<-powerStressMin(dis,kappa=1,lambda=1,nu=-2,weightmat=dis)
res3

## ------------------------------------------------------------------------
res4<-powerStressMin(dis,kappa=2,lambda=2)
res4

## ------------------------------------------------------------------------
res5<-powerStressMin(dis,kappa=2,lambda=1)
res5

## ------------------------------------------------------------------------
res6<-powerStressMin(dis,kappa=2,lambda=1.5)
res6

## ------------------------------------------------------------------------
res7<-powerStressMin(dis,kappa=2,lambda=1.5,nu=-1,weightmat=dis)
res7

## ------------------------------------------------------------------------
res8<-powerStressMin(dis,kappa=2,lambda=1.5,nu=-2,weightmat=dis)
res8

## ------------------------------------------------------------------------
res9<-powerStressMin(dis,kappa=2,lambda=1.5,nu=-1.5,weightmat=2*1-diag(nrow(dis)))
res9
summary(res9)

## ----eval=FALSE,fig.show='hold',fig.width=8, fig.height = 8--------------
#  plot(res9)
#  plot(res9,"transplot")
#  plot(res9,"Shepard")
#  plot(res9,"resplot")
#  plot(res9,"bubbleplot")

## ------------------------------------------------------------------------
resc<-cmdscale(kinshipdelta^3)
resc

## ------------------------------------------------------------------------
summary(resc)

## ----eval=FALSE,fig.show='hold',fig.width=8, fig.height=8----------------
#  summary(resc)
#  plot(resc)

## ------------------------------------------------------------------------
ressm<-stops(kinshipdelta,loss="stress",stressweight=1,structures=c("cclusteredness","clinearity","cdependence"),strucweight=c(-0.33,0.33,-0.33),verbose=0,strucpars=list(c(eps=10,minpts=2),NULL,c(index=1)),type="multiplicative")
ressm

## ----fig.show='hold',fig.width=8,fig.height=8----------------------------
plot(ressm)

## ------------------------------------------------------------------------
ressa<-stops(kinshipdelta,loss="stress",stressweight=1,structures=c("cclusteredness","clinearity","cdependence"),strucweight=c(-0.33,0.33,-0.33),verbose=0,strucpars=list(c(eps=10,minpts=2),NULL,c(index=1)),type="additive")
ressa

## ----fig.show='hold',fig.width=8,fig.height=8----------------------------
plot(ressa)

## ----eval=FALSE----------------------------------------------------------
#  ressa<-stops(kinshipdelta,structure=c("clinearity"),strucweight=0,loss="stress",verbose=0)

## ------------------------------------------------------------------------
resc<-cops(kinshipdelta,loss="strain")
resc
summary(resc)

## ----fig.show='hold',fig.width=8,fig.height=8----------------------------
plot(resc,"confplot")
plot(resc,"Shepard")
plot(resc,"transplot")
plot(resc,"reachplot")

## ------------------------------------------------------------------------
resca<-cops(kinshipdelta,cordweight=0,loss="strain")
resca

## ------------------------------------------------------------------------
rescb<-cops(kinshipdelta,cordweight=20,loss="strain")
rescb

## ----fig.show='hold',fig.width=8,fig.height=8----------------------------
plot(resca,main="with cordweight=0")
plot(rescb,main="with cordweight=20")

## ------------------------------------------------------------------------
data(iris)
res<-optics(iris[,1:4],minpts=2,epsilon=1000)
print(res)
summary(res)

## ----fig.show='hold',fig.width=8,fig.height=8----------------------------
plot(res,withlabels=TRUE)

## ----fig.show='hold',fig.width=8,fig.height=8----------------------------
cres<-cordillera(iris[,1:4],minpts=2,epsilon=1000,scale=FALSE)
cres
summary(cres)

## ----fig.show='hold',fig.width=8,fig.height=8----------------------------
plot(cres)

## ------------------------------------------------------------------------
set.seed(210485)
fwild <- function (x) 10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80
res2<-ljoptim(50, fwild,lower=-50,upper=50,adaptive=FALSE,accd=1e-16,acc=1e-16)
res2

## ----fig.show='hold',fig.width=8,fig.height=8----------------------------
plot(fwild, -50, 50, n = 1000, main = "ljoptim() minimising 'wild function'")
points(res2$par,res2$value,col="red",pch=19)

## ----eval=FALSE----------------------------------------------------------
#  conf_adjust(conf1,conf2)


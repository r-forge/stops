context("Powerstress")

library(stops)
library(smacof)
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

e2<-powerStressMin(ekman,kappa=2*1,lambda=1)
e1<-powerStressMin(ekman,kappa=2*0.5,lambda=1)
e05<-powerStressMin(ekman,kappa=2*0.25,lambda=1)

test_that("Power stress for Ekman data works",{
#these are r stresses 
expect_equal_to_reference(e05<-powerStressMin(ekman,kappa=2*0.25,lambda=1))
expect_that(e05,is_a("smacofP"))
})

test_that("PowerStress equals smacofSym stress to the third digit",{
s1 <- smacofSym(ekman)
expect_that(round(e1$stress.en1,3),equals(round(s1$stress,3)))
})


test_that("summary of power stress works",{
expect_that(summary(e05),is_a("summary.smacofP"))
expect_that(str(summary(e05)),prints_text("List of 2"))
expect_that(str(summary(e05)),prints_text("$ conf",fixed=TRUE))
expect_that(str(summary(e05)),prints_text("$ sppmat",fixed=TRUE))
})


test_that("output is verbose",{
              expect_that(e2 <- powerStressMin(ekman,kappa=1,lambda=1,verbose=0),not(prints_text("Minimizing powerStress with kappa= 1 lambda= 1")))
              expect_that(e2 <- powerStressMin(ekman,kappa=1,lambda=1,verbose=1),prints_text("Minimizing powerStress with kappa= 1 lambda= 1"))
              expect_that(e2 <- powerStressMin(ekman,kappa=1,lambda=1,verbose=2),prints_text("Minimizing powerStress with kappa= 1 lambda= 1"))
              expect_that(e2 <- powerStressMin(ekman,kappa=1,lambda=1,verbose=2),prints_text("raw stress:",fixed=TRUE))
          })
              

test_that("Transformation and Shepard Plots work correctly",{
              expect_that(cor(e05$confdiss,e05$delta),equals(0.9905080545))
              expect_that(cor(e05$confdiss,e05$delta),is_more_than(cor(e1$confdiss,e1$delta)))
              expect_that(cor(e05$confdiss,e05$delta),is_more_than(cor(e2$confdiss,e2$delta)))
            #  expect_that(expect_that(plot(e06,plot.type="NLShepard"),throws_error()),throws_error())
              plot(e1,"transplot")
              plot(e2,"transplot")
              plot(e05,"Shepard")
              plot(e1,"Shepard")
              plot(e2,"Shepard")
          })


test_that("Power stresses are similar when nu is used versus direct weighting",{
              e2<-powerStressMin(ekman,kappa=2*1,lambda=1,nu=1)
              e1<-powerStressMin(ekman,kappa=2*1,lambda=1)
              expect_equal(e1$conf,e2$conf)
              e2<-powerStressMin(ekman,kappa=2*1,lambda=2,nu=-1,weightmat=ekman)
              ekinv <- 1/ekman
              diag(ekinv) <- 1
              e1<-powerStressMin(ekman,kappa=2*1,lambda=2,weightmat=ekinv)
              expect_equal(e1$conf,e2$conf)
              e2<-powerStressMin(ekman,kappa=2*1,lambda=2,nu=-2,weightmat=ekman)
              e1<-powerStressMin(ekman,kappa=2*1,lambda=2,weighmat=1/ekman^2)
              expect_equal(e1$conf,e2$conf)
          })

test_that("Power stresses are equal when nu is used as 0 with identity matrix for w",{
              e2<-powerStressMin(ekman,kappa=2*1,lambda=1,nu=0)
              e1<-powerStressMin(ekman,kappa=2*1,lambda=1,nu=1,weightmat=1-diag(14))
              expect_equal(e1$conf,e2$conf)
          })



#these are r stresses 
e05<-powerStressMin(ekman,kappa=2*0.05,lambda=1,verbose=2,eps=1e-8)
e05
shep(e05)
plot(e05$conf,type="n")
text(e05$conf,labels=row.names(ekman))
e10<-powerStressMin(ekman,lambda=1,kappa=2*.10,verbose=2,eps=1e-10)
e10
shep(e10)
e25<-powerStressMin(ekman,kappa=2*0.25,lambda=1,verbose=2)
e25

shep(e25)

e34n<-powerStressMin(ekman,lambda=1,kappa=2*.75,verbose=1)


t0.05 <- powerStressMin(ekman,lambda=1,kappa=2*0.05,verbose=1,eps=1e-10)
shep(t0.05)
t0.1 <-powerStressMin(ekman,lambda=1,kappa=2*.10,verbose=1,eps=1e-10)
shep(t0.1)
t0.125 <-powerStressMin(ekman,lambda=1,kappa=2*0.125,verbose=1,eps=1e-10)
shep(t0.125)
t0.25 <-powerStressMin(ekman,lambda=1,kappa=2*0.25,verbose=1,eps=1e-10)
shep(t0.25)
t0.5 <-powerStressMin(ekman,lambda=1,kappa=2*0.5,verbose=1,eps=1e-10)
shep(t0.5)
t0.75 <-powerStressMin(ekman,lambda=1,kappa=2*.75,verbose=1,eps=1e-10)
shep(t0.75)
t1 <-powerStressMin(ekman,lambda=1,kappa=2*1,verbose=1,eps=1e-10)
shep(t1)
t2 <-powerStressMin(ekman,lambda=1,kappa=2*2,verbose=1,eps=1e-10)
shep(t2)
t3 <-powerStressMin(ekman,lambda=1,kappa=2*3,verbose=1,eps=1e-10)
shep(t3)


t0.25.0.75 <-powerStressMin(ekman,lambda=0.75,kappa=2*.255,verbose=1,eps=1e-10)
shep(t0.25.0.75)
t0.25.1 <-powerStressMin(ekman,lambda=1,kappa=2*.255,verbose=1,eps=1e-10)
shep(t0.25.1)
t0.25.1.25 <-powerStressMin(ekman,lambda=1.25,kappa=2*.255,verbose=1,eps=1e-10)
shep(t0.25.1.25)
t0.25.1.5 <-powerStressMin(ekman,lambda=1.5,kappa=2*.255,verbose=1,eps=1e-10)
shep(t0.25.1.5)
t0.25.2 <-powerStressMin(ekman,lambda=2,kappa=2*.255,verbose=1,eps=1e-10)
shep(t0.25.2)
t0.25.3 <-powerStressMin(ekman,lambda=3,kappa=2*.255,verbose=1,eps=1e-10)
shep(t0.25.3)
t0.25.4 <-powerStressMin(ekman,lambda=4,kappa=2*.255,verbose=1,eps=1e-10)
shep(t0.25.4)



shep <- function(obj,orig=ekman)
    {
        dout <- (2*sqrt(sqdist(obj$conf)))^obj$pars[1]
        plot(orig^obj$pars[2],dout,main=paste("stressn:",round(obj$stress.m,4),"stresse",round(obj$stress.e,4)))
    }

shep(e25,orig=ekman)

max((2*sqrt(sqdist(e25$conf)))^0.25)

e34n
#normal stress


data(kinshipdelta)
delta <- kinshipdelta

delta <- ekman
weightmat <- 1-diag(nrow(delta))
delta <- ekman/enorm(ekman,weightmat)
deltaold <- delta*enorm(ekman)
eone1<-powerStressMin(deltaold,lambda=3,kappa=1,verbose=1,eps=1e-16)


eone2<-powerStressMin(deltaold,lambda=2,kappa=1,verbose=1,eps=1e-16)
eone2

library(smacof)
smac1 <- smacofSym(deltaold^2)
smac1


smac1$stress
(smac1$stress^2)*sum(smac1$confdiss^2)
(smac1$stress^2)*sum(smac1$confdiss^2)/sum(smac1$obsdiss^2)

sx <- smac1$conf
px <- eone2$conf

sx/px

smadelta <- as.matrix(smac1$obsdiss)
powerStressMin(smadelta)

smadout <- as.matrix(smac1$confdiss)

smadelta <- deltaold^2
smadout <- as.matrix(smac1$confdiss)
sum((smadelta-smadout)^2) #raw stress
sum((smadelta-smadout)^2)/sum(smadelta^2) #enorm stress
sqrt(sum((smadelta-smadout)^2)/sum(smadout^2)) #enorm stress


w1 <- delta/smadelta
diag(w1) <- mean(w1[!is.na(w1)])
w1

pdelt <- w1*smadelta

dout <- as.matrix(dist(eone2$conf))
w2 <- dout/smadout
diag(w2) <- mean(w2[!is.na(w2)])

pdout <- w2*smadout
pdout <- w2*smadout

sum((smadelta-smadout*w2)^2) #raw stress
sum((smadelta-smadout)^2)/sum(smadelta^2) #enorm stress
sqrt(sum((smadelta-smadout)^2)/sum(smadout^2)) #enorm stress

sum((pdelt-pdout)^2) #raw stress
sum((pdelt/w1-pdout/w2)^2) #raw stress
sum((smadelta*w1-smadout*w2)^2) #raw stress

sum((pdelt-pdout)^2) #raw stress
sum((pdelt-pdout)^2)/sum(pdelt^2) #enorm stress
sqrt(sum((pdelt-pdout)^2)/sum(pdout^2)) #enorm stress



pdelt <- smadelta
pdout <- smadout

delta <- pdelt
xnew <- eone2$conf
dout <- 2*sqrt(sqdist(xnew))  #added euclidean distance
stressr <- sum(weightmat*(dout^kappa-delta)^2) #explicitly normed stress
stressr
stressn <- stressr/(sum(weightmat*delta^2)) #normalized to the maximum stress
stressn

stressr <- stresse*(sum(weightmat*deltaold^2)) #raw stress as explicitly normed stress times \sum w deltaold^2
     stress1 <- sqrt(stressr/sum(weightmat*(dout^(2*kappa)))) #implicitly noremd stress
stress1
stressn <- stressr/(sum(weightmat*deltaold^2)) #normalized to the maximum stress
stressn


deltaold <- delta

etwo<-powerStressMin(ekman,lambda=1,kappa=2,verbose=1)
etwo
efour<-powerStressMin(ekman,lambda=1,kappa=4, verbose=1)
efour
esix <- powerStressMin (ekman, lambda=1, kappa = 6, verbose = 0)
esix
#powerstress
epart <- powerStressMin (ekman, lambda = 2, verbose = 1)
epart <- powerStressMin (ekman, kappa = 3, verbose = TRUE)
epart <- powerStressMin (ekman, kappa = 3.4, lambda=0.666, verbose = TRUE)

#test w argument
w <- 1-diag(nrow(ekman)) 
w[c(1,2),c(1,2)] <- 0 #without distance between 1 and 2
epart <- powerStressMin (ekman, verbose = TRUE)
epart0 <- powerStressMin(ekman, lambda=1,kappa=1,w=w)
if(!isTRUE(all.equal(epart$conf,epart0$conf))) cat("PASSED \n")
plot(epart$conf,type="n")
text(epart$conf,labels=row.names(ekman))
plot(epart$conf,type="n")
text(epart0$conf,labels=row.names(ekman))

#test p argument (dimensions)
epart <- powerStressMin (ekman, p=3, verbose = TRUE)
scatterplot3d(epart$conf)
epart <- powerStressMin (ekman, p=3, kappa=0.75, verbose = TRUE)
scatterplot3d(epart$conf)
epart <- powerStressMin (ekman, p=3, kappa=0.5, verbose = TRUE) #not passed!
scatterplot3d(epart$conf)

#test init argument
epart <- powerStressMin (ekman, verbose = 1)
epart
epart2 <- powerStressMin (ekman, verbose = 1, kappa=1.2, init=epart$conf)
epart1 <- powerStressMin (ekman, verbose = 1, kappa=1.2)
if(epart2$itel<epart1$itel) cat("PASSED \n")

test_that("compare obsdiss and so on",{
library(smacof)
library(stops)
data(kinshipdelta)
resp <- powerStressMin(kinshipdelta)
plot(resp$confdiss,resp$obsdiss)

plot(unclass(resp$confdiss),unclass(resp$obsdiss),ylim=c(0,max(resp$confdiss,resp$obsidss)),xlim=c(0,max(resp$confdiss,resp$obsidss)))
})

##testing rstress new version
library(devtools)
load_all()
dis <- smacof::kinshipdelta
resnew <- rStressMin(dis,type="ordinal",r=1)
resold <- rStressMinOld(dis,type="ordinal")
##Ratio works great with or without r. Maybe not show untransformed points as there is not transformation here
## check out the smacof Shepard and resplots
## change the transplot
resnew
resold

plot(Procrustes(resnew$conf,resold$conf))
plot(resnew)


summary(resnew)
summary(resold)

plot.smacofP(resnew,plot.type="transplot")
plotOLD.smacofP(resold,plot.type="transplot")
plot(resnew,plot.type="Shepard")
plotOLD.smacofP(resold,plot.type="Shepard")
smacof:::plot.smacof(resnew,plot.type="Shepard")

plot(resnew,plot.type="resplot")
plot(resnew,plot.type="stressplot")
plot(resnew,plot.type="bubbleplot")
plotOLD.smacofP(resold,plot.type="bubbleplot")


resnew1 <- rStressMin(dis,type="ordinal",r=1)
resnew0.5 <- rStressMin(dis,type="ordinal",r=0.5)
resnew1.5 <- rStressMin(dis,type="ordinal",r=1.5)
resnew0.3 <- rStressMin(dis,type="ordinal",r=0.3)
resnew4 <- rStressMin(dis,type="ordinal",r=4)
resnew0.1 <- rStressMin(dis,type="ordinal",r=0.1)


resnew1 <- rStressMin(dis,type="interval",r=1)
resnew0.5 <- rStressMin(dis,type="interval",r=0.5)
resnew1.5 <- rStressMin(dis,type="interval",r=1.5)
resnew0.3 <- rStressMin(dis,type="interval",r=0.3)

resnew1 <- rStressMin(dis,type="ratio",r=1)
resnew0.5 <- rStressMin(dis,type="ratio",r=0.5)
resnew1.5 <- rStressMin(dis,type="ratio",r=1.5)
resnew0.3 <- rStressMin(dis,type="ratio",r=0.3)

dev.new()
par(mfrow=c(2,3))
plot(resnew0.3,plot.type="Shepard")
plot(resnew0.5,plot.type="Shepard")
plot(resnew1,plot.type="Shepard")
plot(resnew1.5,plot.type="Shepard")
plot(resnew4,plot.type="Shepard")
plot(resnew0.1,plot.type="Shepard")

plot(resnew,plot.type="Shepard")


mean(resnew0.4$confdist/resnew0.5$confdist)

mean(resnew1$confdist/resnew0.5$confdist)

mean(resnew1.5$confdist/resnew0.5$confdist)

mean(resnew1.5$confdist/resnew0.5$confdist)

resnew$dhat/resnew$confdist

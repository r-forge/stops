######## power stress: last passed 2.12.2014
source("rstress.R")
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
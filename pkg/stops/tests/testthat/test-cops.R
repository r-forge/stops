content("COPS tests")

test_that("cop smacof",{
              data(smacof::kinshipdelta)
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
          })


test_that("cop sammon",{
              data(smacof::kinshipdelta)
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
          }

test_that("cop cmdscale",{
              data(smacof::kinshipdelta)
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
              plot(res2$fit$conf)
          })

test_that("cop rstress",{
data(smacof::kinshipdelta)
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
})
          
test_that("cop powerstress",{
data(smacof::kinshipdelta)
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
})

          
test_that("coploss",{
#coploss
res1 <- cop_powerstress(dis,c(1,1),ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=FALSE)
fit <- res1$fit
coploss(fit,a=a,q=q,normed=FALSE,minpts=2,epsilon=10,rang=rang,verbose=1,plot=FALSE,scale=TRUE)
coploss(fit,a=a,q=q,normed=TRUE,minpts=2,epsilon=10,rang=rang,verbose=1,plot=FALSE,scale=TRUE)
coploss(fit,a=a,q=q,normed=TRUE,minpts=3,epsilon=10,rang=rang,verbose=1,plot=FALSE,scale=TRUE)
coploss(fit,a=a,q=q,normed=TRUE,minpts=2,epsilon=0.5,rang=rang,verbose=1,plot=TRUE,scale=TRUE)
coploss(fit,a=a,q=q,normed=TRUE,minpts=2,epsilon=10,rang=c(0,0.5),verbose=1,plot=FALSE,scale=TRUE)
coploss(fit,a=a,q=q,normed=TRUE,minpts=2,epsilon=10,rang=rang,verbose=1,plot=FALSE,scale=FALSE)
})

          
test_that("cops",{
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
plot(optsstress)
lot(optrstress)
plot(optpowerstress)
plot(optpowersammon)
plot(optpowerelastic)
})
          



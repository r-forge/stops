library(stops)

context("COPS tests")

#TODO
#test plots
#test results
#test match call
#test cops arguments

dis <- kinshipdelta
dis <- as.matrix(dis)
test1 <- cops(dis,verbose=4)

test_that("cops runs in default mode",{
#cops
expect_equal_to_reference(test1 <- cops(dis))
})

test_that("cops loss argument work right",{
        test2 <- cops(dis,loss="strain")
        expect_that(test2$fit,is_a("cmdscale"))
        test3 <- cops(dis,loss="stress")
        expect_that(test3$fit,is_a("smacofB"))
        test4 <- cops(dis,loss="powerelastic")
        expect_that(test4$fit,is_a("smacofP"))
        expect_that(test4$fit$weightmat,equals(dis^(-2*test4$par[2])))
        test5 <- cops(dis,loss="rstress")
        expect_that(test5$fit,is_a("smacofP"))            
        test6 <- cops(dis,loss="sammon")
        expect_that(test6$fit,is_a("cmdscale"))            
        test7 <- cops(dis,loss="powersammon",verbose=4)
        expect_that(test6$fit,is_a("smacofP"))
        expact_that(test6$fit$weightmat,equals(dis^(-test6$par[2])))
        test8 <- cops(dis,loss="smacofSphere")
        expect_that(test8$fit,is_a("smacofS"))
        test9 <- cops(dis,loss="elastic")
        expect_that(test9$fit,is_a("smacofB"))
        test10 <- cops(dis,loss="smacofSym")
        expect_that(test10$fit,is_a("smacofB"))
        test11 <- cops(dis,loss="sstress")
        expect_that(test11$fit,is_a("smacofP"))
        test12 <- cops(dis,loss="powermds")
        expect_that(test12$fit,is_a("smacofP"))
        test13 <- cops(dis,loss="powerstress",verbose=4)
        expect_that(test13$fit,is_a("smacofP"))
    })

test_that("cops ndim argument",{
        #ndim      
        test12 <- cops(dis,loss="strain",ndim=3)
        expect_that(dim(test12$fit$points)[2],equals(3))
    })


        
test_that("cops q minpts ndim arguments",{ 
        test13 <- cops(dis,loss="stress",ndim=3,q=2,minpts=3)
        expect_that(test13$OC$minpts,equals(3))
        expect_that(test13$OC$q,equals(2))
        expect_that(dim(test13$fit$conf)[2]),equals(3))
})

test_that("cops epsilon parameter"{
        #epsilon  
        test14 <- cops(dis,loss="strain",ndim=2,q=2,epsilon=2)
        expect_that(test13$OC$epsilon),equals(2))
})
        #verbose 
test_that("cops verbose parameter",{
        expect_that(cops(dis,loss="strain",verbose=0),!prints_text("*"))
        expect_output(cops(dis,loss="strain",verbose=1))
})


test_that("cops plots parameter",{
        #plot
        test12 <- cops(dis,loss="strain",ndim=2,q=2,epsilon=2,plot=TRUE)
})

test_that("cops scale parameter",{
        #scale
        test13 <- cops(dis,loss="strain",ndim=2,q=2,epsilon=2,scale=TRUE)
        test13a <- cops(dis,loss="strain",ndim=2,q=2,epsilon=2,scale=FALSE)
        expect_that(test13,!identical(test13a))
        expect_that(test13$fit$conf,equals(scale(test13a$fit$conf)))
    })
      

test_that("Initsol for cops",{
               initsol <- powerStressMin(dis,kappa=1,lambda=1,verbose=0)
               expect_equal_to_reference(initsol)        
            })

test_that("Powerstress COPS Kinship",{
              teso <- cops(dis,loss="powerstress",q=q,minpts=minpts,epsilon=eps,rang=rang,verbose=0,plot=FALSE,scale=scale,lower=c(0.75,0.75,0.5),upper=c(4,4,2),itmax)
             #expect class
             #expect a value
              expect_equal_to_reference(teso)
          })

test_that("Finding cordweight",{
         dis <- as.matrix(smacof::kinshipdelta)
         rang <- c(0,1.2) #fixed range 
         minpts <- 2
         eps <- 10
         q <- 1
         initcorrd <- cordillera(initsol$conf,rang=rang,q=q,minpts=minpts,plot=TRUE,scale=scale)
         a <- initsol$stress.m/initcorrd$normed
         expect_equal(teso$cordweight,a)
     })

test_that("cops transplot had side effect",{
     testp <- cops(dis,loss="stress",theta=c(1,2),itmax=1)
     testp2 <- smacofSym(dis^2)
     expect_that(testp$fit[1:6],equals(testp2[1:6])) 

     plot(testp,"Shepard") #no side effect
     expect_that(testp$fit[1:6],equals(testp2[1:6]))

     plot(testp,"transplot") #had the side effect of changing the smacof object in testp$fit
     expect_that(testp$fit[1:6],equals(testp2[1:6])) #was not fullfilled before bug was fixed

     testp <- cops(dis,loss="stress",theta=c(2,2),itmax=1)
     testp2 <- powerstressMin(dis,kappa=2,lambda=2)
     expect_that(testp$fit[1:6],equals(testp2[1:6]))

     plot(testp,"Shepard") #no side effect
     expect_that(testp$fit[1:6],equals(testp2[1:6]))

     plot(testp,"transplot") #had the side effect of changing the smacof object in testp$fit
     expect_that(testp$fit[1:6],equals(testp2[1:6])) #was not fullfilled before bug was fixed
  
)}

test_that("cops transplot equals powerstress transplot",{
     testp <- cops(dis,loss="stress",theta=c(2,2),itmax=1)
     testp3 <- powerstressMin(dis,kappa=2,lambda=2) 
     plot(testp,"transplot") #how to test if two plots are the same?
     plot(testp3,"transplot")
})

## library(MASS)
## library(stops)
## data(BankingCrisesDistances)
## set.seed(210485)
## opto <- cops(BankingCrisesDistances[,1:69],theta=1,loss="sammon",verbose=2,acc=1e-16,accd=1e-12,cordweight=0.5,lower=0.5,upper=5)

## opto <- cops(BankingCrisesDistances[,1:69],theta=1,loss="sammon",verbose=2,acc=1e-12,accd=1e-12,lower=0.5,upper=5,scale=FALSE)

## library(stops)
## data(BankingCrisesDistances)
## set.seed(210485)
## optstrain <- cops(BankingCrisesDistances[,1:69],loss="strain",verbose=1)
## optstress <- cops(BankingCrisesDistances[,1:69],loss="stress",verbose=3)
## optsammon <- cops(BankingCrisesDistances[,1:69],loss="sammon",verbose=3)
## optelastic <- cops(BankingCrisesDistances[,1:69],loss="elastic",verbose=3)
## optsstress <- cops(BankingCrisesDistances[,1:69],loss="sstress",verbose=3)
## optrstress <- cops(BankingCrisesDistances[,1:69],loss="rstress",verbose=3)
## optpowerstress <- cops(BankingCrisesDistances[,1:69],loss="powerstress",verbose=3)
## optpowersammon <- cops(BankingCrisesDistances[,1:69],loss="powersammon",verbose=3)
## optpowerelastic <- cops(BankingCrisesDistances[,1:69],loss="powerelastic",verbose=3)

## optstrain
## optstress
## optsammon
## optelastic
## optsstress
## optrstress
## optpowerstress
## optpowersammon
## optpowerelastic



## plot(optstrain)
## plot(optstress)
## plot(optsammon)
## plot(optelastic)
## plot(optsstress)
## lot(optrstress)
## plot(optpowerstress)
## plot(optpowersammon)
## plot(optpowerelastic)
## })


## test_that("cop smacof",{
##               data(smacof::kinshipdelta)
##               dis <- kinshipdelta
##               #defaults
##               res1 <- cop_smacof(dis)
##               #verbose
##               res1 <- cop_smacof(dis,verbose=1)
##               res1 <- cop_smacof(dis,verbose=2)
##                                         #plot and single theta argument, ndim argument
##               res2 <- cop_smacof(dis,3,ndim=3,a=0.5,q=1,minpts=2,epsilon=10,plot=TRUE)
##                                         #vector theta argument
##               res3 <- cop_smacof(dis,c(2,3),ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE) #the kappa value gets ignored as it should be
##               res3
##                                         #cordillera meta parameters: a, q, minpts, epsilon and rang
##               res2 <- cop_smacof(dis,3,ndim=2,a=0.1,q=2,minpts=3,epsilon=1,verbose=0,plot=TRUE)
##               res2
##               res2a <- cop_smacof(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE)
##               res2a
##               res2b <- cop_smacof(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=c(0,1),plot=TRUE)
##               res2b
##               plot(res2a$fit)
##               plot(res2$fit)
##                                         #weightmat and init
##               w <- 1-diag(nrow(dis))
##               w[c(1,2),c(1,2)] <- 0
##               res <- cop_smacof(dis,3)
##               res2 <- cop_smacof(dis,3,weightmat=w)
##               res
##               res2
##               res2a <- cop_smacof(dis,3,init=res2$fit$conf,weightmat=w)
##               res2a

##               #in optimization
##               rang <- c(0,1.45)
##               smacop <- ljoptim(1,function(lambda) cop_smacof(dis,lambda,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)$copstress,lower=c(0.5,0.5),upper=c(5,5))
##               smacop
##               cop_smacof(dis,smacop$par,ndim=2,a=0.5,q=1,minpts=2,eps=10,verbose=0,rang=rang,plot=TRUE)
##               cop_smacof(dis,4.9,ndim=2,a=0.5,q=1,minpts=2,eps=10,verbose=0,rang=rang,plot=TRUE)
##           })


## test_that("cop sammon",{
##               data(smacof::kinshipdelta)
##               dis <- kinshipdelta
##                                         #defaults
##               res1 <- cop_sammon(dis)
##               res1
##                                         #verbose
##               res1 <- cop_sammon(dis,verbose=1)
##               res1 <- cop_sammon(dis,verbose=2)
##                                         #plot and single theta argument, ndim argument
##               res2 <- cop_sammon(dis,3,ndim=3,a=0.5,q=1,minpts=2,epsilon=10,plot=TRUE)
##                                         #vector theta argument
##               res3 <- cop_sammon(dis,c(2,3),ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE) #the kappa value gets ignored as it should be
##               res3
##                                         #cordillera meta parameters: a, q, minpts, epsilon and rang
##               res2 <- cop_sammon(dis,3,ndim=2,a=0.1,q=2,minpts=3,epsilon=1,verbose=0,plot=TRUE)
##               res2
##               res2a <- cop_sammon(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE)
##               res2a
##               res2b <- cop_sammon(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=c(0,1),plot=TRUE)
##               res2b
##               plot(res2a$fit$conf)
##               plot(res2$fit$conf)
##                                         #weightmat and init
##               res2 <- cop_sammon(dis,3,verbose=1)
##               res2
##               res2a <- cop_sammon(dis,3,init=res2$fit$conf,verbose=1)
##               res2a
          

##                                         #in optimization
##               rang <- c(0,1.45)
##               samcop <- ljoptim(1,function(lambda) cop_sammon(dis,lambda,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)$copstress,lower=c(0.5,0.5),upper=c(5,5))
##               samcop
##               cop_sammon(dis,smacop$par,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)
##               cop_sammon(dis,1,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)
##           }

## test_that("cop cmdscale",{
##               data(smacof::kinshipdelta)
##               dis <- kinshipdelta
##                                        #defaults
##               res1 <- cop_cmdscale(dis)
##               res1
##                                         #verbose
##               res1 <- cop_cmdscale(dis,verbose=1)
##               res1 <- cop_cmdscale(dis,verbose=2)
##                                         #plot and single theta argument, ndim argument
##               res2 <- cop_cmdscale(dis,3,ndim=3,a=0.5,q=1,minpts=2,epsilon=10,plot=TRUE)
##                                         #vector theta argument
##               res3 <- cop_cmdscale(dis,c(2,3),ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE) #the kappa value gets ignored as it should be
##               res3
## #cordillera meta parameters: a, q, minpts, epsilon and rang
##               res2 <- cop_cmdscale(dis,3,ndim=2,a=0.1,q=2,minpts=3,epsilon=1,verbose=0,plot=TRUE)
##               res2
##               res2a <- cop_cmdscale(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE)
##               res2a
##               res2b <- cop_cmdscale(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=c(0,1),plot=TRUE)
##               res2b
##               plot(res2a$fit$conf)
##               plot(res2$fit$conf)
##                                         #weightmat and init
##               res2 <- cop_cmdscale(dis,3,verbose=1)
##               res2
##               res2a <- cop_cmdscale(dis,3,init=res2$fit$conf,verbose=1)
##               res2a

##                                         #in optimization
##               rang <- c(0,1.45)
##               cmdscalecop <- ljoptim(1,function(lambda) cop_cmdscale(dis,lambda,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)$copstress,lower=c(0.5,0.5),upper=c(5,5))
##               cmdscalecop
##               res1 <- cop_cmdscale(dis,cmdscalecop$par,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)
##               res2 <- cop_cmdscale(dis,1,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)
##               plot(res1$fit$conf)
##               plot(res2$fit$conf)
##           })

## test_that("cop rstress",{
## data(smacof::kinshipdelta)
## dis <- kinshipdelta
## #defaults
## res1 <- cop_rstress(dis)
## #verbose
## res1 <- cop_rstress(dis,verbose=1)
## res1 <- cop_rstress(dis,verbose=2)
## #plot and single theta argument, ndim argument
## res2 <- cop_rstress(dis,3,ndim=3,a=0.5,q=1,minpts=2,epsilon=10,plot=TRUE)
## #vector theta argument
## res3 <- cop_rstress(dis,c(2,3),ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE) #the lambda value gets ignored as it should be
## res3
## #cordillera meta parameters: a, q, minpts, epsilon and rang 
## res2 <- cop_rstress(dis,3,ndim=2,a=0.1,q=2,minpts=3,epsilon=1,verbose=0,plot=TRUE)
## res2
## res2a <- cop_rstress(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE)
## res2a
## res2b <- cop_rstress(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=c(0,1),plot=TRUE)
## res2b
## plot(res2x$fit$conf)
## dev.new()
## plot(res2$fit$conf)
## res2 <- cop_rstress(dis,3,ndim=2,a=0.1,q=1,minpts=4,epsilon=10,verbose=0,plot=TRUE)
## res2
## res2x <- fpowerStress(dis,3,ndim=2,a=0.1,q=1,minpts=4,epsilon=10,verbose=0,plot=TRUE)
## #weightmat and init
## w <- 1-diag(nrow(dis))
## w[c(1,2),c(1,2)] <- 0
## res <- cop_rstress(dis,3)
## res2 <- cop_rstress(dis,3,weightmat=w)
## res
## res2
## res2a <- cop_rstress(dis,3,init=res2$fit$conf,weightmat=w) 
## res2a

## #in optimization
## rang <- c(0,1.45)
## powercop <- ljoptim(c(1,1),function(theta) cop_rstress(dis,theta,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)$copstress,lower=c(1,0.5),upper=c(5,5))
## powercop
## resopt <- cop_rstress(dis,powercop$par,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)
## res1 <- cop_rstress(dis,c(1,1),ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)

## res1
## resopt
## plot(resopt$fit$conf)
## plot(res1$fit$conf)
## })
          
## test_that("cop powerstress",{
## data(smacof::kinshipdelta)
## dis <- kinshipdelta
## #defaults
## res1 <- cop_powerstress(dis)
## #verbose
## res1 <- cop_powerstress(dis,verbose=1)
## res1 <- cop_powerstress(dis,verbose=2)
## #plot and single theta argument, ndim argument
## res2 <- cop_powerstress(dis,3,ndim=3,a=0.5,q=1,minpts=2,epsilon=10,plot=TRUE)
## #vector theta argument
## res3 <- cop_powerstress(dis,c(2,3),ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE) #the kappa value gets ignored as it should be
## res3
## #cordillera meta parameters: a, q, minpts, epsilon and rang 
## res2 <- cop_powerstress(dis,3,ndim=2,a=0.1,q=2,minpts=3,epsilon=1,verbose=0,plot=TRUE)
## res2
## res2a <- cop_powerstress(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,plot=TRUE)
## res2a
## res2b <- cop_powerstress(dis,3,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=c(0,1),plot=TRUE)
## res2b
## plot(res2a$fit$conf)
## plot(res2$fit$conf)
## #weightmat and init
## w <- 1-diag(nrow(dis))
## w[c(1,2),c(1,2)] <- 0
## res <- cop_powerstress(dis,3)
## res2 <- cop_powerstress(dis,3,weightmat=w)
## res
## res2
## res2a <- cop_powerstress(dis,3,init=res2$fit$conf,weightmat=w) 
## res2a

## #in optimization
## rang <- c(0,1.45)
## powercop <- ljoptim(c(1,1),function(theta) cop_powerstress(dis,theta,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)$copstress,lower=c(1,0.5),upper=c(5,5))
## powercop
## resopt <- cop_powerstress(dis,powercop$par,ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)
## res1 <- cop_powerstress(dis,c(1,1),ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=TRUE)

## res1
## resopt
## plot(resopt$fit$conf)
## plot(res1$fit$conf)
## })

          
## test_that("coploss",{
## #coploss
## res1 <- cop_powerstress(dis,c(1,1),ndim=2,a=0.5,q=1,minpts=2,epsilon=10,verbose=0,rang=rang,plot=FALSE)
## fit <- res1$fit
## coploss(fit,a=a,q=q,normed=FALSE,minpts=2,epsilon=10,rang=rang,verbose=1,plot=FALSE,scale=TRUE)
## coploss(fit,a=a,q=q,normed=TRUE,minpts=2,epsilon=10,rang=rang,verbose=1,plot=FALSE,scale=TRUE)
## coploss(fit,a=a,q=q,normed=TRUE,minpts=3,epsilon=10,rang=rang,verbose=1,plot=FALSE,scale=TRUE)
## coploss(fit,a=a,q=q,normed=TRUE,minpts=2,epsilon=0.5,rang=rang,verbose=1,plot=TRUE,scale=TRUE)
## coploss(fit,a=a,q=q,normed=TRUE,minpts=2,epsilon=10,rang=c(0,0.5),verbose=1,plot=FALSE,scale=TRUE)
## coploss(fit,a=a,q=q,normed=TRUE,minpts=2,epsilon=10,rang=rang,verbose=1,plot=FALSE,scale=FALSE)
## })


#'\donttest{
#'# From De Leuuw et al (2016) example 7.2.
#'#They look at different rstress versions and compare how clustered the configuration is
#'#where stress is minimal and that stress is a monotonically increasing function of r;
#' dats <- c(5.63,5.27, 6.72,4.60, 5.64, 5.46,4.80, 6.22, 4.97, 3.20,7.54 ,5.12, 8.13, 7.84 ,7.80, 6.73 ,4.59 ,7.55, 6.73, 7.08, 4.08, 7.18 ,7.22 ,6.90 ,7.28 ,6.96 ,6.34 ,6.88, 6.17, 5.47, 4.67, 6.13, 6.04 ,7.42, 6.36, 7.36)
#'num_cols <- (1 + sqrt(1 + 8*length(dats)))/2 - 1
#'mat <- matrix(0, num_cols, num_cols)
#'mat[row(mat) <= col(mat)] <- dats
#'mat <- t(mat)
#'mat <- rbind(0, mat)
#'mat <- cbind(mat, 0)
#'colnames(mat) <- rownames(mat) <- c(" KVP", "PvdA" , "VVD" , "ARP" , "CHU" , "CPN" , "PSP" ,  "BP", "D66")
#'dobj <- as.dist(mat)
#'dobj
#'#We can do this in one go by setting cordweight to 0 and find that stress is minimal (0.0033) around r~=0.17 (kappa~=0.34)
#'#and that stress appears thus not monotonically increasing in r
#' set.seed(210485)
#' m1 <- pcops(dobj,loss="rstress",lower=c(0.05,1,1),upper=c(5,1,1),verbose=3,cordweight=0,stressweight=1)
#' m1
#'# They observe increasing clustering for larger r which we can again do systematically:
#'# When only clusteredness is of interest, we use cordweight=1 stressweight=0 and try clusters of at least k=2 and k=3 observations
#' set.seed(210485)
#' m2 <- pcops(dobj,loss="rstress",minpts=2,lower=c(0.05,1,1),upper=c(5,1,1),verbose=3,cordweight=1,stressweight=0) 
#' m3 <- pcops(dobj,loss="rstress",minpts=3,lower=c(0.05,1,1),upper=c(5,1,1),verbose=3,cordweight=1,stressweight=0)
#' m2   #r~=1.24
#' m3   #r~=1.39
#'
#'# It is generally better to trade off clusteredness and fit
#' set.seed(210485)
#' m2t <- pcops(dobj,loss="rstress",minpts=2,theta=c(m1$par[1],1,1),lower=c(0.05,1,1),upper=c(5,1,1),verbose=3,cordweight=1/3,stressweight=2/3)
#' m3t <- pcops(dobj,loss="rstress",minpts=3,theta=c(m1$par[1],1,1),lower=c(0.05,1,1),upper=c(5,1,1),verbose=3,cordweight=1/3,stressweight=2/3)
#' m2t #r~=0.08
#' m4t #r~=1.39
#'}


delta <- kinshipdelta
disobj <- smacof::transPrep(as.dist(delta), trans = "none", spline.intKnots = 2, spline.degree = 2)
r <-0.5
n <- dim(delta)[1]
ndim <- 2
weightmat <- 1-diag(nrow(delta))
stressweight <- 1
cordweight <- 0
q <- 1
minpts <- 2
epsilon <- 1.2
rang <- c(0,5)
scalo <- "std"
normed <- TRUE
init <- NULL
kappa <- 1
lambda <- 1
nu <- 1
xold  <- init <- cops::powerStressFast(delta,kappa=kappa,lambda=lambda,nu=nu,ndim=ndim)$conf
stresstype <- "stress-1"
optimized1 <- dfoptim::hjk(xold,function(par) copsf(par,delta=delta,disobj=disobj,r=r,n=n,ndim=ndim,weightmat=weightmat,stressweight=stressweight,cordweight=cordweight,q=q,minpts=minpts,epsilon=epsilon,rang=rang,scalo=scalo,normed=normed,init=init),control=list(maxfeval=itmax)) #,...)


library(nloptr)
?newuoa

 fr <- function(x) {   ## Rosenbrock Banana function
         100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
 }

itel <- 5
itelmax <- 10
iteli <- itelmax-itel

(S <- newuoa(c(1, 2), fr,control=list(maxeval=itelmax-itel)))


res1<-copstressMin(dis,stressweight=0.95,cordweight=0.05,itmax=10000) #use higher itmax about 10000
res1


##### Trying out all losses in pcops with the correct for theta and upper and lower
library(cops)
dis <- smacof::kinshipdelta #df
dis <- cops::matchphi #matrix
dis <- smacof::morse #dist object
#strain
p1 <- pcops(dis,loss="strain",type="ratio",theta=1,lower=0,upper=5,verbose=3)
p1
summary(p1)
plot(p1,"confplot")
plot(p1,"reachplot")
plot(p1,"Shepard")
plot(p1,"transplot")
#stress
p1 <- pcops(dis,loss="stress",type="ratio",theta=1,lower=0,upper=5,verbose=3)
p1
summary(p1)
plot(p1,"confplot")
plot(p1,"reachplot")
plot(p1,"Shepard") 
plot(p1,"transplot") 
#sammon 
p1 <- pcops(dis,loss="sammon",type="ratio",theta=1,lower=0,upper=5,verbose=3)
p1
summary(p1)
plot(p1,"confplot")
plot(p1,"reachplot")
plot(p1,"Shepard")
plot(p1,"transplot")
#sammon2
p1 <- pcops(dis,loss="sammon2",type="interval",theta=1,lower=0,upper=4,verbose=3)
p1
summary(p1)
plot(p1,"confplot")
plot(p1,"reachplot")
plot(p1,"Shepard")
plot(p1,"transplot")
#rstress
p1 <- pcops(dis,loss="rstress",type="ordinal",theta=1,lower=0,upper=5,verbose=3)
p1
summary(p1)
plot(p1,"confplot")
plot(p1,"reachplot")
plot(p1,"Shepard")
plot(p1,"transplot")
#elastic ##CHECK again why such a large confdist
p1 <- pcops(dis,loss="elastic",type="interval",theta=1,lower=0.5,upper=3,verbose=3)
p1
summary(p1)
plot(p1,"confplot")
plot(p1,"reachplot")
plot(p1,"Shepard")
plot(p1,"transplot")
#powerstrain #same as strain
p1 <- pcops(dis,loss="powerstrain",type="ratio",theta=1,lower=0,upper=5,verbose=3)
p1
summary(p1)
plot(p1,"confplot")
plot(p1,"reachplot")
plot(p1,"Shepard")
plot(p1,"transplot")
#sstress
p1 <- pcops(dis,loss="sstress",type="ratio",theta=1,lower=0,upper=5,verbose=3)
p1
summary(p1)
plot(p1,"confplot")
plot(p1,"reachplot")
plot(p1,"Shepard")
plot(p1,"transplot")
#powermds
p1 <- pcops(dis,loss="powermds",type="interval",theta=c(1,1),lower=c(0,0),upper=c(5,5))
p1
summary(p1)
plot(p1,"confplot")
plot(p1,"reachplot")
plot(p1,"Shepard")
plot(p1,"transplot")
#powersammon
p1 <- pcops(dis,loss="powersammon",type="interval",theta=c(1,1),lower=c(0,0),upper=c(5,5))
p1
summary(p1)
plot(p1,"confplot")
plot(p1,"reachplot")
plot(p1,"Shepard")
plot(p1,"transplot")
#powerelastic
p1 <- pcops(dis,loss="powerelastic",type="interval",theta=c(1,1),lower=c(0.1,0.1),upper=c(3,3),verbose=3)
p1
summary(p1)
plot(p1,"confplot")
plot(p1,"reachplot")
plot(p1,"Shepard")
plot(p1,"transplot")
#rpowerstress
p1 <- pcops(dis,loss="rpowerstress",type="interval",theta=c(1,1),lower=c(0,0),upper=c(5,5))
p1
summary(p1)
plot(p1,"confplot")
plot(p1,"reachplot")
plot(p1,"Shepard")
plot(p1,"transplot")
#apstress
p1 <- pcops(dis,loss="apstress",type="ratio",theta=c(1,1,1),lower=c(1,1,1),upper=c(3,3,3),verbose=3)
p1
summary(p1)
plot(p1,"confplot")
plot(p1,"reachplot")
plot(p1,"Shepard")
plot(p1,"transplot")
#powerstress
p1 <- pcops(dis,loss="powerstress",type="interval",theta=c(1,1,1),lower=c(0,0,0),upper=c(5,5,5),verbose=3)
p1
summary(p1)
plot(p1,"confplot")
plot(p1,"reachplot")
plot(p1,"Shepard")
plot(p1,"transplot")




                                        #try all 2nd kinship 3rd matchphi
losses <- c("stress",#1 ok ok
            "smacofSym",#1 ok ok
            "strain",#1 ok ok
            "rstress", #1 ok ok
          #  "elastic",#1 ok ok
            "apstress", #3  ok
losses <- c(           "sammon2",#1 ok ok
      "powerstrain",#1 ok ok
      "sammon", #1 ok ok
      "powerelastic",#2 ok ok wrong
      "powerstress",#3 ok ok
      "sstress", #1 ok ok
      #)  
#works with 3 works with 1 but sometimes issues with 2
#why only one value by the wrong losses?  
#losses <-     c(
    "powermds", #2 wrong, worked, wrong, worked ok     
      "powersammon",#2 wrong, worked, worked ok
"rpowerstress" #2 wrong, worked ok
  )
      
 
  

dis <- as.matrix(dis)
for(i in losses)
{
 cat("Testing:",i,"\n")   
 p1 <- pcops(dis,loss=i,type="ratio",lower=c(0.1,0.1,0.1),upper=c(5,5,5),verbose=2,itmaxi=5000)
 print(p1)
 summary(p1)
 par(mfrow=c(2,1))
 plot(p1,"confplot",label.conf=list(pos=5))
 plot(p1,"reachplot") 
# plot(p1,"Shepard")
 #plot(p1,"transplot") transplot for stress
 par(mfrow=c(1,1))
}

## elastic has issues, apstress has issues (related to the upper agument)
## what is sstress?
## transplot for smacofB issues, sheaprd issues for apstress

#rstessmin and powermds has only one parameter! what gives?

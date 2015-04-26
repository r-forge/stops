content("OPTICS cordillera")

              q <- 1
              eps <- 2
             
test_that("cordillera works",{
          data(iris)
          expect_equal_to_reference(res<-cordillera(iris[,1:4],minpts=2,epsilon=100))
          expect_that(print(res),not(throws_error()))
        #  expect_that(print(res),matches("        raw      normed  16.51312758  0.07134626",fixed=TRUE))
          expect_that(summary(res),not(throws_error()))
          #expect_that(summary(res),matches("OPTICS cordillera values with minpts= 2 and epsilon= 100",fixed=TRUE)) 
          expect_that(plot(res),not(throws_error()))
          expect_that(plot(res,withlabels=TRUE),gives_warning())
  }

 test_that("cordillera arguments work",{
          c1 <- cordillera(confs)
          c1
          c2 <- cordillera(confs,epsilon=eps)
          c2
          c3 <- cordillera(confs,q=q,plot=TRUE)
          c3
          c4 <- cordillera(confs,minpts=3)
          c4
          c5 <- cordillera(confs,rang=c(0,10),minpts=3,plot=TRUE)
          c5
          c6 <- cordillera(confs,rang=c(0,10),minpts=3,ylim=c(0,1.5),plot=FALSE)
          c6
          c7 <- cordillera(confs,digits=1,minpts=3,plot=FALSE)
          c7
          c8 <- cordillera(confs,digits=1,minpts=3,path="~")
          c8
      })         

q <- 1
eps <- 2
range <- c(0,0.9225246)
test_that("cordillera calculates correctly",{
              #from the paper
library(stops)
#row 1
x <- c(-0.375,-0.125,0.125,0.375,-0.375,-0.125,0.125,0.375)
x <- x*2
y <- c(-0.125,-0.125,-0.125,-0.125,0.125,0.125,0.125,0.125)
y <- y*2
c1 <- cordillera(cbind(x,y),epsilon=eps,q=q,rang=range,scale=FALSE)
expect_that(c1$raw,equals(0))
expect_that(c1$normed,equals(0))
expect_that(c1$normfac,equals(8))
expect_that(c1$dmax,equals(max(range)))

#row 2
obj <- 8
set.seed(210485)
mat <- matrix(rep(1,obj^2),ncol=obj)+rnorm(64,sd=0.0001)
diag(mat) <- 0
mds <- smacofSym(mat)
mdsmat <-as.matrix(mds$conf)
confdiss <- as.matrix(mds$confdiss)
c2 <- cordillera(mdsmat,epsilon=eps,q=q,rang=range,scale=FALSE)
expect_that(c2$raw,equals(0.16612214))
expect_that(c2$normed,equals(0.0257247712))
expect_that(c2$normfac,equals(8))
expect_that(c2$dmax,equals(max(range)))

#row 3
mdsmat <- cbind(rep(c(0,0,-0.5,0.5),2),rep(c(0.5,-0.5,0,0),2))
set.seed(123)
mdsmat <- jitter(mdsmat,factor=3.5) 
c3 <- cordillera(mdsmat,epsilon=eps,q=q,rang=range,scale=FALSE)
expect_that(c3$raw,equals(1.46845612))
expect_that(c3$normed,equals(0.22739712926277))
expect_that(c3$normfac,equals(8))
expect_that(c3$dmax,equals(max(range)))

#row 4
mdsmat <- cbind(rep(c(0,0,-0.5,0.5),2),rep(c(0.5,-0.5,0,0),2))
set.seed(123)
mdsmat <- jitter(mdsmat,factor=1.7) 
c4 <- cordillera(mdsmat,epsilon=eps,q=q,rang=range,scale=FALSE)
expect_that(c4$raw,equals(2.35424102))
expect_that(c4$normed,equals(0.3645649619688))
expect_that(c4$normfac,equals(8))
expect_that(c4$dmax,equals(max(range)))

#row 5
mdsmat <- cbind(rep(c(0,0,-0.5,0.5),2),rep(c(0.5,-0.5,0,0),2))
set.seed(210485)
mdsmat <- jitter(mdsmat,factor=1) 
c5 <- cordillera(mdsmat,epsilon=eps,q=q,rang=range,scale=FALSE)
expect_that(c5$raw,equals(3.92984969))
expect_that(c5$normed,equals(0.60855515242783))
expect_that(c5$normfac,equals(8))
expect_that(c5$dmax,equals(max(range)))


#row 6
ind <- matrix(c(1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,
                0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,
                0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1),ncol=8)
fake <- matrix(1,nrow=8,ncol=8)
dissi <- fake*(1-ind)
exdiss <- as.dist(dissi)
exmds <- smacofSym(exdiss)
mdsmat <- exmds$conf
c6 <- cordillera(mdsmat,q=q,rang=range,epsilon=eps,scale=FALSE)
expect_that(c6$raw,equals(6.4530180))
expect_that(c6$normed,equals(0.999279275))
expect_that(c6$normfac,equals(8))
expect_that(c6$dmax,equals(max(range)))

#use the original length; in cordillera there is some rounding, so perhaps undo the rounding?  
range2 <- c(0,0.92136105)
c6 <- cordillera(mdsmat,q=q,rang=range2,epsilon=eps,scale=FALSE)
expect_that(c6$normed,equals(1))
})


test_that("some comparisons",{
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
})

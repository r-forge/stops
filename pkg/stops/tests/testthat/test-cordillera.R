content("OPTICS cordillera")

test_that("cordillera works",{
              library(MASS)
              par(mfrow=c(3,2))
              q <- 1
              eps <- 2
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
          })


test_that("cordillera tests",{
          data(smacof::kinshipdelta)
          dis <- kinshipdelta
          rescoy <- smacofSym(dis^lamb, itmax = 5000)
          confs <- rescoy$conf #make matrix from fitted distances; the diagonal is zero
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

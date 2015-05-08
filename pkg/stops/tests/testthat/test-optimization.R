context("Optimization")

fbana <- function(x) {
 x1 <- x[1]
 x2 <- x[2]
 100 * (x2 - x1 * x1)^2 + (1 - x1)^2
 }

fwild <- function (x) 10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80

test_that("ALJ works",{
              set.seed(210485)
              res1<-ljoptim(c(-1.2,1),fbana,lower=-5,upper=5,accd=1e-16,acc=1e-16)
              res1
              set.seed(210485)
              res2<-ljoptim(50, fwild,lower=-50,upper=50,adaptive=FALSE,accd=1e-16,acc=1e-16)
              res2  
          })

test_that("ALJ does not stop for one d of length 0",{
       res1 <- ljoptim(c(-1.2,1),fbana,lower=c(-5,1),upper=c(5,1),accd=1e-16,acc=1e-16)
       expect_that(res1$counts[1],not(is_less_than(2)))
       })



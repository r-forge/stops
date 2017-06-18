context("optics")

test_that("OPTICS works",{
         data(smacof::kinshipdelta)
         dis <- kinshipdelta
         rescoy <- smacof::smacofSym(dis, itmax = 5000)
         confs <- rescoy$conf #make matrix from fitted distances; the diagonal is zero
         test1 <- optics(confs,10,2)
         test1
         test2 <- optics(confs,minpts=3)
         test2
         test3 <- optics(confs,minpts=3,keep=TRUE)
         test3
     })

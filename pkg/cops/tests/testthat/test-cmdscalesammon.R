context("cmdscale and sammon")

library(stops)

initsam <- sammon(dis^3)
expect_that(initsam,is_a("sammon"))
expect_that(initsam,is_a("cmdscale"))

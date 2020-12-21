library(testthat)
library(SRRMLR)

test_that("input data obeys our standard",{
  expect_equal(class(HAM$hamimgs),"matrix")
  expect_equal(class(HAM$hamlab),"matrix")
  expect_equal(dim(HAM$hamimgs)[1],dim(HAM$hamlab)[1])
  expect_equal(rep(1,dim(HAM$hamlab)[1]), rowSums(HAM$hamlab))
})

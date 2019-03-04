context("SBM test")
library(greed)

test_that("SBM sim", {
  N = 500
  K = 10
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))
  sbm = rsbm(N,pi,mu)
  expect_equal(dim(sbm$x)[1], N)
  expect_equal(dim(sbm$x)[2], N)
  expect_equal(length(sbm$cl),N)
  expect_gte(min(sbm$cl), 1)
  expect_lte(max(sbm$cl), K)
})

context("DIAGGMM test")
library(greed)
library(ggplot2)
set.seed(1234)




test_that("DIAGGMM hybrid", {
  N=150
  X = rbind(MASS::mvrnorm(N/3,c(-5,0),diag(2)),MASS::mvrnorm(N/3,c(0,5),diag(2)),MASS::mvrnorm(N/3,c(5,0),diag(2)))
  mod = DiagGmm(mu=apply(X,2,mean))
  sol=greed(X,model=mod)
  expect_equal(sol@K, 3)
  
  co=coef(sol)
  expect_equal(nrow(do.call(rbind,co$muk)),3)
  expect_equal(ncol(do.call(rbind,co$muk)),2)
  expect_equal(sum(co$pi),1)
  expect_equal(length(co$pi),3)
  
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(sol,type='tree')))
  expect_true(is.ggplot(plot(sol,type='path')))
  expect_true(is.ggplot(plot(sol,type='front')))
  expect_true(is(plot(sol),"gtable"))
  expect_true(is(plot(sol,type="violins"),"gtable"))
})

test_that("DiagGMM seed", {
  N=150
  X = rbind(MASS::mvrnorm(N/3,c(-5,0),diag(2)),MASS::mvrnorm(N/3,c(0,5),diag(2)),MASS::mvrnorm(N/3,c(5,0),diag(2)))
  mod = DiagGmm()
  sol=greed(X,model=mod,alg=Seed())
  expect_equal(sol@K, 3)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(sol,type='tree')))
  expect_true(is.ggplot(plot(sol,type='path')))
  expect_true(is.ggplot(plot(sol,type='front')))
  expect_true(is(plot(sol),"gtable"))
  expect_true(is(plot(sol,type="violins"),"gtable"))
})



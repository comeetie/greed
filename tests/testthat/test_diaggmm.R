context("DIAGGMM test")
library(greed)
library(ggplot2)
set.seed(1234)




test_that("DIAGGMM hybrid", {
  N=150
  X = rbind(MASS::mvrnorm(N/3,c(-5,0),diag(2)),MASS::mvrnorm(N/3,c(0,5),diag(2)),MASS::mvrnorm(N/3,c(5,0),diag(2)))
  mod = new("diaggmm",mu=apply(X,2,mean))
  sol=greed(X,model=mod)
  expect_equal(sol@K, 3)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(sol,type='tree')))
  expect_true(is.ggplot(plot(sol,type='path')))
  expect_true(is.ggplot(plot(sol,type='front')))
})

test_that("GMM seed", {
  N=150
  X = rbind(MASS::mvrnorm(N/3,c(-5,0),diag(2)),MASS::mvrnorm(N/3,c(0,5),diag(2)),MASS::mvrnorm(N/3,c(5,0),diag(2)))
  mod = new("diaggmm",mu=apply(X,2,mean))
  sol=greed(X,model=mod,alg=new("seed"))
  expect_equal(sol@K, 3)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(sol,type='tree')))
  expect_true(is.ggplot(plot(sol,type='path')))
  expect_true(is.ggplot(plot(sol,type='front')))
})


test_that("GMM multistart", {
  N=150
  X = rbind(MASS::mvrnorm(N/3,c(-5,0),diag(2)),MASS::mvrnorm(N/3,c(0,5),diag(2)),MASS::mvrnorm(N/3,c(5,0),diag(2)))
  mod = new("diaggmm",mu=apply(X,2,mean))
  sol=greed(X,model=mod,alg=new("multistarts"))
  expect_gte(sol@K, 1)
  expect_lte(sol@K, 4)
  expect_true(is.ggplot(plot(sol,type='tree')))
  expect_true(is.ggplot(plot(sol,type='path')))
  expect_true(is.ggplot(plot(sol,type='front')))
})

context("GMM test")
library(greed)
library(ggplot2)
set.seed(1234)



test_that("GMM hybrid", {
  N=600
  X = rbind(MASS::mvrnorm(N/3,c(-5,0),diag(2)),MASS::mvrnorm(N/3,c(0,5),diag(2)),MASS::mvrnorm(N/3,c(5,0),diag(2)))
  sol=greed(X)
  expect_equal(sol@K, 3)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
})

test_that("GMM seed", {
  N=600
  X = rbind(MASS::mvrnorm(N/3,c(-5,0),diag(2)),MASS::mvrnorm(N/3,c(0,5),diag(2)),MASS::mvrnorm(N/3,c(5,0),diag(2)))
  sol=greed(X,alg=new("seed"))
  expect_equal(sol@K, 3)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
})


test_that("GMM multistart", {
  N=600
  X = rbind(MASS::mvrnorm(N/3,c(-5,0),diag(2)),MASS::mvrnorm(N/3,c(0,5),diag(2)),MASS::mvrnorm(N/3,c(5,0),diag(2)))
  sol=greed(X,alg=new("multistarts"))
  expect_equal(sol@K, 3)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
})

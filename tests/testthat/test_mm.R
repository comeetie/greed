context("MM test")
library(greed)
library(ggplot2)
set.seed(1234)

test_that("MM sim", {
  N = 500
  K = 10
  pi = rep(1/K,K)
  mu = cbind(diag(rep(5,K)),matrix(0,K,20))+matrix(runif(K*(20+K)),K,20+K)
  mm = rmm(N,pi,mu,15)
  expect_equal(dim(mm$x)[1], N)
  expect_equal(dim(mm$x)[2], 20+K)
  expect_equal(length(mm$cl),N)
  expect_gte(min(mm$cl), 1)
  expect_lte(max(mm$cl), K)
})


test_that("MM hybrid", {
  N = 200
  K = 4
  pi = rep(1/K,K)
  mu = cbind(diag(rep(5,K)),matrix(0,K,20))+matrix(runif(K*(20+K)),K,20+K)
  mm = rmm(N,pi,mu,30)
  sol=greed(mm$x,model=new("mm"))
  expect_equal(sol@K, K)
  solc = cut(sol,3)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
})

test_that("MM seed", {
  N = 200
  K = 4
  pi = rep(1/K,K)
  mu = cbind(diag(rep(5,K)),matrix(0,K,20))+matrix(runif(K*(20+K)),K,20+K)
  mm = rmm(N,pi,mu,15)
  sol=greed(mm$x,alg=new("seed"),model=new("mm"))
  expect_gte(sol@K, K-2)
  expect_lte(sol@K, K+2)
  solc = cut(sol,3)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
})


test_that("MM multistart", {
  N = 200
  K = 4
  pi = rep(1/K,K)
  mu = cbind(diag(rep(5,K)),matrix(0,K,20))+matrix(runif(K*(20+K)),K,20+K)
  mm = rmm(N,pi,mu,15)
  sol=greed(mm$x,alg=new("multistarts"),model=new("mm"))
  expect_gte(sol@K, K-2)
  expect_lte(sol@K, K+2)
  solc = cut(sol,3)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
})

context("DCSBM test")
library(greed)
library(ggplot2)
set.seed(1234)

test_that("DCSBM sim", {
  N = 100
  K = 5
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rdcsbm(N,pi,mu,rep(15,N),rep(15,N))
  expect_equal(dim(sbm$x)[1], N)
  expect_equal(dim(sbm$x)[2], N)
  expect_equal(length(sbm$cl),N)
  expect_gte(min(sbm$cl), 1)
  expect_lte(max(sbm$cl), K)
})


test_that("DCSBM hybrid", {
  N = 100
  K = 5
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rdcsbm(N,pi,mu,rep(15,N),rep(15,N))
  sol=greed(sbm$x)
  expect_equal(sol@K, K)
  solc = cut(sol,4)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})

test_that("DCSBM seed", {
  N = 100
  K = 5
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rdcsbm(N,pi,mu,rep(15,N),rep(15,N))
  sol=greed(sbm$x,alg=new("seed"))
  expect_gte(sol@K, K-2)
  expect_lte(sol@K, K+2)
  solc = cut(sol,4)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})


test_that("DCSBM multistart", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rdcsbm(N,pi,mu,rep(15,N),rep(15,N))
  sol=greed(sbm$x,alg=new("multistarts"))
  expect_gte(sol@K, K-2)
  expect_lte(sol@K, K+2)
  expect_true(is.ggplot(plot(sol,type='tree')))
  expect_true(is.ggplot(plot(sol,type='path')))
  expect_true(is.ggplot(plot(sol,type='front')))
  expect_true(is.ggplot(plot(sol,type='blocks')))
  expect_true(is.ggplot(plot(sol,type='nodelink')))
})

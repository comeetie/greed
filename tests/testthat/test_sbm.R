context("SBM test")
library(greed)
library(ggplot2)
set.seed(1234)

test_that("SBM sim", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))
  sbm = rsbm(N,pi,mu)
  expect_equal(dim(sbm$x)[1], N)
  expect_equal(dim(sbm$x)[2], N)
  expect_equal(length(sbm$cl),N)
  expect_gte(min(sbm$cl), 1)
  expect_lte(max(sbm$cl), K)
  
})


test_that("SBM hybrid", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rsbm(N,pi,mu)
  sol=greed(sbm$x,model=new('sbm'))
  expect_equal(sol@K, K)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})

test_that("SBM seed", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rsbm(N,pi,mu)
  sol=greed(sbm$x,model=new('sbm'),alg=new("seed"))
  expect_gte(sol@K, K-2)
  expect_lte(sol@K, K+2)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})


test_that("SBM multitstart", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rsbm(N,pi,mu)
  sol=greed(sbm$x,model=new('sbm'),alg=new("multistarts"))
  expect_gte(sol@K, K-2)
  expect_lte(sol@K, K+2)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})

test_that("SBM genetic", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rsbm(N,pi,mu)
  sol=greed(sbm$x,model=new('sbm'),alg=new("genetic"))
  expect_gte(sol@K, K-2)
  expect_lte(sol@K, K+2)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})

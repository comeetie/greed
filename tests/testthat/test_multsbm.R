context("MULTSBM test")
library(greed)
library(ggplot2)
set.seed(1234)

test_that("MULTSBM sim", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = array(dim=c(K,K,3))
  mu[,,1] = diag(rep(1/5,K))+runif(K^2)*0.005
  mu[1,1,1]=runif(1)*0.005
  mu[,,2] = diag(rep(1/5,K))+runif(K^2)*0.005
  mu[2,2,1]=runif(1)*0.005
  mu[,,3] = 1- mu[,,1]-mu[,,2]
  lambda = 10
  multsbm = rmultsbm(N,pi,mu,10)
  expect_equal(dim(multsbm$x)[1], N)
  expect_equal(dim(multsbm$x)[2], N)
  expect_equal(length(multsbm$cl),N)
  expect_gte(min(multsbm$cl), 1)
  expect_lte(max(multsbm$cl), K)
  
})


test_that("MULTSBM hybrid", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = array(dim=c(K,K,3))
  mu[,,1] = diag(rep(1/5,K))+runif(K^2)*0.005
  mu[1,1,1]=runif(1)*0.005
  mu[,,2] = diag(rep(1/5,K))+runif(K^2)*0.005
  mu[2,2,1]=runif(1)*0.005
  mu[,,3] = 1- mu[,,1]-mu[,,2]
  lambda = 10
  multsbm = rmultsbm(N,pi,mu,10)
  sol=greed(multsbm$x,model=new('multsbm'))
  expect_equal(sol@K, K)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})

test_that("MULTSBM seed", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = array(dim=c(K,K,3))
  mu[,,1] = diag(rep(1/5,K))+runif(K^2)*0.005
  mu[1,1,1]=runif(1)*0.005
  mu[,,2] = diag(rep(1/5,K))+runif(K^2)*0.005
  mu[2,2,1]=runif(1)*0.005
  mu[,,3] = 1- mu[,,1]-mu[,,2]
  lambda = 10
  multsbm = rmultsbm(N,pi,mu,10)
  sol=greed(multsbm$x,model=new('multsbm'),alg=new("seed"))
  expect_gte(sol@K, K-2)
  expect_lte(sol@K, K+2)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})


test_that("MULTSBM multitstart", {
  
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = array(dim=c(K,K,3))
  mu[,,1] = diag(rep(1/5,K))+runif(K^2)*0.005
  mu[1,1,1]=runif(1)*0.005
  mu[,,2] = diag(rep(1/5,K))+runif(K^2)*0.005
  mu[2,2,1]=runif(1)*0.005
  mu[,,3] = 1- mu[,,1]-mu[,,2]
  lambda = 10
  multsbm = rmultsbm(N,pi,mu,10)
  sol=greed(multsbm$x,model=new('multsbm'),alg=new("multistarts"))
  expect_gte(sol@K, K-2)
  expect_lte(sol@K, K+2)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})

test_that("MULTSBM genetic", {
  
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = array(dim=c(K,K,3))
  mu[,,1] = diag(rep(1/5,K))+runif(K^2)*0.005
  mu[1,1,1]=runif(1)*0.005
  mu[,,2] = diag(rep(1/5,K))+runif(K^2)*0.005
  mu[2,2,1]=runif(1)*0.005
  mu[,,3] = 1- mu[,,1]-mu[,,2]
  lambda = 10
  multsbm = rmultsbm(N,pi,mu,10)
  sol=greed(multsbm$x,model=new('multsbm'),alg=new("genetic"))
  expect_gte(sol@K, K-2)
  expect_lte(sol@K, K+2)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})

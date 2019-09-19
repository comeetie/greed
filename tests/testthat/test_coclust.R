context("COSBM test")
library(greed)
library(ggplot2)
set.seed(1234)



test_that("COSBM hybrid", {
  N = 500
  K = 10
  pi = rep(1/K,K)
  mu = cbind(diag(rep(5,K)),matrix(0,K,20))+matrix(runif(K*(20+K)),K,20+K)
  mm = rmm(N,pi,mu,40)
  sol=greed(mm$x,model=new('co_dcsbm'))
  expect_equal(sol@Krow, K)
  solc = cut(sol,21)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})

test_that("COSBM seed", {
  N = 500
  K = 10
  pi = rep(1/K,K)
  mu = lapply(1:5,function(x){cbind(diag(rep(10,K)),matrix(0,K,20))+matrix(runif(K*(20+K)),K,20+K)})
  mu = do.call(cbind,mu)
  mm = rmm(N,pi,mu,15)
  sol=greed(mm$x,model=new("co_dcsbm"),alg=new("seed"),K=40)
  expect_gte(sol@Krow, K-2)
  expect_lte(sol@Krow, K+2)
  solc = cut(sol,15)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})


test_that("COSBM multistart", {
  N = 500
  K = 10
  pi = rep(1/K,K)
  mu = lapply(1:5,function(x){cbind(diag(rep(10,K)),matrix(0,K,20))+matrix(runif(K*(20+K)),K,20+K)})
  mu = do.call(cbind,mu)
  mm = rmm(N,pi,mu,15)
  sol=greed(mm$x,model=new('co_dcsbm'),alg=new("multistarts"),K=40)
  expect_gte(sol@Krow, K-2)
  expect_lte(sol@Krow, K+2)
  solc = cut(sol,15)
  expect_true(is.ggplot(plot(sol,type='tree')))
  expect_true(is.ggplot(plot(sol,type='path')))
  expect_true(is.ggplot(plot(sol,type='front')))
  expect_true(is.ggplot(plot(sol,type='blocks')))
  expect_true(is.ggplot(plot(sol,type='nodelink')))
})

context("MISSSBM test")
library(greed)
library(ggplot2)
set.seed(1234)


test_that("MISSSBM hybrid", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rsbm(N,pi,mu)
  sbm$x[cbind(base::sample(1:N,100),base::sample(1:N,100))]=NA
  sol=greed(sbm$x,model=new('misssbm'))
  expect_equal(sol@K, K)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})

test_that("MISSSBM seed", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rsbm(N,pi,mu)
  sbm$x[cbind(base::sample(1:N,100),base::sample(1:N,100))]=NA
  sol=greed(sbm$x,model=new('misssbm'),alg=new("seed"))
  expect_gte(sol@K, K-2)
  expect_lte(sol@K, K+2)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})


test_that("MISSSBM multitstart", {
  N = 150
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rsbm(N,pi,mu)
  sbm$x[cbind(base::sample(1:N,100),base::sample(1:N,100))]=NA
  sol=greed(sbm$x,model=new('misssbm'),alg=new("multistarts"))
  expect_gte(sol@K, K-2)
  expect_lte(sol@K, K+2)
  expect_true(is.ggplot(plot(sol,type='tree')))
  expect_true(is.ggplot(plot(sol,type='path')))
  expect_true(is.ggplot(plot(sol,type='front')))
  expect_true(is.ggplot(plot(sol,type='blocks')))
  expect_true(is.ggplot(plot(sol,type='nodelink')))
})

test_that("MISSSBM genetic", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rsbm(N,pi,mu)
  sbm$x[cbind(base::sample(1:N,100),base::sample(1:N,100))]=NA
  sol=greed(sbm$x,model=new('misssbm'),alg=new("genetic"))
  expect_gte(sol@K, K-2)
  expect_lte(sol@K, K+2)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})

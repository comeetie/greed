context("SBM test")
library(greed)
library(ggplot2)
library(Matrix)
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
  co=coef(sol)
  expect_true(all(dim(co$thetakl)==c(3,3)))
  expect_equal(sum(co$pi),1)
  expect_equal(length(co$pi),3)
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
  co=coef(sol)
  expect_true(all(dim(co$thetakl)==c(3,3)))
  expect_equal(sum(co$pi),1)
  expect_equal(length(co$pi),3)
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
  co=coef(sol)
  expect_true(all(dim(co$thetakl)==c(3,3)))
  expect_equal(sum(co$pi),1)
  expect_equal(length(co$pi),3)
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
  co=coef(sol)
  expect_true(all(dim(co$thetakl)==c(3,3)))
  expect_equal(sum(co$pi),1)
  expect_equal(length(co$pi),3)
  expect_gte(sol@K, K-2)
  expect_lte(sol@K, K+2)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})


test_that("SBM hybrid unidrected", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rsbm(N,pi,mu)
  x=tril(sbm$x)+t(tril(sbm$x))
  diag(x)=0
  sol=greed(x,model=new('sbm',type="undirected"))
  co=coef(sol)
  expect_true(all(dim(co$thetakl)==c(3,3)))
  expect_equal(sum(co$pi),1)
  expect_equal(length(co$pi),3)
  expect_equal(sol@K, K)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})

test_that("SBM seed undirected", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rsbm(N,pi,mu)
  x=tril(sbm$x)+t(tril(sbm$x))
  diag(x)=0
  sol=greed(x,model=new('sbm',type="undirected"),alg=new("seed"))
  co=coef(sol)
  expect_true(all(dim(co$thetakl)==c(3,3)))
  expect_equal(sum(co$pi),1)
  expect_equal(length(co$pi),3)
  expect_gte(sol@K, K-2)
  expect_lte(sol@K, K+2)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})

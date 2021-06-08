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

test_that("MISSSBM seed", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rsbm(N,pi,mu)
  sbm$x[cbind(base::sample(1:N,100),base::sample(1:N,100))]=NA
  sol=greed(sbm$x,model=new('misssbm'),alg=new("seed"))
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


test_that("MISSSBM hybrid mar", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rsbm(N,pi,mu)
  sbm$x[cbind(base::sample(1:N,100),base::sample(1:N,100))]=NA
  sol=greed(sbm$x,model=new('misssbm',sampling="dyad"))
  co=coef(sol)
  expect_true(all(dim(co$thetakl)==c(3,3)))
  expect_equal(sum(co$pi),1)
  expect_equal(length(co$pi),3)
  expect_equal(sol@K, K)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(sol,type='tree')))
  expect_true(is.ggplot(plot(sol,type='path')))
  expect_true(is.ggplot(plot(sol,type='front')))
  expect_true(is.ggplot(plot(sol,type='blocks')))
  expect_true(is.ggplot(plot(sol,type='nodelink')))
})

test_that("MISSSBM seed mar", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rsbm(N,pi,mu)
  sbm$x[cbind(base::sample(1:N,100),base::sample(1:N,100))]=NA
  sol=greed(sbm$x,model=new('misssbm',sampling="dyad"),alg=new("seed"))
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





test_that("MISSSBM hybrid undirected", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rsbm(N,pi,mu)
  sbm$x[cbind(base::sample(1:N,100),base::sample(1:N,100))]=NA
  x=tril(sbm$x)+t(tril(sbm$x))
  diag(x)=0
  sol=greed(x,model=new('misssbm',type="undirected"))
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

test_that("MISSSBM seed undirected", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rsbm(N,pi,mu)
  sbm$x[cbind(base::sample(1:N,100),base::sample(1:N,100))]=NA
  x=tril(sbm$x)+t(tril(sbm$x))
  diag(x)=0
  sol=greed(x,model=new('misssbm',type="undirected"),alg=new("seed"))
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


test_that("MISSSBM hybrid mar undirected", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rsbm(N,pi,mu)
  sbm$x[cbind(base::sample(1:N,100),base::sample(1:N,100))]=NA
  x=tril(sbm$x)+t(tril(sbm$x))
  diag(x)=0
  sol=greed(x,model=new('misssbm',type="undirected",sampling="dyad"))
  co=coef(sol)
  expect_true(all(dim(co$thetakl)==c(3,3)))
  expect_equal(sum(co$pi),1)
  expect_equal(length(co$pi),3)
  expect_equal(sol@K, K)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(sol,type='tree')))
  expect_true(is.ggplot(plot(sol,type='path')))
  expect_true(is.ggplot(plot(sol,type='front')))
  expect_true(is.ggplot(plot(sol,type='blocks')))
  expect_true(is.ggplot(plot(sol,type='nodelink')))
})

test_that("MISSSBM seed block node", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rsbm(N,pi,mu)
  no_obs_node = sample(1:100,20)
  sbm$x[no_obs_node,no_obs_node]=NA
  sol=greed(sbm$x,model=new('misssbm',sampling="block-node"),alg=new("seed"))
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
  expect_true(is.ggplot(plot(sol,type='blocks')))
  expect_true(is.ggplot(plot(sol,type='nodelink')))
})

test_that("MISSSBM seed block node undirected", {
  N = 100
  K = 3
  pi = rep(1/K,K)
  mu = diag(rep(1/5,K))+runif(K*K)*0.01
  sbm = rsbm(N,pi,mu)
  no_obs_node = sample(1:100,20)
  x=tril(sbm$x)+t(tril(sbm$x))
  diag(x)=0
  x[no_obs_node,no_obs_node]=NA
  sol=greed(x,model=new('misssbm',sampling="block-node",type="undirected"),alg=new("seed"))
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
  expect_true(is.ggplot(plot(sol,type='blocks')))
  expect_true(is.ggplot(plot(sol,type='nodelink')))
})






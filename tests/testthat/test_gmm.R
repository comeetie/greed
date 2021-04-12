context("GMM test")
library(greed)
library(ggplot2)
set.seed(1234)



test_that("GMM marginal",{
  N=200
  X1 = MASS::mvrnorm(N/2,c(-5,0),diag(2))
  X2 = MASS::mvrnorm(N/2,c(0,5),diag(2))
  R1 = greed:::gmm_marginal(X1,0.1,3,diag(2),c(0,0))
  R2 = greed:::gmm_marginal(X2,0.1,3,diag(2),c(0,0))
  Rm = greed:::gmm_marginal_merge(R1,R2,0.1,3,diag(2),c(0,0))
  R = greed:::gmm_marginal(rbind(X1,X2),0.1,3,diag(2),c(0,0))
  
  expect_lte(abs(Rm$log_evidence-R$log_evidence),10^-6)
  
  R1m = greed:::gmm_marginal(X1[1:(N/2-1),],0.1,3,diag(2),c(0,0))
  R1mo = greed:::gmm_marginal_del1(R1,X1[N/2,],0.1,3,diag(2),c(0,0))
  expect_lte(sum(abs(R1m$S-R1mo$S)),10^-6)
  expect_lte(sum(abs(R1m$m-R1mo$m)),10^-6)
  expect_lte(abs(R1m$log_evidence-R1mo$log_evidence),10^-6)
  R1a = greed:::gmm_marginal_add1(R1m,X1[N/2,],0.1,3,diag(2),c(0,0))
  expect_lte(sum(abs(R1a$S-R1$S)),10^-6)
  expect_lte(sum(abs(R1a$m-R1$m)),10^-6)
  expect_lte(abs(R1a$log_evidence-R1$log_evidence),10^-6)
})

test_that("GMM hybrid", {
  N=150
  X = rbind(MASS::mvrnorm(N/3,c(-5,0),diag(2)),MASS::mvrnorm(N/3,c(0,5),diag(2)),MASS::mvrnorm(N/3,c(5,0),diag(2)))
  sol=greed(X)
  expect_equal(sol@K, 3)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
})

test_that("GMM seed", {
  N=150
  X = rbind(MASS::mvrnorm(N/3,c(-5,0),diag(2)),MASS::mvrnorm(N/3,c(0,5),diag(2)),MASS::mvrnorm(N/3,c(5,0),diag(2)))
  sol=greed(X,alg=new("seed"))
  expect_equal(sol@K, 3)
  solc = cut(sol,2)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
})


test_that("GMM multistart", {
  N=150
  X = rbind(MASS::mvrnorm(N/3,c(-5,0),diag(2)),MASS::mvrnorm(N/3,c(0,5),diag(2)),MASS::mvrnorm(N/3,c(5,0),diag(2)))
  sol=greed(X,alg=new("multistarts"))
  expect_gte(sol@K, 1)
  expect_lte(sol@K, 4)
  expect_true(is.ggplot(plot(sol,type='tree')))
  expect_true(is.ggplot(plot(sol,type='path')))
  expect_true(is.ggplot(plot(sol,type='front')))
})

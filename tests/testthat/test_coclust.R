context("COSBM test")
library(greed)
library(ggplot2)
set.seed(1234)



test_that("COSBM hybrid", {
  mu=cbind(lower.tri(matrix(1,4,4)),upper.tri(matrix(1,4,4)))*0.5+0.01
  mu[2,4]=0.51
  mu[3,5]=0.51
  mm=rlbm(150,300,rep(1/4,4),rep(1/8,8),mu)
  sol=greed(mm$x,model=new('co_dcsbm'))
  expect_equal(sol@Krow, 4)
  expect_equal(sol@Kcol, 8)
  solc = cut(sol,8)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})

test_that("COSBM seed", {
  mu=cbind(lower.tri(matrix(1,4,4)),upper.tri(matrix(1,4,4)))*0.5+0.01
  mu[2,4]=0.51
  mu[3,5]=0.51
  mm=rlbm(150,300,rep(1/4,4),rep(1/8,8),mu)
  sol=greed(mm$x,model=new("co_dcsbm"),alg=new("seed"),K=40)
  expect_gte(sol@K, 12-2)
  expect_lte(sol@K, 12+2)
  solc = cut(sol,8)
  expect_true(is.ggplot(plot(solc,type='tree')))
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
})


test_that("COSBM multistart", {
  mu=cbind(lower.tri(matrix(1,4,4)),upper.tri(matrix(1,4,4)))*0.5+0.01
  mu[2,4]=0.51
  mu[3,5]=0.51
  mm=rlbm(150,300,rep(1/4,4),rep(1/8,8),mu)
  sol=greed(mm$x,model=new('co_dcsbm'),alg=new("multistarts"),K=40)
  expect_gte(sol@K, 12-2)
  expect_lte(sol@K, 12+2)
  solc = cut(sol,8)
  expect_true(is.ggplot(plot(sol,type='tree')))
  expect_true(is.ggplot(plot(sol,type='path')))
  expect_true(is.ggplot(plot(sol,type='front')))
  expect_true(is.ggplot(plot(sol,type='blocks')))
  expect_true(is.ggplot(plot(sol,type='nodelink')))
})

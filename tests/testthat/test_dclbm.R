context("COSBM test")
library(greed)
library(ggplot2)
set.seed(1234)



test_that("DcLBM hybrid", {
  mu=cbind(lower.tri(matrix(1,4,4)),upper.tri(matrix(1,4,4)))*0.2+0.01
  mu[2,4]=0.21
  mu[3,5]=0.21
  mm=rlbm(1000,800,rep(1/4,4),rep(1/8,8),mu)
  model = DcLbm()
  data = greed:::preprocess(model,mm$x)
  i=sample(200,1)
  oldcl=mm$clr[i]
  newcl = sample(setdiff(1:mm$Kr,oldcl),1)
  cl = c(mm$clr,mm$clc+max(mm$clr))
  expect_lte(greed:::test_swap(model,data,cl,i,newcl),10^-6)
  expect_lte(greed:::test_merge(model,data,cl,oldcl,newcl),10^-6)

})

test_that("COSBM seed", {
  mu=cbind(lower.tri(matrix(1,4,4)),upper.tri(matrix(1,4,4)))*0.2+0.01
  mu[2,4]=0.21
  mu[3,5]=0.21
  mm=rlbm(200,400,rep(1/4,4),rep(1/8,8),mu)
  model = DcLbm()
  sol=greed(mm$x)
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
  mu=cbind(lower.tri(matrix(1,4,4)),upper.tri(matrix(1,4,4)))*0.2+0.01
  mu[2,4]=0.21
  mu[3,5]=0.21
  mm=rlbm(200,400,rep(1/4,4),rep(1/8,8),mu)
  sol=greed(mm$x,model=DcLbm(),alg=Multistarts(),K=40)
  expect_gte(sol@K, 12-2)
  expect_lte(sol@K, 12+2)
  solc = cut(sol,8)
  expect_true(is.ggplot(plot(sol,type='tree')))
  expect_true(is.ggplot(plot(sol,type='path')))
  expect_true(is.ggplot(plot(sol,type='front')))
  expect_true(is.ggplot(plot(sol,type='blocks')))
  expect_true(is.ggplot(plot(sol,type='nodelink')))
})

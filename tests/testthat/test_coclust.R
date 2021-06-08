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
  expect_true(is.ggplot(plot(solc,type='path')))
  expect_true(is.ggplot(plot(solc,type='front')))
  expect_true(is.ggplot(plot(solc,type='blocks')))
  expect_true(is.ggplot(plot(solc,type='nodelink')))
  co = coef(sol)
  expect_equal(length(co$picols), 8)
  expect_equal(sum(co$picols), 1)
  
  expect_equal(length(co$pirows), 4)
  expect_equal(sum(co$pirows), 1)
  
  expect_equal(nrow(co$thetakl), 4)
  expect_equal(ncol(co$thetakl), 8)
  
  expect_equal(sum(co$gammarows),4)
  expect_equal(sum(co$gammacols),8)
})



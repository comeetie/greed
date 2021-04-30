context("GMM test")
library(greed)
library(ggplot2)
set.seed(1234)




test_that("MVMREG hybrid", {
  regs=rmreg(200,rep(1/3,3),mu=cbind(c(5,1,-125),c(1,20,1),c(15,100,200)),sigma=1)
  sol=greed_cond(regs$X,as.matrix(regs$y))
  expect_equal(sol@K, 3)
  expect_true(is.ggplot(plot(sol,type='tree')))
  expect_true(is.ggplot(plot(sol,type='path')))
  expect_true(is.ggplot(plot(sol,type='front')))
})

test_that("GMM seed", {
  regs=rmreg(500,rep(1/3,3),mu=cbind(c(5,1,-125),c(1,20,1),c(15,100,200)),sigma=1)
  sol=greed_cond(regs$X,as.matrix(regs$y),alg=new("seed"))
  expect_equal(sol@K, 3)
})


test_that("GMM multistart", {
  regs=rmreg(500,rep(1/3,3),mu=cbind(c(5,1,-125),c(1,20,1),c(15,100,200)),sigma=1)
  sol=greed_cond(regs$X,as.matrix(regs$y),alg=new("multistarts"))
  expect_gte(sol@K, 1)
  expect_lte(sol@K, 4)
  expect_true(is.ggplot(plot(sol,type='tree')))
  expect_true(is.ggplot(plot(sol,type='path')))
  expect_true(is.ggplot(plot(sol,type='front')))
})


context("MOR test")
library(greed)
library(ggplot2)
set.seed(1234)




test_that("MOR hybrid", {
  regs=rmreg(2000,rep(1/3,3),mu=cbind(c(5,1,-125),c(1,20,1),c(15,100,200)),sigma=1)
  df = data.frame(x1=regs$X[,1],x2=regs$X[,2],y=regs$y)
  sol=greed(df,MoR(y ~ x1 + x2))
  expect_equal(sol@K, 3)
  expect_true(is.ggplot(plot(sol,type='tree')))
  expect_true(is.ggplot(plot(sol,type='path')))
  expect_true(is.ggplot(plot(sol,type='front')))
})

test_that("MVMREG seed", {
  regs=rmreg(500,rep(1/3,3),mu=cbind(c(5,1,-125),c(1,20,1),c(15,100,200)),sigma=1)
  sol=greed_cond(regs$X,as.matrix(regs$y),alg=new("seed"))
  expect_equal(sol@K, 3)
})

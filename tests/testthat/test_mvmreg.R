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
  co=coef(sol)
  expect_true(all(sapply(co$muk,length)==3))
  expect_true(all(sapply(co$Sigmak,nrow)==1))
  expect_equal(sum(co$pi),1)
  expect_equal(length(co$pi),3)
})

test_that("GMM seed", {
  regs=rmreg(500,rep(1/3,3),mu=cbind(c(5,1,-125),c(1,20,1),c(15,100,200)),sigma=1)
  sol=greed_cond(regs$X,as.matrix(regs$y),alg=new("seed"))
  expect_equal(sol@K, 3)
})



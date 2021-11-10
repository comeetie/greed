context("Mixed mixture model test (factor and numeric only)")
library(greed)
library(ggplot2)
set.seed(1234)




test_that("[Mixed mixture model preprocess] Character input",{
  n = 10
  numerical = data.frame(num1 = rnorm(n), num2 = rnorm(n, mean=5))
  cats = data.frame(f1 = sample(c('a', 'b'), size = n, replace = T), 
                    f2 = sample(c('10', '20', 'NoNe'), size = n, replace = T),
                    f3 = factor(sample(c(1,4), size = n, replace = T)))
  X = cbind(numerical, cats)
  
  model = new("mmm")
  data=greed:::preprocess(model, X)
  expect_equal(ncol(numerical), ncol(data$Xnum))
  expect_equal(ncol(cats), ncol(data$Xcat))
  expect_equal(apply(data$Xcat, 2, function(col) dplyr::n_distinct(col)),
               sapply(cats, function(col) dplyr::n_distinct(col)))
})


test_that("[Mixed mixture model preprocess] Drop missing factors",{
  # add unused levels that should be dropped by preprocess(new('mmm), ...)
  n = 10
  numerical = data.frame(num1 = rnorm(n), num2 = rnorm(n, mean=5))
  cats = data.frame(f1 = sample(c('a', 'b'), size = n, replace = T), 
                    f2 = sample(c('10', '20', 'NoNe'), size = n, replace = T),
                    f3 = factor(sample(c(1,4), size = n, replace = T)), 
                    stringsAsFactors = T)
  X = cbind(numerical, cats)
  Xaug = X[!(as.numeric(X$f2) == 2),] # remove second factor of $f2 column
  
  model = new("mmm")
  data=greed:::preprocess(model, Xaug)
  expect_equal(apply(data$Xcat, 2, function(col) dplyr::n_distinct(col)), apply(data$Xcat, 2, function(col) max(col) + 1))
})

test_that("mmm optim",{
  N=500
  theta = list(matrix(c(0.1,0.9,0.9,0.1,0.8,0.2,0.05,0.95),ncol=2,byrow=TRUE),
               matrix(c(0.95,0.05,0.3,0.7,0.05,0.95,0.05,0.95),ncol=2,byrow=TRUE),
               matrix(c(0.95,0.04,0.01,0.9,0.09,0.01,0.01,0.01,0.98,0.9,0.05,0.05),ncol=3,byrow=TRUE),
               matrix(c(1,0,0,1,1,0,0,1),ncol=2,byrow=TRUE))
  lca.data = rlca(N,rep(1/4,4),theta)
  mu = c(-5,-2,2,5)
  X = lca.data$x
  X$num=sapply(lca.data$cl,function(cli){rnorm(1,mu[cli],0.7)})
  sol=greed(X,model=new("mmm"))
  expect_equal(sol@K,4)
  expect_equal(NMI(sol@cl,lca.data$cl),1)
  
})

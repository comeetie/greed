context("Mixed mixture model test (factor and numeric only)")
library(greed)
library(ggplot2)
set.seed(1234)




testthat::test_that("[Mixed mixture model preprocess] Character input",{
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


testthat::test_that("[Mixed mixture model preprocess] Drop missing factors",{
  # add unused levels that should be dropped by preprocess(new('mmm), ...)
  n = 10
  numerical = data.frame(num1 = rnorm(n), num2 = rnorm(n, mean=5))
  cats = data.frame(f1 = sample(c('a', 'b'), size = n, replace = T), 
                    f2 = sample(c('10', '20', 'NoNe'), size = n, replace = T),
                    f3 = factor(sample(c(1,4), size = n, replace = T)), 
                    stringsAsFactors = T)
  X = cbind(numerical, cats)
  Xaug = X[!(as.numeric(Xaug$f2) == 2),] # remove second factor of $f2 column
  
  model = new("mmm")
  data=greed:::preprocess(model, Xaug)
  expect_equal(apply(data$Xcat, 2, function(col) dplyr::n_distinct(col)), apply(data$Xcat, 2, function(col) max(col) + 1))
})
context("LCA test")
library(greed)
library(ggplot2)
set.seed(1234)



test_that("LCA sim",{
  N=1000
  theta = list(matrix(c(0.1,0.9,0.9,0.1,0.8,0.2,0.05,0.95),ncol=2,byrow=TRUE),
                matrix(c(0.95,0.05,0.3,0.7,0.05,0.95,0.05,0.95),ncol=2,byrow=TRUE),
                matrix(c(0.95,0.04,0.01,0.9,0.09,0.01,0.01,0.01,0.98,0.9,0.05,0.05),ncol=3,byrow=TRUE),
               matrix(c(0.7,0,0,0.7,0.7,0,0,0.7),ncol=2,byrow=TRUE))
  lca.data = rlca(N,rep(1/4,4),theta)
  
  
  sol=greed(lca.data$x,model=new("lca"),K=10,verbose = TRUE)
  sol@icl  
  model=new("lca")
  data=greed:::preprocess(model,lca.data$x)

  temp=greed:::fit_greed(model,data,lca.data$cl,"none")
  temp@icl
  sol@icl
  
  
  sol@obs_stats$counts
  sol@obs_stats$x_counts
  table(sol@cl,lca.data$cl)
  model=new("lca",beta=5)
  data=greed:::preprocess(model,lca.data$x)
  sol=greed:::fit_greed(model,data,lca.data$cl,type = "none")
  sol@obs_stats[[2]]$x_counts
  sol@icl
  sol=greed:::fit_greed(model,data,lca.data$cl,type = "both")
  sol@obs_stats[[2]]$x_counts[[3]]
  sol@icl
  table(sol@cl,lca.data$cl)  
  
  theta = list(matrix(c(0.01,0.99,0.99,0.01),ncol=2,byrow=TRUE),
               matrix(c(0.01,0.99,0.99,0.01),ncol=2,byrow=TRUE),
               matrix(c(0.95,0.05,0.05,0.95),ncol=2,byrow=TRUE))
  

  
  sol=greed:::fit_greed(model,data,sample(1:10,N, replace = TRUE),type = "both")
  table(sol@cl,lca.data$cl)
  
})



test_that("LCA icl opt",{
N=500
theta = list(matrix(c(0.1,0.9,0.9,0.1,0.8,0.2,0.05,0.95),ncol=2,byrow=TRUE),
               matrix(c(0.95,0.05,0.3,0.7,0.05,0.95,0.05,0.95),ncol=2,byrow=TRUE),
               matrix(c(0.95,0.04,0.01,0.9,0.09,0.01,0.01,0.01,0.98,0.9,0.05,0.05),ncol=3,byrow=TRUE),
               matrix(c(1,0,0,1,1,0,0,1),ncol=2,byrow=TRUE))
lca.data = rlca(N,rep(1/4,4),theta)
model= new("lca")
data=greed:::preprocess(model,lca.data$x)
cl =lca.data$cl
i = sample(500,1)
newcl = sample(setdiff(1:lca.data$K,cl[i]),1)
expect_lte(greed:::test_swap(model,data,cl,i,newcl), 10^-6)
expect_lte(greed:::test_merge(model,data,cl,1,4), 10^-6)
})


testthat::test_that("[LCA preprocess] Character input",{
  X = data.frame(x=c('a','b', 'a'), 
                 y= c('ctrl', 'tab', 'ctrl'), 
                 z = c('left', 'right', 'None'))
  model = new("lca")
  data=greed:::preprocess(model, X)
  expect_equal(apply(data$X, 2, function(col) dplyr::n_distinct(col)),sapply(X, function(col) dplyr::n_distinct(col)))
})


testthat::test_that("[LCA preprocess] Drop missing factors",{
  N=4
  theta = list(matrix(c(0.1,0.9,0.9,0.1,0.8,0.2,0.05,0.95),ncol=2,byrow=TRUE),
               matrix(c(0.95,0.05,0.3,0.7,0.05,0.95,0.05,0.95),ncol=2,byrow=TRUE),
               matrix(c(0.95,0.04,0.01,0.9,0.09,0.01,0.01,0.01,0.98,0.9,0.05,0.05),ncol=3,byrow=TRUE),
               matrix(c(0.7,0,0,0.7,0.7,0,0,0.7),ncol=2,byrow=TRUE))
  lca.data = rlca(N,rep(1/4,4),theta)
  model = new("lca")
  data=greed:::preprocess(model, lca.data$x)
  expect_equal(apply(data$X, 2, function(col) dplyr::n_distinct(col)), apply(data$X, 2, function(col) max(col) + 1))
})

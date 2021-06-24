context("LCA test")
library(greed)
library(ggplot2)
set.seed(1234)



test_that("LCA sim",{
  N=500
  theta = list(matrix(c(0.1,0.9,0.9,0.1,0.8,0.2,0.05,0.95),ncol=2,byrow=TRUE),
                matrix(c(0.95,0.05,0.3,0.7,0.05,0.95,0.05,0.95),ncol=2,byrow=TRUE),
                matrix(c(0.95,0.04,0.01,0.9,0.09,0.01,0.01,0.01,0.98,0.9,0.05,0.05),ncol=3,byrow=TRUE),
               matrix(c(1,0,0,1,1,0,0,1),ncol=2,byrow=TRUE))
  lca.data = rlca(N,rep(1/4,4),theta)
  
  
  sol=greed(lca.data$x,model=new("lca"),alg=new("hybrid",pop_size=100),K=10)
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

context("MREG test")
library(greed)
library(ggplot2)

test_that("MReg sim evidence and fit", {
  mreg_simu=rmreg(1000,c(0.5,0.5),cbind(c(5,5),c(-5,-5)),0.8)
  
  reg=0.1
  sK21 = greed:::lm_post(mreg_simu$X[mreg_simu$cl==1,],mreg_simu$y[mreg_simu$cl==1],reg,1,1)
  sK22 = greed:::lm_post(mreg_simu$X[mreg_simu$cl==2,],mreg_simu$y[mreg_simu$cl==2],reg,1,1)
  sKm  = greed:::lm_post_merge(sK21,sK22,reg,1,1)
  sK1  = greed:::lm_post(mreg_simu$X,mreg_simu$y,reg,1,1)
  expect_equal(abs(sK1$log_evidence-sKm$log_evidence),0)
  
  expect_lte(sK1$log_evidence,  sK21$log_evidence+sK22$log_evidence)

  
  skadd = greed:::lm_post_add1(sK21,mreg_simu$X[1,],mreg_simu$y[1],reg,1,1)
  skdel = greed:::lm_post_del1(skadd,mreg_simu$X[1,],mreg_simu$y[1],reg,1,1)
  
  expect_equal(sK21$log_evidence,  skdel$log_evidence)
  
  
  sol=greed:::fit_greed(new("mreg"),list(X=mreg_simu$X,y=mreg_simu$y),mreg_simu$cl,type = "both")
  tcf=table(sol@cl,mreg_simu$cl)
  
  
  nb_malcl = sum(apply(tcf,1,min))
  expect_lte(nb_malcl,  100)
  
  
  sol=greed:::fit_greed(new("mreg",alpha=1,a0=1,b0=1,reg=0.1),list(X=mreg_simu$X,y=mreg_simu$y),sample(1:5,length(mreg_simu$y),replace = TRUE))
  tcf= table(sol@cl,mreg_simu$cl)
  nb_malcl = sum(apply(tcf,1,min))
  expect_lte(nb_malcl,  100)
  
})




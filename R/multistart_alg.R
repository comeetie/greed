#' @useDynLib greed
#' @importFrom Rcpp sourceCpp
#' @importFrom future %<-%
#' @name %<-%
NULL

#' @include models_classes.R fit_classes.R cleanpath.R
#' @import Matrix
NULL

multistart = function(model,alg,data,K,verbose=FALSE){
  
  solutions = listenv::listenv()
  for (i in 1:alg@nb_start){
    solutions[[i]] %<-% fit_greed(model,data,sample_cl(model,data,K))
  }
  solutions = as.list(solutions)
  icls = sapply(solutions,function(s){s@icl})
  
  res = solutions[[order(icls,decreasing = TRUE)[1]]]
  path = fit_greed_path(data,res)
  path = cleanpathopt(path)
  path@train_hist = data.frame(icl=icls,K= sapply(solutions,function(s){max(s@cl)}))
  
  path
}


multi_swap = function(model,alg,data,K,verbose=FALSE){
  
  solutions = listenv::listenv()
  for (i in 1:alg@nb_start){
    solutions[[i]] %<-% fit_greed(model,data,sample_cl(model,data,K))
  }
  solutions = as.list(solutions)
  solutions
}
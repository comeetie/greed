#' @useDynLib greed
#' @importFrom Rcpp sourceCpp
#' @importFrom future %<-%
#' @name %<-%
NULL

#' @include models_classes.R fit_classes.R cleanpath.R
#' @import Matrix
NULL

multistart = function(greed_f,path_f,alg,verbose=FALSE){
  
  solutions = listenv::listenv()
  for (i in 1:alg@nb_start){
    solutions[[i]] %<-% greed_f()
  }
  solutions = as.list(solutions)
  icls = sapply(solutions,function(s){s@icl})
  
  res = solutions[[order(icls,decreasing = TRUE)[1]]]
  path = path_f(res)
  cleanpath(path)
  path@train_hist = data.frame(icl=icls,K= sapply(solutions,function(s){max(s@cl)}))
  
  path
}
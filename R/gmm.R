#' @include models_classes.R fit_classes.R mvmreg.R
NULL


#' @title Clustering with a gaussian mixture model description class
#' 
#' @description 
#' An S4 class to represent a multivariate mixture of regression model, extend \code{\link{icl_model-class}}.
#' @slot name name of the model
#' @slot alpha Dirichlet over cluster proportions prior parameter (default to 1)
#' @slot tau Prior parameter (inverse variance) default 0.01 
#' @slot N0 Prior parameter (pseudo count) defulat to 10 ! should be > number of features
#' @slot epsilon Prior parameter covartiance matrix prior
#' @slot mu mean prior
#' @examples
#' new("gmm")
#' new("gmm",alpha=1,tau=0.1,N0=15)
#' @export
setClass("gmm", representation = list(tau = "numeric",mu="numeric",epsilon="matrix",N0="numeric"),
         contains = "icl_model",
         prototype(name="gmm",tau=0.1,N0=10,mu=1,epsilon=matrix(1,1,1),alpha=1))


#' @title Clustering with a gaussian mixture model fit results class
#' 
#' @description An S4 class to represent a fit of a multivariate mixture of regression model, extend \code{\link{icl_fit-class}}.
#' @slot model a \code{\link{mvmreg-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and clolumns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item regs: list of size $K$ with statistics for each clusters
#' }
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history infromation (details depends on the training procedure)
#' @export 
setClass("gmm_fit",slots = list(model="gmm"),contains="icl_fit")


#' @title Clustering with a multivariate mixture of regression model path extraction results class
#' 
#' 
#' @description An S4 class to represent a hierarchical fit of a gaussian mixture model, extend \code{\link{icl_path-class}}.
#' @slot model a \code{\link{gmm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and clolumns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item regs: list of size $K$ with statistics for each clusters
#' }
#' @slot path a list of size K-1 with each part of the path described by:
#' \itemize{
#' \item icl1: icl value reach with this solution for alpha=1 
#' \item logalpha: log(alpha) value were this solution is better than its parent
#' \item K: number of clusters
#' \item cl: vector of cluster indexes
#' \item k,l: index of the cluster that were merged at this step
#' \item merge_mat: lower triangular matrix of delta icl values 
#' \item obs_stats: a list with the same elements
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy ploting with gggplot
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father  
#' @slot train_hist  data.frame with training history infromation (details depends on the training procedure)
#' @export 
setClass("gmm_path",contains=c("icl_path","gmm_fit"))



#' @title plot a \code{\link{gmm_path-class}} object
#' 
#' 
#' @param x a \code{\link{gmm_path-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'front'}: plot the extracted front ICL, log(alpha)
#' \item \code{'path'}: plot the veolution of ICL with repsect to K
#' \item \code{'tree'}: plot the associated dendogram
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export 
setMethod(f = "plot", 
          signature = signature("gmm_path","missing"),
          definition = function(x,type='tree'){
            switch(type,tree = {
              dendo(x)
            },
            path ={
              lapath(x)
            },
            front = {
              plot_front(x)
            })
          })



setMethod(f = "seed", 
          signature = signature("gmm","list","numeric"), 
          definition = function(model,data, K){
            km=stats::kmeans(zscore(data$X),K)
            km$cluster
          })



setMethod(f = "preprocess", 
          signature = signature("gmm"), 
          definition = function(model, data,K){
            list(X=as.matrix(data),N=nrow(data),moves=as.sparse(matrix(1,K,K)))
          })

reorder_gmm = function(obs_stats,or){
  obs_stats$counts = obs_stats$counts[or]
  obs_stats$regs = obs_stats$regs[or]
  obs_stats
}


setMethod(f = "reorder", 
          signature = signature("gmm", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_gmm(obs_stats,order)
          })

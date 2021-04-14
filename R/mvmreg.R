#' @include models_classes.R fit_classes.R
NULL


#' @title Multivariate mixture of regression model description class
#' 
#' @description 
#' An S4 class to represent a multivariate mixture of regression model, extends \code{\link{icl_model-class}}.
#' #' The model follow [minka-linear](https://tminka.github.io/papers/minka-linear.pdf).
#' @slot name name of the model
#' @slot alpha Dirichlet over cluster proportions prior parameter (default to 1)
#' @slot beta Prior parameter (inverse variance) default 0.01 
#' @slot N0 Prior parameter (pseudo count) default to 10 ! should be > number of features
#' @examples
#' new("mvmreg")
#' new("mvmreg",alpha=1,beta=0.1,N0=15)
#' @md
#' @export
setClass("mvmreg",
         representation = list(beta = "numeric",N0="numeric"),
         contains = "icl_model",
         prototype(name="mvmreg",beta=0.01,N0=10,alpha=1))


#' @title Clustering with a multivariate mixture of regression model fit results class
#' 
#' @description An S4 class to represent a fit of a multivariate mixture of regression model, extend \code{\link{icl_fit-class}}.
#' @slot model a \code{\link{mvmreg-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item regs: list of size $K$ with statistics for each clusters
#' }
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details depends on the training procedure)
#' @export 
setClass("mvmreg_fit",slots = list(model="mvmreg"),contains="icl_fit")


#' @title Multivariate mixture of regression model hierarchical fit results class
#' 
#' 
#' @description An S4 class to represent a hierarchical fit of a multivariate mixture of regression model, extend \code{\link{icl_path-class}}.
#' @slot model a \code{\link{mvmreg-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
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
#' \item obs_stats: a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item regs: list of size $K$ with statistics for each clusters
#' }
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy plotting with \code{ggplot2}
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father  
#' @slot train_hist  data.frame with training history information (details depends on the training procedure)
#' @export 
setClass("mvmreg_path",contains=c("icl_path","mvmreg_fit"))



#' @title plot a \code{\link{mvmreg_path-class}} object
#' 
#' 
#' @param x a \code{\link{mvmreg_path-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'front'}: plot the extracted front ICL, log(alpha)
#' \item \code{'path'}: plot the evolution of ICL with respect to K
#' \item \code{'tree'}: plot the associated dendrogram
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export 
setMethod(f = "plot", 
          signature = signature("mvmreg_path","missing"),
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

reorder_mvmreg = function(obs_stats,or){
  obs_stats$counts = obs_stats$counts[or]
  obs_stats$regs = obs_stats$regs[or]
  obs_stats
}


setMethod(f = "reorder", 
          signature = signature("mvmreg", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_mvmreg(obs_stats,order)
          })


setMethod(f = "seed", 
          signature = signature("mvmreg","list","numeric"), 
          definition = function(model,data, K){
            X=cbind(data$X,data$Y)
            sds=apply(X,2,stats::sd)
            X=X[,sds!=0]
            X=t(t(X)/sds[sds!=0])
            km=stats::kmeans(X,K)
            km$cluster
          })



setMethod(f = "preprocess", 
          signature = signature("mvmreg"), 
          definition = function(model, data,K){
            list(Y=as.matrix(data),X=matrix(1,ncol=1,nrow=nrow(data)),N=nrow(data),moves=as.sparse(matrix(1,K,K)))
          })




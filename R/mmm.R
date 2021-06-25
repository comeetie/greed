#' @include models_classes.R fit_classes.R
NULL

#' @title Mixed mixture model description class
#' 
#' @description 
#' An S4 class to represent a Multinomial model model, extends \code{\link{icl_model-class}}.
#' Such model can be used to cluster a data matrix \eqn{X} with the following generative model :  
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i  \sim \mathcal{M}(1,\pi)}
#' \deqn{ \theta_{kv} \sim Dirichlet(\beta)}
#' ....
#' @slot name name of the model
#' @slot alpha Dirichlet over cluster proportions prior parameter (default to 1)
#' @slot beta Dirichlet over vocabulary prior parameter (default to 1)
#' @examples
#' new("mmm")
#' new("mmm",alpha=1,beta=1)
#' @export
setClass("mmm",
         representation = list(beta = "numeric",tau = "numeric",mu="numeric",epsilon="matrix",N0="numeric"),
         contains = "icl_model",
         prototype(name="mmm",beta=1,alpha=1,tau=0.001,N0=NaN,mu=NaN,epsilon=as.matrix(NaN)))

#' @title Mixed mixture model fit results class
#' 
#' @description
#'  An S4 class to represent a fit of a  extend \code{\link{icl_fit-class}}.
#' @slot model a \code{\link{lca-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size K*D with the number of occurrences of each modality for each clusters
#' }
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details depends on the training procedure)
#' @export 
setClass("mmm_fit",slots = list(model="mmm"),contains="icl_fit")


#' @title Mixed mixture model fit results class
#' 
#' 
#' @description An S4 class to represent a fit of a ...., extend \code{\link{icl_path-class}}.
#' @slot model a \code{\link{mm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size K*D with the number of occurrence of modality word in each clusters
#' }
#' @slot path a list of size K-1 with each part of the path described by:
#' \itemize{
#' \item icl1: icl value reach with this solution for alpha=1 
#' \item logalpha: log(alpha) value were this solution is better than its parent
#' \item K: number of clusters
#' \item cl: vector of cluster indexes
#' \item k,l: index of the cluster that were merged at this step
#' \item merge_mat: lower triangular matrix of delta icl values 
#' \item obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size K*D with the number of occurrence of modality word in each clusters
#' }
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy plotting with \code{ggplot2}
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father  
#' @slot train_hist  data.frame with training history information (details depends on the training procedure)
#' @export 
setClass("mmm_path",contains=c("icl_path","mmm_fit"))


#' @title plot a \code{\link{mmm_fit-class}} object
#' 
#' 
#' @param x a \code{\link{mmm_fit-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing connections between clusters
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export 
setMethod(f = "plot", 
          signature = signature("mmm_fit","missing"),
          definition = function(x,type='blocks'){
            ggplot2::ggplot()      
          });


#' @title plot a \code{\link{mmm_path-class}} object
#' 
#' @param x an \code{\link{mmm_path-class}} object
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'front'}: plot the extracted front in the plane ICL, log(alpha)
#' \item \code{'path'}: plot the evolution of ICL with respect to K
#' \item \code{'tree'}: plot the associated dendrogram
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export 
setMethod(f = "plot", 
          signature = signature("mmm_path","missing"),
          definition = function(x,type='blocks'){
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


#' @title Extract parameters from an \code{\link{mmm_fit-class}} object
#' 
#' @param object a \code{\link{mm_fit-class}}
#' @return a list with the model parameters estimates (MAP), the fields are:
#' \itemize{
#' \item \code{'pi'}: cluster proportions 
#' \item \code{'thetav'}: cluster profile probabilites (list of matrix of size K x Dv), 
#' }
#' @export 
setMethod(f = "coef", 
          signature = signature(object = "mmm_fit"),
          definition = function(object){
            list()
          })

reorder_lca = function(obs_stats,or){
  obs_stats[[2]]$counts = obs_stats[[2]]$counts[or]
  for(v in 1:length(obs_stats[[2]]$x_counts)){
    obs_stats[[2]]$x_counts[[v]] = obs_stats[[2]]$x_counts[[v]][or,]
  }
  obs_stats
}


setMethod(f = "reorder", 
          signature = signature("mmm", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_lca(obs_stats,order)
          })


setMethod(f = "seed", 
          signature = signature("mmm","list","numeric"), 
          definition = function(model,data, K){
            km=stats::kmeans(as.matrix(data$X),K)
            km$cluster
          })


setMethod(f = "sample_cl", 
          signature = signature("mmm","list","numeric"), 
          definition = function(model,data,K){
            sample(1:K,data$N,replace = TRUE)
          })

setMethod(f = "preprocess", 
          signature = signature("mmm"), 
          definition = function(model, data){
            if(!methods::is(data,"data.frame")){
              stop("An mmm model expect a data.frame.",call. = FALSE)
            }
            
            if(length(model@alpha)>1){
              stop("Model prior misspecification, alpha must be of length 1.",call. = FALSE)
            }
            if(is.na(model@alpha)){
              stop("Model prior misspecification, alpha is NA.",call. = FALSE)
            }
            if(model@alpha<=0){
              stop("Model prior misspecification, alpha must be positive.",call. = FALSE)
            }
            
            if(length(model@beta)>1){
              stop("Model prior misspecification, beta must be of length 1.",call. = FALSE)
            }
            if(is.na(model@beta)){
              stop("Model prior misspecification, beta is NA.",call. = FALSE)
            }
            if(model@beta<=0){
              stop("Model prior misspecification, beta must be positive.",call. = FALSE)
            }
            facts = sapply(1:ncol(data), function(col){is.factor(data[,col])})
            nums = sapply(1:ncol(data), function(col){is.numeric(data[,col])})
            
            if(!all(facts | nums)){
              stop("An mmm model expect a data.frame whith factors and numeric columns.",call. = FALSE)
            }
            print(list(Xcat=sapply(data[,facts],unclass)-1,Xnum =as.matrix(data[,nums]), N=nrow(data)))
            list(Xcat=sapply(data[,facts],unclass)-1,Xnum =as.matrix(data[,nums]), N=nrow(data))
          })

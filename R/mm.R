#' @include models_classes.R fit_classes.R
NULL

#' @title Clustering with a mixture of multinomial model description class
#' 
#' @description 
#' An S4 class to represent a multinomial model model, extend \code{\link{icl_model-class}}.
#' @slot name name of the model
#' @slot alpha Dirichlet over cluster proportions prior parameter (default to 1)
#' @slot beta Dirchlet over vocabulary prior parameter (default to 1)
#' @examples
#' new("mm")
#' new("mm",alpha=1,beta=1)
#' @export
setClass("mm",
         representation = list(beta = "numeric"),
         contains = "icl_model",
         prototype(name="mm",beta=1,alpha=1))

#' @title Clustering with a degree correted stochastick block model fit results class
#' 
#' @description
#'  An S4 class to represent a fit of a degree corrected stochastick block model for co_clustering, extend \code{\link{icl_fit-class}}.
#' @slot model a \code{\link{mm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and clolumns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size K*D with the number of occurence of modality word in each clusters
#' }
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history infromation (details depends on the training procedure)
#' @export 
setClass("mm_fit",slots = list(model="mm"),contains="icl_fit")


#' @title Clustering with a  mixture of multinomial model path extraction results class
#' 
#' 
#' @description An S4 class to represent a fit of a stochastick block model, extend \code{\link{icl_path-class}}.
#' @slot model a \code{\link{mm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and clolumns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size K*D with the number of occurence of modality word in each clusters
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
setClass("mm_path",contains=c("icl_path","mm_fit"))


#' @title plot a \code{\link{mm_fit-class}} object
#' 
#' 
#' @param x a \code{\link{mm_fit-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing connections between clusters
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export 
setMethod(f = "plot", 
          signature = signature("mm_fit","missing"),
          definition = function(x,type='blocks'){
            mat_blocks(x)      
          });


#' @title plot a \code{\link{mm_path-class}} object
#' 
#' @param x an \code{\link{mm_path-class}} object
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing 
#' \item \code{'nodelink'}: plot a nodelink diagram 
#' \item \code{'front'}: plot the extracted front in the plane ICL, log(alpha)
#' \item \code{'path'}: plot the veolution of ICL with repsect to K
#' \item \code{'tree'}: plot the associated dendogram
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export 
setMethod(f = "plot", 
          signature = signature("mm_path","missing"),
          definition = function(x,type='blocks'){
            switch(type,tree = {
              dendo(x)
            },
            path ={
              lapath(x)
            },
            front = {
              plot_front(x)
            },
            blocks ={
              methods::callNextMethod()
            },
            nodelink={
              methods::callNextMethod()
            })   
          })

reorder_mm = function(obs_stats,or){
  obs_stats$counts = obs_stats$counts[or]
  obs_stats$x_counts = obs_stats$x_counts[,or]
  obs_stats
}


setMethod(f = "reorder", 
          signature = signature("mm", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_mm(obs_stats,order)
          })


setMethod(f = "seed", 
          signature = signature("mm","list","numeric"), 
          definition = function(model,data, K){
            km=stats::kmeans(as.matrix(data$X),K)
            km$cluster
          })


setMethod(f = "sample_cl", 
          signature = signature("mm","list","numeric"), 
          definition = function(model,data,K){
            sample(1:K,data$N,replace = TRUE)
          })


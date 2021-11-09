#' @include models_classes.R fit_classes.R
NULL

#' @title Mixed Mixture model class
#' 
#' @description 
#' An S4 class to represent a Latent Class Analysis model, extends \code{\link{icl_model-class}}.
#' Such model can be used to cluster a data matrix \eqn{X} with the following generative model :  
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i  \sim \mathcal{M}(1,\pi)}
#' \deqn{ \theta_{kv} \sim Dirichlet(\beta)}
#' ....
#' @slot name name of the model
#' @slot alpha Dirichlet over cluster proportions prior parameter (default to 1)
#' @slot beta Dirichlet over vocabulary prior parameter (default to 1)
#' @examples
#' new("lca")
#' new("lca",alpha=1,beta=1)
#' @export
setClass("lca",
         representation = list(beta = "numeric"),
         contains = "icl_model",
         prototype(name="lca",beta=1,alpha=1))

#' @title Latent Class  Analysis fit results class
#'
#' @description An S4 class to represent a fit of a Latent Class Analysis model
#'   for categorical data clustering, extend \code{\link{icl_fit-class}}. The
#'   original data must be an n x p matrix where p is the number of variables
#'   and each variable is encoded as a factor (integer-valued).
#'
#' @slot model a \code{\link{lca-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with cluster indexes
#' @slot obs_stats a list with the following elements: \itemize{ \item counts:
#'   numeric vector of size K with number of elements in each clusters \item
#'   x_counts: matrix of size K*D with the number of occurrences of each
#'   modality for each clusters }
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details
#'   depends on the training procedure)
#' @export
setClass("lca_fit",slots = list(model="lca"),contains="icl_fit")


#' @title Latent Class Analysis hierarchical fit results class
#' 
#' 
#' @description An S4 class to represent a fit of a Latent Class Analysis model, extend \code{\link{icl_path-class}}.
#' @slot model a \code{\link{lca-class}} object to store the model fitted
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
setClass("lca_path",contains=c("icl_path","lca_fit"))


#' @title plot a \code{\link{lca_fit-class}} object
#' 
#' 
#' @param x a \code{\link{lca_fit-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing connections between clusters
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export 
setMethod(f = "plot", 
          signature = signature("lca_fit","missing"),
          definition = function(x,type='marginals'){
            pl=greed:::block_lca(x)
            grid::grid.newpage()
            gpl = grid::grid.draw(pl)
            invisible(pl)
          });


#' @title plot a \code{\link{lca_path-class}} object
#' 
#' @param x an \code{\link{lca_path-class}} object
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'front'}: plot the extracted front in the plane ICL, log(alpha)
#' \item \code{'path'}: plot the evolution of ICL with respect to K
#' \item \code{'tree'}: plot the associated dendrogram
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export 
setMethod(f = "plot", 
          signature = signature("lca_path","missing"),
          definition = function(x,type='marginals'){
            switch(type,tree = {
              dendo(x)
            },
            path ={
              lapath(x)
            },
            front = {
              plot_front(x)
            },invisible(methods::callNextMethod()))
            
          })


#' @title Extract parameters from an \code{\link{lca_fit-class}} object
#' 
#' @param object a \code{\link{mm_fit-class}}
#' @return a list with the model parameters estimates (MAP), the fields are:
#' \itemize{
#' \item \code{'pi'}: cluster proportions 
#' \item \code{'thetav'}: cluster profile probabilites (list of matrix of size K x Dv), 
#' }
#' @export 
setMethod(f = "coef", 
          signature = signature(object = "lca_fit"),
          definition = function(object){
            sol=object
            pi=(sol@obs_stats$counts+sol@model@alpha-1)/sum(sol@obs_stats$counts+sol@model@alpha-1)
            Thetak = lapply(sol@obs_stats$lca,function(mat){(mat+sol@model@beta-1)/rowSums(mat+sol@model@beta-1)})
            list(pi=pi,Thetak=Thetak)
})

reorder_lca = function(obs_stats,or){
  obs_stats[[2]]$counts = obs_stats[[2]]$counts[or]
  for(v in 1:length(obs_stats[[2]]$x_counts)){
    obs_stats[[2]]$x_counts[[v]] = obs_stats[[2]]$x_counts[[v]][or,]
  }
  obs_stats$counts = obs_stats$counts[or] 
  obs_stats
}

setMethod(f = "postprocess", 
          signature = signature("lca_path"), 
          definition = function(path,data,X,Y=NULL){
            path = clean_obs_stats_lca(path)
            path = name_obs_stats_lca(path,X)
            path
          }
          )

clean_obs_stats_lca = function(path) {
  # clean path@obs_stats to have clearer and non-redundant slots
  path@obs_stats = list(counts = path@obs_stats$counts,
                        lca = path@obs_stats[[2]]$x_counts)
  
  # Do the same thing for all submodels in the hierarchy path@path
  for(k in seq_along(path@path)) {
    path@path[[k]]$obs_stats = list(counts = path@path[[k]]$obs_stats$counts,
                                    lca = path@path[[k]]$obs_stats[[2]]$x_counts)
  }
  path
}


name_obs_stats_lca=function(path,X){
  cat_names = colnames(data.frame(X))
  for(v in 1:length(path@obs_stats$lca)){
    path@obs_stats$lca[[v]]=as.matrix(path@obs_stats$lca[[v]])
    colnames(path@obs_stats$lca[[v]])=levels(X[[v]])
    rownames(path@obs_stats$lca[[v]])=paste0("cluster",1:path@K)
  }
  names(path@obs_stats$lca)=cat_names
  
  for(p in 1:length(path@path)){
    for(v in 1:length(path@obs_stats$lca)){
      path@path[[p]]$obs_stats$lca[[v]]=matrix(path@path[[p]]$obs_stats$lca[[v]],nrow = path@path[[p]]$K)
      colnames(path@path[[p]]$obs_stats$lca[[v]])=levels(X[[v]])
      rownames(path@path[[p]]$obs_stats$lca[[v]])=paste0("cluster",1:path@path[[p]]$K)
    }
    names(path@path[[p]]$obs_stats$lca)=cat_names
  }
  path
}


setMethod(f = "reorder", 
          signature = signature("lca", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_lca(obs_stats,order)
          })


setMethod(f = "seed", 
          signature = signature("lca","list","numeric"), 
          definition = function(model,data, K){
            km=stats::kmeans(as.matrix(data$X),K)
            km$cluster
          })


setMethod(f = "sample_cl", 
          signature = signature("lca","list","numeric"), 
          definition = function(model,data,K){
            sample(1:K,data$N,replace = TRUE)
          })

setMethod(f = "preprocess", 
          signature = signature("lca"), 
          definition = function(model, data){
            if(!methods::is(data,"data.frame")){
              stop("An lca model expect a data.frame.",call. = FALSE)
            }
            if(!all(sapply(data,is.factor)) & !all(sapply(data, is.character))){
              stop("An lca model expect a data.frame with only factors or characters.",call. = FALSE)
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
            
            if(all(sapply(data, is.character))) {
              # Convert characters to factors
              data = data.frame(lapply(data, factor))
            }
            data = droplevels(data)
            list(X=sapply(data,unclass)-1,N=nrow(data))
          })

# For now, seed is kmeans on raw table for categorical
# May be improved
setMethod(f = "seed", 
          signature = signature("lca","list","numeric"), 
          definition = function(model,data, K){
            km=stats::kmeans(as.matrix(data$X),K)
            km$cluster
          })
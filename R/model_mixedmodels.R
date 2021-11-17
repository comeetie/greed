#' @include models_classes.R fit_classes.R
NULL



#' @title Mixed Models classes
#'
#' @description
#' An S4 class to represent a mixed clustering models.
#' Such model can be used to cluster graph vertex, and model a square adjacency matrix \eqn{X} with the following generative model :
#' @slot alpha Dirichlet over cluster proportions prior parameter (default to 10)
#' @slot p Exponential prior parameter (default to NaN, in this case p will be estimated from data as the mean connection probability)
#' @slot type define the type of networks (either "directed", "undirected" or "guess", default to "guess")
#' @family DlvmModels
#' @export
setClass("MixedModels",
  representation = list(models = "list"),
  contains = "DlvmPrior",
  prototype(alpha=1,models=list())
)


#' @describeIn MixedModels-class MixedModels class constructor
#' @examples
#' MixedModels(models = list(continuous=Gmm(),discrete=Lca()))
#' @export
MixedModels <- function(alpha = 1, models) {
  methods::new("MixedModels", alpha = alpha,models=models)
}



#' @title Mixed Models fit results class
#'
#' @description
#'  An S4 class to represent a fit of a degree corrected stochastic block model for co_clustering, extend \code{\link{IclFit-class}}.
#' @slot model a \code{\link{DcSbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details depends on the training procedure)
#' @export
setClass("MixedModelsFit", slots = list(model = "MixedModels"), contains = "IclFit")





#' @title Mixed Models hierarchical fit results class
#'
#'
#' @description An S4 class to represent a hierarchical fit of a degree corrected stochastic block model, extend \code{\link{IclPath-class}}.
#' @slot model a \code{\link{DcSbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' @slot path a list of size K-1 with each part of the path described by:
#' \itemize{
#' \item icl1: icl value reach with this solution for alpha=1
#' \item logalpha: log(alpha) value were this solution is better than its parent
#' \item K: number of clusters
#' \item cl: vector of cluster indexes
#' \item k,l: index of the cluster that were merged at this step
#' \item merge_mat: lower triangular matrix of delta icl values
#' \item obs_stats: a list with the elements:
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy plotting with \code{ggplot2}
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father
#' @slot train_hist  data.frame with training history information (details depends on the training procedure)
#' @export
setClass("MixedModelsPath", contains = c("IclPath", "MixedModelsFit"))





#' @title plot a \code{\link{MixedModelsPath-class}} object
#'
#' @param x an \code{\link{MixedModelsPath-class}} object
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'front'}: plot the extracted front in the plane ICL, log(alpha)
#' \item \code{'path'}: plot the evolution of ICL with respect to K
#' \item \code{'tree'}: plot the associated dendrogram
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export
setMethod(
  f = "plot",
  signature = signature("MixedModelsPath", "missing"),
  definition = function(x, type = "tree") {
    switch(type,
      tree = {
        dendo(x)
      },
      path = {
        lapath(x)
      },
      front = {
        plot_front(x)
      }
    )
  }
)


setMethod(
  f = "reorder",
  signature = signature("MixedModels", "list", "integer"),
  definition = function(model, obs_stats, order) {
    mnames=names(model@models)
    new_obs_stats = lapply(mnames,function(current_name){
        reorder(model@models[[current_name]],obs_stats[[current_name]],order)
      })
    names(new_obs_stats)=mnames
    new_obs_stats[["counts"]]=obs_stats$counts[order]
    new_obs_stats
  }
)


setMethod(
  f = "cleanObsStats",
  signature = signature("MixedModels", "list"),
  definition = function(model, obs_stats, data) {
    mnames=names(model@models)
    new_obs_stats = lapply(mnames,function(current_name){
      cleanObsStats(model@models[[current_name]],obs_stats[[current_name]],data[[current_name]])
    })
    names(new_obs_stats)=mnames
    new_obs_stats[["counts"]]=obs_stats$counts
    new_obs_stats
  }
)




setMethod(
  f = "preprocess",
  signature = signature("MixedModels"),
  definition = function(model, data) {
    
    mnames= names(model@models)
    if("counts" %in% mnames){
      stop("Prohibited models name counts is reserved, please use another model name.",.call=FALSE)
    }
    
    if(!(all(names(data) %in% mnames) & all(mnames %in% names(data)))){
      stop("Models names do notch match datasets names, please check the model and dataset list.",.call=FALSE)
    }
    data_prep = lapply(mnames,function(current_name){greed:::preprocess(model@models[[current_name]],data[[current_name]])})
    names(data_prep)=mnames
    Ns = sapply(data_prep,function(x){x$N})
    if(!all(Ns==Ns[1])){
      stop("All datasets must have the same number of elements.",.call=FALSE)
    }
    
    data_prep$N=Ns[1]
    
    data_prep
  }
)

setGeneric("extractSubModel", function(sol, sub_model_name, ...) standardGeneric("extractSubModel"))


#' @title extract a part of a \code{\link{MixedModelsPath-class}} object
#'
#' @param x an \code{\link{MixedModelsPath-class}} object
#' @param sub_model_name a string which specify the part of the model to extract
#' @return a \code{\link{IclFit}} object of the relevant class
#' @export
setMethod(
  f = "extractSubModel",
  signature = signature("MixedModelsPath", "character"),
  definition = function(sol, sub_model_name) {
    model <- sol@model@models[[sub_model_name]]
    model_type <- gsub("Prior", "", as.character(class(model)))
    model <- as(model, model_type)
    model@alpha <- sol@model@alpha
    sol <- as(as(sol, "IclFit"), paste0(model_type, "Fit"))
    sol@model <- model
    sol@obs_stats[[model_type]] <- sol@obs_stats[[sub_model_name]]
    sol@obs_stats <- sol@obs_stats[c("counts", model_type)]
    sol
  }
)

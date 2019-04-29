#' @title Clustering models classes
#' @aliases icl_model icl_model-class sbm sbm-class mm mm-class
#' @name models-classes
NULL


#' @rdname models-classes
#' @title icl_model
#' 
#' An S4 class to represent an abstract clustering model
#' \itemize{
#' \item slots : \code{name,alpha}
#' }
#' @slot name a character vector
#' @slot alpha a numeric vector of length 1 which define the parameters of the dirichlet over the cluster proportions (default to 1)
#' @export
setClass("icl_model",slots = list(name = "character",alpha = "numeric"))






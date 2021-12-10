#' @include models_classes.R
NULL


#' @title Abstract class to represent a clustering result
#'
#' @description
#' An S4 abstract class to represent an icl fit of a clustering model.
#'
#' @slot K a numeric vector of length 1 which correspond to the number of clusters
#' @slot icl a numeric vector of length 1 which store the the icl value
#' @slot cl a numeric vector of length N which store the clusters labels
#' @slot obs_stats a list to store the observed statistics of the model needed to compute ICL.
#' @slot obs_stats_cst a list to store the observed statistics of the model that do not depend on the clustering.
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist a data.frame to store training history (format depends on the used algorithm used).
#' @slot name generative model name
#' @export
setClass("IclFit", slots = list(name = "character", K = "numeric", obs_stats = "list", obs_stats_cst = "list", icl = "numeric", cl = "numeric", train_hist = "data.frame", move_mat = "dgCMatrix"))




#' @title  Abstract class to represent a hierarchical clustering result
#'
#' @description
#' An S4 class to represent a hierarchical path of solution.
#'
#' @slot path a list of merge moves describing the hierarchy of merge followed to complete totally the merge path.
#' @slot tree a tree representation of the merges.
#' @slot ggtree a data.frame for easy plotting of the dendrogram
#' @slot logalpha a numeric value which corresponds to the starting value of log(alpha).
#' @export
setClass("IclPath", slots = list(path = "list", tree = "numeric", ggtree = "data.frame", logalpha = "numeric"))

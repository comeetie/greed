#' @include models_classes.R
NULL

#' @title Clustering solutions classes
#' 
#' @name fits-classes
NULL

#' @rdname fits-classes
#' @title icl_fit
#' 
#' An S4 abstract class to represent an icl fit of a clustering model.
#' \itemize{
#' \item slots : \code{name,K,icl,cl,count}
#' }
#' @slot name of the fit
#' @slot K a numeric vector of length 1 which correspond to the number of clusters
#' @slot icl a numeric vector of length 1 which store the the icl value
#' @slot cl a numeric vector of length N which store the clusters labels
#' @slot obs_stats a list to store the observed statistics of the model needed to compute ICL.
#' @slot train_hist a data.frame to store training history (format depends on the used algorithm used).
#' @export 
setClass("icl_fit",slots = list(name="character",K="numeric",obs_stats="list",icl="numeric",cl="numeric",train_hist="data.frame",move_mat = "dgCMatrix"))






#' @rdname fits-classes
#' @title icl_path
#' 
#' An S4 class to represent a hierachical path of solution.
#' \itemize{
#' \item slots : \code{path,tree,ggtree,logalpha}
#' }
#' @slot path a list of merge moves describing the hierachie of merge followed to complete totaly the merge path.
#' @slot tree a tree representation of the merges.
#' @slot ggtree a data.frame for easy ploting of the dendogram
#' @slot logalpha a numeric value which corresponds to the starting value of log(alpha).
#' @export
setClass("icl_path",slots=list(path="list",tree="numeric",ggtree="data.frame",logalpha="numeric"))




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
#' @slot counts a numeric vector of length K which store the counts for each cluster
#' @export 
setClass("icl_fit",slots = list(name="character",K="numeric",counts="matrix",icl="numeric",cl="matrix"))

#' @rdname fits-classes
#' @title sbm_fit
#' 
#' An S4 class to represent a fit of a stochastick block model that extend \code{icl_fit}.
#' \itemize{
#' \item slots : \code{name,K,icl,cl,counts,x_counts,model}
#' }
#' @slot x_counts a numeric matrix to store the observed counts in the data for each cluster
#' @slot model an \code{\link{icl_model}} to store the model fitted
#' @export 
setClass("sbm_fit",slots = list(x_counts="matrix",model="sbm"),contains="icl_fit")


#' @rdname fits-classes
#' @title mm_fit
#' 
#' An S4 class to represent an icl fit of a mixture of multinomials model that extend \code{icl_fit}.
#' \itemize{
#' \item slots : \code{name,K,icl,cl,counts,x_counts,model}
#' }
#' @export 
setClass("mm_fit",slots = list(x_counts="matrix",model="mm"),contains="icl_fit")

#' @rdname fits-classes
#' @title mm_path
#' 
#' An S4 class to represent a hierachical path of stochastick block models that extend \code{icl_fit}.
#' \itemize{
#' \item slots : \code{name,K,icl,cl,counts,x_counts,path}
#' }
#' @slot path a list of merge moves describing the hierachie of merge followed to complete totaly the merge path.
#' @export
setClass("mm_path",slots = list(x_counts="matrix",path="list"),contains="icl_fit")

#' @rdname fits-classes
#' @title sbm_path
#' 
#' An S4 class to represent a hierachical path pf mixture of multinomials models that extend \code{icl_fit}.
#' \itemize{
#' \item slots : \code{name,K,icl,cl,counts,x_counts,path}
#' }
#' @export 
setClass("sbm_path",slots = list(x_counts="matrix",path="list"),contains="icl_fit")
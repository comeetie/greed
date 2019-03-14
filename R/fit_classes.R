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
setClass("icl_fit",slots = list(name="character",K="numeric",obs_stats="list",icl="numeric",cl="numeric",train_hist="data.frame"))

#' @rdname fits-classes
#' @title sbm_fit
#' 
#' An S4 class to represent a fit of a stochastick block model that extend \code{icl_fit}.
#' \itemize{
#' \item slots : \code{name,K,icl,cl,obs_stats,model}
#' }
#' @slot model an \code{\link{icl_model}} to store the model fitted
#' @export 
setClass("sbm_fit",slots = list(model="sbm"),contains="icl_fit")


#' @rdname fits-classes
#' @title dcsbm_fit
#' 
#' An S4 class to represent a fit of a stochastick block model that extend \code{icl_fit}.
#' \itemize{
#' \item slots : \code{name,K,icl,cl,obs_stats,model}
#' }
#' @slot model an \code{\link{icl_model}} to store the model fitted
#' @export 
setClass("dcsbm_fit",slots = list(model="dcsbm"),contains="icl_fit")




#' @rdname fits-classes
#' @title mm_fit
#' 
#' An S4 class to represent an icl fit of a mixture of multinomials model that extend \code{icl_fit}.
#' \itemize{
#' \item slots : \code{name,K,icl,cl,obs_stats,model}
#' }
#' @export 
setClass("mm_fit",slots = list(model="mm"),contains="icl_fit")


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

#' @rdname fits-classes
#' @title mm_path
#' 
#' An S4 class to represent a hierachical path of solutions for a mixture of mutinomials model that extend \code{mm_fit-class} and \code{icl_path-class}.
#' \itemize{
#' \item slots : \code{name,K,icl,cl,obs_stats, model, path, tree, ggtree, logalpha}
#' }
#' @export
setClass("mm_path",contains=c("icl_path","mm_fit"))

#' @rdname fits-classes
#' @title sbm_path
#' 
#' An S4 class to represent a hierachical path of solutions for a SBM model that extend \code{sbm_fit-class} and \code{icl_path-class}.
#' \itemize{
#' \item slots : \code{name,K,icl,cl,obs_stats, model, path, tree, ggtree, logalpha}
#' }
#' @export 
setClass("sbm_path",contains=c("icl_path","sbm_fit"))


#' @rdname fits-classes
#' @title dcsbm_path
#' 
#' An S4 class to represent a hierachical path of solutions for a DC-SBM model that extend \code{dcsbm_fit-class} and \code{icl_path-class}.
#' \itemize{
#' \item slots : \code{name,K,icl,cl,obs_stats, model, path, tree, ggtree, logalpha}
#' }
#' @export
setClass("dcsbm_path",contains=c("icl_path","dcsbm_fit"))

#' @include models_classes.R fit_classes.R
NULL


#' @title sbm
#' 
#' An S4 class to represent a stochastick block model that extends \code{icl_model} class.
#' @slot alpha dirichlet parameter for the prior over cluyster proportions (default to 1) 
#' @slot a0 a numeric vector of length 1 which define the parameters of the beta prior over the edges (default to 1)
#' @slot b0 a numeric vector of length 1 which define the parameters of the beta prior over the non-edges (default to 1)
#' @examples 
#' new("sbm")
#' new("sbm",a0=0.5,b0=0.5,alpha=1)
#' @export
setClass("sbm",
         representation = list(a0 = "numeric",b0="numeric"),
         contains = "icl_model",
         prototype(name="sbm",a0=1,b0=1,alpha=1))



#' @title sbm_fit
#' 
#' An S4 class to represent a fit of a stochastick block model that extend \code{icl_fit}.
#' @slot model an \code{\link{sbm}} to store the model fitted
#' @export 
setClass("sbm_fit",slots = list(model="sbm"),contains="icl_fit")


#' @title sbm_path
#' 
#' An S4 class to represent a hierachical path of solutions for a SBM model that extend \code{sbm_fit-class} and \code{icl_path-class}.
#' @slot name name of the model
#' @slot K number of extracted cluster
#' @slot icl value of the icl criterion
#' @slot cl vector of cluster label
#' @slot obs_stats list of model observed statistics (counts: number of nodes per clusters, x_counts: matrix with number of edges between each cluster)
#' @slot model \code{\link{sbm}} model used with it's prior parameters
#' @slot path list of clustering solution build using hierachical greedy merge
#' @slot tree merge tree 
#' @slot ggtree data.frame with merge tree meta information
#' @slot logalpha value of log(alpha) were alpha is the dirichlet prior parameter on cluster proportions
#' @slot train_hist data.frame with training history information
#' @slot move_mat sparse.matrix with eventually the constraints used for merge and swap moves 
#' @export 
setClass("sbm_path",contains=c("icl_path","sbm_fit"))


reorder_sbm = function(obs_stats,or){
  obs_stats$counts = obs_stats$counts[or]
  obs_stats$x_counts = obs_stats$x_counts[or,or]
  obs_stats
}


setMethod(f = "reorder", 
          signature = signature("sbm", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_sbm(obs_stats,order)
          })

setMethod(f = "seed", 
          signature = signature("sbm","list","numeric"), 
          definition = function(model,data, K){
            spectral(data$X,K)
          })


#' @rdname plot
#' @export
setMethod(f = "plot", 
          signature = signature("sbm_fit","missing"),
          definition = function(x,type="blocks"){
            switch(type,blocks=graph_blocks(x),nodelink=nodelink(x))
          });

#' @rdname plot
#' @export
setMethod(f = "plot", 
          signature = signature("sbm_path","missing"),
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

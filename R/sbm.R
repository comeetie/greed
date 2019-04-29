#' @include models_classes.R fit_classes.R
NULL


#' @rdname models-classes
#' @title sbm
#' 
#' An S4 class to represent a stochastick block model that extends \code{icl_model} class.
#' \itemize{
#' \item slots : \code{name,alpha,a0,b0}
#' }
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
#' @title sbm_path
#' 
#' An S4 class to represent a hierachical path of solutions for a SBM model that extend \code{sbm_fit-class} and \code{icl_path-class}.
#' \itemize{
#' \item slots : \code{name,K,icl,cl,obs_stats, model, path, tree, ggtree, logalpha}
#' }
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
          signature = signature("sbm","list","integer"), 
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

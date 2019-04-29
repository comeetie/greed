#' @include models_classes.R fit_classes.R
NULL


#' @rdname models-classes
#' @title dcsbm
#' 
#' An S4 class to represent a stochastick block model that extends \code{icl_model} class.
#' \itemize{
#' \item slots : \code{name,alpha,a0,b0}
#' }
#' @slot a0 a numeric vector of length 1 which define the parameters of the beta prior over the edges (default to 1)
#' @slot b0 a numeric vector of length 1 which define the parameters of the beta prior over the non-edges (default to 1)
#' @examples 
#' new("dcsbm")
#' @export
setClass("dcsbm",
         representation = list(),
         contains = "icl_model",
         prototype(name="dcsbm",alpha=1))


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
#' @title dcsbm_path
#' 
#' An S4 class to represent a hierachical path of solutions for a DC-SBM model that extend \code{dcsbm_fit-class} and \code{icl_path-class}.
#' \itemize{
#' \item slots : \code{name,K,icl,cl,obs_stats, model, path, tree, ggtree, logalpha}
#' }
#' @export
setClass("dcsbm_path",contains=c("icl_path","dcsbm_fit"))

reorder_dcsbm = function(obs_stats,or){
  obs_stats$counts = obs_stats$counts[or]
  obs_stats$din = obs_stats$din[or]
  obs_stats$dout = obs_stats$dout[or]
  obs_stats$x_counts = obs_stats$x_counts[or,or]
  obs_stats
}




setMethod(f = "reorder", 
          signature = signature("dcsbm", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_dcsbm(obs_stats,order)
          })


setMethod(f = "seed", 
          signature = signature("dcsbm","list","integer"), 
          definition = function(model,data, K){
            spectral(data$X,K)
          })


#' @rdname plot
#' @export
setMethod(f = "plot", 
          signature = signature("dcsbm_fit","missing"),
          definition = function(x,type="blocks"){
            switch(type,blocks=graph_blocks(x),nodelink=nodelink(x))
          })

#' @rdname plot
#' @export
setMethod(f = "plot", 
          signature = signature("dcsbm_path","missing"),
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

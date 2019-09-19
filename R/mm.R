#' @include models_classes.R fit_classes.R
NULL


#' @rdname models-classes
#' @title mm
#' 
#' An S4 class to represent a mixture of multinomial also known has mixture of unigrams that extends \code{icl_model} class.
#' \itemize{
#' \item slots : \code{name,alpha,beta}
#' }
#' @slot beta a numeric vector of length 1 which define the parameters of the beta prior over the counts (default to 1)
#' @examples
#' new("mm")
#' new("mm",alpha=1,beta=1)
#' @export
setClass("mm",
         representation = list(beta = "numeric"),
         contains = "icl_model",
         prototype(name="mm",beta=1,alpha=1))

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
#' @title mm_path
#' 
#' An S4 class to represent a hierachical path of solutions for a mixture of mutinomials model that extend \code{mm_fit-class} and \code{icl_path-class}.
#' \itemize{
#' \item slots : \code{name,K,icl,cl,obs_stats, model, path, tree, ggtree, logalpha}
#' }
#' @export
setClass("mm_path",contains=c("icl_path","mm_fit"))


reorder_mm = function(obs_stats,or){
  obs_stats$counts = obs_stats$counts[or]
  obs_stats$x_counts = obs_stats$x_counts[or,]
  obs_stats
}


setMethod(f = "reorder", 
          signature = signature("mm", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_mm(obs_stats,order)
          })


setMethod(f = "seed", 
          signature = signature("mm","list","numeric"), 
          definition = function(model,data, K){
            km=stats::kmeans(as.matrix(data$X),K)
            km$cluster
          })


setMethod(f = "sample_cl", 
          signature = signature("mm","list","numeric"), 
          definition = function(model,data,K){
            sample(1:K,data$N,replace = TRUE)
          })

#' @rdname plot
#' @param x \code{\link{icl_fit-class}} object to be ploted
#' @param type type of desired graphics : tree,pathy, blocks, nodelink, front
#' @return a ggplot2 object to visualize the results
#' @export
setMethod(f = "plot", 
          signature = signature("mm_fit","missing"),
          definition = function(x,type='blocks'){
            mat_blocks(x)      
          });


#' @rdname plot
#' @export
setMethod(f = "plot", 
          signature = signature("mm_path","missing"),
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
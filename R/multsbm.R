#' @include models_classes.R fit_classes.R
NULL


#' @title Multinomial Stochastic Block Model class
#' 
#' @description 
#' An S4 class to represent a Multinomial Stochastic Block Model, extends \code{\link{icl_model-class}}. Such model can be used to cluster multilayer-graph vertex, and model a square adjacency cube \eqn{X} of size NxNxM with the following generative model :  
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i  \sim \mathcal{M}(1,\pi)}
#' \deqn{ \theta_{kl} \sim Dirichlet(\beta)}
#' \deqn{ X_{ij.}|Z_{ik}Z_{jl}=1 \sim \mathcal{M}(L_{ij},\theta_{kl})}
#' With \eqn{L_{ij}=\sum_{m=1}^MX_{ijm}}. This class mainly store the prior parameters value \eqn{\alpha,\beta} of this generative model in the following slots:
#' @slot name name of the model
#' @slot alpha Dirichlet over cluster proportions prior parameter
#' @slot beta Dirichlet prior parameter over Multinomial links
#' @slot type define the type of networks (either "directed" or "undirected", default to "directed")
#' @seealso \code{\link{multsbm_fit-class}}, \code{\link{multsbm_path-class}}
#' @export 
setClass("multsbm",
         representation = list(beta = "numeric",type="character"),
         contains = "icl_model",
         prototype(name="multsbm",alpha=1,beta=1,type="directed"))



#' @title Multinomial Stochastic Block Model fit results class
#' 
#' @description An S4 class to represent a fit of a Multinomial Stochastic Block Model, extend \code{\link{icl_fit-class}}.
#' @slot model a \code{\link{multsbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: cube of size KxKxM with the number of links between each pair of clusters 
#' }
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details depends on the training procedure)
#' @export 
setClass("multsbm_fit",slots = list(model="multsbm"),contains="icl_fit")


#' @title Multinomial Stochastic Block Model hierachical fit results class
#' 
#' 
#' @description An S4 class to represent a hierarchical fit of a Multinomial Stochastic Block Model, extend \code{\link{icl_path-class}}.
#' @slot model a \code{\link{multsbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size KxKxM with the number of links between each pair of clusters 
#' }
#' @slot path a list of size K-1 with each part of the path described by:
#' \itemize{
#' \item icl1: icl value reach with this solution for alpha=1 
#' \item logalpha: log(alpha) value were this solution is better than its parent
#' \item K: number of clusters
#' \item cl: vector of cluster indexes
#' \item k,l: index of the cluster that were merged at this step
#' \item merge_mat: lower triangular matrix of delta icl values 
#' \item obs_stats: a list with the elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size KxKxM with the number of links between each pair of clusters 
#' }
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy plotting with \code{ggplot2}
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father  
#' @slot train_hist  data.frame with training history information (details depends on the training procedure)
#' @export 
setClass("multsbm_path",contains=c("icl_path","multsbm_fit"))

#' @title plot a \code{\link{multsbm_fit-class}} object
#' 
#' 
#' @param x a \code{\link{multsbm_fit-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing connections between clusters
#' \item \code{'nodelink'}: plot a nodelink diagram of the graph summarizing connections between clusters
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export 
setMethod(f = "plot", 
          signature = signature("multsbm_fit","missing"),
          definition = function(x,type="blocks"){
            switch(type,blocks=graph_blocks_cube(x),nodelink=nodelink_cube(x))
          });


#' @title plot a \code{\link{sbm_path-class}} object
#' 
#' @param x an \code{\link{sbm_path-class}} object
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing connections between clusters
#' \item \code{'nodelink'}: plot a nodelink diagram of the bipartite graph summarizing connections between clusters
#' \item \code{'front'}: plot the extracted front ICL, log(alpha)
#' \item \code{'path'}: plot the evolution of ICL with respect to K
#' \item \code{'tree'}: plot the associated dendrogram
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export 
setMethod(f = "plot", 
          signature = signature("multsbm_path","missing"),
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

reorder_multsbm = function(obs_stats,or){
  obs_stats$counts = obs_stats$counts[or]
  obs_stats$x_counts = obs_stats$x_counts[or,or,]
  obs_stats
}


setMethod(f = "reorder", 
          signature = signature("multsbm", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_multsbm(obs_stats,order)
          })

setMethod(f = "seed", 
          signature = signature("multsbm","list","numeric"), 
          definition = function(model,data, K){
            # pas terrible a réflechir deplier a droite sur les slices ?
            km=stats::kmeans(data$X[,,1],K)
            km$cluster
          })


setMethod(f = "preprocess", 
          signature = signature("multsbm"), 
          definition = function(model, data){
            if(!methods::is(data,"array") ){
              stop("A multsbm model expect an array with 3 dimensions.",call. = FALSE)
            }
            if(length(dim(data))!=3){
              stop("A multsbm model expect an array with 3 dimensions.",call. = FALSE)
            }
            if(dim(data)[1]!=dim(data)[2]){
              stop("A multsbm model expect an array of 3 dimensions with as many rows as columns.",call. = FALSE)
            }
            if(!all(round(data)==data) || min(data)<0){
              stop("A multsbm model expect an array of 3 dimensions with as many rows as columns filled with postive integers.",call. = FALSE)
            }
            issym = all(sapply(1:dim(data)[3],function(d){ isSymmetric(data[,,d])}))
            if(model@type=="undirected" & !issym){
              stop("An undirected multsbm expect a symmetric array.",call. = FALSE)
            }
            selfloops = sapply(1:dim(data)[3],function(d){ sum(diag(data[,,d]))})
            if(model@type=="undirected" & !all(selfloops==0)) {
              for (d in 1:dim(data)[3]){
                diag(data[,,d])=0
              }
              warning("An undirected multsbm model does not allow self loops, self loops were removed from the graph.",call. = FALSE)
            }
            list(X=data,N=nrow(data))
          })


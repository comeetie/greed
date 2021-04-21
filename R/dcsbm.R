#' @include models_classes.R fit_classes.R
NULL



#' @title Degree Corrected Stochastic Block Model class
#' 
#' @description 
#' An S4 class to represent a degree corrected stochastic block model, extend \code{\link{icl_model-class}}.
#' Such model can be used to cluster graph vertex, and model a square adjacency matrix \eqn{X} with the following generative model :  
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i  \sim \mathcal{M}(1,\pi)}
#' \deqn{ \theta_{kl} \sim Exponential(p)}
#' \deqn{ \gamma_i^+,\gamma_i^- \sim \mathcal{U}(S_k)}
#' \deqn{ X_{ij}|Z_{ik}Z_{jl}=1 \sim \mathcal{P}(\gamma_i^+\theta_{kl}\gamma_j^-)}
#' The individuals parameters \eqn{\gamma_i^+,\gamma_i^-} allow to take into account the node degree heterogeneity. 
#' These parameters have uniform priors over simplex \eqn{S_k}. This class mainly store the prior parameters value \eqn{\alpha} of this generative model in the following slots (the prior parameter \eqn{p} is estimated from the data as the global average probability of connection between two nodes):
#' @slot name name of the model
#' @slot alpha Dirichlet over cluster proportions prior parameter
#' @slot type define the type of networks (either "directed" or "undirected", default to "directed")
#' @export 
setClass("dcsbm",
         representation = list(type="character"),
         contains = "icl_model",
         prototype(name="dcsbm",alpha=1,type="directed"))



#' @title Degree Corrected Stochastic Block Model fit results class
#' 
#' @description
#'  An S4 class to represent a fit of a degree corrected stochastic block model for co_clustering, extend \code{\link{icl_fit-class}}.
#' @slot model a \code{\link{dcsbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item din: numeric vector of size K which store the sums of in-degrees for each clusters
#' \item dout: numeric vector of size K which store the sums of out-degrees for each clusters 
#' \item x_counts: matrix of size K*K with the number of links between each pair of clusters 
#' }
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details depends on the training procedure)
#' @export 
setClass("dcsbm_fit",slots = list(model="dcsbm"),contains="icl_fit")





#' @title Degree Corrected Stochastic Block Model hierarchical fit results class
#' 
#' 
#' @description An S4 class to represent a hierarchical fit of a degree corrected stochastic block model, extend \code{\link{icl_path-class}}.
#' @slot model a \code{\link{dcsbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item din: numeric vector of size K which store the sums of in-degrees for each clusters
#' \item dout: numeric vector of size K which store the sums of out-degrees for each clusters 
#' \item x_counts: matrix of size K*K with the number of links between each pair of clusters 
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
#' \item din: numeric vector of size K which store the sums of in-degrees for each clusters
#' \item dout: numeric vector of size K which store the sums of out-degrees for each clusters 
#' \item x_counts: matrix of size K*K with the number of links between each pair of clusters 
#' }
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy plotting with \code{ggplot2}
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father  
#' @slot train_hist  data.frame with training history information (details depends on the training procedure)
#' @export 
setClass("dcsbm_path",contains=c("icl_path","dcsbm_fit"))


#' @title plot a \code{\link{sbm_fit-class}} object
#' 
#' 
#' @param x a \code{\link{sbm_fit-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing connections between clusters
#' \item \code{'nodelink'}: plot a nodelink diagram of the graph summarizing connections between clusters
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export 
setMethod(f = "plot", 
          signature = signature("dcsbm_fit","missing"),
          definition = function(x,type="blocks"){
            switch(type,blocks=graph_blocks(x),nodelink=nodelink(x))
          })


#' @title plot a \code{\link{sbm_path-class}} object
#' 
#' @param x an \code{\link{sbm_path-class}} object
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing connections between clusters
#' \item \code{'nodelink'}: plot a nodelink diagram of the graph summarizing connections between clusters
#' \item \code{'front'}: plot the extracted front in the plane ICL, log(alpha)
#' \item \code{'path'}: plot the evolution of ICL with respect to K
#' \item \code{'tree'}: plot the associated dendrogram
#' }
#' @return a \code{\link{ggplot2}} graphic
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
          signature = signature("dcsbm","list","numeric"), 
          definition = function(model,data, K){
            spectral(data$X,K)
          })

setMethod(f = "preprocess", 
          signature = signature("dcsbm"), 
          definition = function(model, data){
            if(!(methods::is(data,"dgCMatrix") | methods::is(data,"matrix"))){
              stop("An dcsbm model expect a matrix or a sparse (dgCMatrix) matrix.",call. = FALSE)
            }
            if(nrow(data)!=ncol(data)){
              stop("A dcsbm model expect a square matrix.",call. = FALSE)
            }
            if(!all(round(data)==data) || min(data)<0){
              stop("An dcsbm model expect an integer matrix with postive values.",call. = FALSE)
            }
            if(model@type=="undirected" & !isSymmetric(data)){
              stop("An undirected dcsbm model expect a symmetric matrix.",call. = FALSE)
            }
            if(model@type=="undirected" & sum(diag(data))!=0){
              diag(data)=0
              warning("An undirected dcsbm model does not allow self loops, self loops were removed from the graph.",call. = FALSE)
            }
            list(X=as.sparse(data),N=nrow(data))
          })




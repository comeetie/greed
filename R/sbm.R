#' @include models_classes.R fit_classes.R
NULL


#' @title Stochastic Block Model class
#' 
#' @description 
#' An S4 class to represent a Stochastic Block Model, extends \code{\link{icl_model-class}}. 
#' Such model can be used to cluster graph vertex, and model a square adjacency matrix \eqn{X} with the following generative model :  
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i  \sim \mathcal{M}(1,\pi)}
#' \deqn{ \theta_{kl} \sim Beta(a_0,b_0)}
#' \deqn{ X_{ij}|Z_{ik}Z_{jl}=1 \sim \mathcal{B}(\theta_{kl})}
#' This class mainly store the prior parameters value \eqn{\alpha,a_0,b_0} of this generative model in the following slots:
#' @slot name name of the model
#' @slot alpha Dirichlet over cluster proportions prior parameter (default to 1)
#' @slot a0 Beta prior parameter over links (default to 1)
#' @slot b0 Beta prior parameter over no-links (default to 1)
#' @slot type define the type of networks (either "directed" or "undirected", default to "directed"), for undirected graphs the adjacency matrix is supposed to be symmetric.
#' @seealso \code{\link{sbm_fit-class}},\code{\link{sbm_path-class}}
#' @seealso \code{\link{greed}}
#' @examples 
#' new("sbm")
#' new("sbm",a0=0.5, b0= 0.5,alpha=0.5)
#' sbm = rsbm(100,c(0.5,0.5),diag(2)*0.1+0.01)
#' sol = greed(sbm$x,model=new("sbm",a0=0.5, b0= 0.5,alpha=0.5))
#' @references Nowicki, Krzysztof and Tom A B Snijders (2001). “Estimation and prediction for stochastic block structures”. In:Journal of the American statistical association 96.455, pp. 1077–1087
#' @export 
setClass("sbm",
         representation = list(a0 = "numeric",b0="numeric",type="character"),
         contains = "icl_model",
         prototype(name="sbm",a0=1,b0=1,alpha=1,type="directed"))



#' @title Stochastic Block Model fit results class
#' 
#' @description An S4 class to represent a fit of a Stochastic Block Model, extend \code{\link{icl_fit-class}}.
#' @slot model a \code{\link{sbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over rows and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size K*K with the number of links between each pair of clusters 
#' }
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details depends on the training procedure)
#' @export 
setClass("sbm_fit",slots = list(model="sbm"),contains="icl_fit")


#' @title Stochastic Block Model hierarchical fit results class
#' 
#' @description An S4 class to represent a hierarchical fit of a stochastic block model, extend \code{\link{icl_path-class}}.
#' @slot model a \code{\link{sbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size K*K with the number of links between each pair of clusters 
#' }
#' @slot path a list of size K-1 with that store all the solutions along the path. Each element is a list with the following fields:
#' \itemize{
#' \item icl1: icl value reach with this solution for alpha=1 
#' \item logalpha: log(alpha) value were this solution is better than its parent
#' \item K: number of clusters
#' \item cl: vector of cluster indexes
#' \item k,l: index of the cluster that were merged at this step
#' \item merge_mat: lower triangular matrix of delta icl values 
#' \item obs_stats: a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size K*K with the number of links between each pair of clusters 
#' }
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy plotting with \code{ggplot2}
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father  
#' @slot train_hist  data.frame with training history information (details depends on the training procedure)
#' @export 
setClass("sbm_path",contains=c("icl_path","sbm_fit"))

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
#' @rdname plot-sbm_fit
#' @export 
setMethod(f = "plot", 
          signature = signature("sbm_fit","missing"),
          definition = function(x,type="blocks"){
            switch(type,blocks=graph_blocks(x),nodelink=nodelink(x))
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

#' @title Extract parameters from an \code{\link{sbm_fit-class}} object
#' 
#' @param object a \code{\link{sbm_fit-class}}
#' @return a list with the model parameters estimates (MAP), the fields are:
#' \itemize{
#' \item \code{'pi'}: cluster proportions 
#' \item \code{'thetakl'}: between clusters connections probabilities (matrix of size K x K)
#' }
#' @export 
setMethod(f = "coef", 
          signature = signature(object = "sbm_fit"),
          definition = function(object){
            sol=object
            pi=(sol@obs_stats$counts+sol@model@alpha-1)/sum(sol@obs_stats$counts+sol@model@alpha-1)
            thetakl=(sol@obs_stats$x_counts+sol@model@a0-1)/(t(t(sol@obs_stats$counts))%*%sol@obs_stats$counts+sol@model@a0+sol@model@b0-2)
            if(sol@model@type=="undirected"){
              diag(thethakl)=(diag(sol@obs_stats$x_counts)/2+sol@model@a0-1)/(sol@obs_stats$counts*(sol@obs_stats$counts-1)/2+sol@model@a0+sol@model@b0-2)
            }
            list(pi=pi,thetakl=thetakl)
          })


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

setMethod(f = "preprocess", 
          signature = signature("sbm"), 
          definition = function(model, data){
            if(!(methods::is(data,"dgCMatrix") | methods::is(data,"matrix")| methods::is(data,"data.frame"))){
              stop("An sbm model expect a data.frame, a matrix or a sparse (dgCMatrix) matrix.",call. = FALSE)
            }
            if(methods::is(data,"data.frame")){
              data=as.matrix(data)
            }
            if(nrow(data)!=ncol(data)){
              stop("An sbm model expect a square matrix.",call. = FALSE)
            }
            if(!all(round(data)==data) || min(data)!=0 || max(data)!=1){
              stop("An sbm model expect a binary matrix.",call. = FALSE)
            }
            if(model@type=="undirected" & !isSymmetric(data)){
              stop("An undirected sbm model expect a symmetric matrix.",call. = FALSE)
            }
            if(model@type=="undirected" & sum(diag(data))!=0){
              diag(data)=0
              warning("An undirected sbm model does not allow self loops, self loops were removed from the graph.",call. = FALSE)
            }
            if(length(model@alpha)>1){
              stop("Model prior misspecification, alpha must be of length 1.",call. = FALSE)
            }
            if(is.na(model@alpha)){
              stop("Model prior misspecification, alpha is NA.",call. = FALSE)
            }
            if(model@alpha<=0){
              stop("Model prior misspecification, alpha must be positive.",call. = FALSE)
            }
            if(length(model@a0)>1){
              stop("Model prior misspecification, a0 must be of length 1.",call. = FALSE)
            }
            if(model@a0<=0){
              stop("Model prior misspecification, a0 must be positive.",call. = FALSE)
            }
            if(length(model@b0)>1){
              stop("Model prior misspecification, b0 must be of length 1.",call. = FALSE)
            }
            if(model@b0<=0){
              stop("Model prior misspecification, b0 must be positive.",call. = FALSE)
            }
            
            if(!(model@type %in% c("directed","undirected"))){
              stop("Model prior misspecification, model type must directed or undirected.",call. = FALSE)
            }
            list(X=as.sparse(data),N=nrow(data))
          })



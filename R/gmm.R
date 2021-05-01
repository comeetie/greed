#' @include models_classes.R fit_classes.R mvmreg.R
NULL


#' @title Gaussian mixture model description class
#' 
#' @description 
#' An S4 class to represent a multivariate Gaussian mixture  model, extend \code{\link{icl_model-class}}. 
#' The model corresponds to the following generative model:
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i  \sim \mathcal{M}(1,\pi)}
#' \deqn{ V_k \sim \mathcal{W}(\varepsilon^{-1},n_0)}
#' \deqn{ \mu_k \sim \mathcal{N}(\mu,(\tau V_k)^{-1})}
#' \deqn{ X_{i}|Z_{ik}=1 \sim \mathcal{N}(\mu_k,V_{k}^{-1})}
#' with \eqn{\mathcal{W}(\varepsilon^{-1},n_0)} the Whishart distribution. 
#' @slot name name of the model
#' @slot alpha Dirichlet over cluster proportions prior parameter (default to 1)
#' @slot tau Prior parameter (inverse variance) default 0.01 
#' @slot N0 Prior parameter (pseudo count) should be > number of features (default to NaN, in this case it will be estimated from data as the number of columns of X)
#' @slot epsilon Prior parameter co-variance matrix prior (matrix of size D x D), (default to a matrix of NaN, in this case epsilon will be estimated from data and will corresponds to 0.1 times a diagonal matrix with the variances of the X columns)  
#' @slot mu Prior parameters for the means (vector of size D), (default to NaN, in this case mu will be estimated from the data and will be equal to the mean of X) 
#' @examples
#' new("gmm")
#' new("gmm",alpha=1,tau=0.1,N0=15)
#' @md
#' @references Bertoletti, Marco & Friel, Nial & Rastelli, Riccardo. (2014). Choosing the number of clusters in a finite mixture model using an exact Integrated Completed Likelihood criterion. METRON. 73. 10.1007/s40300-015-0064-5. #' 
#' @export
setClass("gmm", representation = list(tau = "numeric",mu="numeric",epsilon="matrix",N0="numeric"),
         contains = "icl_model",
         prototype(name="gmm",tau=0.001,N0=NaN,mu=NaN,epsilon=as.matrix(NaN),alpha=1))


#' @title Gaussian mixture model fit results class
#' 
#' @description An S4 class to represent a fit of a multivariate mixture of regression model, extend \code{\link{icl_fit-class}}.
#' @slot model a \code{\link{gmm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item regs: list of size $K$ with statistics for each clusters
#' }
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details depends on the training procedure)
#' @export 
setClass("gmm_fit",slots = list(model="gmm"),contains="icl_fit")


#' @title  Gaussian mixture model hierarchical fit results class
#' 
#' 
#' @description An S4 class to represent a hierarchical fit of a gaussian mixture model, extend \code{\link{icl_path-class}}.
#' @slot model a \code{\link{gmm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and clolumns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item regs: list of size $K$ with statistics for each clusters
#' }
#' @slot path a list of size K-1 with each part of the path described by:
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
#' \item regs: list of size $K$ with statistics for each clusters
#' } 
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy plotting with \code{ggplot2}
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father  
#' @slot train_hist  data.frame with training history information (details depends on the training procedure)
#' @export 
setClass("gmm_path",contains=c("icl_path","gmm_fit"))



#' @title plot a \code{\link{gmm_path-class}} object
#' 
#' 
#' @param x a \code{\link{gmm_path-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'front'}: plot the extracted front ICL, log(alpha)
#' \item \code{'path'}: plot the evolution of ICL with respect to K
#' \item \code{'tree'}: plot the associated dendrogram
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export 
setMethod(f = "plot", 
          signature = signature("gmm_path","missing"),
          definition = function(x,type='tree'){
            switch(type,tree = {
              dendo(x)
            },
            path ={
              lapath(x)
            },
            front = {
              plot_front(x)
            })
          })

#' @title Extract mixture parameters from \code{\link{gmm_fit-class}} object
#' 
#' @param object a \code{\link{gmm_fit-class}}
#' @return a list with the mixture parameters estimates (MAP), the fields are:
#' \itemize{
#' \item \code{'pi'}: cluster proportions 
#' \item \code{'muk'}: cluster means
#' \item \code{'Sigmak'}: cluster co-variance matrices
#' }
#' @export 
setMethod(f = "coef", 
          signature = signature(object = "gmm_fit"),
          definition = function(object){
            sol=object
            pi=(sol@obs_stats$counts+sol@model@alpha-1)/sum(sol@obs_stats$counts+sol@model@alpha-1)
            muk = lapply(sol@obs_stats$regs, function(r){(sol@model@tau*sol@model@mu+r$ng*r$m)/(sol@model@tau+r$ng)})
            Sigmak = lapply(sol@obs_stats$regs, function(r){
              mu = (sol@model@tau*sol@model@mu+r$ng*r$m)/(sol@model@tau+r$ng)
              Sc= r$S+t(r$m)%*%r$m-t(mu)%*%mu
              S = (Sc+sol@model@tau*t(mu-sol@model@mu)%*%(mu-sol@model@mu)+sol@model@epsilon)/(r$ng+sol@model@N0-length(mu))
              S
            })
            list(pi=pi,muk=muk,Sigmak=Sigmak)
          })
          




setMethod(f = "seed", 
          signature = signature("gmm","list","numeric"), 
          definition = function(model,data, K){
            km=stats::kmeans(zscore(data$X),K)
            km$cluster
          })



setMethod(f = "preprocess", 
          signature = signature("gmm"), 
          definition = function(model, data){
            if(methods::is(data,"matrix") | methods::is(data,"data.frame") | methods::is(data,"dgCMatrix")){
              X=as.matrix(data)  
            }else{
              stop(paste0("Unsupported data type: ", class(X) ," use a data.frame, a matrix, a sparse dgCMatrix."),call. = FALSE)
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
            if(length(model@tau)>1){
              stop("Model prior misspecification, tau must be of length 1.",call. = FALSE)
            }
            if(is.na(model@tau)){
              stop("Model prior misspecification, tau is NA.",call. = FALSE)
            }
            if(model@tau<=0){
              stop("Model prior misspecification, tau must be positive.",call. = FALSE)
            }
            
            if(length(model@N0)>1){
              stop("Model prior misspecification, N0 must be of length 1.",call. = FALSE)
            }
            if(!is.na(model@N0) & model@N0<ncol(X)){
              stop("Model prior misspecification, N0 must be > ncol(X).",call. = FALSE)
            }
            
            if(prod(dim(model@epsilon))!=1 | !all(is.nan(model@epsilon))){
              if(dim(model@epsilon)[1]!=ncol(X)|| dim(model@epsilon)[2]!=ncol(X)){
                stop("Model prior misspecification, the dimensions of epsilon are not compatible with the data.",call. = FALSE)
              }
            }
            if(!all(is.nan(model@mu)) && length(model@mu)!=ncol(X)){
              stop("Model prior misspecification, mu length is incompatible with the data.",call. = FALSE)
            }
            
            
            list(X=X,N=nrow(X))
          })

reorder_gmm = function(obs_stats,or){
  obs_stats$counts = obs_stats$counts[or]
  obs_stats$regs = obs_stats$regs[or]
  obs_stats
}


setMethod(f = "reorder", 
          signature = signature("gmm", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_gmm(obs_stats,order)
          })


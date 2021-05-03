#' @include models_classes.R fit_classes.R 
NULL


#' @title Diagonal Gaussian mixture model description class
#' 
#' @description 
#' An S4 class to represent a multivariate diagonal Gaussian mixture model, extend \code{\link{icl_model-class}}. 
#' The model corresponds to the following generative model:
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i  \sim \mathcal{M}(1,\pi)}
#' \deqn{ \lambda_k^{(d)} \sim \mathcal{G}(\kappa,\beta)}
#' \deqn{ \mu_k^{(d)} \sim \mathcal{N}(\mu,(\tau \lambda_k)^{-1})}
#' \deqn{ X_{i.}|Z_{ik}=1 \sim \mathcal{N}(\mu_k,\lambda_{k}^{-1})}
#' with \eqn{\mathcal{G}(\kappa,\beta)} the Gamma distribution with shape parameter \eqn{\kappa} and rate parameter \eqn{\beta}. 
#' @slot name name of the model
#' @slot alpha Dirichlet over cluster proportions prior parameter (default to 1)
#' @slot tau Prior parameter (inverse variance), (default 0.01) 
#' @slot kappa Prior parameter (gamma shape), (default to 1)
#' @slot beta Prior parameter (gamma rate), (default to NaN, in this case beta will be estimated from data as 0.1 time the mean of X columns variances) 
#' @slot mu Prior for the means (vector of size D), (default to NaN, in this case mu will be estimated from data as the mean of X)
#' @examples
#' new("diaggmm")
#' new("diaggmm",alpha=1,tau=0.1,beta=0.1)
#' @md
#' @references Bertoletti, Marco & Friel, Nial & Rastelli, Riccardo. (2014). Choosing the number of clusters in a finite mixture model using an exact Integrated Completed Likelihood criterion. METRON. 73. 10.1007/s40300-015-0064-5. #' 
#' @export
setClass("diaggmm", representation = list(tau = "numeric",kappa="numeric",beta="numeric",mu="numeric"),
         contains = "icl_model",
         prototype(name="diaggmm",tau=0.01,kappa=1,beta=NaN,mu=NaN,alpha=1))


#' @title Diagonal Gaussian mixture model fit results class
#' 
#' @description An S4 class to represent a fit of a multivariate diagonal Gaussian mixture model, extend \code{\link{icl_fit-class}}.
#' @slot model a \code{\link{diaggmm-class}} object to store the model fitted
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
setClass("diaggmm_fit",slots = list(model="diaggmm"),contains="icl_fit")


#' @title  Diagonal Gaussian mixture model hierarchical fit results class
#' 
#' 
#' @description An S4 class to represent a hierarchical fit of a diagonal gaussian mixture model, extend \code{\link{icl_path-class}}.
#' @slot model a \code{\link{diaggmm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
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
setClass("diaggmm_path",contains=c("icl_path","diaggmm_fit"))



#' @title plot a \code{\link{diaggmm_path-class}} object
#' 
#' 
#' @param x a \code{\link{diaggmm_path-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'front'}: plot the extracted front ICL, log(alpha)
#' \item \code{'path'}: plot the evolution of ICL with respect to K
#' \item \code{'tree'}: plot the associated dendrogram
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export 
setMethod(f = "plot", 
          signature = signature("diaggmm_path","missing"),
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

#' @title Extract mixture parameters from \code{\link{diaggmm_fit-class}} object
#' 
#' @param object a \code{\link{diaggmm_fit-class}}
#' @return a list with the mixture parameters estimates (MAP), the fields are:
#' \itemize{
#' \item \code{'pi'}: cluster proportions 
#' \item \code{'muk'}: cluster means
#' \item \code{'Sigmak'}: cluster co-variance matrices
#' }
#' @export 
setMethod(f = "coef", 
          signature = signature(object = "diaggmm_fit"),
          definition = function(object){
            sol=object
            pi=(sol@obs_stats$counts+sol@model@alpha-1)/sum(sol@obs_stats$counts+sol@model@alpha-1)
            muk = lapply(sol@obs_stats$regs, function(r){(sol@model@tau*sol@model@mu+r$ng*r$m)/(sol@model@tau+r$ng)})
            Sigmak = lapply(sol@obs_stats$regs, function(r){
              betan = sol@model@beta +0.5*r$S+(sol@model@tau*r$ng*(r$m-sol@model@mu)^2)/(2*sol@model@tau+r$ng)
              alphan = sol@model@kappa+r$ng/2
              dd=as.vector(betan/(alphan-1))
              if(length(dd)>1){
                mode =  diag(dd)
              }else{
                mode = as.matrix(dd)
              }
            })
            list(pi=pi,muk=muk,Sigmak=Sigmak)
          })


setMethod(f = "seed", 
          signature = signature("diaggmm","list","numeric"), 
          definition = function(model,data, K){
            km=stats::kmeans(zscore(data$X),K)
            km$cluster
          })



setMethod(f = "preprocess", 
          signature = signature("diaggmm"), 
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
            
            if(length(model@kappa)>1){
              stop("Model prior misspecification, kappa must be of length 1.",call. = FALSE)
            }
            if(is.na(model@kappa)){
              stop("Model prior misspecification, kappa is NA.",call. = FALSE)
            }
            if(model@kappa<=0){
              stop("Model prior misspecification, kappa must be positive.",call. = FALSE)
            }
            
            if(length(model@beta)>1){
              stop("Model prior misspecification, beta must be of length 1.",call. = FALSE)
            }
            if(!is.nan(model@beta) && model@beta<=0){
              stop("Model prior misspecification, beta must be positive.",call. = FALSE)
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
          signature = signature("diaggmm", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_gmm(obs_stats,order)
          })


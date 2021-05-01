#' @include models_classes.R fit_classes.R
NULL


#' @title Multivariate mixture of regression model description class
#' 
#' @description 
#' An S4 class to represent a multivariate mixture of regression model, extends \code{\link{icl_model-class}}.
#' The model follows [minka-linear](https://tminka.github.io/papers/minka-linear.pdf) .
#' The model corresponds to the following generative model:
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i  \sim \mathcal{M}(1,\pi)}
#' \deqn{ V_k \sim \mathcal{W}(\varepsilon^{-1},n_0)}
#' \deqn{ A_k \sim \mathcal{MN}(0,(V_k)^{-1},\tau X^{t}X)}
#' \deqn{ Y_{i.}|X_{i.}Z_{ik}=1 \sim \mathcal{N}(A_kx_{i.},V_{k}^{-1})}
#' with \eqn{\mathcal{W}(\epsilon^{-1},n_0)} the Whishart distribution and \eqn{\mathcal{MN}} the matrix-normal distribution. 
#' @slot name name of the model
#' @slot alpha Dirichlet over cluster proportions prior parameter (default to 1)
#' @slot tau Prior parameter (inverse variance) default 0.01 
#' @slot epsilon Covariance matrix prior parameter (default to NaN, in this case epsilon will be fixed to a diagonal variance matrix equal to 0.1 time the variance of the regression residuals with only one cluster.) 
#' @slot N0 Prior parameter (default to NaN, in this case N0 will be fixed equal to the number of columns of Y.)
#' @examples
#' new("mvmreg")
#' new("mvmreg",alpha=1,tau=0.1,N0=15)
#' @md
#' @export
setClass("mvmreg",
         representation = list(tau = "numeric",N0="numeric",epsilon ="matrix"),
         contains = "icl_model",
         prototype(name="mvmreg",tau=0.01,N0=NaN,epsilon=as.matrix(NaN),alpha=1))


#' @title Clustering with a multivariate mixture of regression model fit results class
#' 
#' @description An S4 class to represent a fit of a multivariate mixture of regression model, extend \code{\link{icl_fit-class}}.
#' @slot model a \code{\link{mvmreg-class}} object to store the model fitted
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
setClass("mvmreg_fit",slots = list(model="mvmreg"),contains="icl_fit")


#' @title Multivariate mixture of regression model hierarchical fit results class
#' 
#' 
#' @description An S4 class to represent a hierarchical fit of a multivariate mixture of regression model, extend \code{\link{icl_path-class}}.
#' @slot model a \code{\link{mvmreg-class}} object to store the model fitted
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
setClass("mvmreg_path",contains=c("icl_path","mvmreg_fit"))



#' @title plot a \code{\link{mvmreg_path-class}} object
#' 
#' 
#' @param x a \code{\link{mvmreg_path-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'front'}: plot the extracted front ICL, log(alpha)
#' \item \code{'path'}: plot the evolution of ICL with respect to K
#' \item \code{'tree'}: plot the associated dendrogram
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export 
setMethod(f = "plot", 
          signature = signature("mvmreg_path","missing"),
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

#' @title Extract mixture parameters from \code{\link{mvmreg_fit-class}} object
#' 
#' @param object a \code{\link{mvmreg_fit-class}}
#' @return a list with the mixture parameters estimates (MAP), the fields are:
#' \itemize{
#' \item \code{'pi'}: cluster proportions 
#' \item \code{'A'}: cluster regression matrix
#' \item \code{'Sigmak'}: cluster noise co-variance matrices
#' }
#' @export 
setMethod(f = "coef", 
          signature = signature(object = "mvmreg_fit"),
          definition = function(object){
            sol=object
            pi=(sol@obs_stats$counts+sol@model@alpha-1)/sum(sol@obs_stats$counts+sol@model@alpha-1)
            muk = lapply(sol@obs_stats$regs, function(r){1/(1+sol@model@tau)*r$mu})
            Sigmak = lapply(sol@obs_stats$regs, function(r){
              S = (r$Syx+sol@model@epsilon)/(r$n+sol@model@N0)
            })
            list(pi=pi,muk=muk,Sigmak=Sigmak)
          })

reorder_mvmreg = function(obs_stats,or){
  obs_stats$counts = obs_stats$counts[or]
  obs_stats$regs = obs_stats$regs[or]
  obs_stats
}


setMethod(f = "reorder", 
          signature = signature("mvmreg", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_mvmreg(obs_stats,order)
          })


setMethod(f = "seed", 
          signature = signature("mvmreg","list","numeric"), 
          definition = function(model,data, K){
            X=cbind(data$X,data$Y)
            sds=apply(X,2,stats::sd)
            X=X[,sds!=0]
            X=t(t(X)/sds[sds!=0])
            km=stats::kmeans(X,K)
            km$cluster
          })



setMethod(f = "preprocess", 
          signature = signature("mvmreg"), 
          definition = function(model, data){
            if(!any(class(data$Y)%in%c("data.frame","numeric","matrix"))){
              stop("Y must be a data.frame a numeric vector or a matrix.",call.=FALSE)
            }
            
            if(!any(class(data$X)%in%c("data.frame","numeric","matrix"))){
              stop("X must be a data.frame a numeric vector or a matrix.",call.=FALSE)
            }
            if(methods::is(data$Y,"numeric")|methods::is(data$Y,"data.frame")){
              data$Y=as.matrix(data$Y)
            }
            if(methods::is(data$X,"numeric")|methods::is(data$X,"data.frame")){
              data$X=as.matrix(data$X)
            }
            if(nrow(data$X)!=nrow(data$Y)){
              stop("Incomptible sizes between X and Y.")
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
            if(!is.na(model@N0) & model@N0<ncol(data$Y)){
              stop("Model prior misspecification, N0 must be > ncol(Y).",call. = FALSE)
            }
            
            
            if(prod(dim(model@epsilon))!=1 | !all(is.nan(model@epsilon))){
              if(dim(model@epsilon)[1]!=ncol(data$Y)|| dim(model@epsilon)[2]!=ncol(data$Y)){
                stop("Model prior misspecification, the dimensions of epsilon are not compatible with the data.",call. = FALSE)
              }
            }
            list(Y=as.matrix(data$Y),X=as.matrix(data$X),N=nrow(data$X))
          })



#' @include models_classes.R fit_classes.R
NULL


#' @rdname models-classes
#' @title mreg
#' 
#' An S4 class to represent a mixture of multinomial also known has mixture of unigrams that extends \code{icl_model} class.
#' \itemize{
#' \item slots : \code{name,alpha,beta,N0}
#' }
#' @slot beta a numeric vector of length 1 which define the variance parameter of the normal prior over the regression parameters (default to 0.1)
#' @slot N0 a numeric vector of length 1 which define the parameter a0 of the inverse gamma over the regression noise variance parameters (default to 1)
#' @examples
#' new("mreg")
#' new("mreg",alpha=1,reg=5,a0=0.5,b0=0.5)
#' @export
setClass("mvmreg",
         representation = list(beta = "numeric",N0="numeric"),
         contains = "icl_model",
         prototype(name="mvmreg",beta=0.01,N0=10,alpha=1))


#' @rdname fits-classes
#' @title mreg_fit
#' 
#' An S4 class to represent an icl fit of a mixture of multinomials model that extend \code{icl_fit}.
#' \itemize{
#' \item slots : \code{name,K,icl,cl,obs_stats,model}
#' }
#' @export 
setClass("mvmreg_fit",slots = list(model="mvmreg"),contains="icl_fit")


#' @rdname fits-classes
#' @title mreg_path
#' 
#' An S4 class to represent a hierachical path of solutions for a Mixture of Regression model that extend \code{mreg_fit-class} and \code{icl_path-class}.
#' \itemize{
#' \item slots : \code{name,K,icl,cl,obs_stats, model, path, tree, ggtree, logalpha}
#' }
#' @export
setClass("mvmreg_path",contains=c("icl_path","mvmreg_fit"))


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
          definition = function(model, data,K){
            list(Y=as.matrix(data),X=matrix(1,ncol=1,nrow=nrow(data)),N=nrow(data),moves=as.sparse(matrix(1,K,K)))
          })



#' @rdname plot
#' @export
setMethod(f = "plot", 
          signature = signature("mvmreg_path","missing"),
          definition = function(x,type='blocks'){
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

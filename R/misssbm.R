#' @include models_classes.R fit_classes.R sbm.R
NULL


#' @title Stochastic Block Model with sampling scheme class
#' 
#' @description 
#' An S4 class to represent a Stochastic Block Model with a sampling scheme for missing data, extend \code{\link{icl_model-class}}. 
#' Such model can be used to cluster graph vertex, and model a square adjacency matrix \eqn{X} with the following generative model :  
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i  \sim \mathcal{M}(1,\pi)}
#' \deqn{ \theta_{kl} \sim Beta(a_0,b_0)}
#' \deqn{ X_{ij}|Z_{ik}Z_{jl}=1 \sim \mathcal{B}(\theta_{kl})}
#' Missing value are supposed to be generated afterwards, the observation process correspond to a binary matrix \eqn{R} of the same size as \eqn{X}, 
#' with \eqn{R_{ij}=1} for observed entries and \eqn{R_{ij}=0} for missing ones. \eqn{R} may be supposed to be MAR:
#' \deqn{ \epsilon \sim Beta(a_{0obs},b_{0obs})}
#' \deqn{ R_{ij} \sim \mathcal{B}(\epsilon)}
#' this correspond to the "dyad" sampling scheme. But the sampling scheme can also be NMAR:
#' \deqn{ \epsilon_{kl} \sim Beta(a_{0obs},b_{0obs})}
#' \deqn{ R_{ij}|Z_{ik}Z_{jl}=1 \sim \mathcal{B}(\epsilon_{kl})}
#' this correspond to the "block-dyad" sampling scheme.
#' This class mainly store the prior parameters value \eqn{\alpha,a_0,b_0,a_{0obs},b_{0obs}} of this generative model in the following slots:
#' @slot name name of the model
#' @slot alpha Dirichlet over cluster proportions prior parameter (default to 1)
#' @slot a0 Beta prior parameter over links (default to 1)
#' @slot b0 Beta prior parameter over no-links (default to 1)
#' @slot type define the type of networks (either "directed" or "undirected", default to "directed")
#' @slot sampling define the sampling process (either "dyad" or "block-dyad" )
#' @slot sampling_priors define the sampling process priors parameters (list with \code{a0obs} and \code{b0obs} fields.)
#' @seealso \code{\link{misssbm_fit-class}},\code{\link{misssbm_path-class}}  
#' @seealso \code{\link{greed}}
#' @examples 
#' new("misssbm")
#' new("misssbm",a0=0.5, b0= 0.5,alpha=0.5,sampling="dyad",sampling_priors=list(a0obs = 2,b0obs = 1))
#' sbm = rsbm(100,c(0.5,0.5),diag(2)*0.1+0.01)
#' sbm$x[cbind(base::sample(1:100,50),base::sample(1:100,50))]=NA
#' sol = greed(sbm$x,model=new("misssbm",sampling="dyad"))
#' @references Nowicki, Krzysztof and Tom A B Snijders (2001). “Estimation and prediction for stochastic block structures”. In:Journal of the American statistical association 96.455, pp. 1077–1087
#' @export 
setClass("misssbm",
         representation = list(sampling="character",sampling_priors="list"),
         contains = "sbm",
         prototype(name="misssbm",a0=1,b0=1,alpha=1,type="directed",sampling="block-dyad",sampling_priors = list(a0obs=1,b0obs=1)))



#' @title Stochastic Block Model with sampling scheme fit results class
#' 
#' @description An S4 class to represent a fit of a Stochastic Block Model, extend \code{\link{icl_fit-class}}.
#' @slot model a \code{\link{misssbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over rows and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size K*K with the number of observed links between each pair of clusters 
#' \item x_counts_obs: matrix of size K*K with the number of observed dyads between each pair of clusters 
#' }
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details depends on the training procedure)
#' @export 
setClass("misssbm_fit",slots = list(model="misssbm"),contains="icl_fit")


#' @title Stochastic Block Model with sampling scheme hierarchical fit results class 
#' 
#' 
#' @description An S4 class to represent a hierarchical fit of a stochastic block model, extend \code{\link{icl_path-class}}.
#' @slot model a \code{\link{misssbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size K*K with the number of observed links between each pair of clusters 
#' \item x_counts_obs: matrix of size K*K with the number of observed dyads between each pair of clusters 
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
#' \item x_counts: matrix of size K*K with the number of observed links between each pair of clusters 
#' \item x_counts_obs: matrix of size K*K with the number of observed dyads between each pair of clusters 
#' }
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy plotting with \code{ggplot2}
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father  
#' @slot train_hist  data.frame with training history information (details depends on the training procedure)
#' @export 
setClass("misssbm_path",contains=c("icl_path","misssbm_fit"))

#' @title plot a \code{\link{misssbm_fit-class}} object
#' 
#' @param x a \code{\link{misssbm_fit-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing connections between clusters
#' \item \code{'nodelink'}: plot a nodelink diagram of the graph summarizing connections between clusters
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export 
setMethod(f = "plot", 
          signature = signature("misssbm_fit","missing"),
          definition = function(x,type="blocks"){
            switch(type,blocks=graph_blocks(x),nodelink=nodelink(x))
          });


#' @title plot a \code{\link{misssbm_path-class}} object
#' 
#' @param x an \code{\link{misssbm_path-class}} object
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing connections between clusters
#' \item \code{'nodelink'}: plot a nodelink diagram of the bipartite graph summarizing connections between clusters
#' \item \code{'front'}: plot the extracted front ICL, log(alpha)
#' \item \code{'path'}: plot the evolution of ICL with respect to K
#' \item \code{'tree'}: plot the dendrogram
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export 
setMethod(f = "plot", 
          signature = signature("misssbm_path","missing"),
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


#' @title Extract parameters from an \code{\link{misssbm_fit-class}} object
#' 
#' @param object a \code{\link{misssbm_fit-class}}
#' @return a list with the model parameters estimates (MAP), the fields are:
#' \itemize{
#' \item \code{'pi'}: cluster proportions 
#' \item \code{'thetakl'}: between clusters connection probabilities (array of size K x K),
#' \item \code{'epsilonkl'}: between clusters dyad observation probabilities (array of size K x K) for block-dyad sampling and double for dyad sampling,
#' }
#' @export 
setMethod(f = "coef", 
          signature = signature(object = "misssbm_fit"),
          definition = function(object){
            sol=object
            pi=(sol@obs_stats$counts+sol@model@alpha-1)/sum(sol@obs_stats$counts+sol@model@alpha-1)
            if(sol@mode@type=="undirected"){
              x_counts=sol@obs_stats$x_counts
              diag(x_counts)=diag(x_counts)/2
              x_counts_obs=sol@obs_stats$x_counts_obs
              diag(x_counts_obs)=diag(x_counts_obs)/2
            }else{
              x_counts=sol@obs_stats$x_counts
            }
            thetakl=x_counts+sol@model@a0-1
            thetakl = thetakl/(x_counts_obs+sol@model@a0+sol@model@b0-2)
            if(sol@mode@type=="directed"){
              if(sol@model@sampling=="dyad"){
                epsilonkl=(sum(x_counts_obs)+sol@model@sampling_priors$a0obs-1)/(sum(sol@obs_stats$counts)^2+sol@model@sampling_priors$a0obs+sol@model@sampling_priors$b0obs-2)
              }else{
                epsilonkl=(x_counts_obs+sol@model@sampling_priors$a0obs-1)/(t(t(sol@obs_stats$counts))%*%sol@obs_stats$counts+sol@model@sampling_priors$a0obs+sol@model@sampling_priors$b0obs-2)
              }
            }else{
              if(sol@model@sampling=="dyad"){
                nbdyads = sum(sol@obs_stats$counts)*(sum(sol@obs_stats$counts)-1)/2
                epsilonkl=(sum(x_counts_obs[lower.tri(x_counts_obs,1)])+sol@model@sampling_priors$a0obs-1)/(nbdyads+sol@model@sampling_priors$a0obs+sol@model@sampling_priors$b0obs-2)
              }else{
                nbdyads = t(t(sol@obs_stats$counts))%*%sol@obs_stats$counts
                diag(nbdyads)=(diag(nbdyads)-sol@obs_stats$counts)/2
                epsilonkl=(x_counts_obs+sol@model@sampling_priors$a0obs-1)/(nbdyads+sol@model@sampling_priors$a0obs+sol@model@sampling_priors$b0obs-2)
              }
            }            
            list(pi=pi,thetakl=thetakl,epislonkl=epsilonkl)
          })

reorder_misssbm = function(obs_stats,or){
  obs_stats$counts = obs_stats$counts[or]
  obs_stats$x_counts = obs_stats$x_counts[or,or]
  obs_stats$x_counts_obs = obs_stats$x_counts_obs[or,or]
  obs_stats
}


setMethod(f = "reorder", 
          signature = signature("misssbm", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_misssbm(obs_stats,order)
          })

setMethod(f = "seed", 
          signature = signature("misssbm","list","numeric"), 
          definition = function(model,data, K){
            spectral(data$X,K)
          })

setMethod(f = "preprocess", 
          signature = signature("misssbm"), 
          definition = function(model, data){
            if(!(methods::is(data,"dgCMatrix") | methods::is(data,"matrix"))){
              stop("An sbm model expect a matrix or a sparse (dgCMatrix) matrix.",call. = FALSE)
            }
            if(nrow(data)!=ncol(data)){
              stop("An sbm model expect a square matrix.",call. = FALSE)
            }
            if(!all(round(data)==data,na.rm = TRUE) || min(data,na.rm = TRUE)!=0 || max(data,na.rm = TRUE)!=1){
              stop("An sbm model expect a binary matrix.",call. = FALSE)
            }
            if(model@type=="undirected" & !isSymmetric(data)){
              stop("An undirected sbm model expect a symmetric matrix.",call. = FALSE)
            }
            
            X=data
            X[is.na(data)]=0
            ij = Matrix::which(!is.na(data),arr.ind = TRUE)
            Xobs=Matrix::sparseMatrix(ij[,1],ij[,2],x = 1,dims = c(nrow(data),nrow(data)))
            if(model@type=="undirected"){
              diag(Xobs)=0
              diag(X)=0
            }
            if(model@type=="undirected" & sum(diag(data))!=0){
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
            
            if(!(model@sampling %in% c("dyad","block-dyad"))){
              stop("Model prior misspecification, only dyad and block-dyad sampling are supported.",call. = FALSE)
            }
            spnames= sort(names(model@sampling_priors))
            if(length(spnames)!=2 | spnames[1]!="a0obs" | spnames[2]!="b0obs" ){
              stop("Model prior misspecification, sampling prior slot must be a list with a0obs and b0obs fields.",call. = FALSE)
            }
            if(length(model@sampling_priors$a0obs)>1){
              stop("Model prior misspecification, sampling prior a0obs must be of length 1.",call. = FALSE)
            }
            if(model@sampling_priors$a0obs<=0){
              stop("Model prior misspecification, sampling prior a0obs must be positive.",call. = FALSE)
            }
            if(length(model@sampling_priors$b0obs)>1){
              stop("Model prior misspecification, sampling prior b0obs must be of length 1.",call. = FALSE)
            }
            if(model@sampling_priors$b0obs<=0){
              stop("Model prior misspecification, sampling prior b0obs must be positive.",call. = FALSE)
            }
            
            if(!(model@type %in% c("directed","undirected"))){
              stop("Model prior misspecification, model type must directed or undirected.",call. = FALSE)
            }
            list(X=as.sparse(X),Xobs=Xobs,N=nrow(data))
          })



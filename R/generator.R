#' Generate graph adjacency matrix using a SBM
#'
#' \code{rsbm} returns the adjacency matrix and the cluster labels generated randomly unsing a Stochastick Block Model.
#'
#' This function take graph size, cluster proportions and connectivity matrix as input and sample a graph accordingly together with the clusters labels.
#'
#' @param N A numeric value the size of the graph to generate
#' @param pi A numeric vector of length K with clusters proportions. Must sum up to 1.
#' @param mu A numeric matrix of dim K x K with the connectivity pattern to generate. elements in [0,1].
#' @return A list with fields:
#' \itemize{
#' \item x: the graph adjacency matrix as a \code{dgCMatrix}
#' \item K: number of generated clusters
#' \item N: number of vertex
#' \item cl: vector of clusters labels
#' \item pi: clusters proportions
#' \item mu: connectivuty matrix
#' }
#' @examples
#' simu = rsbm(100,rep(1/5,5),diag(rep(0.1,5))+0.001)
#' x  = simu$x
#' xl = simu$cl
#' @export
rsbm = function (N,pi,mu){
  K  = length(pi)
  cl = sample(1:K,N,replace=TRUE,prob = pi)
  x  = matrix(stats::rbinom(N*N,1,mu[cbind(rep(cl,N),rep(cl,each=N))]),N,N)
  links = Matrix::which(x==1,arr.ind = TRUE)
  list(cl=cl, x = Matrix::sparseMatrix(links[,1],links[,2], x = rep(1,nrow(links))), K=K,N=N,pi=pi,mu=mu)
}


#' Generate graph adjacency matrix using a Multinomial Mixture
#'
#' \code{rmm} returns a count matrix and the cluster labels generated randomly unsig a Mixture of Multinomial model.
#'
#' It take the sample size, cluster proportions and emission matrix, and  as input and sample a graph accordingly together with the clusters labels.
#'
#' @param N A numeric value the size of the graph to generate
#' @param pi A numeric vector of length K with clusters proportions. Must sum up to 1.
#' @param mu A numeric matrix of dim k x D with the clusters patterns to generate, all elements in [0,1].
#' @param lambda A numeric value which specify the expectation for the row sums.
#' @return A list with fields:
#' \itemize{
#' \item x: the count matrix as a \code{dgCMatrix}
#' \item K: number of generated clusters
#' \item N: number of vertex
#' \item cl: vector of clusters labels
#' \item pi: clusters proportions
#' \item mu: connectivity matrix
#' \item lambda: expectation of row sums
#' }
#' @export
rmm = function (N,pi,mu,lambda){
  K   = length(pi)
  nbv = dim(mu)[2]
  cl  = sample(1:K,N,replace=TRUE,prob = pi)
  X   = matrix(0,N,nbv) 
  for (k in 1:K){
    X[cl==k]  = t(stats::rmultinom(sum(cl==k),stats::rpois(N,lambda),mu[k,]))
  }
  links = Matrix::which(X>0,arr.ind = TRUE)
  list(cl=cl, x = Matrix::sparseMatrix(links[,1],links[,2], x = X[links]), K=K,N=N,pi=pi,mu=mu, lambda=lambda)
}

#' Generate graph adjacency matrix using a degree corrected SBM
#'
#' \code{rmm} returns a count matrix and the cluster labels generated randomly unsig a Mixture of Multinomial model.
#'
#' It take the sample size, cluster proportions and emission matrix, and  as input and sample a graph accordingly together with the clusters labels.
#'
#' @param N A numeric value the size of the graph to generate
#' @param pi A numeric vector of length K with clusters proportions. Must sum up to 1.
#' @param mu A numeric matrix of dim K x K with the connectivity pattern to generate, elements in [0,1].
#' @param betain A numeric vector of length N which specify the degree-correction will be normalized per cluster during the generation.
#' @return A list with fields:
#' \itemize{
#' \item x: the count matrix as a \code{dgCMatrix}
#' \item K: number of generated clusters
#' \item N: number of vertex
#' \item cl: vector of clusters labels
#' \item pi: clusters proportions
#' \item mu: connectivity matrix
#' \item betain: normalized in-degree parameters 
#' \item betaout: normalized out-degree parameters
#' }
#' @export
rdcsbm = function (N,pi,mu,betain,betaout){
  K  = length(pi)
  cl = sample(1:K,N,replace=TRUE,prob = pi)
  betain_cl = aggregate(betain,list(cl),mean)$x
  betain = betain/betain_cl[cl]
  betaout_cl = aggregate(betaout,list(cl),mean)$x
  betaout = betaout/betaout_cl[cl]
  x  = matrix(stats::rpois(N*N,mu[cbind(rep(cl,N),rep(cl,each=N))]*rep(betain,N)*rep(betaout,each=N)),N,N)
  links = Matrix::which(x>0,arr.ind = TRUE)
  list(cl=cl, x = Matrix::sparseMatrix(links[,1],links[,2], x = x[links]), K=K,N=N,pi=pi,mu=mu)
}


#' Generate X and y with a mixture of regression model
#'
#' \code{rmreg} returns an X matrix, a y vector and the cluster labels generated randomly unsig a Mixture of regression model.
#'
#' It take the sample size, cluster proportions and regression parameters matrix and variance  as input accordingly
#'
#' @param N A numeric value the size of the graph to generate
#' @param pi A numeric vector of length K with clusters proportions. Must sum up to 1.
#' @param mu A numeric matrix of dim K x d with the regression parameters.
#' @param sigma A numeric of length 1 with the target conditional variance.
#' @return A list with fields:
#' \itemize{
#' \item X: the covariate matrix 
#' \item y: the target feature
#' \item K: number of generated clusters
#' \item N: sample size
#' \item cl: vector of clusters labels
#' \item pi: clusters proportions
#' \item mu: regression parameters
#' \item sigma: conditional variance
#' }
#' @export
rmreg = function (N,pi,mu,sigma,X=cbind(matrix(rnorm(N*(ncol(mu)-1)),N,ncol(mu)-1),rep(1,N))){
  K  = length(pi)
  cl = sample(1:K,N,replace=TRUE,prob = pi)
  yt = X%*%mu+rnorm(N,0,sigma)
  y  = yt[cbind(1:N,cl)]
  list(cl=cl, X = X,y=y, K=K,N=N,pi=pi,mu=mu,sigma=sigma)
}


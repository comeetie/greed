#' Compute the entropy of a discrete sample
#'
#' @param cl vector of discrete labels
#' @return the entropy of the sample   
#' @examples
#' cl =sample(2,500,replace=TRUE)
#' H(cl)
#' @export
H = function(cl){
  p=table(cl)/length(cl)
  p=p[p!=0]
  -sum(p*log(p))
}

#' Compute the mutual information of two discrete samples
#'
#' @param cl1 vector of discrete labels
#' @param cl2 vector of discrete labels
#' @return the mutual information between the two discrete samples   
#' @examples
#' cl1 =sample(2,500,replace=TRUE)
#' cl2 =sample(2,500,replace=TRUE)
#' MI(cl1,cl2)
#' @export
MI = function(cl1,cl2){
  pj=table(cl1,cl2)/length(cl1)
  pi=as.matrix(table(cl1),15,1)%*%table(cl2)/length(cl1)^2
  sum(pj[pj!=0]*log(pj[pj!=0])-pj[pj!=0]*log(pi[pj!=0]))
}

#' Compute the normalized mutual information of two discrete samples
#'
#' @param cl1 vector of discrete labels
#' @param cl2 vector of discrete labels
#' @return the normalized mutual information between the two discrete samples   
#' @examples
#' cl1 =sample(2,500,replace=TRUE)
#' cl2 =sample(2,500,replace=TRUE)
#' NMI(cl1,cl2)
#' @export
NMI = function(cl1,cl2){
  MI(cl1,cl2)/max(c(H(cl1),H(cl2)))
}




#' @title Regularized spectral clustering 
#' 
#' @description 
#' performs regularized spectral clustering of a sparse adjacency matrix
#' @references Tai Qin, Karl Rohe. Regularized Spectral Clustering under the Degree-Corrected Stochastic Block Model. Nips 2013.
#' @param X An adjacency matrix in sparse format (see the \code{Matrix} package)
#' @param K Desired number of cluster
#' @return cl Vector of cluster labels
#' @export
spectral= function(X,K){
  X     = X+Matrix::t(X)
  ij    = Matrix::which(X>0,arr.ind = T)
  X[ij] = 1
  nu    = sum(X)/dim(X)[1]
  D     = (Matrix::rowSums(X)+nu)^-0.5
  Ln    = Matrix::sparseMatrix(ij[,1],ij[,2],x=D[ij[,1]]*D[ij[,2]])
  V     = RSpectra::eigs(Ln,K)
  Xp    = V$vectors
  Xpn   = Xp/sqrt(rowSums(Xp)^2)
  Xpn[rowSums(Xp)==0,]=0
  km    = stats::kmeans(Xp,K)
  km$cluster
}


zscore = function(X){
  m=apply(X,2,mean)
  X=t(t(X)-m)
  s=apply(X,2,stats::sd)
  X=t(t(X)/s)
}



as.sparse = function(X){
  S = X
  if(methods::is(X,"matrix")){
    ij= which(X!=0,arr.ind=TRUE)
    S = Matrix::sparseMatrix(ij[,1],ij[,2],x = X[ij]) 
  }
  S
}

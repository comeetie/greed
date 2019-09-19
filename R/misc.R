#' Compute the entropy of a discrete sample
#'
#' @param cl vector of discrete labels
#' @return the entropie of the sample   
#' @examples
#' cl =sample(2,500,replace=TRUE)
#' H(cl)
#' @export
H = function(cl){
  p=table(cl)/length(cl)
  p=p[p!=0]
  -sum(p*log(p))
}

#' Compute the mutual information of two discretes samples
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

#' Compute the normalized mutual information of two discretes samples
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

zscore = function(X){
  m=apply(X,2,mean)
  X=t(t(X)-m)
  s=apply(X,2,sd)
  X=t(t(X)/s)
}

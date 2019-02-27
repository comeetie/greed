

rsbm = function (N,pi,mu){
  K  = length(pi)
  cl = sample(1:K,N,replace=TRUE,prob = pi)
  x  = matrix(rbinom(N*N,1,mu[cbind(rep(cl,N),rep(cl,each=N))]),N,N)
  links = Matrix::which(x==1,arr.ind = TRUE)
  list(cl=cl, x = Matrix::sparseMatrix(links[,1],links[,2], x = rep(1,nrow(links))), K=K,N=N,pi=pi,mu=mu)
}


plot_adjmat = function(x){
  links_mat = Matrix::which(x==1,arr.ind = TRUE)
  links = data.frame(from=links_mat[,1],to=links_mat[,2])
  ggplot2::ggplot(links) + ggplot2::geom_point(ggplot2::aes(y=from,x=to))
}





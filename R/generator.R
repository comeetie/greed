

rsbm = function (N,pi,mu){
  K  = length(pi)
  cl = sample(1:K,N,replace=TRUE,prob = pi)
  x  = matrix(rbinom(N*N,1,mu[cbind(rep(cl,N),rep(cl,each=N))]),N,N)
  links = Matrix::which(x==1,arr.ind = TRUE)
  list(cl=cl, x = Matrix::sparseMatrix(links[,1],links[,2], x = rep(1,nrow(links))), K=K,N=N,pi=pi,mu=mu)
}


rmm = function (N,pi,mu,lambda){
  K   = length(pi)
  nbv = dim(mu)[2]
  cl  = sample(1:K,N,replace=TRUE,prob = pi)
  X   = matrix(0,N,nbv) 
  for (k in 1:K){
    X[cl==k]  = t(rmultinom(sum(cl==k),rpois(N,lambda),mu[k,]))
  }
  links = Matrix::which(X>0,arr.ind = TRUE)
  list(cl=cl, x = Matrix::sparseMatrix(links[,1],links[,2], x = X[links]), K=K,N=N,pi=pi,mu=mu, lambda=lambda)
}

spectral = function(x,k){
  L = diag(1/sqrt(rowSums(x)))%*%(diag(rowSums(x))-sbm$x)%*%diag(1/sqrt(rowSums(x)))
  E = eigen(L,TRUE)
  kmeans(E$vectors[,1:k],k)
}



plot_adjmat = function(x){
  links_mat = Matrix::which(x==1,arr.ind = TRUE)
  links = data.frame(from=links_mat[,1],to=links_mat[,2])
  ggplot2::ggplot(links) + ggplot2::geom_point(ggplot2::aes(y=from,x=to))
}


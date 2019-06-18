library(greed)
library(Matrix)
load("~/Projets/greed/data-raw/X20news.rda")
ij=which(X20news>0,arr.ind = TRUE)
di=dim(X20news)
N=sum(di)
X =  sparseMatrix(i=c(ij[,1],ij[,2]+di[1]),j=c(ij[,2]+di[1],ij[,1]),x = c(X20news[ij],X20news[ij]))

system.time(greed:::fit_greed(new("dcsbm"),list(X=X),sample(1:60,N,replace = TRUE),nb_max_pass=10,verbose = TRUE))


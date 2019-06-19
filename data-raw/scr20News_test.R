library(greed)
library(Matrix)
load("~/Projets/greed/data-raw/X20news.rda")
ij=which(X20news>0,arr.ind = TRUE)
di=dim(X20news)
N=sum(di)
X =  sparseMatrix(i=c(ij[,1],ij[,2]+di[1]),j=c(ij[,2]+di[1],ij[,1]),x = c(X20news[ij],X20news[ij]))

sol=greed:::fit_greed(new("dcsbm"),list(X=X),sample(1:60,N,replace = TRUE),nb_max_pass=10,verbose = TRUE)
system.time(greed:::fit_greed(new("dcsbm"),list(X=X),sample(1:60,N,replace = TRUE),nb_max_pass=10,verbose = TRUE))

a = sparseMatrix(i=sample(1:2500,2000),j=rep(1,2000),x = sample(200,2000,replace = TRUE))
b = sparseMatrix(i=sample(1:25000,10000),j=rep(1,10000),x = sample(200,10000,replace = TRUE))
c=greed:::add_sppat(a,b)

icc=greed:::gsum_mm(cl-1,t(X20news),60)

cl=sample(1:10,2500,replace = TRUE)

sol1=greed:::fit_greed(new("mm"),list(X=X20news[1:2500,1:500]),sample(1:10,2500,replace = TRUE),verbose = TRUE,type = "swap")
sol2=greed:::fit_greed(new("mm"),list(X=X20news[1:2500,1:500]),sample(1:10,2500,replace = TRUE),verbose = TRUE,type = "swap")

ij=which(table(sol1@cl,sol2@cl)>0,arr.ind = TRUE)
M=matrix(0,nrow(ij),nrow(ij))
for(k in 1:sol1@K){
  M[ij[,1]==k,ij[,1]==k]=1
}
for(k in 1:sol2@K){
  M[ij[,2]==k,ij[,2]==k]=1
}
ij=which(tril(M,-1)>0,arr.ind = TRUE)
sM=sparseMatrix(ij[,1],ij[,2],x=rep(1,nrow(ij)))


lapply(which(diff(ij[,2])!=0)


sol@K

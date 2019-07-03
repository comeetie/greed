library(greed)
library(future)
plan(multisession)
library(Matrix)
load("~/Projets/greed/data-raw/X20news.rda")

sol=greed(X20news,K=40,alg=new("hybrid",Kmax=400,nb_max_gen=6))

deb=as.POSIXct(Sys.time())
model=new("co_dcsbm")
data=greed:::preprocess(model,X20news)
solutions=greed:::multi_swap(model,new("multistarts",nb_start=10),data,30,verbose=TRUE)
verbose=TRUE
Kmax=400
fimerge = function(ncl,merge_graph){ greed:::merge_cstr(model,data,ncl,merge_graph,verbose)}
fiswap = function(ncl,move_mat){ greed:::swap_cstr(model,data,ncl,move_mat,verbose = verbose)}


sol1=solutions[[1]]
sol2=solutions[[2]]


sol_12=fco(solutions[[1]],solutions[[2]],fimerge,fiswap,0,Kmax)

sol_12%<-%fco(solutions[[1]],solutions[[2]],fimerge,fiswap,0,Kmax)
sol_34%<-%fco(solutions[[3]],solutions[[4]],fimerge,fiswap,0,Kmax)
sol_56%<-%fco(solutions[[5]],solutions[[6]],fimerge,fiswap,0,Kmax)
sol_78%<-%fco(solutions[[7]],solutions[[8]],fimerge,fiswap,0,Kmax)
fin=as.POSIXct(Sys.time())

deb=as.POSIXct(Sys.time())

Kmax=1000
sol_1234=fco(sol_12,sol_34,fimerge,fiswap,0,Kmax)
sol_5678=fco(sol_56,sol_78,fimerge,fiswap,0,Kmax)
fin=as.POSIXct(Sys.time())



fco = function(sol1,sol2,fimerge,fiswap,pmutation,Kmax){
  # cartesian product on the z of the two solution
  #ncl = unclass(factor(paste(sol1@cl,sol2@cl)))
  sol = sol1
  perm=sample(1:sol2@K,sol2@K)
  cl2 = perm[sol2@cl]
  
  K2=sol2@K
  print(Sys.time())
  
  Ta=table(sol@cl,cl2)
  print(sum(Ta>0))
  nbclust = sol@K-colSums(t(apply(Ta,1,function(r){cumsum(r>0)}))==matrix(rowSums(Ta>0),nrow=nrow(Ta),ncol=ncol(Ta)))+cumsum(colSums(Ta>0))
  if(max(nbclust)>Kmax){
    K2 = which(nbclust>Kmax)[1]
    cl2[cl2>K2]=K2
    # K2 = K2-1 merge sol1 cluster residuals ?
  }
  ij=which(table(sol@cl,cl2)>0,arr.ind = TRUE)
  ncl = as.numeric(factor(paste(sol@cl,"_",cl2,sep=""),levels=paste(ij[,1],"_",ij[,2],sep="")))
  print(max(ncl))
  
  
  M=matrix(0,nrow(ij),nrow(ij))
  for(k in 1:sol1@K){
    M[ij[,1]==k,ij[,1]==k]=1
  }
  for(k in 1:K2){
    M[ij[,2]==k,ij[,2]==k]=1
  }
  ijAm=which(M>0,arr.ind = TRUE)
  ijAm=ijAm[ijAm[,1]!=ijAm[,2],]
  Am=Matrix::sparseMatrix(ijAm[,1],ijAm[,2],x=rep(1,nrow(ijAm)),dims = c(max(ncl),max(ncl)))
  move_mat = Am
  if(nrow(ijAm)>0){
    sol=fimerge(ncl,Matrix::tril(Am))
    move_mat =sol@move_mat+t(sol@move_mat);
    ncl = sol@cl
  }
  print(Sys.time())
  if(runif(1)<pmutation){
    
    sp_cl=sample(max(ncl),1)
    nclold=ncl
    ncl[ncl==sp_cl]=sample(c(sp_cl,max(ncl)+1),sum(ncl==sp_cl),replace=TRUE)
    
    if(max(ncl)>nrow(move_mat) & sum(ncl==sp_cl)>0){
      move_mat = cbind(move_mat,move_mat[,sp_cl])
      move_mat = rbind(move_mat,move_mat[sp_cl,])
      move_mat[sp_cl,max(ncl)]=1
      move_mat[max(ncl),sp_cl]=1
    }else{
      ncl=nclold
    }
    
    
  }
  move_mat =sol@move_mat+t(sol@move_mat);
  for(r in 1:nrow(move_mat)){
    if(sum(move_mat[r,]!=0)>10){
      merges = which(move_mat[r,]!=0)
      best_merges_row = order(move_mat[r,merges],decreasing = TRUE)[1:10]
      move_mat[r,setdiff(merges,merges[best_merges_row])]=0
    }
  }
  
  sol_swap= fiswap(ncl,move_mat)
  }
  print(Sys.time())
  sol
}







library(ggplot2)
ggplot(sol_comp@train_hist)+geom_boxplot(aes(x=generation,group=generation,y=icl))


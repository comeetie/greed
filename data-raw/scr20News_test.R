library(greed)

library(future)
plan(multisession)
library(Matrix)
load("~/Projets/greed/data-raw/X20news.rda")
data("Xvlegislature")
N=10000
K=90
pi=rep(1/K,K)
lambda  = 0.1
lambda_o = 0.005
Ks=10
mu = bdiag(lapply(1:(K/Ks), function(k){matrix(lambda_o,Ks,Ks)+diag(rep(lambda,Ks))}))+0.0001
sbm = rsbm(N,pi,mu)
image(sbm$mu)




deb=as.POSIXct(Sys.time())
mod=new("co_dcsbm")
data=greed:::preprocess(mod,X20news)
sols=greed:::multi_swap(mod,new("multistarts",nb_start=10),data,30,verbose=TRUE)
fin=as.POSIXct(Sys.time())
fin-deb
solutions=sols

fimerge = function(ncl,merge_graph){ greed:::merge_cstr(mod,data,ncl,merge_graph,TRUE)}
fiswap = function(ncl,ws,iclust){ greed:::fit_greed_cstr(mod,data,ncl,ws,iclust,type="swap",verbose = TRUE)}

deb=as.POSIXct(Sys.time())
sol_12=fco(solutions[[1]],solutions[[2]],fimerge,fiswap,0)
sol_34=fco(solutions[[3]],solutions[[4]],fimerge,fiswap,0)
sol_56=fco(solutions[[5]],solutions[[6]],fimerge,fiswap,0)
sol_78=fco(solutions[[7]],solutions[[8]],fimerge,fiswap,0)
sol_78=fco(solutions[[7]],solutions[[8]],fimerge,fiswap,0)
fin=as.POSIXct(Sys.time())

deb=as.POSIXct(Sys.time())
sol_1234=fco(sol_12,sol_34,fimerge,fiswap,0)
sol_5678=fco(sol_56,sol_78,fimerge,fiswap,0)
fin=as.POSIXct(Sys.time())


fco= function(sol1,sol2,fimerge,fiswap,pmutation){
  # cartesian product on the z of the two solution
  #ncl = unclass(factor(paste(sol1@cl,sol2@cl)))
  
  # matrix of possible merge
  ij=which(table(sol1@cl,sol2@cl)>0,arr.ind = TRUE)
  ncl = as.numeric(factor(paste(sol1@cl,"_",sol2@cl,sep=""),levels=paste(ij[,1],"_",ij[,2],sep="")))
  M=matrix(0,nrow(ij),nrow(ij))
  for(k in 1:sol1@K){
    M[ij[,1]==k,ij[,1]==k]=1
  }
  for(k in 1:sol2@K){
    M[ij[,2]==k,ij[,2]==k]=1
  }
  ijAm=which(M>0,arr.ind = TRUE)
  ijAm=ijAm[ijAm[,1]!=ijAm[,2],]
  Am=Matrix::sparseMatrix(ijAm[,1],ijAm[,2],x=rep(1,nrow(ijAm)),dims = c(max(ncl),max(ncl)))
  if(nrow(ijAm)>0){
    sol=fimerge(ncl,Am)
    move_mat =sol@move_mat;
  }else{
    sol=sol1;
    move_mat=Am
  }
  ncl = sol@cl
  iclust = 1:max(ncl)
  if(runif(1)<pmutation){
    
    sp_cl=sample(max(ncl),1)
    
    ncl[ncl==sp_cl]=sample(c(sp_cl,max(ncl)+1),sum(ncl==sp_cl),replace=TRUE)
    
    if(max(ncl)>nrow(move_mat) & sum(ncl==sp_cl)>0){
      move_mat = cbind(move_mat,move_mat[,sp_cl])
      move_mat = rbind(move_mat,move_mat[sp_cl,])
      move_mat[sp_cl,max(ncl)]=1
      move_mat[max(ncl),sp_cl]=1
    }else{
      ncl=sol@cl
    }
    
    
  }
  
  #sol= fiswap(ncl,move_mat)
  sol
}

library(ggplot2)
ggplot(sol_comp@train_hist)+geom_boxplot(aes(x=generation,group=generation,y=icl))


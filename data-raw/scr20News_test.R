library(greed)
library(Matrix)
load("~/Projets/greed/data-raw/X20news.rda")
data("Xvlegislature")
X20news=Xvlegislature$X
ij=which(X20news>0,arr.ind = TRUE)
di=dim(X20news)
N=sum(di)
X =  sparseMatrix(i=c(ij[,1],ij[,2]+di[1]),j=c(ij[,2]+di[1],ij[,1]),x = c(X20news[ij],X20news[ij]))

library(future)
plan("multiprocess")

sol=greed(Xvlegislature$X,K=30,alg=new("hybrid",pop_size=60))

deb=as.POSIXct(Sys.time())
sol=greed(X,K=30,alg=new("hybrid",pop_size=60),verbose = TRUE)
fin=as.POSIXct(Sys.time())

system.time({

})
library(future)
plan(multisession)
library(Matrix)

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
data=greed:::preprocess(mod,X)
sols=greed:::multi_swap(mod,new("multistarts",nb_start=10),data,30,verbose=TRUE)
fin=as.POSIXct(Sys.time())
fin-deb

solutions=sols
new_solutions = listenv::listenv()
deb=as.POSIXct(Sys.time())
pop_size=20
icls = sapply(solutions,function(s){s@icl})
icl_order = order(icls,decreasing = TRUE)
selected  = icl_order[1:(pop_size/2)]

selected_couples = matrix(selected[sample(1:length(selected),length(selected)*2,replace = TRUE)],ncol=2)

fimerge = function(ncl,merge_graph){ soltemp = greed:::merge_cstr(mod,data,ncl,merge_graph,TRUE);

greed:::fit_greed(mod,data,soltemp@cl,type="merge",verbose = TRUE)
  
  }
fiswap = function(ncl,ws,iclust){ greed:::fit_greed_cstr(mod,data,ncl,ws,iclust,type="swap",verbose = TRUE)}
for (i in 1:nrow(selected_couples)){
  new_solutions[[i]] %<-% fco(solutions[[selected_couples[i,1]]],solutions[[selected_couples[i,2]]],fimerge,fiswap,0.1)
}
solutions = c(solutions[selected],as.list(new_solutions))
icls = sapply(solutions,function(s){s@icl})
fin=as.POSIXct(Sys.time())

fin-deb





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

sol_F=fco(sol_1234,sol_5678,fimerge,fiswap,0)
fco = function(sol1,sol2,fimerge,fiswap,pmutation){
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
  ijAm=which(tril(M,-1)>0,arr.ind = TRUE)
  Am=Matrix::sparseMatrix(ijAm[,1],ijAm[,2],x=rep(1,nrow(ijAm)),dims = c(max(ncl),max(ncl)))
  
  print("Prep...")
  
  sol=fimerge(ncl,Am)
  
  print("swap")
  #if(runif(1)<pmutation){
  #  ncl = sol@cl
  #  sp_cl=sample(max(ncl),1)
  #  ws = as.numeric(ncl==sp_cl)
  
  #  ncl[ncl==sp_cl]=sample(c(sp_cl,max(ncl)+1),sum(ncl==sp_cl),replace=TRUE)
  
  #  iclust  = c(sp_cl,max(ncl))
  ncl=sol@cl
  iclust=1:max(ncl)
  #sol= fiswap(ncl,as.numeric(ncl%in%iclust),iclust)
  #}
  
  sol
}

library(ggplot2)
ggplot(sol_comp@train_hist)+geom_boxplot(aes(x=generation,group=generation,y=icl))


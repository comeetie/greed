library(future)
library(dplyr)
library(purrr)
data=greed:::preprocess(new("co_dcsbm"),Xvlegislature$X,20)
solg=greed(Xvlegislature$X)

solbi=biclust::biclust(as.matrix(Xvlegislature$X),biclust::BCSpectral())

clrow=rep(0,nrow(solbi@RowxNumber))
clcol=rep(0,ncol(solbi@NumberxCol))
ir=1
ic=1
for (i in 1:ncol(solbi@RowxNumber)){
  if(sum(clrow[solbi@RowxNumber[,i]])==0){
    clrow[solbi@RowxNumber[,i]]=ir
    ir=ir+1
  } 
  if(sum(clcol[solbi@NumberxCol[i,]])==0){
    clcol[solbi@NumberxCol[i,]]=ic
    ic=ic+1
  } 
}

cl_i=c(clrow,clcol+max(clrow))
sol_bi=fit_greed(new("co_dcsbm"),data,cl_i,"none")
sol_bi@icl


library(future)
plan(multisession)
SOLCO  = listenv::listenv()
for (nbr in 1:20){
  for (nbc in 1:20){
    SOLCO[[nbr+(nbc-1)*20]]%<-%cocluster(as.matrix(Xvlegislature$X),'contingency',nbcocluster = c(nbr,nbc))
  }
}


lr=lapply(SOLCO, function(a){c(a@ICLvalue,max(a@rowclass+1),max(a@colclass+1))}) %>% purrr::keep(function(a){!is.infinite(a[1])})
RES=do.call(rbind,lr)
il = RES[which.max(RES[,1]),2]+(RES[which.max(RES[,1]),3]-1)*20
solco=SOLCO[[il]]
cl_i=c(solco@rowclass+1,solco@colclass+max(solco@rowclass)+2)
sol_co=fit_greed(new("co_dcsbm"),data,cl_i,"none")
sol_co@icl

library(readr);
library(dplyr);
mnist_raw <- read_csv("https://pjreddie.com/media/files/mnist_train.csv", col_names = FALSE)
cl=mnist_raw[,1]
X=as.matrix(mnist_raw[,-1])
Xbin=X>125
ij=which(Xbin,arr.ind = TRUE)
Xbins=Matrix::sparseMatrix(ij[,1],ij[,2],x=1)
library(greed)
library(Matrix)

model=MoM()
N=2000
Kr=40
ii=sample(60000,N)
Xt=Xbins[ii,]
Xt[Xt==1]=30
Xt[Xt==0]=1
data=greed:::preprocess(model,Xt)
clr=sample(Kr,N,replace = TRUE)
res=microbenchmark::microbenchmark({greed:::fit_greed(model,data,clr,type = "swap",verbose = TRUE,nb_max_pass = 1)},times = 1L)
median(res$time/10^9)
time.df=tibble(time=c(6,14.67,32.8,53.1,77.68),N=c(500,1000,2000,3000,4000))
library(ggplot2)
ggplot(time.df)+geom_line(aes(x=N,y=time))+geom_point(aes(x=N,y=time))+xlim(c(0,4000))

sol=greed(Xbins[ii,],model=MoM(),alg=Seed(),K=80)

co=coef(sol)
lp=lapply(1:K(sol),\(k){
  mat=matrix(c(co$thetak[k,],rep(0,5)),28,28)
  ij=which(mat>0,arr.ind = TRUE)
  df=tibble(i=ij[,1],j=ij[,2],v=unlist(mat[ij]))
  gg =ggplot(df)+geom_tile(aes(x=i,y=-j,fill=v))+
    coord_equal()+
    xlim(c(0,29))+
    ylim(c(-29,0))+
    scale_fill_gradient(high="#ffffff",low="#000000",guide="none")+
    theme_void()+
    theme(plot.background = element_rect(fill = "#000000"))
  ggplot2::ggplotGrob(gg)
  })

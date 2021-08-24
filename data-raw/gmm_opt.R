library(future)
library(ggplot2)
library(microbenchmark)

N = 3000
K = 20
X = rbind(MASS::mvrnorm(N/3,c(-3,0),diag(2)),MASS::mvrnorm(N/3,c(0,3),diag(2)),MASS::mvrnorm(N/3,c(3,3),diag(2)))
nbvn = 1
Xn = cbind(X,MASS::mvrnorm(N,rep(0,nbvn),diag(rep(1,nbvn))) %*% matrix(runif(nbvn*nbvn),nbvn,nbvn))*5
model=new("diaggmm")
data = greed:::preprocess(model,Xn)
cli = sample(1:K,N,replace = TRUE)


sol=greed:::fit_greed(model,data,cli,type="swap")
sol=greed:::fit_greed(model,data,rep(1:30,each=N/30),type="merge")
sol=greed:::fit_greed(model,data,rep(c(1,2,3),each=N/3),type="none")


table(sol@cl,rep(c(1,2,3),each=N/3))

timming = microbenchmark({greed:::fit_greed(model,data,cli,type="swap",nb_max_pass = 1)},times = 2,unit = "s")
timming

K = 4
xy = cbind(rep(1:K,K),rep(1:K,each=K))
xy = cbind(runif(K^2)*10,runif(K^2)*10)
Nk= 400
X = do.call(rbind,lapply(1:nrow(xy), \(i){MASS::mvrnorm(Nk,xy[i,],runif(1)*0.3*diag(2))}))
plan(multisession)
sol=greed(X,model=new("gmm"),alg=new("seed"),K=60)
sol=greed(X,model=new("gmm"),alg=new("hybrid",pop_size=50,prob_mutation=0.5,nb_max_gen=20),K=40)
centers = do.call(rbind,lapply(sol@obs_stats[[2]], \(x){x$m}))
ggplot(data.frame(x=X[,1],y=X[,2]))+geom_point(aes(x=x,y=y,color=factor(sol@cl)),alpha=0.3)+
  geom_density2d(aes(x=x,y=y))+
  geom_point(data=data.frame(x=centers[,1],y=centers[,2]),aes(x=x,y=y),color="red",size=3)+
  geom_point(data=data.frame(x=xy[,1],y=xy[,2]),aes(x=x,y=y),color="blue")

model=new("gmm")
data = greed:::preprocess(model,X)
cli = sample(1:100,nrow(X),replace = TRUE)
sol=greed:::fit_greed(model,data,cli,type="swap")
ggplot(data.frame(x=X[,1],y=X[,2]))+geom_point(aes(x=x,y=y,color=factor(sol@cl)),alpha=0.3)
#lq      mean    median        uq       max neval
#0.2720953 0.2770488 0.2755549 0.2802737 0.2870223    10

#lq      mean    median        uq       max neval
#0.2143715 0.2180028 0.2170492 0.2216261 0.2238532    10

# lq      mean   median        uq       max neval
# 0.2765024 0.2781535 0.278051 0.2784919 0.2844987    10

# lq      mean    median        uq       max neval
# 0.307677 0.3111059 0.3103494 0.3126276 0.3210271    10

#        lq      mean    median        uq      max neval
# 0.7320863 0.7486045 0.7354881 0.7391555 0.871886    10


library(greed)
library(Matrix)
N  <- 2000           # Number of node
K  <- 15            # Number of cluster
pi <- rep(1/K,K)    # Clusters proportions 
lambda   <- 0.1     # Building the connectivity matrix template
lambda_o <- 0.01
Ks <- 3
mu <- bdiag(lapply(1:(K/Ks), function(k){
matrix(lambda_o,Ks,Ks)+diag(rep(lambda,Ks))}))+0.001
sbm <- rsbm(N,pi,mu) # Simulation
model=Sbm()
data=greed:::preprocess(model,sbm$x)


sol1=greed:::fit_greed(model,data,sample(20,N,replace=TRUE),type="swap")
table(sbm$cl,clustering(sol1))
plot(sol1,type="blocks")
sol2=greed:::fit_greed(model,data,sample(20,N,replace=TRUE),type="swap")
table(sbm$cl,clustering(sol2))
plot(sol2,type="blocks")
table(clustering(sol1),clustering(sol2))
cl_merge = as.numeric(factor(paste(clustering(sol1),clustering(sol2))))
sol_merge = greed:::fit_greed(model,data,cl_merge,type="none",verbose = TRUE)
sol_merge
plot(sol_merge,type="blocks")
table(sbm$cl,clustering(sol_merge))
sol_merge = greed:::fit_greed(model,data,cl_merge,type="merge",verbose = TRUE)
plot(sol_merge,type="blocks")
sol_merge
table(sbm$cl,clustering(sol_merge))

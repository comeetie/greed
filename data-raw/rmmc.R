rmmc = function(pi,As,Is,Ls,N){
  cl <- sample(1:K, N, replace = TRUE, prob = pi)
  chains = list()
  for (i in 1:length(cl)){
    k=cl[i]
    l = rpois(1,Ls[[k]])+2
    si=sample(1:length(Is[[k]]), 1, replace = TRUE, prob = Is[[k]])
    chain = c(si)
    for(ll in 2:l){
      pm=As[[k]][chain[ll-1],]
      si=sample(1:length(pm), 1, replace = TRUE, prob = pm)
      chain = c(chain,si)
    }
    chains = c(chains,list(chain))
  }
  list(chains=chains,cl=cl,pi=pi,As=As,Is=Is,Ls=Ls)
}


K  = 10
pi = rep(1/K,K)
As = lapply(1:K,\(k){
  A=matrix(runif(100),10,10)
  A[A>0.6]=15
  A/rowSums(A)
  })
Is = lapply(1:K,\(k){
  I=runif(10)
  I[I>0.6]=15
  I/sum(I)
})
Ls = rep(4,K)
N=1000

mmc= rmmc(pi,As,Is,Ls,N)


Xchains = mmc$chains
Xi = data.frame(sti = factor(sapply(Xchains,\(ch){ch[1]}),1:10))

Xtr=lapply(Xchains,\(chain){
ij=cbind(chain[1:(length(chain)-1)],chain[2:length(chain)])
At=table(factor(ij[,1],1:10),factor(ij[,2],1:10))
})

X = lapply(1:10,\(row){
  do.call(rbind,lapply(Xtr,\(At){
    At[row,]
  }))
})
names(X)=paste0("tr",1:10)
X$sti=Xi

library(greed)

models=lapply(1:10,\(v){MoMPrior()})
names(models)=paste0("tr",1:10)
models$sti=LcaPrior()

model=CombinedModels(models = models)
sol=greed(X,model=model)
table(sol@cl,mmc$cl)



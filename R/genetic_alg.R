#' @useDynLib greed
#' @importFrom Rcpp sourceCpp
#' @importFrom future %<-%
#' @name %<-%
NULL

#' @include models_classes.R fit_classes.R cleanpath.R
#' @import Matrix
NULL

genetic = function(model,alg,data,K,verbose=FALSE){
  
  
  init_f =  function(C){fit_greed(model,data,C,"none")}
  train.hist = data.frame(generation=c(),icl=c(),K=c())
  
  # multi-start in //
  #future::plan(future::multiprocess)
  
  solutions = listenv::listenv()
  # first generation of solutions
  pop_size = alg@pop_size
  for (i in 1:pop_size){
    cli = sample_cl(model,data,K)
    cli=as.numeric(factor(cli))
    solutions[[i]] %<-% fit_greed(model,data,cli,type="none")
  }
  solutions = as.list(solutions)
  icls  = sapply(solutions,function(s){s@icl})
  # check for errors 
  solutions=solutions[!is.nan(icls)]
  icls=icls[!is.nan(icls)]
  
  nbgen = 1
  
  
  while(max(icls)-min(icls)>0 && nbgen < alg@nb_max_gen){
    
    train.hist=rbind(train.hist,data.frame(generation=nbgen,icl=icls,K=sapply(solutions,function(s){max(s@cl)})))
    
    ii = order(icls)
    Nsel = round(alg@pop_size*alg@sel_frac)
    ii=ii[(alg@pop_size-Nsel):alg@pop_size]
    icls =icls[ii]
    solutions = solutions[ii]
    bres = solutions[[order(icls,decreasing = TRUE)[1]]]
    new_solutions = listenv::listenv()
    children = listenv::listenv()
    
    for (i in 1:(alg@pop_size-1)){
      ip = 1:Nsel
      i1 = sample(ip,1,prob=ip)
      i2 = sample(ip[-i1],1,prob=ip[-i1])
      children[[i]] = cross_over(solutions[[i1]],solutions[[i2]],init_f)
    }
    children = as.list(children)
    for (i in 1:(alg@pop_size-1)){
      if(stats::runif(1)<alg@prob_mutation){
        new_solutions[[i]] %<-% fit_greed(model,data,children[[i]]@cl,"swap",1)
      }else{
        new_solutions[[i]] = children[[i]]
      }
    }
    solutions = c(bres,as.list(new_solutions))
    
    icls = sapply(solutions,function(s){s@icl})
    nbgen = nbgen + 1;
  }
  
  train.hist=rbind(train.hist,data.frame(generation=nbgen,icl=icls,K=sapply(solutions,function(s){max(s@cl)})))
  
  # best solution
  res = solutions[[order(icls,decreasing = TRUE)[1]]]
  # compute merge path
  
  res = fit_greed(model,data,res@cl,"both")
  path = fit_greed_path(data,res)
  # clean the resuts (compute, merge tree,...)
  
  path = cleanpath(path)
  # store train history
  path@train_hist = train.hist
  path
}


cross_over = function(sol1,sol2,init_f){
  cl1 = sol1@cl
  K1 = sol1@K
  lk1 = 1:K1
  cl2 = sol2@cl
  K2=sol2@K
  lk2=1:K2
  C=rep(0,length(cl1))
  K=0
  while(sum(C==0)>0){
    if(K1>0  && stats::runif(1)>0.5){
      kt = sample(K1,1,1)
      k = lk1[kt]
      K1=K1-1
      lk1=lk1[-kt]
      ink = (cl1==k & C==0)
      if(sum(ink)>0){
        K=K+1
        C[ink]=K
      }
    }else if(K2>0){
      kt = sample(K2,1,1)
      k = lk2[kt]
      K2=K2-1
      lk2=lk2[-kt]
      ink = (cl2==k & C==0)
      if(sum(ink)>0){
        K=K+1
        C[ink]=K
      }
    }
    
  }
  init_f(C)
}






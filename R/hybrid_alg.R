#' @useDynLib greed
#' @importFrom Rcpp sourceCpp
#' @importFrom future %<-%
#' @name %<-%
NULL

#' @include models_classes.R fit_classes.R cleanpath.R
#' @import Matrix
NULL

hybrid = function(model,alg,data,K, verbose=FALSE){
            
            fimerge = function(ncl,merge_graph){ greed:::merge_cstr(model,data,ncl,merge_graph,verbose)}
            fiswap = function(ncl,ws,iclust){ fit_greed_cstr(model,data,ncl,ws,iclust,type="swap",verbose = verbose)}
            train.hist = data.frame(generation=c(),icl=c(),K=c())

            # multi-start in //
            #future::plan(future::multiprocess)

            solutions = listenv::listenv()
            # first generation of solutions
            pop_size = alg@pop_size
            for (i in 1:pop_size){
              cli = sample_cl(model,data,K)
              solutions[[i]] %<-% fit_greed(model,data,cli,verbose = verbose) %globals% c("model","data","cli","verbose","fit_greed")
            }
            solutions = as.list(solutions)
            icls  = sapply(solutions,function(s){s@icl})
            # check for errors 
            solutions=solutions[!is.nan(icls)]
            icls=icls[!is.nan(icls)]
            old_best = -Inf
            best_icl = max(icls)
            nbgen = 1
            # while maximum number of generation // all solutions are equals // no improvements
            print(icls)
            pmut = alg@prob_mutation
            while((max(icls)-min(icls))>1  & nbgen < alg@nb_max_gen){
              

              train.hist=rbind(train.hist,data.frame(generation=nbgen,icl=icls,K=sapply(solutions,function(s){max(s@cl)})))
              # selection keep the top half solutions
              ii = order(icls)
              Nsel = round(alg@pop_size*0.6)
              ii=ii[(alg@pop_size-Nsel):alg@pop_size]
              icls =icls[ii]
              solutions = solutions[ii]
              bres = solutions[[order(icls,decreasing = TRUE)[1]]]
              new_solutions = listenv::listenv()
              for (i in 1:(alg@pop_size-1)){
                ip = 1:min(c(Nsel,length(solutions)))
                i1 = sample(ip,1,prob=ip)
                i2 = sample(ip[-i1],1,prob=ip[-i1])
                s1 = solutions[[i1]]
                s2 = solutions[[i2]]
                new_solutions[[i]] %<-% full_cross_over(s1,s2,fimerge,fiswap,pmut) %globals% c("s1","s2","fimerge","fiswap","pmut","full_cross_over")
              }
              solutions = c(bres,as.list(new_solutions))
              icls = sapply(solutions,function(s){s@icl})
              solutions=solutions[!is.nan(icls)]
              icls=icls[!is.nan(icls)]
              print(icls)
              old_best=best_icl
              best_icl = max(icls)
              nbgen = nbgen + 1;
              
              #print("#################")
              #print(paste0("Generation N",nbgen, " best solution with an ICL of ",round(solutions[[which.max(icls)]]@icl)," and ",solutions[[which.max(icls)]]@K," clusters."))
              #print("#################")
            }
            

            train.hist=rbind(train.hist,data.frame(generation=nbgen,icl=icls,K=sapply(solutions,function(s){max(s@cl)})))
            #parallel::stopCluster(cl)
            # best solution
            res = solutions[[order(icls,decreasing = TRUE)[1]]]
            # compute merge path
            path = fit_greed_path(data,res)
            # clean the resuts (compute, merge tree,...)
            path = greed:::cleanpath(path)
            # store train history
            path@train_hist = train.hist
            # stop future plan
            #oplan <- future::plan()
            #on.exit(future::plan(oplan), add = TRUE)
            path
          }


full_cross_over = function(sol1,sol2,fimerge,fiswap,pmutation){
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
  if(nrow(ijAm)>0){
    Am=Matrix::sparseMatrix(ijAm[,1],ijAm[,2],x=rep(1,nrow(ijAm)),dims = c(max(ncl),max(ncl)))
    sol=fimerge(ncl,Am)
  }else{
    sol=sol1;
  }
  ncl = sol@cl
  iclust = 1:max(ncl)
  if(runif(1)<pmutation){
    
    sp_cl=sample(max(ncl),1)
    ws = as.numeric(ncl==sp_cl)

    ncl[ncl==sp_cl]=sample(c(sp_cl,max(ncl)+1),sum(ncl==sp_cl),replace=TRUE)
    if(max(ncl)>10){
      iclust  = unique(c(sp_cl,max(ncl),sample(max(ncl),8)))  
    }else{
      iclust = 1:max(ncl)
    }
    
    sol= fiswap(ncl,as.numeric(ncl%in% iclust),iclust)
  }

  sol
}





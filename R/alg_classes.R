#' @useDynLib greed
#' @importFrom Rcpp sourceCpp
#' @importFrom future %<-%
#' @name %<-%
NULL

#' @include models_classes.R fit_classes.R
#' @import Matrix
NULL

#' @title Optimization algorithm classes
#' 
#' @name algs-classes 
NULL



#' @rdname algs-classes
#' @title alg 
#' An S4 class to represent an abstract optimisation algorithm.
#' @slot name Name of the algorithm
#' @export
setClass("alg",slots = list(name = "character"))


#' @rdname algs-classes 
#' @title gree 
#' An S4 class to represent a greedy algorithm extends \code{alg} class with multiple start.
#' @slot nb_start number of random starts (default to 10)
#' @export
setClass("greed",
         contains = "alg",
         representation =  list(nb_start="numeric"),
         prototype(name="greed",nb_start=10))


#' @rdname algs-classes 
#' @title km 
#' An S4 class to represent a greedy algorithm extends \code{alg} class with initialization with spectral clustering and or k-means.
#' @export
setClass("seed",
         contains = "alg",
         representation =  list(),
         prototype(name="seed"))


#' @rdname algs-classes
#' @title genetic
#' An S4 class to represent a hybrid genetic/greedy algorithm extends \code{alg} class.
#' @slot pop_size size of the solutions populations (default to 10)
#' @slot nb_max_gen maximal number of generation to produce (default to 4) 
#' @export
setClass("hybrid",
         contains = "alg",
         representation =  list(pop_size = "numeric",nb_max_gen = "numeric"),
         prototype(name="hybrid",pop_size=10, nb_max_gen = 10))


#' @rdname algs-classes
#' @title genetic
#' An S4 class to represent a hybrid genetic/greedy algorithm extends \code{alg} class.
#' @slot pop_size size of the solutions populations (default to 10)
#' @slot nb_max_gen maximal number of generation to produce (default to 4) 
#' @export
setClass("genetic",
         contains = "alg",
         representation =  list(pop_size = "numeric",nb_max_gen = "numeric",prob_mut="numeric",sel_frac="numeric"),
         prototype(name="geno",pop_size=100, nb_max_gen = 20,prob_mut=0.1,sel_frac=0.75))


#' @describeIn fit 
#' @title Fit a clustering model
#' 
#' @param x A sparse Matrix as a \code{dgCMatrix}
#' @param K An initial guess of the maximal number of cluster
#' @param model An \code{\link{IclModel-class}} such as \code{\link{sbm-class}}, \code{\link{dcsbm-class}}, ...
#' @param alg An optimization algorithm such as \code{\link{greed-class}}, \code{\link{genetic-class}} or \code{\link{km-class}}
#' @export
setGeneric("fit", function(x,K,model,alg,...) standardGeneric("fit")) 


#' @describeIn fit 
#' @export
setMethod(f = "fit", 
          signature = signature("dgCMatrix", "numeric","missing","missing"),
          definition = function(x,K,verbose=FALSE){
            # only a sparseMatrix and a number check dim to choose between graph models and mm
            if(dim(x)[1]==dim(x)[2]){
              fit(x,K,new("sbm"),new("hybrid"),verbose=verbose)  
            }else{
              fit(x,K,new("mm"),new("hybrid"),verbose=verbose)
            }
          });

#' @describeIn fit 
#' @export
setMethod(f = "fit", 
          signature = signature("dgCMatrix", "numeric","icl_model","missing"),
          definition = function(x,K,model,verbose=FALSE){
            fit(x,K,model,new("hybrid"),verbose=verbose)
          });

#' @describeIn fit 
#' @export
setMethod(f = "fit", 
          signature = signature("dgCMatrix", "numeric","icl_model","hybrid"), 
          definition = function(x, K,model,alg,verbose=FALSE){

            
            train.hist = data.frame(generation=c(),icl=c(),K=c())
            
            # multi-start in //
            #future::plan(future::multiprocess)
            
            solutions = listenv::listenv()
            # first generation of solutions
            pop_size = alg@pop_size
            for (i in 1:pop_size){
              solutions[[i]] %<-% fit_greed(model,x,K,verbose=verbose)
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

            while((max(icls)-min(icls))>1 & (best_icl > old_best) & nbgen < alg@nb_max_gen){
              
              
              train.hist=rbind(train.hist,data.frame(generation=nbgen,icl=icls,K=sapply(solutions,function(s){max(s@cl)})))
              # selection keep the top half solutions
              icl_order = order(icls,decreasing = TRUE)
              selected  = icl_order[1:(pop_size/2)]
              # cross_over between the kept solution
              new_solutions = listenv::listenv()
              selected_couples = matrix(selected[sample(1:length(selected),length(selected)*2,replace = TRUE)],ncol=2)
              for (i in 1:nrow(selected_couples)){
                new_solutions[[i]] %<-% full_cross_over(solutions[[selected_couples[i,1]]],solutions[[selected_couples[i,2]]],model,x,verbose)
              }
              new_solutions = as.list(new_solutions)
              solutions = c(solutions[selected],new_solutions)
              icls = sapply(solutions,function(s){s@icl})
              old_best=best_icl
              best_icl = max(icls)
              nbgen = nbgen + 1;
            }
            train.hist=rbind(train.hist,data.frame(generation=nbgen,icl=icls,K=sapply(solutions,function(s){max(s@cl)})))
            #parallel::stopCluster(cl)
            # best solution
            res = solutions[[order(icls,decreasing = TRUE)[1]]]
            # compute merge path
            path = fit_greed_path(x,res)
            # clean the resuts (compute, merge tree,...)
            path = cleanpath(path)
            # store train history
            path@train_hist = train.hist
            # stop future plan
            #oplan <- future::plan()
            #on.exit(future::plan(oplan), add = TRUE)
            path
          })


full_cross_over = function(sol1,sol2,model,x,verbose){
  # cartesian product on the z of the two solution
  ncl = unclass(factor(paste(sol1@cl,sol2@cl)))
  # greedy merge
  fit_greed_init(model,x,ncl,"merge",verbose=verbose)
}


#' @describeIn fit 
#' @export
setMethod(f = "fit", 
          signature = signature("dgCMatrix", "numeric","icl_model","genetic"), 
          definition = function(x, K,model,alg,verbose=FALSE){
            
            train.hist = data.frame(generation=c(),icl=c(),K=c())
            
            # multi-start in //
            #future::plan(future::multiprocess)
            
            solutions = listenv::listenv()
            # first generation of solutions
            pop_size = alg@pop_size
            for (i in 1:pop_size){
              Kc = sample(2:K,1)
              cl = sample(Kc,nrow(x),replace = TRUE)
              solutions[[i]] %<-% init(model,x,cl)
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
                children[[i]] = cross_over(solutions[[i1]],solutions[[i2]],model,x,verbose)
              }
              children = as.list(children)
              for (i in 1:(alg@pop_size-1)){
                if(runif(1)<alg@prob_mut){
                  new_solutions[[i]] %<-% fit_greed_init(model,x,children[[i]]@cl,"swap",1,verbose)
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

            res = fit_greed_init(model,x,res@cl,"both",verbose=verbose)
            path = fit_greed_path(x,res)
            # clean the resuts (compute, merge tree,...)

            path = cleanpath(path)
            # store train history
            path@train_hist = train.hist
            path
          })

cross_over = function(sol1,sol2,model,x,verbose){
  cl1 = sol1@cl
  K1 = sol1@K
  lk1 = 1:K1
  cl2 = sol2@cl
  K2=sol2@K
  lk2=1:K2
  C=rep(0,length(cl1))
  K=0
  while(sum(C==0)>0){
    if(K1>0  && runif(1)>0.5){
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
  init(model,x,C)
}

#' @describeIn fit 
#' @export
setMethod(f = "fit", 
          signature = signature("dgCMatrix", "numeric","icl_model","greed"), 
          definition = function(x, K,model,alg,verbose=FALSE){


            solutions = listenv::listenv()
            for (i in 1:alg@nb_start){
              solutions[[i]] %<-% fit_greed(model,x,K,verbose=verbose)
            }
            solutions = as.list(solutions)
            icls = sapply(solutions,function(s){s@icl})

            res = solutions[[order(icls,decreasing = TRUE)[1]]]
            path = fit_greed_path(x,res)
            cleanpath(path)
            path@train_hist = data.frame(icl=icls,K= sapply(solutions,function(s){max(s@cl)}))

            path
          })

#' @describeIn fit 
#' @export
setMethod(f = "fit", 
          signature = signature("dgCMatrix", "numeric","icl_model","seed"), 
          definition = function(x, K,model,alg,verbose=FALSE){
            if(class(model)=="dcsbm" | class(model)=="sbm"){
              cl = spectral(x,K)  
            }
            if(class(model)=="mm"){
              cl = kmeans(x,K)  
            }
            res = fit_greed_init(model,x,cl,"both",verbose=verbose)
            path = fit_greed_path(x,res)
            p=cleanpath(path)   
            

          })


# clean the merge path 
cleanpath = function(pathsol){
  K=pathsol@K
  pathsol@logalpha = 0
  path=pathsol@path
  tree =c(0)
  xtree=c(0)
  cn  = 1
  lab = c(1)
  xpos = c(0)
  H = rep(0,2*K-1)
  x = 0
  w = 0.5
  K=1
  for (lev in seq(length(path),1)){
    pl = length(path)-lev
    ord = order(xpos)
    path[[lev]]$obs_stats = reorder(pathsol@model,path[[lev]]$obs_stats,ord)
    path[[lev]]$cl  = order(ord)[path[[lev]]$cl] 
    k=path[[lev]]$k
    l=path[[lev]]$l
    tree=c(tree,lab[l],lab[l])
    
    H[lab[l]]=-path[[lev]]$logalpha
    if(tree[lab[l]]!=0 && H[lab[l]]>H[tree[lab[l]]]){
      H[lab[l]]=H[tree[lab[l]]]
    }

    lab[l]=cn+1
    fpos = xpos[l]
    xtree=c(xtree,fpos-w^pl,fpos+w^pl)
    # choisir + ou - en fonctionde la taille ?
    xpos[l]=fpos-w^pl
    if(k>K){
      xpos = c(xpos,fpos+w^pl)
      lab=c(lab,cn+2)  
    }else{
      xpos = c(xpos[1:(k-1)],fpos+w^pl,xpos[k:length(lab)])
      lab=c(lab[1:(k-1)],cn+2,lab[k:length(lab)])
    }
    K  = K+1
    cn = cn + 2 
  }
  ggtree=data.frame(H=H,tree=tree,x=xtree,node=1:length(tree),xmin=0,xmax=0)
  # recalculer les x
  leafs = which(ggtree$H==0)
  or = order(ggtree[leafs,"x"])
  ggtree$x[leafs[or]]=seq(-1,1,length.out = length(leafs))
  others = ggtree$node[ggtree$H!=0]
  for(n in others[seq(length(others),1)]){
    sons=which(ggtree$tree==n)
    ggtree$x[n]=mean(ggtree$x[sons])
    ggtree$xmin[n] = min(ggtree$x[sons])
    ggtree$xmax[n] = max(ggtree$x[sons])
  }



  ggtree$Hend = c(-1,ggtree$H[ggtree$tree])
  ggtree$xend = c(-1,ggtree$x[ggtree$tree])

  or = order(xpos)
  pathsol@obs_stats = reorder(pathsol@model,pathsol@obs_stats,or)
  pathsol@cl=order(or)[pathsol@cl] 
  pathsol@path = path
  pathsol@tree = tree
  pathsol@ggtree = ggtree 
  pathsol
} 

reorder_sbm = function(obs_stats,or){
  obs_stats$counts = obs_stats$counts[or]
  obs_stats$x_counts = obs_stats$x_counts[or,or]
  obs_stats
}

reorder_dcsbm = function(obs_stats,or){
  obs_stats$counts = obs_stats$counts[or]
  obs_stats$din = obs_stats$din[or]
  obs_stats$dout = obs_stats$dout[or]
  obs_stats$x_counts = obs_stats$x_counts[or,or]
  obs_stats
}

reorder_mm = function(obs_stats,or){
  obs_stats$counts = obs_stats$counts[or]
  obs_stats$x_counts = obs_stats$x_counts[or,]
  obs_stats
}

setGeneric("reorder", function(model, obs_stats,order) standardGeneric("reorder")) 

setMethod(f = "reorder", 
          signature = signature("sbm", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_sbm(obs_stats,order)
          })
setMethod(f = "reorder", 
          signature = signature("mm", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_mm(obs_stats,order)
          })
setMethod(f = "reorder", 
          signature = signature("dcsbm", "list","integer"), 
          definition = function(model, obs_stats,order){
            reorder_dcsbm(obs_stats,order)
          })

#' @describeIn cut
#' @title Cut a path to a desired number of cluster 
#' 
#' @param x A an \code{icl_path} solution 
#' @param K Desired number of cluster
#' @describeIn cut
#' @export
setMethod(f = "cut", 
          signature = signature("icl_path"), 
          definition = function(x, K){
            i = which(sapply(x@path,function(p){p$K})==K)
            x@tree = x@tree[1:(2*K-1)]
            x@ggtree = x@ggtree[1:(2*K-1),]
            x@K = K
            x@logalpha=x@path[[i]]$logalpha
            x@icl = x@path[[i]]$icl
            x@cl = as.vector(x@path[[i]]$cl)
            for(st in names(x@obs_stats)){
              x@obs_stats[st] = x@path[[i]]$obs_stats[st]
            }

            x@path=x@path[(i+1):length(x@path)]
            x
            
})

#' @title spectral
#' Regularized spectral clustering nips paper 2013
#' @param x An adjacency matrix in sparse format
#' @param K Desired number of cluster
#' @value cl Vector of clsuter labels
#' @export
spectral= function(X,K){
  X = X+Matrix::t(X)
  ij=Matrix::which(X>0,arr.ind = T)
  X[ij]=1
  nu = sum(X)/dim(X)[1]
  D  = (Matrix::rowSums(X)+nu)^-0.5
  Ln = Matrix::sparseMatrix(ij[,1],ij[,2],x=D[ij[,1]]*D[ij[,2]])
  V = rARPACK::eigs(Ln,K)
  Xp = V$vectors
  Xpn = Xp/sqrt(rowSums(Xp)^2)
  Xpn[rowSums(Xp)==0,]=0
  km = stats::kmeans(Xp,K)
  km$cluster
}



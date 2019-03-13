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
#' @title greedy 
#' An S4 class to represent a greedy algorithm extends \code{alg} class.
#' @slot name greed
#' @slot nb_start number of random starts (default to 10)
#' @export
setClass("greed",
         contains = "alg",
         representation =  list(nb_start="numeric"),
         prototype(name="greed",nb_start=10))


#' @rdname algs-classes 
#' @title greedy 
#' An S4 class to represent a greedy algorithm extends \code{alg} class.
#' @slot name greed
#' @slot nb_start number of random starts (default to 10)
#' @export
setClass("km",
         contains = "alg",
         representation =  list(),
         prototype(name="km"))


#' @rdname algs-classes
#' @title genetic
#' An S4 class to represent a hybrid genetic/greedy algorithm extends \code{alg} class.
#' @slot name genetic
#' @slot pop_size size of the solutions populations (default to 10)
#' @slot nb_max_gen maximal number of generation to produce (default to 4) 
#' @export
setClass("genetic",
         contains = "alg",
         representation =  list(pop_size = "numeric",nb_max_gen = "numeric"),
         prototype(name="genetic",pop_size=10, nb_max_gen = 10))


#' @describeIn fit 
#' @title Fit a clustering model
#' 
#' @param x A sparse Matrix dgCMatrix
#' @param K An initial guess of the maximal number of cluster
#' @param model An IclModel
#' @param alg An optimization algorithm
#' @export
setGeneric("fit", function(x,K,model,alg) standardGeneric("fit")) 

#' @describeIn fit 
#' @title Fit a clustering model
#' 
#' @param x A sparse Matrix dgCMatrix
#' @param K An initial guess of the maximal number of cluster
setMethod(f = "fit", 
          signature = signature("dgCMatrix", "numeric","missing","missing"),
          definition = function(x,K){
            # only a sparseMatrix and a number check dim to choose between graph models and 
            if(dim(x)[1]==dim(x)[2]){
              fit(x,K,new("sbm"),new("genetic"))  
            }else{
              fit(x,K,new("mm"),new("genetic"))
            }
          });

#' @describeIn fit 
#' @title Fit a clustering model
#' 
#' @param x A sparse Matrix dgCMatrix
#' @param K An initial guess of the maximal number of cluster
#' @param model An IclModel
#' @export
setMethod(f = "fit", 
          signature = signature("dgCMatrix", "numeric","icl_model","missing"),
          definition = function(x,K,model){
            fit(x,K,model,new("genetic"))
          });

#' @describeIn fit 
#' @title Fit a clustering model
#' 
#' @param x A sparse Matrix dgCMatrix
#' @param K An initial guess of the maximal number of cluster
#' @param model An IclModel
#' @param alg An optimization algorithm
setMethod(f = "fit", 
          signature = signature("dgCMatrix", "numeric","icl_model","genetic"), 
          definition = function(x, K,model,alg){
            #cl <- parallel::makeCluster(K, timeout = 60)
            #future::plan(future::cluster,workers = cl)
            future::plan(future::multiprocess)
            solutions = listenv::listenv()
            # premiere generation
            pop_size = alg@pop_size
            for (i in 1:pop_size){
     
              solutions[[i]] %<-% fit_greed(model,x,K)
            }
   
            solutions = as.list(solutions)
            icls  = sapply(solutions,function(s){s@icl})
            solutions=solutions[!is.nan(icls)]
            icls=icls[!is.nan(icls)]
            print(icls)
            nbgen = 1
            # tout le monde a converger vers la meme solution
            while((max(icls)-min(icls))>1 & nbgen < alg@nb_max_gen){
              # selections 
              print(paste0("GEN: ",nbgen ))
              icl_order = order(icls,decreasing = TRUE)
              selected  = icl_order[1:(pop_size/2)]
              # cross_over
              new_solutions = listenv::listenv()
              selected_couples = matrix(selected[sample(1:length(selected),length(selected)*2,replace = TRUE)],ncol=2)
              
              for (i in 1:nrow(selected_couples)){
                
                new_solutions[[i]] %<-% cross_over(solutions[[selected_couples[i,1]]],solutions[[selected_couples[i,2]]],model,x)
              }
              new_solutions = as.list(new_solutions)
              solutions = c(solutions[selected],new_solutions)
              icls = sapply(solutions,function(s){s@icl})
              print(icls)
              nbgen = nbgen + 1;
            }
            
            #parallel::stopCluster(cl)
            
            res = solutions[[order(icls,decreasing = TRUE)[1]]]
            pathsol = fit_greed_path(x,res)
            cleanpath(pathsol)   
          })


cross_over = function(sol1,sol2,model,x){
  ncl = unclass(factor(paste(sol1@cl,sol2@cl)))
  fit_greed_init(model,x,ncl,"merge")
}

#' @describeIn fit 
#' @title Fit a clustering model
#' 
#' @param x A sparse Matrix dgCMatrix
#' @param K An initial guess of the maximal number of cluster
#' @param model An IclModel
#' @param alg An optimization algorithm
#' @export
setMethod(f = "fit", 
          signature = signature("dgCMatrix", "numeric","icl_model","greed"), 
          definition = function(x, K,model,alg){
            #cl <- parallel::makeCluster(K, timeout = 60)
            #future::plan(future::cluster,workers = cl)
            future::plan(future::multiprocess)
            solutions = listenv::listenv()
            for (i in 1:alg@nb_start){
              solutions[[i]] %<-% fit_greed(model,x,K)
            }
            solutions = as.list(solutions)
            icls = sapply(solutions,function(s){s@icl})
            #parallel::stopCluster(cl)
            res = solutions[[order(icls,decreasing = TRUE)[1]]]
            pathsol = fit_greed_path(x,res)
            cleanpath(pathsol)   
          })

#' @describeIn fit 
#' @title Fit a clustering model
#' 
#' @param x A sparse Matrix dgCMatrix
#' @param K An initial guess of the maximal number of cluster
#' @param model An IclModel
#' @param alg An optimization algorithm
#' @export
setMethod(f = "fit", 
          signature = signature("dgCMatrix", "numeric","icl_model","km"), 
          definition = function(x, K,model,alg){
            if(class(model)=="dcsbm" | class(model)=="sbm"){
              cl = spectral(x,K)  
            }
            if(class(model)=="mm"){
              cl = kmeans(x,K)  
            }
            res = fit_greed_init(model,x,cl,"both")
            pathsol = fit_greed_path(x,res)
            p=cleanpath(pathsol)   
            

          })


cleanpath = function(pathsol){
  K=length(pathsol@obs_stats$counts)
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
    #path[[lev]]$lab  = lab
    #path[[lev]]$xpos = xpos
    #path[[lev]]$perm = order(xpos)
    path[[lev]] = reorder(pathsol@model,path[[lev]],order(xpos))
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
    print(n)
    print(ggtree$H[n])
    print(sons)
    print(ggtree$x[sons])
    ggtree$x[n]=mean(ggtree$x[sons])
    ggtree$xmin[n] = min(ggtree$x[sons])
    ggtree$xmax[n] = max(ggtree$x[sons])
  }



  ggtree$Hend = c(-1,ggtree$H[ggtree$tree])
  ggtree$xend = c(-1,ggtree$x[ggtree$tree])
  print("ordering...")
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
  if(!is.null(obs_stats$cl)){
    obs_stats$cl = order(or)[obs_stats$cl] 
  }
  obs_stats
}

reorder_dcsbm = function(obs_stats,or){
  obs_stats$counts = obs_stats$counts[or]
  obs_stats$din = obs_stats$din[or]
  obs_stats$dout = obs_stats$dout[or]
  obs_stats$x_counts = obs_stats$x_counts[or,or]
  if(!is.null(obs_stats$cl)){
    obs_stats$cl = order(or)[obs_stats$cl] 
  }
  obs_stats
}

reorder_mm = function(obs_stats,or){
  obs_stats$counts = obs_stats$counts[or]
  obs_stats$x_counts = obs_stats$x_counts[or,]
  if(!is.null(obs_stats$cl)){
    obs_stats$cl = order(or)[obs_stats$cl] 
  }
  obs_stats
}

setGeneric("reorder", function(model, obs_stats,order) standardGeneric("reorder")) 

setMethod(f = "reorder", 
          signature = signature("sbm", "list","numeric"), 
          definition = function(model, obs_stats,order){
            reorder_sbm(obs_stats,order)
          })
setMethod(f = "reorder", 
          signature = signature("mm", "list","numeric"), 
          definition = function(model, obs_stats,order){
            reorder_mm(obs_stats,order)
          })
setMethod(f = "reorder", 
          signature = signature("dcsbm", "list","numeric"), 
          definition = function(model, obs_stats,order){
            reorder_dcsbm(obs_stats,order)
          })


setGeneric("cut", function(fit,K,...) standardGeneric("cut")) 
setMethod(f = "cut", 
          signature = signature("icl_path", "numeric"), 
          definition = function(fit, K){
            i = which(sapply(fit@path,function(p){p$K})==K)
            print(i)
            fit@tree = fit@tree[1:(2*K-1)]
            fit@ggtree = fit@ggtree[1:(2*K-1),]
            fit@K = K
            fit@logalpha=fit@path[[i]]$logalpha
            fit@icl = fit@path[[i]]$icl
            fit@cl = fit@path[[i]]$cl
            for(st in names(fit@obs_stats)){
              fit@obs_stats[st] = fit@path[[i]][st]
            }

            fit@path=fit@path[(i+1):length(fit@path)]
            fit
            
})

spectral= function(X,K){
  X = X+t(X)
  X[X>0]=1
  nu=sum(X)/dim(X)[1]
  L = diag((rowSums(X)+nu)^-0.5)%*%X%*%diag((colSums(X)+nu)^-0.5)
  S = svd(L,K,K)
  Xp=S$u
  #Xp = S$u/rowSums(S$u)
  #Xp[rowSums(S$u)==0,]=0
  km = kmeans(Xp,K)
  km$cluster
}



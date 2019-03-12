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
            # première generation
            pop_size = alg@pop_size
            for (i in 1:pop_size){
              solutions[[i]] %<-% fit_greed(model,x,K)
            }
            solutions = as.list(solutions)
            icls  = sapply(solutions,function(s){s@icl})
            print(icls)
            nbgen = 1
            # tout le monde a converger vers la même solution
            while((max(icls)-min(icls))>1 & nbgen < alg@nb_max_gen){
              # sélections 
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
  fit_greed_init(model,x,max(ncl),ncl-1)
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
            cl <- parallel::makeCluster(K, timeout = 60)
            future::plan(future::cluster,workers = cl)
            #future::plan(future::multiprocess)
            solutions = listenv::listenv()
            for (i in 1:alg@nb_start){
              solutions[[i]] %<-% fit_greed(model,x,K)
            }
            solutions = as.list(solutions)
            icls = sapply(solutions,function(s){s@icl})
            print(icls)
            parallel::stopCluster(cl)
            solutions[[order(icls,decreasing = TRUE)[1]]]  
          })

cleanpath = function(pathsol){
  K=length(pathsol@counts)
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
    print(xpos)
    pl = length(path)-lev
    path[[lev]]$lab  = lab
    path[[lev]]$xpos = xpos
    path[[lev]]$perm = order(xpos)
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
  others = setdiff(1:nrow(ggtree),leafs)
  for(n in others[order(H[others])]){
    n
    sons=which(ggtree$tree==n)
    ggtree$x[n]=mean(ggtree$x[sons])
    ggtree$xmin[n] = min(ggtree$x[sons])
    ggtree$xmax[n] = max(ggtree$x[sons])
  }



  ggtree$Hend = c(-1,ggtree$H[ggtree$tree])
  ggtree$xend = c(-1,ggtree$x[ggtree$tree])
  #ggplot()+geom_segment(data=ggtree[-1,],aes(x=x,y=H,xend=xend,yend=Hend))+geom_point(data=ggtree,aes(x=x,y=H))
  print("ordering...")
  pathsol@x_counts = pathsol@x_counts[order(xpos),order(xpos)]
  pathsol@counts = pathsol@counts[order(xpos)]  
  # pathsol@cl = pathsol@cl[order()]
  
  pathsol@path = path
  pathsol@tree = tree
  pathsol@ggtree = ggtree 
  pathsol
} 



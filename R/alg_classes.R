#' @useDynLib greed
#' @importFrom Rcpp sourceCpp
#' @importFrom future %<-%
#' @name %<-%
NULL

#' @include models_classes.R fit_classes.R hybrid_alg.R genetic_alg.R misc.R
#' @import Matrix
NULL


#' @title Abstract optimization algorithm class
#' 
#' @description
#' An S4 class to represent an abstract optimization algorithm.
#' @slot name algorithm name
#' @export
setClass("alg",slots = list(name = "character"))



#' @title Greedy algorithm with multiple start class
#' 
#' @description 
#' An S4 class to represent a greedy algorithm  with multiple start (extends \code{\link{alg-class}} class).
#' @slot nb_start number of random starts (default to 10)
#' @export
setClass("multistarts",
         contains = "alg",
         representation =  list(nb_start="numeric"),
         prototype(name="greed",nb_start=10))



#' @title Greedy algorithm with seeded initialization
#' 
#' @description
#' An S4 class to represent a greedy algorithm with initialization from spectral clustering and or k-means (extends \code{\link{alg-class}} class ).
#' @export
setClass("seed",
         contains = "alg",
         representation =  list(),
         prototype(name="seed"))



#' @title Hybrid optimization algorithm 
#' 
#' @description 
#' An S4 class to represent an hybrid genetic/greedy algorithm (extends \code{\link{alg-class}} class).
#' @slot pop_size size of the solutions populations (default to 20)
#' @slot nb_max_gen maximal number of generation to produce (default to 10)
#' @slot prob_mutation mutation probability (default to 0.25)
#' @slot Kmax maximum number of clusters (default to 100) 
#' @export
setClass("hybrid",
         contains = "alg",
         representation =  list(pop_size = "numeric",nb_max_gen = "numeric",prob_mutation = "numeric",Kmax="numeric"),
         prototype(name="hybrid",pop_size=20, nb_max_gen = 10,prob_mutation=0.25,Kmax=100))



#' @title Genetic optimization algorithm
#' 
#' @description
#' An S4 class to represent a genetic algorithm (extends \code{\link{alg-class}} class).
#' @slot pop_size size of the solutions populations (default to 10)
#' @slot nb_max_gen maximal number of generation to produce (default to 4) 
#' @slot prob_mutation probability of mutation (default to 0.25)
#' @slot sel_frac fraction of best solutions selected for crossing  (default to 0.75)
#' @export
setClass("genetic",
         contains = "alg",
         representation =  list(pop_size = "numeric",nb_max_gen = "numeric",prob_mutation="numeric",sel_frac="numeric"),
         prototype(name="genetic",pop_size=100, nb_max_gen = 20,prob_mutation=0.25,sel_frac=0.75))





#' @title Method to cut a path solution to a desired number of cluster 
#' 
#' @description This method take a \code{\link{icl_path-class}} object and an integer K and return the solution from the path with K clusters 
#' @param x A an \code{icl_path} solution 
#' @param K Desired number of cluster
#' @return an \code{\link{icl_path-class}} object with the desired number of cluster
#' @export
setMethod(f = "cut", 
          signature = signature("icl_path"), 
          definition = function(x, K){
            if(K<x@K){
              i = which(sapply(x@path,function(p){p$K})==K)
              x@K = K
              x@logalpha=x@path[[i]]$logalpha
              x@icl = x@path[[i]]$icl
              x@cl = as.vector(x@path[[i]]$cl)
              for(st in names(x@obs_stats)){
                x@obs_stats[st] = x@path[[i]]$obs_stats[st]
              }
              
              x@path=x@path[(i+1):length(x@path)]
              x
            }else{
              warning(paste0("This clustering has ",x@K," clusters and you requested ",K ," clusters. Please provide a value for K smaller than ",x@K,"."),call. = FALSE)
            }
            x
            
})

#' @title Model based hierarchical clustering 
#' 
#' @description
#' 
#' @param X data to cluster either a matrix,an array or a \code{\link{dgCMatrix-class}}
#' @param K initial number of cluster
#' @param model a generative model to fit \code{\link{sbm-class}}, \code{\link{dcsbm-class}}, \code{\link{co_dcsbm-class}}, \code{\link{mm-class}},\code{\link{gmm-class}}, \code{\link{diaggmm-class}} or \code{\link{mvmreg-class}}
#' @param alg an optimization algorithm of class \code{\link{hybrid-class}} (default), \code{\link{multistarts-class}}, \code{\link{seed-class}} or \code{\link{genetic-class}}
#' @param verbose Boolean for verbose mode 
#' @return an \code{\link{icl_path-class}} object
#' @export
greed = function(X,K=20,model=find_model(X),alg=methods::new("hybrid"),verbose=FALSE){
  data = preprocess(model,X)
  modelname = toupper(model@name)
  if("type" %in% methods::slotNames(model)){
    modelname = paste(model@type,modelname)
  }
  if("sampling" %in% methods::slotNames(model)){
    modelname = paste(modelname,"with",model@sampling,"sampling")
  }
  cat(paste0("------- ",modelname, " model fitting ------\n"))
  sol = fit(model,alg,data,K,verbose)
  sol = postprocess(sol,data)
  cat("------- Final clustering -------\n")
  print(sol)
  cat("\n")
  sol
}


find_model = function(X){
  if(methods::is(X,"array") && length(dim(X))>2){
    dimensions = dim(X)
    if(dimensions[1]!=dimensions[2]){
      stop(paste0("Multinomial SBM expect a cube with as many rows as columns:",dimensions,collapse = " x "))
    }
    if(length(dimensions)!=3){
      stop(paste0("Multinomial SBM expect a cube found an array od dimensions:",dimensions,collapse = " x "))
    }
    issym = all(sapply(1:dim(X)[3],function(d){ isSymmetric(X[,,d])}))
    if(issym){
      model = methods::new("multsbm",type="undirected")
    }else{
      model = methods::new("multsbm")
    }
    
  }else{
    if(methods::is(X,"data.frame")){
      X=as.matrix(X)
    }
    if(methods::is(X,"dgCMatrix") | methods::is(X,"matrix")){
      if(nrow(X)==ncol(X)){
        if(sum(is.na(X))>0){
          if(isSymmetric(X)){
            model = methods::new("misssbm",type="undirected")
          }else{
            model = methods::new("misssbm")
          } 
        }else{
          if(isSymmetric(X)){
            model = methods::new("dcsbm",type="undirected")
          }else{
            model = methods::new("dcsbm")
          } 
        }
      }else{
        if(all(round(X)==X)){
          model = methods::new("co_dcsbm")  
        }else{
          model = methods::new("gmm",N0=ncol(X),epsilon=0.1*diag(diag(stats::cov(X))),mu=apply(X,2,mean),tau=0.01)
          #model = methods::new("diaggmm",mu=apply(X,2,mean),beta=0.1)
        }
      }
    }else{
      stop(paste0("Unsupported data type: ", class(X) ," use a data.frame, a matrix, a sparse dgCMatrix or an array."),call. = FALSE)
    }
  }
  model
}

#' @title Conditional model based hierarchical clustering
#' 
#' @param X design matrix
#' @param Y target variables
#' @param K Desired number of cluster
#' @param model a conditional generative model \code{\link{mvmreg-class}}
#' @param alg an optimization algorithm of class \code{\link{hybrid-class}} (default), \code{\link{multistarts-class}}, \code{\link{seed-class}} or \code{\link{genetic-class}}
#' @param verbose Boolean for verbose mode 
#' @return an \code{\link{icl_path-class}} object
#' @export
greed_cond = function(X,Y,K=20,model=find_model_cond(X,Y),alg=methods::new("hybrid"),verbose=FALSE){
  data = preprocess(model,list(X=X,Y=Y))
  modelname = toupper(model@name)
  if("type" %in% methods::slotNames(model)){
    modelname = paste(model@type,modelname)
  }
  if("sampling" %in% methods::slotNames(model)){
    modelname = paste(modelname,"with",model@sampling,"sampling")
  }
  cat(paste0("------- ",modelname, " model fitting ------\n"))
  sol = fit(model,alg,data,K,verbose)
  sol = postprocess(sol,data)
  cat("------- Final clustering -------\n")
  print(sol)
  cat("\n")
  sol
}


find_model_cond = function(X,Y){
  methods::new("mvmreg",N0=ncol(X)+1)
}
setGeneric("reorder", function(model, obs_stats,order) standardGeneric("reorder")) 

setGeneric("fit",function(model,alg,...) standardGeneric("fit")) 


setMethod(f = "fit", 
          signature = signature("icl_model","hybrid"), 
          definition = function(model,alg,data, K,verbose=FALSE){
            hybrid(model,alg,data,K,verbose)
          })

setMethod(f = "fit", 
          signature = signature("icl_model","genetic"), 
          definition = function(model,alg,data,k,verbose=FALSE){
            genetic(model,alg,data,k,verbose=verbose)
          })


setMethod(f = "fit", 
          signature = signature("icl_model","multistarts"), 
          definition = function(model,alg,data, K=20,verbose=FALSE){
            multistart(model,alg,data,K,verbose)
          })

setGeneric("seed", function(model, data,K) standardGeneric("seed")) 

setMethod(f = "fit", 
          signature = signature("icl_model","seed"), 
          definition = function(model,alg,data, K=20,verbose=FALSE){
            cl = seed(model,data,K)  
            res = fit_greed(model,data,cl,"both",verbose=verbose)
            path = fit_greed_path(data,res)
            p=cleanpathopt(path)   
          })

setGeneric("preprocess", function(model, ...) standardGeneric("preprocess")) 

setMethod(f = "preprocess", 
          signature = signature("icl_model"), 
          definition = function(model, data){
            list(X=as.sparse(data),N=nrow(data))
})

setGeneric("postprocess", function(path, ...) standardGeneric("postprocess")) 

setMethod(f = "postprocess", 
          signature = signature("icl_path"), 
          definition = function(path,data=NULL){
            path    
})


setGeneric("sample_cl", function(model, data,K) standardGeneric("sample_cl")) 

setMethod(f = "sample_cl", 
          signature = signature("icl_model","list","numeric"), 
          definition = function(model,data,K){
            sample(1:K,data$N,replace = TRUE)
          })






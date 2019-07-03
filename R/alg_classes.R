#' @useDynLib greed
#' @importFrom Rcpp sourceCpp
#' @importFrom future %<-%
#' @name %<-%
NULL

#' @include models_classes.R fit_classes.R hybrid_alg.R genetic_alg.R
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
setClass("multistarts",
         contains = "alg",
         representation =  list(nb_start="numeric"),
         prototype(name="greed",nb_start=10))


#' @rdname algs-classes 
#' @title seed 
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
         representation =  list(pop_size = "numeric",nb_max_gen = "numeric",prob_mutation = "numeric",Kmax="numeric"),
         prototype(name="hybrid",pop_size=20, nb_max_gen = 10,prob_mutation=0.25,Kmax=300))


#' @rdname algs-classes
#' @title genetic
#' An S4 class to represent a hybrid genetic/greedy algorithm extends \code{alg} class.
#' @slot pop_size size of the solutions populations (default to 10)
#' @slot nb_max_gen maximal number of generation to produce (default to 4) 
#' @export
setClass("genetic",
         contains = "alg",
         representation =  list(pop_size = "numeric",nb_max_gen = "numeric",prob_mut="numeric",sel_frac="numeric"),
         prototype(name="genetic",pop_size=100, nb_max_gen = 20,prob_mut=0.1,sel_frac=0.75))



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
            p=cleanpath(path)   
})



setGeneric("reorder", function(model, obs_stats,order) standardGeneric("reorder")) 



#' @describeIn cut method to cut a fit to a desired number of cluster 
#' @title Cut
#' cut a path to a desired number of cluster 
#' @param x A an \code{icl_path} solution 
#' @param K Desired number of cluster
#' @return an icl_path obejct with the desired number of cluster
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
#' @param X An adjacency matrix in sparse format
#' @param K Desired number of cluster
#' @return cl Vector of clsuter labels
#' @export
spectral= function(X,K){
  X = X+Matrix::t(X)
  ij=Matrix::which(X>0,arr.ind = T)
  X[ij]=1
  nu = sum(X)/dim(X)[1]
  D  = (Matrix::rowSums(X)+nu)^-0.5
  Ln = Matrix::sparseMatrix(ij[,1],ij[,2],x=D[ij[,1]]*D[ij[,2]])
  V = RSpectra::eigs(Ln,K)
  Xp = V$vectors
  Xpn = Xp/sqrt(rowSums(Xp)^2)
  Xpn[rowSums(Xp)==0,]=0
  km = stats::kmeans(Xp,K)
  km$cluster
}



#' @title greed
#' @param X data to cluster 
#' @param K Desired number of cluster
#' @param model a dcsbm, sbm, mm or mreg model
#' @param alg an optimisation algorithm hybrid, mutlistarts, seed or genetic
#' @param verbose boolean for verbose mode 
#' @return an icl_path object
#' @export
greed = function(X,K=20,model=find_model(X),alg=methods::new("hybrid"),verbose=FALSE){
  data = preprocess(model,X,K)
  cat(paste0("------- Fitting a ",model@name, " model ------\n"))
  sol=fit(model,alg,data,K,verbose)
  sol=postprocess(sol,data)
  sol
}

find_model = function(X){
  if(class(X)=="dgCMatrix" | class(X)=="matrix"){
    if(nrow(X)==ncol(X)){
      model = methods::new("dcsbm")
    }else{
      if(all(round(X)==X)){
        model = methods::new("co_dcsbm")  
      }else{
        model = methods::new("mvmreg")
      }
    }
  }else{
    stop(paste0("Unsupported data type :", class(X) ," use matrix or sparse dgCMatrix."))
  }
  model
}


as.sparse = function(X){
  S = X
  if(class(X)=="matrix"){
    ij= which(X!=0,arr.ind=TRUE)
    S = Matrix::sparseMatrix(ij[,1],ij[,2],x = X[ij]) 
  }
  S
}

setGeneric("preprocess", function(model, ...) standardGeneric("preprocess")) 

#' @title preprocess
#' @param model an icl_model
#' @param X input data to prepare
setMethod(f = "preprocess", 
          signature = signature("icl_model"), 
          definition = function(model, data,K){
            list(X=as.sparse(data),N=nrow(data),moves=as.sparse(matrix(1,K,K)))
})

setGeneric("postprocess", function(path, ...) standardGeneric("postprocess")) 

#' @title postprocess
#' @param path an icl_path object
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


#' @title greed_cond
#' @param X covariable data 
#' @param y target variable
#' @param K Desired number of cluster
#' @param model an mreg model
#' @param alg an optimisation algorithm hybrid, mutlistarts, seed or genetic
#' @param verbose boolean for verbose mode 
#' @return an icl_path object
#' @export
greed_cond = function(X,y,K=20,model=find_model_cond(X,y),alg=methods::new("hybrid"),verbose=FALSE){
  fit(model,alg,list(X=X,y=y,N=dim(X)[1]),K,verbose)
}



find_model_cond = function(X,y){
  methods::new("mreg")
}

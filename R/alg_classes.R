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
         prototype(name="genetic",pop_size=100, nb_max_gen = 20,prob_mut=0.1,sel_frac=0.75))


#' @describeIn fit 
#' @title Fit a clustering model
#' 
#' @param x A sparse Matrix as a \code{dgCMatrix}
#' @param K An initial guess of the maximal number of cluster
#' @param model An \code{\link{IclModel-class}} such as \code{\link{sbm-class}}, \code{\link{dcsbm-class}}, ...
#' @param alg An optimization algorithm such as \code{\link{greed-class}}, \code{\link{genetic-class}} or \code{\link{km-class}}
#' @export
setGeneric("fit",function(model,alg,...) standardGeneric("fit")) 



#' @describeIn fit 
#' @export
setMethod(f = "fit", 
          signature = signature("icl_model","hybrid"), 
          definition = function(model,alg,data, K,verbose=FALSE){
            hybrid(model,alg,data,K,verbose)
})

#' @describeIn fit 
#' @export
setMethod(f = "fit", 
          signature = signature("icl_model","genetic"), 
          definition = function(model,alg,data,k,verbose=FALSE){
            genetic(model,alg,data,k,verbose=verbose)
          })


#' @describeIn fit 
#' @export
setMethod(f = "fit", 
          signature = signature("icl_model","multistarts"), 
          definition = function(model,alg,data, K=20,verbose=FALSE){
            multistart(model,alg,data,K,verbose)
          })

#' @describeIn fit 
#' @export
setMethod(f = "fit", 
          signature = signature("icl_model","seed"), 
          definition = function(model,alg,data, K=20,verbose=FALSE){
            cl = seed(model,data,K)  
            res = fit_greed_init(model,x,cl,"both",verbose=verbose)
            path = fit_greed_path(x,res)
            p=cleanpath(path)   
})



setGeneric("reorder", function(model, obs_stats,order) standardGeneric("reorder")) 



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



#' @title greed
#' @param x An adjacency matrix in sparse format
#' @param K Desired number of cluster
#' @value cl Vector of clsuter labels
#' @export
greed = function(X,K=20,model=find_model(X),alg=new("hybrid"),verbose=FALSE){
  fit(model,alg,list(X=X,N=dim(X)[1]),K,verbose)
}


#' @title greed_cond
#' @param x An adjacency matrix in sparse format
#' @param K Desired number of cluster
#' @value cl Vector of clsuter labels
#' @export
greed_cond = function(X,y,K=20,model=find_model_cond(X,y),alg=new("hybrid"),verbose=FALSE){
  fit(model,alg,list(X=X,y=y,N=dim(X)[1]),K,verbose)
}

find_model = function(X){
  if(class(X)=="dgCMatrix"){
    if(nrow(X)==ncol(X)){
      model = new("dcsbm")
    }else{
      model = new("mm")
    }
  }
  model
}

find_model_cond = function(X,y){
  new("mreg")
}

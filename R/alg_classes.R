#' @useDynLib greed
#' @importFrom Rcpp sourceCpp
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
         prototype(name="genetic",pop_size=10, nb_max_gen = 4))


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
            future::plan(future::multisession)
            solutions = listenv::listenv()
            # première generation
            pop_size = alg@pop_size
            for (i in 1:pop_size){
              s = future::future(fit_greed(model,x,K))
              solutions[[i]] = future::value(s) 
            }
            solutions = as.list(solutions)
            icls  = sapply(solutions,function(s){s@icl})
            print(icls)
            nbgen = 1
            # tout le monde a converger vers la même solution
            while((max(icls)-min(icls))>1 & nbgen < alg@nb_max_gen){
              # sélections 
              print(paste0("GEN :",nbgen ))
              icl_order = order(icls,decreasing = TRUE)
              selected  = icl_order[1:(pop_size/2)]
              # cross_over
              new_solutions = listenv::listenv()
              selected_couples = matrix(selected[sample(1:length(selected),length(selected)*2,replace = TRUE)],ncol=2)
              
              for (i in 1:nrow(selected_couples)){
                
                cvo = future::future(cross_over(solutions[[selected_couples[i,1]]],solutions[[selected_couples[i,2]]],model,x))
                new_solutions[[i]] = future::value(cvo)
              }
              new_solutions = as.list(new_solutions)
              solutions = c(solutions[selected],new_solutions)
              icls = sapply(solutions,function(s){s@icl})
              print(icls)
              nbgen = nbgen + 1;
            }
            
            solutions[[order(icls,decreasing = TRUE)[1]]]
            
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
            future::plan(future::multiprocess)
            solutions = listenv::listenv()
            for (i in 1:alg@nb_start){
              val = future::future(fit_greed(model,x,K))
              solutions[[i]] = future::value(val)  
            }
            solutions = as.list(solutions)
            icls = sapply(solutions,function(s){s@icl})
            print(icls)
            solutions[[order(icls,decreasing = TRUE)[1]]]  
          })

cleanpath = function(path){
  
} 
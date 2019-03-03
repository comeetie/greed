#' @useDynLib gicl
#' @importFrom Rcpp sourceCpp

#' An S4 class to represent a model
#'
#' @slot name a character vector
#' @slot alpha a numeric vector of length 1 which define the parameters of the dirichlet over the cluster proportions (default to 1)
setClass("icl_model",slots = list(name = "character",alpha = "numeric"))

#' An S4 class to represent a stochastick block model that extends \code{icl_model} class.
#'
#' @slot a0 a numeric vector of length 1 which define the parameters of the beta prior over the edges (default to 1)
#' @slot b0 a numeric vector of length 1 which define the parameters of the beta prior over the edges (default to 1)
setClass("sbm",
         representation = list(a0 = "numeric",b0="numeric"),
         contains = "icl_model",
         prototype(name="sbm",a0=1,b0=1,alpha=1))


#' An S4 class to represent a mixture of multinomial also known has mixture of unigrams that extends \code{icl_model} class.
#'
#' @slot a0 a numeric vector of length 1 which define the parameters of the beta prior over the edges (default to 1)
setClass("mm",
         representation = list(beta = "numeric"),
         contains = "icl_model",
         prototype(name="sbm",beta=1,alpha=1))

#' An S4 class to represent a degree corrected stochastick block model that extends \code{icl_model} class.
#'
#' @slot name a character vector
#' @slot p a numeric vector of length 1 which define the parameters of exponential prior over the edges counts (default to 1)
setClass("sbm_dg", representation = list(p="numeric"), contains = "icl_model",prototype(name="sbm_dg",p=1,alpha=1))


#' An S4 class to represent an icl fit of a clustering model.
#'
#' @slot icl_model an \code{icl_model} describing the fitted model
#' @slot K a numeric vector of length 1 which correspond to the number of clusters
#' @slot count a numeric vector of length K which store the counts for each cluster
setClass("icl_fit",slots = list(icl_model="icl_model",K="numeric",counts="matrix"))

check_sbm_fit <- function(object) {
  errors = character()
  length_count = length(object@counts)
  if (length(errors) == 0) TRUE else errors
}

#' An S4 class to represent an icl fit of stochastick block model that extend \code{icl_fit}.
#'
#' @slot RA a matrix of dimension K by K which store the edges counts for pairs of clusters
setClass("sbm_fit",slots = list(x_counts="matrix",cl="matrix"),contains="icl_fit",validity = check_sbm_fit)




#' An S4 class to represent a mixture of multinomial also known has mixture of unigrams that extends \code{icl_model} class.
#'
setClass("alg",slots = list(name = "character"))



#' An S4 class  extends \code{alg} class.
#'
#' @slot a0 a numeric vector of length 1 which define the parameters of the beta prior over the edges (default to 1)
setClass("greed",
         contains = "alg",
         representation =  list(cl="vector",nb_start="numeric"),
         prototype(name="greed",nb_start=10))

setClass("genetic",
         contains = "alg",
         representation =  list(pop_size = "numeric"),
         prototype(name="genetic",pop_size=10))

setGeneric("fit", function(model,x,K,alg) standardGeneric("fit")) 
setMethod(f = "fit", 
          signature = signature("sbm","dgCMatrix", "numeric","greed"), 
          definition = function(model, x, K,alg){
            future::plan(multisession)
            solutions = listenv::listenv()
            for (i in 1:alg@nb_start){
              solutions[[i]]%<-%fit_greed_sbm(model,x,K)  
            }
            solutions = as.list(solutions)
            icls = sapply(solutions,function(s){s$icl})
            solutions[[order(icls,decreasing = TRUE)[1]]]  
          })

setMethod(f = "fit", 
          signature = signature("mm","dgCMatrix", "numeric","greed"), 
          definition = function(model, x, K,alg){
            future::plan(multisession)
            solutions = listenv::listenv()
            for (i in 1:alg@nb_start){
              solutions[[i]]%<-%fit_greed_mm(model,x,K)  
            }
            solutions = as.list(solutions)
            icls = sapply(solutions,function(s){s$icl})
            solutions[[order(icls,decreasing = TRUE)[1]]]
          })

setMethod(f = "fit", 
          signature = signature("sbm","dgCMatrix", "numeric","genetic"), 
          definition = function(model, x, K,alg){
            future::plan(multisession)
            solutions = listenv::listenv()
            # première generation
            pop_size = alg@pop_size
            for (i in 1:pop_size){
                solutions[[i]] %<-% fit_greed_sbm(model,x,K)
            }
            solutions = as.list(solutions)
            icls = sapply(solutions,function(s){s$icl})
            while((max(icls)-min(icls))>1){
              # sélections 
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
              icls = sapply(solutions,function(s){s$icl})
              
            }
            solutions[[order(icls,decreasing = TRUE)[1]]]
            
})


cross_over = function(sol1,sol2,model,x){
  ncl = unclass(factor(paste(sol1$cl,sol2$cl)))
  fit_icl_init(model,x,max(ncl),ncl-1)
}
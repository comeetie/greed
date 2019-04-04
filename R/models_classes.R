#' @title Clustering models classes
#' @aliases icl_model icl_model-class sbm sbm-class mm mm-class
#' @name models-classes
NULL


#' @rdname models-classes
#' @title icl_model
#' 
#' An S4 class to represent an abstract clustering model
#' \itemize{
#' \item slots : \code{name,alpha}
#' }
#' @slot name a character vector
#' @slot alpha a numeric vector of length 1 which define the parameters of the dirichlet over the cluster proportions (default to 1)
#' @export
setClass("icl_model",slots = list(name = "character",alpha = "numeric"))

#' @rdname models-classes
#' @title sbm
#' 
#' An S4 class to represent a stochastick block model that extends \code{icl_model} class.
#' \itemize{
#' \item slots : \code{name,alpha,a0,b0}
#' }
#' @slot a0 a numeric vector of length 1 which define the parameters of the beta prior over the edges (default to 1)
#' @slot b0 a numeric vector of length 1 which define the parameters of the beta prior over the non-edges (default to 1)
#' @examples 
#' new("sbm")
#' new("sbm",a0=0.5,b0=0.5,alpha=1)
#' @export
setClass("sbm",
         representation = list(a0 = "numeric",b0="numeric"),
         contains = "icl_model",
         prototype(name="sbm",a0=1,b0=1,alpha=1))


#' @rdname models-classes
#' @title sbm
#' 
#' An S4 class to represent a stochastick block model that extends \code{icl_model} class.
#' \itemize{
#' \item slots : \code{name,alpha,a0,b0}
#' }
#' @slot a0 a numeric vector of length 1 which define the parameters of the beta prior over the edges (default to 1)
#' @slot b0 a numeric vector of length 1 which define the parameters of the beta prior over the non-edges (default to 1)
#' @examples 
#' new("sbm")
#' new("sbm",a0=0.5,b0=0.5,alpha=1)
#' @export
setClass("dcsbm",
         representation = list(),
         contains = "icl_model",
         prototype(name="dcsbm",alpha=1))



#' @rdname models-classes
#' @title mm
#' 
#' An S4 class to represent a mixture of multinomial also known has mixture of unigrams that extends \code{icl_model} class.
#' \itemize{
#' \item slots : \code{name,alpha,beta}
#' }
#' @slot beta a numeric vector of length 1 which define the parameters of the beta prior over the counts (default to 1)
#' @examples
#' new("mm")
#' new("mm",alpha=1,beta=1)
#' @export
setClass("mm",
         representation = list(beta = "numeric"),
         contains = "icl_model",
         prototype(name="mm",beta=1,alpha=1))



#' @rdname models-classes
#' @title mm
#' 
#' An S4 class to represent a mixture of multinomial also known has mixture of unigrams that extends \code{icl_model} class.
#' \itemize{
#' \item slots : \code{name,alpha,reg,a0,b0}
#' }
#' @slot reg a numeric vector of length 1 which define the variance parameter of the normal prior over the regression parameters (default to 0.1)
#' @slot a0 a numeric vector of length 1 which define the parameter a0 of the inverse gamma over the regression noise variance parameters (default to 1)
#' @slot b0 a numeric vector of length 1 which define the parameter b0 of the inverse gamma prior over the regression noise variance parameters (default to 1)
#' @examples
#' new("mreg")
#' new("mreg",alpha=1,reg=0.8,a0=0.5,b0=0.5)
#' @export
setClass("mreg",
         representation = list(reg = "numeric",a0="numeric",b0="numeric"),
         contains = "icl_model",
         prototype(name="mreg",reg=0.1,a0=1,b0=1,alpha=1))
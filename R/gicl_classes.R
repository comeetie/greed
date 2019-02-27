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


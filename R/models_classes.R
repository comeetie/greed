

#' @title abstract class to represent a generative model
#' 
#' An S4 class to represent an abstract generative model
#' @slot name a character vector
#' @slot alpha a numeric vector of length 1 which define the parameters of the Dirichlet over the cluster proportions (default to 1)
#' @seealso  \code{\link{Lca-class}}, \code{\link{Gmm-class}}, \code{\link{Sbm-class}}, \code{\link{DcSbm-class}},...
#' @export
setClass("DlvmPrior",slots = list(alpha = "numeric"))






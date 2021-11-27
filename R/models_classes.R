

#' @title abstract class to represent a generative model for clustering
#' 
#' An S4 class to represent an abstract generative model
#' @slot alpha a numeric vector of length 1 which define the parameters of the Dirichlet over the cluster proportions (default to 1)
#' @family DlvmModels
#' @export
setClass("DlvmPrior",slots = list(alpha = "numeric"))


#' @title abstract class to represent a generative model for co-clustering
#' 
#' An S4 class to represent an abstract generative model
#' @slot alpha a numeric vector of length 1 which define the parameters of the Dirichlet over the cluster proportions (default to 1)
#' @family DlvmCoModels
#' @export
setClass("DlvmCoPrior",contains="DlvmPrior")






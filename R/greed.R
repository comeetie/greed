#' greed: clustering and model selection for graphs and counts data
#'
#' Greed enables the clustering of networks and count data matrix with different models. 
#' Model selection and clustering are performed in
#' combination by optimizing the Integrated Classification Likelihood.
#' Optimization is performed thanks to a combination of greedy local search and
#' a genetic algorithm. The main entry point is the \code{\link{greed}} function
#' to perform the clustering, which is documented below. The package also
#' provides sampling functions for all the implemented DLVMs.
#'
#' @family DlvmModels
#'
#' @docType package
#' @name greed
NULL
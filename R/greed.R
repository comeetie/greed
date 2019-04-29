#' greed: clustering and model selection for graphs and counts data
#'
#' Greed enable the clustering of networks and counts data such as document/term matrix with different model. 
#' Model selection and clustering is performed in combination by optimizing the Integrated Classification Likelihood 
#' (which is equivalent to minimizing the description length). 
#' Their are four models availables : 
#' \itemize{
#' \item Stochastic Block Model (directed)
#' \item Degree corected Stochastic Block Model (directed)
#' \item Mixture of Multinomials
#' \item Multivariate mixture of poissons
#' }
#' The optimization is performed thanks to a combination of greedy local search and a genetic algorithm. 
#' The main entry point is the \code{\link{greed}} function to perfom the clustering.
#' 
#' 
#'
#' @docType package
#' @name greed
NULL
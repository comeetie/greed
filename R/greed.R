#' greed: clustering and model selection for graphs and counts data
#'
#' Greed enable the clustering of networks and counts data such as document/term matrix with different model. 
#' Model selection and clustering is performed in combination by optimizing the Integrated Classification Likelihood 
#' (which is equivalent to minimizing the description length). 
#' Their following generative models are available : 
#' \itemize{
#' \item Stochastic Block Model (directed)
#' \item Degree corrected Stochastic Block Model (directed)
#' \item Stochastic Block Model with Mutlinomial observations (experimental)
#' \item Mixture of Multinomial
#' \item Degree corrected Latent Block Model (directed)
#' \item Gaussian Mixture Model (experimental)
#' \item Mixture of Multivariate Regressions (experimental)
#' }
#' The optimization is performed thanks to a combination of greedy local search and a genetic algorithm. 
#' The main entry point is the \code{\link{greed}} function to perform the clustering.
#' 
#' 
#'
#' @docType package
#' @name greed
NULL
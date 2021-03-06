% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alg_classes.R, R/greed.R
\docType{package}
\name{greed}
\alias{greed}
\title{Model based hierarchical clustering}
\usage{
greed(
  X,
  K = 20,
  model = find_model(X),
  alg = methods::new("hybrid"),
  verbose = FALSE
)
}
\arguments{
\item{X}{data to cluster either a matrix,an array or a \code{\link{dgCMatrix-class}}}

\item{K}{initial number of cluster}

\item{model}{a generative model to fit \code{\link{sbm-class}}, \code{\link{dcsbm-class}}, \code{\link{co_dcsbm-class}}, \code{\link{mm-class}},\code{\link{gmm-class}}, \code{\link{diaggmm-class}} or \code{\link{mvmreg-class}}}

\item{alg}{an optimization algorithm of class \code{\link{hybrid-class}} (default), \code{\link{multistarts-class}}, \code{\link{seed-class}} or \code{\link{genetic-class}}}

\item{verbose}{Boolean for verbose mode}
}
\value{
an \code{\link{icl_path-class}} object
}
\description{


Greed enables the clustering of networks and count data such as document-term
matrix with different models. Model selection and clustering are performed in
combination by optimizing the Integrated Classification Likelihood.
Optimization is performed thanks to a combination of greedy local search and
a genetic algorithm. The main entry point is the \code{\link{greed}} function
to perform the clustering, which is documented below. The package also
provides sampling functions for all the implemented DLVMs: Mixture of
Multinomials (\code{\link{rmm}}), Stochastic Block Model (\code{\link{rsbm}},
\code{\link{rdcsbm}}) and Latent Block Model (\code{\link{rlbm}}).
}

#' Generate a graph adjacency matrix using a Stochastic Block Model
#'
#' \code{rsbm} returns the adjacency matrix and the cluster labels generated randomly with a Stochastic Block Model.
#'
#' This function takes the desired graph size, cluster proportions and connectivity matrix as input and sample a graph accordingly together with the clusters labels.
#'
#' @param N The size of the graph to generate
#' @param pi A numeric vector of length K with clusters proportions (will be normalized to sum up to 1).
#' @param mu A numeric matrix of dim K x K with the connectivity pattern to generate. elements in [0,1].
#' @return A list with fields:
#' \itemize{
#' \item x: the graph adjacency matrix as a \code{dgCMatrix}
#' \item K: number of generated clusters
#' \item N: number of vertex
#' \item cl: vector of clusters labels
#' \item pi: clusters proportions
#' \item mu: connectivity matrix
#' }
#' @examples
#' simu <- rsbm(100, rep(1 / 5, 5), diag(rep(0.1, 5)) + 0.001)
#' @export
rsbm <- function(N, pi, mu) {
  K <- length(pi)
  cl <- sample(1:K, N, replace = TRUE, prob = pi)
  x <- matrix(stats::rbinom(N * N, 1, mu[cbind(rep(cl, N), rep(cl, each = N))]), N, N)
  links <- Matrix::which(x == 1, arr.ind = TRUE)
  list(cl = cl, x = Matrix::sparseMatrix(links[, 1], links[, 2], x = rep(1, nrow(links))), K = K, N = N, pi = pi, mu = mu)
}

#' Generate a data matrix using a Latent Block Model
#'
#' \code{rlbm} returns the adjacency matrix and the cluster labels generated randomly with a Latent Block Model.
#'
#' This function takes the desired graph size, cluster proportions and connectivity matrix as input and sample a graph accordingly together with the clusters labels.
#'
#' @param Nr desired Number of rows
#' @param Nc desired Number of column
#' @param pir A numeric vector of length Kr with rows clusters proportions (will be normalized to sum up to 1).
#' @param pic A numeric vector of length Kc with columns clusters proportions (will be normalized to sum up to 1).
#' @param mu A numeric matrix of dim Kr x Kc with the connectivity pattern to generate. elements in [0,1].
#' @return A list with fields:
#' \itemize{
#' \item x: the generated data matrix as a \code{dgCMatrix}
#' \item clr: vector of row clusters labels
#' \item clc: vector of column clusters labels
#' \item Kr: number of generated row clusters
#' \item Kc: number of generated column clusters
#' \item Nr: number of rows
#' \item Nc: number of column
#' \item pir: row clusters proportions
#' \item pic: column clusters proportions
#' \item mu: connectivity matrix
#' }
#' @examples
#' simu <- rlbm(500, 1000, rep(1 / 5, 5), rep(1 / 10, 10), matrix(runif(50), 5, 10))
#' @export
rlbm <- function(Nr, Nc, pir, pic, mu) {
  Kr <- length(pir)
  Kc <- length(pic)
  clc <- sample(1:Kc, Nc, replace = TRUE, prob = pic)
  clr <- sample(1:Kr, Nr, replace = TRUE, prob = pir)
  if (dim(mu)[1] != Kr | dim(mu)[1] != Kr) {
    stop("incompatible argument sizes")
  }
  x <- matrix(stats::rbinom(Nr * Nc, 1, mu[cbind(rep(clr, Nc), rep(clc, each = Nr))]), Nr, Nc)
  links <- Matrix::which(x == 1, arr.ind = TRUE)
  list(clr = clr, clc = clc, x = Matrix::sparseMatrix(links[, 1], links[, 2], x = rep(1, nrow(links))), Kr = Kr, Kc = Kc, Nr = Nr, Nc = Nc, pir = pir, pic = pic, mu = mu)
}



#' Generate data using a Multinomial Mixture
#'
#' \code{rmm} returns a count matrix and the cluster labels generated randomly with a Mixture of Multinomial model.
#'
#' It takes the sample size, cluster proportions and emission matrix, and  as input and sample a graph accordingly together with the clusters labels.
#'
#' @param N A numeric value the size of the graph to generate
#' @param pi A numeric vector of length K with clusters proportions. Must sum up to 1.
#' @param mu A numeric matrix of dim k x D with the clusters patterns to generate, all elements in [0,1].
#' @param lambda A numeric value which specify the expectation for the row sums.
#' @return A list with fields:
#' \itemize{
#' \item x: the count matrix as a \code{dgCMatrix}
#' \item K: number of generated clusters
#' \item N: number of vertex
#' \item cl: vector of clusters labels
#' \item pi: clusters proportions
#' \item mu: connectivity matrix
#' \item lambda: expectation of row sums
#' }
#' @export
rmm <- function(N, pi, mu, lambda) {
  K <- length(pi)
  nbv <- dim(mu)[2]
  cl <- sample(1:K, N, replace = TRUE, prob = pi)
  X <- matrix(0, N, nbv)
  if (length(lambda) == 1 | length(lambda) != N) {
    Nr <- stats::rpois(N, lambda)
  } else {
    Nr <- lambda
  }

  for (k in 1:K) {
    X[cl == k] <- t(stats::rmultinom(sum(cl == k), Nr[cl == k], mu[k, ]))
  }
  links <- Matrix::which(X > 0, arr.ind = TRUE)
  list(cl = cl, x = Matrix::sparseMatrix(links[, 1], links[, 2], x = X[links]), K = K, N = N, pi = pi, mu = mu, lambda = lambda)
}

#' Generates graph adjacency matrix using a degree corrected SBM
#'
#' \code{rdcsbm} returns an adjacency matrix and the cluster labels generated randomly using a Degree Corrected Stochastic Block Model.
#'
#' It takes the sample size, cluster proportions and emission matrix, and   as input and sample a graph accordingly together with the clusters labels.
#'
#' @param N A numeric value the size of the graph to generate
#' @param pi A numeric vector of length K with clusters proportions. Must sum up to 1.
#' @param mu A numeric matrix of dim K x K with the connectivity pattern to generate, elements in [0,1].
#' @param betain A numeric vector of length N which specify the in-degree correction will be normalized per cluster during the generation.
#' @param betaout A numeric vector of length N which specify the out-degree correction will be normalized per cluster during the generation.
#' @return A list with fields:
#' \itemize{
#' \item x: the count matrix as a \code{dgCMatrix}
#' \item K: number of generated clusters
#' \item N: number of vertex
#' \item cl: vector of clusters labels
#' \item pi: clusters proportions
#' \item mu: connectivity matrix
#' \item betain: normalized in-degree parameters
#' \item betaout: normalized out-degree parameters
#' }
#' @export
rdcsbm <- function(N, pi, mu, betain, betaout) {
  K <- length(pi)
  cl <- sample(1:K, N, replace = TRUE, prob = pi)
  betain_cl <- stats::aggregate(betain, list(cl), mean)$x
  betain <- betain / betain_cl[cl]
  betaout_cl <- stats::aggregate(betaout, list(cl), mean)$x
  betaout <- betaout / betaout_cl[cl]
  x <- matrix(stats::rpois(N * N, mu[cbind(rep(cl, N), rep(cl, each = N))] * rep(betain, N) * rep(betaout, each = N)), N, N)
  links <- Matrix::which(x > 0, arr.ind = TRUE)
  list(cl = cl, x = Matrix::sparseMatrix(links[, 1], links[, 2], x = x[links],dims = c(N,N)), K = K, N = N, pi = pi, mu = mu)
}


#' Generate data from a mixture of regression model
#'
#' \code{rmreg} returns an X matrix, a y vector and the cluster labels generated randomly with a Mixture of regression model.
#'
#' It takes the sample size, cluster proportions and regression parameters matrix and variance  as input accordingly
#'
#' @param N A numeric value the size of the graph to generate
#' @param pi A numeric vector of length K with clusters proportions (must sum up to 1)
#' @param A A numeric matrix of dim K x d with the regression coefficient
#' @param sigma A numeric of length 1 with the target conditional variance
#' @param X A matrix of covariate
#' @return A list with fields:
#' \itemize{
#' \item X: the covariate matrix
#' \item y: the target feature
#' \item K: number of generated clusters
#' \item N: sample size
#' \item cl: vector of clusters labels
#' \item pi: clusters proportions
#' \item A: regression coefficients used in the simulation
#' \item sigma: conditional variance
#' }
#' @export
rmreg <- function(N, pi, A, sigma, X = cbind(rep(1, N), matrix(stats::rnorm(N * (ncol(A) - 1)), N, ncol(A) - 1))) {
  K <- length(pi)
  cl <- sample(1:K, N, replace = TRUE, prob = pi)
  yt <- X %*% t(A) + stats::rnorm(N, 0, sigma)
  y <- yt[cbind(1:N, cl)]
  list(cl = cl, X = X, y = y, K = K, N = N, pi = pi, A = A, sigma = sigma)
}

#' Generate a graph adjacency matrix using a Stochastic Block Model
#'
#' \code{rmultsbm} returns the multi-graph adjacency matrix and the cluster labels generated randomly with a Multinomial Stochastic Block Model.
#'
#' This function takes the desired graph size, cluster proportions and connectivity matrix as input and sample a graph accordingly together with the clusters labels.
#'
#' @param N The size of the graph to generate
#' @param pi A numeric vector of length K with clusters proportions (will be normalized to sum up to 1).
#' @param mu A numeric array of dim K x K x M with the connectivity pattern to generate. elements in [0,1].
#' @param lambda A double with the Poisson intensity to generate the total counts
#' @return A list with fields:
#' \itemize{
#' \item x: the multi-graph adjacency matrix as an \code{array}
#' \item K: number of generated clusters
#' \item N: number of vertex
#' \item cl: vector of clusters labels
#' \item pi: clusters proportions
#' \item mu: connectivity matrix
#' \item lambda:
#' }
#' @examples
#' simu <- rsbm(100, rep(1 / 5, 5), diag(rep(0.1, 5)) + 0.001)
#' @export
rmultsbm <- function(N, pi, mu, lambda) {
  K <- length(pi)
  cl <- sample(1:K, N, replace = TRUE, prob = pi)
  x <- array(dim = c(N, N, dim(mu)[3]))
  for (i in 1:N) {
    for (j in 1:N) {
      x[i, j, ] <- stats::rmultinom(1, stats::rpois(1, lambda), mu[cl[i], cl[j], ])
    }
  }
  list(cl = cl, x = x, K = K, N = N, pi = pi, mu = mu)
}


#' Generate data from lca model
#'
#' \code{rlca} returns a data.frame with factor sampled from an lca model
#'
#' This function takes the desired graph size, cluster proportions and connectivity matrix as input and sample a graph accordingly together with the clusters labels.
#'
#' @param N The size of the graph to generate
#' @param pi A numeric vector of length K with clusters proportions (will be normalized to sum up to 1).
#' @param theta A list of size V
#' @return A list with fields:
#' \itemize{
#' \item x: the multi-graph adjacency matrix as an \code{array}
#' \item K: number of generated clusters
#' \item N: number of vertex
#' \item cl: vector of clusters labels
#' \item pi: clusters proportions
#' \item theta:
#' }
#' @examples
#' theta <- list(
#'   matrix(c(0.1, 0.9, 0.9, 0.1, 0.5, 0.5, 0.3, 0.7), ncol = 2, byrow = TRUE),
#'   matrix(c(0.5, 0.5, 0.3, 0.7, 0.05, 0.95, 0.3, 0.7), ncol = 2, byrow = TRUE),
#'   matrix(c(0.5, 0.5, 0.9, 0.1, 0.5, 0.5, 0.1, 0.9), ncol = 2, byrow = TRUE)
#' )
#' lca.data <- rlca(100, rep(1 / 4, 4), theta)
#' @export
rlca <- function(N, pi, theta) {
  K <- length(pi)
  cl <- sample(1:K, N, replace = TRUE, prob = pi)
  V <- length(theta)
  x <- data.frame(matrix(NA, N, V))
  for (v in 1:V) {
    for (k in 1:K) {
      x[cl == k, v] <- sample(1:ncol(theta[[v]]), sum(cl == k), prob = theta[[v]][k, ], replace = TRUE)
    }
    x[, v] <- factor(x[, v], levels = 1:ncol(theta[[v]]))
  }
  list(cl = cl, x = x, K = K, N = N, pi = pi, theta = theta)
}

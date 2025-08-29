context("SBM test")
library(greed)
library(ggplot2)
library(Matrix)
set.seed(1234)

test_that("Combined models sbm and gmm", {
  N <- 500
  K <- 10
  pi <- rep(1 / K, K)
  mu <- diag(rep(1 / 5, K))
  sbm <- rsbm(N, pi, mu)
  gmm <- do.call(cbind, lapply(1:K, function(x) {
    rnorm(N, 20 * runif(1) - 10)
  }))
  Xnodes <- as.matrix(gmm[cbind(1:N, sbm$cl)])

  Xinput <- list(graph = sbm$x, Xnodes = Xnodes)
  Mtt <- CombinedModels(models = list(graph = SbmPrior(), Xnodes = GmmPrior()))
  sol <- greed(Xinput, model = Mtt)
  expect_equal(sol@K, K)
  solc <- cut(sol, 8)
})


# test_that("Combined models multsbm and gmm", {
#   N <- 100
#   K <- 3
#   pi <- rep(1 / K, K)
#   mu <- array(dim = c(K, K, 3))
#   mu[, , 1] <- diag(rep(1 / 5, K)) + runif(K^2) * 0.005
#   mu[1, 1, 1] <- runif(1) * 0.005
#   mu[, , 2] <- diag(rep(1 / 5, K)) + runif(K^2) * 0.005
#   mu[2, 2, 1] <- runif(1) * 0.005
#   mu[, , 3] <- 1 - mu[, , 1] - mu[, , 2]
#   lambda <- 10
#   multsbm <- rmultsbm(N, pi, mu, 10)
#   gmm <- do.call(cbind, lapply(1:K, function(x) {
#     rnorm(N, 20 * runif(1) - 10)
#   }))
#   Xnodes <- as.matrix(gmm[cbind(1:N, multsbm$cl)])
# 
#   Xinput <- list(graph = multsbm$x, Xnodes = Xnodes)
#   Mtt <- CombinedModels(models = list(graph = MultSbmPrior(), Xnodes = GmmPrior()))
#   sol <- greed(Xinput, model = Mtt)
#   expect_equal(sol@K, K)
#   solc <- cut(sol, 2)
# })
# 
# test_that("Combined models mom and gmm", {
#   N <- 200
#   K <- 4
#   pi <- rep(1 / K, K)
#   mu <- cbind(diag(rep(5, K)), matrix(0, K, 20)) + matrix(runif(K * (20 + K)), K, 20 + K)
#   mm <- rmm(N, pi, mu, 30)
#   gmm <- do.call(cbind, lapply(1:K, function(x) {
#     rnorm(N, 20 * runif(1) - 10)
#   }))
#   Xnodes <- as.matrix(gmm[cbind(1:N, mm$cl)])
#   Xinput <- list(mom = mm$x, Xnodes = Xnodes)
#   Mtt <- CombinedModels(models = list(mom = MoMPrior(), Xnodes = GmmPrior()))
#   sol <- greed(Xinput, model = Mtt)
#   expect_equal(sol@K, K)
#   solc <- cut(sol, 2)
# })
# 
# test_that("Combined models lca and gmm", {
#   N <- 500
#   theta <- list(
#     matrix(c(0.1, 0.9, 0.9, 0.1, 0.8, 0.2, 0.05, 0.95), ncol = 2, byrow = TRUE),
#     matrix(c(0.95, 0.05, 0.3, 0.7, 0.05, 0.95, 0.05, 0.95), ncol = 2, byrow = TRUE),
#     matrix(c(0.95, 0.04, 0.01, 0.9, 0.09, 0.01, 0.01, 0.01, 0.98, 0.9, 0.05, 0.05), ncol = 3, byrow = TRUE),
#     matrix(c(1, 0, 0, 1, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
#   )
#   lca.data <- rlca(N, rep(1 / 4, 4), theta)
#   K <- 4
#   gmm <- do.call(cbind, lapply(1:K, function(x) {
#     rnorm(N, 20 * runif(1) - 10)
#   }))
#   Xnodes <- as.matrix(gmm[cbind(1:N, lca.data$cl)])
#   Xinput <- list(lca = lca.data$x, Xnodes = Xnodes)
#   Mtt <- CombinedModels(models = list(lca = LcaPrior(), Xnodes = GmmPrior()))
#   sol <- greed(Xinput, model = Mtt)
#   solc <- cut(sol, 2)
# })

# test_that("Combined models lca and diaggmm", {
#   N <- 500
#   theta <- list(
#     matrix(c(0.1, 0.9, 0.9, 0.1, 0.8, 0.2, 0.05, 0.95), ncol = 2, byrow = TRUE),
#     matrix(c(0.95, 0.05, 0.3, 0.7, 0.05, 0.95, 0.05, 0.95), ncol = 2, byrow = TRUE),
#     matrix(c(0.95, 0.04, 0.01, 0.9, 0.09, 0.01, 0.01, 0.01, 0.98, 0.9, 0.05, 0.05), ncol = 3, byrow = TRUE),
#     matrix(c(1, 0, 0, 1, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
#   )
#   lca.data <- rlca(N, rep(1 / 4, 4), theta)
#   K <- 4
#   gmm <- do.call(cbind, lapply(1:K, function(x) {
#     rnorm(N, 20 * runif(1) - 10)
#   }))
#   Xnodes <- as.matrix(gmm[cbind(1:N, lca.data$cl)])
#   Xinput <- list(lca = lca.data$x, Xnodes = Xnodes)
#   Mtt <- CombinedModels(models = list(lca = LcaPrior(), Xnodes = DiagGmmPrior()))
#   sol <- greed(Xinput, model = Mtt)
#   expect_equal(sol@K, K)
#   solc <- cut(sol, 2)
# })

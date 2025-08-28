context("MULTSBM test")
library(greed)
library(ggplot2)
set.seed(1234)

test_that("MULTSBM sim", {
  N <- 100
  K <- 3
  pi <- rep(1 / K, K)
  mu <- array(dim = c(K, K, 3))
  mu[, , 1] <- diag(rep(1 / 5, K)) + runif(K^2) * 0.005
  mu[1, 1, 1] <- runif(1) * 0.005
  mu[, , 2] <- diag(rep(1 / 5, K)) + runif(K^2) * 0.005
  mu[2, 2, 1] <- runif(1) * 0.005
  mu[, , 3] <- 1 - mu[, , 1] - mu[, , 2]
  lambda <- 10
  multsbm <- rmultsbm(N, pi, mu, 10)
  expect_equal(dim(multsbm$x)[1], N)
  expect_equal(dim(multsbm$x)[2], N)
  expect_equal(length(multsbm$cl), N)
  expect_gte(min(multsbm$cl), 1)
  expect_lte(max(multsbm$cl), K)
  model <- MultSbm()
  data <- greed:::preprocess(model, multsbm$x)
  i <- sample(100, 1)
  oldcl <- multsbm$cl[i]
  newcl <- sample(setdiff(1:K, oldcl), 1)
  expect_lte(greed:::test_swap(model, data, multsbm$cl, i, newcl), 10^-6)
  expect_lte(greed:::test_merge(model, data, multsbm$cl, oldcl, newcl), 10^-6)
  expect_lte(max(abs(greed:::test_merge_correction(model, data, multsbm$cl, oldcl, newcl))), 10^-6)
  # undirected 
  multsbm$x[, , 1] <- matrix(tril(multsbm$x[, , 1]) + t(tril(multsbm$x[, , 1])))
  diag(multsbm$x[, , 1]) <- 0
  multsbm$x[, , 2] <- matrix(tril(multsbm$x[, , 2]) + t(tril(multsbm$x[, , 2])))
  diag(multsbm$x[, , 2]) <- 0
  multsbm$x[, , 3] <- matrix(tril(multsbm$x[, , 3]) + t(tril(multsbm$x[, , 3])))
  diag(multsbm$x[, , 3]) <- 0
  model <- MultSbm(type="undirected")
  data <- greed:::preprocess(model, multsbm$x)
  expect_lte(greed:::test_swap(model, data, multsbm$cl, i, newcl), 10^-6)
  expect_lte(greed:::test_merge(model, data, multsbm$cl, oldcl, newcl), 10^-6)
  expect_lte(max(abs(greed:::test_merge_correction(model, data, multsbm$cl, oldcl, newcl))), 10^-6)
})


test_that("MULTSBM hybrid directed", {
  N <- 100
  K <- 3
  pi <- rep(1 / K, K)
  mu <- array(dim = c(K, K, 3))
  mu[, , 1] <- diag(rep(1 / 5, K)) + runif(K^2) * 0.005
  mu[1, 1, 1] <- runif(1) * 0.005
  mu[, , 2] <- diag(rep(1 / 5, K)) + runif(K^2) * 0.005
  mu[2, 2, 1] <- runif(1) * 0.005
  mu[, , 3] <- 1 - mu[, , 1] - mu[, , 2]
  lambda <- 10
  multsbm <- rmultsbm(N, pi, mu, 10)
  sol <- greed(multsbm$x, model = MultSbm())
  expect_equal(sol@K, K)
  solc <- cut(sol, 2)
  # expect_true(is.ggplot(plot(solc, type = "tree")))
  # expect_true(is.ggplot(plot(solc, type = "path")))
  # expect_true(is.ggplot(plot(solc, type = "front")))
  # expect_true(is.ggplot(plot(solc, type = "blocks")))
  # expect_true(is.ggplot(plot(solc, type = "nodelink")))
  co <- coef(sol)
  expect_true(all(dim(co$thetakl) == 3))
  expect_equal(length(co$pi), 3)
  expect_equal(sum(co$pi), 1)
  expect_true(all(rowSums(co$thetakl, 2) == 3))
})

test_that("MULTSBM seed directed", {
  N <- 100
  K <- 3
  pi <- rep(1 / K, K)
  mu <- array(dim = c(K, K, 3))
  mu[, , 1] <- diag(rep(1 / 5, K)) + runif(K^2) * 0.005
  mu[1, 1, 1] <- runif(1) * 0.005
  mu[, , 2] <- diag(rep(1 / 5, K)) + runif(K^2) * 0.005
  mu[2, 2, 1] <- runif(1) * 0.005
  mu[, , 3] <- 1 - mu[, , 1] - mu[, , 2]
  lambda <- 10
  multsbm <- rmultsbm(N, pi, mu, 10)
  sol <- greed(multsbm$x, model = MultSbm(), alg = Seed())
  expect_gte(sol@K, K - 2)
  expect_lte(sol@K, K + 2)
  solc <- cut(sol, 2)
  # expect_true(is.ggplot(plot(solc, type = "tree")))
  # expect_true(is.ggplot(plot(solc, type = "path")))
  # expect_true(is.ggplot(plot(solc, type = "front")))
  # expect_true(is.ggplot(plot(solc, type = "blocks")))
  # expect_true(is.ggplot(plot(solc, type = "nodelink")))
})




test_that("MULTSBM sim", {
  N <- 100
  K <- 3
  pi <- rep(1 / K, K)
  mu <- array(dim = c(K, K, 3))
  mu[, , 1] <- diag(rep(1 / 5, K)) + runif(K^2) * 0.005
  mu[1, 1, 1] <- runif(1) * 0.005
  mu[, , 2] <- diag(rep(1 / 5, K)) + runif(K^2) * 0.005
  mu[2, 2, 1] <- runif(1) * 0.005
  mu[, , 3] <- 1 - mu[, , 1] - mu[, , 2]
  lambda <- 10
  multsbm <- rmultsbm(N, pi, mu, 10)
  expect_equal(dim(multsbm$x)[1], N)
  expect_equal(dim(multsbm$x)[2], N)
  expect_equal(length(multsbm$cl), N)
  expect_gte(min(multsbm$cl), 1)
  expect_lte(max(multsbm$cl), K)
})


test_that("MULTSBM hybrid undirected", {
  N <- 100
  K <- 3
  pi <- rep(1 / K, K)
  mu <- array(dim = c(K, K, 3))
  mu[, , 1] <- diag(rep(1 / 5, K)) + runif(K^2) * 0.005
  mu[1, 1, 1] <- runif(1) * 0.005
  mu[, , 2] <- diag(rep(1 / 5, K)) + runif(K^2) * 0.005
  mu[2, 2, 1] <- runif(1) * 0.005
  mu[, , 3] <- 1 - mu[, , 1] - mu[, , 2]
  lambda <- 10
  multsbm <- rmultsbm(N, pi, mu, 10)
  multsbm$x[, , 1] <- matrix(tril(multsbm$x[, , 1]) + t(tril(multsbm$x[, , 1])))
  diag(multsbm$x[, , 1]) <- 0
  multsbm$x[, , 2] <- matrix(tril(multsbm$x[, , 2]) + t(tril(multsbm$x[, , 2])))
  diag(multsbm$x[, , 2]) <- 0
  multsbm$x[, , 3] <- matrix(tril(multsbm$x[, , 3]) + t(tril(multsbm$x[, , 3])))
  diag(multsbm$x[, , 3]) <- 0
  sol <- greed(multsbm$x, model = MultSbm(type = "undirected"))
  expect_equal(sol@K, K)
  solc <- cut(sol, 2)
  # expect_true(is.ggplot(plot(solc, type = "tree")))
  # expect_true(is.ggplot(plot(solc, type = "path")))
  # expect_true(is.ggplot(plot(solc, type = "front")))
  # expect_true(is.ggplot(plot(solc, type = "blocks")))
  # expect_true(is.ggplot(plot(solc, type = "nodelink")))
  co <- coef(sol)
  expect_true(all(dim(co$thetakl) == 3))
  expect_equal(length(co$pi), 3)
  expect_equal(sum(co$pi), 1)
  expect_true(all(rowSums(co$thetakl, 2) == 3))
})

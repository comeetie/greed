context("SBM test")
library(greed)
library(ggplot2)
library(Matrix)
set.seed(1234)

test_that("SBM sim", {
  N <- 200
  K <- 6
  pi <- rep(1 / K, K)
  mu <- diag(rep(1 / 5, K))
  sbm <- rsbm(N, pi, mu)
  expect_equal(dim(sbm$x)[1], N)
  expect_equal(dim(sbm$x)[2], N)
  expect_equal(length(sbm$cl), N)
  expect_gte(min(sbm$cl), 1)
  expect_lte(max(sbm$cl), K)
  model <- Sbm()
  data <- greed:::preprocess(model, sbm$x)
  i <- sample(100, 1)
  oldcl <- sbm$cl[i]
  newcl <- sample(setdiff(1:K, oldcl), 1)
  expect_lte(greed:::test_swap(model, data, sbm$cl, i, newcl), 10^-6)
  expect_lte(greed:::test_merge(model, data, sbm$cl, oldcl, newcl), 10^-6)
  expect_lte(max(abs(greed:::test_merge_correction(model, data, sbm$cl, oldcl, newcl))), 10^-6)
  x <- tril(sbm$x) + t(tril(sbm$x))
  diag(x) <- 0
  model <- Sbm(type="undirected")
  data <- greed:::preprocess(model,x)
  expect_lte(greed:::test_swap(model, data, sbm$cl, i, newcl), 10^-6)
  expect_lte(greed:::test_merge(model, data, sbm$cl, oldcl, newcl), 10^-6)
  expect_lte(max(abs(greed:::test_merge_correction(model, data, sbm$cl, oldcl, newcl))), 10^-6)
  
  
})

test_that("SBM hybrid", {
  N <- 100
  K <- 3
  pi <- rep(1 / K, K)
  mu <- diag(rep(1 / 5, K)) + runif(K * K) * 0.01
  sbm <- rsbm(N, pi, mu)
  sol <- greed(sbm$x, model = Sbm())
  co <- coef(sol)
  expect_true(all(dim(co$thetakl) == c(3, 3)))
  expect_equal(sum(co$pi), 1)
  expect_equal(length(co$pi), 3)
  expect_equal(sol@K, K)
  solc <- cut(sol, 2)
  # expect_true(is.ggplot(plot(solc, type = "tree")))
  # expect_true(is.ggplot(plot(solc, type = "path")))
  # expect_true(is.ggplot(plot(solc, type = "front")))
  # expect_true(is.ggplot(plot(solc, type = "blocks")))
  # expect_true(is.ggplot(plot(solc, type = "nodelink")))
})

test_that("SBM seed", {
  N <- 100
  K <- 3
  pi <- rep(1 / K, K)
  mu <- diag(rep(1 / 5, K)) + runif(K * K) * 0.01
  sbm <- rsbm(N, pi, mu)
  sol <- greed(sbm$x, model = Sbm(), alg = Seed())
  co <- coef(sol)
  expect_true(all(dim(co$thetakl) == c(3, 3)))
  expect_equal(sum(co$pi), 1)
  expect_equal(length(co$pi), 3)
  expect_gte(sol@K, K - 2)
  expect_lte(sol@K, K + 2)
  solc <- cut(sol, 2)
  # expect_true(is.ggplot(plot(solc, type = "tree")))
  # expect_true(is.ggplot(plot(solc, type = "path")))
  # expect_true(is.ggplot(plot(solc, type = "front")))
  # expect_true(is.ggplot(plot(solc, type = "blocks")))
  # expect_true(is.ggplot(plot(solc, type = "nodelink")))
})


test_that("SBM multistart", {
  N <- 101
  K <- 3
  pi <- rep(1 / K, K)
  mu <- diag(rep(1 / 5, K)) + runif(K * K) * 0.01
  sbm <- rsbm(N, pi, mu)
  sol <- greed(sbm$x, model = Sbm(), alg = Multistarts())
  co <- coef(sol)
  expect_true(all(dim(co$thetakl) == c(3, 3)))
  expect_equal(sum(co$pi), 1)
  expect_equal(length(co$pi), 3)
  expect_gte(sol@K, K - 2)
  expect_lte(sol@K, K + 2)
  solc <- cut(sol, 2)
  # expect_true(is.ggplot(plot(solc, type = "tree")))
  # expect_true(is.ggplot(plot(solc, type = "path")))
  # expect_true(is.ggplot(plot(solc, type = "front")))
  # expect_true(is.ggplot(plot(solc, type = "blocks")))
  # expect_true(is.ggplot(plot(solc, type = "nodelink")))
})

test_that("SBM genetic", {
  N <- 100
  K <- 3
  pi <- rep(1 / K, K)
  mu <- diag(rep(1 / 5, K)) + runif(K * K) * 0.01
  sbm <- rsbm(N, pi, mu)
  sol <- greed(sbm$x, model = Sbm(), alg = Genetic())
  co <- coef(sol)
  expect_true(all(dim(co$thetakl) == c(3, 3)))
  expect_equal(sum(co$pi), 1)
  expect_equal(length(co$pi), 3)
  expect_gte(sol@K, K - 2)
  expect_lte(sol@K, K + 2)
  solc <- cut(sol, 2)
  # expect_true(is.ggplot(plot(solc, type = "tree")))
  # expect_true(is.ggplot(plot(solc, type = "path")))
  # expect_true(is.ggplot(plot(solc, type = "front")))
  # expect_true(is.ggplot(plot(solc, type = "blocks")))
  # expect_true(is.ggplot(plot(solc, type = "nodelink")))
})


test_that("SBM hybrid undirected", {
  N <- 100
  K <- 3
  pi <- rep(1 / K, K)
  mu <- diag(rep(1 / 5, K)) + runif(K * K) * 0.01
  sbm <- rsbm(N, pi, mu)
  x <- tril(sbm$x) + t(tril(sbm$x))
  diag(x) <- 0
  sol <- greed(x, model = Sbm(type = "undirected"))
  co <- coef(sol)
  expect_true(all(dim(co$thetakl) == c(3, 3)))
  expect_equal(sum(co$pi), 1)
  expect_equal(length(co$pi), 3)
  expect_equal(sol@K, K)
  solc <- cut(sol, 2)
  # expect_true(is.ggplot(plot(solc, type = "tree")))
  # expect_true(is.ggplot(plot(solc, type = "path")))
  # expect_true(is.ggplot(plot(solc, type = "front")))
  # expect_true(is.ggplot(plot(solc, type = "blocks")))
  # expect_true(is.ggplot(plot(solc, type = "nodelink")))
})

test_that("SBM seed undirected", {
  N <- 100
  K <- 3
  pi <- rep(1 / K, K)
  mu <- diag(rep(1 / 5, K)) + runif(K * K) * 0.01
  sbm <- rsbm(N, pi, mu)
  x <- tril(sbm$x) + t(tril(sbm$x))
  diag(x) <- 0
  sol <- greed(x, model = Sbm(type = "undirected"), alg = Seed())
  co <- coef(sol)
  expect_true(all(dim(co$thetakl) == c(3, 3)))
  expect_equal(sum(co$pi), 1)
  expect_equal(length(co$pi), 3)
  expect_gte(sol@K, K - 2)
  expect_lte(sol@K, K + 2)
  solc <- cut(sol, 2)
  # expect_true(is.ggplot(plot(solc, type = "tree")))
  # expect_true(is.ggplot(plot(solc, type = "path")))
  # expect_true(is.ggplot(plot(solc, type = "front")))
  # expect_true(is.ggplot(plot(solc, type = "blocks")))
  # expect_true(is.ggplot(plot(solc, type = "nodelink")))
})

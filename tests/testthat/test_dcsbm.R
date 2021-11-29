context("DCSBM test")
library(greed)
library(ggplot2)
library(Matrix)
set.seed(1234)

test_that("DCSBM sim", {
  N <- 100
  K <- 5
  pi <- rep(1 / K, K)
  mu <- diag(rep(1 / 5, K)) + runif(K * K) * 0.01
  sbm <- rdcsbm(N, pi, mu, rep(15, N), rep(15, N))
  expect_equal(dim(sbm$x)[1], N)
  expect_equal(dim(sbm$x)[2], N)
  expect_equal(length(sbm$cl), N)
  expect_gte(min(sbm$cl), 1)
  expect_lte(max(sbm$cl), K)
  model <- DcSbm()
  data <- greed:::preprocess(model, sbm$x)
  i <- sample(100, 1)
  oldcl <- sbm$cl[i]
  newcl <- sample(setdiff(1:K, oldcl), 1)
  expect_lte(greed:::test_swap(model, data, sbm$cl, i, newcl), 10^-6)
  expect_lte(greed:::test_merge(model, data, sbm$cl, oldcl, newcl), 10^-6)
  expect_lte(max(abs(greed:::test_merge_correction(model, data, sbm$cl, oldcl, newcl))), 10^-6)
  x <- tril(sbm$x) + t(tril(sbm$x))
  diag(x) <- 0
  model <- DcSbm(type="undirected")
  data <- greed:::preprocess(model,x)
  expect_lte(greed:::test_swap(model, data, sbm$cl, i, newcl), 10^-6)
  expect_lte(greed:::test_merge(model, data, sbm$cl, oldcl, newcl), 10^-6)
  expect_lte(max(abs(greed:::test_merge_correction(model, data, sbm$cl, oldcl, newcl))), 10^-6)
})


test_that("DCSBM hybrid", {
  N <- 1000
  K <- 5
  pi <- rep(1 / K, K)
  mu <- diag(rep(1 / 5, K)) + runif(K * K) * 0.01
  sbm <- rdcsbm(N, pi, mu, rep(15, N), rep(15, N))
  sol <- greed(sbm$x)
  expect_equal(sol@K, K)
  co <- coef(sol)
  expect_true(all(dim(co$thetakl) == c(5, 5)))
  expect_equal(sum(co$pi), 1)
  expect_equal(length(co$pi), 5)
  solc <- cut(sol, 4)
  expect_true(is.ggplot(plot(solc, type = "tree")))
  expect_true(is.ggplot(plot(solc, type = "path")))
  expect_true(is.ggplot(plot(solc, type = "front")))
  expect_true(is.ggplot(plot(solc, type = "blocks")))
  expect_true(is.ggplot(plot(solc, type = "nodelink")))
})

test_that("DCSBM seed", {
  N <- 100
  K <- 5
  pi <- rep(1 / K, K)
  mu <- diag(rep(1 / 5, K)) + runif(K * K) * 0.01
  sbm <- rdcsbm(N, pi, mu, rep(15, N), rep(15, N))
  sol <- greed(sbm$x, alg = Seed())
  expect_gte(sol@K, K - 2)
  expect_lte(sol@K, K + 2)
  solc <- cut(sol, 4)
  expect_true(is.ggplot(plot(solc, type = "tree")))
  expect_true(is.ggplot(plot(solc, type = "path")))
  expect_true(is.ggplot(plot(solc, type = "front")))
  expect_true(is.ggplot(plot(solc, type = "blocks")))
  expect_true(is.ggplot(plot(solc, type = "nodelink")))
})

test_that("DCSBM hybrid undirected", {
  N <- 200
  K <- 5
  pi <- rep(1 / K, K)
  mu <- diag(rep(1 / 5, K)) + runif(K * K) * 0.01
  sbm <- rdcsbm(N, pi, mu, rep(15, N), rep(15, N))
  x <- tril(sbm$x) + t(tril(sbm$x))
  diag(x) <- 0
  sol <- greed(x)
  expect_equal(sol@K, K)
  co <- coef(sol)
  expect_true(all(dim(co$thetakl) == c(5, 5)))
  expect_equal(sum(co$pi), 1)
  expect_equal(length(co$pi), 5)
  solc <- cut(sol, 4)
  expect_true(is.ggplot(plot(solc, type = "tree")))
  expect_true(is.ggplot(plot(solc, type = "path")))
  expect_true(is.ggplot(plot(solc, type = "front")))
  expect_true(is.ggplot(plot(solc, type = "blocks")))
  expect_true(is.ggplot(plot(solc, type = "nodelink")))
})

test_that("DCSBM seed undirected", {
  N <- 200
  K <- 5
  pi <- rep(1 / K, K)
  mu <- diag(rep(1 / 5, K)) + runif(K * K) * 0.01
  sbm <- rdcsbm(N, pi, mu, rep(15, N), rep(15, N))
  x <- tril(sbm$x) + t(tril(sbm$x))
  diag(x) <- 0

  sol <- greed(x, alg = Seed())
  expect_gte(sol@K, K - 2)
  expect_lte(sol@K, K + 2)
  solc <- cut(sol, 4)
  expect_true(is.ggplot(plot(solc, type = "tree")))
  expect_true(is.ggplot(plot(solc, type = "path")))
  expect_true(is.ggplot(plot(solc, type = "front")))
  expect_true(is.ggplot(plot(solc, type = "blocks")))
  expect_true(is.ggplot(plot(solc, type = "nodelink")))
})

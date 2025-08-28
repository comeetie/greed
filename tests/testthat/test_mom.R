context("MM test")
library(greed)
library(ggplot2)
set.seed(1234)

test_that("MM sim", {
  N <- 200
  K <- 10
  pi <- rep(1 / K, K)
  mu <- cbind(diag(rep(5, K)), matrix(0, K, 20)) + matrix(runif(K * (20 + K)), K, 20 + K)
  mm <- rmm(N, pi, mu, 15)
  expect_equal(dim(mm$x)[1], N)
  expect_equal(dim(mm$x)[2], 20 + K)
  expect_equal(length(mm$cl), N)
  expect_gte(min(mm$cl), 1)
  expect_lte(max(mm$cl), K)
  model <- MoM()
  data <- greed:::preprocess(model, mm$x)
  i <- sample(200, 1)
  oldcl <- mm$cl[i]
  newcl <- sample(setdiff(1:K, oldcl), 1)
  expect_lte(greed:::test_swap(model, data, mm$cl, i, newcl), 10^-6)
  expect_lte(greed:::test_merge(model, data, mm$cl, oldcl, newcl), 10^-6)
  expect_lte(max(abs(greed:::test_merge_correction(model, data, mm$cl, oldcl, newcl))), 10^-6)
})


test_that("MM hybrid", {
  N <- 200
  K <- 4
  pi <- rep(1 / K, K)
  mu <- cbind(diag(rep(5, K)), matrix(0, K, 20)) + matrix(runif(K * (20 + K)), K, 20 + K)
  mm <- rmm(N, pi, mu, 30)
  sol <- greed(mm$x, model = MoM())
  expect_equal(sol@K, K)
  solc <- cut(sol, 3)
  # expect_true(is.ggplot(plot(solc, type = "tree")))
  # expect_true(is.ggplot(plot(solc, type = "path")))
  # expect_true(is.ggplot(plot(solc, type = "front")))
  # expect_true(is.ggplot(plot(solc, type = "blocks")))
  co <- coef(sol)
  expect_equal(sum(co$pi), 1)
  expect_equal(length(co$pi), 4)
  expect_equal(nrow(co$thetak), 4)
  expect_equal(ncol(co$thetak), 24)
})

test_that("MM seed", {
  N <- 200
  K <- 4
  pi <- rep(1 / K, K)
  mu <- cbind(diag(rep(5, K)), matrix(0, K, 20)) + matrix(runif(K * (20 + K)), K, 20 + K)
  mm <- rmm(N, pi, mu, 15)
  sol <- greed(mm$x, alg = Seed(), model = MoM())
  expect_gte(sol@K, K - 2)
  expect_lte(sol@K, K + 2)
  solc <- cut(sol, 3)
  # expect_true(is.ggplot(plot(solc, type = "tree")))
  # expect_true(is.ggplot(plot(solc, type = "path")))
  # expect_true(is.ggplot(plot(solc, type = "front")))
  # expect_true(is.ggplot(plot(solc, type = "blocks")))
})

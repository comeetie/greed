context("LCA test")
library(greed)
library(ggplot2)
set.seed(1234)



test_that("LCA icl opt", {
  N <- 500
  theta <- list(
    matrix(c(0.1, 0.9, 0.9, 0.1, 0.8, 0.2, 0.05, 0.95), ncol = 2, byrow = TRUE),
    matrix(c(0.95, 0.05, 0.3, 0.7, 0.05, 0.95, 0.05, 0.95), ncol = 2, byrow = TRUE),
    matrix(c(0.95, 0.04, 0.01, 0.9, 0.09, 0.01, 0.01, 0.01, 0.98, 0.9, 0.05, 0.05), ncol = 3, byrow = TRUE),
    matrix(c(1, 0, 0, 1, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  )
  lca.data <- rlca(N, rep(1 / 4, 4), theta)
  model <- Lca()
  data <- greed:::preprocess(model, lca.data$x)
  cl <- lca.data$cl
  i <- sample(500, 1)
  newcl <- sample(setdiff(1:lca.data$K, cl[i]), 1)
  expect_lte(greed:::test_swap(model, data, cl, i, newcl), 10^-6)
  expect_lte(greed:::test_merge(model, data, cl, 1, 4), 10^-6)
  expect_lte(max(abs(greed:::test_merge_correction(model, data, cl, 1, 4))), 10^-6)
})


test_that("LCA icl opt", {
  N <- 200
  theta <- list(
    matrix(c(0.1, 0.9, 0.9, 0.1, 0.8, 0.2, 0.05, 0.95), ncol = 2, byrow = TRUE),
    matrix(c(0.95, 0.05, 0.3, 0.7, 0.05, 0.95, 0.05, 0.95), ncol = 2, byrow = TRUE),
    matrix(c(0.95, 0.04, 0.01, 0.9, 0.09, 0.01, 0.01, 0.01, 0.98, 0.9, 0.05, 0.05), ncol = 3, byrow = TRUE),
    matrix(c(0.9, 0.1, 0.1, 0.9, 0.9, 0.1, 0.1, 0.9), ncol = 2, byrow = TRUE)
  )
  lca.data <- rlca(N, rep(1 / 4, 4), theta)
  sol <- greed(lca.data$x, model = Lca())
  expect_gte(sol@K, 3)
  co <- coef(sol)
  expect_equal(nrow(co$Thetak[[1]]), sol@K)
  # expect_true(is.ggplot(plot(sol, type = "tree")))
  # expect_true(is.ggplot(plot(sol, type = "path")))
  # expect_true(is.ggplot(plot(sol, type = "front")))
  # expect_true(methods::is(plot(sol, type = "marginals"), "gtable"))
})

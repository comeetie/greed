
context("MOR test")
library(greed)
library(ggplot2)
set.seed(1234)


test_that("MoR simulation", {
  # test with K != m (4!=3)
  N=500
  K=4
  pi=rep(1 / K, K)
  A = cbind(c(5, 1, -125, -10), c(1, 20, 1, 120), c(15, 100, 200, 1e3))
  sigma = 1
  X = cbind(rep(1, N), matrix(stats::rnorm(N * (ncol(A) - 1)), N, ncol(A) - 1))
  
  regs <- rmreg(N=N, pi=pi, A = A, sigma=sigma, X=X)
  expect_equal(length(regs$y), N)
  
  regs_defaultX <- rmreg(N=N, pi=pi, A = A, sigma=sigma)
  expect_equal(length(regs$y), N)
})


test_that("MOR internal function", {
  regs <- rmreg(500, rep(1 / 3, 3), A = t(cbind(c(5, 1, -125), c(1, 20, 1), c(15, 100, 200))), sigma = 1)
  df <- data.frame(x1 = regs$X[, 2], x2 = regs$X[, 3], y = regs$y)
  model <- MoR(y ~ x1 + x2)
  data <- greed:::preprocess(model, df)
  sol <- greed:::fit_greed(model, data, sample(1:3, 500, replace = TRUE), type = "none")
  Kp <- sol@model@tau * t(data$X) %*% data$X * 1 / 500
  M <- solve(t(data$X) %*% data$X) %*% t(data$X) %*% data$Y
  regF <- greed:::mvlm_post_comp(data$X, data$Y, matrix(0, nrow = 3, ncol = 1), Kp, sol@model@epsilon, sol@model@N0)
  i <- sample(500, 1)
  regD <- greed:::mvlm_post_del1_comp(regF, data$X[i, ], data$Y[i, ], matrix(0, nrow = 3, ncol = 1), Kp, sol@model@epsilon, sol@model@N0)
  regDA <- greed:::mvlm_post_add1_comp(regD, data$X[i, ], data$Y[i, ], matrix(0, nrow = 3, ncol = 1), Kp, sol@model@epsilon, sol@model@N0)
  expect_lte(abs(regDA$log_evidence - regF$log_evidence), 10^-6)

  regP1 <- greed:::mvlm_post_comp(data$X[1:250, ], matrix(data$Y[1:250], 250, 1), matrix(0, nrow = 3, ncol = 1), Kp, sol@model@epsilon, sol@model@N0)
  regP2 <- greed:::mvlm_post_comp(data$X[251:500, ], matrix(data$Y[251:500], 250, 1), matrix(0, nrow = 3, ncol = 1), Kp, sol@model@epsilon, sol@model@N0)
  regFM <- greed:::mvlm_post_merge_comp(regP1, regP2, matrix(0, nrow = 3, ncol = 1), Kp, sol@model@epsilon, sol@model@N0)
  expect_lte(abs(regF$log_evidence - regFM$log_evidence), 10^-6)
})

test_that("MOR hybrid", {
  m=3
  K=3
  regs <- rmreg(500, rep(1 / 3, 3), A = cbind(c(5, 1, -125), c(1, 20, 1), c(15, 100, 200)), sigma = 1)
  df <- data.frame(x1 = regs$X[, 2], x2 = regs$X[, 3], y = regs$y)
  model <- MoR(y ~ x1 + x2)
  sol <- greed(df, model)
  expect_equal(sol@K, 3)
  expect_true(is.ggplot(plot(sol, type = "tree")))
  expect_true(is.ggplot(plot(sol, type = "path")))
  expect_true(is.ggplot(plot(sol, type = "front")))
  
  param <- coef(sol)
  expect_equal(sum(param$pi), 1, tolerance = 1e-8)
  expect_length(param$A, K(sol))
  expect_length(param$Sigmak, K(sol))
  lapply(param$A, function(mat) expect_true(all.equal(dim(mat), c(m, 1))))
  lapply(param$Sigmak, function(mat) expect_true(all.equal(dim(mat), c(1, 1)))) # scalar variance
  
  i <- sample(500, 1)
  cl <- sol@cl
  newcl <- sample(setdiff(1:3, cl[i]), 1)
  data <- greed:::preprocess(model, df)
  expect_lte(greed:::test_swap(model, data, cl, i, newcl), 10^-6)
  expect_lte(greed:::test_merge(model, data, cl, 1, 2), 10^-6)
  expect_lte(max(abs(greed:::test_merge_correction(model, data, cl, 1, 2))), 10^-6)
})



test_that("MOR seed", {
  regs <- rmreg(500, rep(1 / 3, 3), A = cbind(c(5, 1, -125), c(1, 20, 1), c(15, 100, 200)), sigma = 1)
  df <- data.frame(x1 = regs$X[, 2], x2 = regs$X[, 3], y = regs$y)
  model <- MoR(y ~ x1 + x2)
  sol <- greed(df, model, alg = Seed())
  expect_equal(sol@K, 3)
  expect_true(is.ggplot(plot(sol, type = "tree")))
  expect_true(is.ggplot(plot(sol, type = "path")))
  expect_true(is.ggplot(plot(sol, type = "front")))
  sol2 <- cut(sol, 2)
  co <- coef(sol2)
  expect_equal(nrow(co$A[[1]]), 3)
})

irls_b <- function(A, b, family = poisson, lambda = 0.1, intercept = TRUE, maxit = 25, tol = 1e-08) {
  x <- rep(0, ncol(A))
  if (intercept) {
    lambdas <- c(0, rep(lambda, ncol(A)-1))
  } else {
    lambdas <- rep(lambda, ncol(A))
  }
  for (j in 1:maxit)
  {
    eta <- A %*% x
    g <- family()$linkinv(eta)
    gprime <- family()$mu.eta(eta)
    z <- eta + (b - g) / gprime
    W <- as.vector(gprime^2 / family()$variance(g))
    xold <- x
    H = crossprod(A, W * A) + diag(lambdas)
    iH = solve(H,tol = 2 * .Machine$double.eps)
    x <- iH%*%crossprod(A, W * z)
    if (sqrt(crossprod(x - xold)) < tol) break
  }
  list(coefficients = x, iterations = j,H=H,iH=iH)
}

X <- model.matrix(hp ~ -1 + drat + cyl, data = mtcars)
y <- mtcars$hp

lambda = 0.01
scale = sqrt(1/lambda)
irls_b(X,y,poisson,lambda,FALSE,5000,10^-8)
res_greed_glm  = greed:::glm_fit(X,y,lambda,5000,10^-8)


library(rstanarm)
res_stan_glm <- stan_glm(hp ~ -1 + drat + cyl, 
                      data = mtcars, family = poisson, 
                      prior = normal(0, scale), 
                      iter=20000, algorithm="sampling",
                      seed = 12345)


pairs(res_stan_glm)
res_greed_glm$beta
res_stan_glm$stan_summary

res_greed_glm$iH
res_stan_glm$covmat





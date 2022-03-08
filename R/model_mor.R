#' @include models_classes.R fit_classes.R
NULL


#' @title Multivariate mixture of regression Prior model description class
#'
#' @description
#' An S4 class to represent a multivariate mixture of regression model.
#' The model follows [minka-linear](https://tminka.github.io/papers/minka-linear.pdf) .
#' The model corresponds to the following generative model:
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i  \sim \mathcal{M}(1,\pi)}
#' \deqn{ V_k \sim \mathcal{W}(\varepsilon^{-1},n_0)}
#' \deqn{ A_k \sim \mathcal{MN}(0,(V_k)^{-1},\tau XX^\top)}
#' \deqn{ Y_{i.}|X_{i.}, A_k, Z_{ik}=1 \sim \mathcal{N}(A_k x_{i.},V_{k}^{-1})}
#' with \eqn{\mathcal{W}(\epsilon^{-1},n_0)} the Wishart distribution and \eqn{\mathcal{MN}} the matrix-normal distribution.
#' The \code{MoR-class} must be used when fitting a simple Mixture of Regression whereas the \code{MoRPrior-class} must be used when fitting a \code{\link{CombinedModels-class}}.
#' @name MoR
NULL
#> NULL

#' @rdname MoR
#' @family DlvmModels
#' @export
setClass("MoRPrior",
  representation = list(formula = "formula", tau = "numeric", N0 = "numeric", epsilon = "matrix"),
  prototype(tau = 0.001, N0 = NaN, epsilon = as.matrix(NaN))
)

setValidity("MoRPrior", function(object) {
  if (length(object@tau) > 1) {
    return("MoR model prior misspecification, tau must be of length 1.")
  }
  if (is.na(object@tau)) {
    return("MoR model prior misspecification, tau is NA.")
  }
  if (object@tau <= 0) {
    return("MoR model prior misspecification, tau must be positive.")
  }
  if (length(object@N0) > 1) {
    return("MoR model prior misspecification, N0 must be of length 1.")
  }
  TRUE
})

#' @rdname MoR
#' @export
setClass("MoR",
  contains = c("MoRPrior", "DlvmPrior")
)


#' @rdname MoR
#' @param formula a \code{\link{formula}} that describe the linear model to use
#' @param tau Prior parameter (inverse variance) default 0.001
#' @param epsilon Covariance matrix prior parameter (default to NaN, in this case epsilon will be fixed to a diagonal variance matrix equal to 0.1 time the variance of the regression residuals with only one cluster.)
#' @param N0 Prior parameter (default to NaN, in this case N0 will be fixed equal to the number of columns of Y.)
#' @return a \code{MoRPrior-class} object
#' @seealso \code{\link{MoRFit-class}}, \code{\link{MoRPath-class}}
#' @examples
#' MoRPrior(y ~ x1 + x2)
#' MoRPrior(y ~ x1 + x2, N0 = 100)
#' MoRPrior(cbind(y1, y2) ~ x1 + x2, N0 = 100)
#' @export
MoRPrior <- function(formula, tau = 0.001, N0 = NaN, epsilon = as.matrix(NaN)) {
  methods::new("MoRPrior", formula = stats::as.formula(formula), tau = tau, N0 = N0, epsilon = epsilon)
}


#' @rdname MoR
#' @param alpha Dirichlet prior parameter over the cluster proportions (default to 1)
#' @return a \code{MoR-class} object
#' @examples
#' MoR(y ~ x1 + x2)
#' MoR(y ~ x1 + x2, N0 = 100)
#' MoR(cbind(y1, y2) ~ x1 + x2, N0 = 100)
#' @export
MoR <- function(formula, alpha = 1, tau = 0.1, N0 = NaN, epsilon = as.matrix(NaN)) {
  methods::new("MoR", formula = formula, alpha = alpha, tau = tau, N0 = N0, epsilon = epsilon)
}

#' @title Clustering with a multivariate mixture of regression model fit results class
#'
#' @description An S4 class to represent a fit of a multivariate mixture of regression model, extend \code{\link{IclFit-class}}.
#' @slot model a \code{\link{MoR-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item mvmregs: list of size $K$ with statistics for each clusters
#' }
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details depends on the training procedure)
#' @seealso \code{\link{coef,MoRFit-method}}
#' @export
setClass("MoRFit", slots = list(model = "MoR"), contains = "IclFit")


#' @title Multivariate mixture of regression model hierarchical fit results class
#'
#'
#' @description An S4 class to represent a hierarchical fit of a multivariate mixture of regression model, extend \code{\link{IclPath-class}}.
#' @slot model a \code{\link{MoR-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item mvmregs: list of size $K$ with statistics for each clusters
#' }
#' @slot path a list of size K-1 with each part of the path described by:
#' \itemize{
#' \item icl1: icl value reach with this solution for alpha=1
#' \item logalpha: log(alpha) value were this solution is better than its parent
#' \item K: number of clusters
#' \item k,l: index of the cluster that were merged at this step
#' \item merge_mat: lower triangular matrix of delta icl values
#' \item obs_stats: a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item mvregs: list of size $K$ with statistics for each clusters
#' }
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy plotting with \code{ggplot2}
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father
#' @slot train_hist  data.frame with training history information (details depends on the training procedure)
# not implemented yet: @seealso \code{\link{plot,MoRFit,missing-method}}
#' @export
setClass("MoRPath", contains = c("IclPath", "MoRFit"))





#' @title Extract mixture parameters from \code{\link{MoRFit-class}} object using MAP estimation
#'
#' @param object a \code{\link{MoRFit-class}}
#' @return a list with the mixture parameters estimates (MAP), the fields are:
#' \itemize{
#' \item \code{'pi'}: cluster proportions
#' \item \code{'A'}: cluster regression matrix
#' \item \code{'Sigmak'}: cluster noise co-variance matrices
#' }
#' @export
setMethod(
  f = "coef",
  signature = signature(object = "MoRFit"),
  definition = function(object) {
    sol <- object
    pi <- (sol@obs_stats$counts + sol@model@alpha - 1) / sum(sol@obs_stats$counts + sol@model@alpha - 1)
    A <- lapply(sol@obs_stats$MoR, function(reg) {
      reg$iS %*% reg$Xty / (sol@model@tau + 1)
    })
    Sigmak <- lapply(sol@obs_stats$MoR, function(reg){
      (reg$Syx + diag(rep(sol@model@N0, nrow(reg$Syx)))) / (reg$n + sol@model@N0 + 1)
    })
    list(pi = pi, A = A, Sigmak = Sigmak)
  }
)







setMethod(
  f = "preprocess",
  signature = signature("MoR"),
  definition = function(model, data) {
    if (!methods::is(data, "data.frame")) {
      stop("X must be a data.frame a numeric vector or a matrix.", call. = FALSE)
    }

    mod_frame <- stats::model.frame(model@formula, data)
    X <- stats::model.matrix(model@formula, mod_frame)
    Y <- stats::model.response(mod_frame)

    if (prod(dim(model@epsilon)) != 1 | !all(is.nan(model@epsilon))) {
      if (dim(model@epsilon)[1] != ncol(Y) || dim(model@epsilon)[2] != ncol(Y)) {
        stop("Model prior misspecification, the dimensions of epsilon are not compatible with the data.", call. = FALSE)
      }
    }
    list(Y = as.matrix(Y), X = as.matrix(X), N = nrow(X), x_var_names = colnames(X), y_var_names = all.vars(model@formula[[2]]))
  }
)





setMethod(
  f = "seed",
  signature = signature("MoRPrior", "list", "numeric"),
  definition = function(model, data, K) {
    km <- stats::kmeans(zscore(data$X), K)
    km$cluster
  }
)

setMethod(
  f = "reorder",
  signature = signature("MoRPrior", "list", "integer"),
  definition = function(model, obs_stats, order) {
    obs_stats[order]
  }
)


setMethod(
  f = "reorder",
  signature = signature("MoR", "list", "integer"),
  definition = function(model, obs_stats, order) {
    obs_stats$MoR <- obs_stats$MoR[order]
    obs_stats$counts <- obs_stats$counts[order]
    obs_stats
  }
)


setMethod(
  f = "cleanObsStats",
  signature = signature("MoRPrior", "list"),
  definition = function(model, obs_stats, data) {
    xv_names <- data$x_var_names
    yv_names <- data$y_var_names
    new_obs_stats <- lapply(obs_stats, function(clust_stats) {
      new_clust_stats <- clust_stats
      rownames(new_clust_stats$mu) <- xv_names
      colnames(new_clust_stats$S) <- xv_names
      rownames(new_clust_stats$S) <- xv_names
      colnames(new_clust_stats$iS) <- xv_names
      rownames(new_clust_stats$iS) <- xv_names
      rownames(new_clust_stats$Xty) <- xv_names
      colnames(new_clust_stats$Xty) <- yv_names
      colnames(new_clust_stats$Yty) <- yv_names
      rownames(new_clust_stats$Yty) <- yv_names
      colnames(new_clust_stats$Syx) <- yv_names
      rownames(new_clust_stats$Syx) <- yv_names

      new_clust_stats
    })
    names(new_obs_stats) <- paste0("cluster", seq_len(length(obs_stats)))
    new_obs_stats
  }
)


setMethod(
  f = "cleanObsStats",
  signature = signature("MoR", "list"),
  definition = function(model, obs_stats, data) {
    if (!is.null(obs_stats$MoR)) {
      obs_stats$MoR <- methods::callNextMethod(model, obs_stats$MoR, data)
    }
    obs_stats
  }
)

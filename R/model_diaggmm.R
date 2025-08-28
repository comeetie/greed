#' @include models_classes.R fit_classes.R
NULL


#' @title Diagonal Gaussian Mixture Model Prior description class
#'
#' @description
#' An S4 class to represent a multivariate diagonal Gaussian mixture model.
#' The model corresponds to the following generative model:
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i  \sim \mathcal{M}(1,\pi)}
#' \deqn{ \lambda_k^{(d)} \sim \mathcal{G}(\kappa,\beta)}
#' \deqn{ \mu_k^{(d)} \sim \mathcal{N}(\mu,(\tau \lambda_k)^{-1})}
#' \deqn{ X_{i.}|Z_{ik}=1 \sim \mathcal{N}(\mu_k,\lambda_{k}^{-1})}
#' with \eqn{\mathcal{G}(\kappa,\beta)} the Gamma distribution with shape parameter \eqn{\kappa} and rate parameter \eqn{\beta}.
#' These classes mainly store the prior parameters value (\eqn{\alpha,\tau,\kappa\beta,\mu}) of this generative model.
#' The \code{DiagGmm-class} must be used when fitting a simple Diagonal Gaussian Mixture Model whereas the \code{DiagGmmPrior-class} must be sued when fitting a \code{\link{CombinedModels-class}}.
#' @name DiagGmm
NULL  
#> NULL

#' @rdname DiagGmm
#' @family DlvmModels
#' @md
#' @references Bertoletti, Marco & Friel, Nial & Rastelli, Riccardo. (2014). Choosing the number of clusters in a finite mixture model using an exact Integrated Completed Likelihood criterion. METRON. 73. 10.1007/s40300-015-0064-5. #'
#' @export
setClass("DiagGmmPrior",
  representation = list(tau = "numeric", kappa = "numeric", beta = "numeric", mu = "matrix"),
  prototype(tau = 0.01, kappa = 1, beta = NaN, mu = as.matrix(NaN))
)


setValidity("DiagGmmPrior", function(object) {
  if (length(object@tau) > 1) {
    return("DiagGmm model prior misspecification, tau must be of length 1.")
  }
  if (is.na(object@tau)) {
    return("DiagGmm model prior misspecification, tau is NA.")
  }
  if (object@tau <= 0) {
    return("DiagGmm model prior misspecification, tau must be positive.")
  }

  if (length(object@kappa) > 1) {
    return("DiagGmm model prior misspecification, kappa must be of length 1.")
  }
  if (is.na(object@kappa)) {
    return("DiagGmm model prior misspecification, kappa is NA.")
  }
  if (object@kappa <= 0) {
    return("DiagGmm model prior misspecification, kappa must be positive.")
  }

  if (length(object@beta) > 1) {
    return("DiagGmm model prior misspecification, beta must be of length 1.")
  }
  if (!is.nan(object@beta) && object@beta <= 0) {
    return("DiagGmm model prior misspecification, beta must be positive.")
  }
  TRUE
})

#' @rdname DiagGmm
#' @export
setClass("DiagGmm",
  contains = c("DlvmPrior", "DiagGmmPrior")
)


#' @rdname DiagGmm
#' @param tau Prior parameter (inverse variance), (default 0.01)
#' @param kappa Prior parameter (gamma shape), (default to 1)
#' @param beta Prior parameter (gamma rate), (default to NaN, in this case beta will be estimated from data as 0.1 time the mean of X columns variances)
#' @param mu Prior for the means (vector of size D), (default to NaN, in this case mu will be estimated from data as the mean of X)
#' @return a \code{DiagGmmPrior-class} object
#' @seealso \code{\link{DiagGmmFit-class}}, \code{\link{DiagGmmPath-class}}
#' @examples
#' DiagGmmPrior()
#' DiagGmmPrior(tau = 0.1)
#' @export
DiagGmmPrior <- function(tau = 0.01, kappa = 1, beta = NaN, mu = NaN) {
  methods::new("DiagGmmPrior", tau = tau, kappa = kappa, beta = beta, mu = as.matrix(mu))
}


#' @rdname DiagGmm
#' @param alpha Dirichlet prior parameter over the cluster proportions (default to 1)
#' @return a \code{DiagGmm-class} object
#' @export
#' @examples
#' DiagGmm()
#' DiagGmm(tau = 0.1)
#' @export
DiagGmm <- function(alpha = 1, tau = 0.01, kappa = 1, beta = NaN, mu = NaN) {
  methods::new("DiagGmm", alpha = alpha, tau = tau, kappa = kappa, beta = beta, mu = as.matrix(mu))
}


#' @title Diagonal Gaussian mixture model fit results class
#'
#' @description An S4 class to represent a fit of a multivariate diagonal Gaussian mixture model, extend \code{\link{IclFit-class}}.
#' @slot model a \code{\link{DiagGmm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item regs: list of size $K$ with statistics for each clusters
#' }
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details depends on the training procedure)
#' @seealso \code{\link{coef,DiagGmmFit-method}}
#' @export
setClass("DiagGmmFit", slots = list(model = "DiagGmm"), contains = "IclFit")


#' @title  Diagonal Gaussian mixture model hierarchical fit results class
#'
#'
#' @description An S4 class to represent a hierarchical fit of a diagonal gaussian mixture model, extend \code{\link{IclPath-class}}.
#' @slot model a \code{\link{DiagGmm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item regs: list of size $K$ with statistics for each clusters
#' }
#' @slot path a list of size K-1 with each part of the path described by:
#' \itemize{
#' \item icl1: icl value reach with this solution for alpha=1
#' \item logalpha: log(alpha) value were this solution is better than its parent
#' \item K: number of clusters
#' \item cl: vector of cluster indexes
#' \item k,l: index of the cluster that were merged at this step
#' \item merge_mat: lower triangular matrix of delta icl values
#' \item obs_stats: a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item regs: list of size $K$ with statistics for each clusters
#' }
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy plotting with \code{ggplot2}
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father
#' @slot train_hist  data.frame with training history information (details depends on the training procedure)
#' @seealso \code{\link{plot,DiagGmmFit,missing-method}}
#' @export
setClass("DiagGmmPath", contains = c("IclPath", "DiagGmmFit"))

#' @title Plot a \code{\link{DiagGmmFit-class}} object
#'
#'
#' @param x a \code{\link{DiagGmmFit-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'marginals'}: plot the marginal densities
#' \item \code{'violins'}: make a violin plot for each clusters and features
#' }
#' @return a ggplot graphic
#' @export
setMethod(
  f = "plot",
  signature = signature("DiagGmmFit", "missing"),
  definition = function(x, type = "marginals") {
    switch(type,
      marginals = {
        gg <- block_gmm_marginals(x)
        grid::grid.newpage()
        gpl <- grid::grid.draw(gg)
        invisible(gg)
      },
      violins = {
        gg <- block_gmm_marginals_violin(x)
        grid::grid.newpage()
        gpl <- grid::grid.draw(gg)
        invisible(gg)
      }
    )
  }
)


#' @title Extract mixture parameters from \code{\link{DiagGmmFit-class}} object
#'
#' @param object a \code{\link{DiagGmmFit-class}}
#' @return a list with the mixture parameters estimates (MAP), the fields are:
#' \itemize{
#' \item \code{'pi'}: cluster proportions
#' \item \code{'muk'}: cluster means
#' \item \code{'Sigmak'}: cluster co-variance matrices
#' }
#' @export
setMethod(
  f = "coef",
  signature = signature(object = "DiagGmmFit"),
  definition = function(object) {
    sol <- object
    pi <- (sol@obs_stats$counts + sol@model@alpha - 1) / sum(sol@obs_stats$counts + sol@model@alpha - 1)
    muk <- lapply(sol@obs_stats$DiagGmm, function(r) {
      (sol@model@tau * sol@model@mu + r$ng * r$m) / (sol@model@tau + r$ng)
    })
    Sigmak <- lapply(sol@obs_stats$DiagGmm, function(r) {
      betan <- sol@model@beta + 0.5 * r$S + (sol@model@tau * r$ng * (r$m - sol@model@mu)^2) / (2 * sol@model@tau + r$ng)
      alphan <- sol@model@kappa + r$ng / 2
      dd <- as.vector(betan / (alphan - 1))
      if (length(dd) > 1) {
        mode <- diag(dd)
      } else {
        mode <- as.matrix(dd)
      }
    })
    list(pi = pi, muk = muk, Sigmak = Sigmak)
  }
)


setMethod(
  f = "seed",
  signature = signature("DiagGmm", "list", "numeric"),
  definition = function(model, data, K) {
    km <- stats::kmeans(zscore(data$X), K)
    km$cluster
  }
)



setMethod(
  f = "preprocess",
  signature = signature("DiagGmmPrior"),
  definition = function(model, data) {
    if (methods::is(data, "matrix") | methods::is(data, "data.frame") | methods::is(data, "dgCMatrix")) {
      X <- as.matrix(data)
    } else {
      stop(paste0("Unsupported data type: ", class(X), " use a data.frame, a matrix, a sparse dgCMatrix."), call. = FALSE)
    }

    if (!all(is.nan(model@mu)) && length(model@mu) != ncol(X)) {
      stop("Model prior misspecification, mu length is incompatible with the data.", call. = FALSE)
    }

    list(X = X, N = nrow(X), var_names = colnames(data.frame(data)))
  }
)



setMethod(
  f = "reorder",
  signature = signature("DiagGmmPrior", "list", "integer"),
  definition = function(model, obs_stats, order) {
    obs_stats[order]
  }
)


setMethod(
  f = "reorder",
  signature = signature("DiagGmm", "list", "integer"),
  definition = function(model, obs_stats, order) {
    obs_stats$DiagGmm <- obs_stats$DiagGmm[order]
    obs_stats$counts <- obs_stats$counts[order]
    obs_stats
  }
)



setMethod(
  f = "cleanObsStats",
  signature = signature("DiagGmmPrior", "list"),
  definition = function(model, obs_stats, data) {
    num_names <- data$var_names
    new_obs_stats <- lapply(obs_stats, function(clust_stats) {
      new_clust_stats <- clust_stats[c("m", "S", "ng", "log_evidence")]
      colnames(new_clust_stats$m) <- num_names
      colnames(new_clust_stats$S) <- num_names
      new_clust_stats
    })
    names(new_obs_stats) <- paste0("cluster", seq_len(length(obs_stats)))
    new_obs_stats
  }
)

setMethod(
  f = "cleanObsStats",
  signature = signature("DiagGmm", "list"),
  definition = function(model, obs_stats, data) {
    if (!is.null(obs_stats$DiagGmm)) {
      obs_stats$DiagGmm <- methods::callNextMethod(model, obs_stats$DiagGmm, data)
    }
    obs_stats
  }
)

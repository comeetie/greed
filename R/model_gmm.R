#' @include models_classes.R fit_classes.R
NULL


#' @title Gaussian Mixture Model Prior description class
#'
#' @description
#' An S4 class to represent a multivariate Gaussian mixture model.
#' The model corresponds to the following generative model:
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i  \sim \mathcal{M}(1,\pi)}
#' \deqn{ V_k \sim \mathcal{W}(\varepsilon^{-1},n_0)}
#' \deqn{ \mu_k \sim \mathcal{N}(\mu,(\tau V_k)^{-1})}
#' \deqn{ X_{i}|Z_{ik}=1 \sim \mathcal{N}(\mu_k,V_{k}^{-1})}
#' with \eqn{\mathcal{W}(\varepsilon^{-1},n_0)} the Wishart distribution.
#' The \code{Gmm-class} must be used when fitting a simple Gaussian Mixture Model whereas the \code{GmmPrior-class} must be used when fitting a \code{\link{CombinedModels-class}}.
#' @name Gmm
NULL
#> NULL

#' @rdname Gmm
#' @family DlvmModels
#' @md
#' @references Bertoletti, Marco & Friel, Nial & Rastelli, Riccardo. (2014). Choosing the number of clusters in a finite mixture model using an exact Integrated Completed Likelihood criterion. METRON. 73. 10.1007/s40300-015-0064-5. 
#' @export
setClass("GmmPrior",
  representation = list(tau = "numeric", mu = "matrix", epsilon = "matrix", N0 = "numeric"),
  prototype(tau = 0.01, N0 = NaN, mu = as.matrix(NaN), epsilon = as.matrix(NaN))
)

setValidity("GmmPrior", function(object) {
  if (length(object@tau) > 1) {
    return("GMM model prior misspecification, tau must be of length 1.")
  }
  if (is.na(object@tau)) {
    return("GMM model prior misspecification, tau is NA.")
  }
  if (object@tau <= 0) {
    return("GMM model prior misspecification, tau must be positive.")
  }
  if (length(object@N0) > 1) {
    return("GMM model prior misspecification, N0 must be of length 1.")
  }
  TRUE
})

#' @rdname Gmm
#' @export
setClass("Gmm",
  contains = c("GmmPrior", "DlvmPrior")
)

#' @rdname Gmm
#' @param tau Prior parameter (inverse variance) default 0.01
#' @param N0 Prior parameter (pseudo count) should be > number of features (default to NaN, in this case it will be estimated from data as the number of columns of X)
#' @param epsilon Prior parameter co-variance matrix prior (matrix of size D x D), (default to a matrix of NaN, in this case epsilon will be estimated from data and will corresponds to 0.1 times a diagonal matrix with the variances of the X columns)
#' @param mu Prior parameters for the means (vector of size D), (default to NaN, in this case mu will be estimated from the data and will be equal to the mean of X)
#' @return a \code{GmmPrior-class} object
#' @seealso \code{\link{GmmFit-class}}, \code{\link{GmmPath-class}}
#' @examples
#' GmmPrior()
#' GmmPrior(tau = 0.1)
#' @export
GmmPrior <- function(tau = 0.01, N0 = NaN, mu = NaN, epsilon = NaN) {
  methods::new("GmmPrior", tau = tau, N0 = N0, mu = as.matrix(mu), epsilon = as.matrix(epsilon))
}

#' @rdname Gmm
#' @param alpha Dirichlet prior parameter over the cluster proportions (default to 1)
#' @return a \code{Gmm-class} object
#' @examples
#' Gmm()
#' Gmm(tau = 0.1, alpha = 0.5)
#' @export
Gmm <- function(tau = 0.01, N0 = NaN, mu = NaN, epsilon = NaN, alpha = 1) {
  methods::new("Gmm", alpha = alpha, tau = tau, N0 = N0, mu = as.matrix(mu), epsilon = as.matrix(epsilon))
}


#' @title Gaussian mixture model fit results class
#'
#' @description An S4 class to represent a fit of a multivariate mixture of regression model, extend \code{\link{IclFit-class}}.
#' @slot model a \code{\link{GmmPrior-class}} object to store the model fitted
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
#' @seealso \code{\link{coef,GmmFit-method}}
#' @export
setClass("GmmFit", slots = list(model = "Gmm"), contains = "IclFit")


#' @title  Gaussian mixture model hierarchical fit results class
#'
#'
#' @description An S4 class to represent a hierarchical fit of a gaussian mixture model, extend \code{\link{IclPath-class}}.
#' @slot model a \code{\link{GmmPrior-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item gmm: list of size $K$ with statistics for each clusters
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
#' \item gmm: list of size $K$ with statistics for each clusters
#' }
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy plotting with \code{ggplot2}
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father
#' @slot train_hist  data.frame with training history information (details depends on the training procedure)
#' @seealso \code{\link{plot,GmmFit,missing-method}}
#' @export
setClass("GmmPath", contains = c("IclPath", "GmmFit"))

#' @title Plot a \code{\link{GmmFit-class}} object
#'
#'
#' @param x a \code{\link{GmmFit-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'marginals'}: plot the marginal densities
#' \item \code{'violins'}: make a violin plot for each clusters and features
#' }
#' @return a ggplot graphic
#' @export
setMethod(
  f = "plot",
  signature = signature("GmmFit", "missing"),
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
      },
      stop(paste0("No plot available with type :", type), call. = FALSE)
    )
  }
)



#' @title Extract mixture parameters from \code{\link{GmmFit-class}} object
#'
#' @param object a \code{\link{GmmFit-class}}
#' @return a list with the mixture parameters estimates (MAP), the fields are:
#' \itemize{
#' \item \code{'pi'}: cluster proportions
#' \item \code{'muk'}: cluster means
#' \item \code{'Sigmak'}: cluster co-variance matrices
#' }
#' @export
setMethod(
  f = "coef",
  signature = signature(object = "GmmFit"),
  definition = function(object) {
    sol <- object
    pi <- (sol@obs_stats$counts + sol@model@alpha - 1) / sum(sol@obs_stats$counts + sol@model@alpha - 1)
    muk <- lapply(sol@obs_stats$Gmm, function(r) {
      (sol@model@tau * sol@model@mu + r$ng * r$m) / (sol@model@tau + r$ng)
    })
    Sigmak <- lapply(sol@obs_stats$Gmm, function(r) {
      mu <- (sol@model@tau * sol@model@mu + r$ng * r$m) / (sol@model@tau + r$ng)
      Sc <- r$S + t(r$m) %*% r$m - t(mu) %*% mu
      S <- (Sc + sol@model@tau * t(mu - sol@model@mu) %*% (mu - sol@model@mu) + sol@model@epsilon) / (r$ng + sol@model@N0 - length(mu))
      S
    })
    list(pi = pi, muk = muk, Sigmak = Sigmak)
  }
)





setMethod(
  f = "seed",
  signature = signature("Gmm", "list", "numeric"),
  definition = function(model, data, K) {
    km <- stats::kmeans(zscore(data$X), K)
    km$cluster
  }
)



setMethod(
  f = "preprocess",
  signature = signature("GmmPrior"),
  definition = function(model, data) {
    methods::validObject(model)
    if (methods::is(data, "matrix") | methods::is(data, "data.frame") | methods::is(data, "dgCMatrix")) {
      X <- as.matrix(data)
    } else {
      stop(paste0("Unsupported data type: ", class(X), " use a data.frame, a matrix, a sparse dgCMatrix."), call. = FALSE)
    }

    if (!is.na(model@N0) & model@N0 < ncol(X)) {
      stop("Model prior misspecification, N0 must be > ncol(X).", call. = FALSE)
    }

    if (prod(dim(model@epsilon)) != 1 | !all(is.nan(model@epsilon))) {
      if (dim(model@epsilon)[1] != ncol(X) || dim(model@epsilon)[2] != ncol(X)) {
        stop("Model prior misspecification, the dimensions of epsilon are not compatible with the data.", call. = FALSE)
      }
    }
    if (!all(is.nan(model@mu)) && length(model@mu) != ncol(X)) {
      stop("Model prior misspecification, mu length is incompatible with the data.", call. = FALSE)
    }


    list(X = X, N = nrow(X), var_names = colnames(data.frame(data)))
  }
)




setMethod(
  f = "reorder",
  signature = signature("GmmPrior", "list", "integer"),
  definition = function(model, obs_stats, order) {
    obs_stats[order]
  }
)


setMethod(
  f = "reorder",
  signature = signature("Gmm", "list", "integer"),
  definition = function(model, obs_stats, order) {
    obs_stats$Gmm <- obs_stats$Gmm[order]
    obs_stats$counts <- obs_stats$counts[order]
    obs_stats
  }
)



setMethod(
  f = "cleanObsStats",
  signature = signature("GmmPrior", "list"),
  definition = function(model, obs_stats, data) {
    num_names <- data$var_names
    new_obs_stats <- lapply(obs_stats, function(clust_stats) {
      new_clust_stats <- clust_stats[c("m", "S", "ng", "log_evidence")]
      colnames(new_clust_stats$m) <- num_names
      colnames(new_clust_stats$S) <- num_names
      rownames(new_clust_stats$S) <- num_names
      new_clust_stats
    })
    names(new_obs_stats) <- paste0("cluster", seq_len(length(obs_stats)))
    new_obs_stats
  }
)

setMethod(
  f = "cleanObsStats",
  signature = signature("Gmm", "list"),
  definition = function(model, obs_stats, data) {
    if (!is.null(obs_stats$Gmm)) {
      obs_stats$Gmm <- methods::callNextMethod(model, obs_stats$Gmm, data)
    }
    obs_stats
  }
)

#' @include models_classes.R fit_classes.R
NULL

#' @title Latent Class Analysis Model Prior class
#'
#' @description
#' An S4 class to represent a Latent Class Analysis model
#' Such model can be used to cluster a data.frame \eqn{X} with several columns of factors with the following generative model :
#' \deqn{\pi \sim \textrm{Dirichlet}(\alpha),}
#' \deqn{\forall k, \forall j, \quad \theta_{kj} \sim \textrm{Dirichlet}_{d_j}(\beta),}
#' \deqn{Z_i \sim \mathcal{M}_K(1,\pi),}
#' \deqn{\forall j=1, \ldots, p, \quad X_{ij}|Z_{ik}=1 \sim \mathcal{M}_{d_j}(1, \theta_{kj}),}
#' These classes mainly store the prior parameters value (\eqn{\alpha,\beta}) of this generative model.
#' The \code{Lca-class} must be used when fitting a simple Latent Class Analysis whereas the \code{LcaPrior-class} must be used when fitting a \code{\link{CombinedModels-class}}.
#' @name Lca
NULL
#> NULL

#' @rdname Lca
#' @family DlvmModels
#' @export
setClass("LcaPrior",
  representation = list(beta = "numeric"),
  prototype(beta = 1)
)

setValidity("LcaPrior", function(object) {
  if (length(object@beta) > 1) {
    return("Lca model prior misspecification, beta must be of length 1.")
  }
  if (is.na(object@beta)) {
    return("Lca model prior misspecification, beta is NA.")
  }
  if (object@beta <= 0) {
    return("Lca model prior misspecification, beta must be positive.")
  }
  TRUE
})

#' @rdname Lca
#' @export
setClass("Lca",
  contains = c("DlvmPrior", "LcaPrior")
)


#' @rdname Lca
#' @param beta Dirichlet prior parameter for all the categorical feature (default to 1)
#' @return a \code{LcaPrior-class} object
#' @seealso \code{\link{LcaFit-class}}, \code{\link{LcaPath-class}}
#' @examples
#' LcaPrior()
#' LcaPrior(beta = 0.5)
#' @export
LcaPrior <- function(beta = 1) {
  methods::new("LcaPrior", beta = 1)
}


#' @rdname Lca
#' @param alpha Dirichlet prior parameter over the cluster proportions (default to 1)
#' @return a \code{Lca-class} object
#' @examples
#' Lca()
#' Lca(beta = 0.5)
#' @export
Lca <- function(alpha = 1, beta = 1) {
  methods::new("Lca", alpha = alpha, beta = beta)
}

#' @title Latent Class  Analysis fit results class
#'
#' @description An S4 class to represent a fit of a Latent Class Analysis model
#'   for categorical data clustering, extend \code{\link{IclFit-class}}. The
#'   original data must be an n x p matrix where p is the number of variables
#'   and each variable is encoded as a factor (integer-valued).
#'
#' @slot model a \code{\link{Lca-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with cluster indexes
#' @slot obs_stats a list with the following elements: \itemize{ \item counts:
#'   numeric vector of size K with number of elements in each clusters \item
#'   x_counts: matrix of size K*D with the number of occurrences of each
#'   modality for each clusters }
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details
#'   depends on the training procedure)
#' @seealso \code{\link{coef,LcaFit-method}}
#' @export
setClass("LcaFit", slots = list(model = "Lca"), contains = "IclFit")


#' @title Latent Class Analysis hierarchical fit results class
#'
#'
#' @description An S4 class to represent a fit of a Latent Class Analysis model, extend \code{\link{IclPath-class}}.
#' @slot model a \code{\link{Lca-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size K*D with the number of occurrence of modality word in each clusters
#' }
#' @slot path a list of size K-1 with each part of the path described by:
#' \itemize{
#' \item icl1: icl value reach with this solution for alpha=1
#' \item logalpha: log(alpha) value were this solution is better than its parent
#' \item K: number of clusters
#' \item cl: vector of cluster indexes
#' \item k,l: index of the cluster that were merged at this step
#' \item merge_mat: lower triangular matrix of delta icl values
#' \item obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size K*D with the number of occurrence of modality word in each clusters
#' }
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy plotting with \code{ggplot2}
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father
#' @slot train_hist  data.frame with training history information (details depends on the training procedure)
#' @seealso \code{\link{plot,LcaFit,missing-method}}
#' @export
setClass("LcaPath", contains = c("IclPath", "LcaFit"))


#' @title Plot a \code{\link{LcaFit-class}} object
#'
#'
#' @param x a \code{\link{LcaFit-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing connections between clusters
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export
setMethod(
  f = "plot",
  signature = signature("LcaFit", "missing"),
  definition = function(x, type = "marginals") {
    pl <- greed:::block_lca(x)
    grid::grid.newpage()
    gpl <- grid::grid.draw(pl)
    invisible(pl)
  }
)
#' @title Extract parameters from an \code{\link{LcaFit-class}} object
#'
#' @param object a \code{\link{LcaFit-class}}
#' @return a list with the model parameters estimates (MAP), the fields are:
#' \itemize{
#' \item \code{'pi'}: cluster proportions
#' \item \code{'thetav'}: cluster profile probabilities (list of matrix of size K x Dv),
#' }
#' @export
setMethod(
  f = "coef",
  signature = signature(object = "LcaFit"),
  definition = function(object) {
    sol <- object
    pi <- (sol@obs_stats$counts + sol@model@alpha - 1) / sum(sol@obs_stats$counts + sol@model@alpha - 1)
    Thetak <- lapply(sol@obs_stats$Lca$x_counts, function(mat) {
      (mat + sol@model@beta - 1) / rowSums(mat + sol@model@beta - 1)
    })
    list(pi = pi, Thetak = Thetak)
  }
)




setMethod(
  f = "preprocess",
  signature = signature("LcaPrior"),
  definition = function(model, data) {
    if (!methods::is(data, "data.frame")) {
      stop("An Lca model expect a data.frame.", call. = FALSE)
    }
    if (!all(sapply(data, is.factor)) & !all(sapply(data, is.character))) {
      stop("An Lca model expect a data.frame with only factors or characters.", call. = FALSE)
    }

    if (all(sapply(data, is.character))) {
      # Convert characters to factors
      data <- data.frame(lapply(data, factor))
    }
    data <- droplevels(data)
    list(X = sapply(data, unclass) - 1, N = nrow(data), var_names = colnames(data), levels_names = lapply(data, levels))
  }
)

# For now, seed is kmeans on raw table for categorical
# May be improved
setMethod(
  f = "seed",
  signature = signature("Lca", "list", "numeric"),
  definition = function(model, data, K) {
    km <- stats::kmeans(as.matrix(data$X), K)
    km$cluster
  }
)







setMethod(
  f = "reorder",
  signature = signature("Lca", "list", "integer"),
  definition = function(model, obs_stats, order) {
    obs_stats$Lca$x_counts <- lapply(obs_stats$Lca$x_counts, function(var_counts) {
      var_counts[order, ]
    })
    obs_stats$Lca$counts <- obs_stats$Lca$counts[order]
    obs_stats$counts <- obs_stats$counts[order]
    obs_stats
  }
)



setMethod(
  f = "reorder",
  signature = signature("LcaPrior", "list", "integer"),
  definition = function(model, obs_stats, order) {
    obs_stats$x_counts <- lapply(obs_stats$x_counts, function(var_counts) {
      var_counts[order, ]
    })
    obs_stats$counts <- obs_stats$Lca$counts[order]
    obs_stats
  }
)


setMethod(
  f = "cleanObsStats",
  signature = signature("LcaPrior", "list"),
  definition = function(model, obs_stats, data) {
    if (length(obs_stats$x_counts) > 0) {
      cat_names <- data$var_names
      for (v in seq_len(length(obs_stats$x_counts))) {
        nb_levels <- length(data$levels_names[[v]])
        obs_stats$x_counts[[v]] <- matrix(obs_stats$x_counts[[v]], ncol = nb_levels)
        colnames(obs_stats$x_counts[[v]]) <- data$levels_names[[v]]
        rownames(obs_stats$x_counts[[v]]) <- paste0("cluster", seq_len(nrow(obs_stats$x_counts[[v]])))
      }
      names(obs_stats$x_counts) <- cat_names
    }
    obs_stats
  }
)

setMethod(
  f = "cleanObsStats",
  signature = signature("Lca", "list"),
  definition = function(model, obs_stats, data) {
    if (!is.null(obs_stats$Lca)) {
      obs_stats$Lca <- methods::callNextMethod(model, obs_stats$Lca, data)
    }
    obs_stats
  }
)

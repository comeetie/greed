#' @include models_classes.R fit_classes.R
NULL


#' @title Multinomial Stochastic Block Model Prior class
#'
#' @description
#' An S4 class to represent a Multinomial Stochastic Block Model. Such model can be used to cluster multi-layer graph vertex, and model a square adjacency cube \eqn{X} of size NxNxM with the following generative model :
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i  \sim \mathcal{M}(1,\pi)}
#' \deqn{ \theta_{kl} \sim Dirichlet(\beta)}
#' \deqn{ X_{ij.}|Z_{ik}Z_{jl}=1 \sim \mathcal{M}(L_{ij},\theta_{kl})}
#' With \eqn{L_{ij}=\sum_{m=1}^MX_{ijm}}. These classes mainly store the prior parameters value \eqn{\alpha,\beta} of this generative model.
#' The \code{MultSbm-class} must be used when fitting a simple MultSbm whereas the \code{MultSbmPrior-class} must be sued when fitting a \code{\link{CombinedModels-class}}.
#' @name MultSbm
NULL
#> NULL

#' @rdname MultSbm
#' @family DlvmModels
#' @export
setClass("MultSbmPrior",
  representation = list(beta = "numeric", type = "character"),
  prototype(beta = 1, type = "guess")
)


setValidity("MultSbmPrior", function(object) {
  if (length(object@beta) > 1) {
    return("MultSbm model prior misspecification, beta must be of length 1.")
  }
  if (is.na(object@beta)) {
    return("MultSbm model prior misspecification, beta is NA.")
  }
  if (object@beta <= 0) {
    return("MultSbm model prior misspecification, beta must be positive.")
  }
  if (!(object@type %in% c("directed", "undirected", "guess"))) {
    return("MultSbm model prior misspecification, model type must directed, undirected or guess.")
  }
  TRUE
})

#' @rdname MultSbm
#' @export
setClass("MultSbm",
  contains = c("DlvmPrior", "MultSbmPrior")
)


#' @rdname MultSbm
#' @param beta Dirichlet prior parameter over Multinomial links
#' @param type define the type of networks (either "directed", "undirected" or "guess", default to "guess"), for undirected graphs the adjacency matrix is supposed to be symmetric.
#' @return a \code{MultSbmPrior-class} object
#' @seealso \code{\link{MultSbmFit-class}}, \code{\link{MultSbmPath-class}}
#' @examples
#' MultSbmPrior()
#' MultSbmPrior(type = "undirected")
#' @export
MultSbmPrior <- function(beta = 1, type = "guess") {
  methods::new("MultSbmPrior", beta = beta, type = type)
}


#' @rdname MultSbm
#' @param alpha Dirichlet prior parameter over the cluster proportions (default to 1)
#' @return a \code{MultSbm-class} object
#' @examples
#' MultSbm()
#' MultSbm(type = "undirected")
#' @export
MultSbm <- function(alpha = 1, beta = 1, type = "guess") {
  methods::new("MultSbm", alpha = alpha, beta = beta, type = type)
}


#' @title Multinomial Stochastic Block Model fit results class
#'
#' @description An S4 class to represent a fit of a Multinomial Stochastic Block Model, extend \code{\link{IclFit-class}}.
#' @slot model a \code{\link{MultSbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: cube of size KxKxM with the number of links between each pair of clusters
#' }
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details depends on the training procedure)
#' @seealso \code{\link{coef,MultSbmFit-method}}
#' @export
setClass("MultSbmFit", slots = list(model = "MultSbm"), contains = "IclFit")


#' @title Multinomial Stochastic Block Model hierarchical fit results class
#'
#'
#' @description An S4 class to represent a hierarchical fit of a Multinomial Stochastic Block Model, extend \code{\link{IclPath-class}}.
#' @slot model a \code{\link{MultSbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size KxKxM with the number of links between each pair of clusters
#' }
#' @slot path a list of size K-1 with each part of the path described by:
#' \itemize{
#' \item icl1: icl value reach with this solution for alpha=1
#' \item logalpha: log(alpha) value were this solution is better than its parent
#' \item K: number of clusters
#' \item cl: vector of cluster indexes
#' \item k,l: index of the cluster that were merged at this step
#' \item merge_mat: lower triangular matrix of delta icl values
#' \item obs_stats: a list with the elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size KxKxM with the number of links between each pair of clusters
#' }
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy plotting with \code{ggplot2}
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father
#' @slot train_hist  data.frame with training history information (details depends on the training procedure)
#' @seealso \code{\link{plot,MultSbmFit,missing-method}}
#' @export
setClass("MultSbmPath", contains = c("IclPath", "MultSbmFit"))

#' @title Plot a \code{\link{MultSbmFit-class}} object
#'
#'
#' @param x a \code{\link{MultSbmFit-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing connections between clusters
#' \item \code{'nodelink'}: plot a nodelink diagram of the graph summarizing connections between clusters
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export
setMethod(
  f = "plot",
  signature = signature("MultSbmFit", "missing"),
  definition = function(x, type = "blocks") {
    switch(type,
      blocks = graph_blocks_cube(x),
      nodelink = nodelink_cube(x),
      stop(paste0("No plot available with type :", type), call. = FALSE)
    )
  }
)

#' @title Extract parameters from an \code{\link{MultSbmFit-class}} object
#'
#' @param object a \code{\link{MultSbmFit-class}}
#' @return a list with the model parameters estimates (MAP), the fields are:
#' \itemize{
#' \item \code{'pi'}: cluster proportions
#' \item \code{'thetakl'}: cluster profile probabilities (array of size K x K x D),
#' }
#' @export
setMethod(
  f = "coef",
  signature = signature(object = "MultSbmFit"),
  definition = function(object) {
    sol <- object
    pi <- (sol@obs_stats$counts + sol@model@alpha - 1) / sum(sol@obs_stats$counts + sol@model@alpha - 1)
    if (sol@model@type == "undirected") {
      x_counts <- sol@obs_stats$MultSbm$x_counts
      for (k in 1:dim(x_counts)[1]) {
        for (d in 1:dim(x_counts)[3]) {
          x_counts[k, k, d] <- x_counts[k, k, d] / 2
        }
      }
    } else {
      x_counts <- sol@obs_stats$MultSbm$x_counts
    }
    thetakl <- x_counts + sol@model@beta - 1
    norm <- colSums(aperm(thetakl, c(3, 1, 2)), 2)
    for (d in 1:dim(thetakl)[3]) {
      thetakl[, , d] <- thetakl[, , d] / norm
    }
    list(pi = pi, thetakl = thetakl)
  }
)



setMethod(
  f = "seed",
  signature = signature("MultSbm", "list", "numeric"),
  definition = function(model, data, K) {
    # pas terrible a réflechir deplier a droite sur les slices ?
    km <- stats::kmeans(data$X[, , 1], K)
    km$cluster
  }
)


setMethod(
  f = "preprocess",
  signature = signature("MultSbmPrior"),
  definition = function(model, data) {
    methods::validObject(model)
    if (!methods::is(data, "array")) {
      stop("A multsbm model expect an array with 3 dimensions.", call. = FALSE)
    }
    if (length(dim(data)) != 3) {
      stop("A multsbm model expect an array with 3 dimensions.", call. = FALSE)
    }
    if (dim(data)[1] != dim(data)[2]) {
      stop("A multsbm model expect an array of 3 dimensions with as many rows as columns.", call. = FALSE)
    }
    if (!all(round(data) == data) || min(data) < 0) {
      stop("A multsbm model expect an array of 3 dimensions with as many rows as columns filled with postive integers.", call. = FALSE)
    }
    issym <- all(sapply(1:dim(data)[3], function(d) {
      isSymmetric(data[, , d])
    }))
    if (model@type == "undirected" & !issym) {
      stop("An undirected multsbm expect a symmetric array.", call. = FALSE)
    }
    selfloops <- sapply(1:dim(data)[3], function(d) {
      sum(diag(data[, , d]))
    })
    if (model@type == "undirected" & !all(selfloops == 0)) {
      for (d in 1:dim(data)[3]) {
        diag(data[, , d]) <- 0
      }
      warning("An undirected multsbm model does not allow self loops, self loops were removed from the graph.", call. = FALSE)
    }



    list(X = data, N = nrow(data))
  }
)







setMethod(
  f = "reorder",
  signature = signature("MultSbm", "list", "integer"),
  definition = function(model, obs_stats, order) {
    obs_stats$MultSbm$x_counts <- obs_stats$MultSbm$x_counts[order, order, ]
    obs_stats$counts <- obs_stats$counts[order]
    obs_stats
  }
)

setMethod(
  f = "reorder",
  signature = signature("MultSbmPrior", "list", "integer"),
  definition = function(model, obs_stats, order) {
    obs_stats$x_counts <- obs_stats$x_counts[order, order, ]
    obs_stats
  }
)



setMethod(
  f = "cleanObsStats",
  signature = signature("MultSbmPrior", "list"),
  definition = function(model, obs_stats, data) {
    obs_stats
  }
)

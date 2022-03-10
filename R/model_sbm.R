#' @include models_classes.R fit_classes.R
NULL

#' @title Stochastic Block Model Prior class
#'
#' @description
#' An S4 class to represent a Stochastic Block Model.
#' Such model can be used to cluster graph vertex, and model a square adjacency matrix \eqn{X} with the following generative model :
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i  \sim \mathcal{M}(1,\pi)}
#' \deqn{ \theta_{kl} \sim Beta(a_0,b_0)}
#' \deqn{ X_{ij}|Z_{ik}Z_{jl}=1 \sim \mathcal{B}(\theta_{kl})}
#' These classes mainly store the prior parameters value \eqn{\alpha,a_0,b_0} of this generative model.
#' The \code{Sbm-class} must be used when fitting a simple Sbm whereas the \code{SbmPrior-class} must be used when fitting a \code{\link{CombinedModels-class}}.
#' @seealso \code{\link{greed}}
#' @name Sbm
NULL
#> NULL

#' @rdname Sbm
#' @family DlvmModels
#' @examples
#' Sbm()
#' @references Nowicki, Krzysztof and Tom A B Snijders (2001). “Estimation and prediction for stochastic block structures”. In:Journal of the American statistical association 96.455, pp. 1077–1087
#' @export
setClass("SbmPrior",
  representation = list(a0 = "numeric", b0 = "numeric", type = "character"),
  prototype(a0 = 1, b0 = 1, type = "guess")
)

setValidity("SbmPrior", function(object) {
  if (length(object@a0) > 1) {
    return("SBM model prior misspecification, a0 must be of length 1.")
  }
  if (object@a0 <= 0) {
    return("SBM model prior misspecification, a0 must be positive.")
  }
  if (length(object@b0) > 1) {
    return("SBM model prior misspecification, b0 must be of length 1.")
  }
  if (object@b0 <= 0) {
    return("SBM model prior misspecification, b0 must be positive.")
  }
  if (!(object@type %in% c("directed", "undirected", "guess"))) {
    return("SBM model prior misspecification, model type must directed, undirected or guess.")
  }
  TRUE
})

#' @rdname Sbm
#' @export
setClass("Sbm",
  contains = c("DlvmPrior", "SbmPrior")
)

#' @rdname Sbm
#' @param a0 Beta prior parameter over links (default to 1)
#' @param b0 Beta prior parameter over no-links (default to 1)
#' @param type define the type of networks (either "directed", "undirected" or "guess", default to "guess"), for undirected graphs the adjacency matrix is supposed to be symmetric.
#' @return a \code{SbmPrior-class} object
#' @seealso \code{\link{SbmFit-class}},\code{\link{SbmPath-class}}
#' @examples
#' SbmPrior()
#' SbmPrior(type = "undirected")
#' @export
SbmPrior <- function(a0 = 1, b0 = 1, type = "guess") {
  methods::new("SbmPrior", a0 = a0, b0 = b0, type = type)
}



#' @rdname Sbm
#' @param alpha Dirichlet prior parameter over the cluster proportions (default to 1)
#' @return a \code{Sbm-class} object
#' @examples
#' Sbm()
#' Sbm(type = "undirected")
#' @export
Sbm <- function(alpha = 1, a0 = 1, b0 = 1, type = "guess") {
  methods::new("Sbm", alpha = alpha, a0 = a0, b0 = b0, type = type)
}

#' @title Stochastic Block Model fit results class
#'
#' @description An S4 class to represent a fit of a Stochastic Block Model, extend \code{\link{IclFit-class}}.
#' @slot model a \code{\link{Sbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over rows and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size K*K with the number of links between each pair of clusters
#' }
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details depends on the training procedure)
#' @seealso \code{\link{coef,SbmFit-method}}
#' @export
setClass("SbmFit", slots = list(model = "Sbm"), contains = "IclFit")


#' @title Stochastic Block Model hierarchical fit results class
#'
#' @description An S4 class to represent a hierarchical fit of a stochastic block model, extend \code{\link{IclPath-class}}.
#' @slot model a \code{\link{Sbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size K*K with the number of links between each pair of clusters
#' }
#' @slot path a list of size K-1 with that store all the solutions along the path. Each element is a list with the following fields:
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
#' \item x_counts: matrix of size K*K with the number of links between each pair of clusters
#' }
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy plotting with \code{ggplot2}
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father
#' @slot train_hist  data.frame with training history information (details depends on the training procedure)
#' @seealso \code{\link{plot,SbmFit,missing-method}}
#' @export
setClass("SbmPath", contains = c("IclPath", "SbmFit"))

#' @title Plot a \code{\link{SbmFit-class}} object
#'
#' @param x a \code{\link{SbmFit-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing connections between clusters
#' \item \code{'nodelink'}: plot a nodelink diagram of the graph summarizing connections between clusters
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @seealso \code{\link{plot,IclPath,missing-method}}
#' @export
setMethod(
  f = "plot",
  signature = signature("SbmFit", "missing"),
  definition = function(x, type = "blocks") {
    switch(type,
      blocks = graph_blocks(x),
      nodelink = nodelink(x),
      stop(paste0("No plot available with type :", type), call. = FALSE)
    )
  }
)


#' @title Extract parameters from an \code{\link{SbmFit-class}} object
#'
#' @param object a \code{\link{SbmFit-class}}
#' @return a list with the model parameters estimates (MAP), the fields are:
#' \itemize{
#' \item \code{'pi'}: cluster proportions
#' \item \code{'thetakl'}: between clusters connections probabilities (matrix of size K x K)
#' }
#' @export
setMethod(
  f = "coef",
  signature = signature(object = "SbmFit"),
  definition = function(object) {
    sol <- object
    pi <- (sol@obs_stats$counts + sol@model@alpha - 1) / sum(sol@obs_stats$counts + sol@model@alpha - 1)
    thetakl <- (sol@obs_stats$Sbm$x_counts + sol@model@a0 - 1) / (t(t(sol@obs_stats$counts)) %*% sol@obs_stats$counts + sol@model@a0 + sol@model@b0 - 2)
    if (sol@model@type == "undirected") {
      diag(thetakl) <- (diag(sol@obs_stats$Sbm$x_counts) / 2 + sol@model@a0 - 1) / (sol@obs_stats$counts * (sol@obs_stats$counts - 1) / 2 + sol@model@a0 + sol@model@b0 - 2)
    }
    list(pi = pi, thetakl = thetakl)
  }
)




reorder_sbm <- function(obs_stats, or) {
  obs_stats$counts <- obs_stats$counts[or]
  obs_stats$x_counts <- obs_stats$x_counts[or, or]
  obs_stats
}


setMethod(
  f = "reorder",
  signature = signature("SbmPrior", "list", "integer"),
  definition = function(model, obs_stats, order) {
    reorder_sbm(obs_stats, order)
  }
)

setMethod(
  f = "reorder",
  signature = signature("Sbm", "list", "integer"),
  definition = function(model, obs_stats, order) {
    obs_stats$Sbm <- reorder_sbm(obs_stats$Sbm, order)
    obs_stats$counts <- obs_stats$counts[order]
    obs_stats
  }
)

setMethod(
  f = "cleanObsStats",
  signature = signature("SbmPrior", "list"),
  definition = function(model, obs_stats, data) {
    obs_stats
  }
)


setMethod(
  f = "seed",
  signature = signature("Sbm", "list", "numeric"),
  definition = function(model, data, K) {
    spectral(data$X, K)
  }
)

setMethod(
  f = "preprocess",
  signature = signature("SbmPrior"),
  definition = function(model, data) {
    methods::validObject(model)

    if (!(methods::is(data, "dgCMatrix") | methods::is(data, "matrix") | methods::is(data, "data.frame") | methods::is(data,"igraph"))) {
      stop("An SBM model expect a matrix,  a sparse (dgCMatrix) matrix,  a data.frame, or an igraph/tidygraph object.", call. = FALSE)
    }
    if (requireNamespace("igraph", quietly = TRUE) & methods::is(data, "igraph")) {
      
      if(length(igraph::graph_attr_names(data))){
        warning("Sbm model used with an igraph object. Vertex and nodes attributes were removed, the clustering will only use the adjacency matrix.", call. = FALSE)
      }
      data <- igraph::as_adj(data,type = "both",names=TRUE,sparse=TRUE)
    }
    
    if (methods::is(data, "data.frame")) {
      data <- as.matrix(data)
    }
    if (nrow(data) != ncol(data)) {
      stop("An SBM model expect a square adjacency matrix.", call. = FALSE)
    }
    if (!all(round(data) == data) || min(data) != 0 || max(data) != 1) {
      stop("An SBM model expect a binary adjacency matrix.", call. = FALSE)
    }
    if (model@type == "undirected" & !isSymmetric(data)) {
      stop("An undirected sbm model expect a symmetric matrix.", call. = FALSE)
    }
    if (model@type == "undirected" & sum(diag(data)) != 0) {
      diag(data) <- 0
      warning("An undirected SBM model does not allow self loops, self loops were removed from the graph.", call. = FALSE)
    }

    list(X = as.sparse(data), N = nrow(data))
  }
)

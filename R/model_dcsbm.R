#' @include models_classes.R fit_classes.R
NULL

#' @title Degree Corrected Stochastic Block Model Prior class
#'
#' @description
#' An S4 class to represent a Degree Corrected Stochastic Block Model.
#' Such model can be used to cluster graph vertex, and model a square adjacency matrix \eqn{X} with the following generative model :
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i  \sim \mathcal{M}(1,\pi)}
#' \deqn{ \theta_{kl} \sim Exponential(p)}
#' \deqn{ \gamma_i^+,\gamma_i^- \sim \mathcal{U}(S_k)}
#' \deqn{ X_{ij}|Z_{ik}Z_{jl}=1 \sim \mathcal{P}(\gamma_i^+\theta_{kl}\gamma_j^-)}
#' The individuals parameters \eqn{\gamma_i^+,\gamma_i^-} allow to take into account the node degree heterogeneity.
#' These parameters have uniform priors over the simplex \eqn{S_k} ie. \eqn{\sum_{i:z_{ik}=1}\gamma_i^+=1}.
#' These classes mainly store the prior parameters value \eqn{\alpha,p} of this generative model.
#' The \code{DcSbm-class} must be used when fitting a simple Degree Corrected Stochastic Block Model whereas the \code{DcSbmPrior-class} must be used when fitting a \code{\link{CombinedModels-class}}.
#' @name DcSbm
NULL
#> NULL

#' @rdname DcSbm
#' @family DlvmModels
#' @export
setClass("DcSbmPrior",
  representation = list(type = "character", p = "numeric"),
  prototype(p = NaN, type = "guess")
)

setValidity("DcSbmPrior", function(object) {
  if (length(object@p) > 1) {
    return("DcSbm model prior misspecification, p must be of length 1.")
  }
  if (!is.nan(object@p) && object@p <= 0) {
    return("DcSbm model prior misspecification, p must be positive.")
  }
  if (!(object@type %in% c("directed", "undirected", "guess"))) {
    return("DcSbm model prior misspecification, model type must directed, undirected or guess.")
  }

  TRUE
})

#' @rdname DcSbm
#' @export
setClass("DcSbm",
  contains = c("DlvmPrior", "DcSbmPrior")
)


#' @rdname DcSbm
#' @param p Exponential prior parameter (default to NaN, in this case p will be estimated from data as the mean connection probability)
#' @param type define the type of networks (either "directed", "undirected" or "guess", default to "guess")
#' @return a \code{DcSbmPrior-class} object
#' @seealso \code{\link{DcSbmFit-class}}, \code{\link{DcSbmPath-class}}
#' @examples
#' DcSbmPrior()
#' DcSbmPrior(type = "undirected")
#' @export
DcSbmPrior <- function(p = NaN, type = "guess") {
  methods::new("DcSbmPrior", p = p, type = type)
}


#' @rdname DcSbm
#' @param alpha Dirichlet prior parameter over the cluster proportions (default to 1)
#' @return a \code{DcSbm-class} object
#' @examples
#' DcSbm()
#' DcSbm(type = "undirected")
#' @export
DcSbm <- function(alpha = 1, p = NaN, type = "guess") {
  methods::new("DcSbm", alpha = alpha, p = p, type = type)
}



#' @title Degree Corrected Stochastic Block Model fit results class
#'
#' @description
#'  An S4 class to represent a fit of a degree corrected stochastic block model for co_clustering, extend \code{\link{IclFit-class}}.
#' @slot model a \code{\link{DcSbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item din: numeric vector of size K which store the sums of in-degrees for each clusters
#' \item dout: numeric vector of size K which store the sums of out-degrees for each clusters
#' \item x_counts: matrix of size K*K with the number of links between each pair of clusters
#' }
#' @slot obs_stats_cst a list with the following elements:
#' \itemize{
#' \item din_node: node in-degree, a vector of size N
#' \item dout_node: node in-degree vector of size N
#' }
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details depends on the training procedure)
#' @seealso \code{\link{coef,DcSbmFit-method}}
#' @export
setClass("DcSbmFit", slots = list(model = "DcSbm"), contains = "IclFit")





#' @title Degree Corrected Stochastic Block Model hierarchical fit results class
#'
#'
#' @description An S4 class to represent a hierarchical fit of a degree corrected stochastic block model, extend \code{\link{IclPath-class}}.
#' @slot model a \code{\link{DcSbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item din: numeric vector of size K which store the sums of in-degrees for each clusters
#' \item dout: numeric vector of size K which store the sums of out-degrees for each clusters
#' \item x_counts: matrix of size K*K with the number of links between each pair of clusters
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
#' \item din: numeric vector of size K which store the sums of in-degrees for each clusters
#' \item dout: numeric vector of size K which store the sums of out-degrees for each clusters
#' \item x_counts: matrix of size K*K with the number of links between each pair of clusters
#' }
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy plotting with \code{ggplot2}
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father
#' @slot train_hist  data.frame with training history information (details depends on the training procedure)
#' @seealso \code{\link{plot,DcSbmFit,missing-method}}
#' @export
setClass("DcSbmPath", contains = c("IclPath", "DcSbmFit"))


#' @title Plot a \code{\link{DcSbmFit-class}} object
#'
#' @param x a \code{\link{DcSbmFit-class}}
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
  signature = signature("DcSbmFit", "missing"),
  definition = function(x, type = "blocks") {
    switch(type,
      blocks = graph_blocks(x),
      nodelink = nodelink(x),
      stop(paste0("No plot available with type :", type), call. = FALSE)
    )
  }
)


#' @title Extract parameters from an \code{\link{DcSbmFit-class}} object
#'
#' @param object a \code{\link{DcSbmFit-class}}
#' @return a list with the model parameters estimates (MAP), the fields are the following for "directed" models :
#' \itemize{
#' \item \code{'pi'}: cluster proportions
#' \item \code{'thetakl'}: between cluster normalized connection intensities (matrix of size K x K),
#' \item \code{gammain}: node in-degree correction parameter
#' \item \code{gammaout}: node out-degree correction parameter
#' }
#' And as follow for un-directed models :
#' #' \itemize{
#' \item \code{'pi'}: cluster proportions
#' \item \code{'thetakl'}: between cluster normalized connection intensities (matrix of size K x K),
#' \item \code{gamma}: node degree correction parameter
#' }
#' @details in case of undirected graph
#' @export
setMethod(
  f = "coef",
  signature = signature(object = "DcSbmFit"),
  definition = function(object) {
    sol <- object
    pi <- (sol@obs_stats$counts + sol@model@alpha - 1) / sum(sol@obs_stats$counts + sol@model@alpha - 1)
    if (sol@model@type == "directed") {
      thetakl <- (sol@obs_stats$DcSbm$x_counts) / (t(t(sol@obs_stats$counts)) %*% sol@obs_stats$counts + 1 / sol@model@p)
      gammain <- sol@obs_stats_cst$din_node / sol@obs_stats$DcSbm$din[sol@cl]
      gammaout <- sol@obs_stats_cst$dout_node / sol@obs_stats$DcSbm$dout[sol@cl]
      params <- list(pi = pi, thetakl = thetakl, gammain = gammain, gammaout = gammaout)
    } else {
      thetakl <- (sol@obs_stats$DcSbm$x_counts)
      diag(thetakl) <- diag(thetakl) / 2
      norm <- t(t(sol@obs_stats$counts)) %*% sol@obs_stats$counts
      diag(norm) <- diag(norm) / 2 - sol@obs_stats$counts
      thetakl <- thetakl / (norm + 1 / sol@model@p)
      gamma <- sol@obs_stats_cst$din_node / sol@obs_stats$DcSbm$din[sol@cl]
      params <- list(pi = pi, thetakl = thetakl, gamma = gamma)
    }
    params
  }
)


setMethod(
  f = "seed",
  signature = signature("DcSbm", "list", "numeric"),
  definition = function(model, data, K) {
    spectral(data$X, K)
  }
)

setMethod(
  f = "preprocess",
  signature = signature("DcSbmPrior"),
  definition = function(model, data) {
    if (!(methods::is(data, "dgCMatrix") | methods::is(data, "matrix") | methods::is(data, "data.frame") | methods::is(data, "igraph"))) {
      stop("A DcSbm model expect a data.frame, a matrix or a sparse (dgCMatrix) matrix.", call. = FALSE)
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
      stop("A DcSbm model expect a square matrix.", call. = FALSE)
    }
    if (!all(round(data) == data) || min(data) < 0) {
      stop("A DcSbm model expect an integer matrix with postive values.", call. = FALSE)
    }
    if (model@type == "undirected" & !isSymmetric(data)) {
      stop("A undirected DcSbm model expect a symmetric matrix.", call. = FALSE)
    }
    if (model@type == "undirected" & sum(diag(data)) != 0) {
      diag(data) <- 0
      warning("An undirected DcSbm model does not allow self loops, self loops were removed from the graph.", call. = FALSE)
    }
    list(X = as.sparse(data), N = nrow(data))
  }
)




setMethod(
  f = "cleanObsStats",
  signature = signature("DcSbmPrior", "list"),
  definition = function(model, obs_stats, data) {
    obs_stats
  }
)



setMethod(
  f = "reorder",
  signature = signature("DcSbmPrior", "list", "integer"),
  definition = function(model, obs_stats, order) {
    obs_stats$counts <- obs_stats$counts[order]
    obs_stats$din <- obs_stats$din[order]
    obs_stats$dout <- obs_stats$dout[order]
    obs_stats$x_counts <- obs_stats$x_counts[order, order]
    obs_stats
  }
)



setMethod(
  f = "reorder",
  signature = signature("DcSbm", "list", "integer"),
  definition = function(model, obs_stats, order) {
    obs_stats$DcSbm$counts <- obs_stats$DcSbm$counts[order]
    obs_stats$DcSbm$din <- obs_stats$DcSbm$din[order]
    obs_stats$DcSbm$dout <- obs_stats$DcSbm$dout[order]
    obs_stats$DcSbm$x_counts <- obs_stats$DcSbm$x_counts[order, order]
    obs_stats$counts <- obs_stats$counts[order]
    obs_stats
  }
)

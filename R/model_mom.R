#' @include models_classes.R fit_classes.R
NULL

#' @title Mixture of Multinomial Model Prior description class
#'
#' @description
#' An S4 class to represent a Mixture of Multinomial model.
#' Such model can be used to cluster a data matrix \eqn{X} with the following generative model :
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i  \sim \mathcal{M}(1,\pi)}
#' \deqn{ \theta_{k} \sim Dirichlet(\beta)}
#' \deqn{ X_{i.}|Z_{ik}=1 \sim \mathcal{M}(L_i,\theta_{k})}
#' With \eqn{L_i=\sum_d=1^DX_{id}}. These classes mainly store the prior parameters value (\eqn{\alpha,\beta}) of this generative model.
#' The \code{MoM-class} must be used when fitting a simple Mixture of Multinomials whereas the \code{MoMPrior-class} must be sued when fitting a \code{\link{CombinedModels-class}}.
#' @name MoM 
NULL
#> NULL

#' @rdname MoM 
#' @family DlvmModels
#' @export
setClass("MoMPrior",
  representation = list(beta = "numeric"),
  prototype(beta = 1)
)


setValidity("MoMPrior", function(object) {
  if (length(object@beta) > 1) {
    return("MoM model prior misspecification, beta must be of length 1.")
  }
  if (is.na(object@beta)) {
    return("MoM model prior misspecification, beta is NA.")
  }
  if (object@beta <= 0) {
    return("MoM model prior misspecification, beta must be positive.")
  }
  TRUE
})


#' @rdname MoM 
#' @export
setClass("MoM",
  contains = c("DlvmPrior", "MoMPrior")
)

#' @rdname MoM 
#' @param beta Dirichlet over vocabulary prior parameter (default to 1)
#' @return a \code{MoMPrior-class} object
#' @seealso \code{\link{MoMFit-class}}, \code{\link{MoMPath-class}}
#' @examples
#' MoMPrior()
#' MoMPrior(beta = 0.5)
#' @export
MoMPrior <- function(beta = 1) {
  methods::new("MoMPrior", beta = 1)
}

#' @rdname MoM 
#' @param alpha Dirichlet prior parameter over the cluster proportions (default to 1)
#' @return a \code{MoM-class} object
#' @examples
#' MoM()
#' MoM(beta = 0.5)
#' @export
MoM <- function(alpha = 1, beta = 1) {
  methods::new("MoM", alpha = alpha, beta = beta)
}



#' @title Mixture of Multinomial fit results class
#'
#' @description
#'  An S4 class to represent a fit of a degree corrected stochastic block model for co_clustering, extend \code{\link{IclFit-class}}.
#' @slot model a \code{\link{MoM-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item x_counts: matrix of size K*D with the number of occurrences of each modality for each clusters
#' }
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details depends on the training procedure)
#' @seealso \code{\link{coef,LcaFit-method}}
#' @export
setClass("MoMFit", slots = list(model = "MoM"), contains = "IclFit")


#' @title Mixture of Multinomial hierarchical fit results class
#'
#'
#' @description An S4 class to represent a fit of a stochastic block model, extend \code{\link{IclPath-class}}.
#' @slot model a \code{\link{MoM-class}} object to store the model fitted
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
setClass("MoMPath", contains = c("IclPath", "MoMFit"))


#' @title Plot a \code{\link{MoMFit-class}} object
#'
#'
#' @param x a \code{\link{MoMFit-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing connections between clusters
#' }
#' @return a ggplot graphic
#' @export
setMethod(
  f = "plot",
  signature = signature("MoMFit", "missing"),
  definition = function(x, type = "blocks") {
    mat_blocks(x)
  }
)



#' @title Extract parameters from an \code{\link{MoMFit-class}} object
#'
#' @param object a \code{\link{MoMFit-class}}
#' @return a list with the model parameters estimates (MAP), the fields are:
#' \itemize{
#' \item \code{'pi'}: cluster proportions
#' \item \code{'thetak'}: cluster profile probabilities (matrix of size K x D),
#' }
#' @export
setMethod(
  f = "coef",
  signature = signature(object = "MoMFit"),
  definition = function(object) {
    sol <- object
    pi <- (sol@obs_stats$counts + sol@model@alpha - 1) / sum(sol@obs_stats$counts + sol@model@alpha - 1)
    thetak <- (t(sol@obs_stats$MoM$x_counts) + sol@model@beta - 1)
    thetak <- as.matrix(thetak / rowSums(thetak))
    list(pi = pi, thetak = thetak)
  }
)



setMethod(
  f = "reorder",
  signature = signature("MoM", "list", "integer"),
  definition = function(model, obs_stats, order) {
    obs_stats$counts <- obs_stats$counts[order]
    obs_stats$MoM$x_counts <- obs_stats$MoM$x_counts[, order]
    obs_stats
  }
)
setMethod(
  f = "reorder",
  signature = signature("MoMPrior", "list", "integer"),
  definition = function(model, obs_stats, order) {
    obs_stats$x_counts <- obs_stats$x_counts[, order]
    obs_stats
  }
)



setMethod(
  f = "seed",
  signature = signature("MoM", "list", "numeric"),
  definition = function(model, data, K) {
    km <- stats::kmeans(as.matrix(data$X), K)
    km$cluster
  }
)



setMethod(
  f = "preprocess",
  signature = signature("MoMPrior"),
  definition = function(model, data) {
    methods::validObject(model)
    if (!(methods::is(data, "dgCMatrix") | methods::is(data, "matrix") | methods::is(data, "data.frame"))) {
      stop("A MoM model expect a data.frame, a matrix or a sparse (dgCMatrix) matrix.", call. = FALSE)
    }
    if (methods::is(data, "data.frame")) {
      data <- as.matrix(data)
    }
    if (!all(round(data) == data) || min(data) < 0) {
      stop("A MoM model expect an integer matrix with postive values.", call. = FALSE)
    }
    list(X = as.sparse(data), N = nrow(data))
  }
)

setMethod(
  f = "cleanObsStats",
  signature = signature("MoMPrior", "list"),
  definition = function(model, obs_stats, data) {
    obs_stats
  }
)

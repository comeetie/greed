#' @include models_classes.R fit_classes.R
NULL

#' @title Mixture of Multinomial Model Prior description class
#'
#' @description
#' An S4 class to represent a Multinomial model model.
#' Such model can be used to cluster a data matrix \eqn{X} with the following generative model :
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i  \sim \mathcal{M}(1,\pi)}
#' \deqn{ \theta_{k} \sim Dirichlet(\beta)}
#' \deqn{ X_{i.}|Z_{ik}=1 \sim \mathcal{M}(L_i,\theta_{k})}
#' With \eqn{L_i=\sum_d=1^DX_{id}}. This class mainly store the prior parameters value (\eqn{\alpha,\beta}) of this generative model in the following slots:
#' @slot name name of the model
#' @slot alpha Dirichlet over cluster proportions prior parameter (default to 1)
#' @slot beta Dirichlet over vocabulary prior parameter (default to 1)
#' @family DlvmModels
#' @export
setClass("MoMPrior",
  representation = list(beta = "numeric"),
  prototype(beta = 1)
)

#' @describeIn MoMPrior-class MoMPrior class constructor
#' @examples
#' MoMPrior()
#' MoMPrior(beta = 0.5)
#' @export
MoMPrior <- function(beta = 1) {
  methods::new("MoMPrior", beta = 1)
}

#' @describeIn MoMPrior-class MoM class constructor
setClass("MoM",
  contains = c("DlvmPrior", "MoMPrior")
)

#' @describeIn MoMPrior-class MoM class constructor
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
#' @export
setClass("MoMPath", contains = c("IclPath", "MoMFit"))


#' @title plot a \code{\link{MoMFit-class}} object
#'
#'
#' @param x a \code{\link{MoMFit-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing connections between clusters
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export
setMethod(
  f = "plot",
  signature = signature("MoMFit", "missing"),
  definition = function(x, type = "blocks") {
    mat_blocks(x)
  }
)
#' @title plot a \code{\link{MoMPath-class}} object
#'
#' @param x an \code{\link{MoMPath-class}} object
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing
#' \item \code{'nodelink'}: plot a nodelink diagram
#' \item \code{'front'}: plot the extracted front in the plane ICL, log(alpha)
#' \item \code{'path'}: plot the evolution of ICL with respect to K
#' \item \code{'tree'}: plot the associated dendrogram
#' }
#' @return a \code{\link{ggplot2}} graphic
#' @export
setMethod(
  f = "plot",
  signature = signature("MoMPath", "missing"),
  definition = function(x, type = "blocks") {
    switch(type,
      tree = {
        dendo(x)
      },
      path = {
        lapath(x)
      },
      front = {
        plot_front(x)
      },
      blocks = {
        methods::callNextMethod()
      },
      nodelink = {
        methods::callNextMethod()
      }
    )
  }
)


#' @title Extract parameters from an \code{\link{MoMFit-class}} object
#'
#' @param object a \code{\link{MoMFit-class}}
#' @return a list with the model parameters estimates (MAP), the fields are:
#' \itemize{
#' \item \code{'pi'}: cluster proportions
#' \item \code{'thetak'}: cluster profile probabilites (matrix of size K x D),
#' }
#' @export
setMethod(
  f = "coef",
  signature = signature(object = "MoMFit"),
  definition = function(object) {
    sol <- object
    pi <- (sol@obs_stats$counts + sol@model@alpha - 1) / sum(sol@obs_stats$counts + sol@model@alpha - 1)
    thetak <- (t(sol@obs_stats$x_counts) + sol@model@beta - 1)
    thetak <- as.matrix(thetak / rowSums(thetak))
    list(pi = pi, thetak = thetak)
  }
)

reorder_mm <- function(obs_stats, or) {
  obs_stats$counts <- obs_stats$counts[or]
  obs_stats$x_counts <- obs_stats$x_counts[, or]
  obs_stats
}


setMethod(
  f = "reorder",
  signature = signature("MoM", "list", "integer"),
  definition = function(model, obs_stats, order) {
    reorder_mm(obs_stats, order)
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
  f = "sample_cl",
  signature = signature("MoM", "list", "numeric"),
  definition = function(model, data, K) {
    sample(1:K, data$N, replace = TRUE)
  }
)

setMethod(
  f = "preprocess",
  signature = signature("MoM"),
  definition = function(model, data) {
    if (!(methods::is(data, "dgCMatrix") | methods::is(data, "matrix") | methods::is(data, "data.frame"))) {
      stop("n dcsbm model expect a data.frame, a matrix or a sparse (dgCMatrix) matrix.", call. = FALSE)
    }
    if (methods::is(data, "data.frame")) {
      data <- as.matrix(data)
    }
    if (!all(round(data) == data) || min(data) < 0) {
      stop("A MoM model expect an integer matrix with postive values.", call. = FALSE)
    }
    if (length(model@alpha) > 1) {
      stop("Model prior misspecification, alpha must be of length 1.", call. = FALSE)
    }
    if (is.na(model@alpha)) {
      stop("Model prior misspecification, alpha is NA.", call. = FALSE)
    }
    if (model@alpha <= 0) {
      stop("Model prior misspecification, alpha must be positive.", call. = FALSE)
    }

    if (length(model@beta) > 1) {
      stop("Model prior misspecification, beta must be of length 1.", call. = FALSE)
    }
    if (is.na(model@beta)) {
      stop("Model prior misspecification, beta is NA.", call. = FALSE)
    }
    if (model@beta <= 0) {
      stop("Model prior misspecification, beta must be positive.", call. = FALSE)
    }
    list(X = as.sparse(data), N = nrow(data))
  }
)


setMethod(
  f = "postprocess",
  signature = signature("MoMPath"),
  definition = function(path, data, X, Y = NULL) {
    path@obs_stats <- list(
      counts = path@obs_stats$counts,
      x_counts = path@obs_stats$MoM$x_counts
    )
    for (p in 1:length(path@path)) {
      path@path[[p]]$obs_stats <- list(
        counts = path@path[[p]]$obs_stats$counts,
        x_counts = path@path[[p]]$obs_stats$MoM$x_counts
      )
    }
    path
  }
)

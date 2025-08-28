#' @include models_classes.R fit_classes.R
NULL

#' @title Degree Corrected Latent Block Model for bipartite graph class
#'
#' @description An S4 class to represent a degree corrected stochastic block model for co_clustering of bipartite graph.
#' Such model can be used to cluster graph vertex, and model a bipartite graph adjacency matrix \eqn{X} with the following generative model :
#' \deqn{ \pi \sim Dirichlet(\alpha)}
#' \deqn{ Z_i^r  \sim \mathcal{M}(1,\pi^r)}
#' \deqn{ Z_j^c  \sim \mathcal{M}(1,\pi^c)}
#' \deqn{ \theta_{kl} \sim Exponential(p)}
#' \deqn{ \gamma_i^r\sim \mathcal{U}(S_k)}
#' \deqn{ \gamma_i^c\sim \mathcal{U}(S_l)}
#' \deqn{ X_{ij}|Z_{ik}^cZ_{jl}^r=1 \sim \mathcal{P}(\gamma_i^r\theta_{kl}\gamma_j^c)}
#' The individuals parameters \eqn{\gamma_i^r,\gamma_j^c} allow to take into account the node degree heterogeneity.
#' These parameters have uniform priors over simplex \eqn{S_k}.
#' These classes mainly store the prior parameters value \eqn{\alpha,p} of this generative model.
#' The \code{DcLbm-class} must be used when fitting a simple Diagonal Gaussian Mixture Model whereas the \code{DcLbmPrior-class} must be sued when fitting a \code{\link{CombinedModels-class}}.
#' @name DcLbm
NULL
#> NULL

#' @rdname DcLbm
#' @family DlvmModels
#' @export
setClass("DcLbmPrior",
  representation = list(p = "numeric"),
  prototype(p = NaN)
)

#' @rdname DcLbm
#' @export
setClass("DcLbm",
  contains = c("DlvmCoPrior", "DcLbmPrior")
)

#' @rdname DcLbm
#' @param p Exponential prior parameter (default to Nan, in this case p will be estimated from data as the average intensities of X)
#' @return a \code{DcLbmPrior-class}
#' @seealso \code{\link{DcLbmFit-class}}, \code{\link{DcLbmPath-class}}
#' @examples
#' DcLbmPrior()
#' DcLbmPrior(p = 0.7)
#' @export
DcLbmPrior <- function(p = NaN) {
  methods::new("DcLbmPrior", p = p)
}



#' @rdname DcLbm
#' @param alpha Dirichlet prior parameter over the cluster proportions (default to 1)
#' @return a \code{DcLbm-class} object
#' @examples
#' DcLbm()
#' DcLbm(p = 0.7)
#' @export
DcLbm <- function(alpha = 1, p = NaN) {
  methods::new("DcLbm", alpha = alpha, p = p)
}


#' @title Degree corrected Latent Block Model fit results class
#'
#' @description An S4 class to represent a fit of a degree corrected stochastic block model for co_clustering, extend \code{\link{IclFit-class}}.
#' @slot model a \code{\link{DcLbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot Krow number of extracted row clusters
#' @slot Kcol number of extracted column clusters
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item din: numeric vector of size K which store the sums of in-degrees for each clusters
#' \item dout: numeric vector of size K which store the sums of out-degrees for each clusters
#' \item x_counts: matrix of size K*K with the number of links between each pair of clusters
#' \item co_x_counts: matrix of size Krow*Kcol with the number of links between each pair of row and column cluster
#' }
#' @slot clrow a numeric vector with row cluster indexes
#' @slot clcol a numeric vector with column cluster indexes
#' @slot Nrow number of rows
#' @slot Ncol number of columns
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details depends on the training procedure)
#' @seealso \code{\link{coef,DcLbmFit-method}}
#' @export
setClass("DcLbmFit", slots = list(model = "DcLbm", clrow = "numeric", clcol = "numeric", Krow = "numeric", Kcol = "numeric", Nrow = "numeric", Ncol = "numeric"), contains = "IclFit")




#' @title Degree corrected Latent Block Model hierarchical fit results class
#'
#'
#' @description An S4 class to represent a fit of a degree corrected stochastic block model for co_clustering, extend \code{\link{IclPath-class}}.
#' @slot model a \code{\link{DcLbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot Krow number of extracted row clusters
#' @slot Kcol number of extracted column clusters
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' \itemize{
#' \item counts: numeric vector of size K with number of elements in each clusters
#' \item din: numeric vector of size K which store the sums of in-degrees for each clusters
#' \item dout: numeric vector of size K which store the sums of out-degrees for each clusters
#' \item x_counts: matrix of size K*K with the number of links between each pair of clusters
#' \item co_x_counts: matrix of size Krow*Kcol with the number of links between each pair of row and column cluster
#' }
#' @slot clrow a numeric vector with row cluster indexes
#' @slot clcol a numeric vector with column cluster indexes
#' @slot Nrow number of rows
#' @slot Ncol number of columns
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
#' \item co_x_counts: matrix of size Krow*Kcol with the number of links between each pair of row and column cluster
#' }
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy plotting with \code{ggplot2}
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father
#' @slot ggtreerow data.frame with complete merge tree of row clusters for easy plotting with \code{ggplot2}
#' @slot ggtreecol data.frame with complete merge tree of column clusters for easy plotting with \code{ggplot2}
#' @slot train_hist  data.frame with training history information (details depends on the training procedure)
#' @seealso \code{\link{plot,DcLbmPath,missing-method}}
#' @export
setClass("DcLbmPath", slots = list(ggtreerow = "data.frame", ggtreecol = "data.frame"), contains = c("IclPath", "DcLbmFit"))




#' @title Method to cut a DcLbmPath solution to a desired number of cluster
#'
#' @description This method take a \code{\link{DcLbmPath-class}} and an integer K and return the solution from the path with K clusters
#' @param x A an \code{\link{DcLbmPath-class}} solution
#' @param K Desired number of cluster
#' @return an \code{\link{IclPath-class}} object with the desired number of cluster
#' @export
setMethod(
  f = "cut",
  signature = signature("DcLbmPath"),
  definition = function(x, K) {
    if (K < x@K) {
      i <- which(sapply(x@path, function(p) {
        p$K
      }) == K)
      x@K <- K
      x@logalpha <- x@path[[i]]$logalpha
      x@icl <- x@path[[i]]$icl
      # Old version: x@cl = as.vector(x@path[[i]]$cl)
      # Compute cl with the history of fusions (k,l) at each stage
      for (p in 1:i) {
        # Get fusion (k,l)
        k <- x@path[[p]]$k
        l <- x@path[[p]]$l
        x@cl[x@cl == k] <- l
        # rescale @cl to be 1...K
        x@cl <- as.integer(factor(x@cl))
      }
      for (st in names(x@path[[i]]$obs_stats)) {
        x@obs_stats[st] <- x@path[[i]]$obs_stats[st]
      }

      x@path <- x@path[(i + 1):length(x@path)]
      x <- postprocess(x)
    } else {
      warning(paste0("This clustering has ", x@K, " clusters and you requested ", K, " clusters. Please provide a value for K smaller than ", x@K, "."), call. = FALSE)
    }
    x
  }
)


#' @title Plot a \code{\link{DcLbmPath-class}}
#'
#'
#' @param x a \code{\link{DcLbmPath-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'tree'}: plot a co-dendogram of rows and columns clusters
#' \item \code{'blocks'}: plot a block matrix with summarizing connections between row and column clusters
#' \item \code{'biplot'}: plot a block matrix with summarizing connections between row and column clusters aligned with row and clusters drendograms
#' \item \code{'nodelink'}: plot a nodelink diagram of the bipartite graph summarizing connections between row and column clusters
#' }
#' @return a ggplot graphic
#' @export
setMethod(
  f = "plot",
  signature = signature("DcLbmPath", "missing"),
  definition = function(x, type = "tree") {
    switch(type,
      tree = co_dendo(x),
      blocks = co_blocks(x),
      biplot = bi_plot(x),
      nodelink = co_nodelink(x),
      methods::callNextMethod(x, type = type)
    )
  }
)

#' @title Plot a \code{\link{DcLbmFit-class}}
#'
#'
#' @param x a \code{\link{DcLbmFit-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'blocks'}: plot a block matrix with summarizing connections between row and column clusters
#' \item \code{'nodelink'}: plot a nodelink diagram of the bipartite graph summarizing connections between row and column clusters
#' }
#' @return a ggplot graphic
#' @export
setMethod(
  f = "plot",
  signature = signature("DcLbmFit", "missing"),
  definition = function(x, type = "blocks") {
    switch(type,
           blocks = co_blocks(x),
           nodelink = co_nodelink(x)
    )
  }
)


#' @title Extract parameters from an \code{\link{DcLbmFit-class}} object
#'
#' @param object a \code{\link{DcLbmFit-class}}
#' @return a list with the model parameters estimates (MAP), the fields are:
#' \itemize{
#' \item \code{'pirows'}: row cluster proportions
#' \item \code{'picols'}: row cluster proportions
#' \item \code{'thetakl'}: between clusters connection probabilities (matrix of size Krow x Kcol),
#' \item \code{'gammarows'}: rows degree correction parameters (size Nrows),
#' \item \code{'gammacols'}: cols degree correction parameters (size Ncols),
#' }
#' @export
setMethod(
  f = "coef",
  signature = signature(object = "DcLbmFit"),
  definition = function(object) {
    sol <- object
    pirows <- (sol@obs_stats$rows_counts + sol@model@alpha - 1) / sum(sol@obs_stats$rows_counts + sol@model@alpha - 1)
    picols <- (sol@obs_stats$cols_counts + sol@model@alpha - 1) / sum(sol@obs_stats$cols_counts + sol@model@alpha - 1)
    gammarows <- sol@obs_stats_cst$drow / sol@obs_stats$DcLbm$dr[sol@clrow]
    gammacols <- sol@obs_stats_cst$dcol / sol@obs_stats$DcLbm$dc[sol@clcol]
    thetakl <- (sol@obs_stats$DcLbm$co_x_counts) / (t(t(sol@obs_stats$rows_counts)) %*% sol@obs_stats$cols_counts + 1 / sol@model@p)
    list(pirows = pirows, picols = picols, thetakl = thetakl, gammarows = gammarows, gammacols = gammacols)
  }
)

setMethod(
  f = "preprocess",
  signature = signature("DcLbm"),
  definition = function(model, data) {
    X <- as.sparse(data)
    if (!methods::is(X,"dgCMatrix")) {
      stop(paste0("Unsupported data type :", class(X), "for DcLbm model."), call. = FALSE)
    }
    ij <- which(X > 0, arr.ind = TRUE)
    if (!all(X[ij] == round(X[ij]))) {
      stop("Only integer matrix allowed as input, non integer values found.", call. = FALSE)
    }


    di <- dim(X)
    N <- sum(di)
    list(X = X, N = N, Nrows = di[1], Ncols = di[2])
  }
)

setMethod(
  f = "postprocess",
  signature = signature("DcLbmPath"),
  definition = function(path, data = NULL) {
    sol <- path
    if (!is.null(data)) {
      sol@Nrow <- data$Nrows
      sol@Ncol <- data$Ncols
    }
    clusters_type <- apply(table(sol@cl, c(rep(1, sol@Nrow), rep(2, sol@Ncol))), 1, which.max)
    clust_rows <- which(clusters_type == 1)
    clust_cols <- which(clusters_type == 2)
    icol <- (sol@Nrow + 1):length(sol@cl)
    irow <- 1:sol@Nrow

    sol@clrow <- as.numeric(factor(sol@cl[irow], levels = clust_rows))
    sol@Krow <- max(sol@clrow, na.rm = TRUE)
    sol@clcol <- as.numeric(factor(sol@cl[icol], levels = clust_cols))
    sol@Kcol <- max(sol@clcol, na.rm = TRUE)
    sol@obs_stats$DcLbm$co_x_counts <- sol@obs_stats$DcLbm$x_counts[clust_rows, clust_cols]
    sol@obs_stats$DcLbm$dr <- sol@obs_stats$DcLbm$dr[clust_rows]
    sol@obs_stats$DcLbm$dc <- sol@obs_stats$DcLbm$dc[clust_cols]
    sol@obs_stats$rows_counts <- sol@obs_stats$counts[clust_rows]
    sol@obs_stats$cols_counts <- sol@obs_stats$counts[clust_cols]
    if (!is.null(data)) {
      tree <- sol@ggtree
      root <- tree$node[tree$tree == 0]
      tree <- tree[tree$node != root, ]
      tree$tree[tree$tree == root] <- 0

      xcol <- tree[tree$node %in% clust_cols, ]$x
      coltree <- tree[tree$x >= min(xcol) & tree$x <= max(xcol), ]
      xrow <- tree[tree$node %in% clust_rows, ]$x
      rowtree <- tree[tree$x >= min(xrow) & tree$x <= max(xrow), ]
      sol@ggtreecol <- coltree
      sol@ggtreerow <- rowtree
    }
    sol
  }
)




setMethod(
  f = "sample_cl",
  signature = signature("DcLbm", "list", "numeric"),
  definition = function(model, data, K) {
    c(sample(1:floor(K / 2), data$Nrows, replace = TRUE), sample((floor(K / 2) + 1):K, data$Ncols, replace = TRUE))
  }
)


setMethod(
  f = "reorder",
  signature = signature("DcLbm", "list", "integer"),
  definition = function(model, obs_stats, order) {
    obs_stats$counts <- obs_stats$counts[order]
    obs_stats$DcLbm$counts <- obs_stats$DcLbm$counts[order]
    obs_stats$DcLbm$dr <- obs_stats$DcLbm$dr[order]
    obs_stats$DcLbm$dc <- obs_stats$DcLbm$dc[order]
    obs_stats$DcLbm$x_counts <- obs_stats$DcLbm$x_counts[order, order]
    obs_stats
  }
)



setMethod(
  f = "cleanObsStats",
  signature = signature("DcLbmPrior", "list"),
  definition = function(model, obs_stats, data) {
    obs_stats
  }
)

setMethod(
  f = "seed",
  signature = signature("DcLbm", "list", "numeric"),
  definition = function(model, data, K) {
    kmrow <- stats::kmeans(data$X, floor(K / 2))
    kmcol <- stats::kmeans(t(data$X), floor(K / 2))
    c(kmrow$cluster, kmcol$cluster + max(kmrow$cluster))
  }
)

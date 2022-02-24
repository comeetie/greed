#' @include models_classes.R fit_classes.R
NULL



#' @title Combined Models classes
#'
#' @description
#' An S4 class to represent a combined clustering models, where several models are used to model different datasets. A conditional independence assumption between the view knowing the cluster is made.
#' @details
#' The filed name in the models list must match the name of the list use to provide the datasets to cluster together.
#' @name CombinedModels
NULL
#> NULL

#' @rdname CombinedModels
#' @family DlvmModels
#' @export
setClass("CombinedModels",
  representation = list(models = "list"),
  contains = "DlvmPrior",
  prototype(alpha = 1, models = list())
)

setValidity("CombinedModels", function(object) {
  models_classes <- lapply(object@models, class)
  valid_models <- c("GmmPrior", "DiagGmmPrior", "MoRPrior", "SbmPrior", "DcSbmPrior", "MultSbmPrior", "MoMPrior", "LcaPrior")
  if (!all(models_classes %in% valid_models)) {
    return(paste0("At least one of the provided models to CombinedModels is not of the good classe, only ", valid_models, " may be used with a CombinedModels.", collapse = ", "))
  }
  TRUE
})

#' @rdname CombinedModels
#' @param models a named list of DlvmPrior's object
#' @param alpha Dirichlet prior parameter over the cluster proportions (default to 1)
#' @return a \code{CombinedModels-class} object
#' @seealso \code{\link{CombinedModelsFit-class}}, \code{\link{CombinedModelsPath-class}}
#' @examples
#' CombinedModels(models = list(continuous = GmmPrior(), discrete = LcaPrior()))
#' @export
CombinedModels <- function(models, alpha = 1) {
  methods::new("CombinedModels", alpha = alpha, models = models)
}



#' @title Combined Models fit results class
#'
#' @description
#'  An S4 class to represent a fit of a degree corrected stochastic block model for co_clustering, extend \code{\link{IclFit-class}}.
#' @slot model a \code{\link{DcSbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' @slot move_mat binary matrix which store move constraints
#' @slot train_hist data.frame with training history information (details depends on the training procedure)
#' @seealso \code{\link{extractSubModel,CombinedModelsPath,character-method}}
#' @export
setClass("CombinedModelsFit", slots = list(model = "CombinedModels"), contains = "IclFit")





#' @title Combined Models hierarchical fit results class
#'
#'
#' @description An S4 class to represent a hierarchical fit of a degree corrected stochastic block model, extend \code{\link{IclPath-class}}.
#' @slot model a \code{\link{DcSbm-class}} object to store the model fitted
#' @slot name generative model name
#' @slot icl icl value of the fitted model
#' @slot K number of extracted clusters over row and columns
#' @slot cl a numeric vector with row and columns cluster indexes
#' @slot obs_stats a list with the following elements:
#' @slot path a list of size K-1 with each part of the path described by:
#' \itemize{
#' \item icl1: icl value reach with this solution for alpha=1
#' \item logalpha: log(alpha) value were this solution is better than its parent
#' \item K: number of clusters
#' \item cl: vector of cluster indexes
#' \item k,l: index of the cluster that were merged at this step
#' \item merge_mat: lower triangular matrix of delta icl values
#' \item obs_stats: a list with the elements:
#' }
#' @slot logalpha value of log(alpha)
#' @slot ggtree data.frame with complete merge tree for easy plotting with \code{ggplot2}
#' @slot tree numeric vector with merge tree \code{tree[i]} contains the index of \code{i} father
#' @slot train_hist  data.frame with training history information (details depends on the training procedure)
#' @seealso \code{\link{extractSubModel,CombinedModelsPath,character-method}}
#' @export
setClass("CombinedModelsPath", contains = c("IclPath", "CombinedModelsFit"))





setMethod(
  f = "reorder",
  signature = signature("CombinedModels", "list", "integer"),
  definition = function(model, obs_stats, order) {
    mnames <- names(model@models)
    new_obs_stats <- lapply(mnames, function(current_name) {
      reorder(model@models[[current_name]], obs_stats[[current_name]], order)
    })
    names(new_obs_stats) <- mnames
    new_obs_stats[["counts"]] <- obs_stats$counts[order]
    new_obs_stats
  }
)


setMethod(
  f = "cleanObsStats",
  signature = signature("CombinedModels", "list"),
  definition = function(model, obs_stats, data) {
    mnames <- names(model@models)
    new_obs_stats <- lapply(mnames, function(current_name) {
      cleanObsStats(model@models[[current_name]], obs_stats[[current_name]], data[[current_name]])
    })
    names(new_obs_stats) <- mnames
    new_obs_stats[["counts"]] <- obs_stats$counts
    new_obs_stats
  }
)




setMethod(
  f = "preprocess",
  signature = signature("CombinedModels"),
  definition = function(model, data) {
    mnames <- names(model@models)
    if ("counts" %in% mnames) {
      stop("Prohibited models name counts is reserved, please use another model name.", .call = FALSE)
    }

    if (!(all(names(data) %in% mnames) & all(mnames %in% names(data)))) {
      stop("Models names do notch match datasets names, please check the model and dataset list.", .call = FALSE)
    }
    data_prep <- lapply(mnames, function(current_name) {
      greed:::preprocess(model@models[[current_name]], data[[current_name]])
    })
    names(data_prep) <- mnames
    Ns <- sapply(data_prep, function(x) {
      x$N
    })
    if (!all(Ns == Ns[1])) {
      stop("All datasets must have the same number of elements.", .call = FALSE)
    }

    data_prep$N <- Ns[1]

    data_prep
  }
)

#' @title Extract a part of a \code{\link{CombinedModelsPath-class}} object
#'
#' @param sol an \code{\link{CombinedModelsPath-class}} object
#' @param sub_model_name a string which specify the part of the model to
#'   extract. Note that the name must correspond to the one of the names used in
#'   the list of models during the origin call to \code{\link{greed}}.
#' @return a \code{\link{IclFit-class}} object of the relevant class
#' @export
setGeneric("extractSubModel", function(sol, sub_model_name) standardGeneric("extractSubModel"))


#' @describeIn extractSubModel CombinedModelsPath method
#' @export 
setMethod(
  f = "extractSubModel",
  signature = signature("CombinedModelsPath", "character"),
  definition = function(sol, sub_model_name) {
    model <- sol@model@models[[sub_model_name]]
    model_type <- gsub("Prior", "", as.character(class(model)))
    model <- methods::as(model, model_type)
    model@alpha <- sol@model@alpha
    sol <- methods::as(methods::as(sol, "IclFit"), paste0(model_type, "Fit"))
    sol@model <- model
    sol@obs_stats[[model_type]] <- sol@obs_stats[[sub_model_name]]
    sol@obs_stats <- sol@obs_stats[c("counts", model_type)]
    sol
  }
)

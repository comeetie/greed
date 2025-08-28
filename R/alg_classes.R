#' @useDynLib greed
#' @importFrom Rcpp sourceCpp
#' @importFrom future %<-%
#' @name %<-%
NULL

#' @include models_classes.R fit_classes.R alg_hybrid.R alg_genetic.R tools_misc.R
#' @import Matrix
NULL


#' @title Abstract optimization algorithm class
#'
#' @description
#' An S4 class to represent an abstract optimization algorithm.
#' @export
setClass("Alg")

#' @title Greedy algorithm with multiple start class
#'
#' @description
#' An S4 class to represent a greedy algorithm  with multiple start (extends \code{\link{Alg-class}} class).
#' @slot nb_start number of random starts (default to 10)
#' @export
setClass("Multistarts",
  contains = "Alg",
  representation =  list(nb_start = "numeric"),
  prototype(nb_start = 10)
)

#' @describeIn Multistarts-class Multistarts algorithm class constructor
#' @param nb_start number of random starts (default to 10)
#' @return a \code{Multistarts-class} object
#' @examples
#' Multistarts()
#' Multistarts(15)
#' @export
Multistarts <- function(nb_start = 10) {
  methods::new("Multistarts", nb_start = nb_start)
}

#' @title Greedy algorithm with seeded initialization
#'
#' @description
#' An S4 class to represent a greedy algorithm with initialization from spectral clustering and or k-means (extends \code{\link{Alg-class}} class ).
#' @export
setClass("Seed",
  contains = "Alg",
  representation =  list(),
  prototype()
)

#' @describeIn Seed-class Seed algorithm class constructor
#' @return a \code{Seed-class} object
#' @examples
#' Seed()
#' @export
Seed <- function() {
  methods::new("Seed")
}



#' @title Hybrid optimization algorithm
#'
#' @description
#' An S4 class to represent an hybrid genetic/greedy algorithm (extends \code{\link{Alg-class}} class).
#' @slot pop_size size of the solutions populations (default to 20)
#' @slot nb_max_gen maximal number of generation to produce (default to 10)
#' @slot prob_mutation mutation probability (default to 0.25)
#' @slot Kmax maximum number of clusters (default to 100)
#' @export
setClass("Hybrid",
  contains = "Alg",
  representation =  list(pop_size = "numeric", nb_max_gen = "numeric", prob_mutation = "numeric", Kmax = "numeric"),
  prototype(pop_size = 20, nb_max_gen = 10, prob_mutation = 0.25, Kmax = 100)
)

#' @describeIn Hybrid-class Hybrid algorithm class constructor
#' @param pop_size size of the solutions populations (default to 20)
#' @param nb_max_gen maximal number of generation to produce (default to 10)
#' @param prob_mutation mutation probability (default to 0.25)
#' @param Kmax maximum number of clusters (default to 100)
#' @return a \code{Hybrid-class} object
#' @examples
#' Hybrid()
#' Hybrid(pop_size = 100)
#' @export
Hybrid <- function(pop_size = 20, nb_max_gen = 10, prob_mutation = 0.25, Kmax = 100) {
  methods::new("Hybrid", pop_size = pop_size, nb_max_gen = nb_max_gen, prob_mutation = prob_mutation, Kmax = Kmax)
}



#' @title Genetic optimization algorithm
#'
#' @description
#' An S4 class to represent a genetic algorithm (extends \code{\link{Alg-class}} class).
#' @slot pop_size size of the solutions populations (default to 10)
#' @slot nb_max_gen maximal number of generation to produce (default to 4)
#' @slot prob_mutation probability of mutation (default to 0.25)
#' @slot sel_frac fraction of best solutions selected for crossing  (default to 0.75)
#' @export
setClass("Genetic",
  contains = "Alg",
  representation =  list(pop_size = "numeric", nb_max_gen = "numeric", prob_mutation = "numeric", sel_frac = "numeric"),
  prototype(pop_size = 100, nb_max_gen = 20, prob_mutation = 0.25, sel_frac = 0.75)
)

#' @describeIn Genetic-class Genetic algorithm class constructor
#' @param pop_size size of the solutions populations (default to 10)
#' @param nb_max_gen maximal number of generation to produce (default to 4)
#' @param prob_mutation probability of mutation (default to 0.25)
#' @param sel_frac fraction of best solutions selected for crossing  (default to 0.75)
#' @return a \code{Genetic-class} object
#' @examples
#' Genetic()
#' Genetic(pop_size = 500)
#' @export
Genetic <- function(pop_size = 100, nb_max_gen = 20, prob_mutation = 0.25, sel_frac = 0.75) {
  methods::new("Genetic", pop_size = pop_size, nb_max_gen = nb_max_gen, prob_mutation = prob_mutation, sel_frac = sel_frac)
}



#' @title Method to extract the clustering results from an \code{\link{IclFit-class}} object
#'
#' @description This method take a \code{\link{IclFit-class}} object and return an integer vector with the cluster assignments that were found.
#' @param fit an \code{IclFit} solution
#' @return an integer vector with cluster assignments. Zero indicates noise points.
#' @export
setGeneric("clustering", function(fit) standardGeneric("clustering"))

#' @describeIn clustering IclFit-class method
#' @export
setMethod(
  f = "clustering",
  signature = signature("IclFit"),
  definition = function(fit) {
    as.integer(fit@cl)
  }
)

#' @title Generic method to extract the ICL value from an \code{\link{IclFit-class}} object
#'
#' @description This method take a \code{\link{IclFit-class}} object and return its ICL score.
#' @param fit an \code{IclFit} solution
#' @return The ICL value achieved
#' @export
setGeneric("ICL", function(fit) standardGeneric("ICL"))

#' @describeIn ICL IclFit method
#' @export
setMethod(
  f = "ICL",
  signature = signature("IclFit"),
  definition = function(fit) {
    fit@icl
  }
)

#' @title Generic method to get the number of clusters from an \code{\link{IclFit-class}} object
#' @description This method take a \code{\link{IclFit-class}} object and return its ICL score.
#' @param fit an \code{IclFit} solution
#' @return The number of clusters
#' @export
setGeneric("K", function(fit) standardGeneric("K"))

#' @describeIn K IclFit method
#' @export
setMethod(
  f = "K",
  signature = signature("IclFit"),
  definition = function(fit) {
    fit@K
  }
)


#' @title Generic method to extract the prior used to fit \code{\link{IclFit-class}} object
#' @description This method take a \code{\link{IclFit-class}} object and return the prior used.
#' @param fit an \code{IclFit} solution
#' @return An S4 object describing the prior parameters
#' @export
setGeneric("prior", function(fit) standardGeneric("prior"))

#' @describeIn prior IclFit method
#' @export
setMethod(
  f = "prior",
  signature = signature("IclFit"),
  definition = function(fit) {
    fit@model
  }
)



#' @title Generic method to cut a path solution to a desired number of cluster
#'
#' @description This method take a \code{\link{IclPath-class}} object and an integer K and return the solution from the path with K clusters
#' @param x A an \code{IclPath} solution
#' @param K Desired number of cluster
#' @return an \code{\link{IclPath-class}} object with the desired number of cluster
#' @export
setMethod(
  f = "cut",
  signature = signature("IclPath"),
  definition = function(x, K) {
    if (K < x@K) {
      i <- which(sapply(x@path, function(p) {
        p$K
      }) == K)
      # Old version: x@cl = as.vector(x@path[[i]]$cl)
      # Compute cl with the history of fusions (k,l) at each stage
      for (p in 1:i) {
        # Get fusion (k,l)
        k <- x@path[[p]]$k
        l <- x@path[[p]]$l
        x@cl[x@cl == k] <- l
        # rescale @cl to be 1...K
        x@cl[x@cl > k] <- x@cl[x@cl > k] - 1
      }
      x@K <- K
      x@logalpha <- x@path[[i]]$logalpha
      x@icl <- x@path[[i]]$icl


      for (st in names(x@obs_stats)) {
        x@obs_stats[st] <- x@path[[i]]$obs_stats[st]
      }

      x@path <- x@path[(i + 1):length(x@path)]
      x
    } else {
      warning(paste0("This clustering has ", x@K, " clusters and you requested ", K, " clusters. Please provide a value for K smaller than ", x@K, "."), call. = FALSE)
    }
    x
  }
)

#' @title Plot an \code{\link{IclPath-class}} object
#'
#'
#' @param x a \code{\link{IclPath-class}}
#' @param type a string which specify plot type:
#' \itemize{
#' \item \code{'front'}: plot the extracted front ICL, log(alpha)
#' \item \code{'path'}: plot the evolution of ICL with respect to K
#' \item \code{'tree'}: plot the associated dendrogram
#' }
#' @return a ggplot graphic
#' @export
setMethod(
  f = "plot",
  signature = signature("IclPath", "missing"),
  definition = function(x, type = "tree") {
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
      plot(methods::as(x, gsub("Path", "Fit", class(x))), type = type)
    )
  }
)

#' @title Extract parameters from an \code{\link{IclFit-class}} object
#'
#' @param object a \code{\link{IclFit-class}}
#' @return a list with the model parameters estimates (MAP)
#' @details The results depends of the used model, in case the method is not yet implemented for a model, this generic method will be used. Which will return the \code{obs_stats} slot of the model.
#' @export
setMethod(
  f = "coef",
  signature = signature(object = "IclFit"),
  definition = function(object) {
    object@obs_stats
  }
)



#' @title Model based hierarchical clustering
#'
#' @description This function is the main function for fitting Dlvms with greed. 
#' In the simplest case you may only provide a dataset and greed will find a suitable one. 
#' The accepted classes for \code{X}  depends on the generative used which can be specified with the \code{model} argument. 
#' See the \code{\link{DlvmPrior-class}} and the derived classes for details.
#' 
#'
#' @param X data to cluster either a data.frame, a matrix, an array, ... depending on the used generative model
#' @param K initial number of cluster
#' @param model a generative model to fit such as \code{\link{Gmm}},\code{\link{Sbm}},.. 
#' @param alg an optimization algorithm of class \code{\link{Alg-class}} such as \code{\link{Hybrid-class}} (default), \code{\link{Multistarts-class}}, \code{\link{Seed-class}} or \code{\link{Genetic-class}}
#' @param verbose boolean value for verbose mode
#' @return an \code{\link{IclPath-class}} object
#' @examples 
#' sbm <- rsbm(50, c(0.5, 0.5), diag(2) * 0.1 + 0.01)
#' sol <- greed(sbm$x, model = Sbm())
#' table(sbm$cl,clustering(sol))
#' @export
greed <- function(X, model = find_model(X), K = 20, alg = Hybrid(), verbose = FALSE) {
  data <- preprocess(model, X)
  modelname <- toupper(class(model))
  if ("type" %in% methods::slotNames(model)) {
    modelname <- paste(model@type, modelname)
  }
  if ("sampling" %in% methods::slotNames(model)) {
    modelname <- paste(modelname, "with", model@sampling, "sampling")
  }

  cli::cli_h2("Fitting a {modelname} model")
  sol <- fit(model, alg, data, K, verbose)
  sol@obs_stats <- cleanObsStats(model, sol@obs_stats, data)

  if (length(sol@path) > 0) {
    for (p in seq_len(length(sol@path))) {
      sol@path[[p]]$obs_stats <- cleanObsStats(model, sol@path[[p]]$obs_stats, data)
    }
  }


  sol <- postprocess(sol, data)

  cli::cli_h2("Final clustering")
  cli::cli_h3("Clustering with a {toupper(class(sol@model))} model {length(sol@obs_stats$counts)} clusters and an ICL of {round(sol@icl)}")
  sol
}


find_model <- function(X) {
  if (methods::is(X, "array") && length(dim(X)) > 2) {
    dimensions <- dim(X)
    if (dimensions[1] != dimensions[2]) {
      stop(paste0("Multinomial SBM expect a cube with as many rows as columns:", dimensions, collapse = " x "))
    }
    if (length(dimensions) != 3) {
      stop(paste0("Multinomial SBM expect a cube found an array odd dimensions:", dimensions, collapse = " x "))
    }
    issym <- all(sapply(1:dim(X)[3], function(d) {
      isSymmetric(X[, , d])
    }))
    if (issym) {
      model <- MultSbm(type = "undirected")
    } else {
      model <- MultSbm()
    }
  } else {
    if (methods::is(X, "data.frame") & !all(apply(X, 2, is.numeric))) {
      if (all(sapply(X, is.factor)) | all(sapply(X, is.character))) {
        return(Lca())
      }
    }
    if (methods::is(X, "data.frame") & all(apply(X, 2, is.numeric))) {
      X <- as.matrix(X)
    }
    if (methods::is(X, "dgCMatrix") | methods::is(X, "matrix")) {
      if (nrow(X) == ncol(X)) {
        if (sum(is.na(X)) > 0) {
          stop("No missing value allowed for Sbm models. ", .call = FALSE)
        } else {
          model <- DcSbm()
        }
      } else {
        if (all(round(X) == X)) {
          model <- DcLbm()
        } else {
          model <- Gmm()
          # model = methods::new("diaggmm",mu=apply(X,2,mean),beta=0.1)
        }
      }
    } else {
      stop(paste0("Unsupported data type: ", class(X), " use a data.frame, a matrix, a sparse dgCMatrix or an array."), call. = FALSE)
    }
  }
  model
}


setGeneric("reorder", function(model, obs_stats, order) standardGeneric("reorder"))

setGeneric("fit", function(model, alg, ...) standardGeneric("fit"))


setMethod(
  f = "fit",
  signature = signature("DlvmPrior", "Hybrid"),
  definition = function(model, alg, data, K, verbose = FALSE) {
    hybrid(model, alg, data, K, verbose)
  }
)

setMethod(
  f = "fit",
  signature = signature("DlvmPrior", "Genetic"),
  definition = function(model, alg, data, k, verbose = FALSE) {
    genetic(model, alg, data, k, verbose = verbose)
  }
)


setMethod(
  f = "fit",
  signature = signature("DlvmPrior", "Multistarts"),
  definition = function(model, alg, data, K = 20, verbose = FALSE) {
    multistart(model, alg, data, K, verbose)
  }
)

setGeneric("seed", function(model, data, K) standardGeneric("seed"))

setMethod(
  f = "fit",
  signature = signature("DlvmPrior", "Seed"),
  definition = function(model, alg, data, K = 20, verbose = FALSE) {
    cl <- seed(model, data, K)
    res <- fit_greed(model, data, cl, "both", verbose = verbose)
    path <- fit_greed_path(data, res)
    p <- cleanpathopt(path)
  }
)

setGeneric("preprocess", function(model, ...) standardGeneric("preprocess"))




setGeneric("cleanObsStats", function(model, obs_stats, data) standardGeneric("cleanObsStats"))


setGeneric("postprocess", function(path, ...) standardGeneric("postprocess"))

setMethod(
  f = "postprocess",
  signature = signature("IclPath"),
  definition = function(path, data = NULL) {
    path
  }
)


setGeneric("sample_cl", function(model, data, K) standardGeneric("sample_cl"))

setMethod(
  f = "sample_cl",
  signature = signature("DlvmPrior", "list", "numeric"),
  definition = function(model, data, K) {
    sample(1:K, data$N, replace = TRUE)
  }
)

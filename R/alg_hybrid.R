#' @useDynLib greed
#' @importFrom Rcpp sourceCpp
#' @importFrom future %<-%
#' @name %<-%
NULL

#' @importFrom future %globals%
#' @name %globals%
NULL

#' @importFrom future %seed%
#' @name %seed%
NULL

#' @include models_classes.R fit_classes.R tools_cleanpath.R
#' @import Matrix
NULL

hybrid <- function(model, alg, data, K, verbose = FALSE) {
  fimerge <- function(ncl, merge_graph) {
    merge_cstr(model, data, ncl, merge_graph, verbose)
  }
  fiswap <- function(ncl, move_mat) {
    swap_cstr(model, data, ncl, move_mat, verbose = verbose)
  }

  train.hist <- data.frame(generation = c(), icl = c(), K = c())

  # multi-start in //
  # future::plan(future::multiprocess)

  solutions <- listenv::listenv()
  # first generation of solutions
  pop_size <- alg@pop_size

  sb <- cli::cli_status("{cli::symbol$info} Initializing a population of {pop_size} solutions.")
  
  for (i in 1:pop_size) {
    cli <- sample_cl(model, data, K)
    cli <- as.numeric(factor(cli))
    Ko <- max(cli)
    moves <- as.sparse(matrix(1, Ko, Ko))
    solutions[[i]] %<-% fiswap(cli, moves) %globals% c("model", "data", "cli", "verbose", "fiswap", "moves") %seed% TRUE
  }
  solutions <- as.list(solutions)
  icls <- sapply(solutions, function(s) {
    s@icl
  })
  # check for errors
  solutions <- solutions[!is.nan(icls)]
  icls <- icls[!is.nan(icls)]
  old_best <- -Inf
  best_icl <- max(icls)
  cK <- solutions[[order(icls, decreasing = TRUE)[1]]]@K
  nbgen <- 1
  # while maximum number of generation // all solutions are equals // no improvements
  pmut <- alg@prob_mutation
  Kmax <- alg@Kmax

  
  cli::cli_status_update(sb,"{cli::symbol$info} Generation {nbgen} : best solution with an ICL of {round(solutions[[which.max(icls)]]@icl)} and {solutions[[which.max(icls)]]@K} clusters.")

  while ((max(icls) - min(icls)) > 1 & nbgen < alg@nb_max_gen & best_icl > old_best & cK < Kmax) {
    train.hist <- rbind(train.hist, data.frame(generation = nbgen, icl = icls, K = sapply(solutions, function(s) {
      max(s@cl)
    })))
    # selection keep the top half solutions
    ii <- order(icls)
    Nsel <- round(alg@pop_size * 0.4)
    ii <- ii[(length(ii) - Nsel):length(ii)]
    icls <- icls[ii]
    solutions <- solutions[ii]
    bres <- solutions[[order(icls, decreasing = TRUE)[1]]]
    new_solutions <- listenv::listenv()
    for (i in 1:(alg@pop_size - 1)) {
      ip <- 1:min(c(Nsel, length(solutions)))
      i1 <- sample(ip, 1, prob = ip)
      i2 <- sample(ip[-i1], 1, prob = ip[-i1])
      s1 <- solutions[[i1]]
      s2 <- solutions[[i2]]
      new_solutions[[i]] %<-% full_cross_over(s1, s2, fimerge, fiswap, pmut, Kmax) %globals% c("s1", "s2", "fimerge", "fiswap", "pmut", "full_cross_over", "Kmax") %seed% TRUE
    }
    solutions <- c(bres, as.list(new_solutions))
    icls <- sapply(solutions, function(s) {
      s@icl
    })
    if (sum(is.na(icls)) > 0 | sum(is.nan(icls)) > 0) {
      message("NAN in objective function returning problematic solution")
      return(solutions[[which(is.nan(icls))[1]]])
    }
    solutions <- solutions[!is.nan(icls)]

    icls <- icls[!is.nan(icls)]
    old_best <- best_icl
    best_icl <- max(icls)
    cK <- solutions[[order(icls, decreasing = TRUE)[1]]]@K
    nbgen <- nbgen + 1

    cli::cli_status_update(sb, "{cli::symbol$info} Generation {nbgen} : best solution with an ICL of {round(solutions[[which.max(icls)]]@icl)} and {solutions[[which.max(icls)]]@K} clusters.")
    
  }
  if (cK > Kmax) {
    warning("The number of clusters has reached the upper limit.\n Increase Kmax (see ?Hybrid-class) if you want to explore clusterings with more clusters.")
  }

  train.hist <- rbind(train.hist, data.frame(generation = nbgen, icl = icls, K = sapply(solutions, function(s) {
    max(s@cl)
  })))

  # best solution
  res <- solutions[[order(icls, decreasing = TRUE)[1]]]

  # compute merge path
  path <- fit_greed_path(data, res)

  # clean the resuts (compute, merge tree,...)
  path <- cleanpathopt(path)
  # store train history
  path@train_hist <- train.hist
  path
}


full_cross_over <- function(sol1, sol2, fimerge, fiswap, pmutation, Kmax) {
  if (sol1@K == 1) {
    return(sol2)
  }
  if (sol2@K == 1) {
    return(sol1)
  }


  # cartesian product on the z of the two solution
  # ncl = unclass(factor(paste(sol1@cl,sol2@cl)))
  sol <- sol1
  perm <- sample(1:sol2@K, sol2@K)
  cl2 <- perm[sol2@cl]

  K2 <- sol2@K

  Ta <- table(sol@cl, cl2)
  nbclust <- sol@K - colSums(t(apply(Ta, 1, function(r) {
    cumsum(r > 0)
  })) == matrix(rowSums(Ta > 0), nrow = nrow(Ta), ncol = ncol(Ta))) + cumsum(colSums(Ta > 0))
  #  if(max(nbclust)>Kmax){
  #    K2 = which(nbclust>Kmax)[1]
  #    cl2[cl2>K2]=K2
  #  }
  ij <- which(table(sol@cl, cl2) > 0, arr.ind = TRUE)
  ncl <- as.numeric(factor(paste(sol@cl, "_", cl2, sep = ""), levels = paste(ij[, 1], "_", ij[, 2], sep = "")))


  M <- matrix(0, nrow(ij), nrow(ij))
  for (k in 1:sol1@K) {
    M[ij[, 1] == k, ij[, 1] == k] <- 1
  }
  for (k in 1:K2) {
    M[ij[, 2] == k, ij[, 2] == k] <- 1
  }
  ijAm <- which(M > 0, arr.ind = TRUE)
  ijAm <- ijAm[ijAm[, 1] != ijAm[, 2], ]
  Am <- Matrix::sparseMatrix(ijAm[, 1], ijAm[, 2], x = rep(1, nrow(ijAm)), dims = c(max(ncl), max(ncl)))
  move_mat <- Am
  if (nrow(ijAm) > 0) {
    sol <- fimerge(ncl, Matrix::tril(Am))
    move_mat <- sol@move_mat + Matrix::t(sol@move_mat)
    ncl <- sol@cl
    for (r in seq_len(nrow(move_mat))) {
      if (sum(move_mat[r, ] != 0) > 10) {
        merges <- which(move_mat[r, ] != 0)
        best_merges_row <- order(move_mat[r, merges], decreasing = TRUE)[1:10]
        move_mat[r, setdiff(merges, merges[best_merges_row])] <- 0
      }
    }
  }

  if (stats::runif(1) < pmutation) {
    sp_cl <- sample(max(ncl), 1)
    nclold <- ncl
    ncl[ncl == sp_cl] <- sample(c(sp_cl, max(ncl) + 1), sum(ncl == sp_cl), replace = TRUE)

    if (max(ncl) > nrow(move_mat) & sum(ncl == sp_cl) > 0) {
      move_mat <- cbind(move_mat, move_mat[, sp_cl])
      move_mat <- rbind(move_mat, move_mat[sp_cl, ])
      move_mat[sp_cl, max(ncl)] <- 1
      move_mat[max(ncl), sp_cl] <- 1
    } else {
      ncl <- nclold
    }
  }
  if (nrow(ijAm) > 0) {
    sol <- fiswap(ncl, move_mat)
  }
  sol
}


incremental_cross_over <- function(sol1, sol2, fimerge, fiswap, pmutation, Kmax) {
  # cartesian product on the z of the two solution
  # ncl = unclass(factor(paste(sol1@cl,sol2@cl)))
  sol <- sol1
  ncl <- sol@cl
  ncl_old <- sol@cl
  K2 <- sol2@K
  icl <- sol@icl
  for (k2 in 1:K2) {
    cl2 <- ifelse(sol2@cl == k2, 1, 2)
    ncl_old <- ncl
    ij <- which(table(ncl, cl2) > 0, arr.ind = TRUE)
    ncl <- as.numeric(factor(paste(ncl, "_", cl2, sep = ""), levels = paste(ij[, 1], "_", ij[, 2], sep = "")))
    M <- matrix(0, max(ncl), max(ncl))
    M[ij[, 2] == 1, ij[, 2] == 1] <- 1
    diag(M) <- 0
    ijm <- which(M == 1, arr.ind = TRUE)
    move_mat <- sparseMatrix(i = ijm[, 1], j = ijm[, 2], x = rep(1, nrow(ijm)), dims = c(max(ncl), max(ncl)))

    if (sum(move_mat) > 0) {
      sol <- fimerge(ncl, Matrix::tril(move_mat))
      move_mat <- sol@move_mat + Matrix::t(sol@move_mat)
      if (sol@icl > icl) {
        icl <- sol@icl
        ncl <- sol@cl
      } else {
        ncl <- ncl_old
      }
    }
  }
  M <- matrix(0, max(ncl), max(ncl))
  Tac <- table(ncl, sol1@cl)
  for (k in 1:sol1@K) {
    ibrothers <- which(Tac[, k] != 0)
    M[ibrothers, ibrothers] <- 1
  }
  Tac <- table(ncl, sol2@cl)
  for (k in 1:K2) {
    ibrothers <- which(Tac[, k] != 0)
    M[ibrothers, ibrothers] <- 1
  }
  diag(M) <- 0
  ijm <- which(M == 1, arr.ind = TRUE)
  move_mat <- sparseMatrix(i = ijm[, 1], j = ijm[, 2], x = rep(1, nrow(ijm)), dims = c(max(ncl), max(ncl)))

  sol <- fimerge(ncl, Matrix::tril(move_mat))
  move_mat <- sol@move_mat
  ncl <- sol@cl
  for (r in seq_len(nrow(move_mat))) {
    if (sum(move_mat[r, ] != 0) > 10) {
      merges <- which(move_mat[r, ] != 0)
      best_merges_row <- order(move_mat[r, merges], decreasing = TRUE)[1:10]
      move_mat[r, setdiff(merges, merges[best_merges_row])] <- 0
    }
  }

  if (stats::runif(1) < pmutation) {
    sp_cl <- sample(max(ncl), 1)
    nclold <- ncl
    ncl[ncl == sp_cl] <- sample(c(sp_cl, max(ncl) + 1), sum(ncl == sp_cl), replace = TRUE)

    if (max(ncl) > nrow(move_mat) & sum(ncl == sp_cl) > 0) {
      move_mat <- cbind(move_mat, move_mat[, sp_cl])
      move_mat <- rbind(move_mat, move_mat[sp_cl, ])
      move_mat[sp_cl, max(ncl)] <- 1
      move_mat[max(ncl), sp_cl] <- 1
    } else {
      ncl <- nclold
    }
  }

  sol <- fiswap(ncl, move_mat)

  sol
}

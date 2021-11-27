#' @useDynLib greed
#' @importFrom Rcpp sourceCpp
#' @importFrom future %<-%
#' @name %<-%
NULL

#' @include models_classes.R fit_classes.R
#' @import Matrix
NULL


# extract the pareto front
extract_front_height <- function(sol) {
  if (sol@K == 1) {
    H <- c(0)
    return(H)
  }
  # vector of icls value from root to leaves
  icl <- c(sol@icl, sapply(sol@path, function(v) {
    v$icl1
  }))
  icl <- icl[length(icl):1]
  # K
  K <- 1:length(icl)

  # init H
  H <- rep(0, length(icl))

  # current merge position
  cdi <- Inf
  # current best line
  bestline <- 1
  # vector with indexes of solutions that belong to the pareto front
  Front <- c(1)

  # from root to leaves
  for (l in 2:length(icl)) {

    # merge value with current bestline
    di <- (icl[l] - icl[bestline])
    din <- di / (l - bestline)

    # is their a potential merge ?
    if (di > 0) {

      # if this merge did not occurs after the current one update the front
      while (din > cdi & length(Front) > 1) {

        # remove the last solution from the front
        Front <- Front[-length(Front)]
        H[bestline] <- -1

        # update bestline
        bestline <- Front[length(Front)]

        # update merge position
        di <- (icl[l] - icl[bestline])
        din <- di / (l - bestline)
        # update previous merge position
        if (length(Front) > 1) {
          cdi <- (icl[bestline] - icl[Front[length(Front) - 1]]) / (bestline - Front[length(Front) - 1])
        } else {
          cdi <- H[1]
        }
      }

      # add the extracted solution to the front
      H[Front[length(Front)]] <- din
      cdi <- din
      bestline <- l
      Front <- c(Front, l)
    } else {
      # if solution not in front
      H[l] <- -1
    }
  }

  # copy from left previous value
  for (l in 2:length(icl)) {
    if (H[l] == -1) {
      H[l] <- H[l - 1]
    }
  }
  H
}

# clean the merge path
cleanpathopt <- function(pathsol) {
  pathsol@cl <- as.numeric(pathsol@cl)

  if (length(pathsol@path) > 1) {
    if (is.infinite(pathsol@path[[length(pathsol@path)]]$icl1)) {
      pathsol@path[[length(pathsol@path)]]$icl1 <- pathsol@path[[length(pathsol@path) - 1]]$icl1
      pathsol@path[[length(pathsol@path)]]$logalpha <- pathsol@path[[length(pathsol@path) - 1]]$logalpha
    }

    K <- pathsol@K
    pathsol@logalpha <- 0
    path <- pathsol@path


    # check for possible better solution than init with alpha=1 along the path
    icli <- sapply(path, function(p) {
      p$icl1
    })
    if (max(icli) > pathsol@icl) {
      im <- which.max(icli)
      K <- path[[im]]$K
      pathsol@K <- K
      pathsol@obs_stats <- path[[im]]$obs_stats
      pathsol@icl <- path[[im]]$icl1
      for (p in 1:im) {
        # Get fusion (k,l)
        k <- pathsol@path[[p]]$k
        l <- pathsol@path[[p]]$l
        pathsol@cl[pathsol@cl == k] <- l
        # rescale @cl to be 1...K
        pathsol@cl[pathsol@cl > k] <- pathsol@cl[pathsol@cl > k] - 1
      }
      if ((im + 1) <= length(path)) {
        path <- path[(im + 1):length(path)]
        pathsol@path <- path
      } else {
        pathsol@path <- list()
      }
    }


    # check for non empty path
    if (length(pathsol@path) > 1) {
      # compute the pareto front and extract the height as -log(alpha) of each merge in the front
      Hfront <- extract_front_height(pathsol)
      # initialisation
      # build the merge tree in hclust format
      merge <- c()
      cnodes <- -(1:pathsol@K)
      for (m in 1:length(path)) {
        merge <- rbind(merge, c(cnodes[path[[m]]$k], cnodes[path[[m]]$l]))
        cnodes[path[[m]]$l] <- m
        cnodes <- cnodes[-path[[m]]$k]
      }
      # find optimal leaf ordering

      dm <- -path[[1]]$merge_mat - t(path[[1]]$merge_mat)
      dm[is.infinite(dm)] <- 100 * max(dm[!is.infinite(dm)])
      if (length(path) > 1) {
        leaforder <- cba::order.optimal(stats::as.dist(dm), merge)
      } else {
        leaforder <- list(order = 1:2)
      }




      # ordering of initial solution
      pathsol@obs_stats <- reorder(pathsol@model, pathsol@obs_stats, leaforder$order)
      pathsol@cl <- order(leaforder$order)[pathsol@cl]


      # prepare the data.frame to store the tree
      ggtree <- data.frame(H = rep(0, K), tree = 0, x = seq(-1, 1, length.out = K), node = 1:K, xmin = 0, xmax = 0, K = K)
      tree <- rep(0, 2 * K - 1)
      perm <- leaforder$order
      nodes <- 1:K
      cn <- K + 1
      for (m in 1:length(path)) {
        # update the permutation
        oldperm <- perm
        perm <- perm[perm != path[[m]]$k]
        perm[perm > path[[m]]$k] <- perm[perm > path[[m]]$k] - 1
        # update the stats accordingly
        path[[m]]$obs_stats <- reorder(pathsol@model, path[[m]]$obs_stats, as.integer(perm))
        path[[m]]$merge_mat <- tril(path[[m]]$merge_mat[oldperm, oldperm] + t(path[[m]]$merge_mat[oldperm, oldperm]))
        # and the index of the merged cluster
        nkl <- sort(which(oldperm == path[[m]]$k | oldperm == path[[m]]$l))
        path[[m]]$k <- nkl[2]
        path[[m]]$l <- nkl[1]
        # build the tree
        tree[nodes[path[[m]]$k]] <- cn
        tree[nodes[path[[m]]$l]] <- cn
        xchildren <- ggtree$x[c(nodes[path[[m]]$k], nodes[path[[m]]$l])]
        ggtree <- rbind(ggtree, data.frame(H = Hfront[K - m], tree = 0, x = mean(xchildren), node = cn, xmin = min(xchildren), xmax = max(xchildren), K = K - m))
        # update nodes vector
        nodes <- nodes[-path[[m]]$k]
        nodes[path[[m]]$l] <- cn
        cn <- cn + 1
      }
      # store the tree
      ggtree$tree <- tree
      # store height and xpos of father
      ggtree$Hend <- c(ggtree$H[ggtree$tree], -1)
      ggtree$xend <- c(ggtree$x[ggtree$tree], -1)


      # store upated path and tree
      pathsol@path <- path
      pathsol@tree <- tree
      pathsol@ggtree <- ggtree[nrow(ggtree):1, ]
    } else {
      # deals with empty path
      pathsol@tree <- c(0)
      pathsol@ggtree <- data.frame(H = 0, tree = 0, x = 0, node = 1, xmin = 0, max = 0)
    }
  } else {
    # deals with empty path
    pathsol@tree <- c(0)
    pathsol@ggtree <- data.frame(H = 0, tree = 0, x = 0, node = 1, xmin = 0, max = 0)
  }
  pathsol
}


# clean the merge path
cleanpath <- function(pathsol) {
  K <- pathsol@K
  pathsol@logalpha <- 0
  path <- pathsol@path


  # check for possible better solution than init with alpha=1 along the path
  if (length(path) > 0) {
    icli <- sapply(path, function(p) {
      p$icl1
    })
    if (max(icli) > pathsol@icl) {
      im <- which.max(icli)
      K <- path[[im]]$K
      pathsol@K <- K
      pathsol@obs_stats <- path[[im]]$obs_stats
      pathsol@icl <- path[[im]]$icl1

      path <- path[(im + 1):length(path)]

      pathsol@path <- path
    }



    # check for non empty path
    if (length(path) > 0) {

      # compute the pareto front and extract the height as -log(alpha) of each merge in the front
      Hfront <- extract_front_height(pathsol)
      # initialisation
      # vector with the tree information
      tree <- c(0)
      # x position of the tree nodes the root is at 0
      xtree <- c(0)
      # current node
      cn <- 1

      # vector of length K which cluster name in tree notation
      lab <- c(1)
      # vector of length K with cluster x positions
      xpos <- c(0)
      # height of the nodes
      H <- rep(0, 2 * K - 1)
      Kc <- rep(K, 2 * K - 1)
      x <- 0
      w <- 0.5
      K <- 1
      # go back from the tree root
      best_merge <- length(path)
      # for each merge
      for (lev in seq(length(path), 1)) {

        # merge number starting from root
        pl <- length(path) - lev
        # get the order form the x position of the cluster at this level
        ord <- order(xpos)

        # reorder the cluster accordingly
        path[[lev]]$obs_stats <- reorder(pathsol@model, path[[lev]]$obs_stats, ord)


        # store in the tree the merge of k and l
        k <- path[[lev]]$k
        l <- path[[lev]]$l
        tree <- c(tree, lab[l], lab[l])


        # update height and K vectors
        Kc[lab[l]] <- K
        H[lab[l]] <- Hfront[K]

        # update lab[l] to one of the two new nodes
        lab[l] <- cn + 1

        # current father position
        fpos <- xpos[l]

        # update xtree position
        xtree <- c(xtree, fpos - w^pl, fpos + w^pl)
        # better to take size in account for left/right to avoid random effects ?

        # update xpos and lab for the two new nodes
        xpos[l] <- fpos - w^pl
        if (k > K) {
          xpos <- c(xpos, fpos + w^pl)
          lab <- c(lab, cn + 2)
        } else {
          xpos <- c(xpos[1:(k - 1)], fpos + w^pl, xpos[k:length(lab)])
          lab <- c(lab[1:(k - 1)], cn + 2, lab[k:length(lab)])
        }
        # update K
        K <- K + 1
        # update cn
        cn <- cn + 2
      }


      # prepare the data.frame to store the tree
      ggtree <- data.frame(H = H, tree = tree, x = xtree, node = 1:length(tree), xmin = 0, xmax = 0, K = Kc)
      # recompute the x bottom to top for constant spacing of leafs
      # leafs ordering from x
      leafs <- which(ggtree$H == 0)
      or <- order(ggtree[leafs, "x"])
      # equidistant spacing
      ggtree$x[leafs[or]] <- seq(-1, 1, length.out = length(leafs))
      # update x position from leafs to root is easy since leafs are now ordered
      others <- ggtree$node[ggtree$H != 0]
      for (n in others[seq(length(others), 1)]) {
        sons <- which(ggtree$tree == n)
        ggtree$x[n] <- mean(ggtree$x[sons])
        ggtree$xmin[n] <- min(ggtree$x[sons])
        ggtree$xmax[n] <- max(ggtree$x[sons])
      }

      # store height and xpos of father
      ggtree$Hend <- c(-1, ggtree$H[ggtree$tree])
      ggtree$xend <- c(-1, ggtree$x[ggtree$tree])

      # ordering of initial solution
      or <- order(xpos)
      pathsol@obs_stats <- reorder(pathsol@model, pathsol@obs_stats, or)
      pathsol@cl <- order(or)[pathsol@cl]

      # store upated path and tree
      pathsol@path <- path
      pathsol@tree <- tree
      pathsol@ggtree <- ggtree
    } else {
      # deals with empty path
      pathsol@tree <- c(0)
      pathsol@ggtree <- data.frame(H = 0, tree = 0, x = 0, node = 1, xmin = 0, max = 0)
    }
  } else {
    # deals with empty path
    pathsol@tree <- c(0)
    pathsol@ggtree <- data.frame(H = 0, tree = 0, x = 0, node = 1, xmin = 0, max = 0)
  }
  pathsol
}

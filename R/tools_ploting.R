#' @importFrom graphics plot
#' @include models_classes.R fit_classes.R





#' @title Show an IclPath object
#'
#' @description
#' Print an \code{\link{IclPath-class}} object, model type and number of found clusters are provided.
#' @param object \code{\link{IclPath-class}} object to print
#' @return None (invisible NULL). No return value, called for side effects.
#' @export
setMethod(
  f = "show",
  signature = signature("IclFit"),
  definition = function(object) {
    cli::cli_h2("Clustering with a {toupper(class(object@model))} model {length(object@obs_stats$counts)} clusters and an ICL of {round(object@icl)}")
    cli::cli_alert_info("Generic methods to explore a fit:")
    cli::cli_ul()
    cli::cli_li("?clustering, ?K, ?ICL, ?prior, ?plot, ?cut, ?coef")
    cli::cli_ul()
  }
)



pprint <- function(x, M, l) {
  K <- length(x@obs_stats$counts)
  na <- colnames(M)
  D <- Matrix::rowSums(M)
  for (k in 1:K) {
    ii <- which(x@cl == k)
    topk <- order(D[ii], decreasing = TRUE)[1:l]
    print(na[ii[topk]])
  }
}

spy <- function(x) {
  ij <- Matrix::which(x != 0, arr.ind = TRUE)
  gg <- data.frame(i = ij[, 1], j = ij[, 2], v = x[ij])
  ggplot2::ggplot(gg) +
    ggplot2::geom_point(ggplot2::aes_(y = ~ -i, x = ~j, size = ~v)) +
    ggplot2::scale_x_continuous("", c()) +
    ggplot2::scale_y_continuous("", c()) +
    ggplot2::scale_size_area(max_size = 1, guide = "none")
}

groupspy <- function(x, clust) {
  x <- x[order(clust), order(clust)]
  lims <- c(1, cumsum(table(clust)))
  ij <- Matrix::which(x != 0, arr.ind = TRUE)
  gg <- data.frame(i = ij[, 1], j = ij[, 2], v = x[ij])
  ggplot2::ggplot(gg) +
    ggplot2::geom_point(ggplot2::aes_(y = ~ -i, x = ~j, size = ~v), alpha = 0.5) +
    ggplot2::scale_x_continuous("", breaks = lims, labels = c(), minor_breaks = c(), ) +
    ggplot2::scale_y_continuous("", breaks = -lims, labels = c(), minor_breaks = c()) +
    ggplot2::scale_size_area(max_size = 1, guide = "none")
}


plot_front <- function(sol) {
  if (sol@K < 3) {
    message("The fit contains only less than 3 clusters, an empty plot was produced.")
    return(ggplot2::ggplot())
  }
  icl <- c(sol@icl, sapply(sol@path, function(v) {
    v$icl1
  }))
  ggicl <- data.frame(icl = icl[length(icl):1], K = 1:length(icl))
  # ggfront= sol@ggtree %>% mutate(x=-H) %>% select(x,K) %>% arrange(x) %>% head(sol@K) %>% left_join(ggicl) %>% mutate(xp=lag(x))
  ggfront <- merge(sol@ggtree[, c("H", "K")], ggicl)
  ggfront$x <- -ggfront$H
  ggfront <- ggfront[order(ggfront$x)[1:sol@K], ]
  ggfront$xp <- c(min(ggfront$x) - 0.05 * diff(range(ggfront$x)), ggfront$x[1:(nrow(ggfront) - 1)])
  ggfront <- ggfront[ggfront$x != ggfront$xp, ]
  ggplot2::ggplot() +
    ggplot2::geom_abline(data = ggicl, ggplot2::aes_(intercept = ~icl, slope = ~ K - 1), alpha = 0.2) +
    ggplot2::geom_point(data = ggfront, ggplot2::aes_(x = ~x, y = ~ icl + x * (K - 1))) +
    ggplot2::geom_segment(data = ggfront, ggplot2::aes_(x = ~x, y = ~ icl + x * (K - 1), xend = ~xp, yend = ~ icl + xp * (K - 1))) +
    ggplot2::scale_x_continuous(expression(paste("log(", alpha, ")")), limits = c(min(ggfront$xp), 0)) +
    ggplot2::ylab("ICL") +
    ggplot2::ggtitle(paste0(toupper(class(sol@model)), " model with : ", max(sol@cl), " clusters.")) +
    ggplot2::theme_bw()
}

lapath <- function(x) {
  if (x@K < 3) {
    message("The fit contains less than 3 clusters, an empty plot was produced.")
    return(ggplot2::ggplot())
  }
  gg <- data.frame(k = sapply(x@path, function(p) {
    p$K
  }), logalpha = sapply(x@path, function(p) {
    p$logalpha
  }))
  gg <- rbind(gg, data.frame(k = length(x@obs_stats$counts), logalpha = x@logalpha))
  ggplot2::ggplot(data = gg) +
    ggplot2::geom_line(ggplot2::aes_(x = ~k, y = ~ -logalpha)) +
    ggplot2::geom_point(ggplot2::aes_(x = ~k, y = ~ -logalpha)) +
    ggplot2::ylab(expression(paste("-log(", alpha, ")"))) +
    ggplot2::ggtitle(paste0(toupper(class(x@model)), " ", length(x@obs_stats$counts), " clusters")) +
    ggplot2::theme_bw()
}


iclpath <- function(x) {
  if (x@K == 1) {
    message("The fit contains only one cluster, an empty plot was produced.")
    return(ggplot2::ggplot())
  }
  gg <- data.frame(k = sapply(x@path, function(p) {
    length(p$counts)
  }), icl = sapply(x@path, function(p) {
    p$icl
  }))
  gg <- rbind(gg, data.frame(k = length(x@obs_stats$counts), icl = x@icl))
  ggplot2::ggplot(data = gg) +
    ggplot2::geom_line(ggplot2::aes_(x = ~k, y = ~icl)) +
    ggplot2::geom_point(ggplot2::aes_(x = ~k, y = ~icl)) +
    ggplot2::ylab(expression(paste("ICL"))) +
    ggplot2::ggtitle(paste0(toupper(class(x@model)), " ", length(x@obs_stats$counts), " clusters")) +
    ggplot2::theme_bw()
}


# Node link visualisations

nodelink <- function(sol) {
  if (!(methods::is(sol, "SbmFit") | methods::is(sol, "DcSbmFit"))) {
    stop("Nodes and Links diagrams only available for Sbm and DcSbm models", .call = FALSE)
  }
  if (methods::is(sol, "SbmFit")) {
    x_counts <- sol@obs_stats$Sbm$x_counts
  }

  if (methods::is(sol, "DcSbmFit")) {
    x_counts <- sol@obs_stats$DcSbm$x_counts
  }
  ij <- Matrix::which(x_counts > 0, arr.ind = TRUE)
  ld <- x_counts
  # /(sol@obs_stats$counts%*%t(sol@obs_stats$counts))
  ij <- ij[ij[, 1] != ij[, 2], ]
  gglink <- data.frame(from = ij[, 1], to = ij[, 2], p = ld[ij])
  gglink$y <- ifelse(gglink$from < gglink$to, -0.3, 0.3)
  ggnode <- data.frame(i = 1:length(sol@obs_stats$counts), pi = diag(x_counts))
  gl <- ggplot2::guide_legend()
  ggplot2::ggplot() +
    ggplot2::geom_curve(data = gglink, ggplot2::aes_(x = ~from, xend = ~to, y = ~y, yend = ~y, size = ~p, alpha = ~p), arrow = grid::arrow(length = grid::unit(2, "mm")), curvature = 0.7) +
    ggplot2::scale_x_continuous("", c()) +
    ggplot2::scale_y_continuous("", c(), limits = c(-5, 5)) +
    ggplot2::scale_alpha("Link density:", limits = c(0, max(gglink$p)), guide = "none") +
    ggplot2::scale_size_area("Clusters size:", limits = c(0, max(ggnode$pi)), max_size = 4, guide = "none") +
    ggplot2::geom_point(data = ggnode, ggplot2::aes_(x = ~i, y = ~0, size = ~pi)) +
    ggplot2::ggtitle(paste0(toupper(class(sol@model)), " model with : ", max(sol@cl), " clusters.")) +
    ggplot2::theme_minimal()
}



# nodelinklab
nodelinklab <- function(sol, labels, s = 0) {
  ij <- Matrix::which(sol@obs_stats$x_counts > 0, arr.ind = TRUE)
  ld <- sol@obs_stats$x_counts
  # /(sol@obs_stats$counts%*%t(sol@obs_stats$counts))
  ij <- ij[ij[, 1] != ij[, 2], ]
  gglink <- data.frame(from = ij[, 1], to = ij[, 2], p = ld[ij])
  gglink$y <- ifelse(gglink$from < gglink$to, -0.3, 0.3)
  ggnode <- data.frame(i = 1:length(sol@obs_stats$counts), pi = diag(sol@obs_stats$x_counts))
  gl <- ggplot2::guide_legend()
  ggplot2::ggplot() +
    ggplot2::geom_curve(data = gglink[gglink$p > s, ], ggplot2::aes_(x = ~from, xend = ~to, y = ~y, yend = ~y, size = ~p, alpha = ~p), arrow = grid::arrow(length = grid::unit(2, "mm")), curvature = 0.7) +
    ggplot2::scale_x_continuous("", c()) +
    ggplot2::scale_y_continuous("", c(), limits = c(-6, 6)) +
    ggplot2::scale_alpha("Link density:", limits = c(0, max(gglink$p)), guide = "none") +
    ggplot2::scale_size_area("Clusters size:", limits = c(0, max(c(ggnode$pi, gglink$p))), max_size = 4, guide = "none") +
    ggplot2::geom_point(data = ggnode, ggplot2::aes_(x = ~i, y = ~ -0.1, size = ~pi)) +
    ggplot2::ggtitle(paste0(toupper(class(sol@model)), " model with : ", max(sol@cl), " clusters.")) +
    ggplot2::geom_text(data = data.frame(x = 1:length(labels), label = labels), ggplot2::aes_(x = ~x, label = ~label, y = ~0.05)) +
    ggplot2::theme_minimal()
}

co_nodelink <- function(sol) {
  ij <- Matrix::which(sol@obs_stats$DcLbm$co_x_counts > 0, arr.ind = TRUE)
  ld <- sol@obs_stats$DcLbm$co_x_counts
  cccol <- as.numeric(table(sol@clcol))
  ccrow <- as.numeric(table(sol@clrow))

  gglink <- data.frame(from = ij[, 1], to = ij[, 2], p = ld[ij])
  gglink$nrow <- ccrow[gglink$from]
  gglink$ncol <- cccol[gglink$to]
  gglink$pn <- gglink$p / (gglink$nrow * gglink$ncol)
  gglink <- gglink[order(gglink$pn), ]
  ggnode <- rbind(data.frame(type = "col", n = cccol, i = 1:sol@Kcol), data.frame(type = "row", n = ccrow, i = 1:sol@Krow))


  ggplot2::ggplot() +
    ggplot2::geom_curve(data = gglink, ggplot2::aes_(y = ~from, yend = ~ max(from) + 1, x = -5, xend = ~to, size = ~p, alpha = ~p), curvature = 0.35) +
    ggplot2::geom_point(data = ggnode, ggplot2::aes_(x = ~ ifelse(type == "row", -5, i), y = ~ ifelse(type == "row", i, max(gglink$from) + 1), size = ~ n^2)) +
    ggplot2::theme_minimal() +
    ggplot2::scale_size_area("Clusters size:", limits = c(0, max(ggnode$n)^2), max_size = 7, guide = "none") +
    ggplot2::scale_y_continuous("", c()) +
    ggplot2::scale_x_continuous("", c()) +
    ggplot2::scale_alpha("Link density:", limits = c(0, max(gglink$p)), guide = "none") +
    ggplot2::ggtitle("Poisson, co-clustering")
}

# dendogram visualisation
dendo <- function(x) {
  if (x@K < 3) {
    message("The fit contains less than 3 clusters, an empty plot was produced.")
    return(ggplot2::ggplot())
  }
  ggtree <- x@ggtree
  tree <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = ggtree[ggtree$node %in% ggtree$tree, ], ggplot2::aes_(x = ~xmin, y = ~H, xend = ~xmax, yend = ~H)) +
    ggplot2::geom_segment(data = ggtree[-1, ], ggplot2::aes_(x = ~x, y = ~H, xend = ~x, yend = ~Hend)) +
    ggplot2::scale_x_continuous("", breaks = c()) +
    ggplot2::ylab(expression(paste("-log(", alpha, ")"))) +
    ggplot2::ggtitle(paste0(toupper(class(x@model)), " ", length(x@obs_stats$counts), " clusters, dendogram")) +
    ggplot2::theme_bw()
  if (x@K < max(x@ggtree$K)) {
    hc <- x@ggtree$H[x@ggtree$K == x@K]
    tree <- tree + ggplot2::geom_hline(ggplot2::aes(yintercept = hc), color = "red", size = 0.8, linetype = "dashed", alpha = 0.8)
  }
  tree
}

co_dendo <- function(x) {
  ggtree <- x@ggtreerow
  rowt <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = ggtree[ggtree$node %in% ggtree$tree, ], ggplot2::aes_(x = ~xmin, y = ~H, xend = ~xmax, yend = ~H)) +
    ggplot2::geom_segment(data = ggtree, ggplot2::aes_(x = ~x, y = ~H, xend = ~x, yend = ~Hend)) +
    ggplot2::scale_x_continuous("", breaks = c()) +
    ggplot2::ylab(expression(paste("-log(", alpha, ")"))) +
    ggplot2::ggtitle(paste0(x@Krow, " row clusters")) +
    ggplot2::theme_bw()
  if (x@K < max(x@ggtree$K)) {
    hc <- x@ggtree$H[x@ggtree$K == x@K]
    rowt <- rowt + ggplot2::geom_hline(ggplot2::aes(yintercept = hc), color = "red", size = 0.8, linetype = "dashed", alpha = 0.8)
  }

  ggtree <- x@ggtreecol
  colt <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = ggtree[ggtree$node %in% ggtree$tree, ], ggplot2::aes_(x = ~xmin, y = ~H, xend = ~xmax, yend = ~H)) +
    ggplot2::geom_segment(data = ggtree, ggplot2::aes_(x = ~x, y = ~H, xend = ~x, yend = ~Hend)) +
    ggplot2::scale_x_continuous("", breaks = c()) +
    ggplot2::ylab(expression(paste("-log(", alpha, ")"))) +
    ggplot2::ggtitle(paste0(x@Kcol, " column clusters")) +
    ggplot2::theme_bw()
  if (x@K < max(x@ggtree$K)) {
    hc <- x@ggtree$H[x@ggtree$K == x@K]
    colt <- colt + ggplot2::geom_hline(ggplot2::aes(yintercept = hc), color = "red", size = 0.8, linetype = "dashed", alpha = 0.8)
  }
  gridExtra::grid.arrange(rowt, colt, nrow = 1)
}


# matrice blocks visualisation

graph_blocks <- function(x) {
  if (!(methods::is(x, "SbmFit") | methods::is(x, "DcSbmFit"))) {
    stop("Nodes and Links diagrams only available for Sbm and DcSbm models", .call = FALSE)
  }
  if (methods::is(x, "SbmFit")) {
    x_counts <- x@obs_stats$Sbm$x_counts
  }

  if (methods::is(x, "DcSbmFit")) {
    x_counts <- x@obs_stats$DcSbm$x_counts
  }

  K <- length(x@obs_stats$counts)
  gg <- data.frame(
    kc = rep(cumsum(x@obs_stats$counts), each = K),
    lc = rep(cumsum(x@obs_stats$counts), K),
    sizek = rep(x@obs_stats$counts, each = K),
    sizel = rep(x@obs_stats$counts, K),
    count = as.vector(x_counts)
  )
  ggplot2::ggplot(gg[gg$count > 0, ]) +
    ggplot2::geom_tile(ggplot2::aes_(x = ~ kc - sizek / 2, y = ~ lc - sizel / 2, width = ~sizek, height = ~sizel, fill = ~ count / (sizek * sizel), alpha = ~ count / (sizek * sizel))) +
    ggplot2::scale_fill_distiller("Link density", palette = "YlOrRd", direction = 1, guide = ggplot2::guide_legend(), limits = c(0, max(gg$count / (gg$sizek * gg$sizel)))) +
    ggplot2::scale_alpha("Link density", range = c(0, 1), limits = c(0, max(gg$count / (gg$sizek * gg$sizel)))) +
    ggplot2::ggtitle(paste0(toupper(class(x@model)), " model with : ", max(x@cl), " clusters.")) +
    ggplot2::scale_x_continuous("",
      breaks = cumsum(x@obs_stats$counts),
      labels = ifelse(x@obs_stats$counts / sum(x@obs_stats$counts) > 0.05, paste0(round(100 * x@obs_stats$counts / sum(x@obs_stats$counts)), "%"), ""),
      minor_breaks = NULL, expand = ggplot2::expansion(mult = 0, add = 0)
    ) +
    ggplot2::scale_y_continuous("",
      breaks = cumsum(x@obs_stats$counts),
      labels = ifelse(x@obs_stats$counts / sum(x@obs_stats$counts) > 0.05, paste0(round(100 * x@obs_stats$counts / sum(x@obs_stats$counts)), "%"), ""),
      minor_breaks = NULL, expand = ggplot2::expansion(mult = 0, add = 0)
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw()
}


co_blocks <- function(x) {
  K <- x@Krow
  D <- x@Kcol
  ccrow <- cumsum(table(x@clrow))
  cccol <- cumsum(table(x@clcol))
  gg <- data.frame(
    kc = rep(ccrow, D),
    lc = rep(cccol, each = K),
    sizek = rep(table(x@clrow), D),
    sizel = rep(table(x@clcol), each = K),
    count = as.vector(x@obs_stats$DcLbm$co_x_counts)
  )
  ylab <- round(100 * table(x@clcol) / sum(table(x@clcol)))
  xlab <- round(100 * table(x@clrow) / sum(table(x@clrow)))

  ggplot2::ggplot(gg) +
    ggplot2::geom_tile(ggplot2::aes_(y = ~ kc - sizek / 2, x = ~ lc - sizel / 2, height = ~sizek, width = ~sizel, fill = ~ count / (sizek * sizel), alpha = ~ count / (sizek * sizel))) +
    ggplot2::scale_fill_distiller("E[X]", palette = "YlOrRd", direction = 1, guide = ggplot2::guide_legend(), limits = c(0, max(gg$count / (gg$sizek * gg$sizel)))) +
    ggplot2::scale_alpha("E[X]", range = c(0, 1), limits = c(0, max(gg$count / (gg$sizek * gg$sizel)))) +
    ggplot2::ggtitle(paste0("Co-clustering with : ", max(x@cl), " clusters.")) +
    ggplot2::scale_x_continuous("Col clusters", breaks = cccol, labels = ifelse(ylab > 5, paste0(ylab, "%"), ""), minor_breaks = NULL, expand = ggplot2::expansion(mult = 0, add = 0)) +
    ggplot2::scale_y_continuous("Row clusters", breaks = ccrow, labels = ifelse(xlab > 5, paste0(xlab, "%"), ""), minor_breaks = NULL, expand = ggplot2::expansion(mult = 0, add = 0)) +
    ggplot2::theme_bw()
}

mat_blocks <- function(x) {
  K <- length(x@obs_stats$counts)
  D <- dim(x@obs_stats$MoM$x_counts)[1]
  gg <- data.frame(
    kc = rep(cumsum(x@obs_stats$counts), D),
    lc = rep(1:D, each = K),
    sizek = rep(x@obs_stats$counts, D),
    sizel = rep(1, K * D),
    count = as.vector(Matrix::t(x@obs_stats$MoM$x_counts) / Matrix::rowSums(Matrix::t(x@obs_stats$MoM$x_counts)))
  )

  ggplot2::ggplot(gg) +
    ggplot2::geom_tile(ggplot2::aes_(y = ~ kc - sizek / 2, x = ~ lc - sizel / 2, height = ~sizek, width = ~sizel, fill = ~ log(count), alpha = ~count)) +
    ggplot2::scale_fill_distiller("E[X]", palette = "YlOrRd", direction = 1, guide = ggplot2::guide_legend(), limits = c(1, log(max(gg$count)))) +
    ggplot2::scale_alpha("E[X]", range = c(0, 1), limits = c(0, max(gg$count))) +
    ggplot2::ggtitle(paste0("MM Model with : ", max(x@cl), " clusters.")) +
    ggplot2::scale_x_continuous("Features", breaks = 1:D, labels = rep("", D), minor_breaks = NULL, expand = ggplot2::expansion(mult = 0, add = 0)) +
    ggplot2::scale_y_continuous("Clusters", breaks = cumsum(x@obs_stats$counts), labels = paste0(round(100 * x@obs_stats$counts / sum(x@obs_stats$counts)), "%"), minor_breaks = NULL, expand = ggplot2::expansion(mult = 0, add = 0)) +
    ggplot2::theme_bw()
}



# graph_balance
graph_balance <- function(x) {
  K <- length(x@obs_stats$counts)
  B <- x@obs_stats$x_counts - t(x@obs_stats$x_counts)
  gg <- data.frame(
    kc = rep(cumsum(x@obs_stats$counts), each = K),
    lc = rep(cumsum(x@obs_stats$counts), K),
    sizek = rep(x@obs_stats$counts, each = K),
    sizel = rep(x@obs_stats$counts, K),
    dk = as.vector(x@obs_stats$dout %*% t(x@obs_stats$din)),
    count = as.vector(B)
  )
  vm <- max(abs(gg$count))
  ggplot2::ggplot(gg) +
    ggplot2::geom_tile(ggplot2::aes_(x = ~ kc - sizek / 2, y = ~ lc - sizel / 2, width = ~sizek, height = ~sizel, fill = ~count), alpha = 0.7) +
    ggplot2::scale_fill_distiller("Balance :", direction = 1, palette = "RdBu", limits = c(-vm, vm)) +
    ggplot2::ggtitle(paste0(toupper(class(x@model)), " model with : ", max(x@cl), " clusters.")) +
    ggplot2::scale_x_continuous("", breaks = cumsum(x@obs_stats$counts), labels = ifelse(x@obs_stats$counts / sum(x@obs_stats$counts) > 0.05, paste0(round(100 * x@obs_stats$counts / sum(x@obs_stats$counts)), "%"), ""), minor_breaks = NULL) +
    ggplot2::scale_y_continuous("", breaks = cumsum(x@obs_stats$counts), labels = ifelse(x@obs_stats$counts / sum(x@obs_stats$counts) > 0.05, paste0(round(100 * x@obs_stats$counts / sum(x@obs_stats$counts)), "%"), ""), minor_breaks = NULL) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw()
}




mat_reg <- function(x) {
  K <- length(x@obs_stats$counts)
  D <- length(x@obs_stats$regs[[1]]$mu)
  gg <- data.frame(
    kc = rep(cumsum(x@obs_stats$counts), D),
    lc = rep(1:D, each = K),
    sizek = rep(x@obs_stats$counts, D),
    sizel = rep(1, K * D),
    count = as.vector(sapply(x@obs_stats$regs, function(reg) {
      reg$mu
    }))
  )
  ggplot2::ggplot(gg) +
    ggplot2::geom_tile(ggplot2::aes_(y = ~ kc - sizek / 2, x = ~ lc - sizel / 2, height = ~sizek, width = ~sizel, fill = ~count, alpha = ~count)) +
    ggplot2::scale_fill_distiller(expression(paste(" ", beta, " ")), type = "seq", direction = 1, palette = 2) +
    ggplot2::scale_alpha(expression(paste(" ", beta, " "))) +
    ggplot2::ggtitle(paste0("Mixture of Regression Model with : ", max(x@cl), " clusters.")) +
    ggplot2::scale_x_continuous("Features", breaks = 1:D, labels = rep("", D), minor_breaks = NULL) +
    ggplot2::scale_y_continuous("Clusters", breaks = cumsum(x@obs_stats$counts), labels = paste0(round(100 * x@obs_stats$counts / sum(x@obs_stats$counts)), "%"), minor_breaks = NULL) +
    ggplot2::theme_bw()
}




mat_reg_line <- function(x, Xd, yd) {
  K <- length(x@obs_stats$counts)
  D <- length(x@obs_stats$regs[[1]]$mu)
  ggp <- data.frame(x = Xd[, 1], y = yd, K = factor(x@cl, levels = 1:K))
  gg <- data.frame(
    y = as.vector(cbind(seq(min(ggp$x), max(ggp$x), length.out = 20), rep(1, 20)) %*% sapply(x@obs_stats$regs, function(reg) {
      reg$mu
    })),
    x = rep(seq(min(ggp$x), max(ggp$x), length.out = 20), K), K = rep(1:K, each = 20)
  )
  ggplot2::ggplot() +
    ggplot2::geom_point(data = ggp[, 1:2], ggplot2::aes_(x = ~x, y = ~y), alpha = 0.05) +
    ggplot2::geom_path(data = gg, ggplot2::aes_(x = ~x, y = ~y, group = ~K)) +
    ggplot2::geom_point(data = ggp, ggplot2::aes_(x = ~x, y = ~y, col = ~K)) +
    ggplot2::ggtitle(paste0("Mixture of Regression Model with : ", max(x@cl), " clusters.")) +
    ggplot2::facet_wrap(~K) +
    ggplot2::theme_bw()
}







bi_plot <- function(x) {
  K <- x@Krow
  D <- x@Kcol
  ccrow <- cumsum(table(x@clrow))
  cccol <- cumsum(table(x@clcol))
  gg <- data.frame(
    kc = rep(ccrow, D),
    lc = rep(cccol, each = K),
    sizek = rep(table(x@clrow), D),
    sizel = rep(table(x@clcol), each = K),
    count = as.vector(x@obs_stats$DcLbm$co_x_counts)
  )
  ylab <- round(100 * table(x@clcol) / sum(table(x@clcol)))
  xlab <- round(100 * table(x@clrow) / sum(table(x@clrow)))


  theme_mm <- ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_line(colour = "grey50"),
    panel.grid.major = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_line(size = 0),
    panel.ontop = TRUE,
    panel.background = ggplot2::element_rect(fill = NA)
  )
  theme_nolabY <- ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )
  theme_nolabX <- ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )

  mat <- ggplot2::ggplot(gg) +
    ggplot2::geom_tile(ggplot2::aes_(y = ~ kc - sizek / 2, x = ~ lc - sizel / 2, height = ~sizek, width = ~sizel, fill = ~ count / (sizek * sizel), alpha = ~ count / (sizek * sizel))) +
    ggplot2::scale_fill_distiller("E[X]", palette = "YlOrRd", direction = 1, guide = ggplot2::guide_legend(direction = "horizontal", title.position = "top", label.position = "bottom"), limits = c(0, max(gg$count / (gg$sizek * gg$sizel)))) +
    ggplot2::scale_alpha("E[X]", range = c(0, 1), limits = c(0, max(gg$count / (gg$sizek * gg$sizel))), guide = ggplot2::guide_legend(direction = "horizontal", title.position = "top", label.position = "bottom")) +
    ggplot2::scale_x_continuous("Col clusters", breaks = cccol - x@obs_stats$cols_counts / 2, minor_breaks = c(0, cccol), labels = ifelse(ylab > 5, paste0(ylab, "%"), ""), limits = c(-1, max(cccol) + 1), position = "top", expand = c(0, 0)) +
    ggplot2::scale_y_continuous("Row clusters", breaks = ccrow - x@obs_stats$rows_counts / 2, minor_breaks = c(0, ccrow), labels = ifelse(xlab > 5, paste0(xlab, "%"), ""), limits = c(-1, max(ccrow) + 1), expand = c(0, 0)) +
    ggplot2::theme_bw() +
    theme_mm
  legend <- gtable::gtable_filter(ggplot2::ggplotGrob(mat), "guide-box")




  treerow <- x@ggtreerow
  if (x@K < max(x@ggtree$K)) {
    treerow <- treerow[treerow$K <= x@K, ]
  }
  ggtree <- update_tree_prop(treerow, x@obs_stats$rows_counts)
  rowt <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = ggtree[ggtree$node %in% ggtree$tree, ], ggplot2::aes_(y = ~xmin, x = ~H, yend = ~xmax, xend = ~H)) +
    ggplot2::geom_segment(data = ggtree[ggtree$Hend != -1, ], ggplot2::aes_(x = ~H, y = ~x, yend = ~x, xend = ~Hend)) +
    ggplot2::scale_x_reverse("", position = "top") +
    ggplot2::scale_y_continuous("Row clusters", breaks = ccrow - x@obs_stats$rows_counts / 2, minor_breaks = ccrow, labels = ifelse(xlab > 5, paste0(xlab, "%"), ""), limits = c(0, max(ccrow) + 1), expand = c(0, 0)) +
    # ggplot2::scale_y_continuous("",breaks=c(),limits = c(0,max(ccrow)+1))+
    ggplot2::theme_bw() +
    theme_nolabY

  treecol <- x@ggtreecol
  if (x@K < max(x@ggtree$K)) {
    treecol <- treecol[treecol$K <= x@K, ]
  }
  ggtree <- update_tree_prop(treecol, x@obs_stats$cols_counts)
  colt <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = ggtree[ggtree$node %in% ggtree$tree, ], ggplot2::aes_(x = ~xmin, y = ~H, xend = ~xmax, yend = ~H)) +
    ggplot2::geom_segment(data = ggtree[ggtree$Hend != -1, ], ggplot2::aes_(x = ~x, y = ~H, xend = ~x, yend = ~Hend)) +
    # ggplot2::scale_x_continuous("",breaks=c(),limits = c(0,max(cccol)+1))+
    ggplot2::scale_x_continuous("Col clusters", breaks = cccol - x@obs_stats$cols_counts / 2, minor_breaks = cccol, labels = ifelse(ylab > 5, paste0(ylab, "%"), ""), limits = c(0, max(cccol) + 1), position = "top", expand = c(0, 0)) +
    ggplot2::scale_y_continuous("") +
    ggplot2::theme_bw() +
    theme_nolabX

  growt <- ggplot2::ggplotGrob(rowt)
  gcolt <- ggplot2::ggplotGrob(colt)
  glegend <- ggplot2::ggplotGrob(ggplot2::ggplot() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white")) +
    ggplot2::annotation_custom(legend))
  gmat <- ggplot2::ggplotGrob(mat + ggplot2::theme(legend.position = "hidden"))
  grob.mat <- matrix(list(glegend, growt, gcolt, gmat), nrow = 2)
  gt <- gtable::gtable_matrix("demo", grob.mat, ggplot2::unit(c(0.45, 1), "null"), ggplot2::unit(c(0.45, 1), "null"))
  gt <- gtable::gtable_add_grob(gt, legend, t = 1, l = 2)


  gridExtra::grid.arrange(top = paste0("Co-clustering with : ", max(x@cl), " clusters."), gt)
}

update_tree_prop <- function(tree, counts) {
  cumcounts <- cumsum(counts)
  tree$x[tree$H == 0] <- cumcounts[length(counts):1] - counts[length(counts):1] / 2
  tree$size[tree$H == 0] <- counts[length(counts):1]
  fathers <- unique(tree$tree[tree$H == 0])
  while (any(fathers != 0)) {
    newf <- c()
    for (f in fathers) {
      tree$x[tree$node == f] <- mean(tree$x[tree$tree == f])
      tree$xmin[tree$node == f] <- min(tree$x[tree$tree == f])
      tree$xmax[tree$node == f] <- max(tree$x[tree$tree == f])
      tree$size[tree$node == f] <- sum(tree$size[tree$tree == f])
      newf <- unique(c(newf, tree$tree[tree$node == f]))
    }
    fathers <- newf
  }
  tree
}

graph_blocks_cube <- function(x) {
  K <- length(x@obs_stats$counts)
  M <- dim(x@obs_stats$MultSbm$x_counts)[3]
  ggl <- lapply(1:M, function(m) {
    data.frame(
      kc = rep(cumsum(x@obs_stats$counts), each = K),
      lc = rep(cumsum(x@obs_stats$counts), K),
      sizek = rep(x@obs_stats$counts, each = K),
      sizel = rep(x@obs_stats$counts, K),
      count = as.vector(x@obs_stats$MultSbm$x_counts[, , m]), m = paste("Slice ", m), stringsAsFactors = FALSE
    )
  })

  gg <- do.call(rbind, ggl)

  ggplot2::ggplot(gg[gg$count > 0, ]) +
    ggplot2::geom_tile(ggplot2::aes_(x = ~ kc - sizek / 2, y = ~ lc - sizel / 2, width = ~sizek, height = ~sizel, fill = ~ count / (sizek * sizel), alpha = ~ count / (sizek * sizel))) +
    ggplot2::facet_wrap(~m) +
    ggplot2::scale_fill_distiller("Link density", palette = "YlOrRd", direction = 1, guide = ggplot2::guide_legend(), limits = c(0, max(gg$count / (gg$sizek * gg$sizel)))) +
    ggplot2::scale_alpha("Link density", range = c(0, 1), limits = c(0, max(gg$count / (gg$sizek * gg$sizel)))) +
    ggplot2::ggtitle(paste0(toupper(class(x@model)), " model with : ", max(x@cl), " clusters.")) +
    ggplot2::scale_x_continuous("",
      breaks = cumsum(x@obs_stats$counts),
      labels = ifelse(x@obs_stats$counts / sum(x@obs_stats$counts) > 0.05, paste0(round(100 * x@obs_stats$counts / sum(x@obs_stats$counts)), "%"), ""),
      minor_breaks = NULL, expand = ggplot2::expansion(mult = 0, add = 0)
    ) +
    ggplot2::scale_y_continuous("",
      breaks = cumsum(x@obs_stats$counts),
      labels = ifelse(x@obs_stats$counts / sum(x@obs_stats$counts) > 0.05, paste0(round(100 * x@obs_stats$counts / sum(x@obs_stats$counts)), "%"), ""),
      minor_breaks = NULL, expand = ggplot2::expansion(mult = 0, add = 0)
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw()
}


nodelink_cube <- function(sol) {
  M <- dim(sol@obs_stats$MultSbm$x_counts)[3]
  gglink_list <- lapply(1:M, function(m) {
    ij <- which(sol@obs_stats$MultSbm$x_counts[, , m] > 0, arr.ind = TRUE)
    ld <- sol@obs_stats$MultSbm$x_counts[, , m]
    # /(sol@obs_stats$counts%*%t(sol@obs_stats$counts))
    ij <- ij[ij[, 1] != ij[, 2], ]
    gglink <- data.frame(from = ij[, 1], to = ij[, 2], p = ld[ij])
    gglink$y <- ifelse(gglink$from < gglink$to, -0.3, 0.3)
    gglink$m <- paste("Slice", m)
    gglink
  })
  gglink <- do.call(rbind, gglink_list)
  ggnode_list <- lapply(1:M, function(m) {
    data.frame(i = 1:length(sol@obs_stats$counts), pi = diag(sol@obs_stats$MultSbm$x_counts[, , m]), m = paste("Slice", m), stringsAsFactors = FALSE)
  })
  ggnode <- do.call(rbind, ggnode_list)
  gl <- ggplot2::guide_legend()
  ggplot2::ggplot() +
    ggplot2::geom_curve(data = gglink, ggplot2::aes_(x = ~from, xend = ~to, y = ~y, yend = ~y, size = ~p, alpha = ~p), arrow = grid::arrow(length = grid::unit(2, "mm")), curvature = 0.7) +
    ggplot2::facet_wrap(~m) +
    ggplot2::scale_x_continuous("", c()) +
    ggplot2::scale_y_continuous("", c(), limits = c(-5, 5)) +
    ggplot2::scale_alpha("Link density:", limits = c(0, max(gglink$p)), guide = "none") +
    ggplot2::scale_size_area("Clusters size:", limits = c(0, max(ggnode$pi)), max_size = 4, guide = "none") +
    ggplot2::geom_point(data = ggnode, ggplot2::aes_(x = ~i, y = ~0, size = ~pi)) +
    ggplot2::ggtitle(paste0(toupper(class(sol@model)), " model with : ", max(sol@cl), " clusters.")) +
    ggplot2::theme_minimal()
}

#' @title Make a matrix of plots with a given data and gmm fitted parameters
#'
#' @description
#' Make a matrix of plots with a given data and gmm fitted parameters with ellipses.
#' @param sol a \code{\link{GmmFit-class}} or \code{\link{DiagGmmFit-class}}
#' @param X the data used for the fit a data.frame or matrix.
#' @return a \code{\link{ggplot2}} graphic
#' @export
gmmpairs <- function(sol, X) {
  if (!(methods::is(sol, "GmmFit") | methods::is(sol, "DiagGmmFit"))) {
    stop("Input sol must be a gmm_fit or diaggmm_fit object.", call. = FALSE)
  }
  if (!(methods::is(X, "matrix") | methods::is(X, "data.frame"))) {
    stop("Input X must be a matrix or data.frame.", call. = FALSE)
  }
  if (nrow(X) != length(sol@cl) | ncol(X) != length(coef(sol)$muk[[1]])) {
    stop("Dimension mismatch between the fitted model and the data.", call. = FALSE)
  }
  
  # get current ggplot's theme options (size, etc.)
  curr_t = ggplot2::theme_get()
  
  vnames <- names(X)
  plts.df <- list()
  ii <- 1
  nb_ech <- 500
  X$greedcl <- factor(sol@cl)
  if (nrow(X) > 1000) {
    X <- X[sample(1:nrow(X), 1000), ]
  }
  params <- coef(sol)
  for (i in 1:(ncol(X) - 1)) {
    for (j in 1:(ncol(X) - 1)) {
      if (i == j) {
        # plt = ggally_text(paste("Plot #", i, sep = ""))
        limsi <- c(range(X[, vnames[i]])[1] - diff(range(X[, vnames[i]])) * 0.1, range(X[, i])[2] + diff(range(X[, vnames[i]])) * 0.1)
        x <- seq(limsi[1], limsi[2], length.out = nb_ech)
        pdfs <- lapply(1:sol@K, function(k) {
          params$pi[k] * stats::dnorm(x, mean = params$muk[[k]][j], sd = sqrt(params$Sigmak[[k]][j, j]))
        })
        pdfs.df <- data.frame(pdf = do.call(c, pdfs), cl = rep(1:sol@K, each = nb_ech), x = rep(x, sol@K))
        pdfmix <- data.frame(pdf = colSums(do.call(rbind, pdfs)), x = x)

        plt <- ggplot2::ggplot() +
          ggplot2::geom_line(data = pdfmix, ggplot2::aes_(x = ~x, y = ~pdf)) +
          ggplot2::geom_line(data = pdfs.df, ggplot2::aes_(x = ~x, y = ~pdf, color = ~ factor(cl), group = ~cl)) +
          ggplot2::geom_segment(data = X, ggplot2::aes_(x = as.name(vnames[i]), xend = as.name(vnames[i]), y = 0, yend = max(pdfs.df$pdf) * 0.05, col = ~greedcl)) +
          ggplot2::scale_x_continuous(limits = limsi) +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
      } else {
        elipses.df <- do.call(rbind, lapply(1:sol@K, function(k) {
          el.df <- ellips(params$muk[[k]][c(i, j)], params$Sigmak[[k]][c(i, j), c(i, j)], nb_ech = nb_ech)
          el.df$cl <- k
          el.df
        }))

        limsi <- c(range(X[, vnames[i]])[1] - diff(range(X[, vnames[i]])) * 0.1, range(X[, i])[2] + diff(range(X[, vnames[i]])) * 0.1)
        limsj <- c(range(X[, vnames[j]])[1] - diff(range(X[, vnames[j]])) * 0.1, range(X[, j])[2] + diff(range(X[, vnames[j]])) * 0.1)
        # elipses.df = elipses.df[elipses.df$x >= limsi[1] & elipses.df$x <= limsi[2] & elipses.df$y >= limsj[1] & elipses.df$y <= limsj[2], ]
        plt <- ggplot2::ggplot(X, ggplot2::aes_(as.name(vnames[i]), as.name(vnames[j]), color = ~greedcl)) +
          ggplot2::geom_point() +
          ggplot2::geom_path(data = elipses.df, ggplot2::aes_(x = ~x, y = ~y, color = ~ factor(cl), group = ~cl), na.rm = TRUE) +
          ggplot2::scale_x_continuous(limits = limsi) +
          ggplot2::scale_y_continuous(limits = limsj) +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.title = ggplot2::element_blank()) +
          ggplot2::scale_colour_discrete("Clusters: ", guide = ggplot2::guide_legend(direction = "horizontal"))
      }
      plt + ggplot2::theme(plot.margin = ggplot2::unit(c(0.01, 0.01, 0.01, 0.01), "null"))
      if (j != (ncol(X) - 1)) {
        plt <- plt + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
      }
      if (j == 1) {
        plt <- plt + ggplot2::ggtitle(vnames[i], )
      }

      plt <- plt + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())

      plts.df[[ii]] <- plt + ggplot2::theme(plot.margin = ggplot2::unit(c(0.01, 0.01, 0.01, 0.01), "null"))
      ii <- ii + 1
    }
  }

  g_title = grid::textGrob(paste0(class(sol@model), " clustering with ", max(sol@cl), " clusters.\n"),
                     gp=grid::gpar(fontsize=curr_t$text$size,font=4))
  grob.mat <- matrix(lapply(plts.df, function(x) {
    ggplot2::ggplotGrob(x + ggplot2::theme(legend.position = "hidden"))
  }), nrow = (ncol(X) - 1))
  gt <- gtable::gtable_matrix("demo", grob.mat, ggplot2::unit(rep(1, ncol(X) - 1), "null"), ggplot2::unit(rep(1, ncol(X) - 1), "null"))
  legend <- gtable::gtable_filter(ggplot2::ggplotGrob(plts.df[[2]]), "guide-box")
  glegend <- ggplot2::ggplotGrob(ggplot2::ggplot() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white")) +
    ggplot2::annotation_custom(legend))
  gtl <- gridExtra::arrangeGrob(grobs = list(gt, glegend), nrow = 2, heights = c(10, 1))
  finished_plot <- gridExtra::grid.arrange(top = g_title, gtl)
  invisible(finished_plot)
}



ellips <- function(mu = c(0, 0), Sigma = diag(rep(1, 2)), nb_ech = 100, l = stats::qchisq(.95, df = 2)) {
  t <- seq(0, 2 * pi, len = nb_ech)
  a <- sqrt(l * eigen(Sigma)$values[1])
  b <- sqrt(l * eigen(Sigma)$values[2])
  x <- a * cos(t)
  y <- b * sin(t)
  X <- cbind(x, y)
  R <- eigen(Sigma)$vectors
  Xf <- X %*% t(R)
  data.frame(x = Xf[, 1] + mu[1], y = Xf[, 2] + mu[2])
}

block_lca <- function(sol) {
  params <- coef(sol)
  params$Thetak
  theta.df <- lapply(1:length(params$Thetak), function(v) {
    vname <- names(params$Thetak)[v]
    ccol <- params$Thetak[[v]]
    unfolded <- lapply(1:ncol(ccol), function(x) {
      data.frame(cluster = 1:nrow(ccol), prob = ccol[, x], level = colnames(ccol)[x], feature = vname, row.names = NULL)
    })
    do.call(rbind, unfolded)
  })
  plts <- lapply(theta.df, function(x) {
    ggplot2::ggplot(x) +
      ggplot2::geom_point(ggplot2::aes_(y = ~(sol@K+1-cluster), x = ~level, fill = ~prob, size = ~prob, color = ~ factor(cluster))) +
      ggplot2::scale_size_area(limits = c(0, 1)) +
      ggplot2::scale_y_continuous(limits = c(0.5, sol@K + 0.5)) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank()) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) +
      ggplot2::labs(title = x$feature[1]) +
      ggplot2::theme(legend.position = "hidden")
  })
  list_grobs <- lapply(plts, ggplot2::ggplotGrob)
  nbc <- ceiling(sqrt(length(list_grobs)))
  nbr <- ceiling(length(list_grobs) / nbc)
  nbmiss <- nbr * nbc - length(list_grobs)
  if (nbmiss > 0) {
    miss <- lapply(1:nbmiss, function(i) {
      ggplot2::ggplotGrob(ggplot2::ggplot() +
        ggplot2::theme(panel.background = ggplot2::element_blank()))
    })
    list_grobs <- c(list_grobs, miss)
  }
  grobs_matrix <- matrix(list_grobs, ncol = nbc, nrow = nbr)
  gt <- gtable::gtable_matrix("blocks", grobs_matrix, ggplot2::unit(rep(1, nbc), "null"), ggplot2::unit(rep(1, nbr), "null"))
  gridExtra::arrangeGrob(grobs = list(gt), nrow = 1, ncol = 1, heights = c(10))
}


block_gmm_marginals <- function(sol) {
  params <- coef(sol)
  liminf <- lapply(1:length(params$muk), function(k) {
    params$muk[[k]] - 4.1 * sqrt(diag(params$Sigmak[[k]]))
  })
  liminf <- apply(do.call(rbind, liminf), 2, min)
  limsup <- lapply(1:length(params$muk), function(k) {
    params$muk[[k]] + 4.1 * sqrt(diag(params$Sigmak[[k]]))
  })
  limsup <- apply(do.call(rbind, limsup), 2, max)
  nb_ech <- 2000
  plts <- list()
  for (j in 1:length(liminf)) {
    x <- seq(liminf[j], limsup[j], length.out = nb_ech)
    pdfs <- lapply(1:sol@K, function(k) {
      params$pi[k] * stats::dnorm(x, mean = params$muk[[k]][j], sd = sqrt(params$Sigmak[[k]][j, j]))
    })
    pdfs.u <- lapply(1:sol@K, function(k) {
      stats::dnorm(x, mean = params$muk[[k]][j], sd = sqrt(params$Sigmak[[k]][j, j]))
    })
    pdfs.df <- data.frame(pdf = do.call(c, pdfs), pdfu = do.call(c, pdfs.u), cl = rep(1:sol@K, each = nb_ech), x = rep(x, sol@K))
    pdfmix <- data.frame(pdf = colSums(do.call(rbind, pdfs)), x = x)
    plts[[j]] <- ggplot2::ggplot() +
      ggplot2::geom_line(data = pdfmix, ggplot2::aes_(x = ~x, y = ~pdf)) +
      ggplot2::geom_line(data = pdfs.df, ggplot2::aes_(x = ~x, y = ~pdf, color = ~ factor(cl), group = ~cl)) +
      ggplot2::scale_x_continuous(limits = c(liminf[j], limsup[j])) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank()) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) +
      ggplot2::labs(title = colnames(params$muk[[1]])[j]) +
      ggplot2::scale_color_discrete("Clusters: ", guide = ggplot2::guide_legend(direction = "horizontal", nrow = 1, label.position = "top", label.vjust = -7, keyheight = 1, title.position = "left", title.vjust = 0.2))
  }
  list_grobs <- lapply(plts, function(x) {
    ggplot2::ggplotGrob(x + ggplot2::theme(legend.position = "hidden"))
  })
  nbc <- ceiling(sqrt(length(list_grobs)))
  nbr <- ceiling(length(list_grobs) / nbc)
  nbmiss <- nbr * nbc - length(list_grobs)
  if (nbmiss > 0) {
    miss <- lapply(1:nbmiss, function(i) {
      ggplot2::ggplotGrob(ggplot2::ggplot() +
        ggplot2::theme(panel.background = ggplot2::element_blank()))
    })
    list_grobs <- c(list_grobs, miss)
  }
  grobs_matrix <- matrix(list_grobs, ncol = nbc, nrow = nbr, byrow = TRUE)
  gt <- gtable::gtable_matrix("blocks", grobs_matrix, ggplot2::unit(rep(1, nbc), "null"), ggplot2::unit(rep(1, nbr), "null"))
  legend <- gtable::gtable_filter(ggplot2::ggplotGrob(plts[[1]]), "guide-box")
  glegend <- ggplot2::ggplotGrob(ggplot2::ggplot() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white")) +
    ggplot2::annotation_custom(legend))
  gridExtra::arrangeGrob(grobs = list(gt, glegend), nrow = 2, ncol = 1, heights = c(10, 1))
}

block_gmm_marginals_violin <- function(sol) {
  params <- coef(sol)
  liminf <- lapply(1:length(params$muk), function(k) {
    params$muk[[k]] - 4.1 * sqrt(diag(params$Sigmak[[k]]))
  })
  liminf <- apply(do.call(rbind, liminf), 2, min)
  limsup <- lapply(1:length(params$muk), function(k) {
    params$muk[[k]] + 4.1 * sqrt(diag(params$Sigmak[[k]]))
  })
  limsup <- apply(do.call(rbind, limsup), 2, max)
  nb_ech <- 2000
  plts <- list()
  for (j in 1:length(liminf)) {
    x <- seq(liminf[j], limsup[j], length.out = nb_ech)
    pdfs <- lapply(1:sol@K, function(k) {
      params$pi[k] * stats::dnorm(x, mean = params$muk[[k]][j], sd = sqrt(params$Sigmak[[k]][j, j]))
    })
    pdfs.u <- lapply(1:sol@K, function(k) {
      ldn <- stats::dnorm(x, mean = params$muk[[k]][j], sd = sqrt(params$Sigmak[[k]][j, j]), log = TRUE)
      exp(ldn - max(ldn))
    })
    pdfs.df <- data.frame(pdf = do.call(c, pdfs), pdfu = do.call(c, pdfs.u), cl = rep(1:sol@K, each = nb_ech), x = rep(x, sol@K))
    pdfmix <- data.frame(pdf = colSums(do.call(rbind, pdfs)), x = x)

    mus <- data.frame(cl = 1:length(params$muk), mu = sapply(params$muk, function(p) {
      p[j]
    }))
    plts[[j]] <- ggplot2::ggplot(pdfs.df) +
      ggplot2::geom_ribbon(ggplot2::aes_(ymin = ~ cl - 0.45 * pdfu / (max(pdfu)), x = ~x, ymax = ~ cl + 0.45 * pdfu / (max(pdfu)), group = ~cl, fill = ~ factor(cl)), color = "#000000", size = 0.5) +
      ggplot2::geom_segment(data = mus, ggplot2::aes_(x = ~mu, xend = ~mu, y = ~ cl - 0.45, yend = ~ cl + 0.45)) +
      ggplot2::scale_x_continuous(limits = c(liminf[j], limsup[j])) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank()) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) +
      ggplot2::labs(title = colnames(params$muk[[1]])[j]) +
      ggplot2::scale_fill_discrete("Clusters: ", guide = ggplot2::guide_legend(direction = "horizontal", nrow = 1, label.position = "top", label.vjust = -7, keyheight = 1, title.position = "left", title.vjust = 0.2))
  }
  list_grobs <- lapply(plts, function(x) {
    ggplot2::ggplotGrob(x + ggplot2::theme(legend.position = "hidden"))
  })
  nbc <- ceiling(sqrt(length(list_grobs)))
  nbr <- ceiling(length(list_grobs) / nbc)
  nbmiss <- nbr * nbc - length(list_grobs)
  if (nbmiss > 0) {
    miss <- lapply(1:nbmiss, function(i) {
      ggplot2::ggplotGrob(ggplot2::ggplot() +
        ggplot2::theme(panel.background = ggplot2::element_blank()))
    })
    list_grobs <- c(list_grobs, miss)
  }
  grobs_matrix <- matrix(list_grobs, ncol = nbc, nrow = nbr)
  gt <- gtable::gtable_matrix("blocks", grobs_matrix, ggplot2::unit(rep(1, nbc), "null"), ggplot2::unit(rep(1, nbr), "null"))
  legend <- gtable::gtable_filter(ggplot2::ggplotGrob(plts[[1]]), "guide-box")
  glegend <- ggplot2::ggplotGrob(ggplot2::ggplot() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white")) +
    ggplot2::annotation_custom(legend))
  gridExtra::arrangeGrob(grobs = list(gt, glegend), nrow = 2, ncol = 1, heights = c(10, 1))
}

block_mmm <- function(sol) {
  gt_cat <- block_lca(sol)
  gt_cont <- block_gmm_marginals_violin(sol)
  list_grobs <- c(gt_cat$grobs[[1]]$grobs[1:length(sol@obs_stats$lca)], gt_cont$grobs[[1]]$grobs[1:ncol(sol@obs_stats$gmm$cluster1$m)])
  nbc <- ceiling(sqrt(length(list_grobs)))
  nbr <- ceiling(length(list_grobs) / nbc)
  nbmiss <- nbr * nbc - length(list_grobs)
  if (nbmiss > 0) {
    miss <- lapply(1:nbmiss, function(i) {
      ggplot2::ggplotGrob(ggplot2::ggplot() +
        ggplot2::theme(panel.background = ggplot2::element_blank()))
    })
    list_grobs <- c(list_grobs, miss)
  }
  grobs_matrix <- matrix(list_grobs, ncol = nbc, nrow = nbr, byrow = TRUE)
  gt <- gtable::gtable_matrix("blocks", grobs_matrix, ggplot2::unit(rep(1, nbc), "null"), ggplot2::unit(rep(1, nbr), "null"))
  glegend <- gt_cont$grobs[[2]]
  gridExtra::arrangeGrob(grobs = list(gt, glegend), nrow = 2, ncol = 1, heights = c(10, 1))
}

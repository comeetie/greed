library(readr)
edges <- read_csv("./data-raw/new_guinea_edges.csv")
names(edges) <- c("source", "target", "w")
net <- array(0, c(16, 16, 3))
edges$source <- edges$source + 1
edges$target <- edges$target + 1
edges_neg <- edges[edges$w == -1, ]
edges_pos <- edges[edges$w == 1, ]

net[cbind(edges_neg$source, edges_neg$target, 1)] <- 1
net[cbind(edges_neg$target, edges_neg$source, 1)] <- 1
net[cbind(edges_pos$source, edges_pos$target, 2)] <- 1
net[cbind(edges_pos$target, edges_pos$source, 2)] <- 1
net[, , 3] <- 1 - colSums(aperm(net, c(3, 1, 2)), 1)

library(greed)

sol <- greed(net, K = 5)

mod <- MixedModels(list(pos = DcSbmPrior(), neg = DcSbmPrior()))
data <- list(pos = net[, , 2], neg = net[, , 1])
sol_bidcsbm <- greed(data, model = mod, K = 5)

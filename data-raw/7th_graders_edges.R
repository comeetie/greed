library(readr)
edges <- read_csv("./data-raw/7th_graders.csv")
names(edges) <- c("source", "target", "w", "layer")
edges$source <- edges$source + 1
edges$target <- edges$target + 1
N <- max(edges$target)
net <- array(0, c(N, N, 3))
net[cbind(edges$source, edges$target, edges$layer)] <- 1

library(greed)

sol <- greed(net, K = 5)

mod <- MixedModels(list(class = DcSbmPrior(), friends = DcSbmPrior(), work = DcSbmPrior()))
data <- list(class = net[, , 1], friends = net[, , 2], work = net[, , 3])
sol <- greed(data, model = mod, K = 5)

plot(extractSubModel(sol, "work"), type = "blocks")

colSums(aperm(net, c(3, 1, 2)), 1)

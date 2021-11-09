X = read.csv('./agaricus-lepiota.data', sep=',', header = F)
X = data.frame(lapply(X, function(x) factor(x)))
mushroom = data.frame(edibility = X[,1])
mushroom = cbind(mushroom, X[,-1])

usethis::use_data(mushroom)
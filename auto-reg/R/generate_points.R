library(mvtnorm)
library(lhs)
set.seed(1)
source("R/helper.R")
N = 10000
x = randomLHS(N, 12)

## Priors from the Milner, 2012 paper
## c1->c11, U(-5, 1) on log scale
## c12, U(-12, -6) on log scale

for(i in 1:11) 
  x[,i] = exp(scale_lhs(x[,i], -5, -5, 1))
x[,12] = exp(scale_lhs(x[,12], -12, -12, -6))
x = cbind(1:nrow(x), x)

write.table(x, file="auto-reg/data/prior.csv", col.names=FALSE, row.names=FALSE, sep=",")


## Posterior from Milner paper
# post = read.csv("../data/post.csv", header=TRUE)
# path_G = round(read.csv("../data/pathG.csv", header=TRUE))
# path_I = round(read.csv("../data/pathI.csv", header=TRUE))
# 
# complete = cbind(post, path_G, path_I)
# write.csv(complete, file="../data/samples.csv", col.names=FALSE, row.names=FALSE)
# 
# max(cor(post[,2:13]))
# 
# sim = read.csv("../data/sim.csv", header=FALSE)
# 
# 

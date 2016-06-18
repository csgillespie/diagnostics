library("deSolve")
library("lhs")
library("issb")
library("Rcpp")
library("reshape2")
source("R/helper.R")
source("R/lv_functions.R")

###################################################################
## Construct LHD
###################################################################
set.seed(1)
N = 100
x = randomLHS(N,3)
x[,1] = 10^scale_lhs(x[,1], -6, -6, 0)
x[,2] = 10^scale_lhs(x[,2], -6, -6, 0)
x[,3] = scale_lhs(x[,2], 0, 0, 0.0) + 0.3

saveRDS(x, file="data/lv_lhs.RData")

###################################################################
## Parameter scan over LHD
###################################################################

## Moment closure
sim_mc_10 = lv_simulate(x, maxtime=10, is_mc=TRUE)
sim_mc_30 = lv_simulate(x, maxtime=30, is_mc=TRUE)
sim_mc_100 = lv_simulate(x, maxtime=100, is_mc=TRUE)
saveRDS(sim_mc_10, file="data/lv_mc_10.RData")
saveRDS(sim_mc_30, file="data/lv_mc_30.RData")
saveRDS(sim_mc_100, file="data/lv_mc_100.RData")

sim_lna_10 = lv_simulate(x, maxtime=10, is_mc=FALSE)
sim_lna_30 = lv_simulate(x, maxtime=30, is_mc=FALSE)
sim_lna_100 = lv_simulate(x, maxtime=100, is_mc=FALSE)
saveRDS(sim_lna_10, file="data/lv_lna_10.RData")
saveRDS(sim_lna_30, file="data/lv_lna_30.RData")
saveRDS(sim_lna_100, file="data/lv_lna_100.RData")

###################################################################
## Get Simulations
###################################################################
## Form the model
Rcpp::sourceCpp("src/lv.cpp")
i = 3
maxtime = 100
no_of_sims = 1000
probs = matrix(0, ncol=2, nrow(x))
for(i in 1:nrow(x)) {
  theta = x[i,]; p=0
  for(j in 1:no_of_sims) {
    ext = is_extinct(theta[1], theta[2], theta[3], maxtime)  
    probs[i,] = probs[i,] + ext
  }
  message(i)
}
saveRDS(probs, file=paste0("data/lv_extinction_t", maxtime, ".RData"))







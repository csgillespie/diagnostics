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
sim_mc_30 = lv_simulate(x, maxtime=30, is_mc=TRUE)
sim_mc_100 = lv_simulate(x, maxtime=100, is_mc=TRUE)
saveRDS(sim_mc_30, file="data/lv_mc_30.RData")
saveRDS(sim_mc_100, file="data/lv_mc_100.RData")

sim_lna_30 = lv_simulate(x, maxtime=30, is_mc=FALSE)
sim_lna_100 = lv_simulate(x, maxtime=100, is_mc=FALSE)
saveRDS(sim_lna_30, file="data/lv_lna_30.RData")
saveRDS(sim_lna_100, file="data/lv_lna_100.RData")

#x = x[6,, drop=FALSE]
# mc = lv_simulate(x[6,,drop=F], maxtime=1, is_mc=FALSE)
# mc



###################################################################
## Get Simulations
###################################################################
## Form the model
source("models/lv.R")
Rcpp::sourceCpp("src/lv.cpp")
i = 3
maxtime = 100
no_of_sims = 1000
probs = matrix(0, ncol=3, nrow(x))
for(i in 1:nrow(x)) {
  theta = x[i,]; p=0
  for(j in 1:no_of_sims) {
    ext = is_extinct(theta[1], theta[2], theta[3], maxtime)  
    probs[i,ext] = probs[i,ext] + 1
  }
  message(i)
}
#probs_30 = probs
probs_100 = probs
## 44, 79
#saveRDS(probs_30, file="data/lv_extinction_t30.RData")
ms
?issb::multiple_sims
theta = x[3,]

theta[1:2] = c(10^-6, 10^-6)
model$get_pars(theta)
model$get_initial(c(100, 100))
g = gillespie(model, maxtime=100)
head(g)
setnicepar(mfrow=c(1, 2))
plot(g[,1], g[,2], type="l")
plot(g[,1], g[,3], type="l")

tail(g)


ms
d_ms = as.data.frame(ms)
d_ms_sub = melt(d_ms, c("sim_no", "Time"))
saveRDS(d_ms_sub, file="data/lv_ms.RData")





library(lhs)
library(issb)
source("R/helper.R")

###################################################################
## Load model 
###################################################################
source("models/SIR.R")
(theta = model$get_pars())
(initial = model$get_initial())
maxtime = 50

## Quick Simulation
d = gillespie(model, 50)
plot(d[,1], d[,2], type="l")
plot(d[,1], d[,3], type="l")


###################################################################
## LHS
###################################################################
set.seed(1)
N = 1000
x = randomLHS(N, 4)
x[,1] = scale_lhs(x[,1], 0, 0, 2)
x[,2] = scale_lhs(x[,2], 0, 0, 2)
x[,3] = scale_lhs(x[,3], 0, 0, 2)
x[,4] = scale_lhs(x[,4], 0, 0, 0.04)
saveRDS(x, file="data/sir_lhs.RData")

###################################################################
## Direct method
###################################################################
set.seed(1)
sim_gil = matrix(0, ncol=2, nrow=NROW(x))
for(i in seq_along(x[,1])) {
  theta[1:ncol(x)] = x[i,]
  model$get_initial(initial)
  model$get_pars(theta)
  sim_gil[i,] = tail(gillespie(model, maxtime=maxtime, tstep=maxtime), 1)[2:3]
  message(i)
}
saveRDS(sim_gil, "data/sir_gil.RData")

###################################################################
##  CLE/Diffusion 
###################################################################
set.seed(2)
sim_diff= matrix(0, ncol=2, nrow=NROW(x))
for(i in seq_along(x[,1])) {
  theta[1:ncol(x)] = x[i,]
  model$get_initial(initial)
  model$get_pars(theta)
  sim_diff[i,] = tail(diffusion(model, maxtime=maxtime,  ddt=0.05), 1)[2:3]
  message(i)
}
saveRDS(sim_diff, "data/sir_diff.RData")

###################################################################
## tau-leap
###################################################################
set.seed(3)
sim_tau = matrix(0, ncol=2, nrow=NROW(x))
for(i in seq_along(x[,1])) {
  theta[1:ncol(x)] = x[i,]
  model$get_initial(initial)
  model$get_pars(theta)
  sim_tau[i,] = tail(tau_leap(model, maxtime=maxtime), 1)[2:3]
  message(i)
}
saveRDS(sim_tau, "data/sir_tau.RData")


###################################################################
## Looking at the max difference
###################################################################
set.seed(4)
d = sim_gil[,1] - sim_diff[,1]
(theta = x[which.max(d), ])
model$get_initial(initial)
model$get_pars(theta)
ms = multiple_sims(model, 
              simulator=gillespie, 
              maxtime=maxtime, 
              tstep=0.01, 
              no_sims=50, 
              no_cores=4)
ms = as.data.frame(ms)

ms_diff = multiple_sims(model, 
                       simulator=diffusion, 
                       ddt = 0.05,
                       maxtime=maxtime, 
                       tstep=0.01, 
                       no_sims=50, 
                       no_cores=4)
ms_diff = as.data.frame(ms_diff)

saveRDS(ms, "data/sir_mult_gil.RData")
saveRDS(ms_diff, "data/sir_mult_diff.RData")


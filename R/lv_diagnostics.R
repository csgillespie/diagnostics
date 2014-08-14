library(deSolve)
library(lhs)
library(issb)
source("R/helper.R")

maxtime = 30

###################################################################
## Moment closure LV equations. 
## Additional parameter is_mc used to switch to LNA
###################################################################
lvmodel = function(t, x, parms) {
  with(as.list(c(x, parms)),{
    
    k11 = k11*is_mc
    ##Means - notice the covariance term
    dk10 = c1*k10 - c2*k10*k01 - c2*k11
    dk01 = c2*k10*k01 -  c3*k01 + c2*k11
    
    ##Variances and covariance with the Normal approximation
    dk20 = c2*k11 + 2*k20*(c1 -c2*k01) +
      k10*(c1 -2*c2*k11 + c2*k01)
    
    dk11 = k11*(c1- c2*(1 - k01+ k10) - c3) +
      c2*k20*k01 - c2*k10*(k02+k01)
    
    dk02 = c2*k11  + 2*k02*(c2*k10 - c3)  +
      k01*(2*c2*k11 + c3 + c2*k10)
    
    res = c(dk10, dk01, dk20, dk11, dk02)
    list(res)
  })
}


###################################################################
## Form LHD
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
sim_mc = matrix(0, ncol=6, nrow=NROW(x))
for(i in seq_along(x[,1])) {
  pars = c(c1 = x[i,1], c2 = x[i,2], c3 = x[i,3], is_mc = TRUE)
  yini = c(k10=100, k01=100, k20=0, k11=0, k02=0)
  times = seq(0, maxtime, length.out=2)
  
  ## At extreme points the MC approximation breaks down
  ## Catch error and do something sensible
  out = try(ode(yini, times, lvmodel, pars), silent=TRUE)
  if(inherits(out, "try-error")) 
    sim_mc[i,] = rep(NA_real_, 6)
  else 
    sim_mc[i,] = out[2,]
  message(i)
}
saveRDS(sim_mc, file="data/lv_mc.RData")

## LNA
sim_lna = matrix(0, ncol=6, nrow=NROW(x))
for(i in seq_along(x[,1])) {
  pars = c(c1 = x[i,1], c2 = x[i,2], c3 = x[i,3], is_mc = FALSE)
  yini = c(k10=100, k01=100, k20=0, k11=0, k02=0)
  times <- seq(0, maxtime, length.out=2)
  out   <- try(ode(yini, times, lvmodel, pars), silent=TRUE)
  if(inherits(out, "try-error")) 
    sim_lna[i,] = rep(NA_real_, 6)
  else 
    sim_lna[i,] = out[2,]
  message(i)
}
saveRDS(sim_lna, file="data/lv_lna.RData")

###################################################################
## Get Simulations
###################################################################
## Form the model
source("models/lv.R")
theta = c(10^-4, 0.1, 0.30000)
model$get_pars(theta)
model$get_initial(c(100, 100))

## Do a few simulations
ms = multiple_sims(model, 
                   simulator=gillespie, 
                   maxtime=maxtime, 
                   tstep=0.01, 
                   no_sims=50, 
                   no_cores=6)


d_ms = as.data.frame(ms)
d_ms_sub = melt(d_ms, c("sim_no", "Time"))
saveRDS(d_ms_sub, file="data/lv_ms.RData")





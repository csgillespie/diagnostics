library(lhs)
library(issb)
source("R/helper.R")

N = 10000
if(N == 1000) {dir = "data/"
} else {dir = paste0("data/N_", N, "_")}

###################################################################
## Load model 
###################################################################
source("models/schlogl.R")
maxtime = 5.0
(theta = model$get_pars())
(initial = model$get_initial())

###################################################################
## LHS
###################################################################
set.seed(1)
x = randomLHS(N, 2)
x[,1] = 10^scale_lhs(x[,1], -2, -2, 1) ## This is c4
x[,2] = 10^scale_lhs(x[,2], -4, -4, -2) ## This is c3
saveRDS(x, file=paste0(dir, "schlogl_lhs.RData"))

###################################################################
## Direct method
###################################################################
set.seed(2)
sim_gil = matrix(0, ncol=3, nrow=NROW(x))
for(i in seq_along(x[,1])) {
  set.seed(i)
  theta[4:3]= x[i,]
  model$get_initial(initial)
  model$get_pars(theta)
  sim_gil[i,] = tail(gillespie(model, maxtime=maxtime, tstep=maxtime), 1)[2:4]
  message(i)
}
saveRDS(sim_gil, file=paste0(dir, "schlogl_gil.RData"))

###################################################################
## LNA diagnostics
###################################################################
l = list(NROW(x))
z = model$get_initial()
m = z*0
for(i in seq_along(x[,1])) {
  theta[4:3] = x[i,]
  model$get_initial(initial)
  model$get_pars(theta)
  out = issb:::lnastep(model, z, m, maxtime)
  l[[i]] = list(m = out$z, v = out$V)
  message(i)
}

res = matrix(0, ncol=3, nrow=NROW(x))
mah = matrix(0, ncol=1, nrow=NROW(x))
for(i in seq_along(x[,1])) {
  res[i,] = (sim_gil[i, ] - l[[i]]$m)/sqrt(diag(l[[i]]$v))
  mah[i,] = t(sim_gil[i, ] - l[[i]]$m) %*%
    MASS::ginv(l[[i]]$v) %*% (sim_gil[i, ] - l[[i]]$m)
  message(i)
}
saveRDS(l, file=paste0(dir, "schlogl_lna.RData"))
saveRDS(res, file=paste0(dir, "schlogl_residuals.RData"))

###################################################################
## Get simulations for largest residual
###################################################################
my_res = rev(sort(abs(res[,1])))[1]

(max_res = which.max(abs(res[,1])))
#(max_res = which.min(res[,1]-my_res))

theta[4:3] = x[max_res,]
model$get_initial(initial)
model$get_pars(theta)

ms = multiple_sims(model, 
                   simulator=gillespie, 
                   maxtime=maxtime, 
                   tstep=0.01, 
                   no_sims=50, 
                   no_cores=6)
model$get_stoic()

d_ms = as.data.frame(ms)
saveRDS(d_ms, file=paste0(dir, "schlogl_res_sims.RData"))


theta[4:3] = x[max_res,]
model$get_initial(initial)
model$get_pars(theta)
out = lna(model, maxtime, ddt = 0.001, noise=FALSE)

saveRDS(as.data.frame(out), file=paste0(dir, "schlogl_res_lna_sims.RData"))

###################################################################
## Get simulations for largest error bar
###################################################################
#(max_res = which.max(abs(res[,1])))
theta[4:3] = x[max_res,]
model$get_initial(initial)
model$get_pars(theta)

ms = multiple_sims(model, 
                   simulator=gillespie, 
                   maxtime=maxtime, 
                   tstep=0.01, 
                   no_sims=50, 
                   no_cores=6)
d_ms = as.data.frame(ms)
saveRDS(d_ms, file=paste0(dir, "schlogl_err_sims.RData"))


theta[4:3] = x[max_res,]
model$get_initial(initial)
model$get_pars(theta)

out = lna(model, maxtime, ddt = 0.001, noise=FALSE)
saveRDS(as.data.frame(out), file=paste0(dir, "schlogl_err_lna_sims.RData"))












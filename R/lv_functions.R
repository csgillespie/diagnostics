###################################################################
## Moment closure LV equations. 
## Additional parameter is_mc used to switch to LNA
###################################################################
lvmodel = function(t, x, parms) {
  with(as.list(c(x, parms)),{
    
    ##Means - notice the covariance term
    dk10 = c1*k10 - c2*k10*k01 - c2*k11*is_mc
    dk01 = c2*k10*k01 -  c3*k01 + c2*k11*is_mc
    
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

lv_simulate = function(x, maxtime, is_mc=TRUE){
  results = matrix(0, ncol=6, nrow=NROW(x))
  i = 1
  for(i in seq_along(x[,1])) {
    pars = c(c1 = x[i,1], c2 = x[i,2], c3 = x[i,3], is_mc = is_mc)
    yini = c(k10=100, k01=100, k20=0, k11=0, k02=0)
    times = seq(0, maxtime, length.out=2)
    
    ## At extreme points the MC approximation breaks down
    ## Catch error and do something sensible
    (out = try(ode(yini, times, lvmodel, pars, method="impAdams"), silent=TRUE))
    if(inherits(out, "try-error")) 
      results[i,] = rep(NA_real_, 6)
    else 
      results[i,] = out[2,]
    message(i)
  }
  return(results)
}




# library("deSolve")
# i = 96
#  is_mc = T#FALSE
#  maxtime = 100
#  pars = c(c1 = theta[1], c2 = theta[2], c3 = theta[3], is_mc = is_mc)
# yini = c(k10=100, k01=100, k20=0, k11=0, k02=0)
# times = seq(0, maxtime, length.out=1000)
# which(x[,1] > 0.1 & x[,2] < 10^{-4})
## At extreme points the MC approximation breaks down
## Catch error and do something sensible

library("deSolve")
compiler::enableJIT(3)
###################################################################
## Moment closure and LNA equations. 
## Additional parameter is_mc used to switch to LNA
###################################################################
epi_model = function(t, x, parms) {
  with(as.list(c(x, parms)),{
    #alpha, delta,  mu, beta
    ##Means - notice the covariance term
    dk10 = -beta*(k01*k10 + k11*is_mc) + alpha
    dk01 = beta*(k01*k10+k11*is_mc) - mu*k01 + delta
    
    ##Variances and covariance with the Normal approximation
    dk20 = beta*k01*(k10-2*k20) + beta*k11*(1-2*k10) + alpha
    dk11 = beta*k11*(k10-k01-1) - beta*(k02+k01)*k10 + beta*k01*k20 - mu*k11
    dk02 = beta*k10*(k01 + 2*k02) + beta*k11*(1 + 2*k01) + mu*(k01 - 2*k02) + delta
    
    res = c(dk10, dk01, dk20, dk11, dk02)
    list(res)
  })
}
maxtime=100
times = seq(0, maxtime, length.out=1000)

yini = c(k10=0, k01=0, k20=0, k11=0, k02=0)

pars = c(alpha=10, mu=0.7, delta=1, beta=0.008, is_mc = TRUE)
od_mc = ode(yini, times, epi_model, pars)
pars["is_mc"] = FALSE
od_lna = ode(yini, times, epi_model, pars)

setnicepar(mfrow=c(1,1))
plot(od_mc[,1], od_mc[,2], type="l", col=1)
lines(od_mc[,1], od_mc[,3], col=2)
lines(od_lna[,1], od_lna[,2], col=1, lty=3)
lines(od_lna[,1], od_lna[,3], col=2, lty=3)


N = 10000

x = lhs::randomLHS(N, 2)
x[,1] = 10^scale_lhs(x[,1], -4, -4, 1)
x[,2] = 10^scale_lhs(x[,2], -4, -4, 1)

alpha= x[,1]
beta = x[,2]
dd = data.frame(alpha=alpha, beta=beta, lna10 = numeric(N), lna01 = numeric(N), 
                lna20 = numeric(N), lna11 = numeric(N), lna02 = numeric(N),  
                mc10 = numeric(N), mc01 = numeric(N), mc20 = numeric(N), 
                mc11 = numeric(N), mc02 = numeric(N))
times = seq(0, maxtime, length.out=2)

set.seed(1)

i = 1; j = 1; k = 1
for(i in seq_along(alpha)) {
  for(j in seq_along(beta)) {
    pars = c(alpha=alpha[i], mu=0.7, delta=1, beta=beta[j], is_mc = TRUE)
    od_mc = ode(yini, times, epi_model, pars)
    pars["is_mc"] = FALSE
    od_lna = ode(yini, times, epi_model, pars)
    dd[k,] = c(beta[j], alpha[i], od_lna[2,-1], od_mc[2,-1])
    k = k + 1
  }
  message(i)
}

dd_tmp = dd
dd = dd_tmp
library(ggplot2)

dd$diff = abs(dd$lna10 - dd$mc10)/sd(dd$lna20)
dd$diff[dd$diff > 0.5] = 1
dd$diff[dd$diff < 0.5] = 0
dd[is.nan(dd$diff), ]$diff = 1
ggplot(dd, aes(alpha, beta, colour=factor(diff))) + 
  geom_point() +  scale_y_log10() + scale_x_log10()




setnicepar(mfrow=c(2, 2))
plot(dd$lna10, type="l")
lines(dd$mc10, col=2)
plot(abs(dd$lna10 - dd$mc10)/sd(dd$lna20), type="l")

plot(dd$lna01, type="l")
lines(dd$mc01, col=2)
plot(abs(dd$lna01 - dd$mc01)/sd(dd$lna02), type="l")

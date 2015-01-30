library(ggplot2)
library(grid)
source("R/helper.R")
theme_set(theme_bw())
dir = "graphics/"
mypalette()


###################################################################
## Read in data generate from schlogl_diagnostics.R
###################################################################
x = readRDS(file="data/schlogl_lhs.RData")
l = readRDS(file="data/schlogl_lna.RData")
res = readRDS(file="data/schlogl_residuals.RData")
sim_gil = readRDS(file="data/schlogl_gil.RData")
out = readRDS(file="data/schlogl_res_lna_sims.RData")
d_ms = readRDS(file="data/schlogl_res_sims.RData")



#####################################################
## Schlogl latin hyper cube
#####################################################
g0 = ggplot(as.data.frame(x), aes(V2, V1)) + 
  geom_point(col="black", size=0.7, pch=21, bg=2) + 
  annotation_logticks() + 
  xlab(expression(c[3])) + 
  ylab(expression(c[4])) 


g = g0 + scale_x_continuous(trans="log10", 
                           breaks=c(10^{-4}, 10^{-3}, 10^{-2}),
                           labels=c(expression(10^{-4}), expression(10^{-3}), expression(10^{-2}))) + 
  scale_y_continuous(trans="log10", 
                     breaks=c(10^{-2}, 10^{-1}, 10^{0}, 10^{1}), 
                     labels=c(expression(10^{-2}), expression(10^{-1}), expression(10^{0}), expression(10^{1}))) + 
  annotate("text", x=10^{-4}, 10, label="(a)", size=4, vjust = 1)
#g

###################################################################
## Plot the LNA diagnostics
###################################################################

## Manipulate data for ggplot2
cut_off = data.frame(y = c(-qnorm(0.975), -qnorm(0.9975), qnorm(0.975), qnorm(0.9975)), 
                     threshold = factor(1:2))
dd = data.frame(x = x[,1], y = res[,1])

## Plot
g1 = ggplot(dd, aes(x, y)) + 
  geom_point() + 
  ylim(c(-10, 5)) + 
  ylab(expression(e[i]^"*")) + 
  xlab(expression(c[4]))
#g1
g1 = g1 + 
  geom_hline(data=cut_off,aes(yintercept=y, lty=threshold), 
                     col="grey60") 
g2 = g1 + scale_x_continuous(trans="log10", 
                               breaks=c(10^{-2}, 10^{-1}, 10^{0}, 10^1), limits=c(10^{-2}, 10^1), 
                               labels=c(expression(10^{-2}), expression(10^{-1}), 
                                        expression(10^{0}), expression(10^{1}))) + 
  stat_smooth(se=FALSE) +
  annotate("text", x=10^{-2}, 5, label="(b)", size=4)



################################################################
#Stochastic simulation for large residual
################################################################
(max_res = which.max(abs(res[,1])))
d_ms_sub = d_ms[d_ms$sim_no < 51,]
g6 = ggplot(d_ms_sub) + 
  geom_step(aes(Time, X, group=sim_no), alpha=0.15) + 
  xlab("Time") + 
  ylab("Population") + 
  ylim(c(0, 800))
#g6
g7 = g6 + geom_line(data=out, aes(Time, X), col=4)
g8 = g7 + annotate("text", x=1.2, y = 250, label="LNA approximation", col=4, size=4) + 
  annotate("text", x=0, 800, label="(d)", size=4)



################################################################
## Credible region plot
################################################################


dd_l = data.frame(lower=0, upper=numeric(nrow(x)))
for(i in seq_along(x[,1])) 
  dd_l[i,] = l[[i]]$m[1] + c(-1.96, 1.96)*sqrt(diag(l[[i]]$v))[1]
dd_l$x = x[,2]
dd_l$sim = sim_gil[,1]
dd_l$out = dd_l$sim > dd_l$upper | dd_l$sim < dd_l$lower


g9 = ggplot(dd_l, aes(x=x)) + 
  geom_linerange(aes(ymin=lower, ymax=upper), size=0.1) + 
  geom_point(aes(y=sim, colour=out), size=1) + 
  scale_color_manual(values = c("black", 2), guide=FALSE) + 
  xlab(expression(c[3])) + 
  scale_x_continuous(trans="log10", 
                     breaks=c(10^{-4}, 10^{-3}, 10^{-2}),
                     labels=c(expression(10^{-4}), expression(10^{-3}), expression(10^{-2}))) +
  ylim(c(-50, 1250)) + 
  ylab("Population") + 
  annotation_logticks(sides="b") 
g10 = g9 + annotate("text", x=10^{-4}, 1200, label="(c)", size=4)
  

###################################################################
## Save plots as jpg

fname = paste0(dir, "figure2a.jpg")
jpeg(fname, width=4*resol, height=4*resol, res=resol)
print(g)
dev.off()

fname = paste0(dir, "figure2b.jpg")
jpeg(fname, width=4*resol, height=4*resol, res=resol)
print(g2)
dev.off()

fname = paste0(dir, "figure2c.jpg")
jpeg(fname, width=4*resol, height=4*resol, res=resol)
print(g8)
dev.off()

fname = paste0(dir, "figure2d.jpg")
jpeg(fname, width=4*resol, height=4*resol, res=resol)
print(g10)
dev.off()


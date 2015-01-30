library(gridExtra)
library(ggplot2)
source("R/helper.R")
dir = "graphics/"
theme_set(theme_bw())
mypalette()

###################################################################
## Read in data generate from lv_diagnostics.R
###################################################################
x = readRDS(file="data/sir_lhs.RData")
sim_gil = readRDS(file="data/sir_gil.RData")
sim_tau = readRDS(file="data/sir_tau.RData")
sim_diff = readRDS(file="data/sir_diff.RData")
ms = readRDS("data/sir_mult_gil.RData")
ms_diff = readRDS("data/sir_mult_diff.RData")

###################################################################
#qq-plot
###################################################################
## Tau
par(mfrow=c(2,2))
q = qqplot(sim_gil[,1], sim_tau[,1]);abline(0, 1)
dd = data.frame(x = q$x, y=q$y)
q = qqplot(sim_gil[,2], sim_tau[,2]);abline(0, 1)
tmp = data.frame(x = q$x, y=q$y); dd = rbind(dd, tmp)


## CLE
q = qqplot(sim_gil[,1], sim_diff[,1]);abline(0, 1)
tmp = data.frame(x = q$x, y=q$y); dd = rbind(dd, tmp)
q = qqplot(sim_gil[,2], sim_diff[,2]); abline(0, 1)
tmp = data.frame(x = q$x, y=q$y); dd = rbind(dd, tmp)

dd$type = rep(c("tau-leap", "CLE"), each=2*nrow(x))
dd$species = rep(c("X[1]", "X[2]"), each=nrow(x))

dd_line = data.frame(x=c(0, 150), y=c(0, 150), 
                     type=rep(c("tau-leap", "CLE"), each=4), 
                     species = rep(c("X[1]", "X[2]"), each=2))


g1 = ggplot(dd, aes(x, y)) + 
  geom_point(colour="grey40") +
  facet_grid(species ~ type, labeller=label_parsed) + 
  geom_line(data=dd_line, aes(x, y), size=0.5, alpha=0.5, col=4) + 
  xlim(c(0, 150)) + ylim(c(0, 150)) + 
  xlab("Species level") + 
  ylab("Species level")
#g1


################################################################
# Difference
################################################################
dd_dif = data.frame(x=x[,2], 
                    y= c(sim_gil[,1] - sim_tau[,1], sim_gil[,1] - sim_diff[,1]), 
                    type=rep(c("tau-leap", "CLE"), each=nrow(sim_gil)))
g2 = ggplot(dd_dif, aes(x, y)) + 
  geom_point(size=1) + stat_smooth(se = FALSE, lwd=1.5, col=3) + 
  xlab(expression(c[2])) + 
  ylab("Difference") + 
  ylim(c(-60, 100)) + 
  facet_grid(~type, labeller=label_parsed)
#g2
################################################################
#Stochastic simulation plots
################################################################
mm = rbind(ms, ms_diff)
mm$type = rep(c("Direct", "CLE"), each=nrow(ms))
head(mm)
m1 = mm[,c(1:3, 5)]; colnames(m1)[3] = "pop"; m1$species = "X[1]"
m2 = mm[,c(1:2, 4:5)]; colnames(m2)[3] = "pop"; m2$species = "X[2]"
m = rbind(m1, m2)

m = m[m$sim_no < 6,]
g3 = ggplot(m) + 
  geom_step(aes(Time, pop, group=sim_no), alpha=0.20) + 
  xlab("Time") + ylab("Population") + 
  facet_grid(type ~ species, scale="free_y", labeller=label_parsed) + 
  ylim(c(0, 90))
#g3

################################################################
## Save plots
################################################################

fname = paste0(dir, "sir_f1.pdf")
pdf(fname, width=6, height=6)
print(g1)
dev.off()
crop_plot(fname)

fname = paste0(dir, "sir_f2.pdf")
pdf(fname, width=6, height=4)
print(g2)
dev.off()
crop_plot(fname)

fname = paste0(dir, "sir_f3.pdf")
pdf(fname, width=6, height=6)
print(g3)
dev.off()
crop_plot(fname)



fname = paste0(dir, "sir_f1.jpg")
jpeg(fname, width=6*resol, height=6*resol, res=resol)
print(g1)
dev.off()


fname = paste0(dir, "sir_f2.jpg")
jpeg(fname, width=6*resol, height=4*resol, res=resol)
print(g2)
dev.off()


fname = paste0(dir, "sir_f3.jpg")
jpeg(fname, width=6*resol, height=6*resol, res=resol)
print(g3)
dev.off()




























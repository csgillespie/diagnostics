library(ggplot2)
library(grid)
theme_set(theme_bw())
source("R/helper.R")
mypalette()

## Stochastic simulation data
## Load in 
a = read.csv("auto-reg/data/sim.csv", head=FALSE)
head(a);tail(a)

## Rename columns and manipulate
colnames(a) = c("Time", "r_g", "r_i", "g", "i", "G", "I")
a$Time = a$Time - 1

df = data.frame(values = c(a$r_g, a$r_i, a$g, a$i))
df$Time = a$Time
df$Species = as.factor(rep(c("r_g", "r_i", "g", "i"), each=50))

## Plots of the species
g1 = ggplot(data=df, aes(x=Time, y=values, group=Species)) +
  ylab("Population") +
  xlim(c(0, 50)) + ylim(c(0, 30)) +
  geom_point(size=1.5) +
  theme(axis.title.x = element_text(vjust = -0.5)) 
  
g1c = g1 + geom_line(aes(colour=Species)) + scale_color_manual(values=2:5)
#g1_bw = g1 + geom_line(aes(lty=Species))
#g1

g1b_w = ggplot(data=df, aes(x=Time, y=values, group=Species)) +
  geom_line(aes(lty=Species)) +
  ylab("Population") +
  xlim(c(0, 50)) + ylim(c(0, 30)) +
  geom_point(size=1.5) +
  theme(axis.title.x = element_text(vjust = -0.5))
#g1

g2 = ggplot(data =a, aes(x=Time, y=I)) +
  geom_line() + xlim(c(0, 50)) + ylim(c(0, 15)) +
  geom_point(size=1.5) +
  theme(axis.title.x = element_text(vjust = -0.5))
#g2 

g3 = ggplot(data =a, aes(x=Time, y=G)) +
  geom_line() + xlim(c(0, 50)) +
  geom_point(size=1.5) +
  theme(axis.title.x = element_text(vjust = -0.5))


resol = 1000
jpeg("graphics/auto1.jpg",  width=8*resol, height = 8*resol, res = resol)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(g1c, vp = vplayout(1, 1:2))
print(g2, vp = vplayout(2, 1))
print(g3, vp = vplayout(2, 2))
dev.off()


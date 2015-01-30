library(lhs)
library(ggplot2)
source("R/helper.R")
dir = "graphics/"
theme_set(theme_bw())
mypalette()

###########################################################
#Example LHD - figure 1
###########################################################
## Generate data
set.seed(100)
N = 50
x = as.data.frame(randomLHS(N,2))

## Plot
g = ggplot(x, aes(V1, V2)) + 
  geom_point(col="black", size=4, pch=21, bg=2) + 
  geom_rug() + 
  scale_x_continuous(breaks=seq(0, 1, 0.1), 
                     limits=c(0, 1), labels=c(0, "", "", "", "", 0.5, "", "", "", "", 1.0))+
  scale_y_continuous(breaks=seq(0, 1, 0.1), 
                     limits=c(0, 1), labels=c(0, "", "", "", "", 0.5, "", "", "", "", 1.0)) +
  xlab(expression(c[1])) + ylab(expression(c[2]))
g


###########################################################
## Save graph 
###########################################################
fname = paste0(dir, "lhs_f1.pdf")
pdf(fname, width=5, height=5)
print(g)
dev.off()
crop_plot(fname)


fname = paste0(dir, "lhs_f1.jpg")
jpeg(fname, width=4*resol, height=4*resol, res=resol)
print(g)
dev.off()


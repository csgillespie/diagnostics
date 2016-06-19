#####################################################
## Schlogl latin hyper cube
#####################################################

x = readRDS(file="N_10000_schlogl_lhs.RData")
res = readRDS(file="N_10000_schlogl_residuals.RData")


dd = as.data.frame(x)
dd$res = res[,1]
library(ggplot2)
library(plotly)
theme_set(theme_bw())
alpha=230
palette(c(rgb(200,79,178, alpha=alpha,maxColorValue=255), 
          rgb(105,147,45, alpha=alpha, maxColorValue=255),
          rgb(85,130,169, alpha=alpha, maxColorValue=255),
          rgb(204,74,83, alpha=alpha, maxColorValue=255),
          rgb(183,110,39, alpha=alpha, maxColorValue=255),
          rgb(131,108,192, alpha=alpha, maxColorValue=255),
          rgb(63,142,96, alpha=alpha, maxColorValue=255)))




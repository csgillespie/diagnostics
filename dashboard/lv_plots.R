library("ggplot2")
#source("R/helper.R")
library(plotly)
alpha=230
palette(c(rgb(200,79,178, alpha=alpha,maxColorValue=255), 
          rgb(105,147,45, alpha=alpha, maxColorValue=255),
          rgb(85,130,169, alpha=alpha, maxColorValue=255),
          rgb(204,74,83, alpha=alpha, maxColorValue=255),
          rgb(183,110,39, alpha=alpha, maxColorValue=255),
          rgb(131,108,192, alpha=alpha, maxColorValue=255),
          rgb(63,142,96, alpha=alpha, maxColorValue=255)))


make_plot = function(maxtime) {
  if(maxtime ==10) {
    mc = readRDS(file="lv_mc_10.RData")
    lna = readRDS(file="lv_lna_10.RData")
    probs = readRDS(file="lv_extinction_t10.RData")
  }else if(maxtime == 30) {
    mc = readRDS(file="lv_mc_30.RData")
    lna = readRDS(file="lv_lna_30.RData")
    probs = readRDS(file="lv_extinction_t30.RData")
  } else if(maxtime == 100) {
    mc = readRDS(file="lv_mc_100.RData")  
    lna = readRDS(file="lv_lna_100.RData")
    probs = readRDS(file="lv_extinction_t100.RData")
  }
  x = readRDS(file="lv_lhs.RData")
  theme_set(theme_bw())
  
dd = data.frame(c1 = x[,1], c2 = x[,2], z=(mc[,3]- lna[,3])/lna[,6])
dd$Extreme = abs(dd$z) > 5 | is.na(dd$z)
dd$Extreme = mc[,2] < -1 | mc[,3] < -1| mc[,4] < -1| mc[,6] < -1
#dd$Extreme[lna[,2] < -1 | lna[,3] < -1| lna[,4] < -1| lna[,6] < -1] = 2
#dd$Extreme[!dd$Extreme] = (abs(dd$z[!dd$Extreme]) > 5)*3
dd$Extreme  = factor(dd$Extreme)
dd$prey_death = probs[,1]/1000
dd$pred_death = probs[,2]/1000
#dd = dd[-c(44, 79),]
table(dd$Extreme)

g1 = ggplot(dd, aes(c1, c2))  + 
  scale_x_log10(limits=c(1e-6, 1), 
                breaks=c(10^(-6:0)), 
                labels=c(expression(10^{-6}), "", expression(10^{-4}), "", expression(10^{-2}), "", expression(10^{0}))) + 
  scale_y_log10(limits=c(1e-6, 1), 
                breaks=c(10^(-6:0)), 
                labels=c(expression(10^{-6}), "", expression(10^{-4}), "", expression(10^{-2}), "", expression(10^{0}))) +  
  xlab("c1") + ylab("c2") + 
  guides(colour=FALSE, shape=FALSE, size=FALSE) + 
  # annotation_logticks() + 
  scale_color_manual(values=c(2, 4))
#  annotate("text", 10^{-6}, 10^0, label="(a)", size=4, hjust=-0.5) + 
#theme(legend.position = "top", legend.direction = "horizontal") + 
#scale_size(guide = guide_legend(title = "Extinction")) + 
#theme(legend.text = element_text(size = 8), 
#      legend.title = element_text(size = 9))
g1
}
#make_plot(30)

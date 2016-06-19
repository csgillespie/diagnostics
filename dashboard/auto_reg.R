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

######################
## Prior
######################
set.seed(1)
s = sample(1:10000, size = 2500)
dd_prior = read.csv("prior_analysis.csv", header=FALSE)[s,]
pars_prior = read.csv("prior.csv", header=FALSE)[s,]
x = seq(-4, 4, length.out = 1000)


######################
#Posterior
######################
dd_post = read.csv("post_analysis.csv", header=FALSE)[s,]
pars_post = read.csv("post.csv", header=FALSE)[s,]



#rsconnect::deployApp("diagnostics.Rmd", account="csgillespie")

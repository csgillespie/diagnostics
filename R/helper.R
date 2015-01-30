###################################################################
# Helper functions
###################################################################

## lhs outputs values on the (0, 1) scale
## Rescale
scale_lhs = function(x, c, a, b) {
  c + x*(b-a)
}

resol = 1000

#############################################
## Plotting functions
#############################################
## Sensible margins for base graphics
setnicepar = function(mar=c(3,3,2,1), 
                      mgp=c(2,0.4,0), tck=-.01,
                      cex.axis=0.9, las=1,...) {
  par(mar=mar, 
      mgp=mgp, tck=tck,
      cex.axis=cex.axis, las=las,...)
}
  


## Get pretty colours
mypalette = function(alpha=255) {
      palette(c(rgb(200,79,178, alpha=alpha,maxColorValue=255), 
              rgb(105,147,45, alpha=alpha, maxColorValue=255),
              rgb(85,130,169, alpha=alpha, maxColorValue=255),
              rgb(204,74,83, alpha=alpha, maxColorValue=255),
              rgb(183,110,39, alpha=alpha, maxColorValue=255),
              rgb(131,108,192, alpha=alpha, maxColorValue=255),
              rgb(63,142,96, alpha=alpha, maxColorValue=255)))
}

## Save a ggplot object
pdf_plot = function(g, fname, dir="", height=6, width=6){
  fname = paste0(dir, fname)
  pdf(fname, width=width, height=height)
  print(g)
  dev.off()
  crop_plot(fname)
}


## Use pdfcrop to crop the figures
crop_plot = function(fname) {
  cmd = paste("pdfcrop", fname)
  cmd = paste(cmd, "&& rm -v", fname)
  system(cmd)
}

vplayout = function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

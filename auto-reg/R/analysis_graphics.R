######################
## Prior
######################
dd_prior = read.csv("auto-reg/data/prior_analysis.csv", header=FALSE)
pars_prior = read.csv("auto-reg/data/prior.csv", header=FALSE)
sim = read.csv("auto-reg/data/sim.csv", header=FALSE)
x = seq(-4, 4, length.out = 1000)


######################
#Posterior
######################
dd_post = read.csv("auto-reg/data/post_analysis.csv", header=FALSE)
pars_post = read.csv("auto-reg/data/post.csv", header=FALSE)


## Graphics
## Pretty colours and resolutions
mypalette()
resol = 1000


## Prior graphic
jpeg("graphics/figure5.jpg",  width=12*resol, height = 8*resol, res = resol)
I = 6
setnicepar(mfrow=c(2, 3), cex.axis = 1, cex=1.3)
hist(dd_prior[,I], breaks="fd", freq=F, col="grey70", main=NULL, ylim=c(0, 0.5), xlim=c(-4.5, 4.5), 
     xlab=expression(paste(e[i], "*")));
lines(x, dnorm(x), col=3, lwd=3); text(-4, 0.48, "(a)")

qqnorm(dd_prior[,I], pch=21, bg=2, cex=0.7, xlim=c(-4.5, 4.5), ylim=c(-4.5, 4.5), main=NULL); 
grid(); abline(0, 1, col=3, lwd=3); text(-4, 4, "(b)")

plot(log(pars_prior[dd_prior[,7]+1,3]),dd_prior[,I],type="n", ylim=c(-4.5, 4.5), 
     ylab=expression(paste(e[i], "*")), xlab=expression(log(c[2]))); 
grid()
points(log(pars_prior[dd_prior[,7]+1,3]),dd_prior[,I], pch=21, bg=2, cex=0.5)
text(-4.8, 4.2, "(c)")

abline(h=c(qnorm(0.01/2), -qnorm(0.01/2)), col=3, lty=3, lwd=2)
abline(h=c(qnorm(0.001/2), -qnorm(0.001/2)), col=3, lty=3, lwd=2)
abline(h=c(qnorm(0.0001/2), -qnorm(0.0001/2)), col=3, lty=3, lwd=2)

I = 4
hist(dd_prior[,I], breaks="fd", freq=F, col="grey70", main=NULL, ylim=c(0, 0.8),
     xlab=expression(paste(e[i], "*")));lines(x, dnorm(x), col=3, lwd=3)
text(-11, 0.76, "(d)")
qqnorm(dd_prior[,I], pch=21, bg=2, cex=0.7, xlim=c(-5, 5), ylim=c(-15, 15), main=NULL); 
grid(); abline(0, 1, col=3, lwd=3); text(-4, 13, "(e)")

plot(log(pars_prior[dd_prior[,7]+1,3]),dd_prior[,I],type="n", ylim=c(-15, 15), 
     ylab=expression(paste(e[i], "*")), xlab=expression(log(c[2]))); 
grid(); text(-4.8, 13, "(f)")
points(log(pars_prior[dd_prior[,7]+1,3]),dd_prior[,I], pch=21, bg=2, cex=0.5)

abline(h=c(qnorm(0.01/2), -qnorm(0.01/2)), col=3, lty=3, lwd=2)
abline(h=c(qnorm(0.001/2), -qnorm(0.001/2)), col=3, lty=3, lwd=2)
abline(h=c(qnorm(0.0001/2), -qnorm(0.0001/2)), col=3, lty=3, lwd=2)

dev.off()


## Post graphic
jpeg("graphics/figure6.jpg",  width=12*resol, height = 8*resol, res = resol)
setnicepar(mfrow=c(2, 3), cex.axis = 1, cex=1.3)
I = 6
hist(dd_post[,I], breaks="fd", col="grey70", freq=FALSE, ylim=c(0, 0.5), main=NULL,
     xlim=c(-4.5, 4.5), 
     xlab=expression(paste(e[i], "*")));
lines(x, dnorm(x), col=3, lwd=3);text(-4, 0.48, "(a)")
qqnorm(dd_post[,I], pch=21, bg=2, cex=0.7, xlim=c(-4.5, 4.5), ylim=c(-4.5, 4.5), main=NULL);grid(); 
abline(0, 1, col=3, lwd=3); text(-4, 4, "(b)")

plot(log(pars_post[dd_post[,7]+1,3]),dd_post[,I],type="n", ylim=c(-4.5, 4.5), 
     ylab=expression(paste(e[i], "*")), xlab=expression(log(c[2])));  grid()
points(log(pars_post[dd_post[,7]+1,3]),dd_post[,I], pch=21, bg=2, cex=0.5)
abline(h=c(qnorm(0.01/2), -qnorm(0.01/2)), col=3, lty=3, lwd=2)
abline(h=c(qnorm(0.001/2), -qnorm(0.001/2)), col=3, lty=3, lwd=2)
abline(h=c(qnorm(0.0001/2), -qnorm(0.0001/2)), col=3, lty=3, lwd=2)
text(-1.7, 4.2, "(c)")


I=4
hist(dd_post[,I], breaks="fd", col="grey70", freq=FALSE, ylim=c(0, 0.5), main=NULL, 
     xlim=c(-4.5, 4.5), 
     xlab=expression(paste(e[i], "*")))
lines(x, dnorm(x), col=3, lwd=3);text(-4, 0.48, "(d)")
qqnorm(dd_post[,I], pch=21, bg=2, cex=0.7, xlim=c(-4.5, 4.5), ylim=c(-4.5, 4.5), main=NULL);grid(); 

abline(0, 1, col=3, lwd=3); text(-4, 4, "(e)")
plot(log(pars_post[dd_post[,7]+1,3]),dd_post[,I],type="n", ylim=c(-4.5, 4.5), 
     ylab=expression(paste(e[i], "*")), xlab=expression(log(c[2])));  grid()
points(log(pars_post[dd_post[,7]+1,3]),dd_post[,I], pch=21, bg=2, cex=0.5)
abline(h=c(qnorm(0.01/2), -qnorm(0.01/2)), col=3, lty=3, lwd=2)
abline(h=c(qnorm(0.001/2), -qnorm(0.001/2)), col=3, lty=3, lwd=2)
abline(h=c(qnorm(0.0001/2), -qnorm(0.0001/2)), col=3, lty=3, lwd=2)
text(-1.7, 4.2, "(f)")
dev.off()


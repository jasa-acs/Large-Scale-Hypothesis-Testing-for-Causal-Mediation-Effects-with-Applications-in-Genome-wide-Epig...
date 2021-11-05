###############
## This code is for illustrating normal approximation of standarized product method is not accurate
## and under what kind of condition it will be accurate, of course depends on significance threshold.
## 1. we compare probability density functions using kernel density estimator
## 2. We can Normal QQ-plot to assess normalaity of product testing statistics (Sobel 1982)
####################################
library(mvtnorm)
n.sim = 10000

pdf("Fig2_prod_normal_approx.pdf",width=20,height=10)
#par(mfcol=c(2,3),by)
#png("prod_normal_approx.png",width=800,height=600)
layout(matrix(c(1,2,3,4,5,6),2, 3, byrow = T))

#### Case 1, gamma=0, beta = 0.2, mu_beta = 0.1*sqrt(n)
## sample size 100, R-sq = 0.75
n= 100
mu_beta = 0.1*sqrt(n)
a = rmvnorm(n.sim,c(0,mu_beta),diag(2))
## Sobel test statistic
t = a[,1]*a[,2]/sqrt(a[,1]^2 + a[,2]^2)
plot(density(t),main="Case 1, n = 100",cex.main=3,col="black",xlab=" ", ylab= " ", xlim=c(-4,4),ylim=c(0,0.8),cex=3,lwd=3)
x.norm = seq(-4,4,length = 1000)
lines(x.norm,dnorm(x.norm),lty=2,lwd=3)
#legend("topleft",legend=c(expression(mu[beta]==6),"N(0,1)"),col=c("black","green"),lty=1,cex = 2.3,bty="n")
legend("topleft",legend=c("Sobel's test","N(0,1)"),lty=c(1,2),cex = 2.3,bty="n")

n= 500
mu_beta = 0.1*sqrt(n)
a = rmvnorm(n.sim,c(0,mu_beta),diag(2))
## Sobel test statistic
t = a[,1]*a[,2]/sqrt(a[,1]^2 + a[,2]^2)
plot(density(t),main="Case 1, n = 500",cex.main=3,col="black",xlab=" ", ylab= " ", xlim=c(-4,4),ylim=c(0,0.8),cex=3,lwd=3)
x.norm = seq(-4,4,length = 1000)
lines(x.norm,dnorm(x.norm),lty=2,lwd=3)
#legend("topleft",legend=c(expression(mu[beta]==6),"N(0,1)"),col=c("black","green"),lty=1,cex = 2.3,bty="n")
legend("topleft",legend=c("Sobel's test","N(0,1)"),lty=c(1,2),cex = 2.3,bty="n")

n= 5000
mu_beta = 0.1*sqrt(n)
a = rmvnorm(n.sim,c(0,mu_beta),diag(2))
## Sobel test statistic
t = a[,1]*a[,2]/sqrt(a[,1]^2 + a[,2]^2)
plot(density(t),main="Case 1, n = 5000",cex.main=3,col="black",xlab=" ", ylab=" ", xlim=c(-4,4),ylim=c(0,0.8),cex=3,lwd=3)
x.norm = seq(-4,4,length = 1000)
lines(x.norm,dnorm(x.norm),lty=2,lwd=3)
#legend("topleft",legend=c(expression(mu[beta]==6),"N(0,1)"),col=c("black","green"),lty=1,cex = 2.3,bty="n")
legend("topleft",legend=c("Sobel's test","N(0,1)"),lty=c(1,2),cex = 2.3,bty="n")



### Case 3, gamma=0, beta =0
n= 100
mu_beta = 0*sqrt(n)
a = rmvnorm(n.sim,c(0,mu_beta),diag(2))
## Sobel test statistic
t = a[,1]*a[,2]/sqrt(a[,1]^2 + a[,2]^2)
plot(density(t),main="Case 3, n = 100",cex.main=3,col="black",xlab=" ", ylab= " ", xlim=c(-4,4),ylim=c(0,0.8),cex=3,lwd=3)
x.norm = seq(-4,4,length = 1000)
lines(x.norm,dnorm(x.norm),lty=2,lwd=3)
#legend("topleft",legend=c(expression(mu[beta]==6),"N(0,1)"),col=c("black","green"),lty=1,cex = 2.3,bty="n")
legend("topleft",legend=c("Sobel's test","N(0,1)"),lty=c(1,2),cex = 2.3,bty="n")

n= 500
mu_beta = 0*sqrt(n)
a = rmvnorm(n.sim,c(0,mu_beta),diag(2))
## Sobel test statistic
t = a[,1]*a[,2]/sqrt(a[,1]^2 + a[,2]^2)
plot(density(t),main="Case 3, n = 500",cex.main=3,col="black",xlab=" ", ylab= " ", xlim=c(-4,4),ylim=c(0,0.8),cex=3,lwd=3)
x.norm = seq(-4,4,length = 1000)
lines(x.norm,dnorm(x.norm),lty=2,lwd=3)
#legend("topleft",legend=c(expression(mu[beta]==6),"N(0,1)"),col=c("black","green"),lty=1,cex = 2.3,bty="n")
legend("topleft",legend=c("Sobel's test","N(0,1)"),lty=c(1,2),cex = 2.3,bty="n")


n= 5000
mu_beta = 0*sqrt(n)
a = rmvnorm(n.sim,c(0,mu_beta),diag(2))
## Sobel test statistic
t = a[,1]*a[,2]/sqrt(a[,1]^2 + a[,2]^2)
plot(density(t),main="Case 3, n = 5000",cex.main=3,col="black",xlab=" ", ylab=" ", xlim=c(-4,4),ylim=c(0,0.8),cex=3,lwd=3)
x.norm = seq(-4,4,length = 1000)
lines(x.norm,dnorm(x.norm),lty=2,lwd=3)
#legend("topleft",legend=c(expression(mu[beta]==6),"N(0,1)"),col=c("black","green"),lty=1,cex = 2.3,bty="n")
legend("topleft",legend=c("Sobel's test","N(0,1)"),lty=c(1,2),cex = 2.3,bty="n")



dev.off()


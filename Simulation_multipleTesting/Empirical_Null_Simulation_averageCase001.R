

library(qvalue)
library(qqman)
library(locfdr)

###########
myQQ = function (pmatrix, main="title",...) {
  pvector1 = pmatrix[,1]
  pvector2 = pmatrix[,2]

  o1 = -log10(sort(pvector1, decreasing = FALSE))
  e1= -log10(ppoints(length(pvector1)))

  o2 = -log10(sort(pvector2, decreasing = FALSE))
  e2= -log10(ppoints(length(pvector2)))

  plot(e1,o1,col="black",main=main,cex.main=2, pch=20,xlim = c(0, max(e1)+1), cex.lab =1.5, ylim = c(0,
                                                                                                     max(e1)+1), xlab = expression(Expected ~ ~-log[10](italic(p))),
       ylab = expression(Observed ~ ~-log[10](italic(p))))
  points(e2,o2,col="blue",pch=20)
  abline(0, 1, col = "red")
  legend("topleft",cex=1.5,legend = c("Uncorrected DAST","Corrected DAST "),pch=c(20,20),col = c("black","blue"))
}

## setup parameters
n <- 3e5  ## number of tests
w1 <- 0.01 ## Case 1
w2 <- 0.01 ## Case 2
w3 <- 1 - w1 - w2


# if generate all data under the null
n1 <- n*w1
n2 <- n*w2
n3 <- n - n1 - n2

## data generation
Za <- rnorm(n)
Zb <- rnorm(n)
## add the means
## add one signal
Za.ob <- Za + c(rep(0, n1), rnorm(n2,2,1), rep(0, n3))
Zb.ob <- Zb + c(rnorm(n1,2,1), rep(0, n2),rep(0, n3))

pa.ob <- pnorm(abs(Za.ob), lower.tail = F) * 2
pb.ob <- pnorm(abs(Zb.ob), lower.tail = F) * 2




u3.ob <- pa.ob*I(pa.ob > pb.ob) + pb.ob*I(pa.ob <= pb.ob)
p3.ob <- u3.ob^2
min(p3.ob)
#########


## estimate proportions
qa <- qvalue(pa.ob)
pa.est <- qa$pi0
qb <- qvalue(pb.ob)
pb.est <- qb$pi0

w1.est <- 1 - pb.est
w2.est <- 1 - pa.est
w3.est <- 1- w1.est - w2.est

## original DAT
pDAT <- pa.ob * w1 + pb.ob * w2 + p3.ob * w3
#qqman::qq(pDAT,main = "not corrected")
pvalue <- pDAT
chisq <- qchisq(1-pvalue,1)
(lambda = median(chisq)/qchisq(0.5,1))
######
####### empirical null
z = qnorm(1-pDAT)
max(z)
mean(z)
sd(z)
res <- locfdr(z,nulltype = 1)
res$fp0
mean.emp = res$fp0["mlest","delta"]
sd.emp = res$fp0["mlest","sigma"]
sd.emp
#######################
## empirial p-value
pval.emp = pnorm(z,mean = mean.emp,sd = sd.emp ,lower.tail = F)
#pval.emp = pnorm(z,mean = 0,sd = sd.emp ,lower.tail = F)
qqman::qq(pval.emp,main="empirical")
par(mar=c(4.1,6.1,3.1,2.1))
res_averageCase001 = cbind(pDAT,pval.emp)
save(res_averageCase001,file = "res_averageCase001.Rdata")

myQQ(cbind(pDAT,pval.emp),main = "Average-case Scenario")

pvalue <- pval.emp
chisq <- qchisq(1-pvalue,1)
quantile(chisq,0.75)/qchisq(0.75,1)
(lambda = median(chisq)/qchisq(0.5,1))


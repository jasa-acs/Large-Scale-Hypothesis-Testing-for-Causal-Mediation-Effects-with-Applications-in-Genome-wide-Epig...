#######################
## NAS data analysis 
#####################
library(qqman)
library(locfdr)
library(qvalue)
#########
source("./EstNull.func.R")
source("./epsest.func.R")
source("./myManhattan.R")
source("./myQQNAS.R")

Sobel = function(Z1,Z2){
  T = (Z1*Z2)^2/(Z1^2+Z2^2)
  pchisq(T,df=1,lower.tail = F)
}
##################################
med = read.csv("mval_mediator_sorted.csv",header=T)
outcome = read.csv("mval_outcome_sorted.csv",header=T)
dat = merge(med,outcome, by="probe") ## crucial step
n_cpg = dim(dat)[1]
Z_a = dat$z.x
Z_b = dat$z.y
p_a = dat$Pval.x
p_b = dat$Pval.y
p.mat = cbind(p_a,p_b)
p3 = (apply(p.mat,1,max))^2
dat$p3 = p3
##### Sobel Test p-value
P_sobel = Sobel(Z_a,Z_b)

## estimate proportions
mu.sigma_a = EstNull.func(Z_a)
mu.sigma_b = EstNull.func(Z_b)
pi0a = 1- epsest.func(Z_a,mu.sigma_a$mu,mu.sigma_a$s)
pi0b = 1- epsest.func(Z_b,mu.sigma_b$mu,mu.sigma_b$s)
pi0a
pi0b
#### calculate the weights under the null
wg1 = pi0a*(1-pi0b)
wg2 = (1-pi0a)*pi0b
wg3 = pi0a*pi0b
wg4 = (1-pi0a)*(1-pi0b)
c(wg1,wg2,wg3,wg4)
round(c(wg1,wg2,wg3,wg4),5)

### normalize under the null
wg.sum = wg1 + wg2 + wg3
wg.std = c(wg1,wg2,wg3)/wg.sum
wg.std
p.average = wg.std[1]*p_a + wg.std[2]*p_b + wg.std[3]*p3

dat$nie = dat$coef.x * dat$coef.y
## exact se
dat$nie_se = sqrt(dat$coef.x^2*dat$se_coef.y^2 + dat$coef.y^2*dat$se_coef.x^2 + dat$se_coef.x^2*dat$se_coef.y^2)

dat$p.average = p.average
dat$P_sobel = P_sobel
dat$MaxP = apply(p.mat,1,max)
P_sobel = dat$P_sobel
MaxP = dat$MaxP
p.average = dat$p.average
########################
## Empirical null estimation
########################
Z_DACT = qnorm(p.average,lower.tail = F)
### Efron R package locfdr
efron_res = locfdr(Z_DACT,plot = 2)
efron_res$fp0[3,]
efron.mean = efron_res$fp0[3,1]
efron.sd = efron_res$fp0[3,2]
fdr.cutoff = 0.2
sum(efron_res$fdr < fdr.cutoff)
efron_res$mat[which(efron_res$mat[,"fdr"]<fdr.cutoff)[1],"Fdrright"] ## tail FDR

Efron_correct = pnorm(Z_DACT,lower.tail = F, mean=efron.mean, sd=efron.sd)
#### Jin and Cai
mu.sigma.nas = EstNull.func(Z_DACT)
mu.sigma.nas
pi0_nas = 1- epsest.func(Z_DACT,mu.sigma.nas$mu,mu.sigma.nas$s)
pi0_nas
DACT_correct = pnorm(Z_DACT,lower.tail = F, mean=mu.sigma.nas$mu, sd=mu.sigma.nas$s)

### qvalue approach
localFDR <- lfdr(p = DACT_correct,pi0 = 1)
sum(localFDR < 0.2) ## identical to Efron because it implements efron local
qobj = qvalue(DACT_correct)
sum(qobj$qvalues < 0.03846856 )
### Bonferroni equivalent FDR
sum(qobj$qvalues < 0.005 )
dat$probe[qobj$qvalues < 0.005]
############### plot 
q.grid = seq(0,0.3,length.out = 10)
sig.grid = c()
for(i in 1:length(q.grid)){
  sig.grid[i] =sum(qobj$qvalues < q.grid[i])
}
png("qvalue-cutoff.png")
plot(q.grid,sig.grid,,type="l",xlab = "q-value cut-off",ylab = "number of significant tests")
abline(h=29,col="red")
dev.off()
#################### correction


### output to dat dataframe
length(efron_res$fdr)
dat$Efron_correct = Efron_correct
dat$DACT_correct = DACT_correct
dat$qvalues = qobj$qvalues
dat$lfdr = efron_res$fdr
write.csv(dat,file = "mval_mediation_JC2020.csv",quote = F)

plot(-log10(Efron_correct),-log10(DACT_correct))

png("Efron-JC.png",width = 480,height = 480)
plot(-log10(Efron_correct),-log10(DACT_correct), cex.main=2, main ="Efron's versus JC's Estimation" ,xlab = "Efron's Corrected DACT on the -log10 scale", ylab = "JC's Corrected DACT on the -log10 scale",cex.lab=1)
abline(a=0,b=1,col="red")
dev.off()

### lambda, genomic inflation factor
median(qchisq(1-DACT_correct,df=1))/qchisq(0.5,df=1)
median(qchisq(1-MaxP,df=1))/qchisq(0.5,df=1)
median(qchisq(1-P_sobel,df=1))/qchisq(0.5,df=1)
##################
## QQ plot before correction
##################
library(qqman)
png("QQ_3_methods_JC2020.png",width = 480*3,height = 480)
par(mfrow=c(1,3))
par(mar=c(5.1,5.1,4.1,2.1))
qq(P_sobel,ylim=c(0,17),cex.axis=2,cex.lab=2.5,cex=3,lwd=3)
abline(a=0,b=1,col="red",lwd=2)
title("Sobel",cex.main=3)
qq(MaxP,ylim=c(0,17),cex.axis=2,cex.lab=2.5,cex=3)
abline(a=0,b=1,col="red",lwd=2)
title("MaxP",cex.main=3)
qq(p.average,cex.lab=2,ylim=c(0,17),cex.axis=2,cex.lab=2.5,cex=3)
abline(a=0,b=1,col="red",lwd=2)
title("DACT",cex.main=3)
dev.off()
##################
## QQ plot after correction
##################
png("QQ_3_methods_JC2020_corrected.png",width = 480*3,height = 480)
par(mfrow=c(1,3))
par(mar=c(5.1,5.1,4.1,2.1))
qq(P_sobel,ylim=c(0,17),cex.axis=2,cex.lab=2.5,cex=3,lwd=3)
abline(a=0,b=1,col="red",lwd=2)
title("Sobel",cex.main=3)
qq(MaxP,ylim=c(0,17),cex.axis=2,cex.lab=2.5,cex=3)
abline(a=0,b=1,col="red",lwd=2)
title("MaxP",cex.main=3)
qq(DACT_correct,cex.lab=2,ylim=c(0,17),cex.axis=2,cex.lab=2.5,cex=3)
abline(a=0,b=1,col="red",lwd=2)
title("DACT",cex.main=3)
dev.off()


#####################################
### Histogram of Sobel, MaxP, DACT before transform
#####################################
png("DACT_MaxP_Sobel_hist.png",width = 480*3,height = 480)
par(mfrow=c(1,3))
par(mar=c(5.1,5.1,4.1,2.1))
hist(P_sobel,main = "Sobel",cex.main =3,xlab = "Sobel",cex.lab=1.5,breaks = 30)
hist(MaxP,main = "MaxP",cex.main =3,xlab = "MaxP",cex.lab=1.5,breaks = 30)
hist(p.average,main = "DACT",cex.main =3,xlab = "DACT",cex.lab=1.5,breaks = 30)
dev.off()


png("DACT_hist_compare.png",width = 480*2,height = 360)
par(mfrow=c(1,2))
par(mar=c(5.1,5.1,4.1,2.1))
hist(p.average,main = "DACT",cex.main =2,xlab = "DACT",cex.lab=1.5,breaks = 30)
hist(DACT_correct,main = "Corrected DACT",cex.main =2,xlab = "Corrected DACT",cex.lab=1.5,breaks = 30)
dev.off()
#####################################
### Histogram of Sobel, MaxP, DACT after transform
#####################################
png("DACT_MaxP_Sobel_hist_transform.png",width = 480*3,height = 480)
par(mfrow=c(1,3))
par(mar=c(5.1,5.1,4.1,2.1))
hist(qnorm(P_sobel,lower.tail = F),main = "Transformed Sobel",cex.main =3,xlab = "Transformed Sobel",cex.lab=1.5,breaks = 30)
hist(qnorm(MaxP,lower.tail = F),main = "Transformed MaxP",cex.main =3,xlab = "Transformed MaxP",cex.lab=1.5,breaks = 30)
hist(qnorm(p.average,lower.tail = F),main = "Transformed DACT",cex.main =3,xlab = "Transformed DACT",cex.lab=1.5,breaks = 30)
dev.off()


#####################################
### Histogram of Z_DACT with estimated null curve
#####################################
png("DACT-transform-hist.png",width = 480*1.2,height = 480)
#par(mfrow=c(1,2))
#par(mar=c(5.1,5.1,4.1,2.1))
#hist(qnorm(P_sobel,lower.tail = F),main = "Transformed Sobel",cex.main =3,xlab = "Transformed Sobel",cex.lab=1.5,breaks = 30)
#hist(qnorm(MaxP,lower.tail = F),main = "Transformed MaxP",cex.main =3,xlab = "Transformed MaxP",cex.lab=1.5,breaks = 30)
xlab = paste("delta:", round(mu.sigma.nas$mu, 3), "sigma:", round(mu.sigma.nas$s, 3), "p0:", round(pi0_nas, 3))
hist(Z_DACT,main = "Transformed DACT",cex.main =2,cex.lab=1.5,breaks = 120,xlab=xlab,probability = TRUE)
#lines(density(Z_DACT), col = "green",lwd=3)   
curve(dnorm(x, mean=mu.sigma.nas$mu, sd=mu.sigma.nas$s),col="green",lwd=2.5,add=TRUE)
#qqnorm(Z_DACT,main = "Transformed DACT",cex.main =2 ,cex.lab=1.5)
#qqline(Z_DACT, col = "red") 
#boxplot(Z_DACT,main="Transformed DACT",cex.main =2)
dev.off()

###############################
## DACT verus corrected DACT
##############################
# png("DACT-correct-compare.png",width = 480,height = 480)
# plot(-log10(p.average),-log10(DACT_correct),xlab = "DACT on the -log10 scale", ylab = "Corrected DACT on the -log10 scale",cex.lab=1)
# dev.off()

#####################
## Manhattan plot
#####################



head(dat)
results = dat[,c("probe","nie","DACT_correct")]
names(results) = c("probe","coef","pvalue")
annot_file = read.table("cpg_chr_pos.txt",header=T)
names(results)[1]="SiteID"
dat1 = merge(annot_file,results,by="SiteID")
head(dat)
names(dat)[1]="SiteID"
outall = merge(annot_file,dat,by="SiteID")
head(outall)
write.csv(outall,file="outall.csv",quote = F,row.names = F)
##############
png("QQ_gamma_beta.png",width = 480*2,height = 480)
par(mfrow=c(1,2))
par(mar=c(5.1,5.1,4.1,2.1))
qq(outall$Pval.x,ylim=c(0,17),cex.axis=2,cex.lab=2.5,cex=3,lwd=3,xlab="p-values for A-M")
abline(a=0,b=1,col="red",lwd=2)
title("p-values for A-M",cex.main=3)
qq(outall$Pval.y,ylim=c(0,17),cex.axis=2,cex.lab=2.5,cex=3,xlab="p-values for M-Y")
abline(a=0,b=1,col="red",lwd=2)
title("p-values for M-Y",cex.main=3)
dev.off()
########
hist(outall$Pval.x)
hist(outall$Pval.y)

#############
head(dat1)
dat1$Chr[dat1$Chr== "X"]= "23"
dat1$Chr[dat1$Chr== "Y"]= "24"
res.for.manhattan = dat1[,c("SiteID","Chr","Loc","pvalue")]
names(res.for.manhattan) = c("SNP","CHR","BP","P")
#head(res.for.manhattan)
chr  = res.for.manhattan$CHR
n.cpg = dim(res.for.manhattan)[1]
res.for.manhattan[,"CHR"] = as.numeric(chr)
genome.sig.level = 0.05/n.cpg
head(res.for.manhattan)
ylim = max(-log10(res.for.manhattan$P))+1

########
png("./Manhattan.png",height=480,width=800)
myManhattan(res.for.manhattan,col = c("blue4", "orange3"), ylim=c(0,ylim), suggestiveline = F, genomewideline = F,annotatePval =genome.sig.level )
abline(h=-log10(genome.sig.level),col="red")
dev.off()
#####################
## Volcano plot
#####################
xmax = max(abs(dat1$coef))
png("./Volcano.png")
plot(dat1$coef,-log10(dat1$pvalue),xlim=c(-xmax,xmax),xlab="Indirect Effect",ylab="-log10(P-value)",main = "Volcano Plot")
abline(h=-log10(genome.sig.level),col="red")
dev.off()
#########
### three in one 
round(mu.sigma.nas$mu, 3)
round(mu.sigma.nas$s, 3)
round(pi0_nas, 3)
png("Hist_QQ3in1_Volcano_corrected.png",width = 480*3,height = 480)
par(mfrow=c(1,3))
par(mar=c(6.5,5.5,4.1,2.1))
xlab = expression(paste(delta, ": -0.053, ", sigma,": 0.998, ",  p[0],": 1"))
hist(Z_DACT,main = "Transformed DACT",cex.main = 3,cex.lab=2.3,breaks = 120,xlab=xlab,probability = TRUE)
curve(dnorm(x, mean=mu.sigma.nas$mu, sd=mu.sigma.nas$s),col="green",lwd=2.5,add=TRUE)
## sobel, maxp, DACT
pmatrix = cbind(P_sobel,MaxP,DACT_correct)
myQQNAS(pmatrix = pmatrix,main = "QQ Plot")
plot(dat1$coef,-log10(dat1$pvalue),xlim=c(-xmax,xmax),xlab="Mediation Effect Size",ylab = expression(Observed ~ ~-log[10](italic(p))),main = "Volcano Plot",cex.main=3,cex.lab=2.3)
abline(h=-log10(genome.sig.level),col="red")
dev.off()



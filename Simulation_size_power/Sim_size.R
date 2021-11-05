library(cp4p)
library(limma)
library(qvalue)
library(qqman)
library(mvtnorm)



######################
## Sobel function 
######################
Sobel = function(Z1,Z2){
  T = (Z1*Z2)^2/(Z1^2+Z2^2)
  pchisq(T,df=1,lower.tail = F)
  
}
################
## DACT
################
DACT = function(p1,p2){
  
  #propTrueNull(p1) 
  pi0a = qvalue(p1)$pi0 ## qvalue more conservative, tend to give close to 1 estimate of null proportion
  #propTrueNull(p2) 
  pi0b = qvalue(p2)$pi0 ## qvalue more conservative, tend to give close to 1 estimate of null proportion
  
  p3 = (apply(cbind(p1,p2),1,max))^2
  #### calculate the weights 
  wg1 = pi0a*(1-pi0b)
  wg2 = (1-pi0a)*pi0b
  wg3 = pi0a*pi0b
  wg.sum = wg1 + wg2 + wg3
  wg.std = c(wg1,wg2,wg3)/wg.sum
  p.average = wg.std[1]*p1 + wg.std[2]*p2 + wg.std[3]*p3
  p.average
}

##############
## MT-COMP
#############
int<-10
B<-10000; pdf<-NULL
for (i in 1:B){
  pdfi<-besselK(x=int*i/B, nu=0)
  pdf<-c(pdf, pdfi)
  print(i); flush.console()
}
print(int)

myp<-function(cut){
  select<-(int*1:B/B)>cut
  pdf.sub<-pdf[select]
  pval<-sum(pdf.sub)/sum(pdf)
  return(pval)
}

MT_Comp<-function(a, b){
  ab<-a*b
  pp0<-sapply(abs(ab)/sqrt(1), myp)
  pp1<-sapply(abs(ab)/sd(a), myp)
  pp2<-sapply(abs(ab)/sd(b), myp)
  pp.comp<-pp1+pp2-pp0
  return(pp.comp)
}

###############################
### type I error rate function 
###############################
typeIerror = function(pvalues, cutoff = c(0.05,0.01,1e-3,1e-4)){
  len = length(cutoff)
  err.rate = rep(NA,len)
  for(i in 1:len){
    err.rate[i] = mean(pvalues < cutoff[i])
  }
  err.rate
}

#######################
##  Simulation Set up 
#######################





#######################
n = 1000 # sample size 
sim.num = 5e5
a = 0.2 ## effect a-- exposore on mediator
b = 0 ## effect b-- mediator on outcome 


#################
run_simulation = function(n=1000,sim.num=5e5,a=0.2,b=0){

sigma.M = 1
sigma.Y = 2
############
Z.M = c()
Z.Y = c()
p.M = c()
p.Y = c()
## exposure
A = rbinom(n,size =1,prob = 0.5 )
## covariates 
X1 = rnorm(n,mean = 10,sd = 1)
X2 = rnorm(n,mean = 5, sd = 1)

#############
ptm = proc.time()
for (i in 1:sim.num){
## Mediator 
M = a*A + 0.2*X1 + 0.3*X2 + rnorm(n,mean = 0,sd = sigma.M )
fit.M = lm(M~ A + X1 + X2)
Z.M[i] = summary(fit.M)$coef[2,3]
p.M[i] = summary(fit.M)$coef[2,4]
## outcome 
Y = b*M + A + 0.1*X1 + 0.2*X2 + rnorm(n,mean = 0, sd = sigma.Y)
fit.Y = lm(Y ~ M + A + X1 + X2)
Z.Y[i] = summary(fit.Y)$coef[2,3]
p.Y[i] = summary(fit.Y)$coef[2,4]
}
###################
### running time 
##################
proc.time() -ptm

## Sobel's Test
p.sobel = Sobel(Z.M,Z.Y)
## maxP 
p.MaxP = apply(cbind(p.M,p.Y),1,max)
### MT-COMP
p.MTComp = MT_Comp(Z.M,Z.Y)
## DACT 
p.DACT = DACT(p1=p.M,p2 = p.Y)

####################
## results summary 
###################
# png("./Sim_QQ_3_methods.png",width = 480*3,height = 480)
# par(mfrow=c(1,3))
# par(mar=c(5.1,5.1,4.1,2.1))
# qq(p.sobel,cex.lab=2,lwd=1.5,xlim=c(0,6),ylim=c(0,6),cex.axis=2)
# title("Sobel",cex.main=2)
# qq(p.MaxP,cex.lab=2,lwd=1.5,xlim=c(0,6),ylim=c(0,6),cex.axis=2)
# title("MaxP",cex.main=2)
# qq(p.DACT,cex.lab=2,lwd=1.5,xlim=c(0,6),ylim=c(0,6),cex.axis=2)
# title("DACT",cex.main=2)
# dev.off()

###############################
### type I error rates
###############################
out.sobel = typeIerror(p.sobel)
out.maxp = typeIerror(p.MaxP)
out.mtcomp = typeIerror(p.MTComp)
out.DACT = typeIerror(p.DACT)
out = c(out.sobel,out.maxp,out.mtcomp,out.DACT)
names(out)=c(rep("Sobel",4),rep("MaxP",4),rep("MTComp",4),rep("DACT",4))
out
}


###########
sim.num = 1e5

out1 = run_simulation(n=500,sim.num = sim.num, a=0.2,b=0)
out2 = run_simulation(n=500,sim.num = sim.num, a=0,b=0.2)
out3 = run_simulation(n=500,sim.num = sim.num, a=0,b=0)

out4 = run_simulation(n=1000,sim.num = sim.num, a=0.2,b=0)
out5 = run_simulation(n=1000,sim.num = sim.num, a=0,b=0.2)
out6 = run_simulation(n=1000,sim.num = sim.num, a=0,b=0)

out7 = run_simulation(n=2001,sim.num = sim.num, a=0.2,b=0)
out8 = run_simulation(n=2000,sim.num = sim.num, a=0,b=0.2)
out9 = run_simulation(n=2000,sim.num = sim.num, a=0,b=0)


out = rbind(out1,out2,out3,out4,out5,out6,out7,out8,out9)
write.csv(out,file = "size_output.csv",quote = F,row.names = F)

### power 
power1 = run_simulation(n=2000,sim.num = 3000, a=0.5,b=0.6)
power1

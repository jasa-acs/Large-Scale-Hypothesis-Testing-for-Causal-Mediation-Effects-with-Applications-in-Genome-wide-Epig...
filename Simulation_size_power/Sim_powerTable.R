

#######################
library(cp4p)
library(limma)
library(qvalue)
library(qqman)
library(mvtnorm)

Sobel = function(Z1,Z2){
  T = (Z1*Z2)^2/(Z1^2+Z2^2)
  pchisq(T,df=1,lower.tail = F)
  
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

DACT = function(p1,p2){ 
  pi0a = propTrueNull(p1) 
  #pi0a = qvalue(p1)$pi0 ## qvalue more conservative, tend to give close to 1 estimate of null proportion
  #print(pi0a)
  pi0b = propTrueNull(p2)
  #pi0b = qvalue(p2)$pi0 ## qvalue more conservative, tend to give close to 1 estimate of null proportion
  #print(pi0b)
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

testSize = function(p.sobel = p.sobel, p.MaxP = p.MaxP,p.DACT = p.DACT){
  
  sobel1 = mean(p.sobel < 0.05)
  sobel2 = mean(p.sobel < 0.01)
  sobel3 = mean(p.sobel < 0.001)
  
  maxp1 = mean(p.MaxP < 0.05)
  maxp2 = mean(p.MaxP < 0.01)
  maxp3 = mean(p.MaxP < 0.001)
  
  pDACT1 = mean(p.DACT < 0.05)
  pDACT2 = mean(p.DACT < 0.01)
  pDACT3 = mean(p.DACT < 0.001)
  out = list(Sobel = c(sobel1,sobel2,sobel3),MaxP =c(maxp1,maxp2,maxp3),pDACT = c(pDACT1,pDACT2,pDACT3))
  out
}

########################


mediationSim = function(a=0,b=0,sigma.M=1,sigma.Y=1,n=1000,sim.num = 1e4){
  
  ## containers
  Z.M = c()
  Z.Y = c()
  p.M = c()
  p.Y = c()
  ## exposure
  A = rbinom(n,size =1,prob = 0.5 )
  #sigma.A = sd(A)
  ## covariates 
  X1 = rnorm(n,mean = 10,sd = 1)
  X2 = rnorm(n,mean = 5, sd = 1)
  
  #############
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
  out = list(Z.M = Z.M, Z.Y = Z.Y,p.M =p.M,p.Y=p.Y)
  return(out)
}


##############
##########
### power 
###################
testPower = function(p.sobel = p.sobel, p.MaxP = p.MaxP,p.MTComp,p.DACT = p.DACT){
  
  sobel1 = mean(p.sobel < 0.05)
  maxp1 = mean(p.MaxP < 0.05)
  MTComp1 = mean(p.MTComp < 0.05)
  pDACT1 = mean(p.DACT < 0.05)
  
  out = c(sobel1,maxp1,MTComp1,pDACT1)
  out
}
### sample size grid 

powerSim = function(a,b,N,sim.num=1e4){
  
  res.power = matrix(NA,ncol=4,nrow=length(N))
  for (i in 1:length(N1)){
    ### power 
    res = mediationSim(a=a,b=b,sigma.M=1,sigma.Y=2,n=N[i],sim.num=sim.num)
    
    p.sobel = Sobel(res$Z.M,res$Z.Y)
    p.MaxP = apply(cbind(res$p.M,res$p.Y),1,max)
    p.MTComp = MT_Comp(res$Z.M,res$Z.Y)
    p.DACT = DACT(p1=res$p.M,p2 = res$p.Y)
    tmp = testPower(p.sobel,p.MaxP,p.MTComp,p.DACT)
    
    res.power[i,] = tmp
  }
  res.power = cbind(N,res.power)
  colnames(res.power) = c("N","Sobel","MaxP","MTCOMP","DACT")
  res.power 
}
####################


N1 = c(800,1000,1200)
a1 = 0.3
b1 = 0.133

N2 =  c(800,1000,1200)
a2 = 0.2
b2 = 0.2

N3 =  c(800,1000,1200)
a3 = 0.133
b3 = -0.3

N4 =  c(800,1000,1200)
a4 = 0.1
b4 = 0.4
#####################
power1 = powerSim(a=a1,b=b1,N=N1,sim.num=1e4)
power2 = powerSim(a=a2,b=b2,N=N2,sim.num=1e4)
power3 = powerSim(a=a3,b=b3,N=N3,sim.num=1e4)
power4 = powerSim(a=a4,b=b4,N=N4,sim.num=1e4)

power4

t(power1 )
power2 
power3 

power.out = rbind(power1, power2, power3,power4)
power.out
t(power.out)

write.csv(t(power.out),file="powerout.csv",quote=F)
# #### function to plot power results 
plotPower = function(power_mat,a=0.1,b=0.1){
  pnum = dim(power_mat)[2]
  matplot(power_mat[,1],power_mat[,2:pnum],type="b",pch=1:(pnum-1),col = 1:(pnum-1), ylim=c(0,1),xlab="N",ylab="Power")
  axis(side = 1, at = power_mat[,1])
  title(bquote(paste(alpha==.(a),", ",beta==.(b))))
  legend("topleft", legend = colnames(power_mat)[2:pnum], pch=1:(pnum-1),col=1:(pnum-1),cex=1,bty='n') 
}




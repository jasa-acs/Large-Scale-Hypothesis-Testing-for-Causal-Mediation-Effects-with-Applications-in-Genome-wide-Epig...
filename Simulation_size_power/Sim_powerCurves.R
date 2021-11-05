

#######################
#library(cp4p)
library(limma)
#library(qvalue)
#library(qqman)
#library(mvtnorm)

Sobel = function(Z1,Z2){
  T = (Z1*Z2)^2/(Z1^2+Z2^2)
  pchisq(T,df=1,lower.tail = F)
  
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



###################
### power 
###################
testPower = function(p.sobel = p.sobel, p.MaxP = p.MaxP,p.DACT = p.DACT){
  
  sobel1 = mean(p.sobel < 0.05)
  maxp1 = mean(p.MaxP < 0.05)
 # MTComp1 = mean(p.MTComp < 0.05)
  pDACT1 = mean(p.DACT < 0.05)
  
  out = c(sobel1,maxp1,pDACT1)
  out
}
### sample size grid 

powerSim = function(a1,b1,n,sim.num=1e4){

  res.power = matrix(NA,ncol=3,nrow= length(a1))
    for( i in 1:length(a1)){
    ### power 
    res = mediationSim(a=a1[i],b=b1[i],sigma.M=1,sigma.Y=2,n=n,sim.num=sim.num)
    p.sobel = Sobel(res$Z.M,res$Z.Y)
    p.MaxP = apply(cbind(res$p.M,res$p.Y),1,max)
    p.DACT = DACT(p1=res$p.M,p2 = res$p.Y)
    tmp = testPower(p.sobel,p.MaxP,p.DACT)
    res.power[i,] = tmp
    }
  res.power = rbind(res.power)
  colnames(res.power) = c("Sobel","MaxP","DACT")
  res.power 
}
####################



#####################

a1 = c(seq(0.04,0.2,length.out = 200),seq(0.2,0.5,length.out = 200))
a1
b1 = 0.04/a1
b1
a1/b1
ptm = proc.time()
power = powerSim(a1=a1,b1=b1,n=1000,sim.num=1e2)
proc.time() -ptm
power

write.csv(power, file="powerCurves.csv",quote = F)
# #### function to plot power results 
plotPower = function(power_mat,a=0.1,b=0.1){
  pnum = dim(power_mat)[2]
  matplot(power_mat[,1],power_mat[,2:pnum],type="b",pch=1:(pnum-1),col = 1:(pnum-1), ylim=c(0,1),xlab="N",ylab="Power")
  axis(side = 1, at = power_mat[,1])
  title(bquote(paste(alpha==.(a),", ",beta==.(b))))
  legend("topleft", legend = colnames(power_mat)[2:pnum], pch=1:(pnum-1),col=1:(pnum-1),cex=1,bty='n') 
}

#############
x = a1/b1
par(mfrow=c(1,1))
matplot(x,power,type="l",lty=1:3,col = 1:3,ylab="Power", xlab="Effect-size Ratio" )#plot
abline(v=1)
legend("topright", legend = c("Sobel","MaxP","DACT"), lty=1:3,col=1:3,cex=1.5) # optional legend

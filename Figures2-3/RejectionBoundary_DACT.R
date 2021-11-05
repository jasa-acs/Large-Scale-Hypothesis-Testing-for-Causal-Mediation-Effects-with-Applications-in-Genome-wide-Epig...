############################################################
## Composite Null 
############################################################
par(mfrow=c(1,1))

addPoint = function(){
  points(0,0,pch="+",col="red",cex=2)
}

addLines = function(){
  abline(h=0,col="black",lty=1)
  abline(v=0,col="black",lty=1)
  
}

############################
############################
library(graphics)
pdf("Rejection_Sobel_MaxP_DACT.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
C = 1.96
x = seq(-10,10,length=1000)
y = seq(-10,10,length=1000)
z <- outer(x,y,function(x,y) abs(x*y/sqrt(x^2 + y^2)) )
contour(x,y,z,levels=C,drawlabels=FALSE,col="red",cex.lab=1,xlim=c(-10,10),ylim=c(-10,10),xlab =expression(Z[gamma]),ylab=expression(Z[beta]) )
title( "Rejection Boundary",cex=2)
addPoint()
addLines
segments(1.96,1.96,20,1.96,col = "green",cex=2)
segments(1.96,1.96,1.96,20,col = "green",cex=2)

segments(-1.96,1.96,-20,1.96,col = "green",cex=2)
segments(-1.96,1.96,-1.96,20,col = "green",cex=2)

segments(-1.96,-1.96,-20,-1.96,col = "green",cex=2)
segments(1.96,-1.96,20,-1.96,col = "green",cex=2)

segments(-1.96,-1.96,-1.96,-20,col = "green",cex=2)
segments(1.96,-1.96,1.96,-20,col = "green",cex=2)
legend("topleft",c("Sobel","MaxP"),col = c("red","green"),lty =1,bty="n")


###############
DACT = function(x,y){
  
  w1=0.2
  w2=0.2
  w3=1-w1-w2
  p1 = pchisq(x^2,df=1,lower.tail = F)
  p2 = pchisq(y^2,df=1,lower.tail = F)
  p3 = max(p1,p2)^2
  p1*w1+p2*w2+p3*w3
}

C=0.05
n=1000
x = seq(-10,10,length=n)
y = seq(-10,10,length=n)

z = matrix(NA,ncol=n,nrow=n)
for(i in 1:n){
  for(j in 1:n){
    z[i,j] = DACT(x[i],y[j])
  }
}

contour(x,y,z,levels=C,drawlabels=FALSE,col="blue",cex.lab=1,xlim=c(-10,10),ylim=c(-10,10),xlab =expression(Z[gamma]),ylab=expression(Z[beta]) )
title( "Rejection Boundary",cex=2)
addPoint()
addLines
segments(1.96,1.96,20,1.96,col = "green",cex=2)
segments(1.96,1.96,1.96,20,col = "green",cex=2)

segments(-1.96,1.96,-20,1.96,col = "green",cex=2)
segments(-1.96,1.96,-1.96,20,col = "green",cex=2)

segments(-1.96,-1.96,-20,-1.96,col = "green",cex=2)
segments(1.96,-1.96,20,-1.96,col = "green",cex=2)

segments(-1.96,-1.96,-1.96,-20,col = "green",cex=2)
segments(1.96,-1.96,1.96,-20,col = "green",cex=2)
legend("topleft",c("DACT","MaxP"),col = c("blue","green"),lty =1,bty="n")
dev.off()
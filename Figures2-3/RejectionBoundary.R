############################################################
## Composite Null 
############################################################




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


pdf("Rejection_Sobel_MaxP.pdf")
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
legend("topleft",c("Sobel's Test","MaxP Test"),col = c("red","green"),lty =1)
dev.off()
       
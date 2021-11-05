
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
  legend("topleft",cex=1.5,legend = c("Uncorrected DACT","Corrected DACT "),pch=c(20,20),col = c("black","blue"))
}
#############


load("./res_averageCase001.Rdata")
load("./res_averageCase002.Rdata")
load("./res_worstCase.Rdata")
#############
# pdf("QQ_mixNulls.pdf",width = 12, height =4 )
# par(mfrow=c(1,3))
# par(mar=c(4.1,6.1,3.1,2.1))
# ###  for res 1
# myQQ(res_worstCase,main="w1=0.33, w2=.0.33")
# myQQ(res_averageCase002,main="w1=0.05, w2=0.05")
# myQQ(res_averageCase001,main="w1=0.01, w2=0.01")
# dev.off()


#############
png("QQ_mixNulls.png",width = 1000, height =300 )
par(mfrow=c(1,3))
par(mar=c(4.1,6.1,3.1,2.1))
###  for res 1
myQQ(res_worstCase,main="w1=0.33, w2=.0.33")
myQQ(res_averageCase002,main="w1=0.05, w2=0.05")
myQQ(res_averageCase001,main="w1=0.01, w2=0.01")
dev.off()

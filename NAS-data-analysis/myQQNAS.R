myQQNAS = function (pmatrix, main="title",...) {
  pvector1 = pmatrix[,1]
  pvector2 = pmatrix[,2]
  pvector3 = pmatrix[,3]
  
  o1 = -log10(sort(pvector1, decreasing = FALSE))
  e1= -log10(ppoints(length(pvector1)))
  
  o2 = -log10(sort(pvector2, decreasing = FALSE))
  e2= -log10(ppoints(length(pvector2)))
  
  o3 = -log10(sort(pvector3, decreasing = FALSE))
  e3= -log10(ppoints(length(pvector3)))
  plot(e3,o3,col="black",main=main,cex.main=3, pch=20,xlim = c(0, max(e3)+1), cex.lab =2.3, ylim = c(0, max(o3)+1), xlab = expression(Expected ~ ~-log[10](italic(p))), 
       ylab = expression(Observed ~ ~-log[10](italic(p))))    
  points(e2,o2,col="green",pch=20)
  points(e1,o1,col="blue",pch=20)
  abline(0, 1, col = "red")
  legend("topleft",cex=3,legend = c("DACT","MaxP","Sobel"),pch=c(20,20),col = c("black","green","blue"))
}
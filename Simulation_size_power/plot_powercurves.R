dat = read.csv("./powerCurves.csv")
head(dat)
x = dat[,2]
powerout = dat[,3:5]
pdf("powerCurves.pdf")
matplot(x,powerout,type="l",lty=1,col = 1:3,lwd=3,ylab="Power", xlab="Effect-size Ratio" ,cex.lab=1.5)#plot
abline(v=1,col="blue")
legend("topright", legend = c("Sobel","MaxP","DACT"),col=1:3,cex=1.5,lwd=3) # optional legend
dev.off()
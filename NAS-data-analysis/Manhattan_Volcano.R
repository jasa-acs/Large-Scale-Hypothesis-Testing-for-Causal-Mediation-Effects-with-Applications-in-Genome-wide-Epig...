
dat = read.csv("mval_mediation.csv",header = T)
head(dat)

results = dat[,c("probe","nie","nie_se","p.average")]
names(results) = c("probe","coef","se","pvalue")
save(results,file="mediation.RData")


###################
library(qvalue)
library(gap)
library(qqman)

#' Title
#'
#' @param path 
#' @param inputfile 
#' @param outputFolder 
#'
#' @return
#' @export
#'
#' @examples
EWAS_summary = function(path,inputfile,outputFolder){
  
  if(file.exists(outputFolder)){
    print("The folder already exists")
  }
  else{
    dir.create(outputFolder)
  }
  
  load(paste0(path,inputfile,".RData"))
  annot_file = read.table("cpg_chr_pos.txt",header=T)
  names(results)[1]="SiteID"
  dat = merge(annot_file,results,by="SiteID")
  cat("Annotation data merged!\n")
  write.csv(dat,file=paste0("./",inputfile,"_GeneLocation.csv"))
  res.q5 = qvalue(dat$pvalue,fdr.level = 0.05)
  res.q10 = qvalue(dat$pvalue,fdr.level = 0.1)
  dat$qvalueSig5 = res.q5$significant
  dat$qvalueSig10 = res.q10$significant
  write.csv(dat[dat$qvalueSig10,],file=paste0(outputFolder,"/",inputfile,"_Sig.csv"))
  ## top 50 
  top.sorted = dat[order(dat$pvalue),]
  write.csv(top.sorted[1:50,],file=paste0(outputFolder,"/",inputfile,"_Top50.csv"))
  
  cat("Starting Manhattan plot\n")
  ## Manhattan plot 
  dat = read.csv("mediation_GeneLocation.csv",header = T,as.is = T) ## as.is =T is important, keep CHR as character, not Factor!!
  #head(dat)
  #table(dat$Chr)
  #class(dat$Chr)
  dat$Chr[dat$Chr== "X"]= "23"
  dat$Chr[dat$Chr== "Y"]= "24"
  #table(dat$Chr)
  res.for.manhattan = dat[,c("SiteID","Chr","Loc","pvalue")]
  ## use some GWAS plot function...
  names(res.for.manhattan) = c("SNP","CHR","BP","P")
  #head(res.for.manhattan)
  chr  = res.for.manhattan$CHR
  n.cpg = dim(res.for.manhattan)[1]
  res.for.manhattan[,"CHR"] = as.numeric(chr)
  #table(res.for.manhattan$CHR)
  #res.for.manhattan$CHR = as.numeric(res.for.manhattan$CHR) ## this line of code is WRONG!!!!, needs to check!!!
  #head(res.for.manhattan[order(res.for.manhattan$P),])
  genome.sig.level = 0.05/n.cpg
  jpeg(paste0(outputFolder,"/",inputfile,"Manhattan.jpg"),height=480,width=800)
  manhattan(res.for.manhattan,cex = 2, cex.axis = 1,col = c("blue4", "orange3"), suggestiveline = T, genomewideline = FALSE)
  abline(h=-log10(genome.sig.level),col="red")
  dev.off()
  
  cat("Starting Volcano plot\n")
  
  ### Volcano Plot
  xmax = max(abs(dat$coef))
  jpeg(paste0(outputFolder,"/",inputfile,"Volcano.jpg"))
  plot(dat$coef,-log10(dat$pvalue),xlim=c(-xmax,xmax),xlab="Indirect Effect",ylab="-log10(P-value)",main = "Volcano Plot")
  abline(h=-log10(genome.sig.level),col="red")
  #title(inputfile)
  dev.off()
  
  # cat("Starting QQ plot!\n")
  # ## QQ plot
  # qqplot = function(pvector, main=NULL, ...) {
  #   o = -log10(sort(pvector,decreasing=F))
  #   e = -log10( 1:length(o)/length(o) )
  #   plot(e,o,pch=19,cex=1, main=main, ...,
  #        xlab=expression(Expected~~-log[10](italic(p))),
  #        ylab=expression(Observed~~-log[10](italic(p))),
  #        xlim=c(0,max(e)), ylim=c(0,max(o)))
  #   lines(e,e,col="red")
  # }
  # 
  # lambda.GC = round(estlambda(dat$pvalue)$estimate,2)
  # jpeg(paste0(outputFolder,"/",inputfile,"QQ.jpg"))
  # qqplot(dat$pvalue,main = bquote(lambda == .(lambda.GC)))
  # dev.off()
  # 
}

###############

path = "./"
file.names <- dir(path, pattern =".RData")
for(i in 1:length(file.names)){
  str.len = nchar(file.names[i])
  inputfile = substr(file.names[i],start=1,stop=(str.len-6))
  print(inputfile)
  EWAS_summary(path,inputfile,outputFolder ="output")
}


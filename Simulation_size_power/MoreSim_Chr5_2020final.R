### MoreSimulation using Chr 5

# #######################
# outall = read.csv("./outall.csv",header = T)
# outchr5 = outall[outall$Chr=="5",]
# hist(outchr5$Pval.x)
# hist(outchr5$Pval.y)
# qq(outchr5$Pval.y)
# qq(outchr5$Pval.x)
# qq(outchr5$DACT_correct,main="Corrected DACT")
# topCpG = outchr5[order(outchr5$Pval.y),"SiteID"] ## take top 500 CpG sites 
# ####################
# load("./data/mdata.Rdata")
# load("./data/pdata.Rdata")
# annot_file = read.table("./data/cpg_chr_pos.txt",header=T)
# head(annot_file)
# dim(annot_file[annot_file$Chr=="5",])
# probe.chr = annot_file[annot_file$Chr=="5",]$SiteID
# mybeta <- t(mdata[row.names(mdata) %in%probe.chr,])
# dim(mybeta)
# mybeta.sort = mybeta[,topCpG] ## sort by M-Y signals
# n = dim(mybeta)[1]
# n.cpg = dim(mybeta)[2]
# #########
# sum(pdata$samplename==row.names(mybeta))

#####################
## simulation starts 
#####################

library(qqman)
library(locfdr)
library(qvalue)
source("./EstNull.func.R")
source("epsest.func.R")
load("./coefy.Rdata")
load("./gamma.size.Rdata")
sim_fdr_power = function(){
  ##################
  ###############
  ## binary A
  load("./Mvalue.Rdata")
  n = dim(Mvalue)[1]
  m = dim(Mvalue)[2]
  A = rbinom(n,size = 1,prob = 0.5) 
  # A affects CpG sites by changing their means, not covariance among them
  #######################
  ## gamma effect size
  #####################
  Mvalue.new = Mvalue
  for(i in 1:length(gamma.size)){
    #Mvalue.new[,i] = Mvalue[,i] + A*gamma.size[i] + rnorm(1,mean = 0,sd = 0.1)
    Mvalue.new[,i] = Mvalue[,i] + A*gamma.size[i] 
  }
  ##########
  p.gamma = rep(NA,m)
  for(i in 1:m){
    fit.gamma = lm(Mvalue.new[,i]~A)
    p.gamma[i] = summary(fit.gamma)$coef[2,4]
  }
  
  ##qqman::qq(p.gamma)
  ######################
  ## outcome 
  ######################
  Y.sim = cbind(1,A,Mvalue.new[,1:500])%*%coef.y + rnorm(n,0,1.2)
  p.beta = c()
  for(i in 1:m){
    fit = lm(Y.sim  ~ 1 + A + Mvalue.new[,i])
    p.beta[i]=summary(fit)$coef[3,4]
  }
  
  y.sig.set = which(p.beta[1:50] < 0.05)
  #hist(p.beta)
  #qq(p.beta)
  ###############
  ## DACT
  ################
  Z_a = qnorm(p.gamma,lower.tail = F)
  Z_b  = qnorm(p.beta,lower.tail = F)
  #hist(Z_a)
  #hist(Z_b)
  # estimate proportions
  mu.sigma_a = EstNull.func(Z_a)
  mu.sigma_b = EstNull.func(Z_b)
  pi0a = 1- epsest.func(Z_a,mu.sigma_a$mu,mu.sigma_a$s)
  pi0b = 1- epsest.func(Z_b,mu.sigma_b$mu,mu.sigma_b$s)
  #pi0a
  #pi0b
  #### calculate the weights under the null
  wg1 = pi0a*(1-pi0b)
  wg2 = (1-pi0a)*pi0b
  wg3 = pi0a*pi0b
  wg4 = (1-pi0a)*(1-pi0b)
  #c(wg1,wg2,wg3,wg4)
  #round(c(wg1,wg2,wg3,wg4),5)
  ### normalize under the null
  wg.sum = wg1 + wg2 + wg3
  wg.std = c(wg1,wg2,wg3)/wg.sum
  #wg.std
  p.mat = cbind(p.gamma,p.beta)
  maxp = apply(p.mat,1,max)
  p3 = (maxp)^2
  p.average = wg.std[1]*p.gamma + wg.std[2]*p.beta + wg.std[3]*p3
  #hist(p.average)
  #qq(p.average,main="DACT")
  
  z.dact = qnorm(p.average)
  res = EstNull.func(z.dact)
  p.dact_correct = pnorm(z.dact,mean = res$mu,res$s)
  sum(p.dact_correct < 0.05/m) ## FWER
  
  fdr.thres = 0.05
  Pos.index = which(qvalue(p.dact_correct)$qvalues <fdr.thres)
  Pos.num = sum(qvalue(p.dact_correct)$qvalues < fdr.thres)
  True.pos = length(intersect(which(qvalue(p.dact_correct)$qvalues < fdr.thres),y.sig.set)) # the first 50 are signals
  FDP = 1 - True.pos/Pos.num
  TPP = True.pos/length(y.sig.set)
  #############3
  FWER = sum(which(p.dact_correct < 0.05/m) > 500)
  FWER_TPP = length(intersect(which(p.dact_correct < 0.05/m),y.sig.set))/length(y.sig.set)
  ############
  ## maxp 
  ############
  FWER_maxp = sum(which(maxp < 0.05/m) > 500)
  FWER_TPP_maxp = length(intersect(which(maxp < 0.05/m),y.sig.set))/length(y.sig.set)
  Pos.index_maxp = which(qvalue(maxp)$qvalues <fdr.thres)
  Pos.num_maxp = sum(qvalue(maxp)$qvalues < fdr.thres)
  True.pos_maxp = length(intersect(which(qvalue(maxp)$qvalues < fdr.thres),y.sig.set)) # the first 50 are signals
  FDP_maxp = 1 - True.pos_maxp/Pos.num_maxp
  TPP_maxp = True.pos_maxp/length(y.sig.set)
  ##############
  c(FWER,FWER_TPP,FDP,TPP,FWER_maxp,FWER_TPP_maxp,FDP_maxp,TPP_maxp)
}

###########
t = proc.time()
n.sim = 1000
res = matrix(NA,nrow =n.sim ,ncol = 8)
for(j in 1:n.sim){
  res[j,] = sim_fdr_power()
  cat(res[j,],"\n")
}
proc.time() - t

write.csv(res,file = "simChr5.csv",quote = F)




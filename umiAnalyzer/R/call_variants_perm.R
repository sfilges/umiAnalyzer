library(dplyr)
library(VGAM)
args <- commandArgs(trailingOnly=TRUE)
filename <- args[1] #cons file name from umierrocorrect
outfilename <- args[2] #outfilename e.g. variants.txt

#' Calculate prior distribution.
#' param params Non-negative parameters of the Beta distribution.
#' data consensus.data table of a UMisample or UMIexperiment object.
betaNLL <- function(params,data){
  a<-params[1]
  b<-params[2]
  # negative log likelihood for beta
  return(-sum(dbeta(data,shape1=a, shape2=b, log=TRUE)))
}

#' Calculate p-value using permutation-based testing.
#' param object A UMierrorcorrect object.
calc.perm.pvalue<-function(object){
  a1<-object$Coverage*object$Max.Non.ref.Allele.Frequency
  b1<-object$Coverage
  m<-mean(a1/b1)
  v<-var(a1/b1)
  a0<-m*(m * (1-m) / v-1 )
  b0<-(1-m)*(m * (1-m) / v-1 )
  params0=c(a0,b0)

  fit <- nlm(betaNLL,params0,a0/b0)
  a<-fit$estimate[1]
  b<-fit$estimate[2]
  pval<-NULL

  for (i in 1:length(a1)){
    r1<-rbetabinom.ab(10000,b1[i],shape1=a,shape2=b)
    pval[i]=sum(r1>a1[i])/10000
  }
  return(pval)
}

x <- read.table(filename,sep="\t",header=TRUE)
x2 <- x %>% filter(Consensus.group.size==3 & Coverage >= 100 & !Name=="")
pval<-calc.perm.pvalue(x2)
fdr<-p.adjust(pval,method="fdr")
results<-cbind(x2,pval,fdr)
#cutoff is fdr<0.05
write.table(results[fdr<0.05,],outfilename,sep='\t',row.names=FALSE,quote=FALSE)


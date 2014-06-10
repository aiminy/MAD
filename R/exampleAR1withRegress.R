#####################
# This is a simple example of using "AR1withRegress.R"
#####################
library(regress)
library(rrBLUP)
source("./MADII_layout_function.R")
source("./ar1WithRegress.R")

heritability <- 0.5
noErrCov <- NULL
ysErrCov <- NULL
for (i in 1:10){
  numEntries <- 240
  dummy <- MADIIdgn(num.entries=numEntries, num.rows=6, num.cols=NULL, num.sec.chk=3, designID="dummyExpt", annoy=F)
  # Simple simulation of genotype and error effects with no spatial covariance of error
  entries <- sort(unique(dummyExpt$Entry)) # Note: if there is a Fill, it will also be given an effect
  entryEffects <- rnorm(length(entries))
  names(entryEffects) <- entries
  errVar <- (1 - heritability) / heritability
  fieldObs <- data.frame(Entry=dummyExpt$Entry, phenoVal=entryEffects[dummyExpt$Entry] + rnorm(nrow(dummyExpt), sd=sqrt(errVar)))
  kinOut <- kin.blup(fieldObs, geno="Entry", pheno="phenoVal", K=NULL)
  noErrCov <- rbind(noErrCov, c(kinOut$Vg, kinOut$Ve, cor(kinOut$g[order(names(kinOut$g))], entryEffects)))
  bestReg <- ar1optim(dummyExpt, fieldObs)
  ysErrCov <- rbind(ysErrCov, c(bestReg$regOut$sigma, cor(BLUP(bestReg$regOut, RE="Entry")$Mean, entryEffects), bestReg$phi))
}

print(colMeans(noErrCov))
print(colMeans(ysErrCov))
hist(ysErrCov[,4])
print(sum(ysErrCov[,3] > noErrCov[,3]))

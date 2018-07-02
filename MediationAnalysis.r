## This script performs the mediation analysis of DNA methylation between maternal smoking and birth weight.
## Filepaths and sample IDs have been removed for security reasons, therfore it serves as a guide to the analysis

library(bda)
library(mediation)
install.packages("bda")
install.packages("mediation")


## load methylation data 
setwd("")
load("QCed_Data.rda")

pheno<-read.csv("Register/design_minerva2012.csv", stringsAsFactors = FALSE)
### needs reformating
pheno.cont<-pheno[,c(22,18,grep("cont", colnames(pheno)))]
pheno.case<-pheno[,c(21,18,grep("case", colnames(pheno)))]
pheno.cont<-cbind(pheno.cont, "Con")
pheno.case<-cbind(pheno.case, "Case")
colnames(pheno.cont)<-gsub("_cont", "", colnames(pheno.cont))
colnames(pheno.case)<-gsub("_case", "", colnames(pheno.case))
colnames(pheno.cont)[1]<-"participiant_id"
colnames(pheno.case)[1]<-"participiant_id"
colnames(pheno.cont)[13]<-"CasCon"
colnames(pheno.case)[13]<-"CasCon"
pheno<-rbind(pheno.cont,pheno.case)
pheno$participiant_id<-gsub("MMXII_iPSYCH_", "", pheno$participiant_id)
pheno<-pheno[(match(sampleSheet$Sample_Name, pheno$participiant_id)),]

##reformat smoking data into variables
m_smoke<-rep(1, nrow(sampleSheet))
m_smoke[which(pheno$smoking == 0)]<-0
m_smoke[which(pheno$smoking == 99)]<-NA
m_smoke[which(pheno$smoking == 10)]<-NA
m_smoke[which(pheno$smoking == 11)]<-NA

pheno<-cbind(pheno, m_smoke)

counts<-read.csv("Analysis/HousemanEstimatedCellCounts.csv", stringsAsFactors = FALSE, row.names = 1)
counts<-counts[match(colnames(betas) ,rownames(counts)),]

scores<-read.csv("Analysis/QC/SmokingScores.csv", stringsAsFactors = FALSE, row.names = 1)
scores<-scores[match(colnames(betas), scores$Basename),]

## repeat EWAS of maternal smoking without birth weight covariate

res.smoke<-matrix(data = NA, ncol = 3, nrow = nrow(betas))
rownames(res.smoke)<-rownames(betas)
colnames(res.smoke)<-paste(c("P.value", "RegCoeff", "SE"), "m_smoke", sep = ".")

for(i in 1:nrow(betas)){
  model<-lm(betas[i,] ~ as.factor(m_smoke) + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) 
  + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6])
  res.smoke[i,1:3]<-summary(model)$coefficients["as.factor(m_smoke)1", c(4,1,2)])
}

## repeat EWAS of birth weight without maternal smoking covariate
res.bw<-matrix(data = NA, ncol = 6, nrow = nrow(betas))
rownames(res.bw)<-rownames(betas)
colnames(res.bw)<-as.vector(outer(c("P.value", "RegCoeff", "SE"), c("fvagt","gest_week"), paste, sep="."))

for(i in 1:nrow(betas)){
	model<-lm(betas[i,] ~ pheno[, "fvagt"] + pheno[, "gest_week"] + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) 
  + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + days2sampling$Days2Sampling)
	res.bw[i,1:3]<-summary(model)$coefficients['pheno[, "fvagt"]', c(4,1,2)]
	res.bw[i,4:6]<-summary(model)$coefficients['pheno[, "gest_week"]', c(4,1,2)]
}

res.bw<-res.bw[rownames(res.smoke),]
betas<-betas[rownames(res.smoke),]

out<-matrix(data = NA, nrow = length(which(res.smoke$P.value.m_smoke < 1e-7)), ncol = 17)
colnames(out)<-c("SmokingOnBirthweightP", "SmokingOnBirthweightBeta","SmokingP", "SmokingBeta", "BirthweightP", "BirthweightBeta", "AdjSmokingonBirthweightP", "AdjSmokingonBirthweightBeta", "SobelP", "ACMEP", "ACMEBeta","ADEP", "ADEBeta", "TotalP", "TotalBeta", "PropP", "Prop")
rownames(out)<-rownames(res.smoke)[which(res.smoke$P.value.m_smoke < 1e-7)]

### test mediation effects of each smoking associated CpG on birth weight
model<-lm(pheno$fvagt ~ pheno$m_smoke + as.factor(sampleSheet$gender.x.chr))
out[,1]<-summary(model)$coefficients["pheno$m_smoke",4]
out[,2]<-summary(model)$coefficients["pheno$m_smoke",1]

for(i in which(res.smoke$P.value.m_smoke < 1e-7)){

	## check methylation correlates with birth weight
	out[rownames(res.smoke)[i],3:4]<-unlist(res.smoke[i,c("P.value.m_smoke", "RegCoeff.m_smoke")])
	out[rownames(res.smoke)[i],5:6]<-unlist(res.bw[i,c("P.value.fvagt", "RegCoeff.fvagt")])
	
  ## association decreases upon addition of mediator
	model<-lm(pheno$fvagt ~ pheno$m_smoke + as.factor(sampleSheet$gender.x.chr) + betas[rownames(res.smoke)[i],])
	out[rownames(res.smoke)[i],7:8]<-summary(model)$coefficients["pheno$m_smoke",c(4,1)]
	out[rownames(res.smoke)[i],9]<-mediation.test(betas[rownames(res.smoke)[i],!is.na(pheno$m_smoke)],pheno$fvagt[!is.na(pheno$m_smoke)], pheno$m_smoke[!is.na(pheno$m_smoke)])["p.value", "Sobel"]

	## mediation package
	## ACME - average causal mediation effects; ADE - average direct effects
	meth<-betas[rownames(res.smoke)[i],]
	med.fit<-lm(meth ~ pheno$m_smoke + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6])
	out.fit<-lm(pheno$fvagt ~ meth + pheno$m_smoke + as.factor(sampleSheet$gender.x.chr) + pheno$gest_week)
	med.out <- summary(mediate(med.fit, out.fit, treat = "pheno$m_smoke", mediator = "meth", robustSE = TRUE, sims = 100))
	out[rownames(res.smoke)[i],c(10:17)]<-c(med.out$d.avg.p, med.out$d.avg,  med.out$z.avg.p, med.out$z.avg,  med.out$tau.p, med.out$tau.coef, med.out$n.avg.p,  med.out$n.avg)
}

out<-cbind(out, res.smoke[rownames(out),c("chr", "pos", "UCSC_RefGene_Name")])
write.csv(out, "MediationAnalysisForBirthweightandMaternalSmoking.csv")

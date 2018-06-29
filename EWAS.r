## This script uses the output of the QC.r & methylationDerivedVariables.r scripts from the MinervaASDEWAS repo to perform the primary analysis of a genome-wide analysis to identify dna methylation sites associated with ASD status. 
## Filepaths and sample IDs have been removed for security reasons, therfore it serves as a guide to the analysis

uniqueAnno<-function(row){
if(row != ""){
	return(paste(unique(unlist(strsplit(row, "\\;"))), collapse = ";"))
	} else {
	return(row)
	}
	}
  
library(qqman)

setwd("")
load("QCed_Data.rda")

crosshyb<-read.csv("CrossHybridisingProbesPriceORWeksberg.csv", row.names = 1)
probes<-read.csv("SNPsinProbesAnno.csv", row.names = 1)

## remove cross hybridising probes
remove<-match(crosshyb[,1], rownames(betas))
remove<-remove[which(is.na(remove) != TRUE)]
betas<-betas[-remove,]
## remove SNP probes
probes<-probes[row.names(betas),]
betas<-betas[which(probes$Weksburg_CommonSNP_Af_within10bpSBE == "" | probes$Illumina_CommonSNP_Af_within10bpSBE == ""),]
betas<-betas[-grep("rs", rownames(betas)),]

## load cellular composition variables
counts<-read.csv("Analysis/HousemanEstimatedCellCounts.csv", stringsAsFactors = FALSE, row.names = 1)
counts<-counts[match(colnames(betas) ,rownames(counts)),]

### alternatively use smokingScores
scores<-read.csv("Analysis/QC/SmokingScores.csv", stringsAsFactors = FALSE, row.names = 1)
scores<-scores[match(colnames(betas), scores$Basename),]

days2sampling<-read.csv("", stringsAsFactors = FALSE, row.names = 1)
days2sampling<-days2sampling[match(sampleSheet$Sample_Name, rownames(days2sampling)),]

pheno<-read.csv("", stringsAsFactors = FALSE)
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


## EWAS of birth weight and gestational age simulatenously
res<-matrix(data = NA, ncol = 6, nrow = nrow(betas))
rownames(res)<-rownames(betas)
colnames(res)<-as.vector(outer(c("P-value", "RegCoeff", "SE"), c("fvagt","gest_week"), paste, sep=":"))

for(i in 1:nrow(betas)){
	model<-lm(betas[i,] ~ pheno[, "fvagt"] + pheno[, "gest_week"] 
    + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + scores$scores_combined_A
    + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + days2sampling$Days2Sampling)
	res[i,1:3]<-summary(model)$coefficients['pheno[, "fvagt"]', c(4,1,2)]
	res[i,4:6]<-summary(model)$coefficients['pheno[, "gest_week"]', c(4,1,2)]
}

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations)
probeAnnot<-probeAnnot[rownames(betas),]
res<-cbind(res, probeAnnot)
probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other)
probeAnnot<-probeAnnot[rownames(betas),]
res<-cbind(res, probeAnnot)

write.csv(res, "Analysis/EWAS_BW_GA.csv")

## QQ plots
pdf("Analysis/QQPlot.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
par(mar = c(5,5,4,0.2))
qqman::qq(res[,1], main = "Birth weight", cex.axis = 2, cex.lab = 2, cex.main = 2)
qqman::qq(res[,4], main = "Gestational age", cex.axis = 2, cex.lab = 2, cex.main = 2)
dev.off()

## manhattan plots
res<-res[which(res$chr != "Y")]
res$chr<-gsub("chr", "", res$chr)
res$chr[which(res$chr == "X")]<-23
res$chr<-as.numeric(res$chr)
par(mar = c(5,5,2,1))
manhattan(res, p = "P.value.fvagt", chr = "chr", bp = "pos", suggestiveline = -log10(5e-5), genomewideline = -log10(1e-7), main = "Birth weight")
manhattan(res, p = "gest_weeks.P.value", chr = "chr", bp = "pos", suggestiveline = -log10(5e-5), genomewideline = -log10(1e-7), main = "Gestational age")


## scatterplots of top associations
par(mfrow = c(2,5))
for(i in order(res[,1])[1:10]){
	plot(pheno[, "fvagt"], betas[i,]*100, pch = 16, main = rownames(betas)[i], xlab = "Birth weight", ylab = "% DNA methylation")
}

par(mfrow = c(2,5))
for(i in order(res[,4])[1:10]){
	plot(pheno[, "gest_week"], betas[i,]*100, pch = 16, main = rownames(betas)[i], xlab = "Gestational age (weeks)", ylab = "% DNA methylation")
}

## EWAS of maternal smoking
res<-matrix(data = NA, nrow = nrow(betas), ncol = 3)
colnames(res)<-c("P-value", "RegCoeff", "SE")
rownames(res)<-rownames(betas)

for(i in 1:nrow(betas)){
	model<-lm(betas[i,] ~ as.factor(pheno[,"m_smoke"]) + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) 
  + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + pheno$fvagt + pheno$gest_week)
	res[i,]<-summary(model)$coefficients["as.factor(pheno[, each])1", c(4,1,2)])
}

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations)
probeAnnot<-probeAnnot[rownames(betas),]
res<-cbind(res, probeAnnot)
probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other)
probeAnnot<-probeAnnot[rownames(betas),]
res<-cbind(res, probeAnnot)

write.csv(res, "Analysis/EWAS_MaternalSmoking.csv")

## QQ plot
qqman::qq(res[,1])

## Manhattan plot
manhattan(res, p = "P.value.m_smoke", chr = "chr", bp = "pos", suggestiveline = -log10(5e-5), genomewideline = -log10(1e-7))

## boxplots of top ranked associations
par(mfrow = c(2,5))
for(i in order(res$P.value.m_smoke)[1:10]){
	boxplot(betas[rownames(res)[i],] ~ as.factor(m_smoke), ylab = "DNA methylation (%)", main = paste(rownames(res)[i], "\n", uniqueAnno(res$UCSC_RefGene_Name[i]), sep = ""), cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, names = c("N", "Y"), xlab = "Maternal smoking")
	title(main = paste("p = ", signif(res$P.value.m_smoke[i],3), sep = ""), line = 0, adj = 1)
}

## Interaction with ASD

runIntModel<-function(row, pheno, scores, sampleSheet, counts){
	out<-rep(NA, 18)
	model<-lm(row ~ pheno[, "fvagt"] + pheno[, "gest_week"] + scores$scores_combined_A + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + pheno$CasCon + pheno[, "fvagt"]*pheno$CasCon)
	out[1:3]<-summary(model)$coefficients['pheno[, "fvagt"]', c(4,1,2)]
	out[4:6]<-summary(model)$coefficients['pheno[, "fvagt"]:pheno$CasConCase', c(4,1,2)]
	
	model<-lm(row ~ pheno[, "fvagt"] + pheno[, "gest_week"] + scores$scores_combined_A + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + pheno$CasCon + pheno[, "gest_week"]*pheno$CasCon)
	out[7:9]<-summary(model)$coefficients['pheno[, "gest_week"]', c(4,1,2)]
	out[10:12]<-summary(model)$coefficients['pheno[, "gest_week"]:pheno$CasConCase', c(4,1,2)]
	
	model<-lm(row ~ pheno[, "fvagt"] + pheno[, "gest_week"] + scores$scores_combined_A + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + pheno$CasCon + scores$scores_combined_A*pheno$CasCon)
	out[13:15]<-summary(model)$coefficients['scores$scores_combined_A:pheno$CasConCase', c(4,1,2)]
	out[16:18]<-summary(model)$coefficients['scores$scores_combined_A:pheno$CasConCase', c(4,1,2)]
	
	return(out)
}


library(parallel)
cl <- makeCluster(getOption("cl.cores", 2))


res<-t(parApply(cl, betas, 1, runIntModel, pheno, scores, sampleSheet, counts))
colnames(res)<-as.vector(outer(c("P-value", "RegCoeff", "SE", "Int:P-value", "Int:RegCoeff", "Int:SE"), c("fvagt","gest_week", "matenal_smoking"), paste, sep=":"))

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations)
probeAnnot<-probeAnnot[rownames(betas),]
res<-cbind(res, probeAnnot)
probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other)
probeAnnot<-probeAnnot[rownames(betas),]
res<-cbind(res, probeAnnot)

write.csv(res, "Analysis/EWAS_GA_BW_MS_InteractionwithASD.csv")


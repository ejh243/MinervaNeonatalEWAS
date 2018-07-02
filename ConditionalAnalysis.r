## This script takes the results of the EWAS.r analysis and looks for multiple associations within specifc genomic window. If these exist
## conditional analysis is performed. 
## Filepaths and sample IDs have been removed for security reasons, therfore it serves as a guide to the analysis

clump<-function(dat, window = 5000, thres1 = 1e-7, thres2 = 5e-5, listProbes = TRUE){
#' takes all significant association less than thres1 in order of significance as an index site and identifies all other associations within specified window less than thres2. If a less significant probe is clumped with a more significant probe it is excluded for consideration as an index site. 
#'
#' @param dat A matrix with the columns ids, pval, chr and pos,
#' @param window Distance in bp of region around site to look for other sinificant associations.
#' @param thres1 P value threshold to select significant sites to clump.
#' @param thres2 P value threshold to select additional sites to clump with site under consideration
#' @return a ma 

	tmp<-dat[order(dat$pval),]
	out<-NULL
	while(tmp[1,2] < thres1){
		chrom<-tmp[1,3]
		start<-tmp[1,4]-window
		end<-tmp[1,4]+window
		index<-which(tmp[,3] == chrom & tmp[,4] <= end & tmp[,4] >= start & tmp[,4] != tmp[1,4])	## remove self from consideration

		n.sig<-length(which(tmp[index,2] < thres2))
		if(n.sig > 0){
			probes<-tmp[index,1][which(tmp[index,2] < thres2)]
			if(n.sig > 1){
				probes<-paste(probes, collapse = ";")
				}
		} else {
			probes<-NA
		}
		out<-rbind(out, c(unlist(tmp[1,]), length(index), n.sig, probes))
		## remove these probes from furture consideration
		tmp<-tmp[-c(1,index),]
	}
	colnames(out)<-c("Index", "Pvalue", "Chr", "Position", "nSitesinWindow", "nSigSites", "SigDNAmSites")
	return(out)
	
}

setwd("")

load("QCed_Data.rda")

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
m_smoke<-rep(1, nrow(sampleSheet))
m_smoke[which(pheno$smoking == 0)]<-0
m_smoke[which(pheno$smoking == 99)]<-NA
m_smoke[which(pheno$smoking == 10)]<-NA
m_smoke[which(pheno$smoking == 11)]<-NA

## estimated smoking scores
scores<-read.csv("", stringsAsFactors = FALSE, row.names = 1)
scores<-scores[match(colnames(betas), scores$Basename),]

## estimated cellular composition variables
counts<-read.csv("", stringsAsFactors = FALSE, row.names = 1)
counts<-counts[match(colnames(betas) ,rownames(counts)),]

## load birth weight and gestational age EWAS results
## take each association in turn (order by significance) and identify any others within 5kb
res<-read.csv("Analysis/EWAS_BW_GA.csv", stringsAsFactors = FALSE, row.names = 1)

res<-res[which(res$chr != "chrY"),]
window<-5000

## Gestational age EWAS
dat<-data.frame("ids" = rownames(res), "pval" = res$P.value.gest_week, "chr" = res$chr, "pos" = res$pos, stringsAsFactors = FALSE)
clumpGA<-clump(dat, thres2 = 1e-7)

summary(as.numeric(clumpGA[,6]))
summary(as.numeric(clumpGA[,6])/as.numeric(clumpGA[,5]))

## conditional analysis
out.cond<-NULL
for(i in which(as.numeric(clumpGA[,6]) > 0)){
	probes<-c(clumpGA[i,1], unlist(strsplit(clumpGA[i,7], ";")))
	betas.sub<-betas[probes,]
	## perform conditional analysis including less significant terms
	out.sub<-matrix(data = NA, ncol = 3, nrow = nrow(betas.sub)-1)
	rownames(out.sub)<-probes[-1]
	for(j in 2:nrow(betas.sub)){
		model<-lm(betas.sub[j,] ~ betas.sub[1,] + pheno[, "fvagt"] + pheno[, "gest_week"] + as.factor(sampleSheet$gender.x.chr) 
    + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] 
    + days2sampling$Days2Sampling + scores$scores_combined_A)
		out.sub[j-1,]<-summary(model)$coefficients['pheno[, "gest_week"]', c(4,1,2)]
	}
	out.sub<-cbind(clumpMS[i,1], out.sub)
	rownames(out.sub)<-probes[-1]
	out.cond<-rbind(out.cond, out.sub)
	
} 

## summarize results
length(which(as.numeric(out.cond[,2]) < 5e-5))
length(unique(out.cond[which(as.numeric(out.cond[,2]) < 5e-5), 1]))
length(which(as.numeric(out.cond[,2]) > 0.05))
length(unique(out.cond[which(as.numeric(out.cond[,2]) > 0.05), 1]))
length(which(as.numeric(out.cond[,2]) < 5e-5 & res$P.value.gest_week[match(rownames(out.cond), rownames(res))] < 1e-7))
length(which(as.numeric(out.cond[,2]) > 0.05 & res$P.value.gest_week[match(rownames(out.cond), rownames(res))] < 1e-7))

write.csv(out.cond, "Analysis/ConditionalAnalysis_within5kb_GA.csv")
write.csv(clumpGA, "Analysis/ClumpedLoci_within5kb_GA.csv")

### Birth weight
dat<-data.frame("ids" = rownames(res), "pval" = res$P.value.fvagt, "chr" = res$chr, "pos" = res$pos, stringsAsFactors = FALSE)
clumpBW<-clump(dat)

## as all DMPs located in distinct loci no need for conditional analysis.

### Maternal smoking

res<-read.csv("Analysis/EWAS_MaternalSmoking.csv", stringsAsFactors = FALSE, row.names = 1)
res<-res[which(res$chr != "chrY"),]

dat<-data.frame("ids" = rownames(res), "pval" = res$P.value.m_smoke, "chr" = res$chr, "pos" = res$pos, stringsAsFactors = FALSE)
clumpMS<-clump(dat, thres2 = 1e-7)

## conditional analysis
m_smoke<-as.factor(m_smoke)
out.cond<-NULL
for(i in which(as.numeric(clumpMS[,6]) > 0)){
	probes<-c(clumpMS[i,1], unlist(strsplit(clumpMS[i,7], ";")))
	betas.sub<-betas[probes,]
	## perform conditional analysis including less significant terms
	out.sub<-matrix(data = NA, ncol = 3, nrow = nrow(betas.sub)-1)
	
	for(j in 2:nrow(betas.sub)){
		model<-lm(betas.sub[j,] ~ betas.sub[1,] + pheno[, "fvagt"] + pheno[, "gest_week"] + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + m_smoke)
		out.sub[j-1,]<-summary(model)$coefficients['m_smoke1', c(4,1,2)]
	}
	
	out.sub<-cbind(clumpMS[i,1], out.sub)
	rownames(out.sub)<-probes[-1]
	out.cond<-rbind(out.cond, out.sub)
	
} 

## summarize results
length(which(as.numeric(out.cond[,2]) < 5e-5))
length(unique(out.cond[which(as.numeric(out.cond[,2]) < 5e-5), 1]))
length(which(as.numeric(out.cond[,2]) > 0.05))
length(unique(out.cond[which(as.numeric(out.cond[,2]) > 0.05), 1]))
length(which(as.numeric(out.cond[,2]) < 5e-5 & res$P.value.m_smoke[match(rownames(out.cond), rownames(res))] < 1e-7))
length(which(as.numeric(out.cond[,2]) > 0.05 & res$P.value.m_smoke[match(rownames(out.cond), rownames(res))] < 1e-7))


write.csv(out.cond, "Analysis/ConditionalAnalysis_within5kb_MS.csv")
write.csv(clumpMS, "Analysis/ClumpedLoci_within5kb_MS.csv")

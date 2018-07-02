## For CpGs meeting the crieteria for significant mediation repeat analysis using SIMEX procedure.
## This script was developed using the code provided in the Supplementary Material Section A8 of Valeri L, Reese SL, Zhao S, et al. Misclassified exposure in epigenetic mediation analyses. Does DNA methylation mediate effects of smoking on birthweight? Epigenomics. 2017;9(3):253-265. doi:10.2217/epi-2016-0145.
## Filepaths and sample IDs have been removed for security reasons, therfore it serves as a guide to the analysis

library(parallel)
library(DistributionUtils)
library(simex)

f.mediation.naive.parallel<-function(methcol,M,Astar,Y,C1,C2, C3, C4, C5, C6, C7, C8, C9,nboot=500){
				
	Mi<- M[,methcol]				
	data<- as.data.frame(cbind(Astar,Mi,Y,C1,C2, C3, C4, C5, C6, C7, C8, C9))
	#objects	to	save
	nde_naive<- rep(NA,nboot)
	nie_naive_prod<- rep(NA,nboot)
	pm_naive_prod<- rep(NA,nboot)
								
	for	(i in 1:nboot){
		if(is.wholenumber(i/(nboot/10))){	#	print	progress
			cat("nboot=",i,"	;	",sep="")
		}
		set.seed(1234+i)
		databoot<- data[sample(nrow(data),replace=TRUE),]
		databoot$Astar<- as.factor(databoot$Astar)
												
		#linear	outcome	regression
		regym<-glm(Y~Astar+Mi+C1+C2+C3+C4+C5+C6+C7+C8+C9,data=databoot,x=TRUE,y=TRUE)
		regy<-glm(Y~Astar+C1+C2,data=databoot,x=TRUE,y=TRUE)
		#linear	mediator	regression
		regm<-glm(Mi~Astar+C1+C2+C3+C4+C5+C6+C7+C8+C9,data=databoot,x=TRUE,y=TRUE)
												
		#naive	effects
		nde_naive[i]<-regym$coeff[2]
		nie_naive_prod[i]<-regym$coeff[3]*regm$coeff[2]
		pm_naive_prod[i]<-nie_naive_prod[i]/(nie_naive_prod[i]+nde_naive[i])
	}
	return(list(Naive=data.frame(nde_naive=nde_naive,
			nie_naive_prod=nie_naive_prod,
			pm_naive_prod=pm_naive_prod)))
}

f.CI.quantile<-function(x,simex=FALSE){
	#Naive
	ci_naive<-t(apply(x$Naive,2,quantile,probs=c(0.5,0.025,0.975),na.rm=T))
	#Simex
	if(simex == TRUE){
		ci_simex<-t(apply(x$Simex,2,quantile,probs=c(0.5,0.025,0.975),na.rm=T))
		out<-list(Naive=ci_naive,Simex=ci_simex)
	}else{
		out<-list(Naive=ci_naive)
	}
	return(out)
}
f.mediation.simex.parallel	<-function(methcol,M,Astar,Y,C1,C2,C3,C4,C5,C6,C7, C8, C9,
	sn=0.8,sp=1,nboot=500){
	Mi<-M[,methcol]
	data<-as.data.frame(cbind(Astar,Mi,Y,C1,C2,C3,C4,C5,C6,C7, C8, C9))
				
	#objects	to	save
	nde_naive<-rep(NA,nboot)
	nie_naive_prod<-rep(NA,nboot)
	pm_naive_prod<-rep(NA,nboot)
								
	nde_simex<- matrix(nrow=nboot,ncol=length(sn))
	nie_simex_prod<- matrix(nrow=nboot,ncol=length(sn))
	pm_simex_prod<- matrix(nrow=nboot,ncol=length(sn))
								
	for	(i in 1:nboot){
		if(is.wholenumber(i/(nboot/10)))	{
		# print progress
			cat("nboot=",i,"	;	",sep="")
		}
		set.seed(1234+i)
		databoot<- data[sample(nrow(data),replace=TRUE),]
		databoot$Astar<- as.factor(databoot$Astar)
												
		#linear	outcome	regression
		regym<-glm(Y~Astar+Mi+C1+C2+C3+C4+C5+C6+C7+C8+C9,data=databoot,x=TRUE,y=TRUE)
		regy<-glm(Y~Astar+C1+C2,data=databoot,x=TRUE,y=TRUE)
		#linear	mediator regression
		regm<-glm(Mi~Astar+C1+C2+C3+C4+C5+C6+C7+C8+C9,data=databoot,x=TRUE,y=TRUE)												
		#naive	effects
		nde_naive[i]<-regym$coeff[2]
		nie_naive_prod[i]<-regym$coeff[3]*regm$coeff[2]
		pm_naive_prod[i]<-nie_naive_prod[i]/(nie_naive_prod[i]+nde_naive[i])
		for	(k in 1:length(sn)){	
			###SIMEX PROCEDURE###
			Pi	<- matrix(data=c(sp,1-sp,1-sn[k],sn[k]),nrow = 2,byrow= FALSE)
			dimnames(Pi)<- list(levels(databoot$Astar),	levels(databoot$Astar))
			regy_simex	<- mcsimex(regy,mc.matrix=Pi,SIMEXvariable="Astar")
			regym_simex	<- mcsimex(regym,mc.matrix=Pi,SIMEXvariable="Astar")
			regm_simex	<- mcsimex(regm,mc.matrix=Pi,SIMEXvariable="Astar")
			#corrected	effects
			nde_simex[i,k]<-regym_simex$coeff[2]
			nie_simex_prod[i,k]<-regym_simex$coeff[3]*regm_simex$coeff[2]
			pm_simex_prod[i,k]<-nie_simex_prod[i,k]/(nie_simex_prod[i,k]+nde_simex[i,k])
		}
	}
	colnames(nde_simex)<-paste("nde_simex_",sn,sep="")
	colnames(nie_simex_prod)<-paste("nie_simex_prod_",sn,sep="")
	colnames(pm_simex_prod)<-paste("pm_simex_prod_",sn,sep="")
								
	return(list(Naive=data.frame(nde_naive=nde_naive,nie_naive_prod=nie_naive_prod,pm_naive_prod=pm_naive_prod),Simex=data.frame(nde_simex,nie_simex_prod,pm_simex_prod)))
}
f.descriptives<-function(x,...){
	c(mean=mean(x,...),quantile(x,...))
}
f.CI.quantile<-function(x,simex=FALSE,...){
	# Naive
	ci_naive<-t(apply(x$Naive,2,f.descriptives,...))
	# Simex
	if(simex==TRUE){
		ci_simex<- t(apply(x$Simex,2,f.descriptives,...))
		out<-list(Naive=ci_naive,Simex=ci_simex)
	}else{
		out<-list(Naive=ci_naive)
	}
	return(out)
}

## only apply to those with significant sobel test
setwd("")
sobel<-read.csv("MediationAnalysisForBirthweightandMaternalSmoking.csv", stringsAsFactors = FALSE, row.names = 1)
probes<-rownames(sobel)[which(sobel[,5] < 1e-7 & sobel[,7] > sobel[1,1] & sobel[,9] < 0.05/143)]

## load methylation data
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

#observed	exposure:	
Astar<-as.factor(pheno$m_smoke)	#self-reported	smoking

#mediator
tdat<-t(betas[probes,])	#methylation	beta	values
#outcome
Y<-pheno$fvagt #birth	weight
#confounders
C1<-as.factor(sampleSheet$gender.x.chr)
C2<-pheno$gest_week
C3<-as.factor(sampleSheet$Sentrix_ID)
C4<-counts[,1]
C5<-counts[,2]
C6<-counts[,3]
C7<-counts[,4]
C8<-counts[,5]
C9<-counts[,6]

### Simex
specificity<-1
sensitivity<- seq(0.6,0.9,0.1)
B=100 #	number	of	bootstraps

detectCores()
system.time(
	results	<- mclapply(setNames(seq_len(ncol(tdat)), dimnames(tdat)[[2]]),	
  f.mediation.simex.parallel,M=tdat,Astar,Y,C1,C2,C3,C4,C5,C6,C7, C8, C9,sn=sensitivity,sp=specificity,nboot=B,mc.cores=detectCores())
)

conf_int<-lapply(results,f.CI.quantile,simex=TRUE,na.rm=T,probs=c(0.5,0.025,0.975, 8.33e-3, 0.99166))

extract.naive <- function(Naive) lapply(conf_int, function(y) as.numeric(t((y[[Naive]]))))

## convert to data.frame for saving
df.naive <- do.call("rbind", extract.naive("Naive"))
colnames(df.naive)<-c(paste("nde", colnames(conf_int[[1]]$Naive), sep = ""), paste("nie", colnames(conf_int[[1]]$Naive), sep = ""), paste("pm", colnames(conf_int[[1]]$Naive), sep = ""))
colnames(df.naive)<-paste("Naive:", colnames(df.naive))
df.simex <- do.call("rbind", extract.naive("Simex"))
colnames(df.simex)<-rep(c(paste("nde", colnames(conf_int[[1]]$Naive), sep = ""), paste("nie", colnames(conf_int[[1]]$Naive), sep = ""), paste("pm", colnames(conf_int[[1]]$Naive), sep = "")), 4)
colnames(df.simex)<-paste("SIMEX:", colnames(df.simex))
df<-cbind(df.naive, df.simex)

write.csv(df, "SIMEX_Mediation_Results.csv")

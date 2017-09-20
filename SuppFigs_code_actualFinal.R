#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy v1 -- Upham et al. 2017
###
# Supplemental Figures 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~








#######
######
# Also calculate for FAMILIES
####
library(moments)
library(nlme)
library(ape)
library(picante)
library(phytools)
library(geiger)

setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/Diversification_analyses-Condamine/")
source("DR_functions.R")
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut")

library(foreach);library(doSNOW)
cl = makeCluster(15, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees=100

foreach(i=1:ntrees, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

# which backbone?
bbone<- "NDexp" #"FBD" # 

# load in stats about TIPS
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdGenesSampPC_DR-SUMMARY-expanded_",bbone,".txt",sep=""))
colnames(cladesDR)<-c("tiplabel","gen","fam","ord","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")
head(cladesDR)

mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_nexus.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 

# load in DR calcs
res1<-read.table(file=paste(bbone,"_sample100_",i,"_DRtips.txt",sep=""))

# root age
root=max(node.age(mamPhy)$ages)

# just get the fam tips > 2 species
allFamNames<-names(table(cladesDR$fam))
famNames<-allFamNames[which(table(cladesDR$fam) > 2)]

famTipLabels<-vector("list",length(famNames))
for (j in 1:length(famTipLabels)){
	famTipLabels[[j]]<-cladesDR[which(cladesDR$fam==famNames[j]),"tiplabel"]
}

# get the FAM TREE for each, and also the BACKBONE uniting those fams.
famClades<-vector("list",length(famNames))
for (j in 1:length(famTipLabels)){
	toDrop<-setdiff(mamPhy$tip.label,famTipLabels[[j]])
	famClades[[j]]<-drop.tip(mamPhy,toDrop)
}

toKeepFams<-vector("list",length(famTipLabels))
for (j in 1:length(famTipLabels)){ 
	toKeepFams[[j]]<-unlist(famTipLabels[[j]][1])
}
toKeep<-unlist(toKeepFams)
toDropFams<-setdiff(mamPhy$tip.label,toKeep)
famPhy<-drop.tip(mamPhy,toDropFams)

# write famTree
write.tree(famPhy, file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_FAMS_112gr2.trees",sep=""))

# get node times for tree
btimes<-branching.times(mamPhy)

# yule function
ymle = function(tree){ (.subset2(tree,3)-1L)/sum(.subset2(tree,2)) } # this take the # of number of nodes in a tree (minus 1) / sum of branch lengths.

	# empty data frames to fill
	DR_harm<-data.frame(matrix(NA, nrow = length(famClades), ncol = 1))
	DR_cv<-data.frame(matrix(NA, nrow = length(famClades), ncol = 1))
	DR_skew<-data.frame(matrix(NA, nrow = length(famClades), ncol = 1))
	DR_kurt<-data.frame(matrix(NA, nrow = length(famClades), ncol = 1))
	percentSamp<-data.frame(matrix(NA, nrow = length(famClades), ncol = 1))
	richness<-data.frame(matrix(NA, nrow = length(famClades), ncol = 1))
	MRCA<-data.frame(matrix(NA, nrow = length(famClades), ncol = 1))
	PB_Div<-data.frame(matrix(NA, nrow = length(famClades), ncol = 1))
	BD_Lam<-data.frame(matrix(NA, nrow = length(famClades), ncol = 1))
	BD_Mu<-data.frame(matrix(NA, nrow = length(famClades), ncol = 1))
	BD_Div<-data.frame(matrix(NA, nrow = length(famClades), ncol = 1))
	BD.ms0<-data.frame(matrix(NA, nrow = length(famClades), ncol = 1))
	BD.ms0p5<-data.frame(matrix(NA, nrow = length(famClades), ncol = 1))
	BD.ms0p9<-data.frame(matrix(NA, nrow = length(famClades), ncol = 1))
	BD.ms0_stem<-data.frame(matrix(NA, nrow = length(famClades), ncol = 1))
	BD.ms0p5_stem<-data.frame(matrix(NA, nrow = length(famClades), ncol = 1))
	BD.ms0p9_stem<-data.frame(matrix(NA, nrow = length(famClades), ncol = 1))
	
# do per-family calcs
for(j in 1:length(famClades)){
cladeSet<-famClades[[j]]
	cladeSp<-famClades[[j]]$tip.label
		x<-res1[match(cladeSp,rownames(res1)),"DR"]
		DR_harm[j,] <- 1/(mean(1/x))
		DR_cv[j,] <- (sd(x)/mean(x))*100
		DR_skew[j,] <- skewness(x)
		DR_kurt[j,] <- kurtosis(x)
		percentSamp[j,] <- length(which(cladesDR[match(cladeSp,cladesDR$tiplabel),"samp"]=="sampled"))/length(cladeSp)
		richness[j,] <- length(cladeSp)
		node <- getMRCA(mamPhy, cladeSp)
		MRCA[j,] <- btimes[node-5911] #taking the height of SAMPLED tree
		# Yule model
		PB_Div[j,]<-ymle(cladeSet)
		# BD model
		bd<-birthdeath(cladeSet)
		BD_Lam[j,]<-bd$para[[2]]/(1-bd$para[[1]])
		BD_Mu[j,]<-bd$para[[1]]*(bd$para[[2]]/(1-bd$para[[1]]))
		BD_Div[j,]<-bd$para[[2]]
		# BD Mag and Sand
	    BD.ms0_stem[j,]<-bd.ms(phy=cladeSet, missing=0, epsilon=0, crown=FALSE) # Assuming no extinction
    	BD.ms0p5_stem[j,]<-bd.ms(phy=cladeSet, missing=0, epsilon=0.5, crown=FALSE) # Assuming medium extinction 
     	BD.ms0p9_stem[j,]<-bd.ms(phy=cladeSet, missing=0, epsilon=0.9, crown=FALSE) # Assuming high extinction
		cladeSet$root.edge<-0
	    BD.ms0[j,]<-bd.ms(phy=cladeSet, missing=0, epsilon=0, crown=TRUE) # Assuming no extinction
    	BD.ms0p5[j,]<-bd.ms(phy=cladeSet, missing=0, epsilon=0.5, crown=TRUE) # Assuming medium extinction 
     	BD.ms0p9[j,]<-bd.ms(phy=cladeSet, missing=0, epsilon=0.9, crown=TRUE) # Assuming high extinction
} 
res2<-cbind(DR_harm, DR_cv, DR_skew, DR_kurt, percentSamp, richness, MRCA, PB_Div, BD_Lam, BD_Mu, BD_Div, BD.ms0, BD.ms0p5, BD.ms0p9, BD.ms0_stem, BD.ms0p5_stem, BD.ms0p9_stem, i)
colnames(res2)<-c("DR_harm","DR_cv", "DR_skew", "DR_kurt", "percentSamp", "richness", "MRCA", "PB_Div", "BD_Lam", "BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9", "BD.ms0_stem", "BD.ms0p5_stem", "BD.ms0p9_stem", "tree")
rownames(res2)<-famNames

write.table(res2,paste(bbone,"_sample100_",i,"_FAMS_112gr2_cladeSTATS.txt",sep=""))


}

# do the PGLS - FAMILY-level
#######
# load back in tables
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut")
bbone<- "NDexp" #"FBD" # 
library(nlme); library(MASS); library(AICcmodavg); library(ape); library(phytools); library(geiger); library(dplyr)

numTrees<-100

results<-vector("list",length=numTrees)
for (i in 1:numTrees){
	results[[i]]<-read.table(paste(bbone,"_sample100_",i,"_FAMS_112gr2_cladeSTATS.txt",sep=""))
}
allFamRes<-do.call(rbind,results)

# load back in famPhys
famPhys<-vector("list",length=numTrees)
for (i in 1:numTrees){
	famPhys[[i]]<-read.tree(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_FAMS_112gr2.trees",sep=""))
}

# UNIVARIATE first
uniPGLS_allSlopes<-data.frame(matrix(NA, nrow = length(results), ncol = 13))
colnames(uniPGLS_allSlopes)<-c("MRCA","DR_harm","DR_skew","PB_Div","BD_Lam","BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9", "BD.ms0_stem", "BD.ms0p5_stem", "BD.ms0p9_stem")

uniPGLS_allInts<-data.frame(matrix(NA, nrow = length(results), ncol = 13))
colnames(uniPGLS_allInts)<-c(paste("i_",1:13,sep=""))

uniPGLS_allPs<-data.frame(matrix(NA, nrow = length(results), ncol = 13))
colnames(uniPGLS_allPs)<-c(paste("p_",1:13,sep=""))

# MULTIVARIATE next
partSlopesPGLS_All<-data.frame(matrix(NA, nrow = length(results), ncol = 8))
colnames(partSlopesPGLS_All)<-c("int","MRCA","DR_harm","DR_skew","AIC","Pval1","Pval2","Pval3")

for (j in 1:length(results)){
	rownames(results[[j]])<-famPhys[[j]]$tip.label
	cladeData<-treedata(famPhys[[j]],na.omit(results[[j]]))
	dat<-as.data.frame(cladeData$data)

	form<-(log(richness) ~ MRCA + DR_harm + DR_skew)
	fit<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
	sum<-summary(fit)

	partSlopesPGLS_All[j,1]<-round(sum$coef[[1]],digits=3)
	partSlopesPGLS_All[j,2]<-round(sum$coef[[2]],digits=3)
	partSlopesPGLS_All[j,3]<-round(sum$coef[[3]],digits=3)
	partSlopesPGLS_All[j,4]<-round(sum$coef[[4]],digits=3)
	partSlopesPGLS_All[j,5]<-round(sum$AIC,digits=0)
	partSlopesPGLS_All[j,6]<-round(sum$tTable[14], digits=3)
	partSlopesPGLS_All[j,7]<-round(sum$tTable[15], digits=3)
	partSlopesPGLS_All[j,8]<-round(sum$tTable[16], digits=3)

	vars<-c("MRCA","DR_harm","DR_skew","PB_Div","BD_Lam","BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9", "BD.ms0_stem", "BD.ms0p5_stem", "BD.ms0p9_stem")
	for(k in 1:length(uniPGLS_allSlopes)){
		form<-as.formula(paste("log(richness) ~ ", vars[k],sep=""))
		fit1<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
		sum1<-summary(fit1)
		
		uniPGLS_allInts[j,k]<-round(sum1$coef[[1]],digits=3)
		uniPGLS_allSlopes[j,k]<-round(sum1$coef[[2]],digits=3)
		uniPGLS_allPs[j,k]<-round(sum1$tTable[8], digits=3)
	}	
}
write.table(partSlopesPGLS_All,paste(bbone,"_sample100_PGLSmulti_FAMS_112gr2.txt",sep=""))

uniPGLS_fams<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs)

write.table(uniPGLS_fams,paste(bbone,"_sample100_PGLSuniAll_FAMS_112gr2.txt",sep=""))

# Also with STANDARDIZED -- SCALED data - FAMILIES::
######
# load back in tables
results_Z<-vector("list",length=numTrees)
for (i in 1:numTrees){
	res<-read.table(paste(bbone,"_sample100_",i,"_FAMS_112gr2_cladeSTATS.txt",sep=""))
	varsScale<-scale(res, center=TRUE,scale=TRUE)
	richnessN<-res$richness
	results_Z[[i]]<-cbind(richnessN,varsScale)
}	

uniPGLS_allSlopes<-data.frame(matrix(NA, nrow = length(results_Z), ncol = 13))
colnames(uniPGLS_allSlopes)<-c("MRCA","DR_harm","DR_skew","PB_Div","BD_Lam","BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9", "BD.ms0_stem", "BD.ms0p5_stem", "BD.ms0p9_stem")

uniPGLS_allInts<-data.frame(matrix(NA, nrow = length(results_Z), ncol = 13))
colnames(uniPGLS_allInts)<-c(paste("i_",1:13,sep=""))

uniPGLS_allPs<-data.frame(matrix(NA, nrow = length(results_Z), ncol = 13))
colnames(uniPGLS_allPs)<-c(paste("p_",1:13,sep=""))

# MULTIVARIATE next
partSlopesPGLS_All<-data.frame(matrix(NA, nrow = length(results_Z), ncol = 8))
colnames(partSlopesPGLS_All)<-c("int","MRCA","DR_harm","DR_skew","AIC","Pval1","Pval2","Pval3")

for (j in 1:length(results_Z)){
	rownames(results_Z[[j]])<-famPhys[[j]]$tip.label
	cladeData<-treedata(famPhys[[j]],(results_Z[[j]]))
	dat<-as.data.frame(cladeData$data)

	form<-(log(richnessN) ~ MRCA + DR_harm + DR_skew)
	fit<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
	sum<-summary(fit)

	partSlopesPGLS_All[j,1]<-round(sum$coef[[1]],digits=3)
	partSlopesPGLS_All[j,2]<-round(sum$coef[[2]],digits=3)
	partSlopesPGLS_All[j,3]<-round(sum$coef[[3]],digits=3)
	partSlopesPGLS_All[j,4]<-round(sum$coef[[4]],digits=3)
	partSlopesPGLS_All[j,5]<-round(sum$AIC,digits=0)
	partSlopesPGLS_All[j,6]<-round(sum$tTable[14], digits=3)
	partSlopesPGLS_All[j,7]<-round(sum$tTable[15], digits=3)
	partSlopesPGLS_All[j,8]<-round(sum$tTable[16], digits=3)

	vars<-c("MRCA","DR_harm","DR_skew","PB_Div","BD_Lam","BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9", "BD.ms0_stem", "BD.ms0p5_stem", "BD.ms0p9_stem")
	for(k in 1:length(uniPGLS_allSlopes)){
		form<-as.formula(paste("log(richnessN) ~ ", vars[k],sep=""))
		fit1<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
		sum1<-summary(fit1)
		
		uniPGLS_allInts[j,k]<-round(sum1$coef[[1]],digits=3)
		uniPGLS_allSlopes[j,k]<-round(sum1$coef[[2]],digits=3)
		uniPGLS_allPs[j,k]<-round(sum1$tTable[8], digits=3)
	}	
}
write.table(partSlopesPGLS_All,paste(bbone,"_sample100_PGLSmulti_FAMS_112gr2_SCALED.txt",sep=""))

uniPGLS_fams_SCALED<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs)

write.table(uniPGLS_fams_SCALED,paste(bbone,"_sample100_PGLSuniAll_FAMS_112gr2_SCALED.txt",sep=""))


###
# PLOT to summarize

# *** ACTUAL DATA VALUES:
# PGLS uni as based on ACTUAL values (not z-scores)
results<-vector("list",length=numTrees)
for (i in 1:numTrees){
	results[[i]]<-read.table(paste(bbone,"_sample100_",i,"_FAMS_112gr2_cladeSTATS.txt",sep=""))
}
allFamRes<-do.call(rbind,results)
uniPGLS_fams<-read.table(paste(bbone,"_sample100_PGLSuniAll_FAMS_112gr2.txt",sep=""))


jpeg(file=paste("cladeLevel",bbone,"_scatter_FAMS_PGLSuni_Yaxis_PLOTTED_3part.jpg",sep=""),width=4,height=8, units="in", res=600, quality=100)
par(mfrow = c(3,1),oma = c(5,4,5,0) + 0.1, mar = c(3,1,1,1) + 0.1)
# type="n",
plot((richness)~MRCA,data=allFamRes, log="y",type="n",col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
axis()
mtext(text="Clade crown age (Ma)", side=1, padj=3, font=2,cex=1)

plot((richness)~DR_harm,data=allFamRes, log="y",type="n",col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
mtext(text="Clade ESDR mean (species / Ma)", side=1, padj=3, font=2,cex=1)

plot((richness)~DR_skew,data=allFamRes, log="y", type="n",col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
mtext(text="Clade ESDR skewness", side=1, padj=3, font=2,cex=1)

title(main="", xlab = "",
      ylab = "log(richness)",
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)
dev.off()


#pdf(file=paste("cladeLevel",bbone,"_scatter_FAMS_PGLSuni_ALL_PLOTTED_3part.pdf",sep=""),onefile=TRUE,width=4,height=10)
jpeg(file=paste("cladeLevel",bbone,"_scatter_FAMS_PGLSuni_ALL_PLOTTED_3part.jpg",sep=""),width=4,height=8, units="in", res=600, quality=100)
par(mfrow = c(3,1),oma = c(5,4,5,0) + 0.1, mar = c(3,1,1,1) + 0.1)
# type="n",
plot(log(richness)~MRCA,data=allFamRes, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
for (i in 1:length(uniPGLS_fams[,1])){
	if (uniPGLS_fams$p_1[i] < 0.05){
	abline(a=(uniPGLS_fams$i_1[i]),b=(uniPGLS_fams$MRCA[i]), untf=FALSE,col=grey(0,alpha=0.2), lwd=2,lty=2)
	} else 	abline(a=(uniPGLS_fams$i_1[i]),b=(uniPGLS_fams$MRCA[i]), untf=FALSE,col=rgb(red=1,green=0,blue=0,alpha=0.5), lwd=2,lty=2)
}
mtext(text="Clade crown age (Ma)", side=1, padj=3, font=2,cex=1)

plot(log(richness)~DR_harm,data=allFamRes, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
for (i in 1:length(uniPGLS_fams[,1])){
	if (uniPGLS_fams$p_2[i] < 0.05){
	abline(a=(uniPGLS_fams$i_2[i]),b=(uniPGLS_fams$DR_harm[i]), untf=FALSE,col=grey(0,alpha=0.2), lwd=2,lty=2)
	} else 	abline(a=(uniPGLS_fams$i_2[i]),b=(uniPGLS_fams$DR_harm[i]), untf=FALSE,col=rgb(red=1,green=0,blue=0,alpha=0.5), lwd=2,lty=2)
}
mtext(text="Clade ESDR mean (species / Ma)", side=1, padj=3, font=2,cex=1)

plot(log(richness)~DR_skew,data=allFamRes, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
for (i in 1:length(uniPGLS_fams[,1])){
	if (uniPGLS_fams$p_3[i] < 0.05){
	abline(a=(uniPGLS_fams$i_3[i]),b=(uniPGLS_fams$DR_skew[i]), untf=FALSE,col=grey(0,alpha=0.2), lwd=2,lty=2)
	} else 	abline(a=(uniPGLS_fams$i_2[i]),b=(uniPGLS_fams$DR_harm[i]), untf=FALSE,col=rgb(red=1,green=0,blue=0,alpha=0.5), lwd=2,lty=2)
}
mtext(text="Clade ESDR skewness", side=1, padj=3, font=2,cex=1)

title(main="", xlab = "",
      ylab = "log(richness)",
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)
dev.off()


##
# *** ACTUAL DATA VALUES:
#########################
# and FAMILIES -- actual values and univariate slopes
results<-vector("list",length=numTrees)
for (i in 1:numTrees){
	results[[i]]<-read.table(paste(bbone,"_sample100_",i,"_FAMS_112gr2_cladeSTATS.txt",sep=""))
}
allFamRes<-do.call(rbind,results)
uniPGLS_fams<-read.table(paste(bbone,"_sample100_PGLSuniAll_FAMS_112gr2.txt",sep=""))

jpeg(file=paste("cladeLevel",bbone,"_scatter_FAMS_PGLSuni_Yaxis_PLOTTED_3part.jpg",sep=""),width=4,height=8, units="in", res=600, quality=100)
par(mfrow = c(3,1),oma = c(5,4,5,0) + 0.1, mar = c(3,1,1,1) + 0.1)
# type="n",
plot((richness)~MRCA,data=allFamRes, log="y",type="n",col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
axis()
#mtext(text="Clade crown age (Ma)", side=1, padj=3, font=2,cex=1)

plot((richness)~DR_harm,data=allFamRes, log="y",type="n",col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
#mtext(text="Clade ESDR mean (species / Ma)", side=1, padj=3, font=2,cex=1)

plot((richness)~DR_skew,data=allFamRes, log="y", type="n",col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
#mtext(text="Clade ESDR skewness", side=1, padj=3, font=2,cex=1)

title(main="", xlab = "",
      ylab = "log(richness)",
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)
dev.off()


jpeg(file=paste("cladeLevel",bbone,"_scatter_FAMS_PGLSuni_ALL_PLOTTED_3part_new.jpg",sep=""),width=4,height=8, units="in", res=600, quality=100)
par(mfrow = c(3,1),oma = c(5,4,5,0) + 0.1, mar = c(3,1,1,1) + 0.1)
ltys<-1
lwds<-4
colRaw<-col2rgb("light blue")/255
cols<-rgb(colRaw[1],colRaw[2],colRaw[3],alpha=0.5)
#plot(log(richness)~MRCA,data=dat60, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
smoothScatter(log(allFamRes$richness)~allFamRes$MRCA,xlab="",ylab="",nbin=50,nrpoints=0,colramp=colorRampPalette(c("white", "black")))
for (i in 1:length(uniPGLS_fams[,1])){
	abline(a=(uniPGLS_fams$i_1[i]),b=(uniPGLS_fams$MRCA[i]), untf=FALSE,col=cols, lwd=lwds,lty=ltys)
}
for (i in 1:length(uniPGLS_fams[,1])){
	if (uniPGLS_fams$p_1[i] > 0.05) {
	abline(a=(uniPGLS_fams$i_1[i]),b=(uniPGLS_fams$MRCA[i]), untf=FALSE,col=rgb(1,0,0,alpha=0.3), lwd=lwds,lty=ltys)		
	} else next
}
for (i in 1:length(uniPGLS_fams[,1])){
	abline(a=(uniPGLS_fams$i_1[i]),b=(uniPGLS_fams$MRCA[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}

#mtext(text="Clade crown age (Ma)", side=1, padj=3, font=2,cex=1)

#plot(log(richness)~DR_harm,data=dat60, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
smoothScatter(log(allFamRes$richness)~allFamRes$DR_harm,xlab="",ylab="",nbin=50,nrpoints=0,colramp=colorRampPalette(c("white", "black")))
for (i in 1:length(uniPGLS_fams[,1])){
	abline(a=(uniPGLS_fams$i_2[i]),b=(uniPGLS_fams$DR_harm[i]), untf=FALSE,col=cols, lwd=lwds,lty=ltys)
}
for (i in 1:length(uniPGLS_fams[,1])){
	if (uniPGLS_fams$p_2[i] > 0.05) {
	abline(a=(uniPGLS_fams$i_2[i]),b=(uniPGLS_fams$DR_harm[i]), untf=FALSE,col=rgb(1,0,0,alpha=0.3), lwd=lwds,lty=ltys)		
	} else next
}
for (i in 1:length(uniPGLS_fams[,1])){
	abline(a=(uniPGLS_fams$i_2[i]),b=(uniPGLS_fams$DR_harm[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}
#mtext(text="Clade ESDR mean (species / Ma)", side=1, padj=3, font=2,cex=1)

#plot(log(richness)~DR_skew,data=dat60, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
smoothScatter(log(allFamRes$richness)~allFamRes$DR_skew,xlab="",ylab="",nbin=50,nrpoints=0,colramp=colorRampPalette(c("white", "black")))
for (i in 1:length(uniPGLS_fams[,1])){
	abline(a=(uniPGLS_fams$i_3[i]),b=(uniPGLS_fams$DR_skew[i]), untf=FALSE,col=cols, lwd=lwds,lty=ltys)
}
for (i in 1:length(uniPGLS_fams[,1])){
	if (uniPGLS_fams$p_3[i] > 0.05) {
	abline(a=(uniPGLS_fams$i_3[i]),b=(uniPGLS_fams$DR_skew[i]), untf=FALSE,col=rgb(1,0,0,alpha=0.3), lwd=lwds,lty=ltys)		
	} else next
}
for (i in 1:length(uniPGLS_fams[,1])){
	abline(a=(uniPGLS_fams$i_3[i]),b=(uniPGLS_fams$DR_skew[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}
#mtext(text="Clade ESDR skewness", side=1, padj=3, font=2,cex=1)

title(main="", xlab = "",
      ylab = "",
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)
dev.off()







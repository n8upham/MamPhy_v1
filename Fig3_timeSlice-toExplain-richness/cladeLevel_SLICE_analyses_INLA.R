# Getting DATA to make the multi-level model for time slices
#####
#library(INLA); library(lme4)
library(moments)
library(nlme)
library(ape)
library(picante)
library(phytools)
library(geiger)
#library(RPANDA); library(phyloch); library(phylolm); library(paleotree); 
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")


setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/Diversification_analyses-Condamine/")
source("DR_functions.R")

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
#write.tree(mamPhy,file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""))
tree1=scan(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""), what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

#tree2<-read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""))
#==================
# ES and DR on per-tip basis
#==================

#######
# gives pairwise clade matrix from CAIC function
clade_matrix = readCAIC(tree1)

DR = 1/ES_v2(clade_matrix)
ES = ES_v2(clade_matrix)
res = cbind.data.frame(DR,ES)
res1 = res[order(rownames(res)),]

write.table(res1, file=paste(bbone,"_sample100_",i,"_DRtips.txt",sep=""))

#=======================================
# Make time slices to create the clades
#=======================================
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery

root=max(node.age(mamPhy)$ages)

allCladeSets<-vector("list",length=numSlices)
allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSets[[j]]<-treeSlice(mamPhy, slice=root-(sliceEvery*j), trivial=FALSE)
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

# record clade sizes per slice
lengths<-vector()
for (j in 1:length(allCladeSets)){
	lengths[j]<-length(allCladeSets[[j]])
}
names(lengths)<-allCladeSetNames 
lengths
write.table(lengths, file=paste(bbone,"_sample100_",i,"_timeslice_Lengths.txt",sep=""))

## write clades to per-slice files
#for (j in 1:length(allCladeSets)){
#	trees<-allCladeSets[[j]]
#	write.tree(trees,file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"_newick.trees",sep=""))
#	}
#
## load back in
#allCladeSets<-vector("list",length(allCladeSetNames))
#for (j in 1:length(allCladeSetNames)){
#	allCladeSets[[j]]<-read.tree(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"_newick.trees",sep=""))
#	}

#===================================
# Calculate per-SLICE, per-clade summary values 
#===================================

# get node times for tree
btimes<-branching.times(mamPhy)

# yule function
ymle = function(tree){ (.subset2(tree,3)-1L)/sum(.subset2(tree,2)) } # this take the # of number of nodes in a tree (minus 1) / sum of branch lengths.

# do per-slice, per-clade calcs
for(j in 1:length(allCladeSets)){
cladeSet<-allCladeSets[[j]]

	# empty data frames to fill
	DR_harm<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	DR_cv<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	DR_skew<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	DR_kurt<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	percentSamp<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	richness<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	MRCA<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	PB_Div<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD_Lam<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD_Mu<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD_Div<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD.ms0<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD.ms0p5<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD.ms0p9<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD.ms0_stem<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD.ms0p5_stem<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD.ms0p9_stem<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	
	for (k in 1:length(cladeSet)){
	cladeSp<-cladeSet[[k]]$tip.label
		x<-res1[match(cladeSp,rownames(res1)),"DR"]
		DR_harm[k,] <- 1/(mean(1/x))
		DR_cv[k,] <- (sd(x)/mean(x))*100
		DR_skew[k,] <- skewness(x)
		DR_kurt[k,] <- kurtosis(x)
		percentSamp[k,] <- length(which(cladesDR[match(cladeSp,cladesDR$tiplabel),"samp"]=="sampled"))/length(cladeSp)
		richness[k,] <- length(cladeSp)
		node <- getMRCA(mamPhy, cladeSp)
		MRCA[k,] <- btimes[node-5911] #taking the height of SAMPLED tree
	}
	if (length(cladeSp) > 2) {
	# Yule model
		PB_Div[k,]<-ymle(cladeSet[[k]])
		# BD model
		bd<-birthdeath(cladeSet[[k]])
		BD_Lam[k,]<-bd$para[[2]]/(1-bd$para[[1]])
		BD_Mu[k,]<-bd$para[[1]]*(bd$para[[2]]/(1-bd$para[[1]]))
		BD_Div[k,]<-bd$para[[2]]
		# BD Mag and Sand
	    BD.ms0_stem[k,]<-bd.ms(phy=cladeSet[[k]], missing=0, epsilon=0, crown=FALSE) # Assuming no extinction
    	BD.ms0p5_stem[k,]<-bd.ms(phy=cladeSet[[k]], missing=0, epsilon=0.5, crown=FALSE) # Assuming medium extinction 
     	BD.ms0p9_stem[k,]<-bd.ms(phy=cladeSet[[k]], missing=0, epsilon=0.9, crown=FALSE) # Assuming high extinction
		cladeSet[[k]]$root.edge<-0
	    BD.ms0[k,]<-bd.ms(phy=cladeSet[[k]], missing=0, epsilon=0, crown=TRUE) # Assuming no extinction
    	BD.ms0p5[k,]<-bd.ms(phy=cladeSet[[k]], missing=0, epsilon=0.5, crown=TRUE) # Assuming medium extinction 
     	BD.ms0p9[k,]<-bd.ms(phy=cladeSet[[k]], missing=0, epsilon=0.9, crown=TRUE) # Assuming high extinction
		} else NULL
	}
	res2<-cbind(DR_harm, DR_cv, DR_skew, DR_kurt, percentSamp, richness, MRCA, PB_Div, BD_Lam, BD_Mu, BD_Div, BD.ms0, BD.ms0p5, BD.ms0p9, BD.ms0_stem, BD.ms0p5_stem, BD.ms0p9_stem, i, j*5)
	colnames(res2)<-c("DR_harm","DR_cv", "DR_skew", "DR_kurt", "percentSamp", "richness", "MRCA", "PB_Div", "BD_Lam", "BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9", "BD.ms0_stem", "BD.ms0p5_stem", "BD.ms0p9_stem", "tree", "slice")
	
	rownames(res2)<-paste(i,"_",j,"_",c(1:length(res2[,1])),sep="")

	write.table(res2,paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""))
}

# create slice phys
slicePhys<-vector("list",length(allCladeSets))
for (k in 1:length(allCladeSets)){
cladeReps<-vector()
for (j in 1:length(allCladeSets[[k]])){
	cladeSp<-allCladeSets[[k]][[j]]$tip.label
	cladeReps[j]<-cladeSp[1]
	}
toDrop<-setdiff(mamPhy$tip.label,cladeReps)
slicePhys[[k]]<-drop.tip(mamPhy,toDrop)
slicePhys[[k]]$tip.label<-paste(i,"_",k,"_",c(1:length(slicePhys[[k]]$tip.label)),sep="")
}

# write slice phys
for(j in 1:length(slicePhys)){
	write.tree(slicePhys[[j]], file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_slicePhy-5to70Ma_by5.trees",sep=""), append=TRUE)
}
###########
#######

# do the PGLS
#######
# load back in tables
results<-vector("list",length(allCladeSetNames))
for (j in 1: length(allCladeSetNames)){
	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
	results[[j]]<-res
}

# UNIVARIATE first
sliceTimes<-seq(-5,-70,-5)
uniPGLS_allSlopes<-data.frame(matrix(NA, nrow = length(results), ncol = 13),row.names=sliceTimes)
colnames(uniPGLS_allSlopes)<-c("MRCA","DR_harm","DR_skew","PB_Div","BD_Lam","BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9", "BD.ms0_stem", "BD.ms0p5_stem", "BD.ms0p9_stem")

uniPGLS_allInts<-data.frame(matrix(NA, nrow = length(results), ncol = 13),row.names=sliceTimes)
colnames(uniPGLS_allInts)<-c(paste("i_",1:13,sep=""))

uniPGLS_allPs<-data.frame(matrix(NA, nrow = length(results), ncol = 13),row.names=sliceTimes)
colnames(uniPGLS_allPs)<-c(paste("p_",1:13,sep=""))

# MULTIVARIATE next
sliceTimes<-seq(-5,-70,-5)
partSlopesPGLS_All<-data.frame(matrix(NA, nrow = length(results), ncol = 8),row.names=sliceTimes)
colnames(partSlopesPGLS_All)<-c("int","MRCA","DR_harm","DR_skew","AIC","Pval1","Pval2","Pval3")

for (j in 1:length(results)){
	rownames(results[[j]])<-slicePhys[[j]]$tip.label
	cladeData<-treedata(slicePhys[[j]],na.omit(results[[j]]))
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
write.table(partSlopesPGLS_All,paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlice_5Ma-to-70Ma_partSlopes.txt",sep=""))

write.table(uniPGLS_allInts,paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniINTS.txt",sep=""))
write.table(uniPGLS_allSlopes,paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniSLOPES.txt",sep=""))
write.table(uniPGLS_allPs,paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniPs.txt",sep=""))

}


## ^ cool, that should be working then for FREQUENTIST approach. Boom.
###
# RE-do the PGLS as SCALED-- so that the EFFECT sizes are comparable.

# BASIC SETUP, run in PARALLEL.
library(moments)
library(nlme)
library(ape)
library(picante)
library(phytools)
library(geiger)

setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/Diversification_analyses-Condamine/")
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut")

library(foreach);library(doSNOW)
cl = makeCluster(15, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees=100

# which backbone?
bbone<- "NDexp" #"FBD" # 

sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery

allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

foreach(i=1:ntrees, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

# read back in slice phys
slicePhys<-read.tree(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_slicePhy-5to70Ma_by5.trees",sep=""))

# load back in tables
resultsSCALE<-vector("list",length(allCladeSetNames))
for (j in 1: length(allCladeSetNames)){
	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
	vars<-cbind(res[,7],res[,1],res[,3],res[,6],res[,8:17])
	varsScale<-scale(vars, center=TRUE,scale=TRUE)
	resultsSCALE[[j]]<-cbind(res$richness,varsScale)
	colnames(resultsSCALE[[j]])<-c("richness","MRCA","DR_harm","DR_skew","richnessZ",colnames(res[,8:17]))
}

resultsUN<-vector("list",length(allCladeSetNames))
for (j in 1: length(allCladeSetNames)){
	resultsUN[[j]]<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
}

# testing 2
# want to log-transform the 3 predictor variables and then repeat the analysis...
# comparing actual values vs scaled/transformed...

# load in tables and LOG trans them...
resultsLOG<-vector("list",length(allCladeSetNames))
for (j in 1: length(allCladeSetNames)){
	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
	vars<-cbind(res[,7],res[,1],res[,3],res[,6],res[,8:17])
	varsLog<-log(vars)
	resultsLOG[[j]]<-cbind(res$richness,varsLog)
	colnames(resultsLOG[[j]])<-c("richness","MRCA","DR_harm","DR_skew","richnessLog",colnames(res[,8:17]))
}

resultsRECIP<-vector("list",length(allCladeSetNames))
for (j in 1: length(allCladeSetNames)){
	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
	vars<-cbind(res[,7],res[,1],res[,3],res[,6],res[,8:17])
	vars2<-1/vars
	resultsRECIP[[j]]<-cbind(res$richness,vars2)
	colnames(resultsRECIP[[j]])<-c("richness","MRCA","DR_harm","DR_skew","richnessRECIP",colnames(res[,8:17]))
}

resultsCUBIC<-vector("list",length(allCladeSetNames))
for (j in 1: length(allCladeSetNames)){
	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
	vars<-cbind(res[,7],res[,1],res[,3],res[,6],res[,8:17])
	vars2<-vars^(1/3)
	resultsCUBIC[[j]]<-cbind(res$richness,vars2)
	colnames(resultsCUBIC[[j]])<-c("richness","MRCA","DR_harm","DR_skew","richnessCUBIC",colnames(res[,8:17]))
}

resultsSQROOT<-vector("list",length(allCladeSetNames))
for (j in 1: length(allCladeSetNames)){
	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
	vars<-cbind(res[,7],res[,1],res[,3],res[,6],res[,8:17])
	vars2<-vars^(1/2)
	resultsSQROOT[[j]]<-cbind(res$richness,vars2)
	colnames(resultsSQROOT[[j]])<-c("richness","MRCA","DR_harm","DR_skew","richnessSQROOT",colnames(res[,8:17]))
}

resultsSQ<-vector("list",length(allCladeSetNames))
for (j in 1: length(allCladeSetNames)){
	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
	vars<-cbind(res[,7],res[,1],res[,3],res[,6],res[,8:17])
	vars2<-vars^(2)
	resultsSQ[[j]]<-cbind(res$richness,vars2)
	colnames(resultsSQ[[j]])<-c("richness","MRCA","DR_harm","DR_skew","richnessSQ",colnames(res[,8:17]))
}

resultsUN<-vector("list",length(allCladeSetNames))
for (j in 1: length(allCladeSetNames)){
	resultsUN[[j]]<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
}


## 
# same but focusing BY VARIABLE-- with MRCA for each of SEVERAL transforms, etc...
nums<-c(14,10,6,2)
pdf(file="test_predictors_normality_KS-tests_density_5transforms.pdf",width=25,height=14, onefile=TRUE)

VARS<-c("MRCA","DR_harm","DR_skew")
for (i in 1:3){
par(mfrow = c(4,7),oma = c(3,3,3,3) + 0.1, mar = c(2,1,3,1) + 0.1)
for (j in nums){
	yy<-cbind.data.frame(resultsUN[[j]][,VARS[i]],as.vector(resultsSCALE[[j]][,VARS[i]]),as.vector(resultsLOG[[j]][,VARS[i]]),as.vector(resultsRECIP[[j]][,VARS[i]]),as.vector(resultsCUBIC[[j]][,VARS[i]]),as.vector(resultsSQROOT[[j]][,VARS[i]]),as.vector(resultsSQ[[j]][,VARS[i]]))
	colnames(yy)<-c("Raw","Z-score","Log", "Reciprocal","Cubic-root","Square-root","Square")
	var<-na.omit(yy[,"Raw"])
	plot(density(x=var),main=paste("slice ",j*5," Ma",sep=""),ylab="", xlab="")
	ks<-ks.test(x=var,y=pnorm,alternative="two.sided")
	mtext(paste("P = ",round(ks[[2]],digits=3),sep=""), cex=0.8)
	mtext(side=1, line=2, text="Raw", cex=0.8, font=2)

	var<-na.omit(yy[,"Z-score"])
	plot(density(x=var), main="",ylab="", xlab="")
	ks<-ks.test(x=var,y=pnorm,alternative="two.sided")
	mtext(paste("P = ",round(ks[[2]],digits=3),sep=""), cex=0.8)
	mtext(side=1, line=2, text="Z-score", cex=0.8, font=2)

	var<-na.omit(yy[,"Log"])
	plot(density(x=var),main="",ylab="", xlab="")
	ks<-ks.test(x=var,y=pnorm,alternative="two.sided")
	mtext(paste("P = ",round(ks[[2]],digits=3),sep=""), cex=0.8)
	mtext(side=1, line=2, text="Log", cex=0.8, font=2)

	var<-na.omit(yy[,"Reciprocal"])
	plot(density(x=var),main="",ylab="", xlab="")
	ks<-ks.test(x=var,y=pnorm,alternative="two.sided")
	mtext(paste("P = ",round(ks[[2]],digits=3),sep=""), cex=0.8)
	mtext(side=1, line=2, text="Reciprocal", cex=0.8, font=2)

	var<-na.omit(yy[,"Cubic-root"])
	plot(density(x=var),main="",ylab="", xlab="")
	ks<-ks.test(x=var,y=pnorm,alternative="two.sided")
	mtext(paste("P = ",round(ks[[2]],digits=3),sep=""), cex=0.8)
	mtext(side=1, line=2, text="Cubic-root", cex=0.8, font=2)

	var<-na.omit(yy[,"Square-root"])
	plot(density(x=var),main="",ylab="", xlab="")
	ks<-ks.test(x=var,y=pnorm,alternative="two.sided")
	mtext(paste("P = ",round(ks[[2]],digits=3),sep=""), cex=0.8)
	mtext(side=1, line=2, text="Square-root", cex=0.8, font=2)

	var<-na.omit(yy[,"Square"])
	plot(density(x=var),main="",ylab="", xlab="")
	ks<-ks.test(x=var,y=pnorm,alternative="two.sided")
	mtext(paste("P = ",round(ks[[2]],digits=3),sep=""), cex=0.8)
	mtext(side=1, line=2, text="Square", cex=0.8, font=2)
}
title(main=VARS[i],
      outer = TRUE, line = 1,cex.main=1.5,font.main=2,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)

}

dev.off()

## 
# same but BY VARIABLE-- with MRCA for each of 3 transforms, etc...
nums<-c(14,10,6,2)
pdf(file="test_predictors_normality_KS-tests_Raw-z-log_density.pdf",width=10,height=7, onefile=TRUE)

vars<-c("MRCA","DR_harm","DR_skew")
for (i in 1:3){
par(mfrow = c(4,3),oma = c(3,3,3,3) + 0.1, mar = c(2,1,3,1) + 0.1)
for (j in nums){
	yy<-cbind.data.frame(resultsUN[[j]][,vars[i]],resultsSCALE[[j]][,vars[i]],resultsLOG[[j]][,vars[i]])
	colnames(yy)<-c("Raw","Z-score","Log")
	var<-na.omit(yy[,"Raw"])
	plot(density(x=var),main=paste("slice ",j*5," Ma",sep=""),ylab="", xlab="")
	ks<-ks.test(x=var,y=pnorm,alternative="two.sided")
	mtext(paste("P = ",round(ks[[2]],digits=3),sep=""), cex=0.8)
	mtext(side=1, line=2, text="Raw", cex=0.8, font=2)

	var<-na.omit(yy[,"Z-score"])
	plot(density(x=var), main="",ylab="", xlab="")
	ks<-ks.test(x=var,y=pnorm,alternative="two.sided")
	mtext(paste("P = ",round(ks[[2]],digits=3),sep=""), cex=0.8)
	mtext(side=1, line=2, text="Z-score", cex=0.8, font=2)

	var<-na.omit(yy[,"Log"])
	plot(density(x=var),main="",ylab="", xlab="")
	ks<-ks.test(x=var,y=pnorm,alternative="two.sided")
	mtext(paste("P = ",round(ks[[2]],digits=3),sep=""), cex=0.8)
	mtext(side=1, line=2, text="Log", cex=0.8, font=2)

}
title(main=vars[i],
      outer = TRUE, line = 1,cex.main=1.5,font.main=2,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)

}

dev.off()

nums<-c(14,10,6,2)
pdf(file="test_predictors_normality_KS-tests_Raw-z-log_density.pdf",width=10,height=7, onefile=TRUE)

titles<-c("Raw predictors","Z-score predictors","Log predictors")
res<-list(resultsUN,resultsSCALE,resultsLOG)
for (i in 1:3){
results<-res[[i]]
par(mfrow = c(4,3),oma = c(3,3,3,3) + 0.1, mar = c(2,1,3,1) + 0.1)
for (j in nums){
	xx<-cbind.data.frame(results[[j]][,"MRCA"],results[[j]][,"DR_harm"],results[[j]][,"DR_skew"])
	colnames(xx)<-c("MRCA","DR_harm","DR_skew")
	var<-na.omit(xx[,"MRCA"])
	plot(density(x=var),main=paste("slice ",j*5," Ma",sep=""),ylab="", xlab="")
	ks<-ks.test(x=var,y=pnorm,alternative="two.sided")
	mtext(paste("P = ",round(ks[[2]],digits=3),sep=""), cex=0.8)
	mtext(side=1, line=2, text="MRCA", cex=0.8)

	var<-na.omit(xx[,"DR_harm"])
	plot(density(x=var), main="",ylab="", xlab="")
	ks<-ks.test(x=var,y=pnorm,alternative="two.sided")
	mtext(paste("P = ",round(ks[[2]],digits=3),sep=""), cex=0.8)
	mtext(side=1, line=2, text="DR_harm", cex=0.8)

	var<-na.omit(xx[,"DR_skew"])
	plot(density(x=var),main="",ylab="", xlab="")
	ks<-ks.test(x=var,y=pnorm,alternative="two.sided")
	mtext(paste("P = ",round(ks[[2]],digits=3),sep=""), cex=0.8)
	mtext(side=1, line=2, text="DR_skew", cex=0.8)

}
title(main=titles[i],
      outer = TRUE, line = 1,cex.main=1.5,font.main=2,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)

}

dev.off()




pdf(file="test_PGLS_scaling-vs-actual_MRCA_LOGtrans.pdf",width=8,height=(4*length(nums)))
pdf(file="test_PGLS_scaling-vs-actual_DRmean_LOGtrans.pdf",width=8,height=(4*length(nums)))
pdf(file="test_PGLS_scaling-vs-actual_DRskew_LOGtrans.pdf",width=8,height=(4*length(nums)))

par(mfrow=c(length(nums),2))

for (j in nums){
	xx<-cbind.data.frame(resultsUN[[j]][,"DR_skew"],resultsLOG[[j]][,"DR_skew"],resultsUN[[j]][,"richness"], row.names=rownames(resultsLOG[[j]]))
	colnames(xx)<-c("DR_skew","DR_skew_Log","N")

	if(j==2){
		xx<-xx[which(xx$DR_skew!=0),]
	}

	cladeData<-treedata(slicePhys[[j]],na.omit(xx))#results[[j]]))
	dat<-as.data.frame(cladeData$data)
	form1<-as.formula(log(N) ~ DR_skew)	
	fit1<-gls(form1, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
	sum1<-summary(fit1)
	form2<-as.formula(log(N) ~ DR_skew_Log)
	fit2<-gls(form2, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
	sum2<-summary(fit2)

plot(log(N)~DR_skew,data=xx, cex.axis=1.5,cex.lab=1.5)
abline(sum1$coef[[1]],sum1$coef[[2]])
mtext(paste("slope=",round(sum1$coef[[2]],3),"; P=",round(sum1$tTable[[8]],3),sep=""))
title(main=allCladeSetNames[[j]], cex.main=1.5, font=2)
plot(log(N)~DR_skew_Log,data=xx, cex.axis=1.5,cex.lab=1.5)
abline(sum2$coef[[1]],sum2$coef[[2]])
mtext(paste("slope=",round(sum2$coef[[2]],3),"; P=",round(sum2$tTable[[8]],3),sep=""))

}
dev.off()

setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/")
slicePhys<-read.tree(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_slicePhy-5to70Ma_by5.trees",sep=""))
# load back in tables
resultsUN<-vector("list",length(allCladeSetNames))
for (j in 1: length(allCladeSetNames)){
	resultsUN[[j]]<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
}
resultsSQ<-vector("list",length(allCladeSetNames))
for (j in 1: length(allCladeSetNames)){
	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
	vars<-cbind(res[,7],res[,1],res[,3],res[,6],res[,8:17])
	vars2<-vars^(2)
	resultsSQ[[j]]<-cbind(res$richness,vars2)
	colnames(resultsSQ[[j]])<-c("richness","MRCA_SQ","DR_harm","DR_skew","richnessSQ",colnames(res[,8:17]))
}
resultsSQ_SCALE<-vector("list",length(allCladeSetNames))
for (j in 1: length(allCladeSetNames)){
	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
	vars<-cbind(res[,7],res[,1],res[,3],res[,6],res[,8:17])
	vars2<-vars^(2)
	varsScale<-scale(vars2, center=TRUE,scale=TRUE)
	resultsSQ_SCALE[[j]]<-cbind(res$richness,varsScale)
	colnames(resultsSQ_SCALE[[j]])<-c("richness","MRCA_SQ_Scaled","DR_harm","DR_skew","richnessSQ",colnames(res[,8:17]))
}


# testing 1
nums<-c(14,10,6,2)

pdf(file="test_PGLS_scaling-vs-actual_MRCA-SQ.pdf",width=8,height=(4*length(nums)))
#pdf(file="test_PGLS_scaling-vs-actual_DRmean.pdf",width=8,height=(4*length(nums)))
#pdf(file="test_PGLS_scaling-vs-actual_DRskew.pdf",width=8,height=(4*length(nums)))
par(mfrow=c(length(nums),2))

for (j in nums){
xx<-cbind.data.frame(resultsSQ[[j]][,"MRCA_SQ"],resultsSQ_SCALE[[j]][,"MRCA_SQ_Scaled"],resultsSQ[[j]][,"richness"])
colnames(xx)<-c("MRCA_SQ","MRCA_SQ_Scaled","N")
#xx<-cbind(resultsUN[[j]][,"DR_harm"],results[[j]][,"DR_harm"],resultsUN[[j]][,"richness"])
#colnames(xx)<-c("DR_harm","DR_harmScaled","N")
#xx<-cbind(resultsUN[[j]][,"DR_skew"],results[[j]][,"DR_skew"],resultsUN[[j]][,"richness"])
#colnames(xx)<-c("DR_skew","DR_skewScaled","N")

cladeData<-treedata(slicePhys[[j]],na.omit(xx))#results[[j]]))
dat<-as.data.frame(cladeData$data)
form1<-as.formula(log(N) ~ MRCA_SQ)	
fit1<-gls(form1, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
sum1<-summary(fit1)
form2<-as.formula(log(N) ~ MRCA_SQ_Scaled)
fit2<-gls(form2, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
sum2<-summary(fit2)

plot(log(N)~MRCA_SQ,data=xx)
abline(sum1$coef[[1]],sum1$coef[[2]])
mtext(paste("slope=",round(sum1$coef[[2]],3),"; P=",round(sum1$tTable[[8]],3),sep=""))
title(allCladeSetNames[[j]])
plot(log(N)~MRCA_SQ_Scaled,data=xx)
abline(sum2$coef[[1]],sum2$coef[[2]])
mtext(paste("slope=",round(sum2$coef[[2]],3),"; P=",round(sum2$tTable[[8]],3),sep=""))

}

dev.off()


#### 
# RUNNING of the SCALED analyses...
# RE-do the PGLS as SCALED-- so that the EFFECT sizes are comparable.

# BASIC SETUP, run in PARALLEL.
library(moments)
library(nlme)
library(ape)
library(picante)
library(phytools)
library(geiger)

setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/Diversification_analyses-Condamine/")
#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut")

library(foreach);library(doSNOW)
cl = makeCluster(20, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees=100

# which backbone?
bbone<- "NDexp" #"FBD" # 

sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery

allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

foreach(i=1:ntrees, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

# read back in slice phys
slicePhys<-read.tree(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_slicePhy-5to70Ma_by5.trees",sep=""))

# load back in tables
# Z-score STANDARDIZED from RAW
# resultsSCALE<-vector("list",length(allCladeSetNames))
# for (j in 1: length(allCladeSetNames)){
# 	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
# 	vars<-cbind(res[,7],res[,1],res[,3],res[,6],res[,8:17])
# 	varsScale<-scale(vars, center=TRUE,scale=TRUE)
# 	resultsSCALE[[j]]<-cbind(res$richness,varsScale)
# 	colnames(resultsSCALE[[j]])<-c("richness","MRCA","DR_harm","DR_skew","richnessZ",colnames(res[,8:17]))
# }

# Z-score STANDARDIZED from SQ MRCA and RAW others
# results<-vector("list",length(allCladeSetNames))
# for (j in 1: length(allCladeSetNames)){
# 	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
# 	mrcaSQ<-(res[,"MRCA"])^2
# 	vars<-cbind(mrcaSQ,res[,7],res[,1],res[,3],res[,6],res[,8:17])
# 	varsScale<-scale(vars, center=TRUE,scale=TRUE)
# 	results[[j]]<-cbind(res$richness,varsScale)
# 	colnames(results[[j]])<-c("richness","mrcaSQ","MRCA","DR_harm","DR_skew","richnessZ",colnames(res[,8:17]))
# }
 
# RAW but with SQ MRCA also
results<-vector("list",length(allCladeSetNames))
for (j in 1: length(allCladeSetNames)){
	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
	mrcaSQ<-(res[,"MRCA"])^2
	vars<-cbind(mrcaSQ,res[,7],res[,1],res[,3],res[,6],res[,8:17])
	results[[j]]<-cbind(res$richness,vars)
	colnames(results[[j]])<-c("richness","mrcaSQ","MRCA","DR_harm","DR_skew","richnessZ",colnames(res[,8:17]))
}

# RE_SCALED -- MULTI- and UNI-VARIATE
sliceTimes<-seq(-5,-70,-5)
partSlopesPGLS_All<-data.frame(matrix(NA, nrow = length(results), ncol = 8),row.names=sliceTimes)
colnames(partSlopesPGLS_All)<-c("int","mrcaSQ","DR_harm","DR_skew","AIC","Pval1","Pval2","Pval3")
partSlopesPGLS_12<-data.frame(matrix(NA, nrow = length(results), ncol = 6),row.names=sliceTimes)
colnames(partSlopesPGLS_12)<-c("int","mrcaSQ","DR_harm","AIC","Pval1","Pval2")

uniPGLS_allSlopes<-data.frame(matrix(NA, nrow = length(results), ncol = 11),row.names=sliceTimes)
colnames(uniPGLS_allSlopes)<-c("mrcaSQ","MRCA","DR_harm","DR_skew","PB_Div","BD_Lam","BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9")
uniPGLS_allInts<-data.frame(matrix(NA, nrow = length(results), ncol = 11),row.names=sliceTimes)
colnames(uniPGLS_allInts)<-c(paste("i_",1:11,sep=""))
uniPGLS_allPs<-data.frame(matrix(NA, nrow = length(results), ncol = 11),row.names=sliceTimes)
colnames(uniPGLS_allPs)<-c(paste("p_",1:11,sep=""))

for (j in 1:length(results)){
	rownames(results[[j]])<-slicePhys[[j]]$tip.label
	cladeData<-treedata(slicePhys[[j]],na.omit(results[[j]]))
	dat<-as.data.frame(cladeData$data)

	form<-(log(richness) ~ mrcaSQ + DR_harm + DR_skew)
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

	form<-(log(richness) ~ mrcaSQ + DR_harm)
	fit2<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
	sum<-summary(fit2)

	partSlopesPGLS_12[j,1]<-round(sum$coef[[1]],digits=3)
	partSlopesPGLS_12[j,2]<-round(sum$coef[[2]],digits=3)
	partSlopesPGLS_12[j,3]<-round(sum$coef[[3]],digits=3)
	partSlopesPGLS_12[j,4]<-round(sum$AIC,digits=0)
	partSlopesPGLS_12[j,5]<-round(sum$tTable[11], digits=3)
	partSlopesPGLS_12[j,6]<-round(sum$tTable[12], digits=3)

	vars<-c("mrcaSQ","MRCA","DR_harm","DR_skew","PB_Div","BD_Lam","BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9")
	for(k in 1:length(uniPGLS_allSlopes)){
		form<-as.formula(paste("log(richness) ~ ", vars[k],sep=""))
		fit1<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
		sum1<-summary(fit1)
		
		uniPGLS_allInts[j,k]<-round(sum1$coef[[1]],digits=3)
		uniPGLS_allSlopes[j,k]<-round(sum1$coef[[2]],digits=3)
		uniPGLS_allPs[j,k]<-round(sum1$tTable[8], digits=3)
	}	

}
#write.table(partSlopesPGLS_All,paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlices_mrcaSQ_SCALED_partSlopes123.txt",sep=""))
#write.table(partSlopesPGLS_12,paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlices_mrcaSQ_SCALED_partSlopes12.txt",sep=""))
#
#uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs)
#write.table(uniPGLS,paste(bbone,"_sample100_",i,"_PGLS_timeSlices_mrcaSQ_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))

write.table(partSlopesPGLS_All,paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlices_mrcaSQ_RAW_partSlopes123.txt",sep=""))
write.table(partSlopesPGLS_12,paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlices_mrcaSQ_RAW_partSlopes12.txt",sep=""))

uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs)
write.table(uniPGLS,paste(bbone,"_sample100_",i,"_PGLS_timeSlices_mrcaSQ_RAW_uniINTS_uniSLO_uniP.txt",sep=""))

}



##
# With the RE-arrangement of variables too...
##
# RE_SCALED -- MULTI- and UNI-VARIATE
sliceTimes<-seq(-5,-70,-5)
partSlopesPGLS_All<-data.frame(matrix(NA, nrow = length(results), ncol = 8),row.names=sliceTimes)
colnames(partSlopesPGLS_All)<-c("int","MRCA","DR_harm","DR_skew","AIC","Pval1","Pval2","Pval3")

#partSlopesPGLS_213<-data.frame(matrix(NA, nrow = length(results), ncol = 8),row.names=sliceTimes)
#colnames(partSlopesPGLS_213)<-c("int","DR_harm","MRCA","DR_skew","AIC","Pval1","Pval2","Pval3")
#
#partSlopesPGLS_321<-data.frame(matrix(NA, nrow = length(results), ncol = 8),row.names=sliceTimes)
#colnames(partSlopesPGLS_321)<-c("int","DR_skew","DR_harm","MRCA","AIC","Pval1","Pval2","Pval3")

partSlopesPGLS_12BD<-data.frame(matrix(NA, nrow = length(results), ncol = 8),row.names=sliceTimes)
colnames(partSlopesPGLS_12BD)<-c("int","MRCA","DR_harm","BD_Div","AIC","Pval1","Pval2","Pval3")

uniPGLS_allSlopes<-data.frame(matrix(NA, nrow = length(results), ncol = 13),row.names=sliceTimes)
colnames(uniPGLS_allSlopes)<-c("MRCA","DR_harm","DR_skew","PB_Div","BD_Lam","BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9", "BD.ms0_stem", "BD.ms0p5_stem", "BD.ms0p9_stem")
uniPGLS_allInts<-data.frame(matrix(NA, nrow = length(results), ncol = 13),row.names=sliceTimes)
colnames(uniPGLS_allInts)<-c(paste("i_",1:13,sep=""))
uniPGLS_allPs<-data.frame(matrix(NA, nrow = length(results), ncol = 13),row.names=sliceTimes)
colnames(uniPGLS_allPs)<-c(paste("p_",1:13,sep=""))

for (j in 1:length(results)){
	rownames(results[[j]])<-slicePhys[[j]]$tip.label
	cladeData<-treedata(slicePhys[[j]],na.omit(results[[j]]))
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

	form<-(log(richness) ~ DR_harm + MRCA + DR_skew)
	fit2<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
	sum<-summary(fit2)

	partSlopesPGLS_213[j,1]<-round(sum$coef[[1]],digits=3)
	partSlopesPGLS_213[j,2]<-round(sum$coef[[2]],digits=3)
	partSlopesPGLS_213[j,3]<-round(sum$coef[[3]],digits=3)
	partSlopesPGLS_213[j,4]<-round(sum$coef[[4]],digits=3)
	partSlopesPGLS_213[j,5]<-round(sum$AIC,digits=0)
	partSlopesPGLS_213[j,6]<-round(sum$tTable[14], digits=3)
	partSlopesPGLS_213[j,7]<-round(sum$tTable[15], digits=3)
	partSlopesPGLS_213[j,8]<-round(sum$tTable[16], digits=3)

	form<-(log(richness) ~ DR_skew + DR_harm + MRCA)
	fit3<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
	sum<-summary(fit3)

	partSlopesPGLS_321[j,1]<-round(sum$coef[[1]],digits=3)
	partSlopesPGLS_321[j,2]<-round(sum$coef[[2]],digits=3)
	partSlopesPGLS_321[j,3]<-round(sum$coef[[3]],digits=3)
	partSlopesPGLS_321[j,4]<-round(sum$coef[[4]],digits=3)
	partSlopesPGLS_321[j,5]<-round(sum$AIC,digits=0)
	partSlopesPGLS_321[j,6]<-round(sum$tTable[14], digits=3)
	partSlopesPGLS_321[j,7]<-round(sum$tTable[15], digits=3)
	partSlopesPGLS_321[j,8]<-round(sum$tTable[16], digits=3)

	form<-(log(richness) ~ MRCA + DR_harm + BD_Div)
	fit4<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
	sum<-summary(fit4)

	partSlopesPGLS_12BD[j,1]<-round(sum$coef[[1]],digits=3)
	partSlopesPGLS_12BD[j,2]<-round(sum$coef[[2]],digits=3)
	partSlopesPGLS_12BD[j,3]<-round(sum$coef[[3]],digits=3)
	partSlopesPGLS_12BD[j,4]<-round(sum$coef[[4]],digits=3)
	partSlopesPGLS_12BD[j,5]<-round(sum$AIC,digits=0)
	partSlopesPGLS_12BD[j,6]<-round(sum$tTable[14], digits=3)
	partSlopesPGLS_12BD[j,7]<-round(sum$tTable[15], digits=3)
	partSlopesPGLS_12BD[j,8]<-round(sum$tTable[16], digits=3)

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
write.table(partSlopesPGLS_All,paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlices_SCALED_partSlopes123.txt",sep=""))
write.table(partSlopesPGLS_213,paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlices_SCALED_partSlopes213.txt",sep=""))
write.table(partSlopesPGLS_321,paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlices_SCALED_partSlopes321.txt",sep=""))
write.table(partSlopesPGLS_12BD,paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlices_SCALED_partSlopes12BD.txt",sep=""))

uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs)
write.table(uniPGLS,paste(bbone,"_sample100_",i,"_PGLS_timeSlices_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))

}

###
# Do with including ** percentSamp ** as a major variable too...
##
# written in as:
##
res2<-cbind(DR_harm, DR_cv, DR_skew, DR_kurt, percentSamp, richness, MRCA, PB_Div, BD_Lam, BD_Mu, BD_Div, BD.ms0, BD.ms0p5, BD.ms0p9, BD.ms0_stem, BD.ms0p5_stem, BD.ms0p9_stem, i, j*5)
colnames(res2)<-c("DR_harm","DR_cv", "DR_skew", "DR_kurt", "percentSamp", "richness", "MRCA", "PB_Div", "BD_Lam", "BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9", "BD.ms0_stem", "BD.ms0p5_stem", "BD.ms0p9_stem", "tree", "slice")
rownames(res2)<-paste(i,"_",j,"_",c(1:length(res2[,1])),sep="")
write.table(res2,paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""))

# so now redo the PGLS with that var.
# in PARALLEL.
#### 
# BASIC SETUP, run in PARALLEL.
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/Diversification_analyses-Condamine/")
#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut")
library(moments); library(nlme); library(ape); library(picante); library(phytools); library(geiger)

library(foreach);library(doSNOW)
cl = makeCluster(15, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees=100

# which backbone?
bbone<- "NDexp" #"FBD" # 

sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery

allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

all<-c(1:100)
finished<-c(3, 5, 7, 8, 11, 13, 17, 18, 19, 22, 23, 24, 25, 27, 29, 30, 33, 34, 35, 36, 37, 40, 41, 42, 45, 46, 48, 51, 52, 53, 54, 60, 62, 63, 65, 66, 67, 68, 70, 71, 74, 77, 78, 81, 82, 83, 88, 90, 91, 92, 93, 94, 95, 96, 97, 98, 100)
finished2<-c(1, 4, 5, 6, 7, 10, 12, 16, 17, 18, 20, 21, 22, 23, 26, 28, 29, 30, 32, 33, 34, 36, 37, 38, 39, 40)
finished3<-c(2, 7, 8, 9, 10, 13, 15, 16, 17, 19, 20, 21, 23, 24, 25, 26, 27)
#finished4<-c(4, 5, 6, 7, 10, 12, 13, 14, 16, 17, 18, 20, 21, 22, 23, 24)
#allFin<-c(finished,finished2, finished3)#,finished4)

allFinTotal<-c(1, 1, 1, 3, 3, 3, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 16, 16, 16, 17, 17, 17, 18, 18, 18, 19, 19, 19, 21, 21, 21, 22, 22, 22, 23, 23, 23, 24, 24, 24, 25, 25, 25, 27, 27, 27, 29, 29, 29, 30, 30, 30, 32, 32, 32, 33, 33, 33, 34, 34, 34, 35, 35, 35, 36, 36, 36, 37, 37, 37, 38, 38, 38, 39, 39, 39, 40, 40, 40, 41, 41, 41, 42, 42, 42, 44, 44, 44, 45, 45, 45, 46, 46, 46, 47, 47, 47, 48, 48, 48, 49, 49, 49, 50, 50, 50, 51, 51, 51, 52, 52, 52, 53, 53, 53, 54, 54, 54, 57, 57, 57, 59, 59, 59, 60, 60, 60, 61, 61, 61, 62, 62, 62, 63, 63, 63, 64, 64, 64, 65, 65, 65, 66, 66, 66, 67, 67, 67, 68, 68, 68, 70, 70, 70, 71, 71, 71, 72, 72, 72, 73, 73, 73, 74, 74, 74, 75, 75, 75, 77, 77, 77, 78, 78, 78, 79, 79, 79, 80, 80, 80, 81, 81, 81, 82, 82, 82, 83, 83, 83, 84, 84, 84, 85, 85, 85, 86, 86, 86, 88, 88, 88, 90, 90, 90, 91, 91, 91, 92, 92, 92, 93, 93, 93, 94, 94, 94, 95, 95, 95, 96, 96, 96, 97, 97, 97, 98, 98, 98, 100, 100, 100)
allFin<-unique(allFinTotal)
remaining<-setdiff(all,allFin)

#foreach(i=1:ntrees, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

foreach(i=remaining, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

# read back in slice phys
slicePhys<-read.tree(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_slicePhy-5to70Ma_by5.trees",sep=""))

# multiply edges by 100 for corPagel
modTrees<-vector("list",length(slicePhys))
for (j in 1:length(slicePhys)){
	modTree<-slicePhys[[j]]
	modTree$edge.length<-modTree$edge.length*100
	modTrees[[j]]<-modTree
}

# Z-score STANDARDIZED from RAW
results<-vector("list",length(allCladeSetNames))
for (j in 1: length(allCladeSetNames)){
	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
	vars<-res[,1:14]
	varsScale<-scale(vars, center=TRUE,scale=TRUE)
	results[[j]]<-cbind(res$richness,varsScale)
	colnames(results[[j]])<-c("richness","DR_harm","DR_cv","DR_skew","DR_kurt","percentSamp","richnessZ",colnames(res[,7:14]))
}

# ACTUAL data...
#results<-vector("list",length(allCladeSetNames))
#for (j in 1: length(allCladeSetNames)){
#	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
#	results[[j]]<-res[,1:14]
#}

# MULTI- and UNI-VARIATE
sliceTimes<-seq(-5,-70,-5)
partSlopesPGLS_All<-data.frame(matrix(NA, nrow = length(results), ncol = 16),row.names=sliceTimes)
colnames(partSlopesPGLS_All)<-c("int","MRCA","DR_harm","DR_skew","SE1","SE2","SE3","SE4","Lam_low","Lam_mean","Lam_up","AIC","Pval1","Pval2","Pval3","Pval4")
partSlopesPGLS_All_Per<-data.frame(matrix(NA, nrow = length(results), ncol = 17),row.names=sliceTimes)
colnames(partSlopesPGLS_All_Per)<-c("int","MRCA","DR_harm","DR_skew","percentSamp","SE1","SE2","SE3","SE4","SE5","Lam","AIC","Pval1","Pval2","Pval3","Pval4","Pval5")

uniPGLS_allSlopes<-data.frame(matrix(NA, nrow = length(results), ncol = 11),row.names=sliceTimes)
colnames(uniPGLS_allSlopes)<-c("MRCA","DR_harm","DR_skew","percentSamp","PB_Div","BD_Lam","BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9")
uniPGLS_allSEs<-data.frame(matrix(NA, nrow = length(results), ncol = 11),row.names=sliceTimes)
colnames(uniPGLS_allSEs)<-c(paste("SE_",1:11,sep=""))
uniPGLS_allInts<-data.frame(matrix(NA, nrow = length(results), ncol = 11),row.names=sliceTimes)
colnames(uniPGLS_allInts)<-c(paste("i_",1:11,sep=""))
uniPGLS_allPs<-data.frame(matrix(NA, nrow = length(results), ncol = 11),row.names=sliceTimes)
colnames(uniPGLS_allPs)<-c(paste("p_",1:11,sep=""))
uniPGLS_allLams<-data.frame(matrix(NA, nrow = length(results), ncol = 11),row.names=sliceTimes)
colnames(uniPGLS_allLams)<-c(paste("lam_",1:11,sep=""))

for (j in 1:length(results)){
	#rownames(results[[j]])<-slicePhys[[j]]$tip.label
	#cladeData<-treedata(slicePhys[[j]],na.omit(results[[j]]))
	cladeData<-treedata(modTrees[[j]],na.omit(results[[j]]))
	dat<-as.data.frame(cladeData$data)
	#dat<-as.data.frame(na.omit(results[[j]]))

	form<-(log(richness) ~ MRCA + DR_harm + DR_skew)
#	fit<-gls(form, data=dat, method="ML")
#	fit<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
		for (p in seq(0,1,by=0.01)) {possibleError <- tryCatch(
		      gls(form, correlation=corPagel(value=p,phy=cladeData$phy), data=dat, method="ML"),
		      error=function(e) e)
		if(inherits(possibleError, "gls")) break		
		if(inherits(possibleError, "error")) next}
		fit<-possibleError
	sum<-summary(fit)
	#int95s<-intervals(fit,which="var-cov")

	for(k in 1:8){
	partSlopesPGLS_All[j,k]<-round(sum$tTable[k],digits=3)
	}
	#partSlopesPGLS_All[j,9]<-round(int95s[[1]][1],digits=3)
	partSlopesPGLS_All[j,10]<-round(sum$modelStruct[[1]][1],digits=3)#int95s[[1]][2],digits=3)
	#partSlopesPGLS_All[j,11]<-round(int95s[[1]][3],digits=3)
	partSlopesPGLS_All[j,12]<-round(sum$AIC,digits=0)
	partSlopesPGLS_All[j,13]<-round(sum$tTable[13], digits=3)
	partSlopesPGLS_All[j,14]<-round(sum$tTable[14], digits=3)
	partSlopesPGLS_All[j,15]<-round(sum$tTable[15], digits=3)
	partSlopesPGLS_All[j,16]<-round(sum$tTable[16], digits=3)

	form2<-(log(richness) ~ MRCA + DR_harm + DR_skew + percentSamp)
#	fit2<-gls(form2, data=dat, method="ML")
#	fit2<-gls(form2, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
		for (p in seq(0,1,by=0.01)) {possibleError <- tryCatch(
		      fit2<-gls(form2, correlation=corPagel(value=p,phy=cladeData$phy), data=dat, method="ML"),
		      error=function(e) e)
		if(inherits(possibleError, "gls")) break		
		if(inherits(possibleError, "error")) next}
		fit2<-possibleError
	sum2<-summary(fit2)

	for(k in 1:10){
	partSlopesPGLS_All_Per[j,k]<-round(sum2$tTable[k],digits=3)
	}
	partSlopesPGLS_All_Per[j,11]<-round(sum2$modelStruct[[1]][1],digits=3)
	partSlopesPGLS_All_Per[j,12]<-round(sum2$AIC,digits=0)
	partSlopesPGLS_All_Per[j,13]<-round(sum2$tTable[16], digits=3)
	partSlopesPGLS_All_Per[j,14]<-round(sum2$tTable[17], digits=3)
	partSlopesPGLS_All_Per[j,15]<-round(sum2$tTable[18], digits=3)
	partSlopesPGLS_All_Per[j,16]<-round(sum2$tTable[19], digits=3)
	partSlopesPGLS_All_Per[j,17]<-round(sum2$tTable[20], digits=3)

	vars<-c("MRCA","DR_harm","DR_skew","percentSamp","PB_Div","BD_Lam","BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9")
	for(k in 1:length(uniPGLS_allSlopes)){
		form<-as.formula(paste("log(richness) ~ ", vars[k],sep=""))
#		fit1<-gls(form, data=dat, method="ML")
#		fit1<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
		for (p in seq(0,1,by=0.01)) {possibleError <- tryCatch(
		      gls(form, correlation=corPagel(value=p,phy=cladeData$phy), data=dat, method="ML"),
		      error=function(e) e)
		if(inherits(possibleError, "gls")) break		
		if(inherits(possibleError, "error")) next}
		fit1<-possibleError

		sum1<-summary(fit1)
		
		uniPGLS_allInts[j,k]<-round(sum1$tTable[1], digits=3)
		uniPGLS_allSlopes[j,k]<-round(sum1$tTable[2], digits=3)
		uniPGLS_allSEs[j,k]<-round(sum1$tTable[4], digits=3)
		uniPGLS_allPs[j,k]<-round(sum1$tTable[8], digits=3)
		uniPGLS_allLams[j,k]<-round(sum1$modelStruct[[1]][[1]], digits=3)
	}	
}
#corr<-"BROWNIAN"
corr<-"PAGEL"

#write.table(partSlopesPGLS_All,paste(bbone,"_sample100_",i,"_PGLSmulti_",corr,"_timeSlices_ACTUAL_partSlopes123.txt",sep=""))
#write.table(partSlopesPGLS_All_Per,paste(bbone,"_sample100_",i,"_PGLSmulti_",corr,"_timeSlices_wPercentSamp_ACTUAL_partSlopes123.txt",sep=""))
#uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs,uniPGLS_allSEs,uniPGLS_allLams)
#write.table(uniPGLS,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_timeSlices_wPercentSamp_ACTUAL_uniINTS_uniSLO_uniP.txt",sep=""))

write.table(partSlopesPGLS_All,paste(bbone,"_sample100_",i,"_PGLSmulti_",corr,"_noLamCI_timeSlices_SCALED_partSlopes123.txt",sep=""))
write.table(partSlopesPGLS_All_Per,paste(bbone,"_sample100_",i,"_PGLSmulti_",corr,"_noLamCI_timeSlices_wPercentSamp_SCALED_partSlopes123.txt",sep=""))
uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs,uniPGLS_allSEs,uniPGLS_allLams)
write.table(uniPGLS,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_noLamCI_timeSlices_wPercentSamp_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))

#write.table(partSlopesPGLS_All,paste(bbone,"_sample100_",i,"_PGLSmulti_noTREE_timeSlices_SCALED_partSlopes123.txt",sep=""))
#write.table(partSlopesPGLS_All_Per,paste(bbone,"_sample100_",i,"_PGLSmulti_noTREE_timeSlices_wPercentSamp_SCALED_partSlopes123.txt",sep=""))
#uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs,uniPGLS_allSEs,uniPGLS_allLams)
#write.table(uniPGLS,paste(bbone,"_sample100_",i,"_PGLS_noTREE_timeSlices_wPercentSamp_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))

#write.table(partSlopesPGLS_All,paste(bbone,"_sample100_",i,"_PGLSmulti_noTREE_timeSlices_ACTUAL_partSlopes123.txt",sep=""))
#write.table(partSlopesPGLS_All_Per,paste(bbone,"_sample100_",i,"_PGLSmulti_noTREE_timeSlices_wPercentSamp_ACTUAL_partSlopes123.txt",sep=""))
#uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs,uniPGLS_allSEs,uniPGLS_allLams)
#write.table(uniPGLS,paste(bbone,"_sample100_",i,"_PGLS_noTREE_timeSlices_wPercentSamp_ACTUAL_uniINTS_uniSLO_uniP.txt",sep=""))
}





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





##########
###############
#===================================
# SUMMARIZE
# Load data back in and COMBINE all slice clades per-tree, then across all trees
# Post doing PGLS
#===================================
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut")

setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/Diversification_analyses-Condamine/")
bbone<- "NDexp" #"FBD" # 
library(nlme); library(MASS); library(AICcmodavg); library(ape); library(phytools)

numTrees<-100
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery
sliceTimes<-seq(-5,-70,-5)

allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

# Make empty vectors and frames
allSlopes_per_i<-vector("list",length=numTrees)
slice<-c(1:14)
# load in data as is. 
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	#dat<-read.table(paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlices_SCALED_partSlopes123.txt",sep=""))
	#dat<-read.table(paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlice_5Ma-to-70Ma_partSlopes.txt",sep=""))
	dat<-read.table(paste(bbone,"_sample100_",i,"_PGLSmulti_BROWNIAN_timeSlices_SCALED_partSlopes123.txt",sep=""))
	allSlopes_per_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
allDat_MULTI_scaled<-do.call(rbind,allSlopes_per_i)

allSlopes_per_i<-vector("list",length=numTrees)
slice<-c(1:14)
# load in data as is. 
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	dat<-read.table(paste(bbone,"_sample100_",i,"_PGLSmulti_BROWNIAN_timeSlices_wPercentSamp_SCALED_partSlopes123.txt",sep=""))
	allSlopes_per_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
allDat_MULTI_scaled_4var<-do.call(rbind,allSlopes_per_i)

allDat_UNI_i<-vector("list",length=numTrees)
# load in data as is. 
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	dat<-read.table(paste(bbone,"_sample100_",i,"_PGLS_BROWNIAN_timeSlices_wPercentSamp_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))
	#ints<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniINTS.txt",sep=""))
	#slopes<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniSLOPES.txt",sep=""))
	#Ps<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniPs.txt",sep=""))
	#dat<-cbind(ints,slopes,Ps)
	allDat_UNI_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
allDat_UNI<-do.call(rbind,allDat_UNI_i)


mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 
root=max(nodeHeights(mamPhy)[,2])

# plot through time
#pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_timeSlices_SCALED_partSlopes_PLOTTED_4part.pdf",sep=""),onefile=TRUE, width=4,height=9)
pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_timeSlices_UNscaled-ACTUAL_partSlopes_PLOTTED_4part.pdf",sep=""),onefile=TRUE, width=4,height=9)
dat<-allDat
p1<- dat$Pval1
p2<- dat$Pval2
p3<- dat$Pval3
adjLab<-0.08
Header<-"MULTIVARIATE - PGLS per time slice per tree"

#pdf(file=paste("cladeLevel",bbone,"_PGLS-UNI_timeSlices_SCALED_uniSlopes_PLOTTED_4part.pdf",sep=""),onefile=TRUE, width=4,height=9)
# pdf(file=paste("cladeLevel",bbone,"_PGLS-UNI_timeSlices_UNscaled-ACTUAL_uniSlopes_PLOTTED_4part.pdf",sep=""),onefile=TRUE, width=4,height=9)
# dat<-allDat_UNI #allDat
# p1<- dat$p_1# dat$Pval1
# p2<- dat$p_2# dat$Pval2
# p3<- dat$p_3# dat$Pval3
# adjLab<-0.08 # 0.02
# Header<-"UNIVARIATE - PGLS per time slice per tree" # "MULTIVARIATE - PGLS per time slice per tree"
#quartz(width=10,height=4)
#layout(matrix(c(1:3), 3, 1, byrow = TRUE))
par(mfrow = c(4,1), oma = c(5,4,5,0) + 0.1, mar = c(5,1,1,1) + 0.1)
#par(op)

# for NDexp MCC
plot(mamPhy, edge.width=0.3, cex=0.1, show.tip.label=FALSE, no.margin=TRUE, x.lim=c(110.8,175.9)) #, tip.color = tcol )
#mtext(text="(a) MCC tree", adj=0, font=2)#at=c(-72,0.21)
#axisPhylo(cex.axis=1, pos=-70, mgp=c(0,0.2,0))#, at=seq(5,70,5),side=1,labels=seq(5,70,5))
for(i in 1:numSlices){
	abline(v=root-(5*i),lty=2, lwd=1.5, col="grey")
}

cols<-rep(grey(0.3,alpha=0.5),length(dat[,1]))
cols[which(p1>0.05)]<-rgb(1,0,0,alpha=0.5)
pchs<-rep(1,length(dat[,1]))
pchs[which(p1>0.05)]<-2

#yLims<-c(0,2.5)
#yAxis<-c(0,0.5,1,1.5,NA)
yLims<-NULL
yAxis<-NULL

plot(MRCA ~ sliceTimes, data=dat, col=cols, pch=pchs, ylim=yLims, ylab="",xlab="", yaxt="n",xaxt="n")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
axis(side=2,at=yAxis,labels=TRUE)
mtext(text="(a) Crown age", padj=2, adj=adjLab, font=2)#at=c(-72,0.21)

	#fit1<-lm(MRCA ~ sliceTimes, data=allDat)
	#fit2<-lm(MRCA ~ sliceTimes + I(sliceTimes^2), data=allDat)
	#fit3<-lm(MRCA ~ sliceTimes + I(sliceTimes^2) + I(sliceTimes^3), data=allDat)
	#fits<-list(fit1,fit2,fit3)
	#AICs<-sapply(fits,AICc)

	m1<-lme(fixed = MRCA ~ sliceTimes, random = ~ 1 | tree/slice, data=dat)
	sum1<-summary(m1)
	m2<-lme(fixed = MRCA ~ sliceTimes + I(sliceTimes^2), random = ~ 1 | tree/slice, data=dat)
	sum2<-summary(m2)
	m3<-lme(fixed = MRCA ~ sliceTimes + I(sliceTimes^2) + I(sliceTimes^3), random = ~ 1 | tree/slice, data=dat)
	sum3<-summary(m3)
	models<-list(m1,m2,m3)
	AICs<-sapply(models,AICc)

	m<-models[which(min(AICs)==AICs)]
	sum<-summary(m[[1]])
	fixed<-sum$coef$fixed

	if (length(fixed)==4){
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2+ fixed[[4]]*x^3, add=TRUE, lwd=2, col="red")
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2+ fixed[[4]]*x^3, add=TRUE, lwd=1, lty=2)
	}
	
	if (length(fixed)==3){
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2, add=TRUE, lwd=2, col="red")
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2, add=TRUE, lwd=1, lty=2)
	}

	if (length(fixed)==2){
	curve(fixed[[1]] + fixed[[2]]*x, add=TRUE, lwd=2, col="red")
	curve(fixed[[1]] + fixed[[2]]*x, add=TRUE, lwd=1, lty=2)
	}

for(i in 1:numSlices){
	abline(v=(-5*i),lty=2, lwd=1.5, col=grey(0.6, alpha=0.5))
}
	abline(h=1.0,lty=1, lwd=1, col=grey(0.6, alpha=0.5))

cols<-rep(grey(0.3,alpha=0.5),length(dat[,1]))
cols[which(p2>0.05)]<-rgb(1,0,0,alpha=0.5)
pchs<-rep(1,length(dat[,1]))
pchs[which(p2>0.05)]<-2

plot(DR_harm ~ sliceTimes, data=dat, col=cols, pch=pchs, ylim=yLims, ylab="",xlab="", yaxt="n",xaxt="n")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
axis(side=2,at=yAxis,labels=TRUE)
mtext(text="(b) Tip DR mean", adj=adjLab, padj=2, font=2)#at=c(-72,0.21)

	m1<-lme(fixed = DR_harm ~ sliceTimes, random = ~ 1 | tree/slice, data=dat)
	sum1<-summary(m1)
	m2<-lme(fixed = DR_harm ~ sliceTimes + I(sliceTimes^2), random = ~ 1 | tree/slice, data=dat)
	sum2<-summary(m2)
	m3<-lme(fixed = DR_harm ~ sliceTimes + I(sliceTimes^2) + I(sliceTimes^3), random = ~ 1 | tree/slice, data=dat)
	sum3<-summary(m3)
	models<-list(m1,m2,m3)
	AICs<-sapply(models,AICc)

	m<-models[which(min(AICs)==AICs)]
	sum<-summary(m[[1]])
	fixed<-sum$coef$fixed

	if (length(fixed)==4){
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2+ fixed[[4]]*x^3, add=TRUE, lwd=2, col="red")
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2+ fixed[[4]]*x^3, add=TRUE, lwd=1, lty=2)
	}
	
	if (length(fixed)==3){
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2, add=TRUE, lwd=2, col="red")
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2, add=TRUE, lwd=1, lty=2)
	}

	if (length(fixed)==2){
	curve(fixed[[1]] + fixed[[2]]*x, add=TRUE, lwd=2, col="red")
	curve(fixed[[1]] + fixed[[2]]*x, add=TRUE, lwd=1, lty=2)
	}

for(i in 1:numSlices){
	abline(v=(-5*i),lty=2, lwd=1.5, col=grey(0.6, alpha=0.5))
}
	abline(h=1.0,lty=1, lwd=1, col=grey(0.6, alpha=0.5))

cols<-rep(grey(0.3,alpha=0.5),length(dat[,1]))
cols[which(p3>0.05)]<-rgb(1,0,0,alpha=0.5)
pchs<-rep(1,length(dat[,1]))
pchs[which(p3>0.05)]<-2

plot(DR_skew ~ sliceTimes, data=dat, col=cols, pch=pchs, ylim=yLims, yaxt="n",ylab="",xlab="")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
axis(side=2,at=yAxis,labels=TRUE)
mtext(text="(c) Tip DR skewness", adj=adjLab, padj=2, font=2)#at=c(-72,0.21)

	m1<-lme(fixed = DR_skew ~ sliceTimes, random = ~ 1 | tree/slice, data=dat)
	sum1<-summary(m1)
	m2<-lme(fixed = DR_skew ~ sliceTimes + I(sliceTimes^2), random = ~ 1 | tree/slice, data=dat)
	sum2<-summary(m2)
	m3<-lme(fixed = DR_skew ~ sliceTimes + I(sliceTimes^2) + I(sliceTimes^3), random = ~ 1 | tree/slice, data=dat)
	sum3<-summary(m3)
	models<-list(m1,m2,m3)
	AICs<-sapply(models,AICc)

	m<-models[which(min(AICs)==AICs)]
	sum<-summary(m[[1]])
	fixed<-sum$coef$fixed

	if (length(fixed)==4){
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2+ fixed[[4]]*x^3, add=TRUE, lwd=2, col="red")
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2+ fixed[[4]]*x^3, add=TRUE, lwd=1, lty=2)
	}
	
	if (length(fixed)==3){
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2, add=TRUE, lwd=2, col="red")
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2, add=TRUE, lwd=1, lty=2)
	}

	if (length(fixed)==2){
	curve(fixed[[1]] + fixed[[2]]*x, add=TRUE, lwd=2, col="red")
	curve(fixed[[1]] + fixed[[2]]*x, add=TRUE, lwd=1, lty=2)
	}

for(i in 1:numSlices){
	abline(v=(-5*i),lty=2, lwd=1.5, col=grey(0.6, alpha=0.5))
}
	abline(h=1.0,lty=1, lwd=1, col=grey(0.6, alpha=0.5))

title(main=Header, xlab = "Time slices before present (Ma)",
      #ylab = "Standardized effect on log (clade richness)",
      ylab = "Partial residual on log (clade richness)",
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)

dev.off()



####
# Then ALSO plot each separately...
load(file=paste(bbone,"_sample100__ALL100trees_slice5Ma-to-70Ma_cladeSTATS.Rdata",sep=""))

dat10<-allRes_100trees[which(allRes_100trees$slice==10),] # 65598 obs. of  19 variables
dat60<-allRes_100trees[which(allRes_100trees$slice==60),] # 2888 obs. of  19 variables:


###
# Compare all the UNI effect sizes of the variables. 
# boxplot I think...

allDat_UNI
uniDat10<-allDat_UNI[which(allDat_UNI$sliceTimes==-10),]
uniDat60<-allDat_UNI[which(allDat_UNI$sliceTimes==-60),]

uniPGLS_fams_SCALED<-read.table(paste(bbone,"_sample100_PGLSuniAll_FAMS_112gr2_SCALED.txt",sep=""))


# 10 Ma
dd<-sapply(uniDat10[,15:24],mean)
xx<-14+order(dd,decreasing=TRUE)
colnames(head(uniDat10[,xx]))
datMelt<-melt(uniDat10[,xx])
varNames<-paste(c("ESDR mean", "MS net div e=0.9", "Crown age", "BD ML net div", "MS net div e=0.5", "MS net div e=0", "PB ML net div", "ESDR skewness", "BD ML lambda", "BD ML mu"))
cols<-rep("light grey",10)
vcol<-"red"
cols[1]<-vcol; cols[3]<-vcol; cols[8]<-vcol

datMelt<-melt(uniDat10[,15:24])
cols<-c(rep("red",3),rep("light grey",7))
#varNames<-colnames(head(uniDat10[,15:24]))
varNames<-paste(c("Crown age","Tip DR mean", "Tip DR skewness", "PB net div", "BD lambda", "BD mu", "BD net div", "MS net div e=0","MS net div e=0.5", "MS net div e=0.9"))

pdf(file=paste("cladeLevel",bbone,"_PGLS_UNI_10Ma-timeSlice_VARS_byCat.pdf",sep=""),onefile=TRUE, width=6,height=4)
par(font.lab=2,cex.lab=0.8)#mar=c(6, 4.1, 4.1, 2.1))
boxplot(value ~ variable, data=datMelt, las=2, lwd=0.7,xaxt="n", font=2, ylim=c(-0.5,1.5),col=cols,ylab="Standardized effect on log (clade richness)")
axis(1, at=c(1:10), labels = FALSE)
text(x= c(1:10), y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), cex=0.8,font=2,srt = 45, labels = varNames, xpd = TRUE, adj=c(0.95,0.05))
dev.off()

## Now with PLOT CI (plotrix)
library(plotrix)

#dat<-uniDat10[,15:24]
#dat<-uniDat60[,15:24]
dat<-uniPGLS_fams_SCALED[,14:23]
varNames<-paste(c("Crown age","Tip DR mean", "Tip DR skewness", "PB net div", "BD lambda", "BD mu", "BD net div", "MS net div e=0","MS net div e=0.5", "MS net div e=0.9"))

means<-vector()
lower95s<-vector()
upper95s<-vector()
for (j in 1:length(dat)){
	means[j]<-mean(dat[,j])
	dd<-quantile(dat[,j], c(0.025,0.975))
	lower95s[j]<-dd[[1]]
	upper95s[j]<-dd[[2]]
}
cols<-c(rep("red",3),rep("black",7))

#pdf(file=paste("cladeLevel",bbone,"_PGLS_UNI_10Ma-timeSlice_VARS_byCat_95ci.pdf",sep=""),onefile=TRUE, width=6,height=4)
#pdf(file=paste("cladeLevel",bbone,"_PGLS_UNI_60Ma-timeSlice_VARS_byCat_95ci.pdf",sep=""),onefile=TRUE, width=6,height=4)
pdf(file=paste("cladeLevel",bbone,"_PGLS_UNI_FAMS_VARS_byCat_95ci.pdf",sep=""),onefile=TRUE, width=6,height=4)
plotCI(x=1:10, y=means,ui=upper95s,li=lower95s, ylim=c(-0.5,2),xaxt="n", xlab="", ylab="Standardized effect on log (clade richness)", err="y", lwd=2,col=cols,scol="black",pch=20,font.lab=2,cex.axis=0.95,cex.lab=1)
axis(1, at=c(1:10), labels = FALSE)
text(x= c(1:10), y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), cex=0.9,font=2,srt = 45, labels = varNames, xpd = TRUE, adj=c(0.95,0.05))
dev.off()


# 60 Ma

dd<-sapply(uniDat60[,15:24],mean)
xx<-14+order(dd,decreasing=TRUE)
colnames(head(uniDat60[,xx]))
datMelt<-melt(uniDat60[,xx])
varNames<-paste(c("ESDR mean", "MS net div e=0.9", "Crown age", "BD ML net div", "MS net div e=0.5", "MS net div e=0", "PB ML net div", "ESDR skewness", "BD ML lambda", "BD ML mu"))
cols<-rep("light grey",10)
vcol<-"red"
cols[1]<-vcol; cols[3]<-vcol; cols[8]<-vcol

datMelt<-melt(uniDat60[,15:24])
cols<-c(rep("red",3),rep("light grey",7))
#varNames<-colnames(head(uniDat10[,15:24]))
varNames<-paste(c("Crown age","Tip DR mean", "Tip DR skewness", "PB net div", "BD lambda", "BD mu", "BD net div", "MS net div e=0","MS net div e=0.5", "MS net div e=0.9"))

pdf(file=paste("cladeLevel",bbone,"_PGLS_UNI_60Ma-timeSlice_VARS_byCat.pdf",sep=""),onefile=TRUE, width=6,height=4)
par(font.lab=2,cex.lab=0.8)#mar=c(6, 4.1, 4.1, 2.1))
boxplot(value ~ variable, data=datMelt, las=2, xaxt="n", font=2, col=cols,ylab="Standardized effect on log (clade richness)")
axis(1, at=c(1:10), labels = FALSE)
text(x= c(1:10), y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), cex=0.8,font=2,srt = 45, labels = varNames, xpd = TRUE, adj=c(0.95,0.05))
dev.off()




##
# *** ACTUAL DATA VALUES:
# PGLS uni
allUni<-vector("list",length=numTrees)
for (i in 1:numTrees){
	ints_i<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniINTS.txt",sep=""))
	slopes_i<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniSLOPES.txt",sep=""))
	Ps_i<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniPs.txt",sep=""))
	tree<-rep(i,length(sliceTimes))
	allUni[[i]]<-cbind(sliceTimes,ints_i,slopes_i,Ps_i,tree)
}
allUni_100<-do.call(rbind,allUni)

uniDat_10<-allUni_100[which(allUni_100$sliceTimes==-10),]
uniDat_60<-allUni_100[which(allUni_100$sliceTimes==-60),]

# PGLS multi as based on ACTUAL values (not z-scores)
allSlopes_per_i<-vector("list",length=numTrees)
slice<-c(1:14)
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	dat<-read.table(paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlice_5Ma-to-70Ma_partSlopes.txt",sep=""))
	allSlopes_per_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
allDat<-do.call(rbind,allSlopes_per_i)

allDat_10<-allDat[which(allDat$sliceTimes==-10),]
allDat_60<-allDat[which(allDat$sliceTimes==-60),]

# 10 Ma slice
##
#pdf(file=paste("cladeLevel",bbone,"_scatter_timeSlices10_PGLSuni_Yaxis_PLOTTED_3part.pdf",sep=""),onefile=TRUE,width=4,height=10)
jpeg(file=paste("cladeLevel",bbone,"_scatter_timeSlices10_PGLSuni_Yaxis_PLOTTED_3part.jpg",sep=""),width=4,height=8, units="in", res=600, quality=100)
par(mfrow = c(3,1),oma = c(5,4,5,0) + 0.1, mar = c(3,1,1,1) + 0.1)
# type="n",
plot((richness)~MRCA,data=dat10, log="y",type="n",col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
mtext(text="Clade crown age (Ma)", side=1, padj=3, font=2,cex=1)

plot((richness)~DR_harm,data=dat10, log="y",type="n",col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
mtext(text="Clade ESDR mean (species / Ma)", side=1, padj=3, font=2,cex=1)

plot((richness)~DR_skew,data=dat10, log="y", type="n",col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
mtext(text="Clade ESDR skewness", side=1, padj=3, font=2,cex=1)

title(main="", xlab = "",
      ylab = "log(richness)",
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)
dev.off()


#pdf(file=paste("cladeLevel",bbone,"_scatter_timeSlices10_PGLSuni_ALL_PLOTTED_3part.pdf",sep=""),onefile=TRUE,width=4,height=10)
jpeg(file=paste("cladeLevel",bbone,"_scatter_timeSlices10_PGLSuni_ALL_PLOTTED_3part_new.jpg",sep=""),width=4,height=8, units="in", res=600, quality=100)
par(mfrow = c(3,1),oma = c(5,4,5,0) + 0.1, mar = c(3,1,1,1) + 0.1)
ltys<-1
lwds<-4
colRaw<-col2rgb("light blue")/255
cols<-rgb(colRaw[1],colRaw[2],colRaw[3],alpha=0.5)
#plot(log(richness)~MRCA,data=dat10, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
smoothScatter(log(dat10$richness)~dat10$MRCA,xlab="",ylab="",nbin=50,nrpoints=0,colramp=colorRampPalette(c("white", "black")))
for (i in 1:length(uniDat_10[,1])){
	abline(a=(uniDat_10$i_1[i]),b=(uniDat_10$MRCA[i]), untf=FALSE,col=cols, lwd=lwds,lty=ltys)
}
for (i in 1:length(uniDat_10[,1])){
	abline(a=(uniDat_10$i_1[i]),b=(uniDat_10$MRCA[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}

#mtext(text="Clade crown age (Ma)", side=1, padj=3, font=2,cex=1)

#plot(log(richness)~DR_harm,data=dat10, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
smoothScatter(log(dat10$richness)~dat10$DR_harm,xlab="",ylab="",nbin=50,nrpoints=0,colramp=colorRampPalette(c("white", "black")))
for (i in 1:length(uniDat_10[,1])){
	abline(a=(uniDat_10$i_2[i]),b=(uniDat_10$DR_harm[i]), untf=FALSE,col=cols, lwd=lwds,lty=ltys)
}
for (i in 1:length(uniDat_10[,1])){
	abline(a=(uniDat_10$i_2[i]),b=(uniDat_10$DR_harm[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}
#mtext(text="Clade ESDR mean (species / Ma)", side=1, padj=3, font=2,cex=1)

#plot(log(richness)~DR_skew,data=dat10, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
smoothScatter(log(dat10$richness)~dat10$DR_skew,xlab="",ylab="",nbin=50,nrpoints=0,colramp=colorRampPalette(c("white", "black")))
for (i in 1:length(uniDat_10[,1])){
	abline(a=(uniDat_10$i_3[i]),b=(uniDat_10$DR_skew[i]), untf=FALSE,col=cols, lwd=lwds,lty=ltys)
}
for (i in 1:length(uniDat_10[,1])){
	abline(a=(uniDat_10$i_3[i]),b=(uniDat_10$DR_skew[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}
#mtext(text="Clade ESDR skewness", side=1, padj=3, font=2,cex=1)

title(main="", xlab = "",
      ylab = "",
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)
dev.off()




# 60 Ma slice
##
#pdf(file=paste("cladeLevel",bbone,"_scatter_timeSlices60_PGLSuni_Yaxis_PLOTTED_3part.pdf",sep=""),onefile=TRUE,width=4,height=10)
jpeg(file=paste("cladeLevel",bbone,"_scatter_timeSlices60_PGLSuni_Yaxis_PLOTTED_3part.jpg",sep=""),width=4,height=8, units="in", res=600, quality=100)
par(mfrow = c(3,1),oma = c(5,4,5,0) + 0.1, mar = c(3,1,1,1) + 0.1)
# type="n",
plot((richness)~MRCA,data=dat60, log="y",type="n",col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
#mtext(text="Clade crown age (Ma)", side=1, padj=3, font=2,cex=1)

plot((richness)~DR_harm,data=dat60, log="y",type="n",col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
#mtext(text="Clade ESDR mean (species / Ma)", side=1, padj=3, font=2,cex=1)

plot((richness)~DR_skew,data=dat60, log="y", type="n",col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
#mtext(text="Clade ESDR skewness", side=1, padj=3, font=2,cex=1)

title(main="", xlab = "",
      ylab = "log(richness)",
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)
dev.off()


jpeg(file=paste("cladeLevel",bbone,"_scatter_timeSlices60_PGLSuni_ALL_PLOTTED_3part_new.jpg",sep=""),width=4,height=8, units="in", res=600, quality=100)
par(mfrow = c(3,1),oma = c(5,4,5,0) + 0.1, mar = c(3,1,1,1) + 0.1)
ltys<-1
lwds<-4
colRaw<-col2rgb("light blue")/255
cols<-rgb(colRaw[1],colRaw[2],colRaw[3],alpha=0.5)
#plot(log(richness)~MRCA,data=dat60, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
smoothScatter(log(dat60$richness)~dat60$MRCA,xlab="",ylab="",nbin=50,nrpoints=0,colramp=colorRampPalette(c("white", "black")))
for (i in 1:length(uniDat_60[,1])){
	abline(a=(uniDat_60$i_1[i]),b=(uniDat_60$MRCA[i]), untf=FALSE,col=cols, lwd=lwds,lty=ltys)
}
for (i in 1:length(uniDat_60[,1])){
	if (uniDat_60$p_1[i] > 0.05) {
	abline(a=(uniDat_60$i_1[i]),b=(uniDat_60$MRCA[i]), untf=FALSE,col=rgb(1,0,0,alpha=0.3), lwd=lwds,lty=ltys)		
	} else next
}
for (i in 1:length(uniDat_60[,1])){
	abline(a=(uniDat_60$i_1[i]),b=(uniDat_60$MRCA[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}

#mtext(text="Clade crown age (Ma)", side=1, padj=3, font=2,cex=1)

#plot(log(richness)~DR_harm,data=dat60, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
smoothScatter(log(dat60$richness)~dat60$DR_harm,xlab="",ylab="",nbin=50,nrpoints=0,colramp=colorRampPalette(c("white", "black")))
for (i in 1:length(uniDat_60[,1])){
	abline(a=(uniDat_60$i_2[i]),b=(uniDat_60$DR_harm[i]), untf=FALSE,col=cols, lwd=lwds,lty=ltys)
}
for (i in 1:length(uniDat_60[,1])){
	abline(a=(uniDat_60$i_2[i]),b=(uniDat_60$DR_harm[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}
#mtext(text="Clade ESDR mean (species / Ma)", side=1, padj=3, font=2,cex=1)

#plot(log(richness)~DR_skew,data=dat60, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
smoothScatter(log(dat60$richness)~dat60$DR_skew,xlab="",ylab="",nbin=50,nrpoints=0,colramp=colorRampPalette(c("white", "black")))
for (i in 1:length(uniDat_60[,1])){
	abline(a=(uniDat_60$i_3[i]),b=(uniDat_60$DR_skew[i]), untf=FALSE,col=cols, lwd=lwds,lty=ltys)
}
for (i in 1:length(uniDat_60[,1])){
	abline(a=(uniDat_60$i_3[i]),b=(uniDat_60$DR_skew[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}
#mtext(text="Clade ESDR skewness", side=1, padj=3, font=2,cex=1)

title(main="", xlab = "",
      ylab = "",
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)
dev.off()

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







#===============
# Make SIMULATIONS of the full MamPhy for testing these null expectations.
#===============
library(moments); library(nlme); library(ape); library(picante); library(phytools); library(geiger)

# First I want to get the ML birth and death rates for the whole mamPhy.
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/Diversification_analyses-Condamine/")
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut")
source("DR_functions.R")

library(foreach);library(doSNOW)
#cl = makeCluster(15, type = 'SOCK', outfile="")
cl = makeCluster(5, type = 'SOCK', outfile="")
registerDoSNOW(cl)

#ntrees=100
ntrees=5

foreach(i=1:ntrees, .packages=c('moments', 'nlme', 'ape', 'picante', 'phytools','geiger'), .combine=cbind, .verbose=TRUE) %dopar% {
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_nexus-and-newickTrees")

redos<-c(4, 28, 35, 78, 93)
# which backbone?
bbone<- "NDexp" #"FBD" # 
# tree
#mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_nexus.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 

mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",redos[i],"_nexus.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 
root<-max(node.age(mamPhy)$ages)

bd<-birthdeath(mamPhy)

extFrac<-bd$para[[1]] # d/b
netDiv<-bd$para[[2]] # b-d
lam<-netDiv/(1-extFrac)
mu<-lam-netDiv
rateTable<-cbind.data.frame(extFrac,netDiv,lam,mu)
write.table(rateTable,paste("MamPhy_SIMS_mamPhyE_",bbone,"_tree",redos[i],"_RATES.txt",sep=""))
# sim for mamPhy rates
	for (j in 1:10){
	simMam<-pbtree(b=lam,d=mu,n=5911,t=NULL, scale=root,nsim=1,type="continuous", extant.only=TRUE)
	if (class(simMam) == "NULL"){
		cat("another sim...")
		simMam<-pbtree(b=lam,d=mu,n=5911,t=NULL, scale=root,nsim=1,type="continuous", extant.only=TRUE)
		} else break
		cat("got one, simulating...")
	}
write.tree(simMam, paste("MamPhy_SIMS_mamPhyE_",bbone,"_tree",redos[i],".tre",sep=""))

# sim for e=0.8, same lam
muHiE<-0.8*lam
	for (j in 1:20){
	simHiE<-pbtree(b=lam,d=muHiE,n=5911,t=NULL, scale=root,nsim=1,type="continuous", extant.only=TRUE)
	if (class(simHiE) == "NULL"){
		cat("another sim...")
		simHiE<-pbtree(b=lam,d=muHiE,n=5911,t=NULL, scale=root,nsim=1,type="continuous", extant.only=TRUE)
		} else break
		cat("got one, simulating...")
	}
write.tree(simHiE, paste("MamPhy_SIMS_highE_0p8_",bbone,"_tree",redos[i],".tre",sep=""))

# sim for e=0.2, same net div 
muLoE<-0.2*lam
	simLoE<-pbtree(b=lam,d=muLoE,n=5911,t=NULL, scale=root,nsim=1,type="continuous", extant.only=TRUE)

write.tree(simLoE, paste("MamPhy_SIMS_lowE_0p2_",bbone,"_tree",redos[i],".tre",sep=""))

}

}


####
##############
# ALL trees simulated and WRITTEN
# want to LOAD back in, re-do the calculations
#######
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/Diversification_analyses-Condamine/")
library(moments); library(nlme); library(ape); library(picante); library(phytools); library(geiger)
source("DR_functions.R")

library(foreach);library(doSNOW)
#cl = makeCluster(15, type = 'SOCK', outfile="")
cl = makeCluster(5, type = 'SOCK', outfile="")
registerDoSNOW(cl)

#ntrees=100
ntrees=5

foreach(i=1:ntrees, .packages=c('moments', 'nlme', 'ape', 'picante', 'phytools','geiger'), .combine=cbind, .verbose=TRUE) %dopar% {
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_nexus-and-newickTrees")

redos<-c(4, 28, 35, 78, 93)

# which backbone?
bbone<- "NDexp" #"FBD" # 
# tree sims
simMam<-read.tree(paste("MamPhy_SIMS_mamPhyE_",bbone,"_tree",redos[i],".tre",sep=""))
simHiE<-read.tree(paste("MamPhy_SIMS_highE_0p8_",bbone,"_tree",redos[i],".tre",sep=""))
simLoE<-read.tree(paste("MamPhy_SIMS_lowE_0p2_",bbone,"_tree",redos[i],".tre",sep=""))

# CHANGING code -- the simulation name
sims<-c("mamPhyE", "highE_0p8", "lowE_0p2")
trees<-list(simMam,simHiE,simLoE)

for(q in 1:length(sims)){

# SIM TREE for DR calcs
tree1=scan(paste("MamPhy_SIMS_",sims[q],"_",bbone,"_tree",redos[i],".tre",sep=""), what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

# SIM TREE to work on
mamPhy<-trees[[q]]

# ES and DR on per-tip basis
#==================
# gives pairwise clade matrix from CAIC function
clade_matrix = readCAIC(tree1)
DR = 1/ES_v2(clade_matrix)
ES = ES_v2(clade_matrix)
res = cbind.data.frame(DR,ES)
res1 = res[order(rownames(res)),]
write.table(res1, file=paste("MamPhy_SIMS_",sims[q],"_",bbone,"_tree",redos[i],"_DRtips.txt",sep=""))

# Make time slices to create the clades
#=======================================
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery
rootTime<-max(node.age(mamPhy)$ages)

allCladeSets<-vector("list",length=numSlices)
allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSets[[j]]<-treeSlice(mamPhy, slice=rootTime-(sliceEvery*j), trivial=FALSE)
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

# record clade sizes per slice
lengths<-vector()
for (j in 1:length(allCladeSets)){
	lengths[j]<-length(allCladeSets[[j]])
}
names(lengths)<-allCladeSetNames 
lengths
write.table(lengths, file=paste("MamPhy_SIMS_",sims[q],"_",bbone,"_tree",redos[i],"_timeslice_Lengths.txt",sep=""))

# Calculate per-SLICE, per-clade summary values 
#===================================
# get node times for tree
btimes<-branching.times(mamPhy)

# yule function
ymle = function(tree){ (.subset2(tree,2)-1L)/sum(.subset2(tree,4)) } # this take the # of number of nodes in a tree (minus 1) / sum of branch lengths.

#cladeSetz<-c()
#for(j in 1:length(allCladeSets[[1]])){
#	cladeSetz[j]<-length(allCladeSets[[1]][[j]]$tip.label)
#	}

# do per-slice, per-clade calcs
for(j in 1:length(allCladeSets)){
cladeSet<-allCladeSets[[j]]

	# empty data frames to fill
	DR_harm<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	DR_skew<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	richness<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	MRCA<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	PB_Div<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD_Lam<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD_Mu<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD_Div<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD.ms0<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD.ms0p5<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD.ms0p9<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))

	for (k in 1:length(cladeSet)){
	cladeSp<-cladeSet[[k]]$tip.label
		x<-res1[match(cladeSp,rownames(res1)),"DR"]
		DR_harm[k,] <- 1/(mean(1/x))
		DR_skew[k,] <- skewness(x)
		richness[k,] <- length(cladeSp)
		node <- getMRCA(mamPhy, cladeSp)
		MRCA[k,] <- btimes[node-5911] #taking the height of SAMPLED tree
	if (length(cladeSp) > 2) {
	# Yule model
		PB_Div[k,]<-ymle(cladeSet[[k]])
		# BD model
		bd<-birthdeath(cladeSet[[k]])
		BD_Lam[k,]<-bd$para[[2]]/(1-bd$para[[1]])
		BD_Mu[k,]<-bd$para[[1]]*(bd$para[[2]]/(1-bd$para[[1]]))
		BD_Div[k,]<-bd$para[[2]]
		# BD Mag and Sand
		cladeSet[[k]]$root.edge<-0
	    BD.ms0[k,]<-bd.ms(phy=cladeSet[[k]], missing=0, epsilon=0, crown=TRUE) # Assuming no extinction
    	BD.ms0p5[k,]<-bd.ms(phy=cladeSet[[k]], missing=0, epsilon=0.5, crown=TRUE) # Assuming medium extinction 
     	BD.ms0p9[k,]<-bd.ms(phy=cladeSet[[k]], missing=0, epsilon=0.9, crown=TRUE) # Assuming high extinction
		} else NULL
	}
	res2<-cbind(DR_harm, DR_skew, richness, MRCA, PB_Div, BD_Lam, BD_Mu, BD_Div, BD.ms0, BD.ms0p5, BD.ms0p9, i, j*5)
	colnames(res2)<-c("DR_harm","DR_skew", "richness", "MRCA", "PB_Div", "BD_Lam", "BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9", "tree", "slice")
	
	rownames(res2)<-paste(i,"_",j,"_",c(1:length(res2[,1])),sep="")

	write.table(res2,paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",redos[i],"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""))
}

# create slice phys
slicePhys<-vector("list",length(allCladeSets))
for (k in 1:length(allCladeSets)){
cladeReps<-vector()
for (j in 1:length(allCladeSets[[k]])){
	cladeSp<-allCladeSets[[k]][[j]]$tip.label
	cladeReps[j]<-cladeSp[1]
	}
toDrop<-setdiff(mamPhy$tip.label,cladeReps)
slicePhys[[k]]<-drop.tip(mamPhy,toDrop)
slicePhys[[k]]$tip.label<-paste(i,"_",k,"_",c(1:length(slicePhys[[k]]$tip.label)),sep="")
}

# write slice phys
for(j in 1:length(slicePhys)){
	write.tree(slicePhys[[j]], file=paste("MamPhy_SIMS_",sims[q],"_",bbone,"_tree_",redos[i],"_slicePhy-5to70Ma_by5.trees",sep=""), append=TRUE)
}

########
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/Diversification_analyses-Condamine/")
library(moments); library(nlme); library(ape); library(picante); library(phytools); library(geiger)

library(foreach);library(doSNOW)
cl = makeCluster(30, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees=100

foreach(i=1:ntrees, .packages=c('moments', 'nlme', 'ape', 'picante', 'phytools','geiger'), .combine=cbind, .verbose=TRUE) %dopar% {
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_nexus-and-newickTrees")

# which backbone?
bbone<- "NDexp" #"FBD" # 

# CHANGING code -- the simulation name
sims<-c("mamPhyE", "highE_0p8", "lowE_0p2")

for(q in 1:length(sims)){

# Make time slices to create the clades
#=======================================
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery

allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

# read in slice phys
slicePhys<-read.tree(file=paste("MamPhy_SIMS_",sims[q],"_",bbone,"_tree_",i,"_slicePhy-5to70Ma_by5.trees",sep=""))

# multiply edges by 100 for corPagel
modTrees<-vector("list",length(slicePhys))
for (j in 1:length(slicePhys)){
	modTree<-slicePhys[[j]]
	modTree$edge.length<-modTree$edge.length*100
	modTrees[[j]]<-modTree
}

# do the PGLS as SCALED-- so that the EFFECT sizes are comparable.
#####
# load in tables and RESCALE as STANDARDIZED
#results<-vector("list",length(allCladeSetNames))
#for (j in 1: length(allCladeSetNames)){
#	res<-read.table(paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",redos[i],"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
#	varsScale<-scale(res, center=TRUE,scale=TRUE)
#	results[[j]]<-cbind(res$richness,varsScale)
#	colnames(results[[j]])<-c("richnessN",colnames(res))
#}

# Z-score STANDARDIZED from SQ MRCA and RAW others
# Z-score STANDARDIZED normal MRCA
results<-vector("list",length(allCladeSetNames))
for (j in 1: length(allCladeSetNames)){
	res<-read.table(paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
	vars<-res[,1:11]
	varsScale<-scale(vars, center=TRUE,scale=TRUE)
	results[[j]]<-cbind(res$richness,varsScale)
	colnames(results[[j]])<-c("richnessN",colnames(res[,1:11]))
}
 
# RAW but with SQ MRCA also
# results<-vector("list",length(allCladeSetNames))
# for (j in 1: length(allCladeSetNames)){
# 	res<-read.table(paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
# 	mrcaSQ<-(res[,"MRCA"])^2
# 	vars<-cbind(mrcaSQ, res)
# 	results[[j]]<-cbind(res$richness,vars)
# 	colnames(results[[j]])<-c("richnessN","mrcaSQ",colnames(res[,1:13]))
# }


# RE_SCALED -- MULTI- and UNI-VARIATE
sliceTimes<-seq(-5,-70,-5)
partSlopesPGLS_All<-data.frame(matrix(NA, nrow = length(results), ncol = 16),row.names=sliceTimes)
colnames(partSlopesPGLS_All)<-c("int","MRCA","DR_harm","DR_skew","SE1","SE2","SE3","SE4","Lam_low","Lam_mean","Lam_up","AIC","Pval1","Pval2","Pval3","Pval4")
#partSlopesPGLS_12<-data.frame(matrix(NA, nrow = length(results), ncol = 6),row.names=sliceTimes)
#colnames(partSlopesPGLS_12)<-c("int","MRCA","DR_harm","AIC","Pval1","Pval2")

uniPGLS_allSlopes<-data.frame(matrix(NA, nrow = length(results), ncol = 10),row.names=sliceTimes)
colnames(uniPGLS_allSlopes)<-c("MRCA","DR_harm","DR_skew","PB_Div","BD_Lam","BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9")
uniPGLS_allSEs<-data.frame(matrix(NA, nrow = length(results), ncol = 10),row.names=sliceTimes)
colnames(uniPGLS_allSEs)<-c(paste("SE_",1:10,sep=""))
uniPGLS_allInts<-data.frame(matrix(NA, nrow = length(results), ncol = 10),row.names=sliceTimes)
colnames(uniPGLS_allInts)<-c(paste("i_",1:10,sep=""))
uniPGLS_allPs<-data.frame(matrix(NA, nrow = length(results), ncol = 10),row.names=sliceTimes)
colnames(uniPGLS_allPs)<-c(paste("p_",1:10,sep=""))
uniPGLS_allLams<-data.frame(matrix(NA, nrow = length(results), ncol = 10),row.names=sliceTimes)
colnames(uniPGLS_allLams)<-c(paste("lam_",1:10,sep=""))

for (j in 1:length(results)){
	#rownames(results[[j]])<-slicePhys[[j]]$tip.label
	#res_OK<-results[[j]]
	#cladeData<-treedata(slicePhys[[j]],na.omit(res_OK))
	cladeData<-treedata(modTrees[[j]],na.omit(results[[j]]))
	dat<-as.data.frame(cladeData$data)
	#dat<-as.data.frame(na.omit(results[[j]]))

	form<-(log(richnessN) ~ MRCA + DR_harm + DR_skew)
#	fit<-gls(form, data=dat, method="ML")
#	fit<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
		for (p in seq(0,1,by=0.01)) {possibleError <- tryCatch(
		      gls(form, correlation=corPagel(value=p,phy=cladeData$phy), data=dat, method="ML"),
		      error=function(e) e)
		if(inherits(possibleError, "gls")) break		
		if(inherits(possibleError, "error")) next}
		fit<-possibleError
	sum<-summary(fit)
	#int95s<-intervals(fit,which="var-cov")

	for(k in 1:8){
	partSlopesPGLS_All[j,k]<-round(sum$tTable[k],digits=3)
	}
	#partSlopesPGLS_All[j,9]<-round(int95s[[1]][1],digits=3)
	partSlopesPGLS_All[j,10]<-round(sum$modelStruct[[1]][[1]],digits=3)
	#partSlopesPGLS_All[j,11]<-round(int95s[[1]][3],digits=3)
	partSlopesPGLS_All[j,12]<-round(sum$AIC,digits=0)
	partSlopesPGLS_All[j,13]<-round(sum$tTable[13], digits=3)
	partSlopesPGLS_All[j,14]<-round(sum$tTable[14], digits=3)
	partSlopesPGLS_All[j,15]<-round(sum$tTable[15], digits=3)
	partSlopesPGLS_All[j,16]<-round(sum$tTable[16], digits=3)

	vars<-c("MRCA","DR_harm","DR_skew","PB_Div","BD_Lam","BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9")
	for(k in 1:length(uniPGLS_allSlopes)){
		form<-as.formula(paste("log(richnessN) ~ ", vars[k],sep=""))
#		fit1<-gls(form, data=dat, method="ML")
#		fit1<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
		for (p in seq(0,1,by=0.01)) {possibleError <- tryCatch(
		      gls(form, correlation=corPagel(value=p,phy=cladeData$phy), data=dat, method="ML"),
		      error=function(e) e)
		if(inherits(possibleError, "gls")) break		
		if(inherits(possibleError, "error")) next}
		fit1<-possibleError
		sum1<-summary(fit1)
		
		uniPGLS_allInts[j,k]<-round(sum1$tTable[1], digits=3)
		uniPGLS_allSlopes[j,k]<-round(sum1$tTable[2], digits=3)
		uniPGLS_allSEs[j,k]<-round(sum1$tTable[4], digits=3)
		uniPGLS_allPs[j,k]<-round(sum1$tTable[8], digits=3)
		uniPGLS_allLams[j,k]<-round(sum1$modelStruct[[1]][[1]], digits=3)
	}	

}
corr<-"PAGEL"

write.table(partSlopesPGLS_All,paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLSmulti_",corr,"_noLamCI_timeSlices_SCALED_partSlopes123.txt",sep=""))
uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs,uniPGLS_allSEs,uniPGLS_allLams)#,uniPGLS_allLams)
write.table(uniPGLS,paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLS_",corr,"_noLamCI_timeSlices_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))

#write.table(partSlopesPGLS_All,paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLSmulti_PAGEL_timeSlices_SCALED_partSlopes123.txt",sep=""))
##write.table(partSlopesPGLS_12,paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLSmulti_PAGEL_timeSlices_SCALED_partSlopes12.txt",sep=""))
#
#uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs,uniPGLS_allSEs,uniPGLS_allLams)#,uniPGLS_allLams)
#write.table(uniPGLS,paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLS_PAGEL_timeSlices_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))

#write.table(partSlopesPGLS_All,paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLSmulti_timeSlices_mrcaSQ_RAW_partSlopes123.txt",sep=""))
#write.table(partSlopesPGLS_12,paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLSmulti_timeSlices_mrcaSQ_RAW_partSlopes12.txt",sep=""))
#
#uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs)#,uniPGLS_allLams)
#write.table(uniPGLS,paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLS_timeSlices_mrcaSQ_RAW_uniINTS_uniSLO_uniP.txt",sep=""))

}

}





########
######
# Now also do the SIM PGLS as UN-SCALED -- so as to compare to the real data.
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/Diversification_analyses-Condamine/")
library(moments); library(nlme); library(ape); library(picante); library(phytools); library(geiger)

library(foreach);library(doSNOW)
cl = makeCluster(30, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees=100

foreach(i=1:ntrees, .packages=c('moments', 'nlme', 'ape', 'picante', 'phytools','geiger'), .combine=cbind, .verbose=TRUE) %dopar% {
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_nexus-and-newickTrees")

# which backbone?
bbone<- "NDexp" #"FBD" # 

# CHANGING code -- the simulation name
sims<-c("mamPhyE", "highE_0p8", "lowE_0p2")

for(q in 1:length(sims)){

# Make time slices to create the clades
#=======================================
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery

allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

# read in slice phys
slicePhys<-read.tree(file=paste("MamPhy_SIMS_",sims[q],"_",bbone,"_tree_",i,"_slicePhy-5to70Ma_by5.trees",sep=""))

# do the SIMULATION PGLS as ACTUAL values-- UN-SCALED 
#####
# load in raw tables
results<-vector("list",length(allCladeSetNames))
for (j in 1: length(allCladeSetNames)){
	results[[j]]<-read.table(paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
}

# ACTUAL vals -- MULTI- and UNI-VARIATE
sliceTimes<-seq(-5,-70,-5)
partSlopesPGLS_All<-data.frame(matrix(NA, nrow = length(results), ncol = 8),row.names=sliceTimes)
colnames(partSlopesPGLS_All)<-c("int","MRCA","DR_harm","DR_skew","AIC","Pval1","Pval2","Pval3")

partSlopesPGLS_12<-data.frame(matrix(NA, nrow = length(results), ncol = 6),row.names=sliceTimes)
colnames(partSlopesPGLS_12)<-c("int","MRCA","DR_harm","AIC","Pval1","Pval2")

uniPGLS_allSlopes<-data.frame(matrix(NA, nrow = length(results), ncol = 10),row.names=sliceTimes)
colnames(uniPGLS_allSlopes)<-c("MRCA","DR_harm","DR_skew","PB_Div","BD_Lam","BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9")
uniPGLS_allInts<-data.frame(matrix(NA, nrow = length(results), ncol = 10),row.names=sliceTimes)
colnames(uniPGLS_allInts)<-c(paste("i_",1:10,sep=""))
uniPGLS_allPs<-data.frame(matrix(NA, nrow = length(results), ncol = 10),row.names=sliceTimes)
colnames(uniPGLS_allPs)<-c(paste("p_",1:10,sep=""))

for (j in 1:length(results)){
	#rownames(results[[j]])<-slicePhys[[j]]$tip.label
	res_OK<-results[[j]]
	cladeData<-treedata(slicePhys[[j]],na.omit(res_OK))
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

	form<-(log(richness) ~ MRCA + DR_harm)
	fit2<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
	sum<-summary(fit2)

	partSlopesPGLS_12[j,1]<-round(sum$coef[[1]],digits=3)
	partSlopesPGLS_12[j,2]<-round(sum$coef[[2]],digits=3)
	partSlopesPGLS_12[j,3]<-round(sum$coef[[3]],digits=3)
	partSlopesPGLS_12[j,4]<-round(sum$AIC,digits=0)
	partSlopesPGLS_12[j,5]<-round(sum$tTable[11], digits=3)
	partSlopesPGLS_12[j,6]<-round(sum$tTable[12], digits=3)

	vars<-c("MRCA","DR_harm","DR_skew","PB_Div","BD_Lam","BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9")
	for(k in 1:length(uniPGLS_allSlopes)){
		form<-as.formula(paste("log(richness) ~ ", vars[k],sep=""))
		fit1<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
		sum1<-summary(fit1)
		
		uniPGLS_allInts[j,k]<-round(sum1$coef[[1]],digits=3)
		uniPGLS_allSlopes[j,k]<-round(sum1$coef[[2]],digits=3)
		uniPGLS_allPs[j,k]<-round(sum1$tTable[8], digits=3)
		#uniPGLS_allLams[j,k]<-round(sum1$modelStruct[[1]][[1]], digits=3)
	}	

}
write.table(partSlopesPGLS_All,paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLSmulti_timeSlices_ACTUAL_partSlopes123.txt",sep=""))
write.table(partSlopesPGLS_12,paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLSmulti_timeSlices_ACTUAL_partSlopes12.txt",sep=""))

uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs)#,uniPGLS_allLams)
write.table(uniPGLS,paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLS_timeSlices_ACTUAL_uniINTS_uniSLO_uniP.txt",sep=""))

}

}

######
# Now RELOAD the results files (are they all there?)
library(moments); library(nlme); library(ape); library(picante); library(phytools); library(geiger)

# which backbone?
bbone<- "NDexp" #"FBD" # 

# CHANGING code -- the simulation name
sims<-c("mamPhyE", "highE_0p8", "lowE_0p2")




######
# SUMMARIZE the SIMULATIONS
####
# Load data back in and COMBINE all slice clades per-tree, then across all trees
# Post doing PGLS
#===================================
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/Diversification_analyses-Condamine/")
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_SIMULATIONS_on_NDexp")
#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_SIMULATIONS_on_NDexp")
bbone<- "NDexp" #"FBD" # 
library(nlme); library(MASS); library(AICcmodavg); library(ape); library(phytools)

# summarize exFract of the MAMPHY-E runs...
rates<-vector("list",length=100)
extinctFract<-c()
for (i in 1:100){
	rates[[i]]<-read.table(file=paste("MamPhy_SIMS_mamPhyE_",bbone,"_tree",i,"_RATES.txt",sep=""))
	extinctFract[i]<-rates[[i]][[1]]
}
mean(extinctFract) # [1] 0.6546818

numTrees<-100
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery
sliceTimes<-seq(-5,-70,-5)

allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

# set sim
sims<-c("mamPhyE", "highE_0p8", "lowE_0p2")

#######
# ACTUAL values (non-standardized)
#######
# SIMS - MULTIVARIATE
q<-1
SIMS_allSlopes_per_i<-vector("list",length=numTrees)
slice<-c(1:14)
# load in data as is. 
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	#dat<-read.table(paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLSmulti_timeSlices_ACTUAL_partSlopes123.txt",sep=""))
	dat<-read.table(paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLSmulti_timeSlices_mrcaSQ_RAW_partSlopes123.txt",sep=""))
	SIMS_allSlopes_per_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
SIMS_mamPhyE_allDat_ACTUAL<-do.call(rbind,SIMS_allSlopes_per_i)
write.table(SIMS_mamPhyE_allDat_ACTUAL,file=paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_ALL100-TABLE_PGLSmulti_timeSlices_mrcaSQ_RAW_partSlopes123.txt",sep=""))

q<-2
SIMS_allSlopes_per_i<-vector("list",length=numTrees)
# load in data as is. 
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	dat<-read.table(paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLSmulti_timeSlices_mrcaSQ_RAW_partSlopes123.txt",sep=""))
	SIMS_allSlopes_per_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
SIMS_highE_allDat_ACTUAL<-do.call(rbind,SIMS_allSlopes_per_i)
write.table(SIMS_highE_allDat_ACTUAL,file=paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_ALL100-TABLE_PGLSmulti_timeSlices_mrcaSQ_RAW_partSlopes123.txt",sep=""))

q<-3
SIMS_allSlopes_per_i<-vector("list",length=numTrees)
# load in data as is. 
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	dat<-read.table(paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLSmulti_timeSlices_mrcaSQ_RAW_partSlopes123.txt",sep=""))
	SIMS_allSlopes_per_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
SIMS_lowE_allDat_ACTUAL<-do.call(rbind,SIMS_allSlopes_per_i)
write.table(SIMS_lowE_allDat_ACTUAL,file=paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_ALL100-TABLE_PGLSmulti_timeSlices_mrcaSQ_RAW_partSlopes123.txt",sep=""))

# SIMS - UNI-VARIATE
q<-1
SIMS_allDat_UNI_i<-vector("list",length=numTrees)
# load in data as is. 
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	dat<-read.table(paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLS_timeSlices_mrcaSQ_RAW_uniINTS_uniSLO_uniP.txt",sep=""))
	SIMS_allDat_UNI_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
SIMS_mamPhyE_allDat_UNI_ACTUAL<-do.call(rbind,SIMS_allDat_UNI_i)
write.table(SIMS_mamPhyE_allDat_UNI_ACTUAL,file=paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_ALL100-TABLE_PGLS-UNI_timeSlices_mrcaSQ_RAW_all.txt",sep=""))

q<-2
SIMS_allDat_UNI_i<-vector("list",length=numTrees)
# load in data as is. 
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	dat<-read.table(paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLS_timeSlices_mrcaSQ_RAW_uniINTS_uniSLO_uniP.txt",sep=""))
	SIMS_allDat_UNI_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
SIMS_highE_allDat_UNI_ACTUAL<-do.call(rbind,SIMS_allDat_UNI_i)
write.table(SIMS_highE_allDat_UNI_ACTUAL,file=paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_ALL100-TABLE_PGLS-UNI_timeSlices_mrcaSQ_RAW_all.txt",sep=""))

q<-3
SIMS_allDat_UNI_i<-vector("list",length=numTrees)
# load in data as is. 
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	dat<-read.table(paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLS_timeSlices_mrcaSQ_RAW_uniINTS_uniSLO_uniP.txt",sep=""))
	SIMS_allDat_UNI_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
SIMS_lowE_allDat_UNI_ACTUAL<-do.call(rbind,SIMS_allDat_UNI_i)
write.table(SIMS_lowE_allDat_UNI_ACTUAL,file=paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_ALL100-TABLE_PGLS-UNI_timeSlices_mrcaSQ_RAW_all.txt",sep=""))



#######
# SCALED values (standardized)
#######
# SIMS - MULTIVARIATE
slice<-c(1:14)
numTrees<-100

for (q in 1:3){
	SIMS_allSlopes_per_i<-vector("list",length=numTrees)
	for(i in 1:numTrees){
		tree<-rep(i,length(slice))
		dat<-read.table(paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLSmulti_BROWNIAN_timeSlices_SCALED_partSlopes123.txt",sep=""))
		SIMS_allSlopes_per_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
	}
	SIMS_allDat<-do.call(rbind,SIMS_allSlopes_per_i)
	write.table(SIMS_allDat,file=paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_ALL100-TABLE_PGLSmulti_BROWNIAN_timeSlices_SCALED_partSlopes123.txt",sep=""))
}
for (q in 1:3){
	SIMS_allSlopes_per_i<-vector("list",length=numTrees)
	for(i in 1:numTrees){
		tree<-rep(i,length(slice))
		dat<-read.table(paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLSmulti_noTREE_timeSlices_SCALED_partSlopes123.txt",sep=""))
		SIMS_allSlopes_per_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
	}
	SIMS_allDat_gls<-do.call(rbind,SIMS_allSlopes_per_i)
	write.table(SIMS_allDat_gls,file=paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_ALL100-TABLE_PGLSmulti_noTREE_timeSlices_SCALED_partSlopes123.txt",sep=""))
}

# SIMS - UNI-VARIATE
for (q in 1:3){
	SIMS_allDat_UNI_i<-vector("list",length=numTrees)
	# load in data as is. 
	for(i in 1:numTrees){
		tree<-rep(i,length(slice))
		dat<-read.table(paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLS_BROWNIAN_timeSlices_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))
		SIMS_allDat_UNI_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
	}
	SIMS_allDat_UNI<-do.call(rbind,SIMS_allDat_UNI_i)
	write.table(SIMS_allDat_UNI,file=paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_ALL100-TABLE_PGLS-UNI_BROWNIAN_timeSlices_SCALED_all.txt",sep=""))
}
for (q in 1:3){
	SIMS_allDat_UNI_i<-vector("list",length=numTrees)
	# load in data as is. 
	for(i in 1:numTrees){
		tree<-rep(i,length(slice))
		dat<-read.table(paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLS_noTREE_timeSlices_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))
		SIMS_allDat_UNI_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
	}
	SIMS_allDat_UNI_gls<-do.call(rbind,SIMS_allDat_UNI_i)
	write.table(SIMS_allDat_UNI_gls,file=paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_ALL100-TABLE_PGLS-UNI_noTREE_timeSlices_SCALED_all.txt",sep=""))
}
# ^^ Amazing, that actually summarizes ALL of those sim runs!!
#####
# now RELOAD all the SIMS-- Incl the ACTUALs
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_SIMULATIONS_on_NDexp")
library(nlme); library(MASS); library(AICcmodavg); library(ape); library(phytools)
bbone<- "NDexp" #"FBD" # 
numTrees<-100
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery
sliceTimes<-seq(-5,-70,-5)

# Load sims
sims<-c("mamPhyE", "highE_0p8", "lowE_0p2")

SIMS_allDat_MULTI_actual<-vector("list",length=3)
SIMS_allDat_MULTI_scaled<-vector("list",length=3)
SIMS_allDat_UNI_actual<-vector("list",length=3)
SIMS_allDat_UNI_scaled<-vector("list",length=3)
for (i in 1:length(SIMS_allDat_MULTI_actual)){
	#SIMS_allDat_MULTI_actual[[i]]<-read.table(paste("MamPhy_SIMS_",sims[i],"_",bbone,"_sample100_ALL100-TABLE_PGLSmulti_BROWNIAN_timeSlices_ACTUAL_partSlopes123.txt",sep=""))
	SIMS_allDat_MULTI_scaled[[i]]<-read.table(paste("MamPhy_SIMS_",sims[i],"_",bbone,"_sample100_ALL100-TABLE_PGLSmulti_BROWNIAN_timeSlices_SCALED_partSlopes123.txt",sep=""))
	#SIMS_allDat_UNI_actual[[i]]<-read.table(paste("MamPhy_SIMS_",sims[i],"_",bbone,"_sample100_ALL100-TABLE_PGLS-UNI_BROWNIAN_timeSlices_ACTUAL_all.txt",sep=""))
	SIMS_allDat_UNI_scaled[[i]]<-read.table(paste("MamPhy_SIMS_",sims[i],"_",bbone,"_sample100_ALL100-TABLE_PGLS-UNI_BROWNIAN_timeSlices_SCALED_all.txt",sep=""))

	#SIMS_allDat_MULTI_actual[[i]]<-read.table(paste("MamPhy_SIMS_",sims[i],"_",bbone,"_sample100_ALL100-TABLE_PGLSmulti_timeSlices_mrcaSQ_RAW_partSlopes123.txt",sep=""))
	#SIMS_allDat_MULTI_scaled[[i]]<-read.table(paste("MamPhy_SIMS_",sims[i],"_",bbone,"_sample100_ALL100-TABLE_PGLSmulti_timeSlices_mrcaSQ_SCALED_partSlopes123.txt",sep=""))
	#SIMS_allDat_UNI_actual[[i]]<-read.table(paste("MamPhy_SIMS_",sims[i],"_",bbone,"_sample100_ALL100-TABLE_PGLS-UNI_timeSlices_mrcaSQ_RAW_all.txt",sep=""))
	#SIMS_allDat_UNI_scaled[[i]]<-read.table(paste("MamPhy_SIMS_",sims[i],"_",bbone,"_sample100_ALL100-TABLE_PGLS-UNI_timeSlices_mrcaSQ_SCALED_all.txt",sep=""))
}


ranges1<-cbind(range(SIMS_allDat_MULTI_scaled[[1]][,"MRCA"]),range(SIMS_allDat_MULTI_scaled[[1]][,"DR_harm"]),range(SIMS_allDat_MULTI_scaled[[1]][,"DR_skew"]),range(SIMS_allDat_MULTI_scaled[[1]][,"Pval1"]),range(SIMS_allDat_MULTI_scaled[[1]][,"Pval2"]),range(SIMS_allDat_MULTI_scaled[[1]][,"Pval3"]))
ranges2<-cbind(range(SIMS_allDat_MULTI_scaled[[2]][,"MRCA"]),range(SIMS_allDat_MULTI_scaled[[2]][,"DR_harm"]),range(SIMS_allDat_MULTI_scaled[[2]][,"DR_skew"]),range(SIMS_allDat_MULTI_scaled[[2]][,"Pval1"]),range(SIMS_allDat_MULTI_scaled[[2]][,"Pval2"]),range(SIMS_allDat_MULTI_scaled[[2]][,"Pval3"]))
ranges3<-cbind(range(SIMS_allDat_MULTI_scaled[[3]][,"MRCA"]),range(SIMS_allDat_MULTI_scaled[[3]][,"DR_harm"]),range(SIMS_allDat_MULTI_scaled[[3]][,"DR_skew"]),range(SIMS_allDat_MULTI_scaled[[3]][,"Pval1"]),range(SIMS_allDat_MULTI_scaled[[3]][,"Pval2"]),range(SIMS_allDat_MULTI_scaled[[3]][,"Pval3"]))

allRanges_scaled<-rbind(ranges1,ranges2,ranges3)
colnames(allRanges_scaled)<-c("MRCA","DR_mean","DR_skew","Pval1","Pval2","Pval3")
rownames(allRanges_scaled)<-c("min_0.65","max_0.65","min_0.8","max_0.8","min_0.2","max_0.2")
allRanges_scaled

ranges1<-cbind(range(SIMS_allDat_MULTI_actual[[1]][,"MRCA"]),range(SIMS_allDat_MULTI_actual[[1]][,"DR_harm"]),range(SIMS_allDat_MULTI_actual[[1]][,"DR_skew"]),range(SIMS_allDat_MULTI_actual[[1]][,"Pval1"]),range(SIMS_allDat_MULTI_actual[[1]][,"Pval2"]),range(SIMS_allDat_MULTI_actual[[1]][,"Pval3"]))
ranges2<-cbind(range(SIMS_allDat_MULTI_actual[[2]][,"MRCA"]),range(SIMS_allDat_MULTI_actual[[2]][,"DR_harm"]),range(SIMS_allDat_MULTI_actual[[2]][,"DR_skew"]),range(SIMS_allDat_MULTI_actual[[2]][,"Pval1"]),range(SIMS_allDat_MULTI_actual[[2]][,"Pval2"]),range(SIMS_allDat_MULTI_actual[[2]][,"Pval3"]))
ranges3<-cbind(range(SIMS_allDat_MULTI_actual[[3]][,"MRCA"]),range(SIMS_allDat_MULTI_actual[[3]][,"DR_harm"]),range(SIMS_allDat_MULTI_actual[[3]][,"DR_skew"]),range(SIMS_allDat_MULTI_actual[[3]][,"Pval1"]),range(SIMS_allDat_MULTI_actual[[3]][,"Pval2"]),range(SIMS_allDat_MULTI_actual[[3]][,"Pval3"]))

allRanges_actual<-rbind(ranges1,ranges2,ranges3)
colnames(allRanges_actual)<-c("MRCA","DR_mean","DR_skew","Pval1","Pval2","Pval3")
rownames(allRanges_actual)<-c("min_0.65","max_0.65","min_0.8","max_0.8","min_0.2","max_0.2")
allRanges_actual


# Load EMPIRICAL DATA
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut")
# MCC tree
mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 
rootTime=max(nodeHeights(mamPhy)[,2])

# MULTI -- SCALED
allSlopes_per_i<-vector("list",length=numTrees)
slice<-c(1:14)
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	#dat<-read.table(paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlices_SCALED_partSlopes123.txt",sep=""))
	#dat<-read.table(paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlice_5Ma-to-70Ma_partSlopes.txt",sep=""))
	dat<-read.table(paste(bbone,"_sample100_",i,"_PGLSmulti_BROWNIAN_timeSlices_SCALED_partSlopes123.txt",sep=""))
	allSlopes_per_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
allDat_MULTI_scaled<-do.call(rbind,allSlopes_per_i)
write.table(allDat_MULTI_scaled,file=paste(bbone,"_sample100_ALL100-TABLE_PGLSmulti_BROWNIAN_timeSlices_SCALED_3vars.txt",sep=""))

allSlopes_per_i<-vector("list",length=numTrees)
slice<-c(1:14)
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	#dat<-read.table(paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlices_SCALED_partSlopes123.txt",sep=""))
	#dat<-read.table(paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlice_5Ma-to-70Ma_partSlopes.txt",sep=""))
	dat<-read.table(paste(bbone,"_sample100_",i,"_PGLSmulti_noTREE_timeSlices_SCALED_partSlopes123.txt",sep=""))
	allSlopes_per_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
allDat_MULTI_scaled_gls<-do.call(rbind,allSlopes_per_i)
write.table(allDat_MULTI_scaled_gls,file=paste(bbone,"_sample100_ALL100-TABLE_PGLSmulti_noTREE_timeSlices_SCALED_3vars.txt",sep=""))

allSlopes_per_i<-vector("list",length=numTrees)
slice<-c(1:14)
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	dat<-read.table(paste(bbone,"_sample100_",i,"_PGLSmulti_BROWNIAN_timeSlices_wPercentSamp_SCALED_partSlopes123.txt",sep=""))
	allSlopes_per_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
allDat_MULTI_scaled_4var<-do.call(rbind,allSlopes_per_i)
write.table(allDat_MULTI_scaled_4var,file=paste(bbone,"_sample100_ALL100-TABLE_PGLSmulti_BROWNIAN_timeSlices_wPercentSamp_SCALED_4vars.txt",sep=""))

# UNI -- SCALED
allDat_UNI_i<-vector("list",length=numTrees)
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	dat<-read.table(paste(bbone,"_sample100_",i,"_PGLS_BROWNIAN_timeSlices_wPercentSamp_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))
	#ints<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniINTS.txt",sep=""))
	#slopes<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniSLOPES.txt",sep=""))
	#Ps<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniPs.txt",sep=""))
	#dat<-cbind(ints,slopes,Ps)
	allDat_UNI_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
allDat_UNI_scaled<-do.call(rbind,allDat_UNI_i)
write.table(allDat_UNI_scaled,file=paste(bbone,"_sample100_ALL100-TABLE_PGLS_BROWNIAN_timeSlices_wPercentSamp_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))

allDat_UNI_i<-vector("list",length=numTrees)
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	dat<-read.table(paste(bbone,"_sample100_",i,"_PGLS_noTREE_timeSlices_wPercentSamp_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))
	#ints<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniINTS.txt",sep=""))
	#slopes<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniSLOPES.txt",sep=""))
	#Ps<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniPs.txt",sep=""))
	#dat<-cbind(ints,slopes,Ps)
	allDat_UNI_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
allDat_UNI_scaled_gls<-do.call(rbind,allDat_UNI_i)
write.table(allDat_UNI_scaled_gls,file=paste(bbone,"_sample100_ALL100-TABLE_PGLS_noTREE_timeSlices_wPercentSamp_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))

# MULTI -- ACTUAL
allSlopes_per_i<-vector("list",length=numTrees)
slice<-c(1:14)
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	dat<-read.table(paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlice_5Ma-to-70Ma_partSlopes.txt",sep=""))
	#dat<-read.table(paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlices_mrcaSQ_RAW_partSlopes123.txt",sep=""))
	allSlopes_per_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
allDat_MULTI_actual<-do.call(rbind,allSlopes_per_i)
write.table(allDat_MULTI_actual,file=paste(bbone,"_sample100_ALL100-TABLE_PGLSmulti_BROWNIAN_timeSlices_ACTUAL_3vars.txt",sep=""))

# UNI -- ACTUAL
allDat_UNI_i<-vector("list",length=numTrees)
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	ints<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniINTS.txt",sep=""))
	slopes<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniSLOPES.txt",sep=""))
	Ps<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniPs.txt",sep=""))
	dat<-cbind(ints,slopes,Ps)
	#dat<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlices_mrcaSQ_RAW_uniINTS_uniSLO_uniP.txt",sep=""))	
	allDat_UNI_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
allDat_UNI_actual<-do.call(rbind,allDat_UNI_i)
write.table(allDat_UNI_actual,file=paste(bbone,"_sample100_ALL100-TABLE_PGLSmulti_BROWNIAN_timeSlices_ACTUAL_uniINTS_uniSLO_uniP.txt",sep=""))


# RELOAD THOSE ALL...
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/")
library(nlme); library(MASS); library(AICcmodavg); library(ape); library(phytools)
bbone<- "NDexp" #"FBD" # 
numTrees<-100
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery
sliceTimes<-seq(-5,-70,-5)

# Load sims
sims<-c("mamPhyE", "highE_0p8", "lowE_0p2")
SIMS_allDat_MULTI_scaled<-vector("list",length=3)
SIMS_allDat_UNI_scaled<-vector("list",length=3)
for (i in 1:length(SIMS_allDat_MULTI_scaled)){
	SIMS_allDat_MULTI_scaled[[i]]<-read.table(paste("MamPhy_SIMS_",sims[i],"_",bbone,"_sample100_ALL100-TABLE_PGLSmulti_BROWNIAN_timeSlices_SCALED_partSlopes123.txt",sep=""))
	SIMS_allDat_UNI_scaled[[i]]<-read.table(paste("MamPhy_SIMS_",sims[i],"_",bbone,"_sample100_ALL100-TABLE_PGLS-UNI_BROWNIAN_timeSlices_SCALED_all.txt",sep=""))
}

# load Empirical
allDat_MULTI_scaled<-read.table(paste(bbone,"_sample100_ALL100-TABLE_PGLSmulti_BROWNIAN_timeSlices_SCALED_3vars.txt",sep=""))
allDat_MULTI_scaled_gls<-read.table(file=paste(bbone,"_sample100_ALL100-TABLE_PGLSmulti_noTREE_timeSlices_SCALED_3vars.txt",sep=""))
allDat_MULTI_scaled_4var<-read.table(file=paste(bbone,"_sample100_ALL100-TABLE_PGLSmulti_BROWNIAN_timeSlices_wPercentSamp_SCALED_4vars.txt",sep=""))
allDat_UNI_scaled<-read.table(file=paste(bbone,"_sample100_ALL100-TABLE_PGLS_BROWNIAN_timeSlices_wPercentSamp_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))
allDat_UNI_scaled_gls<-read.table(file=paste(bbone,"_sample100_ALL100-TABLE_PGLS_noTREE_timeSlices_wPercentSamp_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))
allDat_MULTI_actual<-read.table(file=paste(bbone,"_sample100_ALL100-TABLE_PGLSmulti_BROWNIAN_timeSlices_ACTUAL_3vars.txt",sep=""))
allDat_UNI_actual<-read.table(file=paste(bbone,"_sample100_ALL100-TABLE_PGLSmulti_BROWNIAN_timeSlices_ACTUAL_uniINTS_uniSLO_uniP.txt",sep=""))

# MCC tree
mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 
rootTime=max(nodeHeights(mamPhy)[,2])

#save.image(file="cladeSlice_workspace_forPlots_update.Rdata")
load(file="cladeSlice_workspace_forPlots_update.Rdata")


# RE-SCALE the SIM actuals to match empiricals...
## FOR the sims scaled - MULTIVARIATE

simActual_Multis<-vector("list",length=3)
for (i in 1:3){
Sim_multi_scaled<-SIMS_allDat_MULTI_scaled[[i]]
# now do the RATIOS -- as MULTIVARIATE
	simActual_Multi_MRCA<-c()
	for (j in 1:length(allDat_MULTI_actual[,"MRCA"])){
		simActual_Multi_MRCA[j]<-(allDat_MULTI_actual[j,"MRCA"] / allDat_MULTI_scaled[j,"MRCA"])*Sim_multi_scaled[j,"MRCA"]
	#for (j in 1:length(allDat_MULTI_actual[,"mrcaSQ"])){
	#	simActual_Multi_MRCA[j]<-(allDat_MULTI_actual[j,"mrcaSQ"] / allDat_MULTI_scaled[j,"mrcaSQ"])*Sim_multi_scaled[j,"mrcaSQ"]
	}
	simActual_Multi_DR_harm<-c()
	for (j in 1:length(allDat_MULTI_actual[,"DR_harm"])){
		simActual_Multi_DR_harm[j]<-(allDat_MULTI_actual[j,"DR_harm"] / allDat_MULTI_scaled[j,"DR_harm"])*Sim_multi_scaled[j,"DR_harm"]
	}
	simActual_Multi_DR_skew<-c()
	for (j in 1:length(allDat_MULTI_actual[,"DR_skew"])){
		simActual_Multi_DR_skew[j]<-(allDat_MULTI_actual[j,"DR_skew"] / allDat_MULTI_scaled[j,"DR_skew"])*Sim_multi_scaled[j,"DR_skew"]
	}
res<-cbind.data.frame(Sim_multi_scaled$sliceTimes, simActual_Multi_MRCA,simActual_Multi_DR_harm,simActual_Multi_DR_skew, Sim_multi_scaled$Pval1,Sim_multi_scaled$Pval2,Sim_multi_scaled$Pval3,Sim_multi_scaled$slice,Sim_multi_scaled$tree)
colnames(res)<-c("sliceTimes","MRCA","DR_harm","DR_skew","Pval1","Pval2","Pval3","slice","tree")
#colnames(res)<-c("sliceTimes","mrcaSQ","DR_harm","DR_skew","Pval1","Pval2","Pval3","slice","tree")
simActual_Multis[[i]]<-res
}

###
# FOR the sims scaled - UNIVARIATE

simActual_Unis<-vector("list",length=3)
for (i in 1:3){
Sim_uni_scaled<-SIMS_allDat_UNI_scaled[[i]]
# now do the RATIOS -- as UNIVARIATE
	simActual_Uni_MRCA<-c()
	for (j in 1:length(allDat_UNI_actual[,"MRCA"])){
		simActual_Uni_MRCA[j]<-(allDat_UNI_actual[j,"MRCA"] / allDat_UNI_scaled[j,"MRCA"])*Sim_uni_scaled[j,"MRCA"]
	#for (j in 1:length(allDat_UNI_actual[,"mrcaSQ"])){
	#	simActual_Uni_MRCA[j]<-(allDat_UNI_actual[j,"mrcaSQ"] / allDat_UNI_scaled[j,"mrcaSQ"])*Sim_uni_scaled[j,"mrcaSQ"]
	}
	simActual_Uni_DR_harm<-c()
	for (j in 1:length(allDat_UNI_actual[,"DR_harm"])){
		simActual_Uni_DR_harm[j]<-(allDat_UNI_actual[j,"DR_harm"] / allDat_UNI_scaled[j,"DR_harm"])*Sim_uni_scaled[j,"DR_harm"]
	}
	simActual_Uni_DR_skew<-c()
	for (j in 1:length(allDat_UNI_actual[,"DR_skew"])){
		simActual_Uni_DR_skew[j]<-(allDat_UNI_actual[j,"DR_skew"] / allDat_UNI_scaled[j,"DR_skew"])*Sim_uni_scaled[j,"DR_skew"]
	}
res<-cbind.data.frame(Sim_uni_scaled$sliceTimes, simActual_Uni_MRCA,simActual_Uni_DR_harm,simActual_Uni_DR_skew, Sim_uni_scaled$p_1,Sim_uni_scaled$p_2,Sim_uni_scaled$p_3,Sim_uni_scaled$slice,Sim_uni_scaled$tree)
colnames(res)<-c("sliceTimes","MRCA","DR_harm","DR_skew","p_1","p_2","p_3","slice","tree")
#colnames(res)<-c("sliceTimes","mrcaSQ","DR_harm","DR_skew","p_1","p_2","p_3","slice","tree")
simActual_Unis[[i]]<-res
}

# now get the DISTRIBUTION - HIST - of the coefficients...
allDat_MULTI_actual
simActual_Multis

ColorsSIM_full<-c("darkgoldenrod2","deepskyblue2","darkorchid3")

pdf(file="cladeLevel_MRCA_coefs_HISTOS_wSims.pdf",width=12,height=3)
#quartz(width=12,height=4)
par(mfrow = c(4,14), oma = c(2,2,2,2) + 0.1, mar = c(0,0,0,0) + 0.6)
#layout(matrix(c(1:14), 1, 14, byrow = TRUE), widths=1, heights=rep(1,14))
	for(i in 1:14){
		hist(x=allDat_MULTI_actual[which(allDat_MULTI_actual$slice==i),"MRCA"], xlim=c(-0.05,0.21),breaks=10, main=paste("slice",i,sep=""), col="black", border="black")
	}
for (j in 1:3){
SIM<-simActual_Multis[[j]]
	for(i in 1:14){
		hist(x=SIM[which(SIM$slice==i),"MRCA"], xlim=c(-0.05,0.21),breaks=10, main="", col=ColorsSIM_full[j], xaxt="n",yaxt="n",border=ColorsSIM_full[j])
	}
}
#title(main="Clade age", line = 1,cex.main=1,font.main=1)
dev.off()

pdf(file="cladeLevel_DRmean_coefs_HISTOS_wSims.pdf",width=12,height=3)
#quartz(width=12,height=4)
par(mfrow = c(4,14), oma = c(2,2,2,2) + 0.1, mar = c(0,0,0,0) + 0.6)
#layout(matrix(c(1:14), 1, 14, byrow = TRUE), widths=1, heights=rep(1,14))
	for(i in 1:14){
		hist(x=allDat_MULTI_actual[which(allDat_MULTI_actual$slice==i),"DR_harm"], xlim=c(0,24),breaks=10, main=paste("slice",i,sep=""), col="black", border="black")
	}
for (j in 1:3){
SIM<-simActual_Multis[[j]]
	for(i in 1:14){
		hist(x=SIM[which(SIM$slice==i),"DR_harm"], xlim=c(0,24),breaks=10, main="", col=ColorsSIM_full[j], xaxt="n",yaxt="n",border=ColorsSIM_full[j])
	}
}
dev.off()

pdf(file="cladeLevel_DRskew_coefs_HISTOS_wSims.pdf",width=12,height=3)
#quartz(width=12,height=4)
par(mfrow = c(4,14), oma = c(2,2,2,2) + 0.1, mar = c(0,0,0,0) + 0.6)
#layout(matrix(c(1:14), 1, 14, byrow = TRUE), widths=1, heights=rep(1,14))
	for(i in 1:14){
		hist(x=allDat_MULTI_actual[which(allDat_MULTI_actual$slice==i),"DR_skew"], xlim=c(-0.05,2.2),breaks=10, main=paste("slice",i,sep=""), col="black", border="black")
	}
for (j in 1:3){
SIM<-simActual_Multis[[j]]
	for(i in 1:14){
		hist(x=SIM[which(SIM$slice==i),"DR_skew"], xlim=c(-0.05,2.2),breaks=10, main="", col=ColorsSIM_full[j], xaxt="n",yaxt="n",border=ColorsSIM_full[j])
	}
}
dev.off()


library(plotrix)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")
# plot through time
#pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_timeSlices_partSlopes_PLOTTED_5part-wPerSamp_SCALED_withSIMS-values.pdf",sep=""),onefile=TRUE, width=4,height=12)
#pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_BROWNIAN_timeSlices_partSlopes_PLOTTED_5part-wPerSamp_SCALED_withSIMS-values_andSE-BAR.pdf",sep=""),onefile=TRUE, width=4,height=12)
#pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_BROWNIAN_timeSlices_partSlopes_PLOTTED_4part_SCALED_withSIMS-values_andSE-BAR.pdf",sep=""),onefile=TRUE, width=4,height=10)
pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_noTREE_timeSlices_partSlopes_PLOTTED_4part_SCALED_withSIMS-values_andSE-BAR_GLS.pdf",sep=""),onefile=TRUE, width=4,height=10)
#pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_timeSlices_partSlopes_PLOTTED_3part_ACTUAL_withSIMS-values.pdf",sep=""),onefile=TRUE, width=4,height=9)
#pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_timeSlices_partSlopes_PLOTTED_3part_mrcaSQ_RAW_withSIMS-values.pdf",sep=""),onefile=TRUE, width=4,height=9)
#pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_timeSlices_partSlopes_PLOTTED_3part_mrcaSQ_SCALED_withSIMS-values.pdf",sep=""),onefile=TRUE, width=4,height=9)

#dat<-allDat_MULTI_actual
#dat<-allDat_MULTI_scaled
#dat<-allDat_MULTI_scaled_4var
dat<-allDat_MULTI_scaled_gls
#Pval4<- dat$Pval5
#lowHigh95_4<-cbind.data.frame(dat$percentSamp-(1.96*dat$SE5),dat$percentSamp,dat$percentSamp+(1.96*dat$SE5),rep(sliceTimes,100))
#colnames(lowHigh95_4)<-c("low","mean","high","sliceTimes")

Pval1<- dat$Pval2
Pval2<- dat$Pval3
Pval3<- dat$Pval4
adjLab<-0.08
redAlpha<-0.3
greyAlpha<-0.1
vertLwd<-0.5
lowHigh95_1<-cbind.data.frame(dat$MRCA-(1.96*dat$SE2),dat$MRCA,dat$MRCA+(1.96*dat$SE2),rep(sliceTimes,100))
colnames(lowHigh95_1)<-c("low","mean","high","sliceTimes")
lowHigh95_2<-cbind.data.frame(dat$DR_harm-(1.96*dat$SE3),dat$DR_harm,dat$DR_harm+(1.96*dat$SE3),rep(sliceTimes,100))
colnames(lowHigh95_2)<-c("low","mean","high","sliceTimes")
lowHigh95_3<-cbind.data.frame(dat$DR_skew-(1.96*dat$SE4),dat$DR_skew,dat$DR_skew+(1.96*dat$SE4),rep(sliceTimes,100))
colnames(lowHigh95_3)<-c("low","mean","high","sliceTimes")
lwdCI<-1.5

Header<-"MULTIVARIATE - PGLS per time slice per tree"
#Sub<-"Actual data values"
#YlabMain<-"Actual unique effect on (log) clade richness"
#horizLine<-0
Sub<-"Standardized data values"
YlabMain<-"Standardized unique effect on (log) clade richness"
horizLine<-0

#SIMS<-simActual_Multis
#yLims1<-c(-0.01,0.24)
#yLims2<-c(0,35)
#yLims3<-c(0,3.5)
#yAxis<-NULL
SIMS<-SIMS_allDat_MULTI_scaled
yLims1<-c(-0.2,2.1)
yLims2<-c(-0.2,2.1)
yLims3<-c(-0.2,2.1)
yLims4<-c(-1,1)
yAxis<-c(NA,0,0.5,1,1.5,2,NA)
simAlpha<-0.05

#######
##
#pdf(file=paste("cladeLevel",bbone,"_PGLS-UNI_timeSlices_uniSlopes_PLOTTED_3part_SCALED_withSIMS-values.pdf",sep=""),onefile=TRUE, width=4,height=9)
#pdf(file=paste("cladeLevel",bbone,"_PGLS-UNI_timeSlices_uniSlopes_PLOTTED_3part_ACTUAL_withSIMS-values.pdf",sep=""),onefile=TRUE, width=4,height=9)
#pdf(file=paste("cladeLevel",bbone,"_PGLS-UNI_timeSlices_uniSlopes_PLOTTED_3part_mrcaSQ_SCALED_withSIMS-values.pdf",sep=""),onefile=TRUE, width=4,height=9)
#pdf(file=paste("cladeLevel",bbone,"_PGLS-UNI_timeSlices_uniSlopes_PLOTTED_3part_mrcaSQ_RAW_withSIMS-values.pdf",sep=""),onefile=TRUE, width=4,height=9)
#dat<-allDat_UNI_scaled
#dat<-allDat_UNI_actual
#Pval1<- dat$p_1
#Pval2<- dat$p_2
#Pval3<- dat$p_3
#adjLab<-0.08 
#Header<-"UNIVARIATE - PGLS per time slice per tree" 

#Sub<-"Actual data values"
#YlabMain<-"Actual slopes vs. (log) clade richness"
#Sub<-"Standardized data values"
#YlabMain<-"Standardized slopes vs. (log) clade richness"
#
#SIMS<-simActual_Unis
#yAxis<-NULL
#yLims1<-c(-0.01,0.24)
#yLims2<-c(0,35)
#yLims3<-c(0,3.5)
#SIMS<-SIMS_allDat_UNI_scaled
#yLims1<-c(-0.2,2.5)
#yLims2<-c(-0.2,2.5)
#yLims3<-c(-0.2,2.5)
#yAxis<-c(NA,0,0.5,1,1.5,2,NA)

#quartz(width=4,height=9)
#layout(matrix(c(1:4), 4, 1, byrow = TRUE), widths=4, heights=c(4,2,2,2))
#par(oma = c(5,4,5,0) + 0.1, mar = c(5,1,1,1) + 0.1)

#layout(matrix(c(1:5), 5, 1, byrow = TRUE), widths=4, heights=c(4,2.5,2.5,2.5,2.5))#c(3,3,3))
layout(matrix(c(1:4), 4, 1, byrow = TRUE), widths=4, heights=c(4,2.5,2.5,2.5))#c(3,3,3))
par(oma = c(5,4,5,3) + 0.1, mar = rep(0.5,4) + 0.1)
#par(mfrow = c(4,1), oma = c(5,4,5,0) + 0.1, mar = c(5,1,1,1) + 0.1)
#par(op)

#part A
# for NDexp MCC
plot(mamPhy, edge.width=0.45, cex=0.1, show.tip.label=FALSE, no.margin=FALSE, x.lim=c(110.8,175.9)) #, tip.color = tcol )
for(i in 1:numSlices){
	abline(v=rootTime-(5*i),lty=2, lwd=vertLwd, col="grey")
}
mtext(side=3,text=Sub)

# part B
##
Col<-"darkgoldenrod2"
Col1<-rgb((col2rgb(Col)/255)[1],(col2rgb(Col)/255)[2],(col2rgb(Col)/255)[3],alpha=simAlpha)
Col<-"deepskyblue2"
Col2<-rgb((col2rgb(Col)/255)[1],(col2rgb(Col)/255)[2],(col2rgb(Col)/255)[3],alpha=simAlpha)
Col<-"darkorchid3"
Col3<-rgb((col2rgb(Col)/255)[1],(col2rgb(Col)/255)[2],(col2rgb(Col)/255)[3],alpha=simAlpha)
ColorsSIM<-c(Col1,Col2,Col3)
ColorsSIM_full<-c("darkgoldenrod2","deepskyblue2","darkorchid3")

#plot(MRCA ~ sliceTimes, data=dat, type="n",col=cols, pch=pchs, ylim=yLims, ylab="",xlab="", yaxt="n",xaxt="n")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
plot(MRCA ~ sliceTimes, data=dat, type="n",col=cols, pch=pchs, ylim=yLims1, ylab="",xlab="", yaxt="n",xaxt="n")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
axis(side=2,at=yAxis,labels=TRUE)
axis(side=1,at=NULL,labels=FALSE)

mtext(text="(b) Crown age", padj=2, adj=adjLab, font=2)#at=c(-72,0.21)
#mtext(text="(a) Crown age (squared)", padj=2, adj=adjLab, font=2)#at=c(-72,0.21)
#mtext(side=3,text=Sub)

cols<-rep(grey(0.3,alpha=greyAlpha),length(dat[,1]))
cols[which(Pval1>0.05)]<-rgb(1,0,0,alpha=redAlpha)
pchs<-rep(1,length(dat[,1]))
pchs[which(Pval1>0.05)]<-17
#sfracs<-rep(0,1400)
#sfracs[which(lowHigh95_1$high==max(lowHigh95_1$high))]<-3
#cols[which(lowHigh95_1$high==max(lowHigh95_1$high))]<-grey(0.3)

plotCI(x=lowHigh95_1$sliceTimes, add=TRUE, y=lowHigh95_1$mean,ui=lowHigh95_1$high,li=lowHigh95_1$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=cols,scol=grey(0.3,alpha=greyAlpha),pch=pchs,font.lab=2,cex.axis=1.1,cex.lab=1.1)

#points(MRCA ~ sliceTimes, data=dat, col=cols, pch=pchs, ylab="",xlab="", yaxt="n",xaxt="n")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
#	points(low ~ sliceTimes, data=lowHigh95_1, col=cols, pch=pchs, ylab="",xlab="", yaxt="n",xaxt="n")
#	points(high ~ sliceTimes, data=lowHigh95_1, col=cols, pch=pchs, ylab="",xlab="", yaxt="n",xaxt="n")

	m1<-lme(fixed = MRCA ~ sliceTimes, random = ~ 1 | tree/slice, data=dat)
	sum1<-summary(m1)
	m2<-lme(fixed = MRCA ~ sliceTimes + I(sliceTimes^2), random = ~ 1 | tree/slice, data=dat)
	sum2<-summary(m2)
	m3<-lme(fixed = MRCA ~ sliceTimes + I(sliceTimes^2) + I(sliceTimes^3), random = ~ 1 | tree/slice, data=dat)
	sum3<-summary(m3)
	models<-list(m1,m2,m3)
	AICs<-sapply(models,AICc)

	m<-models[which(min(AICs)==AICs)]
	sum<-summary(m[[1]])
	fixed<-sum$coef$fixed
	#ci95<-intervals(m[[1]],level=0.95,which="fixed")
	#lower<-ci95[[1]][,1]
	#upper<-ci95[[1]][,3]

	if (length(fixed)==4){
	#curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2+ fixed[[4]]*x^3, add=TRUE, lwd=3, lty=1, col="red")
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2+ fixed[[4]]*x^3, add=TRUE, lwd=3, lty=1, col="black")
	#curve(lower[[1]] + lower[[2]]*x + lower[[3]]*x^2+ lower[[4]]*x^3, add=TRUE, lwd=3, lty=1, col="black")
	#curve(upper[[1]] + upper[[2]]*x + upper[[3]]*x^2+ upper[[4]]*x^3, add=TRUE, lwd=3, lty=2, col="red")
	}
	
	if (length(fixed)==3){
	#curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2, add=TRUE, lwd=3, lty=1, col="red")
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2, add=TRUE, lwd=3, lty=1, col="black")
	}

	if (length(fixed)==2){
	#curve(fixed[[1]] + fixed[[2]]*x, add=TRUE, lwd=3, lty=1, col="red")
	curve(fixed[[1]] + fixed[[2]]*x, add=TRUE, lwd=3, lty=1, col="black")
	}

for (j in 1:length(SIMS)){
	datSIM<-SIMS[[j]]
	p1<- datSIM$Pval1

	cols<-rep(ColorsSIM[j],length(datSIM[,1]))
	cols[which(p1>0.05)]<-ColorsSIM_full[j]
	pchs<-rep(1,length(datSIM[,1]))
	pchs[which(p1>0.05)]<-17

	lowHigh95_s1<-cbind.data.frame(datSIM$MRCA-(1.96*datSIM$SE2),datSIM$MRCA,datSIM$MRCA+(1.96*datSIM$SE2),rep(sliceTimes,100))
	colnames(lowHigh95_s1)<-c("low","mean","high","sliceTimes")

	plotCI(x=lowHigh95_s1$sliceTimes, add=TRUE, y=lowHigh95_s1$mean,ui=lowHigh95_s1$high,li=lowHigh95_s1$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=cols,scol=cols,pch=pchs,font.lab=2,cex.axis=1.1,cex.lab=1.1)

	#points(MRCA ~ sliceTimes, data=datSIM, col=cols, pch=pchs, ylab="",xlab="", yaxt="n",xaxt="n")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
	#points(low ~ sliceTimes, data=lowHigh95_s1, col=cols, pch=pchs, ylab="",xlab="", yaxt="n",xaxt="n")
	#points(high ~ sliceTimes, data=lowHigh95_s1, col=cols, pch=pchs, ylab="",xlab="", yaxt="n",xaxt="n")

	m1<-lme(fixed = MRCA ~ sliceTimes, random = ~ 1 | tree/slice, data=datSIM)
	sum1<-summary(m1)
		for (i in 1:10) {possibleError <- tryCatch(
		      m2<-lme(fixed = MRCA ~ sliceTimes + I(sliceTimes^2), random = ~ 1 | tree/slice, data=datSIM),
		      error=function(e) e)
		if(inherits(possibleError, "lme")) break		
		if(inherits(possibleError, "error")) next}
		m2<-possibleError
	sum2<-summary(m2)
		for (i in 1:10) {possibleError <- tryCatch(
		      m3<-lme(fixed = MRCA ~ sliceTimes + I(sliceTimes^2) + I(sliceTimes^3), random = ~ 1 | tree/slice, data=datSIM),
		      error=function(e) e)
		if(inherits(possibleError, "lme")) break		
		if(inherits(possibleError, "error")) next}
		m3<-possibleError
	sum3<-summary(m3)
	models<-list(m1,m2,m3)
	AICs<-sapply(models,AICc)

	m<-models[which(min(AICs)==AICs)]
	sum<-summary(m[[1]])
	fixed<-sum$coef$fixed
	ci95<-intervals(m[[1]],level=0.95,which="fixed")
	lower<-ci95[[1]][,1]
	upper<-ci95[[1]][,3]

	if (length(fixed)==4){
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2+ fixed[[4]]*x^3, add=TRUE, lwd=2, lty=2, col=ColorsSIM_full[j])
	#curve(lower[[1]] + lower[[2]]*x + lower[[3]]*x^2+ lower[[4]]*x^3, add=TRUE, lwd=2, lty=2, col=ColorsSIM[j])
	#curve(upper[[1]] + upper[[2]]*x + upper[[3]]*x^2+ upper[[4]]*x^3, add=TRUE, lwd=2, lty=2, col=ColorsSIM[j])
	}
	
	if (length(fixed)==3){
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2, add=TRUE, lwd=2, lty=2, col=ColorsSIM_full[j])
	}

	if (length(fixed)==2){
	curve(fixed[[1]] + fixed[[2]]*x, add=TRUE, lwd=2, lty=2, col=ColorsSIM_full[j])
	}

}

for(i in 1:numSlices){
	abline(v=(-5*i),lty=2, lwd=vertLwd, col=grey(0.6, alpha=0.5))
}
	abline(h=horizLine,lty=1, lwd=1, col=grey(0.6, alpha=0.5))


# part C
plot(DR_harm ~ sliceTimes, data=dat, type="n",col=cols, pch=pchs, ylim=yLims2, ylab="",xlab="", yaxt="n",xaxt="n")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
axis(side=2,at=yAxis,labels=TRUE)
axis(side=1,at=NULL,labels=FALSE)
mtext(text="(c) Tip DR mean", padj=2, adj=adjLab, font=2)#at=c(-72,0.21)

cols<-rep(grey(0.3,alpha=greyAlpha),length(dat[,1]))
cols[which(Pval2>0.05)]<-rgb(1,0,0,alpha=redAlpha)
pchs<-rep(1,length(dat[,1]))
pchs[which(Pval2>0.05)]<-17

plotCI(x=lowHigh95_2$sliceTimes, add=TRUE, y=lowHigh95_2$mean,ui=lowHigh95_2$high,li=lowHigh95_2$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=cols,scol=grey(0.3,alpha=greyAlpha),pch=pchs,font.lab=2,cex.axis=1.1,cex.lab=1.1)

#points(DR_harm ~ sliceTimes, data=dat, col=cols, pch=pchs, ylab="",xlab="", yaxt="n",xaxt="n")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
#	points(low ~ sliceTimes, data=lowHigh95_2, col=cols, pch=pchs, ylab="",xlab="", yaxt="n",xaxt="n")
#	points(high ~ sliceTimes, data=lowHigh95_2, col=cols, pch=pchs, ylab="",xlab="", yaxt="n",xaxt="n")

	m1<-lme(fixed = DR_harm ~ sliceTimes, random = ~ 1 | tree/slice, data=dat)
	sum1<-summary(m1)
	m2<-lme(fixed = DR_harm ~ sliceTimes + I(sliceTimes^2), random = ~ 1 | tree/slice, data=dat)
	sum2<-summary(m2)
	m3<-lme(fixed = DR_harm ~ sliceTimes + I(sliceTimes^2) + I(sliceTimes^3), random = ~ 1 | tree/slice, data=dat)
	sum3<-summary(m3)
	models<-list(m1,m2,m3)
	AICs<-sapply(models,AICc)

	m<-models[which(min(AICs)==AICs)]
	sum<-summary(m[[1]])
	fixed<-sum$coef$fixed
	#ci95<-intervals(m[[1]],level=0.95,which="fixed")
	#lower<-ci95[[1]][,1]
	#upper<-ci95[[1]][,3]

	if (length(fixed)==4){
	#curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2+ fixed[[4]]*x^3, add=TRUE, lwd=3, lty=1, col="red")
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2+ fixed[[4]]*x^3, add=TRUE, lwd=3, lty=1, col="black")
	#curve(lower[[1]] + lower[[2]]*x + lower[[3]]*x^2+ lower[[4]]*x^3, add=TRUE, lwd=3, lty=1, col="black")
	#curve(upper[[1]] + upper[[2]]*x + upper[[3]]*x^2+ upper[[4]]*x^3, add=TRUE, lwd=3, lty=2, col="red")
	}
	
	if (length(fixed)==3){
	#curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2, add=TRUE, lwd=3, lty=1, col="red")
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2, add=TRUE, lwd=3, lty=1, col="black")
	}

	if (length(fixed)==2){
	#curve(fixed[[1]] + fixed[[2]]*x, add=TRUE, lwd=3, lty=1, col="red")
	curve(fixed[[1]] + fixed[[2]]*x, add=TRUE, lwd=3, lty=1, col="black")
	}

for (j in 1:length(SIMS)){
	datSIM<-SIMS[[j]]
	p2<- datSIM$Pval2

	cols<-rep(ColorsSIM[j],length(datSIM[,1]))
	cols[which(p2>0.05)]<-ColorsSIM_full[j]
	pchs<-rep(1,length(datSIM[,1]))
	pchs[which(p2>0.05)]<-17

	lowHigh95_s2<-cbind.data.frame(datSIM$DR_harm-(1.96*datSIM$SE3),datSIM$DR_harm,datSIM$DR_harm+(1.96*datSIM$SE3),rep(sliceTimes,100))
	colnames(lowHigh95_s2)<-c("low","mean","high","sliceTimes")

	plotCI(x=lowHigh95_s2$sliceTimes, add=TRUE, y=lowHigh95_s2$mean,ui=lowHigh95_s2$high,li=lowHigh95_s2$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=cols,scol=cols,pch=pchs,font.lab=2,cex.axis=1.1,cex.lab=1.1)

	#points(DR_harm ~ sliceTimes, data=datSIM, col=cols, pch=pchs, ylab="",xlab="", yaxt="n",xaxt="n")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
	#points(low ~ sliceTimes, data=lowHigh95_s2, col=cols, pch=pchs, ylab="",xlab="", yaxt="n",xaxt="n")
	#points(high ~ sliceTimes, data=lowHigh95_s2, col=cols, pch=pchs, ylab="",xlab="", yaxt="n",xaxt="n")

	m1<-lme(fixed = DR_harm ~ sliceTimes, random = ~ 1 | tree/slice, data=datSIM)
	sum1<-summary(m1)
		for (i in 1:10) {possibleError <- tryCatch(
		      m2<-lme(fixed = DR_harm ~ sliceTimes + I(sliceTimes^2), random = ~ 1 | tree/slice, data=datSIM),
		      error=function(e) e)
		if(inherits(possibleError, "lme")) break		
		if(inherits(possibleError, "error")) next}
		m2<-possibleError
	sum2<-summary(m2)
		for (i in 1:10) {possibleError <- tryCatch(
		      m3<-lme(fixed = DR_harm ~ sliceTimes + I(sliceTimes^2) + I(sliceTimes^3), random = ~ 1 | tree/slice, data=datSIM),
		      error=function(e) e)
		if(inherits(possibleError, "lme")) break		
		if(inherits(possibleError, "error")) next}
		m3<-possibleError
	sum3<-summary(m3)
	models<-list(m1,m2,m3)
	AICs<-sapply(models,AICc)

	m<-models[which(min(AICs)==AICs)]
	sum<-summary(m[[1]])
	fixed<-sum$coef$fixed
	ci95<-intervals(m[[1]],level=0.95,which="fixed")
	lower<-ci95[[1]][,1]
	upper<-ci95[[1]][,3]

	if (length(fixed)==4){
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2+ fixed[[4]]*x^3, add=TRUE, lwd=2, lty=2, col=ColorsSIM_full[j])
	#curve(lower[[1]] + lower[[2]]*x + lower[[3]]*x^2+ lower[[4]]*x^3, add=TRUE, lwd=2, lty=2, col=ColorsSIM[j])
	#curve(upper[[1]] + upper[[2]]*x + upper[[3]]*x^2+ upper[[4]]*x^3, add=TRUE, lwd=2, lty=2, col=ColorsSIM[j])
	}
	
	if (length(fixed)==3){
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2, add=TRUE, lwd=2, lty=2, col=ColorsSIM_full[j])
	}

	if (length(fixed)==2){
	curve(fixed[[1]] + fixed[[2]]*x, add=TRUE, lwd=2, lty=2, col=ColorsSIM_full[j])
	}

}

for(i in 1:numSlices){
	abline(v=(-5*i),lty=2, lwd=vertLwd, col=grey(0.6, alpha=0.5))
}
	abline(h=horizLine,lty=1, lwd=1, col=grey(0.6, alpha=0.5))


# part D
plot(DR_skew ~ sliceTimes, data=dat, type="n",col=cols, pch=pchs, ylim=yLims3, ylab="",xlab="", yaxt="n")#,xaxt="n")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
axis(side=2,at=yAxis,labels=TRUE)
mtext(text="(d) Tip DR skewness", padj=2, adj=adjLab, font=2)#at=c(-72,0.21)

cols<-rep(grey(0.3,alpha=greyAlpha),length(dat[,1]))
cols[which(Pval3>0.05)]<-rgb(1,0,0,alpha=redAlpha)
pchs<-rep(1,length(dat[,1]))
pchs[which(Pval3>0.05)]<-17

plotCI(x=lowHigh95_3$sliceTimes, add=TRUE, y=lowHigh95_3$mean,ui=lowHigh95_3$high,li=lowHigh95_3$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=cols,scol=grey(0.3,alpha=greyAlpha),pch=pchs,font.lab=2,cex.axis=1.1,cex.lab=1.1)

#points(DR_skew ~ sliceTimes, data=dat, col=cols, pch=pchs, ylab="",xlab="", yaxt="n",xaxt="n")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
#	points(low ~ sliceTimes, data=lowHigh95_3, col=cols, pch=pchs, ylab="",xlab="", yaxt="n",xaxt="n")
#	points(high ~ sliceTimes, data=lowHigh95_3, col=cols, pch=pchs, ylab="",xlab="", yaxt="n",xaxt="n")

	m1<-lme(fixed = DR_skew ~ sliceTimes, random = ~ 1 | tree/slice, data=dat)
	sum1<-summary(m1)
	m2<-lme(fixed = DR_skew ~ sliceTimes + I(sliceTimes^2), random = ~ 1 | tree/slice, data=dat)
	sum2<-summary(m2)
	m3<-lme(fixed = DR_skew ~ sliceTimes + I(sliceTimes^2) + I(sliceTimes^3), random = ~ 1 | tree/slice, data=dat)
	sum3<-summary(m3)
	models<-list(m1,m2,m3)
	AICs<-sapply(models,AICc)

	m<-models[which(min(AICs)==AICs)]
	sum<-summary(m[[1]])
	fixed<-sum$coef$fixed
	#ci95<-intervals(m[[1]],level=0.95,which="fixed")
	#lower<-ci95[[1]][,1]
	#upper<-ci95[[1]][,3]

	if (length(fixed)==4){
	#curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2+ fixed[[4]]*x^3, add=TRUE, lwd=3, lty=1, col="red")
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2+ fixed[[4]]*x^3, add=TRUE, lwd=3, lty=1, col="black")
	#curve(lower[[1]] + lower[[2]]*x + lower[[3]]*x^2+ lower[[4]]*x^3, add=TRUE, lwd=3, lty=1, col="black")
	#curve(upper[[1]] + upper[[2]]*x + upper[[3]]*x^2+ upper[[4]]*x^3, add=TRUE, lwd=3, lty=2, col="red")
	}
	
	if (length(fixed)==3){
	#curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2, add=TRUE, lwd=3, lty=1, col="red")
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2, add=TRUE, lwd=3, lty=1, col="black")
	}

	if (length(fixed)==2){
	#curve(fixed[[1]] + fixed[[2]]*x, add=TRUE, lwd=3, lty=1, col="red")
	curve(fixed[[1]] + fixed[[2]]*x, add=TRUE, lwd=3, lty=1, col="black")
	}

for (j in 1:length(SIMS)){
	datSIM<-SIMS[[j]]
	p3<- datSIM$Pval3

	cols<-rep(ColorsSIM[j],length(datSIM[,1]))
	cols[which(p3>0.05)]<-ColorsSIM_full[j]
	pchs<-rep(1,length(datSIM[,1]))
	pchs[which(p3>0.05)]<-17

	lowHigh95_s3<-cbind.data.frame(datSIM$DR_skew-(1.96*datSIM$SE4),datSIM$DR_skew,datSIM$DR_skew+(1.96*datSIM$SE4),rep(sliceTimes,100))
	colnames(lowHigh95_s3)<-c("low","mean","high","sliceTimes")

	plotCI(x=lowHigh95_s3$sliceTimes, add=TRUE, y=lowHigh95_s3$mean,ui=lowHigh95_s3$high,li=lowHigh95_s3$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=cols,scol=cols,pch=pchs,font.lab=2,cex.axis=1.1,cex.lab=1.1)

	#points(DR_skew ~ sliceTimes, data=datSIM, col=cols, pch=pchs, ylab="",xlab="", yaxt="n",xaxt="n")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
	#points(low ~ sliceTimes, data=lowHigh95_s3, col=cols, pch=pchs, ylab="",xlab="", yaxt="n",xaxt="n")
	#points(high ~ sliceTimes, data=lowHigh95_s3, col=cols, pch=pchs, ylab="",xlab="", yaxt="n",xaxt="n")

	m1<-lme(fixed = DR_skew ~ sliceTimes, random = ~ 1 | tree/slice, data=datSIM)
	sum1<-summary(m1)
		for (i in 1:10) {possibleError <- tryCatch(
		      m2<-lme(fixed = DR_skew ~ sliceTimes + I(sliceTimes^2), random = ~ 1 | tree/slice, data=datSIM),
		      error=function(e) e)
		if(inherits(possibleError, "lme")) break		
		if(inherits(possibleError, "error")) next}
		m2<-possibleError
	sum2<-summary(m2)
		for (i in 1:10) {possibleError <- tryCatch(
		      m3<-lme(fixed = DR_skew ~ sliceTimes + I(sliceTimes^2) + I(sliceTimes^3), random = ~ 1 | tree/slice, data=datSIM),
		      error=function(e) e)
		if(inherits(possibleError, "lme")) break		
		if(inherits(possibleError, "error")) next}
		m3<-possibleError
	sum3<-summary(m3)
	models<-list(m1,m2,m3)
	AICs<-sapply(models,AICc)

	m<-models[which(min(AICs)==AICs)]
	sum<-summary(m[[1]])
	fixed<-sum$coef$fixed
	ci95<-intervals(m[[1]],level=0.95,which="fixed")
	lower<-ci95[[1]][,1]
	upper<-ci95[[1]][,3]

	if (length(fixed)==4){
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2+ fixed[[4]]*x^3, add=TRUE, lwd=2, lty=2, col=ColorsSIM_full[j])
	#curve(lower[[1]] + lower[[2]]*x + lower[[3]]*x^2+ lower[[4]]*x^3, add=TRUE, lwd=2, lty=2, col=ColorsSIM[j])
	#curve(upper[[1]] + upper[[2]]*x + upper[[3]]*x^2+ upper[[4]]*x^3, add=TRUE, lwd=2, lty=2, col=ColorsSIM[j])
	}
	
	if (length(fixed)==3){
	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2, add=TRUE, lwd=2, lty=2, col=ColorsSIM_full[j])
	}

	if (length(fixed)==2){
	curve(fixed[[1]] + fixed[[2]]*x, add=TRUE, lwd=2, lty=2, col=ColorsSIM_full[j])
	}

}

for(i in 1:numSlices){
	abline(v=(-5*i),lty=2, lwd=vertLwd, col=grey(0.6, alpha=0.5))
}
	abline(h=horizLine,lty=1, lwd=1, col=grey(0.6, alpha=0.5))

## part E
#plot(percentSamp ~ sliceTimes, data=dat, type="n",col=cols, pch=pchs, ylim=yLims4, ylab="",xlab="", yaxt="n")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
#axis(side=2,at=c(-0.5,0,0.5),labels=TRUE)
#mtext(text="(e) Percent sampled tips", padj=2, adj=adjLab, font=2)#at=c(-72,0.21)
#
#cols<-rep(grey(0.3,alpha=greyAlpha),length(dat[,1]))
#cols[which(Pval4>0.05)]<-rgb(1,0,0,alpha=redAlpha)
#pchs<-rep(1,length(dat[,1]))
#pchs[which(Pval4>0.05)]<-17
#
#plotCI(x=lowHigh95_4$sliceTimes, add=TRUE, y=lowHigh95_4$mean,ui=lowHigh95_4$high,li=lowHigh95_4$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=cols,scol=grey(0.3,alpha=greyAlpha),pch=pchs,font.lab=2,cex.axis=1.1,cex.lab=1.1)
##points(percentSamp ~ sliceTimes, data=dat, col=cols, pch=pchs, ylim=yLims, ylab="",xlab="", yaxt="n",xaxt="n")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
#
#	m1<-lme(fixed = percentSamp ~ sliceTimes, random = ~ 1 | tree/slice, data=dat)
#	sum1<-summary(m1)
#	m2<-lme(fixed = percentSamp ~ sliceTimes + I(sliceTimes^2), random = ~ 1 | tree/slice, data=dat)
#	sum2<-summary(m2)
#	m3<-lme(fixed = percentSamp ~ sliceTimes + I(sliceTimes^2) + I(sliceTimes^3), random = ~ 1 | tree/slice, data=dat)
#	sum3<-summary(m3)
#	models<-list(m1,m2,m3)
#	AICs<-sapply(models,AICc)
#
#	m<-models[which(min(AICs)==AICs)]
#	sum<-summary(m[[1]])
#	fixed<-sum$coef$fixed
#	#ci95<-intervals(m[[1]],level=0.95,which="fixed")
#	#lower<-ci95[[1]][,1]
#	#upper<-ci95[[1]][,3]
#
#	if (length(fixed)==4){
#	#curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2+ fixed[[4]]*x^3, add=TRUE, lwd=3, lty=1, col="red")
#	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2+ fixed[[4]]*x^3, add=TRUE, lwd=3, lty=1, col="black")
#	#curve(lower[[1]] + lower[[2]]*x + lower[[3]]*x^2+ lower[[4]]*x^3, add=TRUE, lwd=3, lty=1, col="black")
#	#curve(upper[[1]] + upper[[2]]*x + upper[[3]]*x^2+ upper[[4]]*x^3, add=TRUE, lwd=3, lty=2, col="red")
#	}
#	
#	if (length(fixed)==3){
#	#curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2, add=TRUE, lwd=3, lty=1, col="red")
#	curve(fixed[[1]] + fixed[[2]]*x + fixed[[3]]*x^2, add=TRUE, lwd=3, lty=1, col="black")
#	}
#
#	if (length(fixed)==2){
#	#curve(fixed[[1]] + fixed[[2]]*x, add=TRUE, lwd=3, lty=1, col="red")
#	curve(fixed[[1]] + fixed[[2]]*x, add=TRUE, lwd=3, lty=1, col="black")
#	}
#
#for(i in 1:numSlices){
#	abline(v=(-5*i),lty=2, lwd=vertLwd, col=grey(0.6, alpha=0.5))
#}
#	abline(h=horizLine,lty=1, lwd=1, col=grey(0.6, alpha=0.5))


title(main=Header, xlab = "Time slices before present (Ma)",
      ylab = YlabMain,
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)

dev.off()

#####
# Compare to PAGEL.

setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_PAGELruns")

allFinTotal<-c(1, 1, 1, 3, 3, 3, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 16, 16, 16, 17, 17, 17, 18, 18, 18, 19, 19, 19, 21, 21, 21, 22, 22, 22, 23, 23, 23, 24, 24, 24, 25, 25, 25, 27, 27, 27, 29, 29, 29, 30, 30, 30, 32, 32, 32, 33, 33, 33, 34, 34, 34, 35, 35, 35, 36, 36, 36, 37, 37, 37, 38, 38, 38, 39, 39, 39, 40, 40, 40, 41, 41, 41, 42, 42, 42, 44, 44, 44, 45, 45, 45, 46, 46, 46, 47, 47, 47, 48, 48, 48, 49, 49, 49, 50, 50, 50, 51, 51, 51, 52, 52, 52, 53, 53, 53, 54, 54, 54, 57, 57, 57, 59, 59, 59, 60, 60, 60, 61, 61, 61, 62, 62, 62, 63, 63, 63, 64, 64, 64, 65, 65, 65, 66, 66, 66, 67, 67, 67, 68, 68, 68, 70, 70, 70, 71, 71, 71, 72, 72, 72, 73, 73, 73, 74, 74, 74, 75, 75, 75, 77, 77, 77, 78, 78, 78, 79, 79, 79, 80, 80, 80, 81, 81, 81, 82, 82, 82, 83, 83, 83, 84, 84, 84, 85, 85, 85, 86, 86, 86, 88, 88, 88, 90, 90, 90, 91, 91, 91, 92, 92, 92, 93, 93, 93, 94, 94, 94, 95, 95, 95, 96, 96, 96, 97, 97, 97, 98, 98, 98, 100, 100, 100)
allFin<-unique(allFinTotal)

allSlopes_per_i<-vector("list",length=numTrees)
slice<-c(1:14)
for(i in allFin){
	tree<-rep(i,length(slice))
	dat<-read.table(paste(bbone,"_sample100_",i,"_PGLSmulti_PAGEL_noLamCI_timeSlices_SCALED_partSlopes123.txt",sep=""))
	#dat<-read.table(paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlices_mrcaSQ_RAW_partSlopes123.txt",sep=""))
	allSlopes_per_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
Pagel_MULTI_scaled<-do.call(rbind,allSlopes_per_i)
write.table(Pagel_MULTI_scaled,file=paste(bbone,"_sample100_Just83trees-TABLE_PGLSmulti_PAGEL_timeSlices_SCALED_3vars.txt",sep=""))

allSlopes_per_i<-vector("list",length=numTrees)
slice<-c(1:14)
for(i in allFin){
	tree<-rep(i,length(slice))
	dat<-read.table(paste(bbone,"_sample100_",i,"_PGLS_PAGEL_noLamCI_timeSlices_wPercentSamp_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))
	#dat<-read.table(paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlices_mrcaSQ_RAW_partSlopes123.txt",sep=""))
	allSlopes_per_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
Pagel_UNI_scaled<-do.call(rbind,allSlopes_per_i)
write.table(Pagel_UNI_scaled,file=paste(bbone,"_sample100_Just83trees-TABLE_PGLS-UNI_PAGEL_timeSlices_SCALED.txt",sep=""))


allSlopes_per_i<-vector("list",length=numTrees)
simFin<-c(7,11,33)
for(i in simFin){
	tree<-rep(i,length(slice))
	dat<-read.table(paste("MamPhy_SIMS_mamPhyE_",bbone,"_sample100_",i,"_PGLSmulti_PAGEL_noLamCI_timeSlices_SCALED_partSlopes123.txt",sep=""))
	#dat<-read.table(paste(bbone,"_sample100_",i,"_PGLSmulti_timeSlices_mrcaSQ_RAW_partSlopes123.txt",sep=""))
	allSlopes_per_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
simPagel_multi<-do.call(rbind,allSlopes_per_i)



pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_PAGEL_timeSlices_partSlopes_PLOTTED_4part_SCALED_w1Sim_83trees.pdf",sep=""),onefile=TRUE, width=4,height=10)

#dat<-allDat_MULTI_actual
#dat<-allDat_MULTI_scaled
#dat<-allDat_MULTI_scaled_4var
dat<-Pagel_MULTI_scaled
#Pval4<- dat$Pval5
#lowHigh95_4<-cbind.data.frame(dat$percentSamp-(1.96*dat$SE5),dat$percentSamp,dat$percentSamp+(1.96*dat$SE5),rep(sliceTimes,100))
#colnames(lowHigh95_4)<-c("low","mean","high","sliceTimes")

Pval1<- dat$Pval2
Pval2<- dat$Pval3
Pval3<- dat$Pval4
adjLab<-0.08
redAlpha<-0.3
greyAlpha<-0.1
vertLwd<-0.5
lowHigh95_1<-cbind.data.frame(dat$MRCA-(1.96*dat$SE2),dat$MRCA,dat$MRCA+(1.96*dat$SE2),rep(sliceTimes,83))
colnames(lowHigh95_1)<-c("low","mean","high","sliceTimes")
lowHigh95_2<-cbind.data.frame(dat$DR_harm-(1.96*dat$SE3),dat$DR_harm,dat$DR_harm+(1.96*dat$SE3),rep(sliceTimes,83))
colnames(lowHigh95_2)<-c("low","mean","high","sliceTimes")
lowHigh95_3<-cbind.data.frame(dat$DR_skew-(1.96*dat$SE4),dat$DR_skew,dat$DR_skew+(1.96*dat$SE4),rep(sliceTimes,83))
colnames(lowHigh95_3)<-c("low","mean","high","sliceTimes")
lwdCI<-1.5

Header<-"MULTIVARIATE - PGLS per time slice per tree"
#Sub<-"Actual data values"
#YlabMain<-"Actual unique effect on (log) clade richness"
#horizLine<-0
Sub<-"Standardized data values"
YlabMain<-"Standardized unique effect on (log) clade richness"
horizLine<-0

#SIMS<-simActual_Multis
#yLims1<-c(-0.01,0.24)
#yLims2<-c(0,35)
#yLims3<-c(0,3.5)
#yAxis<-NULL
datSIM<-simPagel_multi
yLims1<-c(-0.2,2.1)
yLims2<-c(-0.2,2.1)
yLims3<-c(-0.2,2.1)
yLims4<-c(-1,1)
yAxis<-c(NA,0,0.5,1,1.5,2,NA)
simAlpha<-0.05

j<-1

#layout(matrix(c(1:5), 5, 1, byrow = TRUE), widths=4, heights=c(4,2.5,2.5,2.5,2.5))#c(3,3,3))
layout(matrix(c(1:4), 4, 1, byrow = TRUE), widths=4, heights=c(2.5,2.5,2.5,2.5))#c(3,3,3))
par(oma = c(5,4,5,3) + 0.1, mar = rep(0.5,4) + 0.1)
#par(mfrow = c(4,1), oma = c(5,4,5,0) + 0.1, mar = c(5,1,1,1) + 0.1)
#par(op)

#part A
# PAGEL's LAMBDA
cols<-rep(grey(0.3,alpha=0.5),length(dat[,1]))
pchs<-rep(1,length(dat[,1]))

plot(Lam_mean ~ sliceTimes, data=dat,col=cols, pch=pchs, ylim=c(-8,3), ylab="",xlab="", xaxt="n")#,yaxt="n",ylab="Partial residual", xlab="Slice times (Ma)",main="")
#axis(side=2,at=c(2,1,0,-1,-2,-3,-4),labels=TRUE)
axis(side=1,at=NULL,labels=FALSE)

mtext(text="(a) Pagel's lambda", padj=2, adj=adjLab, font=2)#at=c(-72,0.21)
#mtext(text="(a) Crown age (squared)", padj=2, adj=adjLab, font=2)#at=c(-72,0.21)
#mtext(side=3,text=Sub)

for(i in 1:numSlices){
	abline(v=(-5*i),lty=2, lwd=vertLwd, col=grey(0.6, alpha=0.5))
}
	abline(h=horizLine,lty=1, lwd=1, col=grey(0.6, alpha=0.5))

mtext(side=3,text=Sub)

# part B
##
Col<-"darkgoldenrod2"
Col1<-rgb((col2rgb(Col)/255)[1],(col2rgb(Col)/255)[2],(col2rgb(Col)/255)[3],alpha=simAlpha)
Col<-"deepskyblue2"
Col2<-rgb((col2rgb(Col)/255)[1],(col2rgb(Col)/255)[2],(col2rgb(Col)/255)[3],alpha=simAlpha)
Col<-"darkorchid3"
Col3<-rgb((col2rgb(Col)/255)[1],(col2rgb(Col)/255)[2],(col2rgb(Col)/255)[3],alpha=simAlpha)
#ColorsSIM<-c(Col1,Col2,Col3)
ColorsSIM_full<-c("darkgoldenrod2","deepskyblue2","darkorchid3")
ColorsSIM<-ColorsSIM_full

#plot(MRCA ~ sliceTimes, data=dat, type="n",col=cols, pch=pchs, ylim=yLims, ylab="",xlab="", yaxt="n",xaxt="n")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
plot(MRCA ~ sliceTimes, data=dat, type="n",col=cols, pch=pchs, ylim=yLims1, ylab="",xlab="", yaxt="n",xaxt="n")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
axis(side=2,at=yAxis,labels=TRUE)
axis(side=1,at=NULL,labels=FALSE)

mtext(text="(b) Crown age", padj=2, adj=adjLab, font=2)#at=c(-72,0.21)
#mtext(text="(a) Crown age (squared)", padj=2, adj=adjLab, font=2)#at=c(-72,0.21)
#mtext(side=3,text=Sub)

cols<-rep(grey(0.3,alpha=greyAlpha),length(dat[,1]))
cols[which(Pval1>0.05)]<-rgb(1,0,0,alpha=redAlpha)
pchs<-rep(1,length(dat[,1]))
pchs[which(Pval1>0.05)]<-17
#sfracs<-rep(0,1400)
#sfracs[which(lowHigh95_1$high==max(lowHigh95_1$high))]<-3
#cols[which(lowHigh95_1$high==max(lowHigh95_1$high))]<-grey(0.3)

plotCI(x=lowHigh95_1$sliceTimes, add=TRUE, y=lowHigh95_1$mean,ui=lowHigh95_1$high,li=lowHigh95_1$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=cols,scol=grey(0.3,alpha=greyAlpha),pch=pchs,font.lab=2,cex.axis=1.1,cex.lab=1.1)

	p1<- datSIM$Pval1

	cols<-rep(ColorsSIM[j],length(datSIM[,1]))
	cols[which(p1>0.05)]<-ColorsSIM_full[j]
	pchs<-rep(1,length(datSIM[,1]))
	pchs[which(p1>0.05)]<-17

	lowHigh95_s1<-cbind.data.frame(datSIM$MRCA-(1.96*datSIM$SE2),datSIM$MRCA,datSIM$MRCA+(1.96*datSIM$SE2),rep(sliceTimes,3))
	colnames(lowHigh95_s1)<-c("low","mean","high","sliceTimes")

	plotCI(x=lowHigh95_s1$sliceTimes, add=TRUE, y=lowHigh95_s1$mean,ui=lowHigh95_s1$high,li=lowHigh95_s1$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=cols,scol=cols,pch=pchs,font.lab=2,cex.axis=1.1,cex.lab=1.1)


for(i in 1:numSlices){
	abline(v=(-5*i),lty=2, lwd=vertLwd, col=grey(0.6, alpha=0.5))
}
	abline(h=horizLine,lty=1, lwd=1, col=grey(0.6, alpha=0.5))


# part C
plot(DR_harm ~ sliceTimes, data=dat, type="n",col=cols, pch=pchs, ylim=yLims2, ylab="",xlab="", yaxt="n",xaxt="n")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
axis(side=2,at=yAxis,labels=TRUE)
axis(side=1,at=NULL,labels=FALSE)
mtext(text="(c) Tip DR mean", padj=2, adj=adjLab, font=2)#at=c(-72,0.21)

cols<-rep(grey(0.3,alpha=greyAlpha),length(dat[,1]))
cols[which(Pval2>0.05)]<-rgb(1,0,0,alpha=redAlpha)
pchs<-rep(1,length(dat[,1]))
pchs[which(Pval2>0.05)]<-17

plotCI(x=lowHigh95_2$sliceTimes, add=TRUE, y=lowHigh95_2$mean,ui=lowHigh95_2$high,li=lowHigh95_2$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=cols,scol=grey(0.3,alpha=greyAlpha),pch=pchs,font.lab=2,cex.axis=1.1,cex.lab=1.1)

	p2<- datSIM$Pval2

	cols<-rep(ColorsSIM[j],length(datSIM[,1]))
	cols[which(p2>0.05)]<-ColorsSIM_full[j]
	pchs<-rep(1,length(datSIM[,1]))
	pchs[which(p2>0.05)]<-17

	lowHigh95_s2<-cbind.data.frame(datSIM$DR_harm-(1.96*datSIM$SE3),datSIM$DR_harm,datSIM$DR_harm+(1.96*datSIM$SE3),rep(sliceTimes,3))
	colnames(lowHigh95_s2)<-c("low","mean","high","sliceTimes")

	plotCI(x=lowHigh95_s2$sliceTimes, add=TRUE, y=lowHigh95_s2$mean,ui=lowHigh95_s2$high,li=lowHigh95_s2$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=cols,scol=cols,pch=pchs,font.lab=2,cex.axis=1.1,cex.lab=1.1)


for(i in 1:numSlices){
	abline(v=(-5*i),lty=2, lwd=vertLwd, col=grey(0.6, alpha=0.5))
}
	abline(h=horizLine,lty=1, lwd=1, col=grey(0.6, alpha=0.5))


# part D
plot(DR_skew ~ sliceTimes, data=dat, type="n",col=cols, pch=pchs, ylim=yLims3, ylab="",xlab="", yaxt="n")#,xaxt="n")#,ylab="Partial residual", xlab="Slice times (Ma)",main="")
axis(side=2,at=yAxis,labels=TRUE)
mtext(text="(d) Tip DR skewness", padj=2, adj=adjLab, font=2)#at=c(-72,0.21)

cols<-rep(grey(0.3,alpha=greyAlpha),length(dat[,1]))
cols[which(Pval3>0.05)]<-rgb(1,0,0,alpha=redAlpha)
pchs<-rep(1,length(dat[,1]))
pchs[which(Pval3>0.05)]<-17

plotCI(x=lowHigh95_3$sliceTimes, add=TRUE, y=lowHigh95_3$mean,ui=lowHigh95_3$high,li=lowHigh95_3$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=cols,scol=grey(0.3,alpha=greyAlpha),pch=pchs,font.lab=2,cex.axis=1.1,cex.lab=1.1)

	p3<- datSIM$Pval3

	cols<-rep(ColorsSIM[j],length(datSIM[,1]))
	cols[which(p3>0.05)]<-ColorsSIM_full[j]
	pchs<-rep(1,length(datSIM[,1]))
	pchs[which(p3>0.05)]<-17

	lowHigh95_s3<-cbind.data.frame(datSIM$DR_skew-(1.96*datSIM$SE4),datSIM$DR_skew,datSIM$DR_skew+(1.96*datSIM$SE4),rep(sliceTimes,3))
	colnames(lowHigh95_s3)<-c("low","mean","high","sliceTimes")

	plotCI(x=lowHigh95_s3$sliceTimes, add=TRUE, y=lowHigh95_s3$mean,ui=lowHigh95_s3$high,li=lowHigh95_s3$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=cols,scol=cols,pch=pchs,font.lab=2,cex.axis=1.1,cex.lab=1.1)

for(i in 1:numSlices){
	abline(v=(-5*i),lty=2, lwd=vertLwd, col=grey(0.6, alpha=0.5))
}
	abline(h=horizLine,lty=1, lwd=1, col=grey(0.6, alpha=0.5))


title(main=Header, xlab = "Time slices before present (Ma)",
      ylab = YlabMain,
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)

dev.off()




#######################################
# ==============
#  Also get the UNIVARIATE PENDANTS going...
# ===========
#####
# Load in the UNIVARIATE standardized data
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_SIMULATIONS_on_NDexp")
library(nlme); library(MASS); library(AICcmodavg); library(ape); library(phytools)
bbone<- "NDexp" #"FBD" # 
numTrees<-100
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery
sliceTimes<-seq(-5,-70,-5)

setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut")

# Load sims
sims<-c("mamPhyE", "highE_0p8", "lowE_0p2")

SIMS_allDat_UNI_scaled<-vector("list",length=3)
for (i in 1:length(SIMS_allDat_MULTI_actual)){
	SIMS_allDat_UNI_scaled[[i]]<-read.table(paste("MamPhy_SIMS_",sims[i],"_",bbone,"_sample100_ALL100-TABLE_PGLS-UNI_BROWNIAN_timeSlices_SCALED_all.txt",sep=""))
}

# Load EMPIRICAL DATA
allDat_UNI_scaled<-read.table(file=paste(bbone,"_sample100_ALL100-TABLE_PGLS_BROWNIAN_timeSlices_wPercentSamp_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))
allDat_UNI_scaled_gls<-read.table(file=paste(bbone,"_sample100_ALL100-TABLE_PGLS_noTREE_timeSlices_wPercentSamp_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))

## TIMESLICES-- UNI -- SCALED
#allDat_UNI_i<-vector("list",length=numTrees)
#for(i in 1:numTrees){
#	tree<-rep(i,length(slice))
#	dat<-read.table(paste(bbone,"_sample100_",i,"_PGLS_BROWNIAN_timeSlices_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))
#	allDat_UNI_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
#}
#allDat_UNI_scaled<-do.call(rbind,allDat_UNI_i)

uniDat10<-allDat_UNI_scaled[which(allDat_UNI_scaled$sliceTimes==-10),]
uniDat35<-allDat_UNI_scaled[which(allDat_UNI_scaled$sliceTimes==-35),]
uniDat60<-allDat_UNI_scaled[which(allDat_UNI_scaled$sliceTimes==-60),]

# FAMS -- UNI -- SCALED
#uniPGLS_fams_SCALED<-read.table(paste(bbone,"_sample100_PGLSuniAll_FAMS_112gr2_SCALED.txt",sep=""))
#uniPGLS_fams_SCALED

# All the datas to plot. -- 6, across each of 10 vars.
uniDat10
uniDat35
uniDat60
Sim1<-SIMS_allDat_UNI_scaled[[1]]
Sim2<-SIMS_allDat_UNI_scaled[[2]]
Sim3<-SIMS_allDat_UNI_scaled[[3]]

#range<-12:21
range<-c(12:15,18,32:35,38)

Sim1_60<-Sim1[which(Sim1$sliceTimes==-60),range]
Sim2_60<-Sim2[which(Sim2$sliceTimes==-60),range]
Sim3_60<-Sim3[which(Sim3$sliceTimes==-60),range]

Sim1_35<-Sim1[which(Sim1$sliceTimes==-35),range]
Sim2_35<-Sim2[which(Sim2$sliceTimes==-35),range]
Sim3_35<-Sim3[which(Sim3$sliceTimes==-35),range]

Sim1_10<-Sim1[which(Sim1$sliceTimes==-10),range]
Sim2_10<-Sim2[which(Sim2$sliceTimes==-10),range]
Sim3_10<-Sim3[which(Sim3$sliceTimes==-10),range]

Sims_60<-rbind(Sim1_60,Sim2_60,Sim3_60)
Sims_35<-rbind(Sim1_35,Sim2_35,Sim3_35)
Sims_10<-rbind(Sim1_10,Sim2_10,Sim3_10)


##
# Now plot, incl CIs, using plotrix
library(plotrix)

# for means of ALL slices in SIMS
ALL_datasets<-list(uniPGLS_fams_SCALED[,14:23],uniDat60[,15:24], uniDat35[,15:24], uniDat10[,15:24],SIMS_allDat_UNI_scaled[[1]][,12:21],SIMS_allDat_UNI_scaled[[2]][,12:21],SIMS_allDat_UNI_scaled[[3]][,12:21])

# for 60, 35, 10 slices of SIMS -- gives 13 datasets...
# AS SEQUENTIAL of 60,60,60, etc 
ALL_datasets<-list(uniPGLS_fams_SCALED[,14:23],uniDat60[,15:24], Sim1_60, Sim2_60, Sim3_60, uniDat35[,15:24], Sim1_35, Sim2_35, Sim3_35,  uniDat10[,15:24], Sim1_10, Sim2_10, Sim3_10)

# AS SEQUENTIAL, but w MEAN of the 3 sims in RED
ALL_datasets<-list(uniPGLS_fams_SCALED[,14:23],uniDat60[,15:24], Sims_60, uniDat35[,15:24], Sims_35, uniDat10[,15:24], Sims_10)

# AS all Empirical then all SIMS, and with MamPhyE only out of the sims 
ALL_datasets<-list(uniPGLS_fams_SCALED[,14:23],uniDat60[,15:24], uniDat35[,15:24], uniDat10[,15:24], Sim1_60, Sim1_35, Sim1_10)

# AS ALTERNATING -- Empirical then SIM, and with MamPhyE only out of the sims 
ALL_datasets<-list(uniPGLS_fams_SCALED[,14:23],uniDat60[,15:24], Sim1_60, uniDat35[,15:24], Sim1_35, uniDat10[,15:24], Sim1_10)

# AS just the Crown age / Tip DR mean / Tip DR skew for each of the 60,35,10 slices-- ALL 3 sims vs empirical
ALL_datasets<-list(uniDat60[,15:17], Sim1_60, Sim2_60, Sim3_60, uniDat35[,15:17], Sim1_35, Sim2_35, Sim3_35, uniDat10[,15:17], Sim1_10,Sim2_10,Sim3_10)

# ONLY the percentSamp guys for every timeSlice...
ALL_datasets<-vector("list",length=14)
	means<-vector()
	lower95s<-vector()
	upper95s<-vector()
for (i in 1:14){
	ALL_datasets[[i]]<-allDat_UNI_scaled[which(allDat_UNI_scaled$slice==i),"percentSamp"]
		means[i]<-mean(ALL_datasets[[i]])
		dd<-quantile(ALL_datasets[[i]], c(0.025,0.975))
		lower95s[i]<-dd[[1]]
		upper95s[i]<-dd[[2]]
}
pdf(file=paste("cladeLevel",bbone,"_PGLS_UNI_percentSamp_14slices_95ci.pdf",sep=""),onefile=TRUE, width=10,height=5)
plotCI(x=seq(-5,-70,-5), y=means,ui=upper95s,li=lower95s, cex=1.3,sfrac=0, xlab="Time before present (Ma)", ylab="Standardized effect on (log) clade richness", err="y", lwd=Lwd,col="grey",scol="grey",pch=16,font.lab=2,cex.axis=1.1,cex.lab=1.1)
abline(h=0,lty=1, lwd=2, col=grey(0.6, alpha=0.5)) #rgb(1,0,0,alpha=0.5))#
text(x= -65, y = 1.4, cex=1.1,font=2,srt = 0, pos=4,labels = "Percent sampled tips, per clade per slice", xpd = TRUE, adj=c(0.95,0.05))
dev.off()
#####
# The 95% CI around all mean estimates per var (OLD)
ALL_means<-vector("list",length(ALL_datasets))
ALL_lower95s<-vector("list",length(ALL_datasets))
ALL_upper95s<-vector("list",length(ALL_datasets))
for (i in 1:length(ALL_datasets)){
	means<-vector()
	lower95s<-vector()
	upper95s<-vector()
	dat<-ALL_datasets[[i]]
	for (j in 1:length(dat)){
		means[j]<-mean(dat[,j])
		dd<-quantile(dat[,j], c(0.025,0.975))
		lower95s[j]<-dd[[1]]
		upper95s[j]<-dd[[2]]
	}
	ALL_means[[i]]<-means
	ALL_lower95s[[i]]<-lower95s
	ALL_upper95s[[i]]<-upper95s
}

# AS Crown age / Tip DR mean / Tip DR skew / PB Div / BD Div for each of the 60,35,10 slices-- ALL 3 sims vs empirical
# means and SEs to calc 95 CIs...
vars<-c(13:15,17,20,35:37,39,42)
vars_sim<-c(1:10)
ALL_datasets<-list(uniDat60[,vars], Sim1_60[,vars_sim], Sim2_60[,vars_sim], Sim3_60[,vars_sim], uniDat35[,vars], Sim1_35[,vars_sim], Sim2_35[,vars_sim], Sim3_35[,vars_sim], uniDat10[,vars], Sim1_10[,vars_sim],Sim2_10[,vars_sim],Sim3_10[,vars_sim])

# SEs to calc 95 CIs...
#vars<-c(35:37,39,42)
#vars_sim<-c(6:10)
#ALL_SEs<-list(uniDat60[,vars], Sim1_60[,vars_sim], Sim2_60[,vars_sim], Sim3_60[,vars_sim], uniDat35[,vars], Sim1_35[,vars_sim], Sim2_35[,vars_sim], Sim3_35[,vars_sim], uniDat10[,vars], Sim1_10[,vars_sim],Sim2_10[,vars_sim],Sim3_10[,vars_sim])

# The 95% CI for each estimate of each var, plotted together
ALLdata_fiveVarsToPlot<-vector("list",length(ALL_datasets))

for (i in 1:length(ALL_datasets)){
	dat<-ALL_datasets[[i]]
	fiveVarsToPlot_i<-vector("list",length=5)
	for (j in 1:5){
		lowHigh95_j<-cbind.data.frame(dat[,j]-(1.96*dat[,j+5]),dat[,j],dat[,j]+(1.96*dat[,j+5]),rep(j,100))
		colnames(lowHigh95_j)<-c("low","mean","high","xVal")
		fiveVarsToPlot_i[[j]]<-lowHigh95_j
	}
	ALLdata_fiveVarsToPlot[[i]]<-fiveVarsToPlot_i
}


#varNames<-paste(c("Crown age","Tip DR mean", "Tip DR skewness", "PB net div", "BD lambda", "BD mu", "BD net div", "MS net div e=0","MS net div e=0.5", "MS net div e=0.9"))
varNames<-paste(c("Crown age","Tip DR mean", "Tip DR skew.","Clade DR (PB)","Clade DR (BD)"))
#varNames2<-paste(c("(a)","(b)", "(c)",))

####
# PLOT with 5 variables and ERROR BARS-- each comparing Empirical to ALL 3 SIMS -- across 3 time slices...
###


#layout(matrix(c(1:3), 1, 3, byrow = TRUE), widths=c(2,2,2), heights=2)
pdf(file=paste("cladeLevel",bbone,"_PGLS_UNI_wALL3SIMS-vsEmpirical_5VARS_bySlice_full95CI-fromSE.pdf",sep=""),onefile=TRUE, width=10,height=5)
#par(oma = rep(5,4) + 0.1, mar = rep(0.5,4) + 0.1)

LwdCI<-1.5
greyAlpha<-0.2
empirPointCol<-grey(0.3,alpha=0.3)
pchs<-c(rep(1,4),rep(2,4),rep(0,4)) # circle, triangle, square

simAlpha<-0.2
Col<-"darkgoldenrod2"
Col1<-rgb((col2rgb(Col)/255)[1],(col2rgb(Col)/255)[2],(col2rgb(Col)/255)[3],alpha=simAlpha)
Col<-"deepskyblue2"
Col2<-rgb((col2rgb(Col)/255)[1],(col2rgb(Col)/255)[2],(col2rgb(Col)/255)[3],alpha=simAlpha)
Col<-"darkorchid3"
Col3<-rgb((col2rgb(Col)/255)[1],(col2rgb(Col)/255)[2],(col2rgb(Col)/255)[3],alpha=simAlpha)
ColorsSIM<-c(Col1,Col2,Col3)
ColorsSIM_full<-c("darkgoldenrod2","deepskyblue2","darkorchid3")


#quartz(width=10,height=5)
dat<-ALLdata_fiveVarsToPlot[[1]][[1]]
plotCI(x=dat$xVal, y=dat$mean,ui=dat$high,li=dat$low, cex=1.1,xlim=c(1,15.4),ylim=c(-0.5,2.2), sfrac=0,xaxt="n", xlab="", ylab="Standardized effect on (log) clade richness", err="y", lwd=LwdCI,col=empirPointCol,scol=grey(0.3,alpha=greyAlpha),pch=pchs[1],font.lab=2,cex.axis=1.1,cex.lab=1.1)
abline(h=0,lty=1, lwd=2, col=grey(0.6, alpha=0.5)) #rgb(1,0,0,alpha=0.5))#
for(j in 2:5){
	dat<-ALLdata_fiveVarsToPlot[[1]][[j]]
	plotCI(x=dat$xVal, add=TRUE, y=dat$mean,ui=dat$high,li=dat$low, cex=1.1,xlim=c(1,15.4),ylim=c(-0.5,2), sfrac=0,xaxt="n", xlab="", ylab="", err="y", lwd=LwdCI,col=empirPointCol,scol=grey(0.3,alpha=greyAlpha),pch=pchs[1],font.lab=2,cex.axis=1.1,cex.lab=1.1)
}

for(i in 2:4){
simDat<-ALLdata_fiveVarsToPlot[[i]]
for(j in 1:5){
	dat<-simDat[[j]]
	plotCI(x=dat$xVal+0.2*(i-1), add=TRUE, y=dat$mean,ui=dat$high,li=dat$low, cex=1.1,xlim=c(1,15.4),ylim=c(-0.5,2), sfrac=0,xaxt="n", xlab="", ylab="", err="y", lwd=LwdCI,col=ColorsSIM[i-1],scol=ColorsSIM[i-1],pch=pchs[i],font.lab=2,cex.axis=1.1,cex.lab=1.1)
	}
}

i<-5
empirDat<-ALLdata_fiveVarsToPlot[[i]]
for(j in 1:5){
	dat<-empirDat[[j]]
	plotCI(x=(dat$xVal+5), add=TRUE, y=dat$mean,ui=dat$high,li=dat$low, cex=1.1,xlim=c(1,15.4),ylim=c(-0.5,2), sfrac=0,xaxt="n", xlab="", ylab="", err="y", lwd=LwdCI,col=empirPointCol,scol=grey(0.3,alpha=greyAlpha),pch=pchs[i],font.lab=2,cex.axis=1.1,cex.lab=1.1)
}

for(i in 6:8){
simDat<-ALLdata_fiveVarsToPlot[[i]]
for(j in 1:5){
	dat<-simDat[[j]]
	plotCI(x=(dat$xVal+5)+0.2*(i-5), add=TRUE, y=dat$mean,ui=dat$high,li=dat$low, cex=1.1,xlim=c(1,15.4),ylim=c(-0.5,2), sfrac=0,xaxt="n", xlab="", ylab="", err="y", lwd=LwdCI,col=ColorsSIM[i-5],scol=ColorsSIM[i-5],pch=pchs[i],font.lab=2,cex.axis=1.1,cex.lab=1.1)
	}
}

i<-9
empirDat<-ALLdata_fiveVarsToPlot[[i]]
for(j in 1:5){
	dat<-empirDat[[j]]
	plotCI(x=(dat$xVal+10), add=TRUE, y=dat$mean,ui=dat$high,li=dat$low, cex=1.1,xlim=c(1,15.4),ylim=c(-0.5,2), sfrac=0,xaxt="n", xlab="", ylab="", err="y", lwd=LwdCI,col=empirPointCol,scol=grey(0.3,alpha=greyAlpha),pch=pchs[i],font.lab=2,cex.axis=1.1,cex.lab=1.1)
}

for(i in 10:12){
simDat<-ALLdata_fiveVarsToPlot[[i]]
for(j in 1:5){
	dat<-simDat[[j]]
	plotCI(x=(dat$xVal+10)+0.2*(i-9), add=TRUE, y=dat$mean,ui=dat$high,li=dat$low, cex=1.1,xlim=c(1,15.4),ylim=c(-0.5,2), sfrac=0,xaxt="n", xlab="", ylab="", err="y", lwd=LwdCI,col=ColorsSIM[i-9],scol=ColorsSIM[i-9],pch=pchs[i],font.lab=2,cex.axis=1.1,cex.lab=1.1)
	}
}

for(i in c(5,10)){
	abline(v=(i+0.75),lty=1, lwd=2, col="black")
}

axis(1, at=c(1:15)+0.25, labels = FALSE)
text(x= c(1:15)+0.25, y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), cex=1.1,font=2,srt = 45, labels = varNames, xpd = TRUE, adj=c(0.95,0.05))

#text(x= c(1:9)+0.35, y = -0.5, cex=1.1,font=2,srt = 0, labels = varNames2, xpd = TRUE, adj=c(0.95,0.05))

dev.off()




####
# PLOT with just 5 variables, each comparing Empirical to ALL 3 SIMS -- across 3 time slices...
###

cols2<-rep(c("black","darkgoldenrod2","deepskyblue2","darkorchid3"),3)
colPerDataset_bars<-vector("list",length(cols2))
for (i in 1:length(cols2)){
	colPerDataset_bars[[i]]<-rep(cols2[i],5)
}
colPerDataset_bars[[1]]<-c("grey","black","black","grey","grey")
colPerDataset_bars[[9]]<-c("black","black","black","grey","black")


Lwd<-3
pchs<-c(rep(16,4),rep(17,4),rep(15,4)) # circle, triangle, square

#layout(matrix(c(1:3), 1, 3, byrow = TRUE), widths=c(2,2,2), heights=2)
pdf(file=paste("cladeLevel",bbone,"_PGLS_UNI_wALL3SIMS-vsEmpirical_5VARS_bySlice_95ci.pdf",sep=""),onefile=TRUE, width=10,height=5)
#par(oma = rep(5,4) + 0.1, mar = rep(0.5,4) + 0.1)

#quartz(width=10,height=5)
plotCI(x=1:5, y=ALL_means[[1]],ui=ALL_upper95s[[1]],li=ALL_lower95s[[1]], cex=1.3,xlim=c(1,15.4),ylim=c(-0.5,2), sfrac=0,xaxt="n", xlab="", ylab="Standardized effect on (log) clade richness", err="y", lwd=Lwd,col=colPerDataset_bars[[1]],scol=colPerDataset_bars[[1]],pch=pchs[1],font.lab=2,cex.axis=1.1,cex.lab=1.1)
abline(h=0,lty=1, lwd=2, col=grey(0.6, alpha=0.5)) #rgb(1,0,0,alpha=0.5))#
for (i in 1:3){
	plotCI(x=(1:5)+0.14*i, add=TRUE,  y=ALL_means[[i+1]],ui=ALL_upper95s[[i+1]],li=ALL_lower95s[[i+1]], cex=1.3, sfrac=0, err="y", lwd=Lwd,col=colPerDataset_bars[[i+1]],scol=colPerDataset_bars[[i+1]],pch=pchs[i+1])
}

j<-5
plotCI(x=6:10, y=ALL_means[[j]],add=TRUE, ui=ALL_upper95s[[j]],li=ALL_lower95s[[j]], cex=1.3,xlim=c(0.8,12.65),ylim=c(-0.5,2), sfrac=0,xaxt="n", xlab="", ylab="Standardized effect on (log) clade richness", err="y", lwd=Lwd,col=colPerDataset_bars[[j]],scol=colPerDataset_bars[[j]],pch=pchs[j],font.lab=2,cex.axis=1.1,cex.lab=1.1)
for (i in 5:7){
	plotCI(x=(6:10)+0.14*(i-4), add=TRUE,  y=ALL_means[[i+1]],ui=ALL_upper95s[[i+1]],li=ALL_lower95s[[i+1]], cex=1.3, sfrac=0, err="y", lwd=Lwd,col=colPerDataset_bars[[i+1]],scol=colPerDataset_bars[[i+1]],pch=pchs[i+1])
}

j<-9
plotCI(x=11:15, y=ALL_means[[j]],add=TRUE, ui=ALL_upper95s[[j]],li=ALL_lower95s[[j]], cex=1.3,xlim=c(0.8,12.65),ylim=c(-0.5,2), sfrac=0,xaxt="n", xlab="", ylab="Standardized effect on (log) clade richness", err="y", lwd=Lwd,col=colPerDataset_bars[[j]],scol=colPerDataset_bars[[j]],pch=pchs[j],font.lab=2,cex.axis=1.1,cex.lab=1.1)
for (i in 9:11){
	plotCI(x=(11:15)+0.14*(i-8), add=TRUE,  y=ALL_means[[i+1]],ui=ALL_upper95s[[i+1]],li=ALL_lower95s[[i+1]], cex=1.3, sfrac=0, err="y", lwd=Lwd,col=colPerDataset_bars[[i+1]],scol=colPerDataset_bars[[i+1]],pch=pchs[i+1])
}

for(i in c(5,10)){
	abline(v=(i+0.7),lty=1, lwd=2, col="black")
}

axis(1, at=c(1:15)+0.25, labels = FALSE)
text(x= c(1:15)+0.25, y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), cex=1.1,font=2,srt = 45, labels = varNames, xpd = TRUE, adj=c(0.95,0.05))

#text(x= c(1:9)+0.35, y = -0.5, cex=1.1,font=2,srt = 0, labels = varNames2, xpd = TRUE, adj=c(0.95,0.05))

dev.off()

####
# PLOT with just the first 3 variables compared by Empirical - to ALL 3 SIMS for each TIME-SLICE, showing in its native color...
###

cols2<-rep(c("black","darkgoldenrod2","deepskyblue2","darkorchid3"),3)
colPerDataset_bars<-vector("list",length(cols2))
for (i in 1:length(cols2)){
	colPerDataset_bars[[i]]<-rep(cols2[i],3)
}

Lwd<-3
pchs<-c(rep(16,12),rep(17,12),rep(15,12)) # circle, triangle, square

#layout(matrix(c(1:3), 1, 3, byrow = TRUE), widths=c(2,2,2), heights=2)
#par(oma = c(5,4,5,3) + 0.1, mar = rep(0.5,4) + 0.1)
#quartz(width=10,height=6)
pdf(file=paste("cladeLevel",bbone,"_PGLS_UNI_wALL3SIMS-vsEmpirical_VARSbySlice_95ci.pdf",sep=""),onefile=TRUE, width=10,height=5)

plotCI(x=1:3, y=ALL_means[[1]],ui=ALL_upper95s[[1]],li=ALL_lower95s[[1]], cex=1.3,xlim=c(1,9.4),ylim=c(-0.5,2), sfrac=0,xaxt="n", xlab="", ylab="Standardized effect on (log) clade richness", err="y", lwd=Lwd,col=colPerDataset_bars[[1]],scol=colPerDataset_bars[[1]],pch=pchs[1],font.lab=2,cex.axis=1.1,cex.lab=1.1)
abline(h=0,lty=1, lwd=2, col=grey(0.6, alpha=0.5)) #rgb(1,0,0,alpha=0.5))#
for (i in 1:3){
	plotCI(x=(1:3)+0.14*i, add=TRUE,  y=ALL_means[[i+1]],ui=ALL_upper95s[[i+1]],li=ALL_lower95s[[i+1]], cex=1.3, sfrac=0, err="y", lwd=Lwd,col=colPerDataset_bars[[i+1]],scol=colPerDataset_bars[[i+1]],pch=pchs[i+1])
}
j<-5
plotCI(x=4:6, y=ALL_means[[j]],add=TRUE, ui=ALL_upper95s[[j]],li=ALL_lower95s[[j]], cex=1.3,xlim=c(0.8,12.65),ylim=c(-0.5,2), sfrac=0,xaxt="n", xlab="", ylab="Standardized effect on (log) clade richness", err="y", lwd=Lwd,col=colPerDataset_bars[[j]],scol=colPerDataset_bars[[j]],pch=pchs[j],font.lab=2,cex.axis=1.1,cex.lab=1.1)
for (i in 5:7){
	plotCI(x=(4:6)+0.14*(i-4), add=TRUE,  y=ALL_means[[i+1]],ui=ALL_upper95s[[i+1]],li=ALL_lower95s[[i+1]], cex=1.3, sfrac=0, err="y", lwd=Lwd,col=colPerDataset_bars[[i+1]],scol=colPerDataset_bars[[i+1]],pch=pchs[i+1])
}
j<-9
plotCI(x=7:9, y=ALL_means[[j]],add=TRUE, ui=ALL_upper95s[[j]],li=ALL_lower95s[[j]], cex=1.3,xlim=c(0.8,12.65),ylim=c(-0.5,2), sfrac=0,xaxt="n", xlab="", ylab="Standardized effect on (log) clade richness", err="y", lwd=Lwd,col=colPerDataset_bars[[j]],scol=colPerDataset_bars[[j]],pch=pchs[j],font.lab=2,cex.axis=1.1,cex.lab=1.1)
for (i in 9:11){
	plotCI(x=(7:9)+0.14*(i-8), add=TRUE,  y=ALL_means[[i+1]],ui=ALL_upper95s[[i+1]],li=ALL_lower95s[[i+1]], cex=1.3, sfrac=0, err="y", lwd=Lwd,col=colPerDataset_bars[[i+1]],scol=colPerDataset_bars[[i+1]],pch=pchs[i+1])
}

for(i in c(3,6)){
	abline(v=(i+0.7),lty=1, lwd=2, col="black")
}

axis(1, at=c(1:9)+0.25, labels = FALSE)
text(x= c(1:3)+0.25, y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), cex=1.1,font=2,srt = 45, labels = varNames, xpd = TRUE, adj=c(0.95,0.05))
text(x= c(1:9)+0.35, y = -0.5, cex=1.1,font=2,srt = 0, labels = varNames2, xpd = TRUE, adj=c(0.95,0.05))

dev.off()

####
# PLOT with the 0.65 SIMS (MamPhyE) for each TIME-SLICE, showing in its native color...
# but with ALTERNATING mamPhy-to-SIM-to-mamPhy...
###

cols2<-c("black",rep(c("black","darkgoldenrod2"),3))
colPerDataset_bars<-vector("list",length(cols2))
for (i in 1:length(cols2)){
	colPerDataset_bars[[i]]<-rep(cols2[i],10)
}

Lwd<-3
pchs<-c(23,1,1,2,2,0,0)

pdf(file=paste("cladeLevel",bbone,"_PGLS_UNI_all7datasets_wSIMS-1mamPhyE_alternate_VARSbyCat_95ci.pdf",sep=""),onefile=TRUE, width=10,height=5)
#quartz(width=10,height=6)
plotCI(x=1:10, y=ALL_means[[1]],ui=ALL_upper95s[[1]],li=ALL_lower95s[[1]], cex=1.3,xlim=c(1.2,10.65),ylim=c(-0.5,2), sfrac=0,xaxt="n", xlab="", ylab="Standardized effect on (log) clade richness", err="y", lwd=Lwd,col=colPerDataset_bars[[1]],scol=colPerDataset_bars[[1]],pch=pchs[1],font.lab=2,cex.axis=1.1,cex.lab=1.1)
abline(h=0,lty=1, lwd=2, col=grey(0.6, alpha=0.5)) #rgb(1,0,0,alpha=0.5))#
for (i in 1:6){
	plotCI(x=(1:10)+0.14*i, add=TRUE,  y=ALL_means[[i+1]],ui=ALL_upper95s[[i+1]],li=ALL_lower95s[[i+1]], cex=1.3, sfrac=0, err="y", lwd=Lwd,col=colPerDataset_bars[[i+1]],scol=colPerDataset_bars[[i+1]],pch=pchs[i+1])

}
for(i in 1:10){
	abline(v=(i+0.92),lty=1, lwd=2, col="black")
	abline(v=(i+0.06),lty=2, lwd=2, col=grey(0.6, alpha=0.5))
}

axis(1, at=c(1:10)+0.4, labels = FALSE)
text(x= c(1:10)+0.4, y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), cex=1.1,font=2,srt = 45, labels = varNames, xpd = TRUE, adj=c(0.95,0.05))

dev.off()


####
# PLOT with the 0.65 SIMS (MamPhyE) for each TIME-SLICE, showing in its native color...
###

cols2<-c(grey(0.6),rep("black",3), rep("darkgoldenrod2",3))
colPerDataset_bars<-vector("list",length(cols2))
for (i in 1:length(cols2)){
	colPerDataset_bars[[i]]<-rep(cols2[i],10)
}

Lwd<-3
pchs<-c(16,rep(c(16,15,17),2))

pdf(file=paste("cladeLevel",bbone,"_PGLS_UNI_all7datasets_wSIMS-1mamPhyE_VARSbyCat_95ci.pdf",sep=""),onefile=TRUE, width=10,height=5)
#quartz(width=10,height=6)
plotCI(x=1:10, y=ALL_means[[1]],ui=ALL_upper95s[[1]],li=ALL_lower95s[[1]], cex=1.4,xlim=c(1.2,10.65),ylim=c(-0.5,2), sfrac=0,xaxt="n", xlab="", ylab="Standardized effect on (log) clade richness", err="y", lwd=Lwd,col=colPerDataset_bars[[1]],scol=colPerDataset_bars[[1]],pch=pchs[1],font.lab=2,cex.axis=1.3,cex.lab=1.3)
abline(h=0,lty=1, lwd=2, col=grey(0.6, alpha=0.5)) #rgb(1,0,0,alpha=0.5))#
for (i in 1:6){
	plotCI(x=(1:10)+0.14*i, add=TRUE,  y=ALL_means[[i+1]],ui=ALL_upper95s[[i+1]],li=ALL_lower95s[[i+1]], cex=1.4, sfrac=0, err="y", lwd=Lwd,col=colPerDataset_bars[[i+1]],scol=colPerDataset_bars[[i+1]],pch=pchs[i+1])

}
for(i in 1:10){
	abline(v=(i+0.92),lty=1, lwd=2, col="black")
	abline(v=(i+0.06),lty=2, lwd=2, col=grey(0.6, alpha=0.5))
}

axis(1, at=c(1:10)+0.4, labels = FALSE)
text(x= c(1:10)+0.4, y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), cex=1.3,font=2,srt = 45, labels = varNames, xpd = TRUE, adj=c(0.95,0.05))

dev.off()

####
# PLOT with MEAN of the SIMS in one-- for each TIME-SLICE, show the mean of 3 sims as RED.
###
# make the MEAN of the 3 sims...

cols2<-c("black",rep(c("black",grey(0.7)),3))
colPerDataset_bars<-vector("list",length(cols2))
for (i in 1:length(cols2)){
	colPerDataset_bars[[i]]<-rep(cols2[i],10)
}

Lwd<-3
pchs<-c(16,rep(16,2),rep(15,2),rep(17,2))

pdf(file=paste("cladeLevel",bbone,"_PGLS_UNI_all7datasets_wSIMS-1meanEa_VARSbyCat_95ci.pdf",sep=""),onefile=TRUE, width=10,height=6)
#quartz(width=10,height=6)
plotCI(x=1:10, y=ALL_means[[1]],ui=ALL_upper95s[[1]],li=ALL_lower95s[[1]], cex=1.2,xlim=c(1.2,10.65),ylim=c(-0.5,2), sfrac=0,xaxt="n", xlab="", ylab="Standardized effect on (log) clade richness", err="y", lwd=Lwd,col=colPerDataset_bars[[1]],scol=colPerDataset_bars[[1]],pch=pchs[1],font.lab=2,cex.axis=1.5,cex.lab=1.5)
abline(h=0,lty=2, lwd=2, col=rgb(1,0,0,alpha=0.5))#grey(0.6, alpha=0.5))
for (i in 1:6){
	plotCI(x=(1:10)+0.14*i, add=TRUE,  y=ALL_means[[i+1]],ui=ALL_upper95s[[i+1]],li=ALL_lower95s[[i+1]], cex=1.2, sfrac=0, err="y", lwd=Lwd,col=colPerDataset_bars[[i+1]],scol=colPerDataset_bars[[i+1]],pch=pchs[i+1])

}
for(i in 1:9){
	abline(v=(i+0.92),lty=1, lwd=2, col="black")
	abline(v=(i+0.05),lty=2, lwd=1, col=grey(0.6, alpha=0.5))
	abline(v=(i+0.31),lty=2, lwd=1, col=grey(0.6, alpha=0.5))
	abline(v=(i+0.6),lty=2, lwd=1, col=grey(0.6, alpha=0.5))
}

axis(1, at=c(1:10)+0.4, labels = FALSE)
text(x= c(1:10)+0.4, y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), cex=1.5,font=2,srt = 45, labels = varNames, xpd = TRUE, adj=c(0.95,0.05))

dev.off()




###
# PLOT with 1 avg mean per SIM

cols1<-c("white",grey(0),grey(0.6),grey(0.85),"darkgoldenrod2","deepskyblue2","darkorchid3")
colPerDataset_dots<-vector("list",length(cols1))
for (i in 1:length(cols1)){
	colPerDataset_dots[[i]]<-rep(cols1[i],10)
}

cols2<-c("black","black","black","black","darkgoldenrod2","deepskyblue2","darkorchid3")
colPerDataset_bars<-vector("list",length(cols2))
for (i in 1:length(cols2)){
	colPerDataset_bars[[i]]<-rep(cols2[i],10)
}
pchs<-rep(20,7)

pdf(file=paste("cladeLevel",bbone,"_PGLS_UNI_all7datasets_wSIMS-1400mean_VARSbyCat_95ci.pdf",sep=""),onefile=TRUE, width=10,height=6)
#quartz(width=10,height=6)
plotCI(x=1:10, y=ALL_means[[1]],ui=ALL_upper95s[[1]],li=ALL_lower95s[[1]], cex=1.8,xlim=c(1.2,10.65),ylim=c(-0.5,2), sfrac=0,xaxt="n", xlab="", ylab="Standardized effect on (log) clade richness", err="y", lwd=2,col=colPerDataset_dots[[1]],scol=colPerDataset_bars[[1]],pch=pchs[1],font.lab=2,cex.axis=0.95,cex.lab=1)
abline(h=0,lty=2, lwd=1, col=grey(0.6, alpha=0.5))
for (i in 1:6){
	plotCI(x=(1:10)+0.14*i, add=TRUE,  pt.bg=par("bg"), y=ALL_means[[i+1]],ui=ALL_upper95s[[i+1]],cex=1.8,li=ALL_lower95s[[i+1]], sfrac=0, err="y", lwd=2,col=colPerDataset_dots[[i+1]],scol=colPerDataset_bars[[i+1]],pch=pchs[i+1])

}
for(i in 1:9){
	abline(v=(i+0.92),lty=1, lwd=1, col="black")
	#abline(v=(i+0.5),lty=2, lwd=1, col=grey(0.6, alpha=0.5))
}
axis(1, at=c(1:10)+0.4, labels = FALSE)
text(x= c(1:10)+0.4, y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), cex=0.9,font=2,srt = 45, labels = varNames, xpd = TRUE, adj=c(0.95,0.05))

dev.off()


### 
# PLOT with 3 slices per SIM
#cols1<-c(grey(0),grey(0.6),grey(0.85),"white",rep("darkgoldenrod2",3),rep("deepskyblue2",3),rep("darkorchid3",3))
#cols1<-c("white",grey(0),grey(0.6),grey(0.85),rep(c("darkgoldenrod2","deepskyblue2","darkorchid3"),3))
cols1<-c("white",rep(grey(0),4),rep(grey(0.6),4),rep(grey(0.85),4))
colPerDataset_dots<-vector("list",length(cols1))
for (i in 1:length(cols1)){
	colPerDataset_dots[[i]]<-rep(cols1[i],10)
}

cols2<-c("black",rep(c("black","darkgoldenrod2","deepskyblue2","darkorchid3"),3))
colPerDataset_bars<-vector("list",length(cols2))
for (i in 1:length(cols2)){
	colPerDataset_bars[[i]]<-rep(cols2[i],10)
}

Lwd<-3
pchs<-c(16,rep(16,4),rep(15,4),rep(17,4))

pdf(file=paste("cladeLevel",bbone,"_PGLS_UNI_all7datasets_wSIMS-3slicesEa_VARSbyCat_95ci.pdf",sep=""),onefile=TRUE, width=20,height=6)
#quartz(width=10,height=6)
plotCI(x=1:10, y=ALL_means[[1]],ui=ALL_upper95s[[1]],li=ALL_lower95s[[1]], cex=1.2,xlim=c(1.2,10.65),ylim=c(-0.5,2), sfrac=0,xaxt="n", xlab="", ylab="Standardized effect on (log) clade richness", err="y", lwd=Lwd,col=colPerDataset_bars[[1]],scol=colPerDataset_bars[[1]],pch=pchs[1],font.lab=2,cex.axis=1.5,cex.lab=1.5)
abline(h=0,lty=2, lwd=2, col=rgb(1,0,0,alpha=0.5))#grey(0.6, alpha=0.5))
for (i in 1:12){
	plotCI(x=(1:10)+0.07*i, add=TRUE,  y=ALL_means[[i+1]],ui=ALL_upper95s[[i+1]],li=ALL_lower95s[[i+1]], cex=1.2, sfrac=0, err="y", lwd=Lwd,col=colPerDataset_bars[[i+1]],scol=colPerDataset_bars[[i+1]],pch=pchs[i+1])

}
for(i in 1:9){
	abline(v=(i+0.92),lty=1, lwd=2, col="black")
	abline(v=(i+0.05),lty=2, lwd=1, col=grey(0.6, alpha=0.5))
	abline(v=(i+0.31),lty=2, lwd=1, col=grey(0.6, alpha=0.5))
	abline(v=(i+0.6),lty=2, lwd=1, col=grey(0.6, alpha=0.5))
}

axis(1, at=c(1:10)+0.4, labels = FALSE)
text(x= c(1:10)+0.4, y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), cex=1.5,font=2,srt = 45, labels = varNames, xpd = TRUE, adj=c(0.95,0.05))

dev.off()




##############
# Make the SCATTER plots...
####
# Then ALSO plot each separately...
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/Diversification_analyses-Condamine/")
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_SIMULATIONS_on_NDexp")
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_SIMULATIONS_on_NDexp")
bbone<- "NDexp" #"FBD" # 
library(nlme); library(MASS); library(AICcmodavg); library(ape); library(phytools)

numTrees<-100
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery

allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

# load the empirical RAW scatter data...
allRes<-vector("list",length=numTrees)
for(i in 1:numTrees){
results<-vector("list",length(allCladeSetNames))
for (j in 1:length(allCladeSetNames)){
	results[[j]]<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""))
	rownames(results[[j]])<-paste(i,"_",j,"_",c(1:length(results[[j]][,1])),sep="")
}
allRes[[i]]<-do.call(rbind,results) # each is 19 vars * 3000 entries about
#write.table(allRes,file=paste(bbone,"_sample100_",i,"_slice5Ma-to-70Ma_cladeSTATS.txt",sep=""))
}
allRes_100trees<-do.call(rbind,allRes) # 304274 obs. of  19 variables:

save(allRes_100trees, file=paste(bbone,"_sample100__ALL100trees_slice5Ma-to-70Ma_cladeSTATS.Rdata",sep=""))


# load the SIMULATED RAW scatter data...
sims<-c("mamPhyE", "highE_0p8", "lowE_0p2")
q<-1

allRes_SIM1<-vector("list",length=numTrees)
for(i in 1:numTrees){
results<-vector("list",length(allCladeSetNames))
for (j in 1:length(allCladeSetNames)){
	results[[j]]<-read.table(paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""))
	rownames(results[[j]])<-paste(i,"_",j,"_",c(1:length(results[[j]][,1])),sep="")
}
allRes_SIM1[[i]]<-do.call(rbind,results) # each is 19 vars * 3000 entries about
#write.table(allRes,file=paste(bbone,"_sample100_",i,"_slice5Ma-to-70Ma_cladeSTATS.txt",sep=""))
}
allRes_SIM1_100trees<-do.call(rbind,allRes_SIM1) # 836514 obs. of  13 variables... cut some vars... yep. All good, depends on # clades in slices.

setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")

write.table(allRes_SIM1_100trees,file=paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100__ALL100trees_slice5Ma-to-70Ma_cladeSTATS_table.txt",sep=""))
save(allRes_SIM1_100trees, file=paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100__ALL100trees_slice5Ma-to-70Ma_cladeSTATS.Rdata",sep=""))

	# load back in...
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/")
load(file=paste(bbone,"_sample100__ALL100trees_slice5Ma-to-70Ma_cladeSTATS.Rdata",sep=""))
#load(file=paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100__ALL100trees_slice5Ma-to-70Ma_cladeSTATS.Rdata",sep=""))

# subset RAW SCATTER data
dat60<-allRes_100trees[which(allRes_100trees$slice==60),] # 2888 obs. of  19 variables:
dat35<-allRes_100trees[which(allRes_100trees$slice==35),] # 8960 obs. of  19 variables
dat10<-allRes_100trees[which(allRes_100trees$slice==10),] # 65598 obs. of  19 variables
# get a sample from the 10 Ma data...
#library(dplyr)
#dat10<-sample_n(dat10all,10000) # 8960 obs. of  19 variables: awesome.
dats<-list(dat60,dat35,dat10)

# get RAW scatter data from tree 1 ONLY
dat60_t1<-allRes_100trees[which(allRes_100trees$slice==60 & allRes_100trees$tree==1),] # 2888 obs. of  19 variables:
dat35_t1<-allRes_100trees[which(allRes_100trees$slice==35 & allRes_100trees$tree==1),] # 8960 obs. of  19 variables
dat10_t1<-allRes_100trees[which(allRes_100trees$slice==10 & allRes_100trees$tree==1),] # 65598 obs. of  19 variables
# get a sample from the 10 Ma data...
#library(dplyr)
#dat10<-sample_n(dat10all,10000) # 8960 obs. of  19 variables: awesome.
dats_t1<-list(dat60_t1,dat35_t1,dat10_t1)



# subset SIM raw scatter data...
SIMdat60<-allRes_SIM1_100trees[which(allRes_SIM1_100trees$slice==60),] # 16503 obs. of  13 variables:
SIMdat35<-allRes_SIM1_100trees[which(allRes_SIM1_100trees$slice==35),] # 48752 obs. of  13 variables:
SIMdat10<-allRes_SIM1_100trees[which(allRes_SIM1_100trees$slice==10),] # 139349 obs. of  13 variables:
SIMdats<-list(SIMdat60,SIMdat35,SIMdat10)

# load in SLOPES univariate actual...

# UNI -- ACTUAL
allDat_UNI_i<-vector("list",length=numTrees)
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	ints<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniINTS.txt",sep=""))
	slopes<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniSLOPES.txt",sep=""))
	Ps<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlice_5Ma-to-70Ma_uniPs.txt",sep=""))
	dat<-cbind(ints,slopes,Ps)
	allDat_UNI_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
allDat_UNI_actual<-do.call(rbind,allDat_UNI_i)

# you need all vals, incl intercepts and Ps
slopes60<-allDat_UNI_actual[which(allDat_UNI_actual$sliceTimes==-60),] #15:24
slopes35<-allDat_UNI_actual[which(allDat_UNI_actual$sliceTimes==-35),]
slopes10<-allDat_UNI_actual[which(allDat_UNI_actual$sliceTimes==-10),]
slopes<-list(slopes60,slopes35,slopes10)

# you need all vals, incl intercepts and Ps-- for tree 1 ONLY
slopes60_t1<-allDat_UNI_actual[which(allDat_UNI_actual$sliceTimes==-60 & allDat_UNI_actual$tree==1),] #15:24
slopes35_t1<-allDat_UNI_actual[which(allDat_UNI_actual$sliceTimes==-35 & allDat_UNI_actual$tree==1),]
slopes10_t1<-allDat_UNI_actual[which(allDat_UNI_actual$sliceTimes==-10 & allDat_UNI_actual$tree==1),]
slopes_t1<-list(slopes60_t1,slopes35_t1,slopes10_t1)


#SIMS_allDat_UNI_actual<-vector("list",length=3)
#for (i in 1:length(SIMS_allDat_UNI_actual)){
#	SIMS_allDat_UNI_actual[[i]]<-read.table(paste("MamPhy_SIMS_",sims[i],"_",bbone,"_sample100_ALL100-TABLE_PGLS-UNI_timeSlices_ACTUAL_all.txt",sep=""))
#}
#
## you need all vals, incl intercepts and Ps
#sim<-SIMS_allDat_UNI_actual[[1]]
#SIMslopes60<-sim[which(sim$sliceTimes==-60),] #12:21
#SIMslopes35<-sim[which(sim$sliceTimes==-35),]
#SIMslopes10<-sim[which(sim$sliceTimes==-10),]
#SIMslopes<-list(SIMslopes60,SIMslopes35,SIMslopes10)

# Get the SIM actual slopes as a RATIO of the empirical actual vs standard...
# sims scaled
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/")
load(file="cladeSlice_workspace_forPlots.Rdata")
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_SIMULATIONS_on_NDexp")
SIMS_allDat_UNI_scaled<-vector("list",length=3)
for (i in 1:length(SIMS_allDat_MULTI_actual)){
	SIMS_allDat_UNI_scaled[[i]]<-read.table(paste("MamPhy_SIMS_",sims[i],"_",bbone,"_sample100_ALL100-TABLE_PGLS-UNI_timeSlices_SCALED_all.txt",sep=""))
}
Sim1_scaled<-SIMS_allDat_UNI_scaled[[1]]

# empirical scaled
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/")
allDat_UNI_i<-vector("list",length=numTrees)
for(i in 1:numTrees){
	tree<-rep(i,length(slice))
	dat<-read.table(paste(bbone,"_sample100_",i,"_PGLS_timeSlices_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))
	allDat_UNI_i[[i]]<-cbind(sliceTimes,dat,slice,tree)
}
allDat_UNI_scaled<-do.call(rbind,allDat_UNI_i)

# now do the RATIOS
simActual_MRCAs<-vector("list",length=3)
simActual_i_1s<-vector("list",length=3)
for(i in 1:3){
	for (j in 1:length(allDat_UNI_actual[,"MRCA"])){
		simActual_MRCAs[[i]][j]<-(allDat_UNI_actual[j,"MRCA"] / allDat_UNI_scaled[j,"MRCA"])*SIMS_allDat_UNI_scaled[[i]][j,"MRCA"]
		simActual_i_1s[[i]][j]<-(allDat_UNI_actual[j,"i_1"] / allDat_UNI_scaled[j,"i_1"])*SIMS_allDat_UNI_scaled[[i]][j,"i_1"]
	}
}

simActual_DRharms<-vector("list",length=3)
simActual_i_2s<-vector("list",length=3)
for(i in 1:3){
	for (j in 1:length(allDat_UNI_actual[,"DR_harm"])){
		simActual_DRharms[[i]][j]<-(allDat_UNI_actual[j,"DR_harm"] / allDat_UNI_scaled[j,"DR_harm"])*SIMS_allDat_UNI_scaled[[i]][j,"DR_harm"]
		simActual_i_2s[[i]][j]<-(allDat_UNI_actual[j,"i_2"] / allDat_UNI_scaled[j,"i_2"])*SIMS_allDat_UNI_scaled[[i]][j,"i_2"]
	}
}

simActual_DRskews<-vector("list",length=3)
simActual_i_3s<-vector("list",length=3)
for(i in 1:3){
	for (j in 1:length(allDat_UNI_actual[,"DR_skew"])){
		simActual_DRskews[[i]][j]<-(allDat_UNI_actual[j,"DR_skew"] / allDat_UNI_scaled[j,"DR_skew"])*SIMS_allDat_UNI_scaled[[i]][j,"DR_skew"]
		simActual_i_3s[[i]][j]<-(allDat_UNI_actual[j,"i_3"] / allDat_UNI_scaled[j,"i_3"])*SIMS_allDat_UNI_scaled[[i]][j,"i_3"]
	}
}
#simActual_MRCA  compare to -- sim$MRCA
xx<-simActual_MRCA-sim$MRCA
hist(xx) # centered around 0, some pos, some neg

# combine them...
simActualCorrect_all<-vector("list",length=3)
simActualsToPlot_100<-vector("list",length=3)
simActualsToPlot_means<-vector("list",length=3)
for(i in 1:3){
	simActualCorrect_all[[i]]<-cbind.data.frame(SIMS_allDat_UNI_scaled[[i]][,"sliceTimes"],simActual_i_1s[[i]],simActual_i_2s[[i]],simActual_i_3s[[i]],simActual_MRCAs[[i]],simActual_DRharms[[i]],simActual_DRskews[[i]],SIMS_allDat_UNI_scaled[[i]][,"p_1"],SIMS_allDat_UNI_scaled[[i]][,"p_2"],SIMS_allDat_UNI_scaled[[i]][,"p_3"])
	colnames(simActualCorrect_all[[i]])<-c("sliceTimes","i_1","i_2","i_3","MRCA","DR_harm","DR_skew","p_1","p_2","p_3")
	sim<-simActualCorrect_all[[i]]
	simActualsToPlot_100[[i]]<-list(sim[which(sim$sliceTimes==-60),],sim[which(sim$sliceTimes==-35),],sim[which(sim$sliceTimes==-10),])

	means<-data.frame(matrix(NA, nrow = 3, ncol = 10))
	colnames(means)<-colnames(simActualCorrect_all[[i]])
for (j in 1:3){
	dat<-simActualsToPlot_100[[i]][[j]]
	means[j,1:10]<-colMeans(dat)		
	simActualsToPlot_means[[i]]<-means
}
}

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")
#library(plotrix)
# PLOT -- SLICES - 60, 35, 10 Ma -- empirical data for 1 tree ONLY...
##
for (j in 1:3){

n<-c(60,35,10)
dat<-dats_t1[[j]]
slopeDat<-slopes_t1[[j]]
#SIMdat<-SIMdats[[j]]
SIMslopeDat<-simActualsToPlot_means[[j]]

pdf(file=paste("cladeLevel",bbone,"_scatter_timeSlices",n[j],"_PGLSuni_tree1only_PLOTTED_3part_withSIMmeans.pdf",sep=""),width=1.5,height=3.674)

par(mfrow = c(3,1),oma = c(1,1,0,0) + 0.1, mar = c(0.5,1,1,1) + 0.1)
ltys<-1
ltySim<-2
lwds<-2
cexAxis<-0.8
CexPt<-1
LwdPt<-0.6
logs<-c(5,20,50,200,500,2000)
colSims<-c("darkgoldenrod2","deepskyblue2","darkorchid3")
colEmpir<-"black"
colPts<-grey(0.3)

#plot(log(richness)~MRCA,data=dat60, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
#smoothScatter(log(dat60$richness)~dat60$MRCA,xlab="",ylab="",nbin=50,nrpoints=0,colramp=colorRampPalette(c("white", "black")))
#plot((richness)~MRCA,data=dat,type="p",log="y",col="black",xlab="", ylab="log(richness)",cex.axis=cexAxis,cex.lab=1.5,font.lab=2)

#quartz(,width=1.5,height=3.674)
#par(mfrow = c(3,1),oma = c(1,1,0,0) + 0.1, mar = c(0.5,1,1,1) + 0.1)
###
# MRCA
plot(log(richness)~MRCA,data=dat,type="n",col="black",xlab="", ylab="log(richness)",yaxt="n",xaxt="n",cex.axis=cexAxis,cex.lab=1.5,font.lab=2)
axis(side=2,at=log(logs),labels=logs, cex.axis=cexAxis, font=1)
axis(side=1,at=NULL,labels=TRUE, cex.axis=cexAxis, font=1, padj=-1)

# sims
#points(log(richness)~MRCA,data=SIMdat,type="p",col="darkgoldenrod2",xlab="", ylab="",cex.lab=1.5,font.lab=2)

# empirical
points(log(richness)~MRCA,data=dat,type="p",col=colPts,xlab="", ylab="",lwd=LwdPt,cex=CexPt,cex.lab=1.5,font.lab=2)

#lines
for (i in 1:length(SIMslopeDat[,1])){
	a=(SIMslopeDat$i_1[i])
	b=(SIMslopeDat$MRCA[i])
	X=seq(min(dat$MRCA),max(dat$MRCA),length.out=20)
	lines(x=X,y=b*X+a,col=colSims[i], lwd=lwds,lty=ltySim) #x1=min(dat$MRCA),x2=max(dat$MRCA),y1=min(log(dat$richness)),y2=max(log(dat$richness)),untf=FALSE,
}
a=(slopeDat$i_1)
b=(slopeDat$MRCA)
X=seq(min(dat$MRCA),max(dat$MRCA),length.out=20)
lines(x=X,y=b*X+a,col=colEmpir, lwd=lwds,lty=ltys) #ablineclip(a=(slopeDat$i_1),b=(slopeDat$MRCA), x1=min(dat$MRCA),x2=max(dat$MRCA),y1=min(log(dat$richness)),y2=max(log(dat$richness)),untf=FALSE,


###
# DR_mean
plot(log(richness)~DR_harm,data=dat,type="n",col="black",xlab="", ylab="log(richness)",yaxt="n",xaxt="n",cex.axis=cexAxis,cex.lab=1.5,font.lab=2)
axis(side=2,at=log(logs),labels=logs, cex.axis=cexAxis, font=1)
axis(side=1,at=NULL,labels=TRUE, cex.axis=cexAxis, font=1, padj=-1)

# sims
#points(log(richness)~MRCA,data=SIMdat,type="p",col="darkgoldenrod2",xlab="", ylab="",cex.lab=1.5,font.lab=2)

# empirical
points(log(richness)~DR_harm,data=dat,type="p",col=colPts,xlab="", ylab="",lwd=LwdPt,cex=CexPt,cex.lab=1.5,font.lab=2)

#lines
for (i in 1:length(SIMslopeDat[,1])){
	a=(SIMslopeDat$i_2[i])
	b=(SIMslopeDat$DR_harm[i])
	X=seq(min(dat$DR_harm),max(dat$DR_harm),length.out=20)
	lines(x=X,y=b*X+a,col=colSims[i], lwd=lwds,lty=ltySim) #x1=min(dat$MRCA),x2=max(dat$MRCA),y1=min(log(dat$richness)),y2=max(log(dat$richness)),untf=FALSE,
}
a=(slopeDat$i_2)
b=(slopeDat$DR_harm)
X=seq(min(dat$DR_harm),max(dat$DR_harm),length.out=20)
lines(x=X,y=b*X+a,col=colEmpir, lwd=lwds,lty=ltys) #ablineclip(a=(slopeDat$i_1),b=(slopeDat$MRCA), x1=min(dat$MRCA),x2=max(dat$MRCA),y1=min(log(dat$richness)),y2=max(log(dat$richness)),untf=FALSE,


###
# DR_skew
plot(log(richness)~DR_skew,data=dat,type="n",col="black",xlab="", ylab="log(richness)",yaxt="n",xaxt="n",cex.axis=cexAxis,cex.lab=1.5,font.lab=2)
axis(side=2,at=log(logs),labels=logs, cex.axis=cexAxis, font=1)
axis(side=1,at=NULL,labels=TRUE, cex.axis=cexAxis, font=1, padj=-1)

# sims
#points(log(richness)~MRCA,data=SIMdat,type="p",col="darkgoldenrod2",xlab="", ylab="",cex.lab=1.5,font.lab=2)

# empirical
points(log(richness)~DR_skew,data=dat,type="p",col=colPts,xlab="", ylab="",lwd=LwdPt,cex=CexPt,cex.lab=1.5,font.lab=2)

#lines
for (i in 1:length(SIMslopeDat[,1])){
	a=(SIMslopeDat$i_3[i])
	b=(SIMslopeDat$DR_skew[i])
	X=seq(min(na.omit(dat$DR_skew)),max(na.omit(dat$DR_skew)),length.out=20)
	lines(x=X,y=b*X+a,col=colSims[i], lwd=lwds,lty=ltySim) #x1=min(dat$MRCA),x2=max(dat$MRCA),y1=min(log(dat$richness)),y2=max(log(dat$richness)),untf=FALSE,
}
a=(slopeDat$i_3)
b=(slopeDat$DR_skew)
X=seq(min(na.omit(dat$DR_skew)),max(na.omit(dat$DR_skew)),length.out=20)
lines(x=X,y=b*X+a,col=colEmpir, lwd=lwds,lty=ltys) #ablineclip(a=(slopeDat$i_1),b=(slopeDat$MRCA), x1=min(dat$MRCA),x2=max(dat$MRCA),y1=min(log(dat$richness)),y2=max(log(dat$richness)),untf=FALSE,

title(main="", xlab = "",
      ylab = "",
      outer = TRUE, line = 1,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)
dev.off()

}




setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/")
load(file="cladeSlice_workspace_forPlots.Rdata")
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")

# PLOT -- SLICES - 60, 35, 10 Ma -- WITH DATA POINTS, 100 trees
##

for (j in 1:3){

n<-c(60,35,10)
dat<-dats[[j]]
slopeDat<-slopes[[j]]
#SIMdat<-SIMdats[[j]]
SIMslopeDat<-simActualSlopes[[j]]

#jpeg(file=paste("cladeLevel",bbone,"_scatter_timeSlices",n[j],"_PGLSuni_ALL_PLOTTED_3part_withSIM.jpg",sep=""),width=4,height=8, units="in", res=600, quality=100)
pdf(file=paste("cladeLevel",bbone,"_scatter_timeSlices",n[j],"_PGLSuni_ALL_PLOTTED_3part_withSIM.pdf",sep=""),width=1.5,height=3.674)

par(mfrow = c(3,1),oma = c(0,0,0,0) + 0.1, mar = c(0.5,1,1,1) + 0.1)
ltys<-1
lwds<-4
cexAxis<-0.8
CexPt<-0.6
LwdPt<-0.3
colRaw<-col2rgb("light blue")/255
cols<-rgb(colRaw[1],colRaw[2],colRaw[3],alpha=0.5)
colRaw2<-col2rgb("darkgoldenrod2")/255
cols2<-rgb(colRaw2[1],colRaw2[2],colRaw2[3],alpha=0.5)
logs<-c(5,20,50,200,500,2000)
#plot(log(richness)~MRCA,data=dat60, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
#smoothScatter(log(dat60$richness)~dat60$MRCA,xlab="",ylab="",nbin=50,nrpoints=0,colramp=colorRampPalette(c("white", "black")))
#plot((richness)~MRCA,data=dat,type="p",log="y",col="black",xlab="", ylab="log(richness)",cex.axis=cexAxis,cex.lab=1.5,font.lab=2)

###
# MRCA
plot(log(richness)~MRCA,data=dat,type="n",col="black",xlab="", ylab="log(richness)",yaxt="n",xaxt="n",cex.axis=cexAxis,cex.lab=1.5,font.lab=2)
axis(side=2,at=log(logs),labels=logs, cex.axis=cexAxis, font=1)
axis(side=1,at=NULL,labels=TRUE, cex.axis=cexAxis, font=1, padj=-1)

# sims
#points(log(richness)~MRCA,data=SIMdat,type="p",col="darkgoldenrod2",xlab="", ylab="",cex.lab=1.5,font.lab=2)

# empirical
points(log(richness)~MRCA,data=dat,type="p",col=grey(0.6),xlab="", ylab="",lwd=LwdPt,cex=CexPt,cex.lab=1.5,font.lab=2)

#lines
for (i in 1:length(SIMslopeDat[,1])){
	abline(a=(SIMslopeDat$i_1[i]),b=(SIMslopeDat$MRCA[i]), untf=FALSE,col=cols2, lwd=lwds,lty=ltys)
}
for (i in 1:length(SIMslopeDat[,1])){
	if (SIMslopeDat$p_1[i] > 0.05) {
	abline(a=(SIMslopeDat$i_1[i]),b=(SIMslopeDat$MRCA[i]), untf=FALSE,col=rgb(1,0,0,alpha=0.3), lwd=lwds,lty=ltys)		
	} else next
}
for (i in 1:length(SIMslopeDat[,1])){
	abline(a=(SIMslopeDat$i_1[i]),b=(SIMslopeDat$MRCA[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}

# lines
for (i in 1:length(slopeDat[,1])){
	abline(a=(slopeDat$i_1[i]),b=(slopeDat$MRCA[i]), untf=FALSE,col=rgb(0,0,0,alpha=0.2), lwd=lwds,lty=ltys)
}
for (i in 1:length(slopeDat[,1])){
	if (slopeDat$p_1[i] > 0.05) {
	abline(a=(slopeDat$i_1[i]),b=(slopeDat$MRCA[i]), untf=FALSE,col=rgb(1,0,0,alpha=0.3), lwd=lwds,lty=ltys)		
	} else next
}
for (i in 1:length(slopeDat[,1])){
	abline(a=(slopeDat$i_1[i]),b=(slopeDat$MRCA[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}


####
# DR mean
plot(log(richness)~DR_harm,data=dat,type="n",col="black",xlab="", ylab="log(richness)",yaxt="n",xaxt="n",cex.axis=cexAxis,cex.lab=1.5,font.lab=2)
axis(side=2,at=log(logs),labels=logs, cex.axis=cexAxis, font=1)
axis(side=1,at=NULL,labels=TRUE, cex.axis=cexAxis, font=1, padj=-1)

# sims
#points(log(richness)~DR_harm,data=SIMdat,type="p",col="darkgoldenrod2",xlab="", ylab="",cex.lab=1.5,font.lab=2)

# empirical
points(log(richness)~DR_harm,data=dat,type="p",col=grey(0.6),xlab="", ylab="",lwd=LwdPt,cex=CexPt,cex.lab=1.5,font.lab=2)
# lines
for (i in 1:length(SIMslopeDat[,1])){
	abline(a=(SIMslopeDat$i_2[i]),b=(SIMslopeDat$DR_harm[i]), untf=FALSE,col=cols2, lwd=lwds,lty=ltys)
}
for (i in 1:length(SIMslopeDat[,1])){
	if (SIMslopeDat$p_2[i] > 0.05) {
	abline(a=(SIMslopeDat$i_2[i]),b=(SIMslopeDat$DR_harm[i]), untf=FALSE,col=rgb(1,0,0,alpha=0.3), lwd=lwds,lty=ltys)		
	} else next
}
for (i in 1:length(SIMslopeDat[,1])){
	abline(a=(SIMslopeDat$i_2[i]),b=(SIMslopeDat$DR_harm[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}
# lines
for (i in 1:length(slopeDat[,1])){
	abline(a=(slopeDat$i_2[i]),b=(slopeDat$DR_harm[i]), untf=FALSE,col=rgb(0,0,0,alpha=0.2), lwd=lwds,lty=ltys)
}
for (i in 1:length(slopeDat[,1])){
	if (slopeDat$p_2[i] > 0.05) {
	abline(a=(slopeDat$i_2[i]),b=(slopeDat$DR_harm[i]), untf=FALSE,col=rgb(1,0,0,alpha=0.3), lwd=lwds,lty=ltys)		
	} else next
}
for (i in 1:length(slopeDat[,1])){
	abline(a=(slopeDat$i_2[i]),b=(slopeDat$DR_harm[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}


####
# DR skew
plot(log(richness)~DR_skew,data=dat,type="n",col="black",xlab="", ylab="log(richness)",yaxt="n",xaxt="n",cex.axis=cexAxis,cex.lab=1.5,font.lab=2)
axis(side=2,at=log(logs),labels=logs, cex.axis=cexAxis, font=1)
axis(side=1,at=NULL,labels=TRUE, cex.axis=cexAxis, font=1, padj=-1)

# sims
#points(log(richness)~DR_skew,data=SIMdat,type="p",col="darkgoldenrod2",xlab="", ylab="",cex.lab=1.5,font.lab=2)
# empirical
points(log(richness)~DR_skew,data=dat,type="p",col=grey(0.6),xlab="", ylab="",lwd=LwdPt,cex=CexPt,cex.lab=1.5,font.lab=2)

# lines
for (i in 1:length(SIMslopeDat[,1])){
	abline(a=(SIMslopeDat$i_3[i]),b=(SIMslopeDat$DR_skew[i]), untf=FALSE,col=cols2, lwd=lwds,lty=ltys)
}
for (i in 1:length(SIMslopeDat[,1])){
	if (SIMslopeDat$p_3[i] > 0.05) {
	abline(a=(SIMslopeDat$i_3[i]),b=(SIMslopeDat$DR_skew[i]), untf=FALSE,col=rgb(1,0,0,alpha=0.3), lwd=lwds,lty=ltys)		
	} else next
}
for (i in 1:length(SIMslopeDat[,1])){
	abline(a=(SIMslopeDat$i_3[i]),b=(SIMslopeDat$DR_skew[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}# lines
for (i in 1:length(slopeDat[,1])){
	abline(a=(slopeDat$i_3[i]),b=(slopeDat$DR_skew[i]), untf=FALSE,col=rgb(0,0,0,alpha=0.2), lwd=lwds,lty=ltys)
}
for (i in 1:length(slopeDat[,1])){
	if (slopeDat$p_3[i] > 0.05) {
	abline(a=(slopeDat$i_3[i]),b=(slopeDat$DR_skew[i]), untf=FALSE,col=rgb(1,0,0,alpha=0.3), lwd=lwds,lty=ltys)		
	} else next
}
for (i in 1:length(slopeDat[,1])){
	abline(a=(slopeDat$i_3[i]),b=(slopeDat$DR_skew[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}

title(main="", xlab = "",
      ylab = "",
      outer = TRUE, line = 1,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)
dev.off()

}
# ==============
save.image(file="cladeSlice_workspace_forPlots.Rdata")


setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/")
load(file="cladeSlice_workspace_forPlots.Rdata")



#mtext(text="Clade crown age (Ma)", side=1, padj=3, font=2,cex=1)

#plot(log(richness)~DR_harm,data=dat60, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
smoothScatter(log(dat60$richness)~dat60$DR_harm,xlab="",ylab="",nbin=50,nrpoints=0,colramp=colorRampPalette(c("white", "black")))
for (i in 1:length(uniDat_60[,1])){
	abline(a=(uniDat_60$i_2[i]),b=(uniDat_60$DR_harm[i]), untf=FALSE,col=cols, lwd=lwds,lty=ltys)
}
for (i in 1:length(uniDat_60[,1])){
	abline(a=(uniDat_60$i_2[i]),b=(uniDat_60$DR_harm[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}
#mtext(text="Clade ESDR mean (species / Ma)", side=1, padj=3, font=2,cex=1)

#plot(log(richness)~DR_skew,data=dat60, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
smoothScatter(log(dat60$richness)~dat60$DR_skew,xlab="",ylab="",nbin=50,nrpoints=0,colramp=colorRampPalette(c("white", "black")))
for (i in 1:length(uniDat_60[,1])){
	abline(a=(uniDat_60$i_3[i]),b=(uniDat_60$DR_skew[i]), untf=FALSE,col=cols, lwd=lwds,lty=ltys)
}
for (i in 1:length(uniDat_60[,1])){
	abline(a=(uniDat_60$i_3[i]),b=(uniDat_60$DR_skew[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}
#mtext(text="Clade ESDR skewness", side=1, padj=3, font=2,cex=1)

title(main="", xlab = "",
      ylab = "",
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)
dev.off()

}




# 10 Ma slice
##
#pdf(file=paste("cladeLevel",bbone,"_scatter_timeSlices10_PGLSuni_Yaxis_PLOTTED_3part.pdf",sep=""),onefile=TRUE,width=4,height=10)
jpeg(file=paste("cladeLevel",bbone,"_scatter_timeSlices10_PGLSuni_Yaxis_PLOTTED_3part.jpg",sep=""),width=4,height=8, units="in", res=600, quality=100)
par(mfrow = c(3,1),oma = c(5,4,5,0) + 0.1, mar = c(3,1,1,1) + 0.1)
# type="n",
plot((richness)~MRCA,data=dat10, log="y",type="n",col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
mtext(text="Clade crown age (Ma)", side=1, padj=3, font=2,cex=1)

plot((richness)~DR_harm,data=dat10, log="y",type="n",col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
mtext(text="Clade ESDR mean (species / Ma)", side=1, padj=3, font=2,cex=1)

plot((richness)~DR_skew,data=dat10, log="y", type="n",col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
mtext(text="Clade ESDR skewness", side=1, padj=3, font=2,cex=1)

title(main="", xlab = "",
      ylab = "log(richness)",
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)
dev.off()


#pdf(file=paste("cladeLevel",bbone,"_scatter_timeSlices10_PGLSuni_ALL_PLOTTED_3part.pdf",sep=""),onefile=TRUE,width=4,height=10)
jpeg(file=paste("cladeLevel",bbone,"_scatter_timeSlices10_PGLSuni_ALL_PLOTTED_3part_new.jpg",sep=""),width=4,height=8, units="in", res=600, quality=100)
par(mfrow = c(3,1),oma = c(5,4,5,0) + 0.1, mar = c(3,1,1,1) + 0.1)
ltys<-1
lwds<-4
colRaw<-col2rgb("light blue")/255
cols<-rgb(colRaw[1],colRaw[2],colRaw[3],alpha=0.5)
#plot(log(richness)~MRCA,data=dat10, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
smoothScatter(log(dat10$richness)~dat10$MRCA,xlab="",ylab="",nbin=50,nrpoints=0,colramp=colorRampPalette(c("white", "black")))
for (i in 1:length(uniDat_10[,1])){
	abline(a=(uniDat_10$i_1[i]),b=(uniDat_10$MRCA[i]), untf=FALSE,col=cols, lwd=lwds,lty=ltys)
}
for (i in 1:length(uniDat_10[,1])){
	abline(a=(uniDat_10$i_1[i]),b=(uniDat_10$MRCA[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}

#mtext(text="Clade crown age (Ma)", side=1, padj=3, font=2,cex=1)

#plot(log(richness)~DR_harm,data=dat10, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
smoothScatter(log(dat10$richness)~dat10$DR_harm,xlab="",ylab="",nbin=50,nrpoints=0,colramp=colorRampPalette(c("white", "black")))
for (i in 1:length(uniDat_10[,1])){
	abline(a=(uniDat_10$i_2[i]),b=(uniDat_10$DR_harm[i]), untf=FALSE,col=cols, lwd=lwds,lty=ltys)
}
for (i in 1:length(uniDat_10[,1])){
	abline(a=(uniDat_10$i_2[i]),b=(uniDat_10$DR_harm[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}
#mtext(text="Clade ESDR mean (species / Ma)", side=1, padj=3, font=2,cex=1)

#plot(log(richness)~DR_skew,data=dat10, col="light grey",xlab="", ylab="log(richness)",cex.lab=1.5,font.lab=2)
smoothScatter(log(dat10$richness)~dat10$DR_skew,xlab="",ylab="",nbin=50,nrpoints=0,colramp=colorRampPalette(c("white", "black")))
for (i in 1:length(uniDat_10[,1])){
	abline(a=(uniDat_10$i_3[i]),b=(uniDat_10$DR_skew[i]), untf=FALSE,col=cols, lwd=lwds,lty=ltys)
}
for (i in 1:length(uniDat_10[,1])){
	abline(a=(uniDat_10$i_3[i]),b=(uniDat_10$DR_skew[i]), untf=FALSE,col="black", lwd=0.1,lty=1)
}
#mtext(text="Clade ESDR skewness", side=1, padj=3, font=2,cex=1)

title(main="", xlab = "",
      ylab = "",
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)
dev.off()




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

















# ======================

#######################################

quartz(width=10,height=4)
par(mfrow=c(1,3))
plot(mamPhy, show.tip.label=FALSE)
axisPhylo()
plot(ladderize(sim), show.tip.label=FALSE)
axisPhylo()
plot(ladderize(simLoE), show.tip.label=FALSE)
axisPhylo()



setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut")


########################
##############

slicePhys10<-vector("list",length=numTrees)
slicePhys60<-vector("list",length=numTrees)
for (i in 1:numTrees){
	slicePhys<-read.tree(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_slicePhy-5to70Ma_by5.trees",sep=""))
	slicePhys10[[i]]<-slicePhys[[2]]
	slicePhys60[[i]]<-slicePhys[[12]]
}

dat<-allRes_100trees[which(allRes_100trees$slice==10),]

m1<-glmmPQL(fixed = richness ~ MRCA, random = ~ 1 | tree, correlation=corBrownian(value=1,phy=slicePhys10,form=~1|tree),family = poisson, data = dat)
sum1<-summary(m1)
sum1


# make vcv matrices
treeCorMats<-vector("list",length(slicePhys10))
nrows<-vector("list",length(slicePhys10))
labs<-vector("list",length(slicePhys10))
for(j in 1:100){
	treeCorMats[[j]]<-vcv(slicePhys10[[j]],corr=TRUE)
	nrows[[j]]<-length(slicePhys10[[j]]$tip.label)
	labs[[j]]<-slicePhys10[[j]]$tip.label
}
nrows_all<-sum(do.call(rbind,nrows)[,1])
labs_all<-as.data.frame(unlist(labs))

comboCorrs<-data.frame(matrix(0, nrow = nrows_all, ncol = nrows_all),row.names=labs_all[,1])
colnames(comboCorrs)<-labs_all[,1]
head(comboCorrs[,1:10])

# ^ would be 65k * 65k == 4 billion cells!  That aint gonna work...

dims<-vector()
dims[1]<-dim(sliceCorMats[[1]])[1]
for (j in 2:length(slicePhys)){
	dims[j]<-dim(sliceCorMats[[j]])[1]+dims[j-1]
}

comboCorrs[1:dims[1],1:dims[1]]<-sliceCorMats[[1]]
for (j in 2:length(sliceCorMats)){
	comboCorrs[(dims[j-1]+1):dims[j],(dims[j-1]+1):dims[j]]<-sliceCorMats[[j]]
}


fit1<-gls(richness ~ MRCA, correlation=corSymm(comboCorrs[lower.tri(comboCorrs)],fixed=TRUE),data=dat)
#fit2<-gls(richness ~ MRCA, correlation=corSymm(comboCorrs,fixed=TRUE),data=dat) # reproduced the 'wrong dimension' error...
sum<-summary(fit1)
sum$tTable ## THIS WORKS.

xx<-corSymm(comboCorrs[lower.tri(comboCorrs)],fixed=TRUE)
qq<-Initialize(xx,dat) # THIS WORKS.


######
# Effect of age on DR -- but needs to be PHYLOGENETIC in its correction.

fit<-lme(fixed=DR_harm~MRCA,data=dat10, random= ~1 | tree)






log="y",
ylim=c(0,max(dat10$richness))
ylim=c(1,125),


jpeg(file=paste("cladeLevel",bbone,"_scatter_timeSlices10and60_PLOTTED_3part.jpg",sep=""),width=9,height=6,units="in", res=300, quality=100)

par(mfrow = c(2,3))#,oma = c(5,4,5,0) + 0.1, mar = c(2,1,1,1) + 0.1)

plot(richness ~ MRCA, data=dat10, log="y",col=grey(0.3,alpha=0.5), ylab="log(richness)",xlab="")# xlim=rev(range(MRCA)), ,ylab="Partial residual", xlab="Slice times (Ma)",main="")
mtext(text="(a) 10 Ma slice", adj=0, font=2)#at=c(-72,0.21)
plot((richness) ~ DR_harm, data=dat10, log="y",col=grey(0.3,alpha=0.5), ylab="log(richness)",xlab="")# xlim=rev(range(MRCA)), ,ylab="Partial residual", xlab="Slice times (Ma)",main="")
plot((richness) ~ DR_skew, data=dat10, log="y",col=grey(0.3,alpha=0.5), ylab="log(richness)",xlab="")# xlim=rev(range(MRCA)), ,ylab="Partial residual", xlab="Slice times (Ma)",main="")


plot((richness) ~ MRCA, data=dat60, log="y",col=grey(0.3,alpha=0.5), ylab="log(richness)",xlab="Crown age (Ma)")# xlim=rev(range(MRCA)), ,ylab="Partial residual", xlab="Slice times (Ma)",main="")
mtext(text="(b) 60 Ma slice", adj=0, font=2)#at=c(-72,0.21)
plot((richness) ~ DR_harm, data=dat60, log="y",col=grey(0.3,alpha=0.5), ylab="log(richness)",xlab="ESDR mean")# xlim=rev(range(MRCA)), ,ylab="Partial residual", xlab="Slice times (Ma)",main="")
plot((richness) ~ DR_skew, data=dat60, log="y",col=grey(0.3,alpha=0.5), ylab="log(richness)",xlab="ESDR skewness")# xlim=rev(range(MRCA)), ,ylab="Partial residual", xlab="Slice times (Ma)",main="")

#title(main="", xlab = "",     ylab = "log(richness)", outer = TRUE, line = 3,cex.main=1,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)

dev.off()




########

# read back in slice phys
slicePhys<-vector("list",length=14)
slicePhys<-read.tree(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_slicePhy-5to70Ma_by5.trees",sep=""))

# make vcv matrices
sliceCorMats<-vector("list",length(slicePhys))
nrows<-vector("list",length(slicePhys))
labs<-vector("list",length(slicePhys))
for(j in 1:length(slicePhys)){
	sliceCorMats[[j]]<-vcv(slicePhys[[j]],corr=TRUE)
	nrows[[j]]<-length(slicePhys[[j]]$tip.label)
	labs[[j]]<-slicePhys[[j]]$tip.label
}
nrows_all<-sum(do.call(rbind,nrows)[,1])
labs_all<-as.data.frame(unlist(labs))

comboCorrs<-data.frame(matrix(0, nrow = nrows_all, ncol = nrows_all),row.names=labs_all[,1])
colnames(comboCorrs)<-labs_all[,1]
head(comboCorrs[,1:10])

dims<-vector()
dims[1]<-dim(sliceCorMats[[1]])[1]
for (j in 2:length(slicePhys)){
	dims[j]<-dim(sliceCorMats[[j]])[1]+dims[j-1]
}

comboCorrs[1:dims[1],1:dims[1]]<-sliceCorMats[[1]]
for (j in 2:length(sliceCorMats)){
	comboCorrs[(dims[j-1]+1):dims[j],(dims[j-1]+1):dims[j]]<-sliceCorMats[[j]]
}

dat<-allRes_100trees[which(allRes_100trees$tree==1),]
fit1<-gls(richness ~ MRCA, correlation=corSymm(comboCorrs[lower.tri(comboCorrs)],fixed=TRUE),data=dat)
#fit2<-gls(richness ~ MRCA, correlation=corSymm(comboCorrs,fixed=TRUE),data=dat) # reproduced the 'wrong dimension' error...
sum<-summary(fit1)
sum$tTable ## THIS WORKS.

xx<-corSymm(comboCorrs[lower.tri(comboCorrs)],fixed=TRUE)
qq<-Initialize(xx,dat) # THIS WORKS.





#########
#####
# Code to get the WITHIN tree results across slices-- pre PGLS
# Prior to doing this across all 100...
###
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut")
bbone<- "NDexp" #"FBD" # 
library(nlme); library(MASS)

numTrees<-100
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery

allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

allRes<-vector("list",length=numTrees)
for(i in 1:numTrees){
results<-vector("list",length(allCladeSetNames))
for (j in 1:length(allCladeSetNames)){
	results[[j]]<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""))
	rownames(results[[j]])<-paste(i,"_",j,"_",c(1:length(results[[j]][,1])),sep="")
}
allRes[[i]]<-do.call(rbind,results) # each is 19 vars * 3000 entries about
#write.table(allRes,file=paste(bbone,"_sample100_",i,"_slice5Ma-to-70Ma_cladeSTATS.txt",sep=""))
}
allRes_100trees<-do.call(rbind,allRes) # 304274 obs. of  19 variables:

save(allRes_100trees, file=paste(bbone,"_sample100__ALL100trees_slice5Ma-to-70Ma_cladeSTATS.Rdata",sep=""))

load(file=paste(bbone,"_sample100__ALL100trees_slice5Ma-to-70Ma_cladeSTATS.Rdata",sep=""))





m1<-glmmPQL(fixed = richness ~ MRCA, random = ~ 1 | slice | tree, family = poisson, data = allRes_100trees)
sum1<-summary(m1)
sum1

m2<-glmmPQL(fixed = richness ~ MRCA + DR_harm, random = ~ 1 | slice | tree, family = poisson, data = allRes_100trees)
sum2<-summary(m2)
sum2


m2a<-glmmPQL(fixed = richness ~ DR_harm + MRCA, random = ~ 1 | slice | tree, family = poisson, data = allRes_100trees)
sum2a<-summary(m2a)
sum2a


######
#===================
# Now with PHYLOGENETIC correlation structures. Specific to each tree/slice combination.
#==================
# CREATE the slice phys, per tree
library(ape)
library(picante)
library(phytools)
library(geiger)

bbone<- "NDexp" #"FBD" # 
setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_",bbone,"_nexus-and-newickTrees",sep=""))

i=1
# load in tree
mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_nexus.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 


## Make time slices to create the clades
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery

root=max(node.age(mamPhy)$ages)

allCladeSets<-vector("list",length=numSlices)
allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSets[[j]]<-treeSlice(mamPhy, slice=root-(sliceEvery*j), trivial=FALSE, prompt=FALSE)
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

# create slice phys
slicePhys<-vector("list",length(allCladeSets))
for (k in 1:length(allCladeSets)){
cladeReps<-vector()
for (j in 1:length(allCladeSets[[k]])){
	cladeSp<-allCladeSets[[k]][[j]]$tip.label
	cladeReps[j]<-cladeSp[1]
	}
toDrop<-setdiff(mamPhy$tip.label,cladeReps)
slicePhys[[k]]<-drop.tip(mamPhy,toDrop)
}

# write slice phys
for(j in 1:length(slicePhys)){
	slicePhys[[j]]$tip.label<-paste(i,"_",j,"_",c(1:length(slicePhys[[j]]$tip.label)),sep="")
	write.tree(slicePhys[[j]], file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_slicePhy-5to70Ma_by5.trees",sep=""), append=TRUE)
}

# read back in slice phys
slicePhys<-vector("list",length=14)
slicePhys<-read.tree(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_slicePhy-5to70Ma_by5.trees",sep=""))

# make vcv matrices
sliceCorMats<-vector("list",length(slicePhys))
nrows<-vector("list",length(slicePhys))
labs<-vector("list",length(slicePhys))
for(j in 1:length(slicePhys)){
	sliceCorMats[[j]]<-vcv(slicePhys[[j]],corr=TRUE)
	nrows[[j]]<-length(slicePhys[[j]]$tip.label)
	labs[[j]]<-slicePhys[[j]]$tip.label
}
nrows_all<-sum(do.call(rbind,nrows)[,1])
labs_all<-as.data.frame(unlist(labs))

comboCorrs<-data.frame(matrix(0, nrow = nrows_all, ncol = nrows_all),row.names=labs_all[,1])
colnames(comboCorrs)<-labs_all[,1]
head(comboCorrs[,1:10])

dims<-vector()
dims[1]<-dim(sliceCorMats[[1]])[1]
for (j in 2:length(slicePhys)){
	dims[j]<-dim(sliceCorMats[[j]])[1]+dims[j-1]
}

comboCorrs[1:dims[1],1:dims[1]]<-sliceCorMats[[1]]
for (j in 2:length(sliceCorMats)){
	comboCorrs[(dims[j-1]+1):dims[j],(dims[j-1]+1):dims[j]]<-sliceCorMats[[j]]
}

dat<-allRes_100trees[which(allRes_100trees$tree==1),]
fit1<-gls(richness ~ MRCA, correlation=corSymm(comboCorrs[lower.tri(comboCorrs)],fixed=TRUE),data=dat)
#fit2<-gls(richness ~ MRCA, correlation=corSymm(comboCorrs,fixed=TRUE),data=dat) # reproduced the 'wrong dimension' error...
sum<-summary(fit1)
sum$tTable ## THIS WORKS.

xx<-corSymm(comboCorrs[lower.tri(comboCorrs)],fixed=TRUE)
qq<-Initialize(xx,dat) # THIS WORKS.

mats<-vector("list",length(sliceCorMats))
for(j in 1:length(sliceCorMats)){
	mats[[j]]<-sliceCorMats[[j]][lower.tri(sliceCorMats[[j]])]
}

xx1<-corSymm(mats,form= ~1 | slice,fixed=TRUE)
qq1<-Initialize(xx1,dat) # THIS breaks.


xx1<-corSymm(sliceCorMats,form= ~1 | slice,fixed=TRUE)
qq1<-Initialize(xx1,dat) # THIS breaks.


m1_phylo<-glmmPQL(fixed = richness ~ MRCA, random = ~ 1 | slice, family = poisson, correlation=corSymm(form = ~ 1 | slice,fixed=TRUE), data = dat) # WORKS.

m1_phylo<-glmmPQL(fixed = richness ~ MRCA, random = ~ 1 | slice, family = poisson, correlation=corSymm(form = ~ comboCorrs[lower.tri(comboCorrs)] | slice,fixed=TRUE), data = dat) # WORKS.


m1_phylo<-glmmPQL(fixed = richness ~ MRCA, random = ~ 1 | slice, family = poisson, correlation=corSymm(comboCorrs[lower.tri(comboCorrs)],fixed=TRUE), data = dat)


dd<-pdSymm(value=comboCorrs[lower.tri(comboCorrs)],form= ~ 1 | slice, data=dat)

fit<-lme(richness ~ MRCA, random = ~ 1 | slice, correlation=pdSymm(value=comboCorrs[lower.tri(comboCorrs)],form= ~ 1 | slice, data=dat), data = dat)



corSymm(comboCorrs[lower.tri(comboCorrs, diag=TRUE)],form = ~ 1 | slice,fixed=TRUE)

m1_phylo<-
lme(richness ~ slice, random = list(slice = pdSymm(comboCorrs,~slice)), data = dat)

lme(Y~Random, data = DATA,random = list(Random = pdSymm(CovM,~Random))) 

pdSymm(comboCorrs)

	Corr1<-corPagel(value=1,phy=slicePhys[[j]])
	Corr2<-Initialize(object=Corr1,data=allRes_100trees[which(allRes_100trees$tree==i & allRes_100trees$slice==j*5),])
	sliceCorMats[[j]]<-corMatrix(Corr2, corr=TRUE)	


vcv(slicePhys[[14]], model="Brownian", corr=TRUE)

lowerT<-comboCorrs[lower.tri(comboCorrs,diag=TRUE)]
obj<-corSymm(comboCorrs[lower.tri(comboCorrs)],fixed=TRUE)


comboCorrs<-data.frame(matrix(0, nrow = nrows_all, ncol = nrows_all),row.names=labs_all)
colnames(comboCorrs)<-labs_all

# Now try phylo with GLMM...
DATA<-cbind.data.frame(allRes_100trees[which(allRes_100trees$tree==1),"richness"],allRes_100trees[which(allRes_100trees$tree==1),"MRCA"],allRes_100trees[which(allRes_100trees$tree==1),"slice"])
rownames(DATA)<-labs_all
colnames(DATA)<-c("richness","MRCA","slice")

xx<-corSymm(value=comboCorrs[lower.tri(comboCorrs,diag=FALSE)],fixed=TRUE)#,form= ~1 | slice)
qq<-Initialize(xx,data=allRes_100trees[which(allRes_100trees$tree==1),])

m1_phylo<-glmmPQL(fixed = richness ~ MRCA, random = ~ 1 | slice, family = poisson, correlation=qq, data = allRes_100trees[which(allRes_100trees$tree==1),])

# ^^ This appears to work now?? que interesante... 

# https://r-forge.r-project.org/scm/viewvc.php/pkg/nlmeU/R/corStruct.R?view=markup&root=nlmeu&pathrev=53
corSymm <-
  ## Constructor for the corSymm class
  function(value = numeric(0), form = ~ 1, fixed = FALSE)
{
  attr(value, "formula") <- form
  attr(value, "fixed") <- fixed
  class(value) <- c("corSymm", "corStruct")
  value
}



m1_phylo<-glmmPQL(fixed = richness ~ MRCA, random = ~ 1 | slice, family = poisson, correlation=corSymm(value=comboCorrs[lower.tri(comboCorrs)],fixed=TRUE,form= ~1 | slice), data = DATA)

m1_phylo<-glmmPQL(fixed = richness ~ MRCA, random = ~ 1 | slice , family = poisson, correlation=corSymm(value=comboCorrs,fixed=TRUE,form= ~1 | slice), data = allRes_100trees[which(allRes_100trees$tree==1),])
sum1<-summary(m1)
sum1



# http://grokbase.com/p/r/r-sig-phylo/117ntq7bqg/follow-up-for-lambda-transform-query-how-to-use-an-existing-covariance-matrix-directly-in-a-gls-procedure-in-r
data(bird.orders)
dat<-data.frame(y=rnorm(23),x=rnorm(23))
rownames(dat)<-bird.orders$tip.label
mat<-vcv(bird.orders,cor=TRUE)
# following Blomberg to here
lambda<-0.75 # for instance
fit1<-gls(y~x,correlation=corSymm(lambda*mat[lower.tri(mat)],fixed=TRUE),data=dat)
# compare to corPagel
fit2<-gls(y~x,correlation=corPagel(0.75,bird.orders,fixed=TRUE),data=dat)
summary(fit1)
summary(fit2)







# Need the CORRELATION structure in there though... 
# For each slice on each tree, there is an underlying structure of how the clades are related.
# USE the fact that the clades are NUMBERED in the same way that the slice phylogeny is created... RIGHT?? ..?


# Load in slice tree REP data
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")

slicePhys<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_slicePhy-5to70Ma_by5.trees")

# load back in tables
results<-vector("list",length(allCladeSets))
for (i in 1: length(allCladeSets)){
	res<-read.table(paste("MamPhy_5911sp_",bbone,"_cladeLevel_DRstats_trees",allCladeSetNames[i],".txt",sep=""), header=TRUE)
	results[[i]]<-res
}

# example: 1060 clades, richness ~ MRCA + DR_harm
i=1
rownames(results[[i]])<-slicePhys[[i]]$tip.label
cladeData<-treedata(slicePhys[[i]],na.omit(results[[i]]))

form<-(log(richness) ~ MRCA_mcc + DR_harm + DR_skew)
fit<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=as.data.frame(cladeData$data), method="ML")
sum<-summary(fit) #ss<-residuals(fit)

fit2<-gls(form, correlation=corPagel(value=1,phy=cladeData$phy), data=as.data.frame(cladeData$data), method="ML")
sum2<-summary(fit2)






# fixed effects
richness ~ age + dr.mean + dr.skew + time.slice + age:time.slice + dr.mean:time.slice + dr.skew:time.slice

# random effects
~ (age|tree|time.slice) + (dr.mean|tree|time.slice) + (dr.skew|tree|time.slice)





######
# OTHER ways to calculate stuff:

# original, for newick tree: ymle = function(tree) (.subset2(tree,2)-1L)/sum(.subset2(tree,4)) # this take the # of number of nodes in a tree / sum of branch lengths.
ymle = function(tree){ 
	(.subset2(tree,3)-1L)/sum(.subset2(tree,2)) # this take the # of number of nodes in a tree (minus 1) / sum of branch lengths.
}
tree=cladeSet[[1]]
ymle(tree)

# check the results:
yule<-fit_bd(phylo=cladeSet[[k]], tot_time=max(node.age(cladeSet[[k]])$ages), f.lamb=rate_conT, f.mu=rate_conT, lamb_par, mu_par=0, f = 1, meth = "Nelder-Mead", cst.lamb = TRUE, cst.mu = TRUE, expo.lamb = FALSE, expo.mu = FALSE, fix.mu = TRUE, dt=1e-3, cond = "crown")
bd$lamb_par+bd$mu_par	

# BD Kendall and Moran (complete= no extinction)
BD.km[k,]<-bd.km(phy=cladeSet[[k]], missing=0, crown=TRUE) # Assuming no extinction

dd<-birthdeath(cladeSet[[6]])
dd$para
bd$lamb_par
bd$mu_par

#morlon model
bd<-fit_bd(phylo=cladeSet[[k]], tot_time=max(node.age(cladeSet[[k]])$ages), f.lamb=rate_conT, f.mu=rate_conT, lamb_par, mu_par, f = 1, meth = "Nelder-Mead", cst.lamb = TRUE, cst.mu = TRUE, expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE, dt=1e-3, cond = "crown")
# set params for BD model
rate_conT<-function(t,y){y[1]} # BOTH constant rate through time
lamb_par<-c(0.01)
mu_par<-c(0.005)
		BD_Lam[k,]<-bd$lamb_par
		BD_Mu[k,]<-bd$mu_par	
		BD_Div1[k,]<-bd$lamb_par-bd$mu_par
		BD_Div2[k,]<-bd$lamb_par+bd$mu_par




########

i=1
rownames(results[[i]])<-slicePhys[[i]]$tip.label
cladeData<-treedata(slicePhys[[i]],na.omit(results[[i]]))

form<-(log(richness) ~ MRCA_mcc + DR_harm + DR_skew)
fit<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=as.data.frame(cladeData$data), method="ML")
sum<-summary(fit) #ss<-residuals(fit)

fit2<-gls(form, correlation=corPagel(value=1,phy=cladeData$phy), data=as.data.frame(cladeData$data), method="ML")
sum2<-summary(fit2)


glmer
vcvMAM<-vcv(corPagel(1, cladeData$phy))

     library(nlme) # will be loaded automatically if omitted
     summary(glmmPQL(y ~ trt + I(week > 2), random = ~ 1 | ID,
                     family = binomial, data = bacteria))
     

# fixed effects
richness ~ age + dr.mean + dr.skew + time.slice + age:time.slice + dr.mean:time.slice + dr.skew:time.slice

# random effects
~ (age|tree|time.slice) + (dr.mean|tree|time.slice) + (dr.skew|tree|time.slice)








# lme = linear mixed effects = multi-level model


# random intercept

m1 = lmer(dr ~ st.ele + (1|poId), data = sdat)	

form = dr ~ st.ele + f(poId, model = 'iid')


#organize data

n.loc = max(sdat$poId)

sdat$i.int = sdat$poId

sdat$j.int = sdat$poId + n.loc

sdat$k.int = sdat$j.int + n.loc


# uncorrelated random intercept and random slope

m2 = lmer(dr ~ st.ele + (1|poId) + (0 + st.ele|poId), data = sdat)

form = dr ~ st.ele + f(i.int, model = 'iid') + f(j.int, st.ele, model = 'iid')


# correlated random intercept and random slope

m3 = lmer(dr ~ st.ele + (st.ele |poId), data = sdat)

form = dr ~ st.ele + f(i.int, model = 'iid2d', n = 2*n.loc) + f(j.int, st.ele, copy = 'i.int')





# uncorrelated random intercept and random slope of st.ele^2

m4 = lmer(dr ~ st.ele + I(st.ele^2) + (1|poId) + (0 + I(st.ele^2)|poId), data = sdat)

form = dr ~ st.ele + I(st.ele^2) + f(i.int, model = 'iid') + f(j.int, I(st.ele^2), model = 'iid')

	

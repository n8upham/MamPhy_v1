#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy -- Upham et al. 2017
###
# Figure 4 - TIP-level to explain clade diversification rate
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%50
# Perform the tip (species) level analyses, to explain the div rates using ecological variables
######
# - Load in 1 tree of 100
# - Calculate (or load) ES and DR on per-tip basis
# - Subset out trait data comparison
# - Now in parallel, do the PGLS also...


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# packages
library(moments); library(nlme); library(ape); library(geiger); library(picante); library(phytools); library(caper)
library(foreach);library(doSNOW)

# directory and source
#dirname = "ysm-gpfs/home/nu35/project/speciesLevel_explainingTipDR"
##dirname = "/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/timeSlices_CladeDiv"
#setwd(dirname)
#source("DR_functions.R")

# open cluster for parallel processing
cl = makeCluster(100, type = 'MPI', outfile="")
registerDoSNOW(cl)

# start parallel loop
ntrees=100
foreach(i=1:ntrees, .packages=c('geiger','moments', 'nlme', 'ape'), .combine=cbind, .verbose=TRUE) %dopar% {

# which backbone?
bbone<- "NDexp" #"FBD" # 

#==================
# Load in 1 tree of 100
#==================
mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_nexus.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 
#write.tree(mamPhy,file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""))
#tree1=scan(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""), what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

#==================
# Calculate ES and DR on per-tip basis
#==================
# gives pairwise clade matrix from CAIC function
#clade_matrix = readCAIC(tree1)
#
## calculate and write to file
#DR = 1/ES_v2(clade_matrix)
#ES = ES_v2(clade_matrix)
#res = cbind.data.frame(DR,ES)
#res1 = res[order(rownames(res)),]
#
#write.table(res1, file=paste(bbone,"_sample100_",i,"_DRtips.txt",sep=""))
res1<-read.table(file=paste(bbone,"_sample100_",i,"_DRtips.txt",sep=""), header=TRUE)

# tip DR values to use for this tree's comparisons
tipDR_i<-res1$DR
names(tipDR_i)<-rownames(res1)

#==================
# Load in stats about tip TAXONOMY and GENE SAMPLING
#==================
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
head(cladesDR)
sorted<-cladesDR[order(cladesDR$tiplabel),]

taxonomyVars<-c("tiplabel","gen", "fam", "ord", "higher", "genes")
taxFinal<-sorted[,taxonomyVars]

#==================
# Subset by comparison of ecological traits to use
#==================
tipDataAll<-read.table(file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)
	
PRED_LOGvars<-c("homeRange_km2_ext", "trophic123", "geoArea_km2", "DispDistAll_km_ext", "lifemode1234", "GenerationLength_d", "activity123") # "BM_final_kg", 
PRED_NONvars<-c("Lat_centroid", "Lon_centroid")
PRED_names<-c(paste("log",c("HomeRange_km2","Trophic123", "GeoArea_km2", "DispersalIndex", "Lifemode1234","GenLength_d", "Activity123"),sep=""), PRED_NONvars, "genes") # "BodyMass_kg", "DispDist_km", 
#PRED_names<-c(paste("log",c("HomeRange_km2","BodyMass_kg", "DispDist_km"),sep=""), PRED_NONvars,"genes") #"Trophic123", 

catsSelected<-c("CarnOrNot", "HerbOrNot", "OmniOrNot", "AquaOrNot", "ArboOrNot", "FlysOrNot", "SubtOrNot", "TerrOrNot", "NoctOrNot", "CathOrNot", "DiurOrNot", "MarineOrNot")

pred_all<-cbind.data.frame(log(tipDataAll[,PRED_LOGvars[1:3]]),log(tipDataAll[,PRED_LOGvars[4]]*tipDataAll[,PRED_LOGvars[5]]), log(tipDataAll[,PRED_LOGvars[5:7]]),tipDataAll[,PRED_NONvars],(taxFinal$genes))
predScale<-scale(pred_all,center=TRUE, scale=TRUE)

allDatFinal<-cbind.data.frame(taxFinal[,1:5],log(tipDR_i), predScale, tipDataAll[,catsSelected])
colnames(allDatFinal)<-c(taxonomyVars[1:5],"logTipDR",PRED_names,catsSelected)

respVar<-"logTipDR"
# if doing GLOBAL
datFinal<-allDatFinal[,c("tiplabel",respVar,PRED_names)] # just the 10 vars in the MODELS

GROUP<-"global"
z<-1

# if doing SUBSETS
#for(z in 1:length(catsSelected)){
#datFinal<-allDatFinal[which(allDatFinal[,catsSelected[z]]=="1"),]

predictors<-PRED_names
perPredictor<-vector("list",length(predictors))
	for(j in 1:length(predictors)){

	pred<-predictors[j]

	dat1<-cbind.data.frame(datFinal[,c("logTipDR",pred)])
	colnames(dat1)<-c("logTipDR",pred)
	rownames(dat1)<-datFinal$tiplabel

	dat<-na.omit(dat1)
	treeDat<-treedata(mamPhy,dat)
	datPlot<-as.data.frame(treeDat$dat)

	# all predictors as NON-LOG
	form<-as.formula(paste("logTipDR ~ ", pred,sep=""))

#	# some predictors as LOG
#	if(j < 4) {	form<-as.formula(paste("log(tipDR) ~ ", predictors[j],sep=""))
#				} else { 
#				form<-as.formula(paste("log(tipDR) ~ log(", predictors[j], ")",sep=""))}

	#	fit1<-gls(form, data=dat, method="ML")
		for (p in c(0.5,seq(1,0,by=-0.01))) {possibleError <- tryCatch(
		      gls(form, correlation=corPagel(value=p,phy=treeDat$phy), data=datPlot, method="ML"),
		      error=function(e) e)
		if(inherits(possibleError, "gls")) break		
		if(inherits(possibleError, "error")) next}
		fit1<-possibleError

		sum<-summary(fit1)

		a=round(sum$tTable[1], digits=6)
		b=round(sum$tTable[2], digits=6)
		SE=round(sum$tTable[4], digits=6)
		pVal=round(sum$tTable[8], digits=6)
		lambda=round(sum$modelStruct[[1]][[1]], digits=3)
	#	lambda=NA

		perPredictor[[j]]<-c(a,b,SE,pVal,lambda)
	}

allRes_i<-cbind.data.frame(do.call(rbind,perPredictor),rep(i,length(predictors)),rep(GROUP[z],length(predictors)))
colnames(allRes_i)<-c("a","b","SE","pVal","lam","tree","cat")
#toNameRows<-c(predictors[1:3],paste("log(",predictors[4:8],")",sep=""))
#rownames(allRes_i)<-toNameRows
rownames(allRes_i)<-predictors

#write.table(allRes_i,file=paste("speciesLevel_GLS-NoTree_LogTipDR-vs-Log-EcoGeoTraits_",GROUP,"_tree",i,".txt",sep=""))
#write.table(allRes_i,file=paste("speciesLevel_PGLS-PAGEL_LogTipDR-vs-Log-EcoGeoTraits_",GROUP,"_tree",i,".txt",sep=""))
write.table(allRes_i,file=paste("tipLevel_PGLS-PAGEL_LogTipDR-vs-Log-EcoGeoTraits_",GROUP,"_tree",i,"_NEW_10vars.txt",sep=""))

#} # cycle the subset groups

} # cycle 100 trees

stopCluster(cl)

q()

n




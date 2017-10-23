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
library(moments); library(nlme); library(ape); library(geiger) # library(picante); library(phytools)
library(foreach);library(doSNOW)

# directory and source
#dirname = "ysm-gpfs/home/nu35/project/speciesLevel_explainingTipDR"
##dirname = "/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/timeSlices_CladeDiv"
#setwd(dirname)
source("DR_functions.R")

# open cluster for parallel processing
cl = makeCluster(25, type = 'SOCK', outfile="")
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
tree1=scan(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""), what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

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

taxonomyVars<-c("gen", "fam", "ord", "higher", "genes")
taxFinal<-sorted[,taxonomyVars]

#==================
# Subset by comparison of ecological traits to use
#==================
tipDataAll<-read.table(file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)
	
varsSelected<-c("Lat_centroid", "Lon_centroid","geoArea_km2", "BM_final_g", "homeRange_km2_ext", "DispDistAll_km_ext", "GenerationLength_d")
catsSelected<-c("CarnOrNot", "HerbOrNot", "OmniOrNot", "AquaOrNot", "ArboOrNot", "FlysOrNot", "SubtOrNot", "TerrOrNot", "NoctOrNot", "CathOrNot", "DiurOrNot", "MarineOrNot")

allDatFinal<-cbind.data.frame(tipDR_i,taxFinal,tipDataAll[,varsSelected],tipDataAll[,catsSelected])
## if doing GLOBAL
#catsSelected<-"global"; z=1
#datFinal<-allDatFinal 

# if doing TERRESTRIAL-only
catsSelected<-"terrOnly"; z=1
datFinal<-allDatFinal[which(allDatFinal[,"MarineOrNot"]=="0"),] 

## if doing SUBSETS
#for(z in 1:length(catsSelected)){
#datFinal<-allDatFinal[which(allDatFinal[,catsSelected[z]]=="1"),]

varsSelected2<-c("genes",varsSelected)
perPredictor<-vector("list",length(varsSelected2))
for(j in 1:length(varsSelected2)){

var<-varsSelected2[j]

dat1<-cbind.data.frame(datFinal[,c("tipDR_i",var)])
colnames(dat1)<-c("tipDR",var)
rownames(dat1)<-rownames(datFinal)

dat<-na.omit(dat1)
treeDat<-treedata(mamPhy,dat)
datPlot<-as.data.frame(treeDat$dat)

if(j < 4) {	form<-as.formula(paste("log(tipDR) ~ ", varsSelected2[j],sep=""))
			} else { 
			form<-as.formula(paste("log(tipDR) ~ log(", varsSelected2[j], ")",sep=""))}

#	fit1<-gls(form, data=dat, method="ML")
	for (p in c(0.5,seq(1,0,by=-0.01))) {possibleError <- tryCatch(
	      gls(form, correlation=corPagel(value=p,phy=treeDat$phy), data=datPlot, method="ML"),
	      error=function(e) e)
	if(inherits(possibleError, "gls")) break		
	if(inherits(possibleError, "error")) next}
	fit1<-possibleError

	sum<-summary(fit1)

	a=round(sum$tTable[1], digits=3)
	b=round(sum$tTable[2], digits=3)
	SE=round(sum$tTable[4], digits=3)
	pVal=round(sum$tTable[8], digits=3)
	lambda=round(sum$modelStruct[[1]][[1]], digits=3)
#	lambda=NA

	perPredictor[[j]]<-c(a,b,SE,pVal,lambda)
}

allRes_i<-cbind.data.frame(do.call(rbind,perPredictor),rep(i,length(varsSelected2)),rep(catsSelected[z],length(varsSelected2)))
colnames(allRes_i)<-c("a","b","SE","pVal","lam","tree","cat")
toNameRows<-c(varsSelected2[1:3],paste("log(",varsSelected2[4:8],")",sep=""))
rownames(allRes_i)<-toNameRows

GROUP<-catsSelected[z]

#write.table(allRes_i,file=paste("speciesLevel_GLS-NoTree_LogTipDR-vs-Log-EcoGeoTraits_",GROUP,"_tree",i,".txt",sep=""))
write.table(allRes_i,file=paste("speciesLevel_PGLS-PAGEL_LogTipDR-vs-Log-EcoGeoTraits_",GROUP,"_tree",i,".txt",sep=""))

#} # cycle the subset groups

} # cycle 100 trees

q()

n




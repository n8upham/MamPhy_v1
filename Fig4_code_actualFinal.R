#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy -- Upham et al. 2017
###
# Figure 4 - time slices to explain clade diversification rate
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Perform the timeslice analyses, now with summary stats for the div rates plus ecological variables
######
# - Load in 1 tree of 100
# - Calculate ES and DR on per-tip basis
# - Make time slices every 5 Ma to create clades
# - Calculate per-slice, per-clade summary values 
# - Now in parallel, and with doing the PGLS also...



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# packages
library(moments); library(nlme); library(ape); library(picante); library(phytools); library(geiger)
library(foreach);library(doSNOW)

# directory and source
dirname = "/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/timeSlices_CladeDiv"
setwd(dirname)
source("DR_functions.R")

# open cluster for parallel processing
cl = makeCluster(15, type = 'SOCK', outfile="")
registerDoSNOW(cl)

# start parallel loop
ntrees=100
foreach(i=1:ntrees, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

# which backbone?
bbone<- "NDexp" #"FBD" # 

# load in stats about TIP SAMPLING
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
head(cladesDR)

# load in stats about TIP TRAITS
###
tipDataAll<-read.table(file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)

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
clade_matrix = readCAIC(tree1)

# calculate and write to file
DR = 1/ES_v2(clade_matrix)
ES = ES_v2(clade_matrix)
res = cbind.data.frame(DR,ES)
res1 = res[order(rownames(res)),]

write.table(res1, file=paste(bbone,"_sample100_",i,"_DRtips.txt",sep=""))

#=======================================
# Make time slices every 5 Ma to create clades
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

# write clades to per-slice files
for (j in 1:length(allCladeSets)){
	trees<-allCladeSets[[j]]
	write.tree(trees,file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"_newick.trees",sep=""))
	}

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
	BD_Turn<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD.ms0<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD.ms0p5<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD.ms0p9<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	
	logBodyMass_kg<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))	
	logGeoRange_km2<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	logHomeRange_km2<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	rankHomeRange_1to10<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	rankHomeRange_1to10_binsEq<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	logDispDist_km<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	rankDispDist_1to10<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	rankDispDist_1to10_binsEq<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	logGenLength_d<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))

	trophic123<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	CarnOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	HerbOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	OmniOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	perVert<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	perInvert<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	perPlant<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	AquaOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	ArboOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	FlysOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	SubtOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	TerrOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	NoctOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	CathOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	DiurOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))

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
	# TRAIT DATA
	# Continuous
		logBodyMass_kg[k,]<-mean(na.omit(log(tipDataAll[match(cladeSp,cladesDR$tiplabel),"BM_final_kg"])))
		logGeoRange_km2[k,]<-mean(na.omit(log(tipDataAll[match(cladeSp,cladesDR$tiplabel),"geoRange_km2"])))
		logHomeRange_km2[k,]<-mean(na.omit(log(tipDataAll[match(cladeSp,cladesDR$tiplabel),"homeRange_km2_ext"])))
		rankHomeRange_1to10[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"lnHomeRange_rank1to10"]))
		rankHomeRange_1to10_binsEq[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"lnHomeRange_rank1to10_binsEqual"]))
		logDispDist_km[k,]<-mean(na.omit(log(tipDataAll[match(cladeSp,cladesDR$tiplabel),"DispDistAll_km_ext"])))
		rankDispDist_1to10[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"lnDispDistAll_rank1to10"]))
		rankDispDist_1to10_binsEq[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"lnDispDistAll_rank1to10_binsEqual"]))
		logGenLength_d[k,]<-mean(na.omit(log(tipDataAll[match(cladeSp,cladesDR$tiplabel),"GenerationLength_d"])))
	# Categorical
		trophic123[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"trophic123"]))
		CarnOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"CarnOrNot"]))
		HerbOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"HerbOrNot"]))
		OmniOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"OmniOrNot"]))
		perVert[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"perVert"]))
		perInvert[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"perInvert"]))
		perPlant[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"perPlant"]))
		AquaOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"AquaOrNot"]))
		ArboOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"ArboOrNot"]))
		FlysOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"FlysOrNot"]))
		SubtOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"SubtOrNot"]))
		TerrOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"TerrOrNot"]))
		NoctOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"NoctOrNot"]))
		CathOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"CathOrNot"]))
		DiurOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"DiurOrNot"]))

	if (length(cladeSp) > 2) {
	# Yule model
		PB_Div[k,]<-ymle(cladeSet[[k]])
		# BD model
		bd<-birthdeath(cladeSet[[k]])
		BD_Lam[k,]<-bd$para[[2]]/(1-bd$para[[1]])
		BD_Mu[k,]<-bd$para[[1]]*(bd$para[[2]]/(1-bd$para[[1]]))
		BD_Div[k,]<-bd$para[[2]]
		BD_Turn[k,]<-bd$para[[1]]
		# BD Mag and Sand
		cladeSet[[k]]$root.edge<-0
	    BD.ms0[k,]<-bd.ms(phy=cladeSet[[k]], missing=0, epsilon=0, crown=TRUE) # Assuming no extinction
    	BD.ms0p5[k,]<-bd.ms(phy=cladeSet[[k]], missing=0, epsilon=0.5, crown=TRUE) # Assuming medium extinction 
     	BD.ms0p9[k,]<-bd.ms(phy=cladeSet[[k]], missing=0, epsilon=0.9, crown=TRUE) # Assuming high extinction
		} else NULL
	}

	res2<-cbind(DR_harm, DR_cv, DR_skew, DR_kurt, percentSamp, richness, MRCA, PB_Div, BD_Lam, BD_Mu, BD_Div, BD_Turn, BD.ms0, BD.ms0p5, BD.ms0p9, logBodyMass_kg,logGeoRange_km2,logHomeRange_km2,rankHomeRange_1to10,rankHomeRange_1to10_binsEq,logDispDist_km,rankDispDist_1to10,rankDispDist_1to10_binsEq,logGenLength_d,trophic123,CarnOrNot,HerbOrNot,OmniOrNot,perVert,perInvert,perPlant,AquaOrNot,ArboOrNot,FlysOrNot,SubtOrNot,TerrOrNot,NoctOrNot,CathOrNot,DiurOrNot,i, j*5)

	colnames(res2)<-c("DR_harm","DR_cv", "DR_skew", "DR_kurt", "percentSamp", "richness", "MRCA", "PB_Div", "BD_Lam", "BD_Mu", "BD_Div", "BD_Turn", "BD.ms0", "BD.ms0p5", "BD.ms0p9", "logBodyMass_kg", "logGeoRange_km2", "logHomeRange_km2", "rankHomeRange_1to10", "rankHomeRange_1to10_binsEq", "logDispDist_km", "rankDispDist_1to10", "rankDispDist_1to10_binsEq", "logGenLength_d", "trophic123", "CarnOrNot", "HerbOrNot", "OmniOrNot", "perVert","perInvert","perPlant","AquaOrNot", "ArboOrNot", "FlysOrNot", "SubtOrNot", "TerrOrNot", "NoctOrNot", "CathOrNot", "DiurOrNot", "tree", "slice")

	rownames(res2)<-paste(i,"_",j,"_",c(1:length(res2[,1])),sep="")
	
	write.table(res2,paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS_withTraits.txt",sep=""))
} # end 1-tree loop across all 14 slices (every 5 Ma)


# create slice phys
slicePhys<-vector("list",length(allCladeSets))
for (j in 1:length(allCladeSets)){
	cladeSet<-allCladeSets[[j]]

	cladeReps<-vector()
	for (k in 1:length(cladeSet)){
		cladeSp<-cladeSet[[k]]$tip.label
		cladeReps[k]<-cladeSp[1]
		}
	toDrop<-setdiff(mamPhy$tip.label,cladeReps)
	slicePhys[[j]]<-drop.tip(mamPhy,toDrop)
	slicePhys[[j]]$tip.label<-paste(i,"_",j,"_",c(1:length(slicePhys[[j]]$tip.label)),sep="")
}

# write slice phys
for(j in 1:length(slicePhys)){
	write.tree(slicePhys[[j]], file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_slicePhy-5to70Ma_by5.trees",sep=""), append=TRUE)
}
} # end 100-tree loop


#####
# How many slices did we end up with?

lengths<-vector("list", length=100)
for (i in 1:length(lengths)){
	lengths[[i]]<-read.table(file=paste(bbone,"_sample100_",i,"_timeslice_Lengths.txt",sep=""))
}
lengthsALL1<-do.call(cbind,lengths)
colnames(lengthsALL1)<-paste("tree_",1:100,sep="")
lengthsALL<-as.matrix(lengthsALL1)

means<-c()
for(i in 1:length(lengthsALL[,1])){
	means[i]<-mean(lengthsALL[i,])
}
names(means)<-allCladeSetNames 
#   10Ma   15Ma   20Ma   25Ma   30Ma   35Ma   40Ma   45Ma   50Ma   55Ma   60Ma 
# 655.98 379.38 238.00 158.98 115.71  89.60  71.33  57.84  44.30  35.92  28.88 


















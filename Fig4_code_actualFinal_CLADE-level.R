#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy -- Upham et al. 2017
###
# Figure 4 - time slices to explain clade diversification rate
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create timeslice clades, calculate summary stats for div rates and ecological variables
######
# - In parallel, load in 1 tree of 100
# - Calculate ES and DR on per-tip basis
# - Make time slices every 5 Ma to create clades
# - Relabel the backbone "slicePhys" uniting  time slice clades
# - Calculate per-slice, per-clade summary values-- including singletons
# - end 100-tree loop

# Run PGLS -- aim to ** explain variation DIV RATE ** using timeslice clade summary stats 
######
# - In parallel, read back in slice phylo backbones (14 slices at 5 Ma intervals)
# - Load in rate and trait slice summaries for 1 of 100 trees (standardize predictors)
# - Run loop of 4 response vars 
# - Setup empty dataframes to receive PGLS results 
# - Nest the loop of 26 predictor vars


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# packages
library(moments); library(nlme); library(ape); library(picante); library(phytools); library(geiger)
library(foreach);library(doSNOW)

# directory and source
#dirname = "/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/timeSlices_CladeDiv"
#setwd(dirname)
source("DR_functions.R")

# open cluster for parallel processing
cl = makeCluster(100, type = 'MPI', outfile="")
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
#tree1=scan(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""), what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 


##==================
## Calculate ES and DR on per-tip basis
##==================
## gives pairwise clade matrix from CAIC function
#clade_matrix = readCAIC(tree1)
#
## calculate and write to file
#DR = 1/ES_v2(clade_matrix)
#ES = ES_v2(clade_matrix)
#res = cbind.data.frame(DR,ES)
#res1 = res[order(rownames(res)),]
#
#write.table(res1, file=paste(bbone,"_sample100_",i,"_DRtips.txt",sep=""))
res1<-read.table(file=paste(bbone,"_sample100_",i,"_DRtips.txt",sep=""))

#=======================================
# Make time slices every 5 Ma to create clades
#=======================================
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery

root=max(node.age(mamPhy)$ages)

allCladeSets<-vector("list",length=numSlices)
slicePhys<-vector("list",length=numSlices)
allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
	# make tipward slice clades
	allCladeSets[[j]]<-treeSlice(mamPhy, slice=root-(sliceEvery*j), trivial=TRUE) # now keeping trivial (single tips)
		# write clades to per-slice files
		write.tree(allCladeSets[[j]],file=paste(bbone,"_sample100_",i,"_timeslice_cladeTreesWithSingletons_",(sliceEvery*j),"Ma.trees",sep=""))
	# make rootward slice *backbones*
	slicePhys[[j]]<-treeSlice(mamPhy, slice=root-(sliceEvery*j), trivial=TRUE, orientation="rootwards") # toward root, slicePhys
}

# re-label the slicePhys with standard names
# Note: need to match tips with nodes to re-label these properly! 
for (j in 1:numSlices){
	cladeSet<-allCladeSets[[j]]
	newTipNames<-paste(i,"_",j,"_",c(1:length(cladeSet)),sep="")
	
	for (k in 1:length(cladeSet)){
		cladeSp<-cladeSet[[k]]$tip.label
		if(length(cladeSp)==1){ 
			slicePhys[[j]]$tip.label[which(slicePhys[[j]]$tip.label==cladeSp)]<-newTipNames[k]
		} else {
			node <- getMRCA(mamPhy, cladeSp) # find the MRCA node of those species
			slicePhys[[j]]$tip.label[which(slicePhys[[j]]$tip.label==node)]<-newTipNames[k]
		}
	}
}

# write slice phys (multiPhy object)
for(j in 1:length(slicePhys)){
	write.tree(slicePhys[[j]], file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_CORRECT-sliceRootwardSingletons_slicePhy-5to70Ma_by5.trees",sep=""), append=TRUE)
	}

# record clade sizes per slice
allLengths_i<-vector("list",length(allCladeSets))
for (j in 1:length(allCladeSets)){
	lengths<-vector()
	for(k in 1:length(allCladeSets[[j]])){
		lengths[k]<-length(allCladeSets[[j]][[k]]$tip.label)
	}
allLengths_i[[j]]<-lengths
}
names(allLengths_i)<-allCladeSetNames 
save(allLengths_i, file=paste(bbone,"_sample100_",i,"_timeslice_cladeRichnesses_wSingletons.Rda",sep=""))

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
	#DR_cv<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	DR_skew<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	#DR_kurt<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	percentSamp<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	richness<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	MRCA<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	PB_Div<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD_Lam<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD_Mu<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD_Div<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	BD_Turn<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	#BD.ms0<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	#BD.ms0p5<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	#BD.ms0p9<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	
	BodyMass_kg<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))	
	GeoArea_km2<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	Lat_centroid<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	Lon_centroid<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))

	HomeRange_km2<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	rankHomeRange_1to10<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	rankHomeRange_1to10_binsEq<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	DispDist_km<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	rankDispDist_1to10<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	rankDispDist_1to10_binsEq<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	GenLength_d<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	
	perPlant<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	perInvert<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	perVert<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	trophic123<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	HerbOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	OmniOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	CarnOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	
	lifemode1234<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	AquaOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	ArboOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	FlysOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	SubtOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	TerrOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))

	activity123<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	NoctOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	CathOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	DiurOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))

	for (k in 1:length(cladeSet)){
	cladeSp<-cladeSet[[k]]$tip.label
		x<-res1[match(cladeSp,rownames(res1)),"DR"]
		DR_harm[k,] <- 1/(mean(1/x))
		percentSamp[k,] <- length(which(cladesDR[match(cladeSp,cladesDR$tiplabel),"samp"]=="sampled"))/length(cladeSp)
		richness[k,] <- length(cladeSp)
	if(length(cladeSp) > 1){
		#DR_cv[k,] <- (sd(x)/mean(x))*100
		DR_skew[k,] <- skewness(x)
		#DR_kurt[k,] <- kurtosis(x)
		node <- getMRCA(mamPhy, cladeSp)
		MRCA[k,] <- btimes[node-5911] #taking the height of SAMPLED tree
		}
	# TRAIT DATA
	# Continuous
		BodyMass_kg[k,]<-mean(na.omit((tipDataAll[match(cladeSp,tipDataAll$tiplabel),"BM_final_kg"])))
		GeoArea_km2[k,]<-mean(na.omit((tipDataAll[match(cladeSp,tipDataAll$tiplabel),"geoArea_km2"])))
		Lat_centroid[k,]<-mean(na.omit((tipDataAll[match(cladeSp,tipDataAll$tiplabel),"Lat_centroid"])))
		Lon_centroid[k,]<-mean(na.omit((tipDataAll[match(cladeSp,tipDataAll$tiplabel),"Lon_centroid"])))

		HomeRange_km2[k,]<-mean(na.omit((tipDataAll[match(cladeSp,tipDataAll$tiplabel),"homeRange_km2_ext"])))
		rankHomeRange_1to10[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"lnHomeRange_rank1to10"]))
		rankHomeRange_1to10_binsEq[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"lnHomeRange_rank1to10_binsEqual"]))
		DispDist_km[k,]<-mean(na.omit((tipDataAll[match(cladeSp,tipDataAll$tiplabel),"DispDistAll_km_ext"])))
		rankDispDist_1to10[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"lnDispDistAll_rank1to10"]))
		rankDispDist_1to10_binsEq[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"lnDispDistAll_rank1to10_binsEqual"]))
		GenLength_d[k,]<-mean(na.omit((tipDataAll[match(cladeSp,tipDataAll$tiplabel),"GenerationLength_d"])))
	# Categorical
		perPlant[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"perPlant"]))
		perInvert[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"perInvert"]))
		perVert[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"perVert"]))
		trophic123[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"trophic123"]))
		HerbOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"HerbOrNot"]))
		OmniOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"OmniOrNot"]))
		CarnOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"CarnOrNot"]))
		
		lifemode1234[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"lifemode1234"]))
		AquaOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"AquaOrNot"]))
		ArboOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"ArboOrNot"]))
		FlysOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"FlysOrNot"]))
		SubtOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"SubtOrNot"]))
		TerrOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"TerrOrNot"]))

		activity123[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"activity123"]))
		NoctOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"NoctOrNot"]))
		CathOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"CathOrNot"]))
		DiurOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,tipDataAll$tiplabel),"DiurOrNot"]))

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
		#cladeSet[[k]]$root.edge<-0
	    #BD.ms0[k,]<-bd.ms(phy=cladeSet[[k]], missing=0, epsilon=0, crown=TRUE) # Assuming no extinction
    	#BD.ms0p5[k,]<-bd.ms(phy=cladeSet[[k]], missing=0, epsilon=0.5, crown=TRUE) # Assuming medium extinction 
     	#BD.ms0p9[k,]<-bd.ms(phy=cladeSet[[k]], missing=0, epsilon=0.9, crown=TRUE) # Assuming high extinction
		}
	}

	res2<-cbind.data.frame(DR_harm, DR_skew, percentSamp, richness, MRCA, PB_Div, BD_Lam, BD_Mu, BD_Div, BD_Turn, BodyMass_kg,GeoArea_km2, Lat_centroid, Lon_centroid, HomeRange_km2,rankHomeRange_1to10,rankHomeRange_1to10_binsEq,DispDist_km,rankDispDist_1to10,rankDispDist_1to10_binsEq,GenLength_d,perPlant,perInvert,perVert,trophic123,HerbOrNot,OmniOrNot,CarnOrNot,lifemode1234,AquaOrNot,ArboOrNot,FlysOrNot,SubtOrNot,TerrOrNot,activity123,NoctOrNot,CathOrNot,DiurOrNot,i, j*-5)

	colnames(res2)<-c("DR_harm","DR_skew", "percentSamp", "richness", "MRCA", "PB_Div", "BD_Lam", "BD_Mu", "BD_Div", "BD_Turn", "BodyMass_kg", "GeoArea_km2", "Lat_centroid", "Lon_centroid", "HomeRange_km2", "rankHomeRange_1to10", "rankHomeRange_1to10_binsEq", "DispDist_km", "rankDispDist_1to10", "rankDispDist_1to10_binsEq", "GenLength_d", "perPlant","perInvert","perVert","trophic123", "HerbOrNot","OmniOrNot","CarnOrNot","lifemode1234","AquaOrNot","ArboOrNot","FlysOrNot","SubtOrNot","TerrOrNot","activity123","NoctOrNot","CathOrNot","DiurOrNot","tree", "time")

	rownames(res2)<-paste(i,"_",j,"_",c(1:length(cladeSet)),sep="") # same as the slicePhy names
	
	write.table(res2,paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"_cladeSTATS_withEcoTraits_andSingletons.txt",sep=""))
} # end 1-tree loop across all 14 slices (every 5 Ma)


} # end 100-tree loop


# ==============
# do the PGLS-- UNIVARIATE
#######
# intialize
library(ape); library(phytools); library(picante); library(geiger); library(moments); library(nlme)
library(foreach);library(doSNOW)

# which backbone?
bbone<- "NDexp" #"FBD" # 

# get the clade names
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery
allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

# start cluster
cl = makeCluster(100, type = 'MPI')
registerDoSNOW(cl)

#ranTrees<-c(2, 6, 10, 11, 12, 21, 23, 24, 25, 27, 29, 31, 33, 40, 41, 48, 53, 54, 59, 60, 63, 65, 67, 68, 72, 73, 74, 77, 82, 83, 84, 85, 86, 91, 94, 96, 98)
#toRun<-setdiff(1:100,ranTrees)
#foreach(i=toRun, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

ntrees=100
foreach(i=1:ntrees, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {


# read back in slice phys
# ========================
slicePhys<-read.tree(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_CORRECT-sliceRootwardSingletons_slicePhy-5to70Ma_by5.trees",sep=""))

## multiply edges by 100 for corPagel
#modTrees<-vector("list",length(slicePhys))
#for (j in 1:length(slicePhys)){
#	modTree<-slicePhys[[j]]
#	modTree$edge.length<-modTree$edge.length*100
#	modTrees[[j]]<-modTree
#}

# Load in slice clade summaries and STANDARDIZE
# ========================================
results<-vector("list",length(allCladeSetNames))
for (q in 1: length(allCladeSetNames)){
	# Load in
	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[q]],"_cladeSTATS_withEcoTraits_andSingletons.txt",sep=""), header=TRUE)
	# select RESPONSE vars
	RESP<-cbind.data.frame(res[,c("DR_harm","BD_Div","MRCA","DR_skew")])
	# select PREDICTOR vars
	PRED_LOGvars<-c("HomeRange_km2", "GeoArea_km2", "DispDist_km", "lifemode1234", "GenLength_d","trophic123", "activity123") # "BM_final_kg", "", 
	PRED_NONvars<-c("Lat_centroid", "Lon_centroid", "percentSamp")
	PRED_names<-c(paste("log",c("HomeRange_km2","GeoArea_km2", "DispersalIndex", "Lifemode1234","GenLength_d", "Trophic123", "Activity123"),sep=""), PRED_NONvars)  

	PRED_all<-cbind.data.frame(log(res[,PRED_LOGvars[1:2]]),log(res[,PRED_LOGvars[3]]*res[,PRED_LOGvars[4]]),log(res[,PRED_LOGvars[4:7]]),res[,PRED_NONvars])
		# STANDARDIZE predictors as Z-SCORES
		predScale<-scale(PRED_all, center=TRUE,scale=TRUE)

	# Join together RESP and PRED
	results[[q]]<-cbind.data.frame(log(RESP[,1:3]),RESP[,4],predScale)
	#results[[q]]<-cbind.data.frame(log(RESP[,1:3]),RESP[,4],PRED_all)
	colnames(results[[q]])<-c("logTipDR","logBD_Div","logMRCA","DRskew",PRED_names)
}


# Run for ONE of 14 time slices 
# ===============================

for (j in 1:length(results)){
#j=2 # 10 Ma time slice
#j=7 # 35 Ma time slice
#j=9 # 45 Ma time slice
#j=12 # 60 Ma time slice
cladeData<-treedata(slicePhys[[j]],na.omit(results[[j]])) 
dat<-as.data.frame(cladeData$data)


# Setup empty dataframes to receive PGLS results 
# ==============================================
# UNIVARIATE first
sliceTimes<-seq(-5,-70,-5)

predictors<-PRED_names
nCols<-length(predictors)

uniPGLS_allSlopes<-data.frame(matrix(NA, nrow = length(results), ncol = nCols),row.names=sliceTimes)
colnames(uniPGLS_allSlopes)<-predictors
uniPGLS_allSEs<-data.frame(matrix(NA, nrow = length(results), ncol = nCols),row.names=sliceTimes)
colnames(uniPGLS_allSEs)<-c(paste("SE_",1:nCols,sep=""))
uniPGLS_allInts<-data.frame(matrix(NA, nrow = length(results), ncol = nCols),row.names=sliceTimes)
colnames(uniPGLS_allInts)<-c(paste("i_",1:nCols,sep=""))
uniPGLS_allPs<-data.frame(matrix(NA, nrow = length(results), ncol = nCols),row.names=sliceTimes)
colnames(uniPGLS_allPs)<-c(paste("p_",1:nCols,sep=""))
uniPGLS_allLams<-data.frame(matrix(NA, nrow = length(results), ncol = nCols),row.names=sliceTimes)
colnames(uniPGLS_allLams)<-c(paste("lam_",1:nCols,sep=""))


# Run for ONE of 4 response vars 
# ===============================
respVar<-colnames(results[[1]][,c(1:4)]) 

uniPGLS<-data.frame(matrix(NA, nrow = length(respVar), ncol = nCols*5),row.names=respVar)
colnames(uniPGLS)<-c(colnames(uniPGLS_allInts),colnames(uniPGLS_allSlopes),colnames(uniPGLS_allPs),colnames(uniPGLS_allSEs),colnames(uniPGLS_allLams))

for(z in 1:length(respVar)){

# Nest the loop of *nCols* predictor vars
# ===============================

# UNIVARIATE
	for(k in 1:length(predictors)){

		form<-as.formula(paste(respVar[z], " ~ ", predictors[k], sep=""))
#		fit1<-gls(form, data=dat, method="ML")
#		fit1<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
		for (p in c(0.5,seq(0,1,by=0.01))) {possibleError <- tryCatch(
		      gls(form, correlation=corPagel(value=p,phy=cladeData$phy), data=dat, method="ML"),
		      error=function(e) e)
		if(inherits(possibleError, "gls")) break		
		if(inherits(possibleError, "error")) next}
		fit1<-possibleError

		sum1<-summary(fit1)
		
		uniPGLS_allInts[j,k]<-round(sum1$tTable[1], digits=9)
		uniPGLS_allSlopes[j,k]<-round(sum1$tTable[2], digits=9)
		uniPGLS_allSEs[j,k]<-round(sum1$tTable[4], digits=9)
		uniPGLS_allPs[j,k]<-round(sum1$tTable[8], digits=9)
		uniPGLS_allLams[j,k]<-round(sum1$modelStruct[[1]][[1]], digits=3)		

		} # cycles each of *nCols* predictors

		uniPGLS[z,]<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs,uniPGLS_allSEs,uniPGLS_allLams)[j,]

	} # cycles the 4 respVars

	#corr<-"NO_TREE"
	#corr<-"BROWNIAN"
	corr<-"PAGEL"
	sliceN<-allCladeSetNames[[j]] #"ALL" #"45Ma-only" #"10to60Ma" # # #"35Ma-only" ##"all" #"60Ma"#"all" #

	write.table(uniPGLS,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-DR-MRCA_perSlice_",sliceN,"_SCALED_uniVar_Correct-rootwardWithSingle_10preds.txt",sep=""))

	} # cycles each time slice

#uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs,uniPGLS_allSEs,uniPGLS_allLams)
#write.table(uniPGLS,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-",respVar[z],"_timeSlices_",sliceN,"_SCALED_uniVar_Correct_7preds.txt",sep=""))

} # cycles the 100 trees




# ==============
# do the SUMMARIZATION of results
#######
# intialize
library(ape); library(phytools); library(picante); library(geiger); library(moments); library(nlme)
library(plotrix)

# which backbone?
bbone<- "NDexp" #"FBD" # 

# get the clade names
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery
allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}
sliceTimes<-seq(-5,-70,-5)

# set dir
dirname="/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/cladeLevel_SLICES_explainingDR_uniVar"
setwd(dirname)

# load back in the results
# =========================
ntrees=100
	# non-log results (96 trees)
	#whichTrees<-c(100,10,11,12,13,14,15,16,17,18,19,1,20,21,22,23,24,25,26,27,28,29,2,30,31,32,33,34,35,36,38,39,3,40,41,42,43,44,45,46,47,48,49,4,50,51,52,53,54,55,56,57,58,59,5,60,61,62,63,64,65,66,67,69,6,70,71,72,73,75,76,77,78,79,7,80,81,82,83,84,85,86,87,88,8,90,91,92,93,94,95,96,97,98,99,9)
	
	# log results (94 trees)
	whichTrees<-c(100,10,11,12,13,14,15,16,17,18,19,1,20,21,22,23,24,25,26,27,29,2,30,31,32,33,34,35,36,38,39,3,40,41,42,43,44,45,46,47,48,49,4,50,52,53,54,55,56,57,58,59,5,60,61,62,63,64,65,66,67,69,6,70,71,72,73,75,76,77,78,79,7,80,81,82,83,84,85,86,87,88,8,90,91,92,93,94,95,96,97,98,99,9)

fewerTrees<-ntrees-length(whichTrees)
nslices=10 # just analyze 5 Ma to 50 Ma slices
corr="PAGEL"

respVar<-c("logTipDR", "logBD_Div", "logMRCA", "DRskew")

RES_PerResp<-vector("list",length(respVar))
for(z in 1:length(respVar)){

	RES_PerSlices<-vector("list",length=nslices)
	for(j in 1:nslices){
	sliceN<-allCladeSetNames[[j]]

		RES_PerTrees<-vector("list",length=fewerTrees)	

		for(i in whichTrees){
		uniPGLS<-read.table(file=paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-DR-MRCA_perSlice_",sliceN,"_SCALED_uniVar_Correct-rootwardWithSingle_10preds.txt",sep=""), header=TRUE)
		#uniPGLS<-read.table(file=paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-4varDR_perSlice_",sliceN,"_SCALED_uniVar_Correct-rootwardWithSingle_16preds.txt",sep=""), header=TRUE)
		#uniPGLS<-read.table(file=paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-4varDR_perSlice_",sliceN,"_SCALED_uniVar_Correct-rootwardWithSingle_16predsNON-LOG.txt",sep=""), header=TRUE)
		tree<-i
		slice<-sliceTimes[[j]]
		num<-j

		RES_PerTrees[[i]]<-cbind(uniPGLS[z,],slice,num,tree)
		}
		
	RES_PerSlices[[j]]<-do.call(rbind,RES_PerTrees)
	}
RES_PerResp[[z]]<-do.call(rbind,RES_PerSlices)
} # end 4 resp loop


# calculate 95% CIs (and significance) for the slopes per tree & slice
# ====================================================================
#predictors<-colnames(exp_BDdiv[[1]])[8:14]
predictors<-colnames(RES_PerResp[[1]])[11:20]
nPred<-length(predictors)

dataToPlot_allPreds_allSlices_allRespVars<-vector("list",length(RES_PerResp))

for(z in 1:length(RES_PerResp)){
RESP<-RES_PerResp[[z]]

dataToPlot_allPreds_allSlices<-vector("list",length=nslices)

	for(j in 1:10){ # for each slice from 5 to 50 Ma...

	slice<-sliceTimes[[j]]
	uniDat<-RESP[which(RESP$num==j),]

	# SEs to calc 95 CIs...
	slopesToPlot<-uniDat[,predictors]
	pvalsToPlot<-uniDat[,c(paste("p_",(1:nPred),sep=""))]
	SEsToPlot<-uniDat[,c(paste("SE_",(1:nPred),sep=""))]
	lamsToPlot<-uniDat[,c(paste("lam_",(1:nPred),sep=""))]

		# The 95% CI for each estimate of each var, plotted together
		dataToPlot_allPreds_perSlice<-vector("list",length(predictors))

		for (i in 1:(ntrees-fewerTrees)){ # so for a given tree IN a slice... taking each of *nPred* vars and calc 95%CI from SE
		#for (i in 1:ntrees){ 
			slopes_i<-slopesToPlot[i,]
			SEs_i<-SEsToPlot[i,]
			pVals_i<-pvalsToPlot[i,]
			lams_i<-lamsToPlot[i,]
			dataToPlot_k<-vector("list",length(predictors))
			for (k in 1:length(predictors)){
				color<-if(pVals_i[,k] <= 0.05){ "grey" } else { "red" } 
				lowHigh95_k<-cbind.data.frame(slopes_i[,k]-(1.96*SEs_i[,k]),slopes_i[,k],slopes_i[,k]+(1.96*SEs_i[,k]),color,predictors[k],lams_i[,k],sliceTimes[j],j,i)
				colnames(lowHigh95_k)<-c("low","mean","high","color","xVal","lam","time","slice","tree")
				dataToPlot_k[[k]]<-lowHigh95_k
			} # end *nPred* loop
			dataToPlot_allPreds_perSlice[[i]]<-do.call(rbind,dataToPlot_k)
		} # end 100 tree loop
		dataToPlot_allPreds_allSlices[[j]]<-do.call(rbind,dataToPlot_allPreds_perSlice)
	} # end 10 slice loop
dataToPlot_allPreds_allSlices_allRespVars[[z]]<-do.call(rbind,dataToPlot_allPreds_allSlices)
} # end 4 respVar loop


# count the NUMBER of significant runs per time slice per variable across 100 trees
# ==================================================================================
numSignif_perSlice_ALL<-vector("list",length(respVar))
lamMean_perSlice_ALL<-vector("list",length(respVar))
for(z in 1:length(respVar)){
dat<-dataToPlot_allPreds_allSlices_allRespVars[[z]]
	
	numSignif_perSlice<-vector("list",length=nslices)
	lamMean_perSlice<-vector("list",length=nslices)
	for(j in 1:nslices){
	slice<-dat[which(dat$slice==j),]

		numSignif<-c()
		lamMean<-c()
		for (k in 1:length(predictors)){
		predPerSlice<-slice[which(slice$xVal==predictors[k]),]
		num<-length(predPerSlice[which(predPerSlice$color=="grey"),][,1])
		numSignif[k]<-round((num/(ntrees-fewerTrees))*ntrees,0)
		
		lamMean[k]<-round(mean(predPerSlice[,"lam"]),2)
		}
		numSignif_perSlice[[j]]<-c(numSignif,sliceTimes[j],j)
		lamMean_perSlice[[j]]<-c(lamMean,sliceTimes[j],j)
	}
	RES1<-do.call(rbind,numSignif_perSlice)
	colnames(RES1)<-c(predictors,"time","slice")
	rownames(RES1)<-rep(respVar[z],nslices)
	
	RES2<-do.call(rbind,lamMean_perSlice)
	colnames(RES2)<-c(predictors,"time","slice")
	rownames(RES2)<-rep(respVar[z],nslices)

numSignif_perSlice_ALL[[z]]<-RES1
lamMean_perSlice_ALL[[z]]<-RES2
}


# PLOT the 95% CIs per slice for each respVar by predictor comparison across 100 trees
# ====================================================================================
corr="PAGEL"
sliceN="5to50Ma"
predictors<-colnames(RES_PerResp[[1]])[11:20]
nPred<-length(predictors)
respVar<-rownames(uniPGLS)
nResp<-length(respVar)
nTotal<-nPred*nResp
sliceTimes<-seq(-5,-70,-5)[1:10]

pdf(file=paste("cladeLevel",bbone,"_PGLSuni_",corr,"_explaining-DR-MRCA_timeSlices_",sliceN,"_95pCI_SCALED_rootwardWithSingle_10vars.pdf",sep=""),onefile=TRUE, width=(4*nResp),height=(4*nPred))
#pdf(file=paste("cladeLevel",bbone,"_PGLSuni_",corr,"_explaining_4varDR_timeSlices_",sliceN,"_95pCI_7vars.pdf",sep=""),onefile=TRUE, width=(4*nResp),height=(4*nPred))

#quartz(width=(6*4),height=(4*25))
layout(matrix(c(1:nTotal), nPred, nResp, byrow = FALSE), widths=rep(4,nTotal), heights=rep(3,nTotal))
par(oma = c(5,4,5,3) + 0.1, mar = c(4,1,1,1) + 0.1)

vertLwd<-0.5
lwdCI<-1.5
horizLine<-0
empirPointCol<-grey(0.3,alpha=0.5)
redCol<-rgb(1,0,0,alpha=0.3)
yLims1<-c(-1.5,1.5)

for(z in 1:length(respVar)){
uniDat_ALL<-dataToPlot_allPreds_allSlices_allRespVars[[z]]
signifDat_ALL<-numSignif_perSlice_ALL[[z]]
lamMean_ALL<-lamMean_perSlice_ALL[[z]]

	for(k in 1:length(predictors)){
	uniDat<-uniDat_ALL[which(uniDat_ALL$xVal==predictors[k]),]

	# plot base (no points)
	plot(formula(mean ~ time), data=uniDat, ylim=yLims1, ylab="",xlab="", yaxt="n",xaxt="n",type="n")
	axis(side=2,at=NULL,labels=TRUE)
	axis(side=1,at=NULL,labels=TRUE)
#	if(z==1 && k==1){ mtext(side=3, line=3, text="UNIVARIATE - non-log predictors")}
	if(z==1 && k==1){ mtext(side=3, line=3, text="UNIVARIATE - log predictors")}

	# add guide lines
	for(j in 1:nslices){
		abline(v=-5*j,lty=2, lwd=vertLwd, col=grey(0.6, alpha=0.5))
	}
		abline(h=horizLine,lty=1, lwd=1, col=grey(0.6, alpha=0.5))

# plot EMPIRICAL data
	dat_UNsignif<-uniDat[which(uniDat$color=="red"),]
	if(length(dat_UNsignif[,1])>0){
	plotCI(x=dat_UNsignif$time, add=TRUE,y=dat_UNsignif$mean,ui=dat_UNsignif$high,li=dat_UNsignif$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=redCol,scol=redCol,pch=1,font.lab=2,cex.axis=1.1,cex.lab=1.1)
	}
	dat_signif<-uniDat[which(uniDat$color=="grey"),]
	if(length(dat_signif[,1])>0){
	plotCI(x=dat_signif$time, add=TRUE,y=dat_signif$mean,ui=dat_signif$high,li=dat_signif$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=empirPointCol,scol=empirPointCol,pch=1,font.lab=2,cex.axis=1.1,cex.lab=1.1)
	}
	mtext(side=3,text=paste(respVar[z]," ~ "),font=2,adj=0,line=1)
	mtext(side=3,text=predictors[k],font=2,adj=0)
	
	# add numSignif data
	sigDat<-signifDat_ALL[,c(predictors[k],"time")]
	text(x=sigDat[,2], y = yLims1[2], cex=1.2,font=2, labels = sigDat[,1], col="dark grey")

	# add LAMBDA data
	lamDat<-lamMean_ALL[,c(predictors[k],"time")]
	text(x=lamDat[,2], y = yLims1[1], cex=1.2,font=2, labels = lamDat[,1], col="dark grey")

	} # end *nPred* loop

} # end 4 response var loop

dev.off()




### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
###############
# Code a multivariate one-- but as an ANOVA !!
# ======================================



# ==============
# do the PGLS-- MULTIVARIATE ANOVA
#######
# intialize
library(ape); library(phytools); library(picante); library(geiger); library(moments); library(nlme); library(caper)
library(foreach);library(doSNOW)

# which backbone?
bbone<- "NDexp" #"FBD" # 

# get the clade names
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery
allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

# start cluster
cl = makeCluster(100, type = 'MPI')
registerDoSNOW(cl)

#ranTrees<-c(2, 6, 10, 11, 12, 21, 23, 24, 25, 27, 29, 31, 33, 40, 41, 48, 53, 54, 59, 60, 63, 65, 67, 68, 72, 73, 74, 77, 82, 83, 84, 85, 86, 91, 94, 96, 98)
#toRun<-setdiff(1:100,ranTrees)
#foreach(i=toRun, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

ntrees=100
foreach(i=1:ntrees, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools', 'caper'), .combine=cbind, .verbose=TRUE) %dopar% {


# read back in slice phys
# ========================
slicePhys<-read.tree(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_CORRECT-sliceRootwardSingletons_slicePhy-5to70Ma_by5.trees",sep=""))

## multiply edges by 100 for corPagel
#modTrees<-vector("list",length(slicePhys))
#for (j in 1:length(slicePhys)){
#	modTree<-slicePhys[[j]]
#	modTree$edge.length<-modTree$edge.length*100
#	modTrees[[j]]<-modTree
#}

# Load in slice clade summaries and STANDARDIZE
# ========================================
results<-vector("list",length(allCladeSetNames))
for (q in 1: length(allCladeSetNames)){
	# Load in
	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[q]],"_cladeSTATS_withEcoTraits_andSingletons.txt",sep=""), header=TRUE)
	# select RESPONSE vars
	RESP<-cbind.data.frame(res[,c("DR_harm","BD_Div","MRCA","DR_skew")])
	# select PREDICTOR vars
	PRED_LOGvars<-c("HomeRange_km2", "GeoArea_km2", "DispDist_km", "lifemode1234", "GenLength_d") # "BM_final_kg", "", 
	PRED_NONvars<-c("trophic123", "activity123", "Lat_centroid", "Lon_centroid")#, "percentSamp")
	PRED_names<-c(paste("log",c("HomeRange_km2","GeoArea_km2", "DispersalIndex", "GenLength_d"),sep=""), PRED_NONvars)  

	PRED_all<-cbind.data.frame(log(res[,PRED_LOGvars[1:2]]),log(res[,PRED_LOGvars[3]]*res[,PRED_LOGvars[4]]),log(res[,PRED_LOGvars[5]]),res[,PRED_NONvars])
	#PRED_all<-cbind(log(res[,3:36]))
		# STANDARDIZE predictors as Z-SCORES
		#predScale<-scale(PRED_all, center=TRUE,scale=TRUE)

	# Join together RESP and PRED
	#results[[q]]<-cbind(log(RESP[,1:3]),RESP[,4],predScale)
	results[[q]]<-cbind.data.frame(rownames(res),log(RESP[,1:3]),RESP[,4],PRED_all)
	colnames(results[[q]])<-c("names","logTipDR","logBD_Div","logMRCA","DRskew",PRED_names)
}


# Run for ONE of 14 time slices 
# ===============================

for (j in 1:length(results)){
#j=2 # 10 Ma time slice
#j=7 # 35 Ma time slice
#j=9 # 45 Ma time slice
#j=12 # 60 Ma time slice
	#cladeData<-treedata(slicePhys[[j]],na.omit(results[[j]])) 
	#dat<-as.data.frame(cladeData$data)
cladeData<-comparative.data(phy=slicePhys[[j]],data=results[[j]], names.col="names",na.omit=TRUE)

# Setup empty dataframes to receive PGLS results 
# ==============================================
# MULTIVARIATE complex-- just subtracting one var each time and looking at DIFF IN AICc scores.

respVar<-c("logTipDR","logBD_Div","logMRCA","DRskew")
predictors<-PRED_names
nPred<-length(predictors)

	sumFull_params<-vector("list",length(respVar))
	sumFull_RSqs<-vector("list",length(respVar))

	fullRes<-vector("list",length(respVar))
	for(z in 1:length(respVar)){
	#z=1
		form<-as.formula(paste(respVar[z], " ~ ", predictors[1], "+", predictors[2], "+", predictors[3], "+", predictors[4], "+", predictors[5], "+", predictors[6], "+", predictors[7], "+", predictors[8], sep=""))  #"+", predictors[9],
		model_ALL<-pgls(form, data=cladeData, lambda=1.0)

		#form<-as.formula(paste(respVar[z], " ~ ", predictors[1], "+", predictors[2], "+", predictors[3], "+", predictors[4], "+", predictors[5], "+", predictors[6], "+", predictors[7], "+", predictors[8], sep=""))
		#model_9<-pgls(form, data=cladeData, lambda=1.0)

		form<-as.formula(paste(respVar[z], " ~ ", predictors[1], "+", predictors[2], "+", predictors[3], "+", predictors[4], "+", predictors[5], "+", predictors[6], "+", predictors[7], sep=""))
		model_8<-pgls(form, data=cladeData, lambda=1.0)

		form<-as.formula(paste(respVar[z], " ~ ", predictors[1], "+", predictors[2], "+", predictors[3], "+", predictors[4], "+", predictors[5], "+", predictors[6],  "+", predictors[8], sep=""))
		model_7<-pgls(form, data=cladeData, lambda=1.0)

		form<-as.formula(paste(respVar[z], " ~ ", predictors[1], "+", predictors[2], "+", predictors[3], "+", predictors[4], "+", predictors[5], "+", predictors[7], "+", predictors[8], sep=""))
		model_6<-pgls(form, data=cladeData, lambda=1.0)

		form<-as.formula(paste(respVar[z], " ~ ", predictors[1], "+", predictors[2], "+", predictors[3], "+", predictors[4], "+", predictors[6], "+", predictors[7], "+", predictors[8], sep=""))
		model_5<-pgls(form, data=cladeData, lambda=1.0)

		form<-as.formula(paste(respVar[z], " ~ ", predictors[1], "+", predictors[2], "+", predictors[3], "+", predictors[5], "+", predictors[6], "+", predictors[7], "+", predictors[8], sep=""))
		model_4<-pgls(form, data=cladeData, lambda=1.0)

		form<-as.formula(paste(respVar[z], " ~ ", predictors[1], "+", predictors[2], "+", predictors[4], "+", predictors[5], "+", predictors[6], "+", predictors[7], "+", predictors[8], sep=""))
		model_3<-pgls(form, data=cladeData, lambda=1.0)

		form<-as.formula(paste(respVar[z], " ~ ", predictors[1], "+", predictors[3], "+", predictors[4], "+", predictors[5], "+", predictors[6], "+", predictors[7], "+", predictors[8], sep=""))
		model_2<-pgls(form, data=cladeData, lambda=1.0)

		form<-as.formula(paste(respVar[z], " ~ ", predictors[2], "+", predictors[3], "+", predictors[4], "+", predictors[5], "+", predictors[6], "+", predictors[7], "+", predictors[8], sep=""))
		model_1<-pgls(form, data=cladeData, lambda=1.0)

		anovaCompare<-as.data.frame(anova.pgls(model_ALL, model_1, model_2, model_3, model_4, model_5, model_6, model_7, model_8))
		aicCompare<-AIC(model_ALL, model_1, model_2, model_3, model_4, model_5, model_6, model_7, model_8)

		deltaAIC<-(aicCompare[,2]-aicCompare[1,2]) # diff between model with variable removed and full model
		names(deltaAIC)<-c("ALL",predictors)
		ord<-sort(deltaAIC, decreasing=TRUE)
		nums<-1:length(aicCompare[,1])
		
		res<-cbind(ord,nums)
		rank<-res[match(c("ALL",predictors),rownames(res)),2]

			#bestName1<-names(which(res[,"nums"]==1))
			#bestNumber1<-match(bestName1,predictors)
			#sum<-summary(get(paste("model_",bestNumber1,sep="")))
		sumFull<-summary(model_ALL)
			tTable<-sumFull$coef
			RSq<-sumFull$r.squared
			adjRSq<-sumFull$adj.r.squared[1]
		
		respName<-rep(respVar[z],nPred+1)
		sumFull_params[[z]]<-cbind.data.frame(respName,tTable)
		sumFull_RSqs[[z]]<-c(respVar[z],RSq,adjRSq)

		resp<-rep(respVar[z],length(predictors)+1)
		pred<-c("ALL",predictors)
		tree<-i

		fullRes[[z]]<-cbind.data.frame(resp,pred,rank,deltaAIC,anovaCompare,aicCompare,tree)
		rownames(fullRes[[z]])<-1:length(aicCompare[,1])

	} # cycles 4 response vars
	fullResALL<-do.call(rbind,fullRes)

	sumFull_paramsALL<-do.call(rbind,sumFull_params)
	sumFull_RSqsALL<-as.data.frame(do.call(rbind,sumFull_RSqs))


	#corr<-"PAGEL"
	corr<-"BROWNIAN"
	sliceN<-allCladeSetNames[[j]] #"ALL" #"45Ma-only" #"10to60Ma" # # #"35Ma-only" ##"all" #"60Ma"#"all" #

	write.table(fullResALL, file=paste("cladeLevel",bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-DR-MRCA_perSlice_",sliceN,"_multiANOVA_Correct-rootwardWithSingle_8predsNoSamp.txt",sep=""))

	write.table(sumFull_paramsALL, file=paste("cladeLevel",bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-DR-MRCA_perSlice_",sliceN,"_multiANOVA_Correct-rootwardWithSingle_8predsNoSamp_tTable.txt",sep=""))
	write.table(sumFull_RSqsALL, file=paste("cladeLevel",bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-DR-MRCA_perSlice_",sliceN,"_multiANOVA_Correct-rootwardWithSingle_8predsNoSamp_RSqs.txt",sep=""))

	} # cycles each time slice

#uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs,uniPGLS_allSEs,uniPGLS_allLams)
#write.table(uniPGLS,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-",respVar[z],"_timeSlices_",sliceN,"_SCALED_uniVar_Correct_7preds.txt",sep=""))

} # cycles the 100 trees

stopCluster(cl)

q()

n

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# ===================================
# SUMMARIZE those results as TABLES (then visualize...)
# ===================================
# initialize
library(ape); library(phytools); library(picante); library(geiger); library(moments); library(nlme); library(caper)

dirname<-"/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/cladeLevel_SLICES_explainingDR"
setwd(dirname)

bbone<- "NDexp" #"FBD" # 
ntrees<-100

# get the clade names
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery
allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

sliceTimes<-seq(-5,-70,-5)

corr<-"BROWNIAN"
respVar<-c("logTipDR","logBD_Div","logMRCA","DRskew")
	PRED_LOGvars<-c("HomeRange_km2", "trophic123", "GeoArea_km2", "DispDist_km", "lifemode1234","GenLength_d","activity123") # "BM_final_kg", 
	PRED_NONvars<-c("Lat_centroid", "percentSamp") # "Lon_centroid", 
	PRED_names<-c(paste("log",c("HomeRange_km2","Trophic123", "GeoArea_km2", "GenLength_d", "Activity123","DispersalIndex"),sep=""), PRED_NONvars) # "BodyMass_kg", "DispDist_km", 
predictors<-PRED_names
nPred<-length(predictors)

# Load results back in
############
allRes_anova<-vector("list",(numSlices-1))
for(j in 1:(numSlices-1)){
	sliceN<-allCladeSetNames[[j]] #"ALL" #"45Ma-only" #"10to60Ma" # # #"35Ma-only" ##"all" #"60Ma"#"all" #

	res<-vector("list",ntrees)
	for(i in 1:ntrees){
		time<-sliceTimes[j]
		res[[i]]<-cbind(read.table(file=paste("cladeLevel",bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-DR-MRCA_perSlice_",sliceN,"_multiANOVA_Correct-rootwardWithSingle_8preds.txt",sep=""), header=TRUE),time)
	}
	allRes_anova[[j]]<-do.call(rbind,res)
}

allRes_tTable<-vector("list",(numSlices-1))
allRes_RSqs<-vector("list",(numSlices-1))
for(j in 1:(numSlices-1)){
	sliceN<-allCladeSetNames[[j]] #"ALL" #"45Ma-only" #"10to60Ma" # # #"35Ma-only" ##"all" #"60Ma"#"all" #

	resT<-vector("list",ntrees)
	resR<-vector("list",ntrees)
	for(i in 1:ntrees){
		tree<-i
		time<-sliceTimes[j]
		term<-c("INT",predictors)
		resT[[i]]<-cbind(term,read.table(file=paste("cladeLevel",bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-DR-MRCA_perSlice_",sliceN,"_multiANOVA_Correct-rootwardWithSingle_8preds_tTable.txt",sep=""), header=TRUE),time, tree)
		resR[[i]]<-cbind(read.table(file=paste("cladeLevel",bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-DR-MRCA_perSlice_",sliceN,"_multiANOVA_Correct-rootwardWithSingle_8preds_RSqs.txt",sep=""), header=TRUE),time, tree)
	}
	allRes_tTable[[j]]<-do.call(rbind,resT)
	allRes_RSqs[[j]]<-do.call(rbind,resR)

}


# per slice, look at best variables
########

allRes_perResp_perSlice<-vector("list",(numSlices-1))
for(j in 1:(numSlices-1)){

	Anova_all<-allRes_anova[[j]]
	tTable_all<-allRes_tTable[[j]]
	RSqs_all<-allRes_RSqs[[j]]

	allRes_perResp<-vector("list",length(respVar))
	for(z in 1:length(respVar)){
		Anova_z<-Anova_all[which(Anova_all$resp==respVar[z]),]
		tTable_z<-tTable_all[which(tTable_all$respName==respVar[z]),]
		RSqs_z<-RSqs_all[which(RSqs_all$V1==respVar[z]),]

		num1<-c(); num2<-c(); num3<-c()
		delAIC_pred<-vector("list",length(predictors))
		Est_pred<-vector("list",length(predictors))
		for(k in 1:length(predictors)){
			Anova_pred<-Anova_z[which(Anova_z$pred==predictors[k]),]
			tTable_pred<-tTable_z[which(tTable_z$term==predictors[k]),]

			num1[k]<-length(Anova_pred[which(Anova_pred$rank==1),1])
			num2[k]<-length(Anova_pred[which(Anova_pred$rank==2),1])
			num3[k]<-length(Anova_pred[which(Anova_pred$rank==3),1])

			Mean<-mean(Anova_pred[,"deltaAIC"])
			delAIC_pred[[k]]<-c(Mean,quantile(Anova_pred[,"deltaAIC"], c(0.025,0.975)))

			MeanEst<-mean(tTable_pred[,"Estimate"])
			direction<-c()
			for(i in 1:ntrees){
				if(tTable_pred[i,"Estimate"] >0){ direction[i]<- 1} else { direction[i]<- 0}
			}
			numPos<-sum(direction)
			Est_pred[[k]]<-c(MeanEst,quantile(tTable_pred[,"Estimate"], c(0.025,0.975)),numPos)

		}
		delAIC_predALL<-cbind.data.frame(num1,num2,num3,do.call(rbind,delAIC_pred))
		rownames(delAIC_predALL)<-predictors
	
		Est_predALL<-do.call(rbind,Est_pred)
		rownames(Est_predALL)<-predictors

		allRes<-cbind.data.frame(respVar[z],delAIC_predALL,Est_predALL)
		colnames(allRes)<-c("respVar","num1", "num2", "num3", "delAIC_mean", "delAIC_low95", "delAIC_up95", "Est_mean", "Est_low95", "Est_up95", "Est_numPos")

	allRes_perResp[[z]]<-allRes
	}
	allRes_perRespALL<-do.call(rbind,allRes_perResp)

slice<-rep(sliceTimes[j],length(allRes_perRespALL[,1]))
allRes_perResp_perSlice[[j]]<-cbind.data.frame(allRes_perRespALL,slice)
}
allRes_TOTAL<-do.call(rbind,allRes_perResp_perSlice)

# then look at PER RESPONSE - PER SLICE again, since I will be PLOTTING it that way::
# TO PLOT
corr="BROWNIAN"
sliceN="5to65Ma"
nResp<-length(respVar)
nSlices<-numSlices-1
nTotal<-nSlices*nResp

library(viridis)
cols<-viridis(nPred)

# BARPLOT
pdf(file=paste("cladeLevel",bbone,"_PGLS_",corr,"_multiANOVA_explaining-DR-MRCA_timeSlices_",sliceN,"_rootwardWithSingle_8pred_BARplot.pdf",sep=""),onefile=TRUE, width=(4*nSlices),height=(4*nResp))

#quartz(width=(6*4),height=(4*25))
layout(matrix(c(1:nTotal), nResp, nSlices, byrow = TRUE), widths=rep(4,nTotal), heights=rep(3,nTotal))
par(oma = c(5,4,5,3) + 0.1, mar = c(4,1,1,1) + 0.1)

for(z in 1:length(respVar)){
	resResp<-allRes_TOTAL[which(allRes_TOTAL$respVar==respVar[z]),] #c(1:5, 11, 12)]

	for(j in (numSlices-1):1){
		resResp_perSlice<-resResp[which(resResp$slice==sliceTimes[j]),]

		barplot(resResp_perSlice$delAIC_mean, ylim=c(-5,65), col=cols)#, horiz=TRUE)#, names.arg=rownames(resResp_perSlice), cex.names=0.5)
		text(x=seq(from=0.5,to=9,length.out=8),y=0.5,labels=predictors, adj=c(0,1), srt=90, font=2, cex=1.2)
		mtext(side=3, text=paste(respVar[z], " slice ",sliceTimes[j]), font=2)

	} # cycle 13 slices
} # cycle 4 respVars
dev.off()


library(viridis)
# LINE-PLOT
pdf(file=paste("cladeLevel",bbone,"_PGLS_",corr,"_multiANOVA_explaining-DR-MRCA_timeSlices_",sliceN,"_rootwardWithSingle_8pred_LINEplot.pdf",sep=""),onefile=TRUE, width=(4*1),height=(4*nResp))

#quartz(width=(6*4),height=(4*25))
layout(matrix(c(1:nResp), 4, 1 , byrow = TRUE), widths=rep(4,nTotal), heights=rep(3,nTotal))
par(oma = c(5,4,5,3) + 0.1, mar = c(4,1,1,1) + 0.1)

yMax<-65

for(z in 1:length(respVar)){
	resResp<-allRes_TOTAL[which(allRes_TOTAL$respVar==respVar[z]),] #c(1:5, 11, 12)]

	cols<-viridis(nPred)
	for(k in 1:nPred){
		resResp_perPred<-resResp[c(k+((numSlices-2):0)*8),]

		if(k==1) { 	plot(delAIC_mean ~ slice, data=resResp_perPred, type="l", col=cols[k], lwd=3, ylim=c(-5,yMax), xlab="")
					mtext(side=2, line=3, cex=0.7, text="deltaAIC (model decrease without predictor)")
		} else { points(delAIC_mean ~ slice, data=resResp_perPred, type="l", col=cols[k], lwd=3) 
		}
		
	} # cycles 8 preds
	legend(x=-65, y=yMax, legend=predictors, col=cols, lty=1, lwd=3, cex=0.8, title=paste("Predictors of ",respVar[z],sep=""))
	mtext(side=3, text=paste(respVar[z], " all slices",sep=""), font=2)

	if(z==4) {	mtext(side=1, line=3, cex=0.7, text="time before present (Ma)")
		} else { NULL }

} # cycle 4 respVars

dev.off()




## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# ==================
# RAW DATA, scatter plots of clades at different slices.
#####
dirname<-"/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/cladeLevel_SLICES_explainingDR"
setwd(dirname)
# intialize
library(ape); library(phytools); library(picante); library(geiger); library(moments); library(nlme); library(caper)

# which backbone?
bbone<- "NDexp" #"FBD" # 

# get the clade names
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery
allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

ntrees=100

i=1
# Load in slice clade summaries and STANDARDIZE
# ========================================
results<-vector("list",length(allCladeSetNames))
for (q in 1: length(allCladeSetNames)){
	# Load in
	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[q]],"_cladeSTATS_withEcoTraits_andSingletons.txt",sep=""), header=TRUE)
	# select RESPONSE vars
	RESP<-cbind.data.frame(res[,c("DR_harm","BD_Div","MRCA","DR_skew")])
	# select PREDICTOR vars
	PRED_LOGvars<-c("HomeRange_km2", "GeoArea_km2", "DispDist_km", "lifemode1234", "GenLength_d") # "BM_final_kg", "", 
	PRED_NONvars<-c("trophic123", "activity123", "Lat_centroid", "Lon_centroid", "percentSamp")
	PRED_names<-c(paste("log",c("HomeRange_km2","GeoArea_km2", "DispersalIndex", "GenLength_d"),sep=""), PRED_NONvars)  

	PRED_all<-cbind.data.frame(log(res[,PRED_LOGvars[1:2]]),log(res[,PRED_LOGvars[3]]*res[,PRED_LOGvars[4]]),log(res[,PRED_LOGvars[5]]),res[,PRED_NONvars])
	#PRED_all<-cbind(log(res[,3:36]))
		# STANDARDIZE predictors as Z-SCORES
		#predScale<-scale(PRED_all, center=TRUE,scale=TRUE)

	# Join together RESP and PRED
	#results[[q]]<-cbind(log(RESP[,1:3]),RESP[,4],predScale)
	results[[q]]<-cbind.data.frame(rownames(res),log(RESP[,1:3]),RESP[,4],PRED_all)
	colnames(results[[q]])<-c("names","logTipDR","logBD_Div","logMRCA","DRskew",PRED_names)
}

# COMPARE scatterplots across slices
# PLOT
#####
predictors<-PRED_names
nPred<-length(predictors)

respVar<-c("logTipDR","logTipDR","logBD_Div","logBD_Div","logMRCA","logMRCA","DRskew")

nResp<-length(respVar)
#nSlices<-numSlices
sliceN="5to60Ma-by20"
whichSlices<-c(1, 4, 8, 12)
nSlices<-length(whichSlices)
nTotal<-nSlices*nResp

library(viridis)
cols<-viridis(nPred)

yLabels<-c("LOG tip DR mean","tip DR mean","LOG BD_Div", "BD_Div","LOG MRCA","MRCA","DR_skew")


# BARPLOT
pdf(file=paste("cladeLevel",bbone,"_explaining-DR-MRCA_timeSlices_",sliceN,"_SCATTERS_vsHomeRange.pdf",sep=""),onefile=TRUE, width=(4*nSlices),height=(4*nResp))

#quartz(width=(6*4),height=(4*25))
layout(matrix(c(1:nTotal), nResp, nSlices, byrow = TRUE), widths=rep(4,nTotal), heights=rep(3,nTotal))
par(oma = c(5,4,5,3) + 0.1, mar = c(4,1,1,1) + 0.1)

for(z in 1:length(respVar)){

	for(q in rev(whichSlices)){
		res<-results[[q]]

		if(z==2 || z==4 || z==6){
		form<-as.formula(paste("exp(",respVar[z], ") ~ ", "logHomeRange_km2", sep=""))
		plot(form, data=res)
		} else {
		form<-as.formula(paste(respVar[z], " ~ ", "logHomeRange_km2", sep=""))
		plot(form, data=res)
		}
		mtext(side=3, text=allCladeSetNames[[q]], line=-1, font=2)
		
		if(q==12){
		mtext(side=2, text=yLabels[z], line =2, font=2)
		} else { NULL }
	}
}

dev.off()

















########
# TEST CODE

# Work this example at 70 Ma...
root=max(node.age(mamPhy)$ages)

treesTip<-treeSlice(mamPhy, slice=root-70, trivial=TRUE)
treesRoot<-treeSlice(mamPhy, slice=root-70, trivial=TRUE, orientation="rootwards")

	cladeSet<-treesTip
	newTipNames<-paste(i,"_",j,"_",c(1:length(cladeSet)),sep="")

	for (k in 1:length(cladeSet)){
		cladeSp<-cladeSet[[k]]$tip.label
		if(length(cladeSp)==1){ 
			treesRoot$tip.label[which(treesRoot$tip.label==cladeSp)]<-newTipNames[k]
		} else {
		node <- getMRCA(mamPhy, cladeSp) # find the MRCA node of those species
		treesRoot$tip.label[which(treesRoot$tip.label==node)]<-newTipNames[k]
		}
	}

treesRoot_orig<-treeSlice(mamPhy, slice=root-70, trivial=TRUE, orientation="rootwards")

pdf(file="test_treeSlice_tree1_70Ma_rootwards_reLabel.pdf", width=5, height=10, onefile=TRUE)
plot(treesRoot_orig)
axisPhylo()

plot(treesRoot)
axisPhylo()

dev.off()
	# nice, that works.


pdf(file="test_treeSlice_tree1_all_nodes.pdf", width=10, height=40)
plot(mamPhy, show.tip.label=TRUE, cex=0.1)
nodelabels(cex=0.1)
axisPhylo()
dev.off()




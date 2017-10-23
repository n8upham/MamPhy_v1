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



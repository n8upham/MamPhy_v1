#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy v1 -- Upham et al. 2017
###
# Figure 3 - time slices to explain clade richness
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run and visualize time slice analyses
######
# - In parallel, load in 1 tree of 100
# - Calculate ES and DR on per-tip basis
# - Make time slices every 5 Ma to create clades
# - Calculate per-slice, per-clade summary values 
# - Relabel the backbone "slicePhys" uniting  time slice clades
# - Calculate per-slice, per-clade summary values-- including singletons
# - end 100-tree loop


# packages
library(moments); library(nlme); library(ape); library(picante); library(phytools); library(geiger)
library(foreach);library(doSNOW)

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

# load in 1 tree of 100
mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_nexus.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 

# load in DR per species
res1<-read.table(file=paste(bbone,"_sample100_",i,"_DRtips.txt",sep=""))


# load back in sliceClades (pre-calculated)
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery

root=max(node.age(mamPhy)$ages)

allCladeSetNames<-vector("list",length=numSlices)
allCladeSets<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
	# re-load
	allCladeSets[[j]]<-read.tree(file=paste(bbone,"_sample100_",i,"_timeslice_cladeTreesWithSingletons_",(sliceEvery*j),"Ma.trees",sep=""))
	}

#===================================
# Calculate per-SLICE, per-clade summary values 
#===================================

# get node times for tree
btimes<-branching.times(mamPhy)

# yule function
#ymle = function(tree){ (.subset2(tree,3)-1L)/sum(.subset2(tree,2)) } # this take the # of number of nodes in a tree (minus 1) / sum of branch lengths.
ymle = function(tree){ (.subset2(tree,2)-1L)/sum(.subset2(tree,4)) } # this take the # of number of nodes in a tree (minus 1) / sum of branch lengths.

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
	
	for (k in 1:length(cladeSet)){
	cladeSp<-cladeSet[[k]]$tip.label
		x<-res1[match(cladeSp,rownames(res1)),"DR"]
		DR_harm[k,] <- 1/(mean(1/x))
		percentSamp[k,] <- length(which(cladesDR[match(cladeSp,cladesDR$tiplabel),"samp"]=="sampled"))/length(cladeSp)
		richness[k,] <- length(cladeSp)
	if(length(cladeSp) > 1){
		DR_cv[k,] <- (sd(x)/mean(x))*100
		DR_skew[k,] <- skewness(x)
		DR_kurt[k,] <- kurtosis(x)
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
		BD_Turn[k,]<-bd$para[[1]]
		# BD Mag and Sand
		cladeSet[[k]]$root.edge<-0
	    BD.ms0[k,]<-bd.ms(phy=cladeSet[[k]], missing=0, epsilon=0, crown=TRUE) # Assuming no extinction
    	BD.ms0p5[k,]<-bd.ms(phy=cladeSet[[k]], missing=0, epsilon=0.5, crown=TRUE) # Assuming medium extinction 
     	BD.ms0p9[k,]<-bd.ms(phy=cladeSet[[k]], missing=0, epsilon=0.9, crown=TRUE) # Assuming high extinction
		}
	}

	res2<-cbind.data.frame(DR_harm, DR_skew, DR_kurt, DR_cv, percentSamp, richness, MRCA, PB_Div, BD_Lam, BD_Mu, BD_Div, BD_Turn, BD.ms0, BD.ms0p5, BD.ms0p9, i, j*-5)

	colnames(res2)<-c("DR_harm","DR_skew", "DR_kurt", "DR_cv", "percentSamp", "richness", "MRCA", "PB_Div", "BD_Lam", "BD_Mu", "BD_Div", "BD_Turn", "BD.ms0", "BD.ms0p5", "BD.ms0p9", "tree", "time")

	rownames(res2)<-paste(i,"_",j,"_",c(1:length(cladeSet)),sep="") # same as the slicePhy names
	
	write.table(res2,paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"_cladeSTATS_Age-n-Rate_andSingletons.txt",sep=""))
} # end 1-tree loop across all 14 slices (every 5 Ma)


} # end 100-tree loop



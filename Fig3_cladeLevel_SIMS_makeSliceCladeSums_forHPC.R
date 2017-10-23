#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy v1 -- Upham et al. 2017
###
# Figure 3 - time slices to explain clade richness
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# For SIMULATIONS, calculate per-slice, per-clade summary values 
######
# - load back in 1 tree of 100
# - Calculate ES and DR on per-tip basis
# - Calculate per-slice, per-clade summary values-- including singletons


# packages
library(moments); library(nlme); library(ape); library(picante); library(phytools); library(geiger)
library(foreach);library(doSNOW)
source("DR_functions.R")

# open cluster for parallel processing
cl = makeCluster(100, type = 'MPI', outfile="")
registerDoSNOW(cl)

# start parallel loop
ntrees<-100
foreach(i=1:ntrees, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

# loop through 3 categories of simulation
sims<-c("mamPhyE","lowE_0p2","highE_0p8")
bbone<-"NDexp"

for(z in 1:length(sims)){

## read in 1 of 100 sumulated full trees
simPhy<-read.tree(paste("MamPhy_SIMS_",sims[z],"_",bbone,"_tree",i,".tre",sep=""))
tree1=scan(paste("MamPhy_SIMS_",sims[z],"_",bbone,"_tree",i,".tre",sep=""), what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 


# Calculate ES and DR on per-tip basis
# ==================
# gives pairwise clade matrix from CAIC function
clade_matrix = readCAIC(tree1)

# calculate and write to file
DR = 1/ES_v2(clade_matrix)
ES = ES_v2(clade_matrix)
res = cbind.data.frame(DR,ES)
res1 = res[order(rownames(res)),]

write.table(res1, file=paste("MamPhy_SIMS_",sims[z],"_",bbone,"_sample100_",i,"_DRtips.txt",sep=""))
res1<-read.table(file=paste("MamPhy_SIMS_",sims[z],"_",bbone,"_sample100_",i,"_DRtips.txt",sep=""))


# load back in sliceClades (pre-calculated)
# =========================
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery

root=max(node.age(simPhy)$ages)

allCladeSetNames<-vector("list",length=numSlices)
allCladeSets<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
	# re-load
	allCladeSets[[j]]<-read.tree(file=paste("MamPhy_SIMS_",sims[z],"_",bbone,"_sample100_",i,"_timeslice_cladeTrees_",(sliceEvery*j),"Ma.trees",sep=""))
	}


#===================================
# Calculate per-SLICE, per-clade summary values 
#===================================

# get node times for tree
btimes<-branching.times(simPhy)

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
	#percentSamp<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
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
		#percentSamp[k,] <- length(which(cladesDR[match(cladeSp,cladesDR$tiplabel),"samp"]=="sampled"))/length(cladeSp)
		richness[k,] <- length(cladeSp)
	if(length(cladeSp) > 1){
		DR_cv[k,] <- (sd(x)/mean(x))*100
		DR_skew[k,] <- skewness(x)
		DR_kurt[k,] <- kurtosis(x)
		node <- getMRCA(simPhy, cladeSp)
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

	res2<-cbind.data.frame(DR_harm, DR_skew, DR_kurt, DR_cv, richness, MRCA, PB_Div, BD_Lam, BD_Mu, BD_Div, BD_Turn, BD.ms0, BD.ms0p5, BD.ms0p9, i, j*-5)

	colnames(res2)<-c("DR_harm","DR_skew", "DR_kurt", "DR_cv", "richness", "MRCA", "PB_Div", "BD_Lam", "BD_Mu", "BD_Div", "BD_Turn", "BD.ms0", "BD.ms0p5", "BD.ms0p9", "tree", "time")

	rownames(res2)<-paste(i,"_",j,"_",c(1:length(cladeSet)),sep="") # same as the slicePhy names
	
	write.table(res2,paste("MamPhy_SIMS_",sims[z],"_",bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"_cladeSTATS_Age-n-Rate_andSingletons.txt",sep=""))
} # end 1-tree loop across all 14 slices (every 5 Ma)

} # end 3 simulation loop

} # end 100-tree loop






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy v1 -- Upham et al. 2017
###
# Figure 3 - time slices to explain clade richness
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run MamPhy time slice analyses
######
# - In parallel, load in 1 tree of 100
# - Calculate ES and DR on per-tip basis
# - Make time slices every 5 Ma to create clades
# - Relabel the backbone "slicePhys" uniting  time slice clades
# - Calculate per-slice, per-clade summary values-- including singletons
# - end 100-tree loop

# Run MamPhy PGLS -- aim to ** explain variation RICHNESS ** using timeslice clade summary stats 
######
# - In parallel, read back in slice phylo backbones (14 slices at 5 Ma intervals)
# - Load in age and rate slice summaries for 1 of 100 trees (standardize predictors)
# - Setup empty dataframes to receive PGLS results 
# - Nest the loop of *10* predictor vars
# - Do univariate and multivariate analyses


# For SIMULATIONS, do time slice analyses
######
# - Generate simulations of the full MamPhy for testing null expectations.
# - Relabel the backbone "slicePhys" uniting  time slice clades

# For SIMULATIONS, calculate per-slice, per-clade summary values 
######
# - load back in 1 tree of 100
# - Calculate ES and DR on per-tip basis
# - Calculate per-slice, per-clade summary values-- including singletons


# For SIMULATIONS, run PGLS to ** explain variation RICHNESS ** using timeslice clade summary stats 
######
# - In parallel, read back in slice phylo backbones (14 slices at 5 Ma intervals)
# - Load in age and rate slice summaries for 1 of 100 trees (standardize predictors)
# - Setup empty dataframes to receive PGLS results 
# - Nest the loop of *10* predictor vars
# - Do univariate and multivariate analyses



# For SIMS-- summarize and re-load all

# For ALL-- Make the actual data SCATTER plots...
# For ALL-- Also get the UNIVARIATE PENDANTS going...

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
#tipDataAll<-read.table(file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)

#==================
# Load in 1 tree of 100
#==================
mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_nexus.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 
#write.tree(mamPhy,file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""))
tree1=scan(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""), what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 


##==================
## Calculate ES and DR on per-tip basis
##==================
# gives pairwise clade matrix from CAIC function
clade_matrix = readCAIC(tree1)

# calculate and write to file
DR = 1/ES_v2(clade_matrix)
ES = ES_v2(clade_matrix)
res = cbind.data.frame(DR,ES)
res1 = res[order(rownames(res)),]

write.table(res1, file=paste(bbone,"_sample100_",i,"_DRtips.txt",sep=""))
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

# write slice phys
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
allCladeSets<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSets[[j]]<-read.tree(file=paste(bbone,"_sample100_",i,"_timeslice_cladeTreesWithSingletons_",(sliceEvery*j),"Ma.trees",sep=""))
	}

#===================================
# Calculate per-SLICE, per-clade summary values 
#===================================

# get node times for tree
btimes<-branching.times(mamPhy)

# yule function
ymle = function(tree){ (.subset2(tree,3)-1L)/sum(.subset2(tree,2)) } # this take the # of number of nodes in a tree (minus 1) / sum of branch lengths.
# use this one if RE-LOADING the slice clades::
#ymle = function(tree){ (.subset2(tree,2)-1L)/sum(.subset2(tree,4)) } # this take the # of number of nodes in a tree (minus 1) / sum of branch lengths.

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



###########
#######


# do the PGLS - UNIVAR and MULTIVAR
# ==================================
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

ntrees=100
foreach(i=1:ntrees, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

# read back in slice phys
# ========================
slicePhys<-read.tree(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_CORRECT-sliceRootwardSingletons_slicePhy-5to70Ma_by5.trees",sep=""))


# Load in slice clade summaries and STANDARDIZE
# ========================================
results<-vector("list",length(allCladeSetNames))
for (q in 1: length(allCladeSetNames)){
	# Load in
	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[q]],"_cladeSTATS_Age-n-Rate_andSingletons.txt",sep=""), header=TRUE)
	# select RESPONSE vars
	RESP<-log(res[,"richness"])
	# select PREDICTOR vars
	PRED<-res[,c(1:2,4:5,7:12)] # excluding KURTOSIS, excluding the BD.ms vars
		# STANDARDIZE predictors as Z-SCORES
		predScale<-scale(PRED, center=TRUE,scale=TRUE)

	# Join together RESP and standardized PRED
	results[[q]]<-cbind(RESP,predScale[,1:length(PRED)])
	colnames(results[[q]])<-c("logRichness",colnames(PRED))
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

predictors<-colnames(PRED) # focusing on continuous vars, trophic vars
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

uniPGLS<-data.frame(matrix(NA, nrow = 1, ncol = nCols*5),row.names="logRichness")
colnames(uniPGLS)<-c(colnames(uniPGLS_allInts),colnames(uniPGLS_allSlopes),colnames(uniPGLS_allPs),colnames(uniPGLS_allSEs),colnames(uniPGLS_allLams))

# MULTIVARIATE next
sliceTimes<-seq(-5,-70,-5)

multiPGLS<-data.frame(matrix(NA, nrow = 1, ncol = 13),row.names="logRichness")
colnames(multiPGLS)<-c("int","MRCA","DR_harm","DR_skew","SE1","SE2","SE3","SE4","lam","Pval1","Pval2","Pval3","Pval4")
multiPGLS_Per<-data.frame(matrix(NA, nrow = 1, ncol = 16),row.names="logRichness")
colnames(multiPGLS_Per)<-c("int","MRCA","DR_harm","DR_skew","percentSamp","SE1","SE2","SE3","SE4","SE5","lam","Pval1","Pval2","Pval3","Pval4","Pval5")


# Nest the loop of *nCols* predictor vars
# ===============================

# UNIVARIATE
	for(k in 1:length(predictors)){

		form<-as.formula(paste("logRichness", " ~ ", predictors[k], sep=""))
#		fit1<-gls(form, data=dat, method="ML")
#		fit1<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
		for (p in c(0.5,seq(0,1,by=0.01))) {possibleError <- tryCatch(
		      gls(form, correlation=corPagel(value=p,phy=cladeData$phy), data=dat, method="ML"),
		      error=function(e) e)
		if(inherits(possibleError, "gls")) break		
		if(inherits(possibleError, "error")) next}
		fit1<-possibleError

		sum1<-summary(fit1)
		
		uniPGLS_allInts[j,k]<-round(sum1$tTable[1], digits=6)
		uniPGLS_allSlopes[j,k]<-round(sum1$tTable[2], digits=6)
		uniPGLS_allSEs[j,k]<-round(sum1$tTable[4], digits=6)
		uniPGLS_allPs[j,k]<-round(sum1$tTable[8], digits=6)
		uniPGLS_allLams[j,k]<-round(sum1$modelStruct[[1]][[1]], digits=6)		

		} # cycles each of *nCols* predictors

		uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs,uniPGLS_allSEs,uniPGLS_allLams)[j,]

	#corr<-"NO_TREE"
	#corr<-"BROWNIAN"
	corr<-"PAGEL"
	sliceN<-allCladeSetNames[[j]] #"ALL" #"45Ma-only" #"10to60Ma" # # #"35Ma-only" ##"all" #"60Ma"#"all" #

	write.table(uniPGLS,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-Richness_perSlice_",sliceN,"_SCALED_uniVar_Correct-rootwardWithSingle_10preds.txt",sep=""))

# MULTIVARIATE
# 3 var
	form<-(logRichness ~ MRCA + DR_harm + DR_skew)
	#fit1<-gls(form2, data=dat, method="ML")
	#fit1<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
		for (p in c(0.5,seq(0,1,by=0.01))) {possibleError <- tryCatch(
		      gls(form, correlation=corPagel(value=p,phy=cladeData$phy), data=dat, method="ML"),
		      error=function(e) e)
		if(inherits(possibleError, "gls")) break		
		if(inherits(possibleError, "error")) next}
		fit1<-possibleError

		sum1<-summary(fit1)

	for(k in 1:8){
	multiPGLS[k]<-round(sum1$tTable[k],digits=6)
	}
	multiPGLS[9]<-round(sum1$modelStruct[[1]][1],digits=6)
	multiPGLS[10]<-round(sum1$tTable[13], digits=6)
	multiPGLS[11]<-round(sum1$tTable[14], digits=6)
	multiPGLS[12]<-round(sum1$tTable[15], digits=6)
	multiPGLS[13]<-round(sum1$tTable[16], digits=6)

	write.table(multiPGLS,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-Richness_perSlice_",sliceN,"_SCALED_multiVar_Correct-rootwardWithSingle_age-DRmean-DRskew.txt",sep=""))

# MULTIVARIATE
# 3 var + percentSampling
	form2<-(logRichness ~ MRCA + DR_harm + DR_skew + percentSamp)
	#fit2<-gls(form2, data=dat, method="ML")
	#fit2<-gls(form2, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
		for (p in c(0.5,seq(0,1,by=0.01))) {possibleError <- tryCatch(
		      gls(form2, correlation=corPagel(value=p,phy=cladeData$phy), data=dat, method="ML"),
		      error=function(e) e)
		if(inherits(possibleError, "gls")) break		
		if(inherits(possibleError, "error")) next}
		fit2<-possibleError

		sum2<-summary(fit2)

	for(k in 1:10){
	multiPGLS_Per[k]<-round(sum2$tTable[k],digits=6)
	}
	multiPGLS_Per[11]<-round(sum2$modelStruct[[1]][1],digits=6)
	multiPGLS_Per[12]<-round(sum2$tTable[16], digits=6)
	multiPGLS_Per[13]<-round(sum2$tTable[17], digits=6)
	multiPGLS_Per[14]<-round(sum2$tTable[18], digits=6)
	multiPGLS_Per[15]<-round(sum2$tTable[19], digits=6)
	multiPGLS_Per[16]<-round(sum2$tTable[20], digits=6)

	write.table(multiPGLS_Per,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-Richness_perSlice_",sliceN,"_SCALED_multiVar_Correct-rootwardWithSingle_age-DRmean-DRskew-perSamp.txt",sep=""))

	} # cycles each time slice

} # cycles the 100 trees



# ==============
# Visualize 
#======================================
# UNIVARIATE-- adding in the SIMULATIONS too.
#======================================
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
dirname="/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/cladeLevel_SLICES_explainingRichness"
setwd(dirname)

# ~~~~~~~~~~~~~~~~~~~
# EMPIRICAL
# load back in the results
# =========================
ntrees=100
whichTrees<-c(100, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 1, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 2, 30, 32, 33, 34, 35, 36, 37, 38, 39, 3, 40, 41, 42, 44, 45, 46, 47, 48, 49, 4, 50, 51, 52, 53, 55, 56, 57, 58, 59, 5, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 6, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 7, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 8, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 9)
fewerTrees<-ntrees-length(whichTrees)
nslices=14 # analyze all 14 slices, 5 Ma to 70 Ma 
corr="PAGEL"

exp_N<-vector("list",length=nslices)

for(j in 1:nslices){
sliceN<-allCladeSetNames[[j]]

	N<-vector("list",length=ntrees)
	
	for(i in whichTrees){
	#for(i in 1:ntrees){
	#uniPGLS<-read.table(file=paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-4varDR_perSlice_",sliceN,"_SCALED_uniVar_Correct_7preds.txt",sep=""), header=TRUE)
	uniPGLS<-read.table(file=paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-Richness_perSlice_",sliceN,"_SCALED_uniVar_Correct-rootwardWithSingle_10preds.txt",sep=""), header=TRUE)
	tree<-i
	slice<-sliceTimes[[j]]
	num<-j

	N[[i]]<-cbind(uniPGLS[1,],slice,num,tree)
	}

exp_N[[j]]<-do.call(rbind,N)
}

respVar<-"logRichness"


# calculate 95% CIs (and significance) for the slopes per tree & slice
# ====================================================================
#predictors<-colnames(exp_BDdiv[[1]])[8:14]
predictors<-colnames(exp_N[[1]])[11:20]
nPred<-length(predictors)

RESP<-exp_N

dataToPlot_allPreds_allSlices<-vector("list",length=nslices)

for(j in 1:nslices){ # for each slice from 5 to 70 Ma...

slice<-sliceTimes[[j]]
uniDat<-RESP[[j]]

# SEs to calc 95 CIs...
slopesToPlot<-uniDat[,predictors]
pvalsToPlot<-uniDat[,c(paste("p_",(1:nPred),sep=""))]
SEsToPlot<-uniDat[,c(paste("SE_",(1:nPred),sep=""))]
lamsToPlot<-uniDat[,c(paste("lam_",(1:nPred),sep=""))]

	# The 95% CI for each estimate of each var, plotted together
	dataToPlot_allPreds_perSlice<-vector("list",length(predictors))

	for (i in 1:(ntrees-fewerTrees)){ # so for a given tree IN a slice... taking each of *nPred* vars and calc 95%CI
	#for (i in 1:ntrees){ # so for a given tree IN a slice... taking each of *nPred* vars and calc 95%CI
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
dataToPlot_allPreds_allSlices_allRespVars<-do.call(rbind,dataToPlot_allPreds_allSlices)


# count the NUMBER of significant runs per time slice per variable across 100 trees
# ==================================================================================
dat<-dataToPlot_allPreds_allSlices_allRespVars
		
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
	rownames(RES1)<-rep(respVar,nslices)
	
	RES2<-do.call(rbind,lamMean_perSlice)
	colnames(RES2)<-c(predictors,"time","slice")
	rownames(RES2)<-rep(respVar,nslices)

numSignif_perSlice_ALL<-RES1
lamMean_perSlice_ALL<-RES2


# **SIMULATIONS**
# ~~~~~~~~~~~~~~~~
# load back in the results
# =========================
# set NEW dir
dirname="/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/cladeLevel_SLICES_explainingRichness/_res_SIMS_explainingRichness"
setwd(dirname)

ntrees=100
	# for UNI - mamPhyE -- 72 trees
	whichTrees<-c(100,10,11,12,13,18,1,20,21,22,23,26,27,28,29,2,30,32,33,35,36,38,39,3,40,41,43,44,45,46,47,48,4,50,51,52,53,54,55,56,57,58,5,60,61,63,64,68,6,70,71,72,73,75,76,77,79,7,80,81,82,83,84,86,8,90,91,93,95,96,99,9)
	# for UNI - lowE_0p2 -- 50? trees
	# for UNI - highE_0p8 -- 0 trees...

#whichTreesAll<-list(whichTrees_1,whichTrees_2,whichTrees_3)

sims<-c("mamPhyE","lowE_0p2","highE_0p8")

dataToPlot_allPreds_allSlices_allRespVars_SIMS<-vector("list",length(sims))
numSignif_perSlice_ALL_SIMS<-vector("list",length(sims))
lamMean_perSlice_ALL_SIMS<-vector("list",length(sims))

z=1
#for(z in 1:length(sims)){

#whichTrees<-whichTreesAll[[z]]
fewerTrees<-ntrees-length(whichTrees)
nslices=14 # analyze all 14 slices, 5 Ma to 70 Ma 
corr="PAGEL"

exp_N_SIM<-vector("list",length=nslices)

for(j in 1:nslices){
sliceN<-allCladeSetNames[[j]]

	N<-vector("list",length(whichTrees))
	
	for(i in whichTrees){
	#for(i in 1:ntrees){
	uniPGLS<-read.table(file=paste("MamPhy_SIMS_",sims[z],"_",bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-Richness_perSlice_",sliceN,"_SCALED_uniVar_Correct-rootwardWithSingle_9preds.txt",sep=""), header=TRUE)
	tree<-i
	slice<-sliceTimes[[j]]
	num<-j

	N[[i]]<-cbind(uniPGLS[1,],slice,num,tree)
	}

exp_N_SIM[[j]]<-do.call(rbind,N)
}

respVar<-"logRichness"

# calculate 95% CIs (and significance) for the slopes per tree & slice
# ====================================================================
predictors<-colnames(exp_N_SIM[[1]])[10:18]
nPred<-length(predictors)

RESP<-exp_N_SIM

dataToPlot_allPreds_allSlices<-vector("list",length=nslices)

for(j in 1:nslices){ # for each slice from 5 to 70 Ma...

slice<-sliceTimes[[j]]
uniDat<-RESP[[j]]

# SEs to calc 95 CIs...
slopesToPlot<-uniDat[,predictors]
pvalsToPlot<-uniDat[,c(paste("p_",(1:nPred),sep=""))]
SEsToPlot<-uniDat[,c(paste("SE_",(1:nPred),sep=""))]
lamsToPlot<-uniDat[,c(paste("lam_",(1:nPred),sep=""))]

	# The 95% CI for each estimate of each var, plotted together
	dataToPlot_allPreds_perSlice<-vector("list",length(predictors))

	for (i in 1:(ntrees-fewerTrees)){ # so for a given tree IN a slice... taking each of *nPred* vars and calc 95%CI
	#for (i in 1:ntrees){ # so for a given tree IN a slice... taking each of *nPred* vars and calc 95%CI
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
dataToPlot_allPreds_allSlices_allRespVars_SIMS[[z]]<-do.call(rbind,dataToPlot_allPreds_allSlices)


# count the NUMBER of significant runs per time slice per variable across 100 trees
# ==================================================================================
dat<-dataToPlot_allPreds_allSlices_allRespVars_SIMS[[z]]
		
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
	rownames(RES1)<-rep(respVar,nslices)
	
	RES2<-do.call(rbind,lamMean_perSlice)
	colnames(RES2)<-c(predictors,"time","slice")
	rownames(RES2)<-rep(respVar,nslices)

numSignif_perSlice_ALL_SIMS[[z]]<-RES1
lamMean_perSlice_ALL_SIMS[[z]]<-RES2


#} # end 3 sim loop



# UNIVARATE, explainingRichness
# **BOTH-- EMPIRICAL + SIMULATIONS -- PLOTTING**
# PLOT the 95% CIs per slice for each respVar by predictor comparison across 100 trees
# ====================================================================================
corr="PAGEL"
sliceN="5to70Ma"
#predictorsORIG<-colnames(exp_N[[1]])[11:20]
predictors<-c("MRCA", "DR_harm", "DR_skew", "DR_cv", "percentSamp", "PB_Div","BD_Div", "BD_Turn", "BD_Lam", "BD_Mu")
nPred<-length(predictors)
respVar<-"logRichness"
nResp<-length(respVar)
nTotal<-nPred*nResp
sliceTimes<-seq(-5,-70,-5)[1:14]

vertLwd<-0.5
lwdCI<-1.5
horizLine<-0
empirPointCol<-grey(0.3,alpha=0.5)
redCol<-rgb(1,0,0,alpha=0.3)
yLims1<-c(-1.5,2.5)

simAlpha<-0.05
Col<-"darkgoldenrod2"
Col1<-rgb((col2rgb(Col)/255)[1],(col2rgb(Col)/255)[2],(col2rgb(Col)/255)[3],alpha=simAlpha)
Col<-"deepskyblue2"
Col2<-rgb((col2rgb(Col)/255)[1],(col2rgb(Col)/255)[2],(col2rgb(Col)/255)[3],alpha=simAlpha)
Col<-"darkorchid3"
Col3<-rgb((col2rgb(Col)/255)[1],(col2rgb(Col)/255)[2],(col2rgb(Col)/255)[3],alpha=simAlpha)
ColorsSIM<-c(Col1,Col2,Col3)
ColorsSIM_full<-c("darkgoldenrod2","deepskyblue2","darkorchid3")


#pdf(file=paste("cladeLevel",bbone,"_PGLSuni_",corr,"_explaining-Richness_timeSlices_",sliceN,"_95pCI_rootwardWithSingle_10vars.pdf",sep=""),onefile=TRUE, width=(8*nResp),height=(2*nPred))
pdf(file=paste("cladeLevel",bbone,"_PGLSuni_",corr,"_explaining-Richness_timeSlices_",sliceN,"_95pCI_rootwardWithSingle_10vars_WITH-1SIM.pdf",sep=""),onefile=TRUE, width=(8*nResp),height=(2*nPred))

#quartz(width=(8*nResp),height=(2*nPred))
layout(matrix(c(1:nTotal), nPred/2, nResp*2, byrow = TRUE), widths=rep(4,nTotal), heights=rep(3,nTotal))
par(oma = c(5,4,5,3) + 0.1, mar = c(4,1,1,1) + 0.1)



uniDat_ALL<-dataToPlot_allPreds_allSlices_allRespVars
signifDat_ALL<-numSignif_perSlice_ALL

	for(k in 1:length(predictors)){
	uniDat<-uniDat_ALL[which(uniDat_ALL$xVal==predictors[k]),]

	# plot base (no points)
	plot(formula(mean ~ time), data=uniDat, ylim=yLims1, ylab="",xlab="", yaxt="n",xaxt="n",type="n")
	axis(side=2,at=NULL,labels=TRUE)
	axis(side=1,at=NULL,labels=TRUE)

	# add guide lines
	for(j in 1:nslices){
		abline(v=-5*j,lty=2, lwd=vertLwd, col=grey(0.6, alpha=0.5))
	}
		abline(h=horizLine,lty=1, lwd=1, col=grey(0.6, alpha=0.5))

	# plot data
	dat_UNsignif<-uniDat[which(uniDat$color=="red"),]
	if(length(dat_UNsignif[,1])>0){
	plotCI(x=dat_UNsignif$time, add=TRUE,y=dat_UNsignif$mean,ui=dat_UNsignif$high,li=dat_UNsignif$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=redCol,scol=redCol,pch=1,font.lab=2,cex.axis=1.1,cex.lab=1.1)
	}
	dat_signif<-uniDat[which(uniDat$color=="grey"),]
	plotCI(x=dat_signif$time, add=TRUE,y=dat_signif$mean,ui=dat_signif$high,li=dat_signif$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=empirPointCol,scol=empirPointCol,pch=1,font.lab=2,cex.axis=1.1,cex.lab=1.1)

	mtext(side=3,text=paste(respVar," ~ "),font=2,adj=0,line=1)
	mtext(side=3,text=predictors[k],font=2,adj=0)
	
	# add numSignif data
	sigDat<-signifDat_ALL[,c(predictors[k],"time")]
	text(x=sigDat[,2], y = yLims1[2], cex=0.9,font=2, labels = sigDat[,1], col="dark grey")

	} # end *nPred* loop

dev.off()








pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_",corr,"_explaining-Richness_timeSlices_",sliceN,"_95pCI_rootwardWithSingle_3vars_WITH-3-SIMS.pdf",sep=""),onefile=TRUE, width=(4*nResp),height=(4*nPred))
#pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_",corr,"_explaining-Richness_timeSlices_",sliceN,"_95pCI_rootwardWithSingle_3vars.pdf",sep=""),onefile=TRUE, width=(4*nResp),height=(4*nPred))
#pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_",corr,"_explaining-Richness_timeSlices_",sliceN,"_95pCI_rootwardWithSingle_4varsPerSamp.pdf",sep=""),onefile=TRUE, width=(4*nResp),height=(4*nPred))

#quartz(width=(8*nResp),height=(2*nPred))
layout(matrix(c(1:nTotal), nPred, nResp, byrow = TRUE), widths=rep(4,nTotal), heights=rep(3,nTotal))
par(oma = c(5,4,5,3) + 0.1, mar = c(4,1,1,1) + 0.1)


for(k in 1:length(predictors)){
	empiricalDat<-dataToPlot_allPreds_allSlices_allRespVars[which(dataToPlot_allPreds_allSlices_allRespVars$xVal==predictors[k]),]

	# plot base (no points)
	plot(formula(mean ~ time), data=empiricalDat, ylim=yLims1, ylab="",xlab="", yaxt="n",xaxt="n",type="n")
	axis(side=2,at=NULL,labels=TRUE)
	axis(side=1,at=NULL,labels=TRUE)
	if(k==1){ mtext(side=3, line=2, text="MULTIVARIATE")}

	# add guide lines
	for(j in 1:nslices){
		abline(v=-5*j,lty=2, lwd=vertLwd, col=grey(0.6, alpha=0.5))
	}
		abline(h=horizLine,lty=1, lwd=1, col=grey(0.6, alpha=0.5))

# plot EMPIRICAL data
	dat_UNsignif<-empiricalDat[which(empiricalDat$color=="red"),]
	if(length(dat_UNsignif[,1])>0){
	plotCI(x=dat_UNsignif$time, add=TRUE,y=dat_UNsignif$mean,ui=dat_UNsignif$high,li=dat_UNsignif$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=redCol,scol=redCol,pch=1,font.lab=2,cex.axis=1.1,cex.lab=1.1)
	}
	dat_signif<-empiricalDat[which(empiricalDat$color=="grey"),]
	plotCI(x=dat_signif$time, add=TRUE,y=dat_signif$mean,ui=dat_signif$high,li=dat_signif$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=empirPointCol,scol=empirPointCol,pch=1,font.lab=2,cex.axis=1.1,cex.lab=1.1)

	mtext(side=3,text=paste(respVar," ~ "),font=2,adj=0,line=1)
	mtext(side=3,text=predictors[k],font=2,adj=0)
	
	# add numSignif data
	sigDat<-numSignif_perSlice_ALL[,c(predictors[k],"time")]
	text(x=sigDat[,2], y = yLims1[2], cex=0.9,font=2, labels = sigDat[,1], col="dark grey")

	# add LAMBDA data
	lamDat<-lamMean_perSlice_ALL[,c(predictors[k],"time")]
	text(x=lamDat[,2], y = yLims1[1], cex=0.6,font=2, labels = lamDat[,1], col="dark grey")


# plot SIMULATED data
	for(z in 1:length(sims)){

	simDat<-dataToPlot_allPreds_allSlices_allRespVars_SIMS[[z]]
	simColor<-ColorsSIM[[z]]
	simulatedDat<-simDat[which(simDat$xVal==predictors[k]),]

	dat_UNsignif<-simulatedDat[which(simulatedDat$color=="red"),]
	if(length(dat_UNsignif[,1])>0){
	plotCI(x=dat_UNsignif$time, add=TRUE,y=dat_UNsignif$mean,ui=dat_UNsignif$high,li=dat_UNsignif$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=redCol,scol=redCol,pch=1,font.lab=2,cex.axis=1.1,cex.lab=1.1)
	}
	dat_signif<-simulatedDat[which(simulatedDat$color=="grey"),]
	plotCI(x=dat_signif$time, add=TRUE,y=dat_signif$mean,ui=dat_signif$high,li=dat_signif$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=simColor,scol=simColor,pch=1,font.lab=2,cex.axis=1.1,cex.lab=1.1)

	} # end 3 sim loop


} # end *nPred* loop


dev.off()



# ==============
# Visualize 
#======================================
# MULTIVARIATE-- adding in the SIMULATIONS too.
#======================================
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
dirname="/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/cladeLevel_SLICES_explainingRichness"
setwd(dirname)


# **EMPIRICAL**
# ~~~~~~~~~~~~~~~~
# load back in the results
# =========================
ntrees=100
	# for MULTI - 3 var-- 97 trees
	whichTrees<-c(100, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 1, 21, 22, 23, 24, 25, 27, 28, 29, 2, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 3, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 4, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 5, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 6, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 7, 80, 81, 82, 83, 84, 85, 86, 88, 89, 8, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 9)	
	# for MULTI - 4 var (3 + perSamp)-- 98 trees
#	whichTrees<-c(100, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 1, 20, 21, 22, 23, 24, 25, 27, 28, 29, 2, 30, 31, 32, 33, 34, 35, 37, 38, 39, 3, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 4, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 5, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 6, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 7, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 8, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 9)
	
fewerTrees<-ntrees-length(whichTrees)
nslices=14 # analyze all 14 slices, 5 Ma to 70 Ma 
corr="PAGEL"

exp_N<-vector("list",length=nslices)

for(j in 1:nslices){
sliceN<-allCladeSetNames[[j]]

	N<-vector("list",length(whichTrees))
	
	for(i in whichTrees){
	#for(i in 1:ntrees){
	multiPGLS<-read.table(file=paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-Richness_perSlice_",sliceN,"_SCALED_multiVar_Correct-rootwardWithSingle_age-DRmean-DRskew.txt",sep=""), header=TRUE)
#	multiPGLS<-read.table(file=paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-Richness_perSlice_",sliceN,"_SCALED_multiVar_Correct-rootwardWithSingle_age-DRmean-DRskew-perSamp.txt",sep=""), header=TRUE)
	tree<-i
	slice<-sliceTimes[[j]]
	num<-j

	N[[i]]<-cbind(multiPGLS[1,],slice,num,tree)
	}

exp_N[[j]]<-do.call(rbind,N)
}

respVar<-"logRichness"


# calculate 95% CIs (and significance) for the slopes per tree & slice
# ====================================================================
predictors<-colnames(exp_N[[1]])[2:4]
#predictors<-colnames(exp_N[[1]])[2:5]
nPred<-length(predictors)

RESP<-exp_N

dataToPlot_allPreds_allSlices<-vector("list",length=nslices)

for(j in 1:nslices){ # for each slice from 5 to 50 Ma...

slice<-sliceTimes[[j]]
uniDat<-RESP[[j]]

# SEs to calc 95 CIs...
slopesToPlot<-uniDat[,predictors]
pvalsToPlot<-uniDat[,c(paste("Pval",(2:(nPred+1)),sep=""))]
SEsToPlot<-uniDat[,c(paste("SE",(2:(nPred+1)),sep=""))]
lamsToPlot<-uniDat[,"lam"]

	# The 95% CI for each estimate of each var, plotted together
	dataToPlot_allPreds_perSlice<-vector("list",length(predictors))

	for (i in 1:(ntrees-fewerTrees)){ # so for a given tree IN a slice... taking each of *nPred* vars and calc 95%CI
	#for (i in 1:ntrees){ # so for a given tree IN a slice... taking each of *nPred* vars and calc 95%CI
		slopes_i<-slopesToPlot[i,]
		SEs_i<-SEsToPlot[i,]
		pVals_i<-pvalsToPlot[i,]
		lams_i<-lamsToPlot[i]
		dataToPlot_k<-vector("list",length(predictors))
		for (k in 1:length(predictors)){
			color<-if(pVals_i[,k] <= 0.05){ "grey" } else { "red" } 
			lowHigh95_k<-cbind.data.frame(slopes_i[,k]-(1.96*SEs_i[,k]),slopes_i[,k],slopes_i[,k]+(1.96*SEs_i[,k]),color,predictors[k],lams_i,sliceTimes[j],j,i)
			colnames(lowHigh95_k)<-c("low","mean","high","color","xVal","lam","time","slice","tree")
			dataToPlot_k[[k]]<-lowHigh95_k
		} # end *nPred* loop
		dataToPlot_allPreds_perSlice[[i]]<-do.call(rbind,dataToPlot_k)
	} # end 100 tree loop
	dataToPlot_allPreds_allSlices[[j]]<-do.call(rbind,dataToPlot_allPreds_perSlice)
} # end 10 slice loop
dataToPlot_allPreds_allSlices_allRespVars<-do.call(rbind,dataToPlot_allPreds_allSlices)


# count the NUMBER of significant runs per time slice per variable across 100 trees
# ==================================================================================
dat<-dataToPlot_allPreds_allSlices_allRespVars
		
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
	rownames(RES1)<-rep(respVar,nslices)
	
	RES2<-do.call(rbind,lamMean_perSlice)
	colnames(RES2)<-c(predictors,"time","slice")
	rownames(RES2)<-rep(respVar,nslices)

numSignif_perSlice_ALL<-RES1
lamMean_perSlice_ALL<-RES2



# **SIMULATIONS**
# ~~~~~~~~~~~~~~~~
# load back in the results
# =========================
# set NEW dir
dirname="/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/cladeLevel_SLICES_explainingRichness/_res_SIMS_explainingRichness"
setwd(dirname)

ntrees=100
	# for MULTI - mamPhyE -- 97 trees
	whichTrees_1<-c(100,10,11,12,13,15,16,17,18,19,1,20,21,22,23,24,25,26,27,28,29,2,30,31,32,33,34,35,36,37,38,39,3,40,41,42,43,44,45,46,47,48,49,4,50,51,52,53,54,55,56,57,58,59,5,60,61,62,63,64,65,66,67,69,6,70,71,72,73,74,75,76,77,78,79,7,80,81,82,83,84,85,86,87,88,89,8,90,92,93,94,95,96,97,98,99,9)
	# for MULTI - lowE_0p2 -- 97 trees
	whichTrees_2<-c(100,10,11,12,13,15,16,17,18,19,1,20,21,22,23,24,25,26,27,28,2,30,31,32,33,34,35,36,38,39,3,40,41,42,43,44,45,46,47,48,49,4,51,52,53,54,55,56,57,58,59,5,60,61,63,64,65,66,67,69,6,70,71,72,73,74,75,76,77,78,79,7,80,81,82,83,84,85,86,87,88,89,8,90,92,93,94,95,96,97,98,99,9)
	# for MULTI - highE_0p8 -- 97 trees
	whichTrees_3<-c(100,10,11,12,13,15,16,17,18,19,1,20,22,23,24,25,26,27,28,30,31,32,33,34,35,36,38,39,3,40,41,42,43,44,45,46,48,49,4,51,52,53,54,55,57,58,59,5,60,61,63,64,65,66,67,69,6,70,71,72,73,74,75,77,78,79,7,80,81,82,83,84,85,86,87,88,89,92,93,94,95,96,97,98,99,9)

whichTreesAll<-list(whichTrees_1,whichTrees_2,whichTrees_3)

sims<-c("mamPhyE","lowE_0p2","highE_0p8")

dataToPlot_allPreds_allSlices_allRespVars_SIMS<-vector("list",length(sims))
numSignif_perSlice_ALL_SIMS<-vector("list",length(sims))

for(z in 1:length(sims)){

whichTrees<-whichTreesAll[[z]]
fewerTrees<-ntrees-length(whichTrees)
nslices=14 # analyze all 14 slices, 5 Ma to 70 Ma 
corr="PAGEL"

exp_N_SIM<-vector("list",length=nslices)

for(j in 1:nslices){
sliceN<-allCladeSetNames[[j]]

	N<-vector("list",length(whichTrees))
	
	for(i in whichTrees){
	#for(i in 1:ntrees){
	multiPGLS<-read.table(file=paste("MamPhy_SIMS_",sims[z],"_",bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-Richness_perSlice_",sliceN,"_SCALED_multiVar_Correct-rootwardWithSingle_age-DRmean-DRskew.txt",sep=""), header=TRUE)
	tree<-i
	slice<-sliceTimes[[j]]
	num<-j

	N[[i]]<-cbind(multiPGLS[1,],slice,num,tree)
	}

exp_N_SIM[[j]]<-do.call(rbind,N)
}

respVar<-"logRichness"

# calculate 95% CIs (and significance) for the slopes per tree & slice
# ====================================================================
predictors<-colnames(exp_N_SIM[[1]])[2:4]
nPred<-length(predictors)

RESP<-exp_N_SIM

dataToPlot_allPreds_allSlices<-vector("list",length=nslices)

for(j in 1:nslices){ # for each slice from 5 to 50 Ma...

slice<-sliceTimes[[j]]
uniDat<-RESP[[j]]

# SEs to calc 95 CIs...
slopesToPlot<-uniDat[,predictors]
pvalsToPlot<-uniDat[,c(paste("Pval",(2:(nPred+1)),sep=""))]
SEsToPlot<-uniDat[,c(paste("SE",(2:(nPred+1)),sep=""))]

	# The 95% CI for each estimate of each var, plotted together
	dataToPlot_allPreds_perSlice<-vector("list",length(predictors))

	for (i in 1:(ntrees-fewerTrees)){ # so for a given tree IN a slice... taking each of *nPred* vars and calc 95%CI
	#for (i in 1:ntrees){ # so for a given tree IN a slice... taking each of *nPred* vars and calc 95%CI
		slopes_i<-slopesToPlot[i,]
		SEs_i<-SEsToPlot[i,]
		pVals_i<-pvalsToPlot[i,]
		dataToPlot_k<-vector("list",length(predictors))
		for (k in 1:length(predictors)){
			color<-if(pVals_i[,k] <= 0.05){ "grey" } else { "red" } 
			lowHigh95_k<-cbind.data.frame(slopes_i[,k]-(1.96*SEs_i[,k]),slopes_i[,k],slopes_i[,k]+(1.96*SEs_i[,k]),color,predictors[k],sliceTimes[j],j,i)
			colnames(lowHigh95_k)<-c("low","mean","high","color","xVal","time","slice","tree")
			dataToPlot_k[[k]]<-lowHigh95_k
		} # end *nPred* loop
		dataToPlot_allPreds_perSlice[[i]]<-do.call(rbind,dataToPlot_k)
	} # end 100 tree loop
	dataToPlot_allPreds_allSlices[[j]]<-do.call(rbind,dataToPlot_allPreds_perSlice)
} # end 10 slice loop
dataToPlot_allPreds_allSlices_allRespVars_SIMS[[z]]<-do.call(rbind,dataToPlot_allPreds_allSlices)


# count the NUMBER of significant runs per time slice per variable across 100 trees
# ==================================================================================
dat<-dataToPlot_allPreds_allSlices_allRespVars_SIMS[[z]]
	
	numSignif_perSlice<-vector("list",length=nslices)
	for(j in 1:nslices){
	slice<-dat[which(dat$slice==j),]

		numSignif<-c()
		for (k in 1:length(predictors)){
		predPerSlice<-slice[which(slice$xVal==predictors[k]),]
		num<-length(predPerSlice[which(predPerSlice$color=="grey"),][,1])
		numSignif[k]<-round((num/(ntrees-fewerTrees))*ntrees,0)
		}
		numSignif_perSlice[[j]]<-c(numSignif,sliceTimes[j],j)
	}
	RES<-do.call(rbind,numSignif_perSlice)
	colnames(RES)<-c(predictors,"time","slice")
	rownames(RES)<-rep(respVar,nslices)
numSignif_perSlice_ALL_SIMS[[z]]<-RES

} # end 3 sim loop



# MULTIVARIATE
# **BOTH-- EMPIRICAL + SIMULATIONS -- PLOTTING**
# PLOT the 95% CIs per slice for each respVar by predictor comparison across 100 trees
# ====================================================================================
corr="PAGEL"
sliceN="5to70Ma"
#predictorsORIG<-colnames(exp_N[[1]])[11:20]
predictors<-c("MRCA", "DR_harm", "DR_skew")
#predictors<-c("MRCA", "DR_harm", "DR_skew", "percentSamp")
nPred<-length(predictors)
respVar<-"logRichness"
nResp<-length(respVar)
nTotal<-nPred*nResp
sliceTimes<-seq(-5,-70,-5)[1:14]

vertLwd<-0.5
lwdCI<-1.5
horizLine<-0
empirPointCol<-grey(0.3,alpha=0.5)
redCol<-rgb(1,0,0,alpha=0.3)
yLims1<-c(-1.5,2.5)

simAlpha<-0.05
Col<-"darkgoldenrod2"
Col1<-rgb((col2rgb(Col)/255)[1],(col2rgb(Col)/255)[2],(col2rgb(Col)/255)[3],alpha=simAlpha)
Col<-"deepskyblue2"
Col2<-rgb((col2rgb(Col)/255)[1],(col2rgb(Col)/255)[2],(col2rgb(Col)/255)[3],alpha=simAlpha)
Col<-"darkorchid3"
Col3<-rgb((col2rgb(Col)/255)[1],(col2rgb(Col)/255)[2],(col2rgb(Col)/255)[3],alpha=simAlpha)
ColorsSIM<-c(Col1,Col2,Col3)
ColorsSIM_full<-c("darkgoldenrod2","deepskyblue2","darkorchid3")


pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_",corr,"_explaining-Richness_timeSlices_",sliceN,"_95pCI_rootwardWithSingle_3vars_WITH-3-SIMS.pdf",sep=""),onefile=TRUE, width=(4*nResp),height=(4*nPred))
#pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_",corr,"_explaining-Richness_timeSlices_",sliceN,"_95pCI_rootwardWithSingle_3vars.pdf",sep=""),onefile=TRUE, width=(4*nResp),height=(4*nPred))
#pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_",corr,"_explaining-Richness_timeSlices_",sliceN,"_95pCI_rootwardWithSingle_4varsPerSamp.pdf",sep=""),onefile=TRUE, width=(4*nResp),height=(4*nPred))

#quartz(width=(8*nResp),height=(2*nPred))
layout(matrix(c(1:nTotal), nPred, nResp, byrow = TRUE), widths=rep(4,nTotal), heights=rep(3,nTotal))
par(oma = c(5,4,5,3) + 0.1, mar = c(4,1,1,1) + 0.1)


for(k in 1:length(predictors)){
	empiricalDat<-dataToPlot_allPreds_allSlices_allRespVars[which(dataToPlot_allPreds_allSlices_allRespVars$xVal==predictors[k]),]

	# plot base (no points)
	plot(formula(mean ~ time), data=empiricalDat, ylim=yLims1, ylab="",xlab="", yaxt="n",xaxt="n",type="n")
	axis(side=2,at=NULL,labels=TRUE)
	axis(side=1,at=NULL,labels=TRUE)
	if(k==1){ mtext(side=3, line=2, text="MULTIVARIATE")}

	# add guide lines
	for(j in 1:nslices){
		abline(v=-5*j,lty=2, lwd=vertLwd, col=grey(0.6, alpha=0.5))
	}
		abline(h=horizLine,lty=1, lwd=1, col=grey(0.6, alpha=0.5))

# plot EMPIRICAL data
	dat_UNsignif<-empiricalDat[which(empiricalDat$color=="red"),]
	if(length(dat_UNsignif[,1])>0){
	plotCI(x=dat_UNsignif$time, add=TRUE,y=dat_UNsignif$mean,ui=dat_UNsignif$high,li=dat_UNsignif$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=redCol,scol=redCol,pch=1,font.lab=2,cex.axis=1.1,cex.lab=1.1)
	}
	dat_signif<-empiricalDat[which(empiricalDat$color=="grey"),]
	plotCI(x=dat_signif$time, add=TRUE,y=dat_signif$mean,ui=dat_signif$high,li=dat_signif$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=empirPointCol,scol=empirPointCol,pch=1,font.lab=2,cex.axis=1.1,cex.lab=1.1)

	mtext(side=3,text=paste(respVar," ~ "),font=2,adj=0,line=1)
	mtext(side=3,text=predictors[k],font=2,adj=0)
	
	# add numSignif data
	sigDat<-numSignif_perSlice_ALL[,c(predictors[k],"time")]
	text(x=sigDat[,2], y = yLims1[2], cex=0.9,font=2, labels = sigDat[,1], col="dark grey")

	# add LAMBDA data
	lamDat<-lamMean_perSlice_ALL[,c(predictors[k],"time")]
	text(x=lamDat[,2], y = yLims1[1], cex=0.6,font=2, labels = lamDat[,1], col="dark grey")


# plot SIMULATED data
	for(z in 1:length(sims)){

	simDat<-dataToPlot_allPreds_allSlices_allRespVars_SIMS[[z]]
	simColor<-ColorsSIM[[z]]
	simulatedDat<-simDat[which(simDat$xVal==predictors[k]),]

	dat_UNsignif<-simulatedDat[which(simulatedDat$color=="red"),]
	if(length(dat_UNsignif[,1])>0){
	plotCI(x=dat_UNsignif$time, add=TRUE,y=dat_UNsignif$mean,ui=dat_UNsignif$high,li=dat_UNsignif$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=redCol,scol=redCol,pch=1,font.lab=2,cex.axis=1.1,cex.lab=1.1)
	}
	dat_signif<-simulatedDat[which(simulatedDat$color=="grey"),]
	plotCI(x=dat_signif$time, add=TRUE,y=dat_signif$mean,ui=dat_signif$high,li=dat_signif$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=simColor,scol=simColor,pch=1,font.lab=2,cex.axis=1.1,cex.lab=1.1)

	} # end 3 sim loop


} # end *nPred* loop


dev.off()

















#===============
# Make SIMULATIONS of the full MamPhy for testing these null expectations.
#===============
library(moments); library(nlme); library(ape); library(picante); library(phytools); library(geiger)

# First I want to get the ML birth and death rates for the whole mamPhy.
setwd(dirname)
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


#=======================================
# For SIMS-- Make time slices every 5 Ma to create clades-- each of the three SIMS sets.
#=======================================
# Make the time slices at 5 Ma 
###
library(moments); library(nlme); library(ape); library(picante); library(phytools); library(geiger)
library(foreach);library(doSNOW)

cl = makeCluster(100, type = 'MPI', outfile="")
registerDoSNOW(cl)

ntrees<-100
bbone<-"NDexp"
#for(i in 1:ntrees){
foreach(i=1:ntrees, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

sims<-c("mamPhyE","lowE_0p2","highE_0p8")

lengthsALL_SIMS<-vector("list",length(sims))
sumDat_SIMS<-vector("list",length(sims))
for(z in 1:length(sims)){

## read in 1 of 100 sumulated full trees
simPhy<-read.tree(paste("MamPhy_SIMS_",sims[z],"_",bbone,"_tree",i,".tre",sep=""))

# set up the slices
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery

root=max(node.age(simPhy)$ages)

allCladeSets<-vector("list",length=numSlices)
slicePhys<-vector("list",length=numSlices)
allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
	# make tipward slice clades
	allCladeSets[[j]]<-treeSlice(simPhy, slice=root-(sliceEvery*j), trivial=TRUE) # now keeping trivial (single tips)
		# write clades to per-slice files
		write.tree(allCladeSets[[j]],file=paste("MamPhy_SIMS_",sims[z],"_",bbone,"_sample100_",i,"_timeslice_cladeTrees_",(sliceEvery*j),"Ma.trees",sep=""))
	# make rootward slice *backbones*
	slicePhys[[j]]<-treeSlice(simPhy, slice=root-(sliceEvery*j), trivial=TRUE, orientation="rootwards") # toward root, slicePhys
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
			node <- getMRCA(simPhy, cladeSp) # find the MRCA node of those species
			slicePhys[[j]]$tip.label[which(slicePhys[[j]]$tip.label==node)]<-newTipNames[k]
		}
	}
}

# write slice phys
for(j in 1:length(slicePhys)){
	write.tree(slicePhys[[j]], file=paste("MamPhy_SIMS_",sims[z],"_",bbone,"_tree_",i,"_CORRECT_slicePhy-5to70Ma_by5.trees",sep=""), append=TRUE)
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
save(allLengths_i, file=paste("MamPhy_SIMS_",sims[z],"_",bbone,"_sample100_",i,"_timeslice_cladeRichnesses_wSingletons.Rda",sep=""))


} # end 3 sim loop

} # end 100 tree loop



#===============
# For SIMULATIONS--  Calculate per-slice, per-clade summary values 
#===============
#######
# ALL trees simulated and WRITTEN
# want to LOAD back in, re-do the calculations
#######

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




#===============
# For SIMULATIONS-- Now in parallel, and with doing the PGLS also...
#===============
########
setwd(dirname)
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

# Z-score STANDARDIZED normal MRCA
results<-vector("list",length(allCladeSetNames))
for (j in 1: length(allCladeSetNames)){
	res<-read.table(paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS.txt",sep=""), header=TRUE)
	vars<-res[,1:11]
	varsScale<-scale(vars, center=TRUE,scale=TRUE)
	results[[j]]<-cbind(res$richness,varsScale)
	colnames(results[[j]])<-c("richnessN",colnames(res[,1:11]))
}
 

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
	partSlopesPGLS_All[j,k]<-round(sum$tTable[k],digits=6)
	}
	#partSlopesPGLS_All[j,9]<-round(int95s[[1]][1],digits=6)
	partSlopesPGLS_All[j,10]<-round(sum$modelStruct[[1]][[1]],digits=6)
	#partSlopesPGLS_All[j,11]<-round(int95s[[1]][3],digits=6)
	partSlopesPGLS_All[j,12]<-round(sum$AIC,digits=0)
	partSlopesPGLS_All[j,13]<-round(sum$tTable[13], digits=6)
	partSlopesPGLS_All[j,14]<-round(sum$tTable[14], digits=6)
	partSlopesPGLS_All[j,15]<-round(sum$tTable[15], digits=6)
	partSlopesPGLS_All[j,16]<-round(sum$tTable[16], digits=6)

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
		
		uniPGLS_allInts[j,k]<-round(sum1$tTable[1], digits=6)
		uniPGLS_allSlopes[j,k]<-round(sum1$tTable[2], digits=6)
		uniPGLS_allSEs[j,k]<-round(sum1$tTable[4], digits=6)
		uniPGLS_allPs[j,k]<-round(sum1$tTable[8], digits=6)
		uniPGLS_allLams[j,k]<-round(sum1$modelStruct[[1]][[1]], digits=6)
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

	partSlopesPGLS_All[j,1]<-round(sum$coef[[1]],digits=6)
	partSlopesPGLS_All[j,2]<-round(sum$coef[[2]],digits=6)
	partSlopesPGLS_All[j,3]<-round(sum$coef[[3]],digits=6)
	partSlopesPGLS_All[j,4]<-round(sum$coef[[4]],digits=6)
	partSlopesPGLS_All[j,5]<-round(sum$AIC,digits=0)
	partSlopesPGLS_All[j,6]<-round(sum$tTable[14], digits=6)
	partSlopesPGLS_All[j,7]<-round(sum$tTable[15], digits=6)
	partSlopesPGLS_All[j,8]<-round(sum$tTable[16], digits=6)

	form<-(log(richness) ~ MRCA + DR_harm)
	fit2<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
	sum<-summary(fit2)

	partSlopesPGLS_12[j,1]<-round(sum$coef[[1]],digits=6)
	partSlopesPGLS_12[j,2]<-round(sum$coef[[2]],digits=6)
	partSlopesPGLS_12[j,3]<-round(sum$coef[[3]],digits=6)
	partSlopesPGLS_12[j,4]<-round(sum$AIC,digits=0)
	partSlopesPGLS_12[j,5]<-round(sum$tTable[11], digits=6)
	partSlopesPGLS_12[j,6]<-round(sum$tTable[12], digits=6)

	vars<-c("MRCA","DR_harm","DR_skew","PB_Div","BD_Lam","BD_Mu", "BD_Div", "BD.ms0", "BD.ms0p5", "BD.ms0p9")
	for(k in 1:length(uniPGLS_allSlopes)){
		form<-as.formula(paste("log(richness) ~ ", vars[k],sep=""))
		fit1<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
		sum1<-summary(fit1)
		
		uniPGLS_allInts[j,k]<-round(sum1$coef[[1]],digits=6)
		uniPGLS_allSlopes[j,k]<-round(sum1$coef[[2]],digits=6)
		uniPGLS_allPs[j,k]<-round(sum1$tTable[8], digits=6)
		#uniPGLS_allLams[j,k]<-round(sum1$modelStruct[[1]][[1]], digits=6)
	}	

}
write.table(partSlopesPGLS_All,paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLSmulti_timeSlices_ACTUAL_partSlopes123.txt",sep=""))
write.table(partSlopesPGLS_12,paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLSmulti_timeSlices_ACTUAL_partSlopes12.txt",sep=""))

uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs)#,uniPGLS_allLams)
write.table(uniPGLS,paste("MamPhy_SIMS_",sims[q],"_",bbone,"_sample100_",i,"_PGLS_timeSlices_ACTUAL_uniINTS_uniSLO_uniP.txt",sep=""))

}

}


#===========================
# For SIMULATIONS-- summarize and re-load all
#===========================
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


#===========================
# For ALL-- final multivariate plotting 
#===========================
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



# ==============
#  For ALL-- Also get the UNIVARIATE PENDANTS going...
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



#==================
# for ALL-- Make the actual data SCATTER plots...
#==================

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







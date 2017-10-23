#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy -- Upham et al. 2017
###
# Figure 4 - time slices to explain clade diversification rate
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run PGLS on the timeslice clade summary stats -- aim to ** explain variation in div rate and clade age**
######
# - In parallel, read back in slice phylo backbones that unite the clades at each slice (14 slices at 5 Ma intervals)
# - Load in rate and trait slice summaries for 1 of 100 trees (standardize predictors)
# - Run for loop of 2 DR vars + MRCA to explain variation in them using 7 predictors
# - Compare models with ANOVA and AIC
# - Record results and loop across 14 slices, 100 trees

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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













#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy v1 -- Upham et al. 2017
###
# Figure 3 - time slices to explain clade richness
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run PGLS -- aim to ** explain variation RICHNESS ** using timeslice clade summary stats 
######
# - In parallel, read back in slice phylo backbones (14 slices at 5 Ma intervals)
# - Load in age and rate slice summaries for 1 of 100 trees (standardize predictors)
# - Run loop of 4 response vars 
# - Setup empty dataframes to receive PGLS results 
# - Nest the loop of 26 predictor vars



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
	PRED<-res[,c(1:5,7:15)]
		# STANDARDIZE predictors as Z-SCORES
		predScale<-scale(PRED, center=TRUE,scale=TRUE)

	# Join together RESP and standardized PRED
	results[[q]]<-cbind(RESP,predScale[,1:14])
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

## UNIVARIATE
#	for(k in 1:length(predictors)){
#
#		form<-as.formula(paste("logRichness", " ~ ", predictors[k], sep=""))
##		fit1<-gls(form, data=dat, method="ML")
##		fit1<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
#		for (p in c(0.5,seq(0,1,by=0.01))) {possibleError <- tryCatch(
#		      gls(form, correlation=corPagel(value=p,phy=cladeData$phy), data=dat, method="ML"),
#		      error=function(e) e)
#		if(inherits(possibleError, "gls")) break		
#		if(inherits(possibleError, "error")) next}
#		fit1<-possibleError
#
#		sum1<-summary(fit1)
#		
#		uniPGLS_allInts[j,k]<-round(sum1$tTable[1], digits=6)
#		uniPGLS_allSlopes[j,k]<-round(sum1$tTable[2], digits=6)
#		uniPGLS_allSEs[j,k]<-round(sum1$tTable[4], digits=6)
#		uniPGLS_allPs[j,k]<-round(sum1$tTable[8], digits=6)
#		uniPGLS_allLams[j,k]<-round(sum1$modelStruct[[1]][[1]], digits=6)		
#
#		} # cycles each of *nCols* predictors
#
#		uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs,uniPGLS_allSEs,uniPGLS_allLams)[j,]
#
	#corr<-"NO_TREE"
	#corr<-"BROWNIAN"
	corr<-"PAGEL"
	sliceN<-allCladeSetNames[[j]] #"ALL" #"45Ma-only" #"10to60Ma" # # #"35Ma-only" ##"all" #"60Ma"#"all" #
#
#	write.table(uniPGLS,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-Richness_perSlice_",sliceN,"_SCALED_uniVar_Correct-rootwardWithSingle_14preds.txt",sep=""))

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

## MULTIVARIATE
## 3 var + percentSampling
#	form2<-(logRichness ~ MRCA + DR_harm + DR_skew + percentSamp)
#	#fit2<-gls(form2, data=dat, method="ML")
#	#fit2<-gls(form2, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
#		for (p in c(0.5,seq(0,1,by=0.01))) {possibleError <- tryCatch(
#		      gls(form2, correlation=corPagel(value=p,phy=cladeData$phy), data=dat, method="ML"),
#		      error=function(e) e)
#		if(inherits(possibleError, "gls")) break		
#		if(inherits(possibleError, "error")) next}
#		fit2<-possibleError
#
#		sum2<-summary(fit2)
#
#	for(k in 1:10){
#	multiPGLS_Per[k]<-round(sum2$tTable[k],digits=6)
#	}
#	multiPGLS_Per[11]<-round(sum2$modelStruct[[1]][1],digits=6)
#	multiPGLS_Per[12]<-round(sum2$tTable[16], digits=6)
#	multiPGLS_Per[13]<-round(sum2$tTable[17], digits=6)
#	multiPGLS_Per[14]<-round(sum2$tTable[18], digits=6)
#	multiPGLS_Per[15]<-round(sum2$tTable[19], digits=6)
#	multiPGLS_Per[16]<-round(sum2$tTable[20], digits=6)
#
#	write.table(multiPGLS_Per,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-Richness_perSlice_",sliceN,"_SCALED_multiVar_Correct-rootwardWithSingle_age-DRmean-DRskew-perSamp.txt",sep=""))

	} # cycles each time slice

} # cycles the 100 trees

stopCluster(cl)

q()

n





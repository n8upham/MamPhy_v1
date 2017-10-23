#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy -- Upham et al. 2017
###
# Figure 4 - time slices to explain clade diversification rate
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run PGLS on the timeslice clade summary stats -- aim to ** explain variation in div rate **
######
# - In parallel, read back in slice phylo backbones (14 slices at 5 Ma intervals)
# - Load in rate and trait slice summaries for 1 of 100 trees (standardize predictors)
# - Run for ONE of 4 response vars 
# - Setup empty dataframes to receive PGLS results 
# - Nest the loop of 26 predictor vars

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# ==============
# do the PGLS
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
slicePhys<-read.tree(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_CORRECT_slicePhy-5to70Ma_by5.trees",sep=""))

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
	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[q]],"cladeSTATS_withEcoTraits.txt",sep=""), header=TRUE)
	# select RESPONSE vars
	RESP<-cbind.data.frame(res[,c("DR_harm","BD_Div","richness","DR_skew")])
	# select PREDICTOR vars
	PRED_LOGvars<-colnames(head(res,0)[c(11:13,16,19)])
	PRED_NONvars<-colnames(head(res,0)[c(3,5,15,18,20:36)])
	PRED_all<-cbind(log(res[,PRED_LOGvars]), res[,PRED_NONvars])
		# STANDARDIZE predictors as Z-SCORES
		predScale<-scale(PRED_all, center=TRUE,scale=TRUE)

	# Join together RESP and standardized PRED
	results[[q]]<-cbind(log(RESP[,1:3]),RESP[,4],predScale)
	colnames(results[[q]])<-c("logTipDR","logBD_Div","logRichness","DR_skew",paste("log",PRED_LOGvars,sep=""),PRED_NONvars)
}


# Run for ONE of 4 response vars 
# ===============================
respVar<-colnames(results[[1]][,c(2,3,4,1)]) #"logBD_Div"   "logRichness" "DR_skew"     "logTipDR" 
#z=1 #"logBD_Div"   
#z=2 #"logRichness" 
#z=3 #""DR_skew"    
#z=4 #""logTipDR" 

for(z in 1:length(respVar)){

# Setup empty dataframes to receive PGLS results 
# ==============================================
# UNIVARIATE first
sliceTimes<-seq(-5,-70,-5)

predictorsALL<-c(paste("log",PRED_LOGvars,sep=""),PRED_NONvars)
predictors<-predictorsALL[c(1:7)]
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


# Nest the loop of 24 predictor vars
# ===============================
	for (j in 1:length(results)){
	#j=2 # 10 Ma time slice
	#j=7 # 35 Ma time slice
	#j=9 # 45 Ma time slice
	#j=12 # 60 Ma time slice
	cladeData<-treedata(slicePhys[[j]],na.omit(results[[j]])) 
	dat<-as.data.frame(cladeData$data)

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
		
		uniPGLS_allInts[j,k]<-round(sum1$tTable[1], digits=6)
		uniPGLS_allSlopes[j,k]<-round(sum1$tTable[2], digits=6)
		uniPGLS_allSEs[j,k]<-round(sum1$tTable[4], digits=6)
		uniPGLS_allPs[j,k]<-round(sum1$tTable[8], digits=6)
		uniPGLS_allLams[j,k]<-round(sum1$modelStruct[[1]][[1]], digits=3)		

		} # cycles each of 24 predictors

	} # cycles each time slice

#corr<-"NO_TREE"
#corr<-"BROWNIAN"
corr<-"PAGEL"
sliceN<-"ALL" #"45Ma-only" #"10to60Ma" # # #"35Ma-only" ##"all" #"60Ma"#"all" #

uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs,uniPGLS_allSEs,uniPGLS_allLams)
write.table(uniPGLS,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-",respVar[z],"_timeSlices_",sliceN,"_SCALED_uniVar_Correct_7preds.txt",sep=""))

#uniPGLS[z,]<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs,uniPGLS_allSEs,uniPGLS_allLams)[j,]

} # cycles the 4 respVars

#write.table(uniPGLS,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining8_timeSlices_",sliceN,"_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))

} # cycles the 100 trees



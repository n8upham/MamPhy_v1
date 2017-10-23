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
	PRED_LOGvars<-colnames(head(res,0)[c(11:13,16,19)]) #,"richness")
	PRED_NONvars<-colnames(head(res,0)[c(3,5,15,18,20:36)]) # "DR_harm","DR_skew","BD_Div"
	PRED_all<-cbind(log(res[,PRED_LOGvars]), res[,PRED_NONvars])
		# STANDARDIZE predictors as Z-SCORES
		predScale<-scale(PRED_all, center=TRUE,scale=TRUE)

	# Join together RESP and standardized PRED
	results[[q]]<-cbind(log(RESP[,1:3]),RESP[,4],predScale)
	colnames(results[[q]])<-c("logTipDR","logBD_Div","logRichness","DR_skew",paste("log",PRED_LOGvars,sep=""),PRED_NONvars)
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

predictorsALL<-c(paste("log",PRED_LOGvars,sep=""),PRED_NONvars)
predictors<-c(predictorsALL[c(1:16)]) # focusing on continuous vars, trophic vars
nCols<-length(predictors)

#uniPGLS_allSlopes<-data.frame(matrix(NA, nrow = length(results), ncol = nCols),row.names=sliceTimes)
#colnames(uniPGLS_allSlopes)<-predictors
#uniPGLS_allSEs<-data.frame(matrix(NA, nrow = length(results), ncol = nCols),row.names=sliceTimes)
#colnames(uniPGLS_allSEs)<-c(paste("SE_",1:nCols,sep=""))
#uniPGLS_allInts<-data.frame(matrix(NA, nrow = length(results), ncol = nCols),row.names=sliceTimes)
#colnames(uniPGLS_allInts)<-c(paste("i_",1:nCols,sep=""))
#uniPGLS_allPs<-data.frame(matrix(NA, nrow = length(results), ncol = nCols),row.names=sliceTimes)
#colnames(uniPGLS_allPs)<-c(paste("p_",1:nCols,sep=""))
#uniPGLS_allLams<-data.frame(matrix(NA, nrow = length(results), ncol = nCols),row.names=sliceTimes)
#colnames(uniPGLS_allLams)<-c(paste("lam_",1:nCols,sep=""))



# MULTIVARIATE
sliceTimes<-seq(-5,-70,-5)

Cols<-c("int","CarnOrNot", "HerbOrNot", "SE1","SE2","SE3","Lam","AIC","Pval1","Pval2")
partSlopesPGLS_trophic2<-data.frame(matrix(NA, nrow = length(results), ncol = length(Cols)),row.names=sliceTimes)
colnames(partSlopesPGLS_trophic2)<-Cols

Cols<-c("int","DiurOrNot", "NoctOrNot","SE1","SE2","SE3","Lam","AIC","Pval1","Pval2")
partSlopesPGLS_activity2<-data.frame(matrix(NA, nrow = length(results), ncol = length(Cols)),row.names=sliceTimes)
colnames(partSlopesPGLS_activity2)<-Cols

Cols<-c("int","SubtOrNot", "AquaOrNot", "TerrOrNot", "ArboOrNot", "SE1","SE2","SE3","SE4","SE5","Lam","AIC","Pval1","Pval2","Pval3","Pval4")
partSlopesPGLS_lifemode4<-data.frame(matrix(NA, nrow = length(results), ncol = length(Cols)),row.names=sliceTimes)
colnames(partSlopesPGLS_lifemode4)<-Cols

#for (j in 2:(length(results)-2)){

#j=2 # 10 Ma time slice
#j=7 # 35 Ma time slice
j=9 # 45 Ma time slice
#j=12 # 60 Ma time slice
	#rownames(results[[j]])<-slicePhys[[j]]$tip.label
	cladeData<-treedata(slicePhys[[j]],na.omit(results[[j]])) # <<< 60 Ma slice fails here because of NAs in the SubOrNot variable... only in some trees though... 
	dat<-as.data.frame(cladeData$data)

#dat<-as.data.frame(na.omit(results[[j]]))
# UNIVARIATE
	for(k in 1:length(uniPGLS_allSlopes)){
#		respVar<-"tipDR_mean"
#		respVar<-"BD_Div"
#		respVar<-"BD_Lam"
#		respVar<-"BD_Mu"
		if(k==11) { next } else {
#		form<-as.formula(paste(respVar[z], " ~ ", predictors[k], sep=""))
#		form<-as.formula(paste("log(",respVar[z], ") ~ ", predictors[k], sep=""))
		form<-as.formula(paste(predictors[k], " ~ ", respVar[z],sep=""))
#		fit1<-gls(form, data=dat, method="ML")
#		fit1<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
		for (p in seq(0,1,by=0.01)) {possibleError <- tryCatch(
		      gls(form, correlation=corPagel(value=p,phy=cladeData$phy), data=dat, method="ML"),
		      error=function(e) e)
		if(inherits(possibleError, "gls")) break		
		if(inherits(possibleError, "error")) next}
		fit1<-possibleError

		sum1<-summary(fit1)
		
		uniPGLS_allInts[j,k]<-round(sum1$tTable[1], digits=3)
		uniPGLS_allSlopes[j,k]<-round(sum1$tTable[2], digits=3)
		uniPGLS_allSEs[j,k]<-round(sum1$tTable[4], digits=3)
		uniPGLS_allPs[j,k]<-round(sum1$tTable[8], digits=3)
		uniPGLS_allLams[j,k]<-round(sum1$modelStruct[[1]][[1]], digits=3)
		}
	}	
## MULTIVARIATE
#	multiPred<-"CarnOrNot + HerbOrNot"
#	form<-as.formula(paste(respVar, " ~ ", multiPred,sep=""))
#
##	fit<-gls(form, data=dat, method="ML")
##	fit<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
#		for (p in seq(0,1,by=0.01)) {possibleError <- tryCatch(
#		      gls(form, correlation=corPagel(value=p,phy=cladeData$phy), data=dat, method="ML"),
#		      error=function(e) e)
#		if(inherits(possibleError, "gls")) break		
#		if(inherits(possibleError, "error")) next}
#		fit<-possibleError
#	sum<-summary(fit)
#
#	for(k in 1:6){
#	partSlopesPGLS_trophic2[j,k]<-round(sum$tTable[k],digits=3)
#	}
#	partSlopesPGLS_trophic2[j,7]<-round(sum$modelStruct[[1]][1],digits=3)#int95s[[1]][2],digits=3)
#	partSlopesPGLS_trophic2[j,8]<-round(sum$AIC,digits=0)
#	partSlopesPGLS_trophic2[j,9]<-round(sum$tTable[11], digits=3)
#	partSlopesPGLS_trophic2[j,10]<-round(sum$tTable[12], digits=3)
#
#######
#	multiPred<-"DiurOrNot + NoctOrNot"
#	form<-as.formula(paste(respVar, " ~ ", multiPred,sep=""))
#
##	fit<-gls(form, data=dat, method="ML")
##	fit<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
#		for (p in seq(0,1,by=0.01)) {possibleError <- tryCatch(
#		      gls(form, correlation=corPagel(value=p,phy=cladeData$phy), data=dat, method="ML"),
#		      error=function(e) e)
#		if(inherits(possibleError, "gls")) break		
#		if(inherits(possibleError, "error")) next}
#		fit<-possibleError
#	sum<-summary(fit)
#
#	for(k in 1:6){
#	partSlopesPGLS_activity2[j,k]<-round(sum$tTable[k],digits=3)
#	}
#	partSlopesPGLS_activity2[j,7]<-round(sum$modelStruct[[1]][1],digits=3)#int95s[[1]][2],digits=3)
#	partSlopesPGLS_activity2[j,8]<-round(sum$AIC,digits=0)
#	partSlopesPGLS_activity2[j,9]<-round(sum$tTable[11], digits=3)
#	partSlopesPGLS_activity2[j,10]<-round(sum$tTable[12], digits=3)
#
#######
#	multiPred<-"SubtOrNot + AquaOrNot + TerrOrNot + ArboOrNot"
#	form<-as.formula(paste(respVar, " ~ ", multiPred,sep=""))
#
##	fit<-gls(form, data=dat, method="ML")
##	fit<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=dat, method="ML")
#		for (p in seq(0,1,by=0.01)) {possibleError <- tryCatch(
#		      gls(form, correlation=corPagel(value=p,phy=cladeData$phy), data=dat, method="ML"),
#		      error=function(e) e)
#		if(inherits(possibleError, "gls")) break		
#		if(inherits(possibleError, "error")) next}
#		fit<-possibleError
#	sum<-summary(fit)
#
#	for(k in 1:10){
#	partSlopesPGLS_lifemode4[j,k]<-round(sum$tTable[k],digits=3)
#	}
#	partSlopesPGLS_lifemode4[j,11]<-round(sum$modelStruct[[1]][1],digits=3)#int95s[[1]][2],digits=3)
#	partSlopesPGLS_lifemode4[j,12]<-round(sum$AIC,digits=0)
#	partSlopesPGLS_lifemode4[j,13]<-round(sum$tTable[17], digits=3)
#	partSlopesPGLS_lifemode4[j,14]<-round(sum$tTable[18], digits=3)
#	partSlopesPGLS_lifemode4[j,15]<-round(sum$tTable[19], digits=3)
#	partSlopesPGLS_lifemode4[j,16]<-round(sum$tTable[20], digits=3)



# Run for ONE of 4 response vars 
# ===============================
respVar<-colnames(results[[1]][,c(2,3,4,1)]) #"logBD_Div"   "logRichness" "DR_skew"     "logTipDR" 
#z=1 #"logBD_Div"   
#z=2 #"logRichness" 
#z=3 #""DR_skew"    
#z=4 #""logTipDR" 
uniPGLS<-data.frame(matrix(NA, nrow = length(respVar), ncol = nCols*5),row.names=respVar)
colnames(uniPGLS)<-c(colnames(uniPGLS_allInts),colnames(uniPGLS_allSlopes),colnames(uniPGLS_allPs),colnames(uniPGLS_allSEs),colnames(uniPGLS_allLams))

#for(z in 1:length(respVar)){
z=1 # logBDdiv

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
		
		uniPGLS_allInts[j,k]<-round(sum1$tTable[1], digits=6)
		uniPGLS_allSlopes[j,k]<-round(sum1$tTable[2], digits=6)
		uniPGLS_allSEs[j,k]<-round(sum1$tTable[4], digits=6)
		uniPGLS_allPs[j,k]<-round(sum1$tTable[8], digits=6)
		uniPGLS_allLams[j,k]<-round(sum1$modelStruct[[1]][[1]], digits=3)		

		} # cycles each of *nCols* predictors

		uniPGLS[z,]<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs,uniPGLS_allSEs,uniPGLS_allLams)[j,]

	} # cycles the 4 respVars

	#corr<-"NO_TREE"
	#corr<-"BROWNIAN"
	corr<-"PAGEL"
	sliceN<-allCladeSetNames[[j]] #"ALL" #"45Ma-only" #"10to60Ma" # # #"35Ma-only" ##"all" #"60Ma"#"all" #

	write.table(uniPGLS,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-4varDR_perSlice_",sliceN,"_SCALED_uniVar_Correct_16preds.txt",sep=""))

	} # cycles each time slice

#uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs,uniPGLS_allSEs,uniPGLS_allLams)
#write.table(uniPGLS,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-",respVar[z],"_timeSlices_",sliceN,"_SCALED_uniVar_Correct_7preds.txt",sep=""))

} # cycles the 100 trees



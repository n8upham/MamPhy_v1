#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy -- Upham et al. 2017
###
# Figure 4 - time slices to explain clade diversification rate
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run PGLS on the TIP LEVEL -- aim to ** explain variation in div ratee**
######
# - In parallel, read back in mamPhy
# - Load in rate and trait data for SPECIES -- DR for 1 of 100 trees 
# - Run for loop of tipDR vars to explain variation using 7 predictors
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

# start cluster
cl = makeCluster(32, type = 'MPI')
registerDoSNOW(cl)

#ranTrees<-c(2, 6, 10, 11, 12, 21, 23, 24, 25, 27, 29, 31, 33, 40, 41, 48, 53, 54, 59, 60, 63, 65, 67, 68, 72, 73, 74, 77, 82, 83, 84, 85, 86, 91, 94, 96, 98)
#toRun<-setdiff(1:100,ranTrees)
#foreach(i=toRun, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

ntrees=100
foreach(i=1:ntrees, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools', 'caper'), .combine=cbind, .verbose=TRUE) %dopar% {

#==================
# Load in 1 tree of 100
#==================
mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_nexus.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 
#write.tree(mamPhy,file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""))
#tree1=scan(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""), what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

#==================
# Calculate ES and DR on per-tip basis
#==================
# gives pairwise clade matrix from CAIC function
#clade_matrix = readCAIC(tree1)
#
## calculate and write to file
#DR = 1/ES_v2(clade_matrix)
#ES = ES_v2(clade_matrix)
#res = cbind.data.frame(DR,ES)
#res1 = res[order(rownames(res)),]
#
#write.table(res1, file=paste(bbone,"_sample100_",i,"_DRtips.txt",sep=""))
res1<-read.table(file=paste(bbone,"_sample100_",i,"_DRtips.txt",sep=""), header=TRUE)

# tip DR values to use for this tree's comparisons
tipDR_i<-res1$DR
names(tipDR_i)<-rownames(res1)

#==================
# Load in stats about tip TAXONOMY and GENE SAMPLING
#==================
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
head(cladesDR)
sorted<-cladesDR[order(cladesDR$tiplabel),]

taxonomyVars<-c("tiplabel","gen", "fam", "ord", "higher", "genes")
taxFinal<-sorted[,taxonomyVars]

#==================
# Subset by comparison of ecological traits to use
#==================
tipDataAll<-read.table(file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)
	
PRED_LOGvars<-c("homeRange_km2_ext", "trophic123", "geoArea_km2", "GenerationLength_d")#,"lifemode1234", "activity123") # "BM_final_kg", "DispDistAll_km_ext", 
#PRED_LOGvars<-c("homeRange_km2_ext", "trophic123", "geoArea_km2", "GenerationLength_d","DispDistAll_km_ext", "lifemode1234") # "BM_final_kg", 
PRED_NONvars<-c("Lat_centroid") #, "Lon_centroid")
PRED_names<-c(paste("log",c("HomeRange_km2","Trophic123", "GeoArea_km2", "GenLength_d"),sep=""), PRED_NONvars,"genes") # "BodyMass_kg", "DispDist_km", 
#PRED_names<-c(paste("log",c("HomeRange_km2","Trophic123", "GeoArea_km2", "GenLength_d", "DispersalIndex"),sep=""), PRED_NONvars,"genes") # "BodyMass_kg", "DispDist_km", 
#PRED_names<-c(paste("log",c("HomeRange_km2","BodyMass_kg", "DispDist_km"),sep=""), PRED_NONvars,"genes") #"Trophic123", 

catsSelected<-c("CarnOrNot", "HerbOrNot", "OmniOrNot", "AquaOrNot", "ArboOrNot", "FlysOrNot", "SubtOrNot", "TerrOrNot", "NoctOrNot", "CathOrNot", "DiurOrNot", "MarineOrNot")

allDatFinal<-cbind.data.frame(taxFinal[,c(1:5)],log(tipDR_i),log(tipDataAll[,PRED_LOGvars[1:4]]),tipDataAll[,PRED_NONvars],(taxFinal$genes),tipDataAll[,catsSelected])
#allDatFinal<-cbind.data.frame(taxFinal,log(tipDR_i),log(tipDataAll[,PRED_LOGvars[1:4]]),log(tipDataAll[,PRED_LOGvars[5]]*tipDataAll[,PRED_LOGvars[6]]),tipDataAll[,PRED_NONvars],(taxFinal$genes),tipDataAll[,catsSelected])
colnames(allDatFinal)<-c(taxonomyVars[1:5],"logTipDR",PRED_names,catsSelected)

respVar<-"logTipDR"
# if doing GLOBAL
datFinal<-allDatFinal[,c("tiplabel",respVar,PRED_names)] # just the 9 vars in the MODELS

# setup CAPER object
cladeData<-comparative.data(phy=mamPhy,data=datFinal, names.col="tiplabel",na.omit=TRUE)

# Setup empty dataframes to receive PGLS results 
# ==============================================
# MULTIVARIATE complex-- just subtracting one var each time and looking at DIFF IN AICc scores.

predictors<-PRED_names
nPred<-length(PRED_names)
 
#  # how correlated are the predictor variables?
#  res<-c()
#  	for(j in 1:(nPred-1)){
#  		form<-as.formula(paste(" ~ ", predictors[1], "+", predictors[j+1], sep=""))
#  		test<-cor.test(formula=form, data=cladeData$data, method="spearman", alternative="two.sided")
#  		res[j]<-test$estimate[[1]]
#  	}
#  names(res)<-predictors[2:nPred]
# ============
# ## vs HOME RANGE...
# ============
# # logBodyMass_kg  logTrophic123 logGeoArea_km2 logDispDist_km logGenLength_d 
# #     0.99460067    -0.42617973     0.07759230     0.92544866     0.43563327 
# #   Lat_centroid   Lon_centroid          genes 
# #    -0.07995504     0.04806887     0.21767717 
# ============
#    logTrophic123    logGeoArea_km2    logGenLength_d logDispersalIndex 
#      -0.42617973        0.07759230        0.43563327        0.88391716 
#     Lat_centroid      Lon_centroid             genes 
#      -0.07995504        0.04806887        0.21767717 
# ============
#  logTrophic123  logGeoArea_km2  logGenLength_d logLifemode1234  logActivity123 
#    -0.42813083      0.07716077      0.43354361     -0.36527282      0.34414994 
#   Lat_centroid           genes 
#    -0.08442728      0.22009588 


		form<-as.formula(paste(respVar, " ~ ", predictors[1], "+", predictors[2], "+", predictors[3], "+", predictors[4], "+", predictors[5], "+", predictors[6], sep=""))  #"+", predictors[7], "+", predictors[8], 
		model_ALL<-pgls(form, data=cladeData, lambda=1.0)

		#form<-as.formula(paste(respVar, " ~ ", predictors[1], "+", predictors[2], "+", predictors[3], "+", predictors[4], "+", predictors[5], "+", predictors[6], "+", predictors[7], sep=""))
		#model_8<-pgls(form, data=cladeData, lambda=1.0)

		#form<-as.formula(paste(respVar, " ~ ", predictors[1], "+", predictors[2], "+", predictors[3], "+", predictors[4], "+", predictors[5], "+", predictors[6],  sep=""))
		#model_7<-pgls(form, data=cladeData, lambda=1.0)

		form<-as.formula(paste(respVar, " ~ ", predictors[1], "+", predictors[2], "+", predictors[3], "+", predictors[4], "+", predictors[5], sep=""))
		model_6<-pgls(form, data=cladeData, lambda=1.0)

		form<-as.formula(paste(respVar, " ~ ", predictors[1], "+", predictors[2], "+", predictors[3], "+", predictors[4], "+", predictors[6], sep=""))
		model_5<-pgls(form, data=cladeData, lambda=1.0)

		form<-as.formula(paste(respVar, " ~ ", predictors[1], "+", predictors[2], "+", predictors[3], "+", predictors[5], "+", predictors[6], sep=""))
		model_4<-pgls(form, data=cladeData, lambda=1.0)

		form<-as.formula(paste(respVar, " ~ ", predictors[1], "+", predictors[2], "+", predictors[4], "+", predictors[5], "+", predictors[6], sep=""))
		model_3<-pgls(form, data=cladeData, lambda=1.0)

		form<-as.formula(paste(respVar, " ~ ", predictors[1], "+", predictors[3], "+", predictors[4], "+", predictors[5], "+", predictors[6], sep=""))
		model_2<-pgls(form, data=cladeData, lambda=1.0)

		form<-as.formula(paste(respVar, " ~ ", predictors[2], "+", predictors[3], "+", predictors[4], "+", predictors[5], "+", predictors[6], sep=""))
		model_1<-pgls(form, data=cladeData, lambda=1.0)

		anovaCompare<-as.data.frame(anova.pgls(model_ALL, model_1, model_2, model_3, model_4, model_5, model_6))
		aicCompare<-AIC(model_ALL, model_1, model_2, model_3, model_4, model_5, model_6)

		deltaAIC<-(aicCompare[,2]-aicCompare[1,2]) # diff between model with variable removed and full model
		names(deltaAIC)<-c("ALL",predictors)
		ord<-sort(deltaAIC, decreasing=TRUE)
		nums<-1:length(aicCompare[,1])
		
		res<-cbind(ord,nums)
		rank<-res[match(c("ALL",predictors),rownames(res)),2]

		resp<-rep(respVar,length(predictors)+1)
		pred<-c("ALL",predictors)
		tree<-i

	fullResTIPS<-cbind.data.frame(resp,pred,rank,deltaAIC,anovaCompare,aicCompare,tree)
	rownames(fullResTIPS)<-1:length(aicCompare[,1])

	#corr<-"PAGEL"
	corr<-"BROWNIAN"

	write.table(fullResALL, file=paste("tipLevel",bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-tipDR_multiANOVA_6preds.txt",sep=""))

} # cycles the 100 trees

stopCluster(cl)

q()

n



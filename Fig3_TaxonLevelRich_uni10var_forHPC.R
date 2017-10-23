
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy v1 -- Upham et al. 2017
###
# Figure 3 - time slices to explain clade richness
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# For TAXA, calculate summary values per-TAXON category and per taxon (gen, fam, ord)  
######
# - load back in 1 tree of 100
# - Calculate ES and DR on per-tip basis
# - Calculate per-TAXON category and per taxon (gen, fam, ord)  

# for TAXA, run PGLS -- aim to ** explain variation RICHNESS ** 
######
# - In parallel, read back in slice phylo backbones (14 slices at 5 Ma intervals)
# - Load in age and rate slice summaries for 1 of 100 trees (standardize predictors)
# - Run loop of 4 response vars 
# - Setup empty dataframes to receive PGLS results 
# - Nest the loop of 26 predictor vars



# do the PGLS - UNIVAR and MULTIVAR
# ==================================
# set wd
#dirname<-"/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/cladeLevel_TAXA_explainingRichnessDR"
#setwd(dirname)

# intialize
library(ape); library(phytools); library(picante); library(geiger); library(moments); library(nlme)
library(foreach);library(doSNOW)

# which backbone?
bbone<- "NDexp" #"FBD" # 

# start cluster
cl = makeCluster(100, type = 'MPI')
registerDoSNOW(cl)

ntrees=100
foreach(i=1:ntrees, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

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

taxonomyVars<-c("gen", "fam", "ord", "higher", "genes")
taxFinal<-sorted[,taxonomyVars]

#==================
# Subset by comparison of ecological traits to use
#==================
tipDataAll<-read.table(file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)
	
#varsSelected<-c("Lat_centroid", "Lon_centroid","geoArea_km2", "BM_final_kg", "homeRange_km2_ext", "DispDistAll_km_ext", "GenerationLength_d", )
#catsSelected<-c("CarnOrNot", "HerbOrNot", "OmniOrNot", "AquaOrNot", "ArboOrNot", "FlysOrNot", "SubtOrNot", "TerrOrNot", "NoctOrNot", "CathOrNot", "DiurOrNot", "MarineOrNot")

allDatFinal<-cbind.data.frame(tipDR_i,taxFinal,tipDataAll) #[,varsSelected],tipDataAll[,catsSelected])

allGens<-names(table(allDatFinal$gen)[which(as.vector(table(allDatFinal$gen)) >= 4)]) # 392 genera ≥ 4 species
allFams<-names(table(allDatFinal$fam)[which(as.vector(table(allDatFinal$fam)) >= 4)]) # 105 families ≥ 4 species
allOrds<-names(table(allDatFinal$ord)[which(as.vector(table(allDatFinal$ord)) >= 4)]) # 23 orders ≥ 4 species

taxonSubsets<-list(allGens,allFams,allOrds)
taxonCatVars<-c("gen","fam","ord")

#===================================
# Calculate summary values per-TAXON category and per taxon (gen, fam, ord)  
#===================================
# get node times for tree
btimes<-branching.times(mamPhy)

# yule function
ymle = function(tree){ (.subset2(tree,3)-1L)/sum(.subset2(tree,2)) } # this take the # of number of nodes in a tree (minus 1) / sum of branch lengths.

# do per-slice, per-clade calcs
for(j in 1:length(taxonSubsets)){
taxonCat<-taxonSubsets[[j]]
taxon<-taxonCatVars[j]

	# empty data frames to fill
	taxonReps<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))

	DR_harm<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	#DR_cv<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	DR_skew<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	#DR_kurt<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	percentSamp<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	richness<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	MRCA<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	PB_Div<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	BD_Lam<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	BD_Mu<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	BD_Div<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	BD_Turn<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	#BD.ms0<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	#BD.ms0p5<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	#BD.ms0p9<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	
	BodyMass_kg<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))	
	GeoArea_km2<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	HomeRange_km2<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	rankHomeRange_1to10<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	rankHomeRange_1to10_binsEq<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	DispDist_km<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	rankDispDist_1to10<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	rankDispDist_1to10_binsEq<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	GenLength_d<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	
	perPlant<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	perInvert<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	perVert<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	trophic123<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	HerbOrNot<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	OmniOrNot<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	CarnOrNot<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	
	lifemode1234<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	AquaOrNot<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	ArboOrNot<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	FlysOrNot<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	SubtOrNot<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	TerrOrNot<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))

	activity123<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	NoctOrNot<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	CathOrNot<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))
	DiurOrNot<-data.frame(matrix(NA, nrow = length(taxonCat), ncol = 1))

	for (k in 1:length(taxonCat)){
	cladeSpData<-allDatFinal[which(allDatFinal[,taxon]==taxonCat[k]),]
	cladeSp<-rownames(cladeSpData)
	cladeSet<-drop.tip(mamPhy,(setdiff(mamPhy$tip.label,cladeSp)))
	
	taxonReps[k,]<-cladeSp[1]

		x<-cladeSpData[,"tipDR_i"]
		DR_harm[k,] <- 1/(mean(1/x))
		percentSamp[k,] <- length(which(cladeSpData$genes!=0))/length(cladeSp)
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
		BodyMass_kg[k,]<-mean(na.omit(cladeSpData[,"BM_final_kg"]))
		GeoArea_km2[k,]<-mean(na.omit(cladeSpData[,"geoArea_km2"]))
		HomeRange_km2[k,]<-mean(na.omit(cladeSpData[,"homeRange_km2_ext"]))
		rankHomeRange_1to10[k,]<-mean(na.omit(cladeSpData[,"lnHomeRange_rank1to10"]))
		rankHomeRange_1to10_binsEq[k,]<-mean(na.omit(cladeSpData[,"lnHomeRange_rank1to10_binsEqual"]))
		DispDist_km[k,]<-mean(na.omit(cladeSpData[,"DispDistAll_km_ext"]))
		rankDispDist_1to10[k,]<-mean(na.omit(cladeSpData[,"lnDispDistAll_rank1to10"]))
		rankDispDist_1to10_binsEq[k,]<-mean(na.omit(cladeSpData[,"lnDispDistAll_rank1to10_binsEqual"]))
		GenLength_d[k,]<-mean(na.omit(cladeSpData[,"GenerationLength_d"]))
	# Categorical
		perPlant[k,]<-mean(na.omit(cladeSpData[,"perPlant"]))
		perInvert[k,]<-mean(na.omit(cladeSpData[,"perInvert"]))
		perVert[k,]<-mean(na.omit(cladeSpData[,"perVert"]))
		trophic123[k,]<-mean(na.omit(cladeSpData[,"trophic123"]))
		HerbOrNot[k,]<-mean(na.omit(cladeSpData[,"HerbOrNot"]))
		OmniOrNot[k,]<-mean(na.omit(cladeSpData[,"OmniOrNot"]))
		CarnOrNot[k,]<-mean(na.omit(cladeSpData[,"CarnOrNot"]))
		
		lifemode1234[k,]<-mean(na.omit(cladeSpData[,"lifemode1234"]))
		AquaOrNot[k,]<-mean(na.omit(cladeSpData[,"AquaOrNot"]))
		ArboOrNot[k,]<-mean(na.omit(cladeSpData[,"ArboOrNot"]))
		FlysOrNot[k,]<-mean(na.omit(cladeSpData[,"FlysOrNot"]))
		SubtOrNot[k,]<-mean(na.omit(cladeSpData[,"SubtOrNot"]))
		TerrOrNot[k,]<-mean(na.omit(cladeSpData[,"TerrOrNot"]))

		activity123[k,]<-mean(na.omit(cladeSpData[,"activity123"]))
		NoctOrNot[k,]<-mean(na.omit(cladeSpData[,"NoctOrNot"]))
		CathOrNot[k,]<-mean(na.omit(cladeSpData[,"CathOrNot"]))
		DiurOrNot[k,]<-mean(na.omit(cladeSpData[,"DiurOrNot"]))

	if (length(cladeSp) > 2) {
	# Yule model
		PB_Div[k,]<-ymle(cladeSet)
		# BD model
		bd<-birthdeath(cladeSet)
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

	res2<-cbind.data.frame(taxonCat, DR_harm, DR_skew, percentSamp, richness, MRCA, PB_Div, BD_Lam, BD_Mu, BD_Div, BD_Turn, BodyMass_kg,GeoArea_km2,HomeRange_km2,rankHomeRange_1to10,rankHomeRange_1to10_binsEq,DispDist_km,rankDispDist_1to10,rankDispDist_1to10_binsEq,GenLength_d,perPlant,perInvert,perVert,trophic123,HerbOrNot,OmniOrNot,CarnOrNot,lifemode1234,AquaOrNot,ArboOrNot,FlysOrNot,SubtOrNot,TerrOrNot,activity123,NoctOrNot,CathOrNot,DiurOrNot,i)

	colnames(res2)<-c(taxon, "DR_harm","DR_skew", "percentSamp", "richness", "MRCA", "PB_Div", "BD_Lam", "BD_Mu", "BD_Div", "BD_Turn", "BodyMass_kg", "GeoArea_km2", "HomeRange_km2", "rankHomeRange_1to10", "rankHomeRange_1to10_binsEq", "DispDist_km", "rankDispDist_1to10", "rankDispDist_1to10_binsEq", "GenLength_d", "perPlant","perInvert","perVert","trophic123", "HerbOrNot","OmniOrNot","CarnOrNot","lifemode1234","AquaOrNot","ArboOrNot","FlysOrNot","SubtOrNot","TerrOrNot","activity123","NoctOrNot","CathOrNot","DiurOrNot","tree")

	rownames(res2)<-as.vector(taxonReps[,1])

	write.table(res2,paste("taxonLevel_",taxon,"_",bbone,"_sample100_",i,"_cladeSTATS_withEcoTraits.txt",sep=""))


#===================================
# For TAXA, do PGLS 
#===================================

# Load in TAXON summaries and STANDARDIZE
# ========================================
# Explaining log(richness)
# ========================================

# select RESPONSE vars
	RESP<-log(res2[,"richness"])
	# select PREDICTOR vars
	PRED<-res2[,c(2:11)] # age and rate vars
		# STANDARDIZE predictors as Z-SCORES
		predScale<-scale(PRED, center=TRUE,scale=TRUE)

	# Join together RESP and standardized PRED
	results<-cbind(RESP,predScale[,1:length(PRED)])
	colnames(results)<-c("logRichness",colnames(PRED))

# get BACKBONE phylogeny
cladeData<-treedata(mamPhy,na.omit(results)) 
dat<-as.data.frame(cladeData$data)


# Setup empty dataframes to receive PGLS results 
# ==============================================
# UNIVARIATE first
predictors<-colnames(PRED) # focusing on age and rate
nCols<-length(predictors)

uniPGLS_allSlopes<-data.frame(matrix(NA, nrow = 1, ncol = nCols))
colnames(uniPGLS_allSlopes)<-predictors
uniPGLS_allSEs<-data.frame(matrix(NA, nrow = 1, ncol = nCols))
colnames(uniPGLS_allSEs)<-c(paste("SE_",1:nCols,sep=""))
uniPGLS_allInts<-data.frame(matrix(NA, nrow = 1, ncol = nCols))
colnames(uniPGLS_allInts)<-c(paste("i_",1:nCols,sep=""))
uniPGLS_allPs<-data.frame(matrix(NA, nrow = 1, ncol = nCols))
colnames(uniPGLS_allPs)<-c(paste("p_",1:nCols,sep=""))
uniPGLS_allLams<-data.frame(matrix(NA, nrow = 1, ncol = nCols))
colnames(uniPGLS_allLams)<-c(paste("lam_",1:nCols,sep=""))

uniPGLS<-data.frame(matrix(NA, nrow = 1, ncol = nCols*5),row.names="logRichness")
colnames(uniPGLS)<-c(colnames(uniPGLS_allInts),colnames(uniPGLS_allSlopes),colnames(uniPGLS_allPs),colnames(uniPGLS_allSEs),colnames(uniPGLS_allLams))

# MULTIVARIATE next
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

	write.table(uniPGLS,paste("taxonLevel_",taxon,"_",bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-Richness_SCALED_uniVar_10preds.txt",sep=""))

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

	write.table(multiPGLS,paste("taxonLevel_",taxon,"_",bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-Richness_SCALED_multiVar_age-DRmean-DRskew.txt",sep=""))

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

	write.table(multiPGLS_Per,paste("taxonLevel_",taxon,"_",bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-Richness_SCALED_multiVar_age-DRmean-DRskew-perSamp.txt",sep=""))


} # end 1-tree loop across all 3 taxon categories


} # end 100-tree loop


stopCluster(cl)

q()

n

#===================================
# For TAXA, summarize and plot
#===================================
# Explaining log(richness)
# ========================================
# set wd
dirname<-"/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/cladeLevel_TAXA_explainingRichnessDR"
setwd(dirname)

# intialize
library(ape); library(phytools); library(picante); library(geiger); library(moments); library(nlme)
library(plotrix)

# which backbone?
bbone<- "NDexp" #"FBD" # 
ntrees<-100
taxonCatVars<-c("gen","fam","ord")
ntaxa<-length(taxonCatVars)
corr<-"PAGEL"


# load results by taxon categories for 100 trees
uniPGLS_all<-vector("list",ntaxa)
for(j in 1:length(taxonCatVars)){
taxon<-taxonCatVars[j]

	uniPGLS<-vector("list",ntrees)
	for(i in 1:ntrees){
	res<-cbind(read.table(paste("taxonLevel_",taxon,"_",bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-Richness_SCALED_uniVar_10preds.txt",sep=""), header=TRUE), i)
	uniPGLS[[i]]<-res
	}
	uniPGLS_all[[j]]<-do.call(rbind,uniPGLS)

}

respVar<-"logRichness"

# calculate 95% CIs (and significance) for the slopes per tree & slice
# ====================================================================
predictors<-colnames(uniPGLS_all[[1]])[11:20]
nPred<-length(predictors)

RESP<-uniPGLS_all

dataToPlot_allPreds_allTAXA<-vector("list",ntaxa)

for(j in 1:ntaxa){ # for each gen, fam, ord

taxon<-taxonCatVars[j]
uniDat<-RESP[[j]]

# SEs to calc 95 CIs...
slopesToPlot<-uniDat[,predictors]
pvalsToPlot<-uniDat[,c(paste("p_",(1:nPred),sep=""))]
SEsToPlot<-uniDat[,c(paste("SE_",(1:nPred),sep=""))]
lamsToPlot<-uniDat[,c(paste("lam_",(1:nPred),sep=""))]

	# The 95% CI for each estimate of each var, plotted together
	dataToPlot_allPreds_perTaxon<-vector("list",nPred)

	#for (i in 1:(ntrees-fewerTrees)){ # 
	for (i in 1:ntrees){ # so for a given tree IN a TAXON, taking each of *nPred* vars and calc 95%CI
		slopes_i<-slopesToPlot[i,]
		SEs_i<-SEsToPlot[i,]
		pVals_i<-pvalsToPlot[i,]
		lams_i<-lamsToPlot[i,]
		dataToPlot_k<-vector("list",length(predictors))
		for (k in 1:length(predictors)){
			color<-if(pVals_i[,k] <= 0.05){ "grey" } else { "red" } 
			lowHigh95_k<-cbind.data.frame(slopes_i[,k]-(1.96*SEs_i[,k]),slopes_i[,k],slopes_i[,k]+(1.96*SEs_i[,k]),color,predictors[k],lams_i[,k],taxon,j,i)
			colnames(lowHigh95_k)<-c("low","mean","high","color","xVal","lam","taxon","taxNum","tree")
			dataToPlot_k[[k]]<-lowHigh95_k
		} # end *nPred* loop
		dataToPlot_allPreds_perTaxon[[i]]<-do.call(rbind,dataToPlot_k)
	} # end 100 tree loop
	dataToPlot_allPreds_allTAXA[[j]]<-do.call(rbind,dataToPlot_allPreds_perTaxon)
} # end 10 slice loop
dataToPlot_allPreds_allTAXA_join<-do.call(rbind,dataToPlot_allPreds_allTAXA)


# count the NUMBER of significant runs per TAXON per variable across 100 trees
# ==================================================================================
dat<-dataToPlot_allPreds_allTAXA_join
		
	numSignif_perTaxon<-vector("list",length=ntaxa)
	lamMean_perTaxon<-vector("list",length=ntaxa)
	for(j in 1:ntaxa){
	taxonDat<-dat[which(dat$taxon==taxonCatVars[j]),]

		numSignif<-c()
		lamMean<-c()
		for (k in 1:length(predictors)){
		predPerTaxon<-taxonDat[which(taxonDat$xVal==predictors[k]),]
		num<-length(predPerTaxon[which(predPerTaxon$color=="grey"),][,1])
		#numSignif[k]<-round((num/(ntrees-fewerTrees))*ntrees,0)
		numSignif[k]<-num
		
		lamMean[k]<-round(mean(predPerTaxon[,"lam"]),2)
		}
		numSignif_perTaxon[[j]]<-numSignif
		lamMean_perTaxon[[j]]<-lamMean
	}
	RES1<-do.call(rbind,numSignif_perTaxon)
	colnames(RES1)<-predictors
	rownames(RES1)<-taxonCatVars
	
	RES2<-do.call(rbind,lamMean_perTaxon)
	colnames(RES2)<-predictors
	rownames(RES2)<-taxonCatVars

numSignif_perTaxon_ALL<-RES1
lamMean_perTaxon_ALL<-RES2


# PLOT the 95% CIs per TAXON by predictor comparison across 100 trees
# ====================================================================
corr="PAGEL"
#predictorsORIG<-colnames(exp_N[[1]])[11:20]
predictors<-c("MRCA", "DR_harm", "DR_skew", "percentSamp", "richness", "PB_Div","BD_Div", "BD_Turn", "BD_Lam", "BD_Mu")
nPred<-length(predictors)
respVar<-"logRichness"
nResp<-length(respVar)
nTotal<-nPred*nResp

ntrees<-100
taxonCatVars<-c("gen","fam","ord")
ntaxa<-length(taxonCatVars)

vertLwd<-0.5
lwdCI<-1.5
horizLine<-0
empirPointCol<-grey(0.3,alpha=0.5)
redCol<-rgb(1,0,0,alpha=0.3)
yLims1<-c(-1.5,2.5)


#pdf(file=paste("cladeLevel",bbone,"_PGLSuni_",corr,"_explaining-Richness_timeSlices_",sliceN,"_95pCI_rootwardWithSingle_10vars.pdf",sep=""),onefile=TRUE, width=(8*nResp),height=(2*nPred))
#pdf(file=paste("cladeLevel",bbone,"_PGLSuni_",corr,"_explaining-Richness_timeSlices_",sliceN,"_95pCI_rootwardWithSingle_10vars_WITH-1SIM.pdf",sep=""),onefile=TRUE, width=(8*nResp),height=(2*nPred))
pdf(file=paste("taxonLevel",bbone,"_PGLSuni_",corr,"_explaining-Richness_TAXA-gen-fam-ord_95pCI_10vars.pdf",sep=""),onefile=TRUE, height=(8*nResp),width=(2*nPred))

#quartz(width=(8*nResp),height=(2*nPred))
layout(matrix(c(1:nTotal), nResp*2, nPred/2, byrow = TRUE), widths=rep(4,nTotal), heights=rep(3,nTotal))
par(oma = c(5,4,5,3) + 0.1, mar = c(4,1,1,1) + 0.1)

for(k in 1:length(predictors)){
	uniDat<-dataToPlot_allPreds_allTAXA_join[which(dataToPlot_allPreds_allTAXA_join$xVal==predictors[k]),]

	# plot base (no points)
	plot(formula(mean ~ taxNum), data=uniDat, ylim=yLims1, ylab="",xlab="", yaxt="n",xaxt="n",type="n")
	axis(side=2,at=NULL,labels=TRUE)
	axis(side=1,at=c(1,2,3),labels=c("Genera","Families","Orders"), cex=2)
	if(k==1){ mtext(side=3, line=2, text="UNIVARIATE--TAXA")}

	# add guide lines
	abline(h=horizLine,lty=1, lwd=1, col=grey(0.6, alpha=0.5))

# plot EMPIRICAL data
	dat_UNsignif<-uniDat[which(uniDat$color=="red"),]
	if(length(dat_UNsignif[,1])>0){
	plotCI(x=dat_UNsignif$taxNum, add=TRUE,y=dat_UNsignif$mean,ui=dat_UNsignif$high,li=dat_UNsignif$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=redCol,scol=redCol,pch=1,font.lab=2,cex.axis=1.1,cex.lab=1.1)
	}
	dat_signif<-uniDat[which(uniDat$color=="grey"),]
	if(length(dat_signif[,1])>0){
	plotCI(x=dat_signif$taxNum, add=TRUE,y=dat_signif$mean,ui=dat_signif$high,li=dat_signif$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=empirPointCol,scol=empirPointCol,pch=1,font.lab=2,cex.axis=1.1,cex.lab=1.1)
	}	
	mtext(side=3,text=paste(respVar," ~ "),font=2,adj=0,line=1)
	mtext(side=3,text=predictors[k],font=2,adj=0)
	
	# add numSignif data
	sigDat<-numSignif_perTaxon_ALL[,predictors[k]]
	text(x=c(1:3), y = yLims1[2], cex=0.9,font=2, labels = sigDat, col="dark grey")

	# add LAMBDA data
	lamDat<-lamMean_perTaxon_ALL[,predictors[k]]
	text(x=c(1:3), y = yLims1[1], cex=0.9,font=2, labels = lamDat, col="dark grey")

} # end *nPred* loop

dev.off()



>> Now do for multivariate too-- ultimately I want to ADD IN this to the TIME-SLICE analyses...







# What is the mean age of EACH mammalian family?
#####
# set wd
dirname<-"/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/cladeLevel_TAXA_explainingRichnessDR"
setwd(dirname)
bbone<-"NDexp"

taxonCatVars<-c("gen","fam","ord")
sumTaxonMRCA<-vector("list",length(taxonCatVars))
for(q in 1:length(taxonCatVars)){
taxon<-taxonCatVars[q]

	res1<-vector("list",ntrees)
	for(i in 1:ntrees){
		res2<-read.table(paste("taxonLevel_",taxon,"_",bbone,"_sample100_",i,"_cladeSTATS_withEcoTraits.txt",sep=""), header=TRUE)
		res1[[i]]<-res2[,c("MRCA")]
		names(res1[[i]])<-res2[,taxon]
	}
	allRes<-do.call(cbind,res1)

	res<-vector("list",length(allRes[,1]))
	for(j in 1:length(allRes[,1])){
		taxon<-allRes[j,]
		meanAge<-mean(taxon)
		names(meanAge)<-"mean"
		res[[j]]<-c(meanAge,quantile(taxon,c(0.5,0.025,0.975)))
	}
	sumTaxonMRCA[[q]]<-do.call(rbind,res)
	rownames(sumTaxonMRCA[[q]])<-rownames(allRes)

write.table(sumTaxonMRCA[[q]],file=paste("sumTaxonMRCA_",taxonCatVars[q],"_",bbone,"_sample100.txt",sep=""))
}

famMRCA<-sumTaxonMRCA[[2]]

mean(famMRCA[,"mean"])
#[1] 23.62606
#[1] 23.62606
quantile(famMRCA[,"mean"],c(0.025,0.5,0.975))
## if take means of means
#     2.5%       50%     97.5% 
# 4.667611 14.729245 83.812267 

## if take all 100 trees together
#     2.5%       50%     97.5% 
# 4.048956 14.735034 85.522966



allMRCA<-allRes[,c("MRCA","tree")]
rownames(allMRCA)<-paste(allRes[,taxon],"_",allRes[,"tree"],sep="")



#===================================
# For TAXA, do PGLS 
#===================================

# Load in TAXON summaries and STANDARDIZE
# ========================================
# Explaining log(richness)
# ========================================

# select RESPONSE vars
	RESP<-log(res2[,"richness"])
	# select PREDICTOR vars
	PRED<-res2[,c(2:11)] # age and rate vars
		# STANDARDIZE predictors as Z-SCORES
		predScale<-scale(PRED, center=TRUE,scale=TRUE)

	# Join together RESP and standardized PRED
	results<-cbind(RESP,predScale[,1:length(PRED)])
	colnames(results)<-c("logRichness",colnames(PRED))

# get BACKBONE phylogeny
cladeData<-treedata(mamPhy,na.omit(results)) 
dat<-as.data.frame(cladeData$data)






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













	multiPGLS<-read.table(paste("taxonLevel_",taxon,"_",bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-Richness_SCALED_multiVar_age-DRmean-DRskew.txt",sep=""), header=TRUE)
	multiPGLS_Per<-read.table(paste("taxonLevel_",taxon,"_",bbone,"_sample100_",i,"_PGLS_",corr,"_explaining-Richness_SCALED_multiVar_age-DRmean-DRskew-perSamp.txt",sep=""), header=TRUE)

	resPGLS<-list(uniPGLS,multiPGLS,multiPGLS_Per)


















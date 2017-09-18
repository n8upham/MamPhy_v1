library(ape)
library(phytools)
library(phyloch)
library(moments)
library(nlme)
library(geiger)
library(phylolm)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")

#	- DR harmMean
#	- DR CV
#	- species richness (total)
#	- crown age

cladesDR<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_DR-SUMMARY-expanded_FBD.txt")
head(cladesDR)
colnames(cladesDR)<-c("tiplabel","gen","fam","ord","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")

###
# Get a histogram of DR_harm for each ORDER
ordNames<-names(table(cladesDR$ord))

ordNames<-rownames(res[with(res,order(DR_skew)),]) # is the ORDER names by category (skewness)

plot(density(cladesDR[,"harmMeans"], col="grey", breaks=100, main=NULL, sub=NULL, xlab=NULL))


pdf(file="cladeLevelFBD_DRcompare_DRhist-perORD_bySkew.pdf", onefile=TRUE, width=5, heigh=5)
for(i in 1:length(ordNames)){
	hist(cladesDR[which(cladesDR$ord==ordNames[i]),"harmMeans"], breaks=10, main=NULL, xlab="DR_harmonicMean")
	title(main=ordNames[i], cex.sub=1, sub=paste("richness =",length(cladesDR[which(cladesDR$ord==ordNames[i]),"harmMeans"]),"mean =",round(mean(cladesDR[which(cladesDR$ord==ordNames[i]),"harmMeans"]),digits=3),"skewness =",round(skewness(cladesDR[which(cladesDR$ord==ordNames[i]),"harmMeans"]),digits=3),"kurtosis =",round(kurtosis(cladesDR[which(cladesDR$ord==ordNames[i]),"harmMeans"]),digits=3)))
}
dev.off()


#######
## Make those per-clade summaries:

ordNames<-names(table(cladesDR$ord))

DR_harm<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1), row.names=ordNames)
DR_cv<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1), row.names=ordNames)
DR_skew<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1), row.names=ordNames)
DR_kurt<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1), row.names=ordNames)
richness<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1), row.names=ordNames)

for(i in 1:length(ordNames)){
	DR_harm[i,] <- mean(cladesDR[which(cladesDR$ord==ordNames[i]),"harmMeans"])
	DR_cv[i,] <- mean(cladesDR[which(cladesDR$ord==ordNames[i]),"cv"]*100)
	DR_skew[i,] <- skewness(cladesDR[which(cladesDR$ord==ordNames[i]),"harmMeans"])
	DR_kurt[i,] <- kurtosis(cladesDR[which(cladesDR$ord==ordNames[i]),"harmMeans"])
	richness[i,] <- length(cladesDR[which(cladesDR$ord==ordNames[i]),"cv"])
}

res<-cbind(DR_harm, DR_cv, DR_skew, DR_kurt, richness)

colnames(res)<-c("DR_harm","DR_cv", "DR_skew", "DR_kurt", "richness")

write.table(res,"MamPhy_5911sp_FBD_cladeLevel_ORDS_skewKurt.txt")

#####
## and each FAMILY

famNames<-names(table(cladesDR$fam))

famNames<-rownames(res2[with(res2,order(DR_skew)),]) # is the FAM names by category (skewness)

pdf(file="cladeLevelFBD_DRcompare_DRhist-perFAM_bySkew.pdf", onefile=TRUE, width=5, heigh=5)
for(i in 1:length(famNames)){
	hist(cladesDR[which(cladesDR$fam==famNames[i]),"harmMeans"], breaks=10, main=NULL, xlab="DR_harmonicMean")
	title(main=famNames[i], cex.sub=1, sub=paste("richness=",length(cladesDR[which(cladesDR$fam==famNames[i]),"harmMeans"]),"mean=",round(mean(cladesDR[which(cladesDR$fam==famNames[i]),"harmMeans"]),digits=3),"skewness=",round(skewness(cladesDR[which(cladesDR$fam==famNames[i]),"harmMeans"]),digits=3),"kurtosis=",round(kurtosis(cladesDR[which(cladesDR$fam==famNames[i]),"harmMeans"]),digits=3)))
}
dev.off()
#####

##
famNames<-names(table(cladesDR$fam))

DR_harm<-data.frame(matrix(NA, nrow = length(famNames), ncol = 1), row.names=famNames)
DR_cv<-data.frame(matrix(NA, nrow = length(famNames), ncol = 1), row.names=famNames)
DR_skew<-data.frame(matrix(NA, nrow = length(famNames), ncol = 1), row.names=famNames)
DR_kurt<-data.frame(matrix(NA, nrow = length(famNames), ncol = 1), row.names=famNames)
richness<-data.frame(matrix(NA, nrow = length(famNames), ncol = 1), row.names=famNames)

for(i in 1:length(famNames)){
	DR_harm[i,] <- mean(cladesDR[which(cladesDR$fam==famNames[i]),"harmMeans"])
	DR_cv[i,] <- mean(cladesDR[which(cladesDR$fam==famNames[i]),"cv"]*100)
	DR_skew[i,] <- skewness(cladesDR[which(cladesDR$fam==famNames[i]),"harmMeans"])
	DR_kurt[i,] <- kurtosis(cladesDR[which(cladesDR$fam==famNames[i]),"harmMeans"])
	richness[i,] <- length(cladesDR[which(cladesDR$fam==famNames[i]),"cv"])
}

res2<-cbind(DR_harm, DR_cv, DR_skew, DR_kurt, richness)

colnames(res2)<-c("DR_harm","DR_cv", "DR_skew", "DR_kurt", "richness")

write.table(res2,"MamPhy_5911sp_FBD_cladeLevel_FAMS_skewKurt.txt")


####
# PATCHES

patchNames<-names(table(cladesDR$PC))

patchNames<-rownames(res3[with(res3,order(DR_skew)),]) # is the FAM names by category (skewness)

pdf(file="cladeLevelFBD_DRcompare_DRhist-perPC_bySkew.pdf", onefile=TRUE, width=5, heigh=5)
for(i in 1:length(patchNames)){
	hist(cladesDR[which(cladesDR$PC==patchNames[i]),"harmMeans"], breaks=10, main=NULL, xlab="DR_harmonicMean")
	title(main=patchNames[i], cex.sub=1, sub=paste("richness=",length(cladesDR[which(cladesDR$PC==patchNames[i]),"harmMeans"]),"mean=",round(mean(cladesDR[which(cladesDR$PC==patchNames[i]),"harmMeans"]),digits=3),"skewness=",round(skewness(cladesDR[which(cladesDR$PC==patchNames[i]),"harmMeans"]),digits=3),"kurtosis=",round(kurtosis(cladesDR[which(cladesDR$PC==patchNames[i]),"harmMeans"]),digits=3)))
}
dev.off()

#####
##
patchNames<-names(table(cladesDR$PC))

DR_harm<-data.frame(matrix(NA, nrow = length(patchNames), ncol = 1), row.names=patchNames)
DR_cv<-data.frame(matrix(NA, nrow = length(patchNames), ncol = 1), row.names=patchNames)
DR_skew<-data.frame(matrix(NA, nrow = length(patchNames), ncol = 1), row.names=patchNames)
DR_kurt<-data.frame(matrix(NA, nrow = length(patchNames), ncol = 1), row.names=patchNames)
richness<-data.frame(matrix(NA, nrow = length(patchNames), ncol = 1), row.names=patchNames)

for(i in 1:length(patchNames)){
	DR_harm[i,] <- mean(cladesDR[which(cladesDR$PC==patchNames[i]),"harmMeans"])
	DR_cv[i,] <- mean(cladesDR[which(cladesDR$PC==patchNames[i]),"cv"]*100)
	DR_skew[i,] <- skewness(cladesDR[which(cladesDR$PC==patchNames[i]),"harmMeans"])
	DR_kurt[i,] <- kurtosis(cladesDR[which(cladesDR$PC==patchNames[i]),"harmMeans"])
	richness[i,] <- length(cladesDR[which(cladesDR$PC==patchNames[i]),"cv"])
}

res3<-cbind(DR_harm, DR_cv, DR_skew, DR_kurt, richness)

colnames(res3)<-c("DR_harm","DR_cv", "DR_skew", "DR_kurt", "richness")

write.table(res3,"MamPhy_5911sp_FBD_cladeLevel_PCs_skewKurt.txt")




###
# ADDING the MRCA data

ordRes<-read.table("MamPhy_5911sp_FBD_cladeLevel_ORDS_skewKurt.txt")
famRes<-read.table("MamPhy_5911sp_FBD_cladeLevel_FAMS_skewKurt.txt")
patchRes<-read.table("MamPhy_5911sp_FBD_cladeLevel_PCs_skewKurt.txt")

ordNames<-rownames(ordRes)
famNames<-rownames(famRes)
patchNames<-rownames(patchRes)

mamPhy<-read.beast("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target.tre") ## Using the TARGET node heights = from MCC topology tree  
btimes<-branching.times(mamPhy)

# ORDS
MRCA_mean<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1), row.names=ordNames)
MRCA_mcc<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1), row.names=ordNames)

for(i in 1:23){
	node <- getMRCA(mamPhy, as.vector(cladesDR[which(cladesDR$ord==ordNames[i]),"tiplabel"]))
	MRCA_mean[i,] <- mamPhy$height[node-5912] #taking the MEAN height
	MRCA_mcc[i,] <- btimes[node-5912] #taking the MCC height
}
	
ordRes2<-cbind(ordRes, MRCA_mean[with(MRCA_mean,match(rownames(ordRes),rownames(MRCA_mean))),], MRCA_mcc[with(MRCA_mcc,match(rownames(ordRes),rownames(MRCA_mcc))),])

colnames(ordRes2)<-c("DR_harm","DR_cv", "DR_skew", "DR_kurt", "richness","MRCA_mean","MRCA_mcc")

write.table(ordRes2,"MamPhy_5911sp_FBD_cladeLevel_ORDS_skewKurt_withMRCAs.txt")

# FAMS
MRCA_mean<-data.frame(matrix(NA, nrow = length(famNames), ncol = 1), row.names=famNames)
MRCA_mcc<-data.frame(matrix(NA, nrow = length(famNames), ncol = 1), row.names=famNames)

famsGrThan1 <- rownames(famRes[which(famRes$richness!=1),])
#famsEq1 <- rownames(famRes[which(famRes$richness==1),])

for(i in 1:length(famsGrThan1)){
	node <- getMRCA(mamPhy, as.vector(cladesDR[which(cladesDR$fam==famNames[i]),"tiplabel"]))
	MRCA_mean[i,] <- mamPhy$height[node-5912] #taking the MEAN height
	MRCA_mcc[i,] <- btimes[node-5912] #taking the MCC height
}

famRes2<-cbind(famRes, MRCA_mean[with(MRCA_mean,match(rownames(famRes),rownames(MRCA_mean))),], MRCA_mcc[with(MRCA_mcc,match(rownames(famRes),rownames(MRCA_mcc))),])

colnames(famRes2)<-c("DR_harm","DR_cv", "DR_skew", "DR_kurt", "richness","MRCA_mean","MRCA_mcc")

write.table(famRes2,"MamPhy_5911sp_FBD_cladeLevel_FAMS_skewKurt_withMRCAs.txt")

# PATCHES
MRCA_mean<-data.frame(matrix(NA, nrow = length(patchNames), ncol = 1), row.names=patchNames)
MRCA_mcc<-data.frame(matrix(NA, nrow = length(patchNames), ncol = 1), row.names=patchNames)

for(i in 1:length(patchNames)){
	node <- getMRCA(mamPhy, as.vector(cladesDR[which(cladesDR$PC==patchNames[i]),"tiplabel"]))
	MRCA_mean[i,] <- mamPhy$height[node-5912] #taking the MEAN height
	MRCA_mcc[i,] <- btimes[node-5912] #taking the MCC height
}
	
#mamPhy$height[8206-5912] #PC18 MRCA
#mamPhy$height[8205-5912] #PC18-PC16 MRCA = 33.86435
#mamPhy$height[8275-5912] #PC16 MRCA

# So then set those nodes to 33.85
MRCA_mean[8,]<-33.85
MRCA_mean[10,]<-33.85

patchRes2<-cbind(patchRes, MRCA_mean[with(MRCA_mean,match(rownames(patchRes),rownames(MRCA_mean))),], MRCA_mcc[with(MRCA_mcc,match(rownames(patchRes),rownames(MRCA_mcc))),])

colnames(patchRes2)<-c("DR_harm","DR_cv", "DR_skew", "DR_kurt", "richness","MRCA_mean","MRCA_mcc")

write.table(patchRes2,"MamPhy_5911sp_FBD_cladeLevel_PCs_skewKurt_withMRCAs.txt")


###
## Now prune big tree to ORDER reps, FAM reps, PC reps
ordSp<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_1spMostGenes_perORD.txt")
ordSpOrd<-cbind(ordSp[1],ordSp[4])
colnames(ordSpOrd)<-c("tip","ord")

famSp<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_1spMostGenes_perFAM.txt")
famSpFam<-cbind(famSp[1],famSp[3])
colnames(famSpFam)<-c("tip","fam")

patchSp<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_1spMostGenes_perPC.txt")
patchSpPatch<-cbind(patchSp[1],patchSp[7])
colnames(patchSpPatch)<-c("tip","PC")

ordPhy<-drop.tip(mamPhy,setdiff(mamPhy$tip.label,ordSpOrd$tip))
write.tree(ordPhy, "MamPhy_BDvr_pcsFIXED_FBD_MCC_target_27ORDERS_fullName.tre")
ordPhy$tip.label<-as.vector(ordSpOrd[match(ordPhy$tip.label, as.vector(ordSpOrd$tip)),"ord"])
write.tree(ordPhy, "MamPhy_BDvr_pcsFIXED_FBD_MCC_target_27ORDERS.tre")

famPhy<-drop.tip(mamPhy,setdiff(mamPhy$tip.label,famSpFam$tip))
write.tree(famPhy, "MamPhy_BDvr_pcsFIXED_FBD_MCC_target_162FAMILIES_fullName.tre")
famPhy$tip.label<-as.vector(famSpFam[match(famPhy$tip.label, as.vector(famSpFam$tip)),"fam"])
write.tree(famPhy, "MamPhy_BDvr_pcsFIXED_FBD_MCC_target_162FAMILIES.tre")

patchPhy<-drop.tip(mamPhy,setdiff(mamPhy$tip.label,patchSpPatch$tip))
write.tree(patchPhy, "MamPhy_BDvr_pcsFIXED_FBD_MCC_target_28PATCHES_fullName.tre")
patchPhy$tip.label<-as.vector(patchSpPatch[match(patchPhy$tip.label, as.vector(patchSpPatch$tip)),"PC"])
write.tree(patchPhy, "MamPhy_BDvr_pcsFIXED_FBD_MCC_target_28PATCHES.tre")


###
# NON-PHYLO -- GLMS
###
library(ape)
library(geiger)
library(phytools)
library(phylolm)
library(MCMCglmm) #vignette("Overview", "MCMCglmm")
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")

ordRes<-read.table("MamPhy_5911sp_FBD_cladeLevel_ORDS_skewKurt_withMRCAs.txt")
famRes<-read.table("MamPhy_5911sp_FBD_cladeLevel_FAMS_skewKurt_withMRCAs.txt")
patchRes<-read.table("MamPhy_5911sp_FBD_cladeLevel_PCs_skewKurt_withMRCAs.txt")

ordNames<-rownames(ordRes)
famNames<-rownames(famRes)
patchNames<-rownames(patchRes)

ordPhy<-read.tree("MamPhy_BDvr_pcsFIXED_FBD_MCC_target_27ORDERS.tre")
famPhy<-read.tree("MamPhy_BDvr_pcsFIXED_FBD_MCC_target_162FAMILIES.tre")
patchPhy<-read.tree("MamPhy_BDvr_pcsFIXED_FBD_MCC_target_28PATCHES.tre")

##
# GLM with partial residuals
##
# pRes.srrel <- resid(srrel.reg, type='partial') # MARTA's way of getting the partial residuals...
# pRes <- resid(fit, type='partial')
##
# ORDS - start HERE

ordData<-na.omit(ordRes)	
# >> plot the AIC and SLOPE, rather than the Pval-- all the same 

# richness ~ MRCA
form=(richness ~ MRCA_mean)
fit<-glm(form, data=ordData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byORD_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

pdf(file="cladeLevelFBD_GLMs_byORD_withPartialResiduals.pdf", onefile=TRUE, width=5, heigh=5)

termplot(fit, partial.resid=TRUE, cex.main=0.8, main=paste("23 ords; ",ss[2],ss[1],ss[3], "; X=", round(dd$coef[[2]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ DR_harm
form=(richness ~ DR_harm)
fit<-glm(form, data=ordData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byORD_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, cex.main=0.8, main=paste("23 ords; ",ss[2],ss[1],ss[3], "; X=", round(dd$coef[[2]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ DR_skew
form=(richness ~ DR_skew)
fit<-glm(form, data=ordData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byORD_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, cex.main=0.8, main=paste("23 ords; ",ss[2],ss[1],ss[3], "; X=", round(dd$coef[[2]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ MRCA_mean + DR_harm
form=(richness ~ MRCA_mean + DR_harm)
fit<-glm(form, data=ordData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byORD_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, cex.main=0.65, main=paste("23 ords; ",ss[2],ss[1],ss[3], "; X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ MRCA_mean + DR_skew
form=(richness ~ MRCA_mean + DR_skew)
fit<-glm(form, data=ordData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byORD_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, cex.main=0.65, main=paste("23 ords; ",ss[2],ss[1],ss[3], "; X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ MRCA_mean + DR_harm + DR_skew
form=(richness ~ MRCA_mean + DR_harm + DR_skew)
fit<-glm(form, data=ordData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byORD_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, cex.main=0.55, main=paste("23 ords; ",ss[2],ss[1],ss[3], "; X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; X3=", round(dd$coef[[4]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ MRCA_mean + DR_harm + DR_skew + DR_cv
form=(richness ~ MRCA_mean + DR_harm + DR_skew + DR_cv)
fit<-glm(form, data=ordData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byORD_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, cex.main=0.5, main=paste("23 ords; ",ss[2],ss[1],ss[3], "; X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; X3=", round(dd$coef[[4]],digits=3), "; X4=", round(dd$coef[[5]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

dev.off()


##
# FAMS - start HERE
famData<-na.omit(famRes)	

# richness ~ MRCA
form=(richness ~ MRCA_mean)
fit<-glm(form, data=famData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byFAM_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

pdf(file="cladeLevelFBD_GLMs_byFAM_withPartialResiduals.pdf", onefile=TRUE, width=5, heigh=5)

termplot(fit, partial.resid=TRUE, cex.main=0.8, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; X=", round(dd$coef[[2]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ DR_harm
form=(richness ~ DR_harm)
fit<-glm(form, data=famData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byFAM_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, cex.main=0.8, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; X=", round(dd$coef[[2]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ DR_skew
form=(richness ~ DR_skew)
fit<-glm(form, data=famData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byFAM_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, cex.main=0.8, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; X=", round(dd$coef[[2]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ MRCA_mean + DR_harm
form=(richness ~ MRCA_mean + DR_harm)
fit<-glm(form, data=famData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byFAM_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, cex.main=0.65, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ MRCA_mean + DR_skew
form=(richness ~ MRCA_mean + DR_skew)
fit<-glm(form, data=famData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byFAM_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, cex.main=0.65, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ MRCA_mean + DR_harm + DR_skew
form=(richness ~ MRCA_mean + DR_harm + DR_skew)
fit<-glm(form, data=famData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byFAM_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, cex.main=0.55, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; X3=", round(dd$coef[[4]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ MRCA_mean + DR_harm + DR_skew + DR_cv
form=(richness ~ MRCA_mean + DR_harm + DR_skew + DR_cv)
fit<-glm(form, data=famData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byFAM_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, cex.main=0.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; X3=", round(dd$coef[[4]],digits=3), "; X4=", round(dd$coef[[5]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))


dev.off()

##
# PATCHES - start HERE
patchData<-na.omit(patchRes)	

# richness ~ MRCA
form=(richness ~ MRCA_mean)
fit<-glm(form, data=patchData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byPC_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

pdf(file="cladeLevelFBD_GLMs_byPC_withPartialResiduals.pdf", onefile=TRUE, width=5, heigh=5)

termplot(fit, partial.resid=TRUE, cex.main=0.8, main=paste("28 PCs; ",ss[2],ss[1],ss[3], "; X=", round(dd$coef[[2]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ DR_harm
form=(richness ~ DR_harm)
fit<-glm(form, data=patchData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byPC_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, cex.main=0.8, main=paste("28 PCs; ",ss[2],ss[1],ss[3], "; X=", round(dd$coef[[2]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ DR_skew
form=(richness ~ DR_skew)
fit<-glm(form, data=patchData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byPC_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, cex.main=0.8, main=paste("28 PCs; ",ss[2],ss[1],ss[3], "; X=", round(dd$coef[[2]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ MRCA_mean + DR_harm
form=(richness ~ MRCA_mean + DR_harm)
fit<-glm(form, data=patchData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byPC_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, cex.main=0.65, main=paste("28 PCs; ",ss[2],ss[1],ss[3], "; X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ MRCA_mean + DR_skew
form=(richness ~ MRCA_mean + DR_skew)
fit<-glm(form, data=patchData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byPC_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, cex.main=0.65, main=paste("28 PCs; ",ss[2],ss[1],ss[3], "; X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ MRCA_mean + DR_harm + DR_skew
form=(richness ~ MRCA_mean + DR_harm + DR_skew)
fit<-glm(form, data=patchData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byPC_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, cex.main=0.55, main=paste("28 PCs; ",ss[2],ss[1],ss[3], "; X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; X3=", round(dd$coef[[4]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ MRCA_mean + DR_harm + DR_skew + DR_cv
form=(richness ~ MRCA_mean + DR_harm + DR_skew + DR_cv)
fit<-glm(form, data=patchData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink("cladeLevelFBD_GLMs_byPC_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, cex.main=0.5, main=paste("28 PCs; ",ss[2],ss[1],ss[3], "; X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; X3=", round(dd$coef[[4]],digits=3), "; X4=", round(dd$coef[[5]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

dev.off()


##############
# TIMESLICES
###########
#########
####
# NOW, make TIMESLICES through the tree.
##
library(ape)
library(phytools)
library(phyloch)
library(moments)
library(nlme)
library(geiger)
library(phylolm)
library(paleotree)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")

# Working from the MAMPHY !

mamPhy<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target.tre") ## Use the MCC target node heights one... 

root=max(nodeHeights(mamPhy)[,2])

trees5Ma<-treeSlice(mamPhy, slice=root-5, trivial=FALSE, prompt=FALSE)

trees10Ma<-treeSlice(mamPhy, slice=root-10, trivial=FALSE, prompt=FALSE)

trees20Ma<-treeSlice(mamPhy, slice=root-20, trivial=FALSE, prompt=FALSE)

trees30Ma<-treeSlice(mamPhy, slice=root-30, trivial=FALSE, prompt=FALSE)

trees40Ma<-treeSlice(mamPhy, slice=root-40, trivial=FALSE, prompt=FALSE)

trees50Ma<-treeSlice(mamPhy, slice=root-50, trivial=FALSE, prompt=FALSE)

trees60Ma<-treeSlice(mamPhy, slice=root-60, trivial=FALSE, prompt=FALSE)

lengths<-c(length(trees5Ma),length(trees10Ma),length(trees20Ma),length(trees30Ma),length(trees40Ma),length(trees50Ma),length(trees60Ma))

names(lengths)<-c("5Ma","10Ma","20Ma","30Ma","40Ma","50Ma","60Ma")

# 5Ma 10Ma 20Ma 30Ma 40Ma 50Ma 60Ma 
#1060  545  196   99   65   43   30 

allCladeSets<-list(trees5Ma,trees10Ma,trees20Ma,trees30Ma,trees40Ma,trees50Ma,trees60Ma)
allCladeSetNames<-c("5Ma","10Ma","20Ma","30Ma","40Ma","50Ma","60Ma")

for (i in 1:length(allCladeSets)){
	trees<-allCladeSets[[i]]
	write.tree(trees,file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_slice",allCladeSetNames[[i]],"_newick.trees",sep=""))
	}

#cc<-timeSliceTree(famPhy, sliceTime=50, drop.extinct = FALSE, plot = TRUE)
#	# Does the REVERSE of what I'm interested in ! Still interesting later maybe...
#
#cutFams10<-timeslice.phy(famPhy,slice.time=10,plot=TRUE)
#	#old BABST function... actually I think this does the CORRECT thing-- same as PHYTOOLS-- just distance above the root.

# Now want to summarize the PER CLADE values for each of these clades...
#[1] "DR_harm"  "DR_cv"    "DR_skew"  "DR_kurt"  "richness" "MRCA"    

cladesDR<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_DR-SUMMARY-expanded_FBD.txt")

colnames(cladesDR)<-c("tiplabel","gen","fam","ord","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")

#######
## Make those per-clade summaries:

# Have to get the TIPLABELS for each of the taxa in EACH tree of the time slice...
# use mamPhy as above
btimes<-branching.times(mamPhy)

# 60 Ma slice
##
cladeSet<-trees60Ma

DR_harm<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
DR_cv<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
DR_skew<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
DR_kurt<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
richness<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
MRCA_mcc<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))

for (i in 1:length(cladeSet)){
cladeSp<-cladeSet[[i]]$tip.label
	DR_harm[i,] <- mean(na.omit(cladesDR[match(cladesDR$tiplabel,cladeSp),"harmMeans"]))
	DR_cv[i,] <- mean(na.omit(cladesDR[match(cladesDR$tiplabel,cladeSp),"cv"]*100))
	DR_skew[i,] <- skewness(na.omit(cladesDR[match(cladesDR$tiplabel,cladeSp),"harmMeans"]))
	DR_kurt[i,] <- kurtosis(na.omit(cladesDR[match(cladesDR$tiplabel,cladeSp),"harmMeans"]))
	richness[i,] <- length(cladeSp)
	node <- getMRCA(mamPhy, cladeSp)
	MRCA_mcc[i,] <- btimes[node-5912] #taking the MCC height
}

res<-cbind(DR_harm, DR_cv, DR_skew, DR_kurt, richness, MRCA_mcc)

colnames(res)<-c("DR_harm","DR_cv", "DR_skew", "DR_kurt", "richness", "MRCA_mcc")

write.table(res,"MamPhy_5911sp_FBD_cladeLevel_DRstats_trees60Ma.txt")

##
# Now LOOPING across all these SLICES
##
allCladeSets<-list(trees5Ma,trees10Ma,trees20Ma,trees30Ma,trees40Ma,trees50Ma,trees60Ma)
allCladeSetNames<-c("5Ma","10Ma","20Ma","30Ma","40Ma","50Ma","60Ma")

for(j in 1:length(allCladeSets)){
cladeSet<-allCladeSets[[j]]

	DR_harm<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	DR_cv<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	DR_skew<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	DR_kurt<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	richness<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	MRCA_mcc<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))

	for (i in 1:length(cladeSet)){
	cladeSp<-cladeSet[[i]]$tip.label
		DR_harm[i,] <- mean(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"])
		DR_cv[i,] <- mean(cladesDR[match(cladeSp,cladesDR$tiplabel),"cv"]*100)
		DR_skew[i,] <- skewness(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"])
		DR_kurt[i,] <- kurtosis(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"])
		richness[i,] <- length(cladeSp)
		node <- getMRCA(mamPhy, cladeSp)
		MRCA_mcc[i,] <- btimes[node-5912] #taking the MCC height
	}

	res<-cbind(DR_harm, DR_cv, DR_skew, DR_kurt, richness, MRCA_mcc)

	colnames(res)<-c("DR_harm","DR_cv", "DR_skew", "DR_kurt", "richness", "MRCA_mcc")

	write.table(res,paste("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees", allCladeSetNames[j], ".txt", sep=""))

}

## Now LOAD them back in !

res60Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees60Ma.txt", header=TRUE)
res50Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees50Ma.txt", header=TRUE)
res40Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees40Ma.txt", header=TRUE)
res30Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees30Ma.txt", header=TRUE)
res20Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees20Ma.txt", header=TRUE)
res10Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees10Ma.txt", header=TRUE)
res5Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees5Ma.txt", header=TRUE)

results<-list(res5Ma, res10Ma, res20Ma, res30Ma, res40Ma, res50Ma, res60Ma)

# Per-clade time slice histograms
##

results[[1]][match(order(results[[1]]$DR_skew),rownames(results[[1]])),]

allCladeSets[[7]]
allCladeSetNames[[7]]
results

for (k in 1:length(results)){

#cladeNames<-as.numeric(rownames(results[[k]][with(results[[k]],order(results[[k]]$DR_skew)),]))
cladeNames<-order(results[[k]]$DR_skew)

pdf(file=paste("cladeLevelFBD_DRcompare_DRhist_bySkew_trees", allCladeSetNames[[k]], ".pdf",sep=""), onefile=TRUE, width=5, height=5)

for(j in 1:length(cladeNames)){
	i<-cladeNames[j]

	cladeSp<-allCladeSets[[k]][[i]]$tip.label

	hist(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"], breaks=10, main=NULL, xlab="DR_harmonicMean")
	title(cex.main=1.2,main=paste("slice ",allCladeSetNames[[k]],", clade ",i,": ",cladesDR[match(cladeSp,cladesDR$tiplabel),"gen"][1]," to ",cladesDR[match(cladeSp,cladesDR$tiplabel),"gen"][length(cladesDR[match(cladeSp,cladesDR$tiplabel),"gen"])],sep=""), cex.sub=1, sub=paste("richness =",length(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"]),"mean =",round(mean(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"]),digits=3),"skewness =",round(skewness(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"]),digits=3),"kurtosis =",round(kurtosis(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"]),digits=3)))
	mtext(cex=0.9,paste("MRCA = ",round(results[[k]][i,]$MRCA_mcc,digits=2)," Ma, ",length(which(cladesDR[match(cladeSp,cladesDR$tiplabel),"samp"]=="sampled")),"/",length(cladesDR[match(cladeSp,cladesDR$tiplabel),"samp"]=="sampled")," spp samp; ", cladesDR[match(cladeSp,cladesDR$tiplabel),"fam"][1],", ",cladesDR[match(cladeSp,cladesDR$tiplabel),"ord"][1],sep=""))
}

dev.off()

}

## NOW as HIST with the x-axis all the SAME...
mean=mean(cladesDR$harmMeans)
Q1=quantile(cladesDR$harmMeans,c(0.025,0.975))[[1]]
Q2=quantile(cladesDR$harmMeans,c(0.025,0.975))[[2]]

for (k in 1:length(results)){

cladeNames<-order(results[[k]]$DR_skew)

pdf(file=paste("cladeLevelFBD_DRcompare_DRhist_bySkew_Xsame_trees", allCladeSetNames[[k]], ".pdf",sep=""), onefile=TRUE, width=5, height=5)

for(j in 1:length(cladeNames)){
	i<-cladeNames[j]

	cladeSp<-allCladeSets[[k]][[i]]$tip.label

	hist(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"], breaks=10, main=NULL, xlab="DR_harmonicMean", xlim=c(0.01,1.07))
	abline(v=mean, lty=1, lwd=1, col="blue")
	abline(v=Q1, lty=2, lwd=1, col="black")
	abline(v=Q2, lty=2, lwd=1, col="black")

	title(cex.main=1.2,main=paste("slice ",allCladeSetNames[[k]],", clade ",i,": ",cladesDR[match(cladeSp,cladesDR$tiplabel),"gen"][1]," to ",cladesDR[match(cladeSp,cladesDR$tiplabel),"gen"][length(cladesDR[match(cladeSp,cladesDR$tiplabel),"gen"])],sep=""), cex.sub=1, sub=paste("richness =",length(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"]),"mean =",round(mean(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"]),digits=3),"skewness =",round(skewness(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"]),digits=3),"kurtosis =",round(kurtosis(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"]),digits=3)))
	mtext(cex=0.9,paste("MRCA = ",round(results[[k]][i,]$MRCA_mcc,digits=2)," Ma, ",length(which(cladesDR[match(cladeSp,cladesDR$tiplabel),"samp"]=="sampled")),"/",length(cladesDR[match(cladeSp,cladesDR$tiplabel),"samp"]=="sampled")," spp samp; ", cladesDR[match(cladeSp,cladesDR$tiplabel),"fam"][1],", ",cladesDR[match(cladeSp,cladesDR$tiplabel),"ord"][1],sep=""))
}

dev.off()

}

## And as DENSITY plots with MEAN of MAMALIA in the background...
mean=mean(cladesDR$harmMeans)
Q1=quantile(cladesDR$harmMeans,c(0.025,0.975))[[1]]
Q2=quantile(cladesDR$harmMeans,c(0.025,0.975))[[2]]

for (k in 1:length(results)){

cladeNames<-order(results[[k]]$DR_skew)

pdf(file=paste("cladeLevelFBD_DRcompare_DRhist_bySkew_Xsame-DensMamMean_trees", allCladeSetNames[[k]], ".pdf",sep=""), onefile=TRUE, width=5, height=5)

for(j in 1:length(cladeNames)){
	i<-cladeNames[j]

	cladeSp<-allCladeSets[[k]][[i]]$tip.label

	plot(density(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"]), xlab="DR_harmonicMean", xlim=range(cladesDR$harmMeans), cex.main=1.2,main=paste("slice ",allCladeSetNames[[k]],", clade ",i,": ",cladesDR[match(cladeSp,cladesDR$tiplabel),"gen"][1]," to ",cladesDR[match(cladeSp,cladesDR$tiplabel),"gen"][length(cladesDR[match(cladeSp,cladesDR$tiplabel),"gen"])],sep=""))
	polygon(density(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"]), col="light grey", border="black", bty="n", xlab="")
	abline(v=mean, lty=1, lwd=1, col="blue")
	abline(v=Q1, lty=2, lwd=1, col="black")
	abline(v=Q2, lty=2, lwd=1, col="black")
	title(cex.sub=1, sub=paste("richness =",length(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"]),"mean =",round(mean(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"]),digits=3),"skewness =",round(skewness(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"]),digits=3),"kurtosis =",round(kurtosis(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"]),digits=3)))
	mtext(cex=0.9,paste("MRCA = ",round(results[[k]][i,]$MRCA_mcc,digits=2)," Ma, ",length(which(cladesDR[match(cladeSp,cladesDR$tiplabel),"samp"]=="sampled")),"/",length(cladesDR[match(cladeSp,cladesDR$tiplabel),"samp"]=="sampled")," spp samp; ", cladesDR[match(cladeSp,cladesDR$tiplabel),"fam"][1],", ",cladesDR[match(cladeSp,cladesDR$tiplabel),"ord"][1],sep=""))
}

dev.off()

}

## And as DENSITY plots with DISTRIBUTION of MAMALIA in the background...
x.tick <- as.vector(quantile(cladesDR$harmMeans, c(0.01,0.5,0.99,1)))
dens.rate <- density(cladesDR$harmMeans)$y

for (k in 1:length(results)){

cladeNames<-order(results[[k]]$DR_skew)

pdf(file=paste("cladeLevelFBD_DRcompare_DRhist_bySkew_Xsame-DensMamDist_trees", allCladeSetNames[[k]], ".pdf",sep=""), onefile=TRUE, width=5, height=5)

for(j in 1:length(cladeNames)){
	i<-cladeNames[j]

	cladeSp<-allCladeSets[[k]][[i]]$tip.label

	plot(density(cladesDR$harmMeans), col="light grey", main="", bty="n", axes=F, xlim=range(cladesDR$harmMeans), ylim=range(dens.rate), xlab="", ylab="")
	polygon(density(cladesDR$harmMeans), col="light grey", border="light grey", bty="n", xlab="")
	axis(at=c(0,x.tick), labels=c(0,round(x.tick,2)), side=1, line=1.3, cex=1, lwd=1, tck=-0.05, cex.axis=0.8, mgp=c(-1,-1,-1))
	axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=1, las=1, lwd=1, cex.axis=1, tck=-0.05, mgp=c(1,1,0))
	
	par(new=TRUE)
	richness=length(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"])
	plot(density(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"]), xlab="DR_harmonicMean", xlim=range(cladesDR$harmMeans), ylim=c(0,(5911/richness)*max(dens.rate)),cex.main=1.2,main=paste("slice ",allCladeSetNames[[k]],", clade ",i,": ",cladesDR[match(cladeSp,cladesDR$tiplabel),"gen"][1]," to ",cladesDR[match(cladeSp,cladesDR$tiplabel),"gen"][length(cladesDR[match(cladeSp,cladesDR$tiplabel),"gen"])],sep=""), col="black", lty=2, lwd=3, bty="n", axes=FALSE)
	
	title(cex.sub=1, sub=paste("richness =",richness,"mean =",round(mean(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"]),digits=3),"skewness =",round(skewness(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"]),digits=3),"kurtosis =",round(kurtosis(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"]),digits=3)))
	mtext(cex=0.9,paste("MRCA = ",round(results[[k]][i,]$MRCA_mcc,digits=2)," Ma, ",length(which(cladesDR[match(cladeSp,cladesDR$tiplabel),"samp"]=="sampled")),"/",length(cladesDR[match(cladeSp,cladesDR$tiplabel),"samp"]=="sampled")," spp samp; ", cladesDR[match(cladeSp,cladesDR$tiplabel),"fam"][1],", ",cladesDR[match(cladeSp,cladesDR$tiplabel),"ord"][1],sep=""))
	}

	dev.off()

	}

# >> histos DONE. OKAY, now to the GLM analyses::
##########
#####
# 
res60Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees60Ma.txt", header=TRUE)
res50Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees50Ma.txt", header=TRUE)
res40Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees40Ma.txt", header=TRUE)
res30Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees30Ma.txt", header=TRUE)
res20Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees20Ma.txt", header=TRUE)
res10Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees10Ma.txt", header=TRUE)
res5Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees5Ma.txt", header=TRUE)

results<-list(res5Ma, res10Ma, res20Ma, res30Ma, res40Ma, res50Ma, res60Ma)


allCladeSetNames<-c("5Ma","10Ma","20Ma","30Ma","40Ma","50Ma","60Ma")

for (i in 1:length(allCladeSetNames))
    {
    assign(paste("trees",allCladeSetNames[[i]],sep=""), read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_slice",allCladeSetNames[[i]],"_newick.trees",sep=""))) 
    }

allCladeSets<-list(trees5Ma,trees10Ma,trees20Ma,trees30Ma,trees40Ma,trees50Ma,trees60Ma)

##
# GLMS
# Loop through each of these for the TIME-SLICES
##


# PATCHES - start HERE
for (i in 1:length(results)){

cladeData<-na.omit(results[[i]])	

# richness ~ MRCA
form=(richness ~ MRCA_mcc)
fit<-glm(form, data=cladeData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink(file=paste("cladeLevelFBD_GLMs_bySLICE_", allCladeSetNames[i],"_RESULTS.txt",sep=""), append=TRUE)
print(form)
print(dd)
sink()
ss<-as.character(form)

pdf(file=paste("cladeLevelFBD_GLMs_bySLICE_", allCladeSetNames[i],"_withPartialResiduals.pdf",sep=""), onefile=TRUE, width=5, heigh=5)

termplot(fit, partial.resid=TRUE, rug=TRUE, cex.main=1, font.sub=2, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1]),"; ",ss[2],ss[1],ss[3], sep=""), sub=paste("X1=", round(dd$coef[[2]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ DR_harm
form=(richness ~ DR_harm)
fit<-glm(form, data=cladeData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink(file=paste("cladeLevelFBD_GLMs_bySLICE_", allCladeSetNames[i],"_RESULTS.txt",sep=""), append=TRUE)
print(form)
print(dd)
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, rug=TRUE, cex.main=1, font.sub=2, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1]),"; ",ss[2],ss[1],ss[3], sep=""), sub=paste("X1=", round(dd$coef[[2]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ DR_skew
form=(richness ~ DR_skew)
fit<-glm(form, data=cladeData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink(file=paste("cladeLevelFBD_GLMs_bySLICE_", allCladeSetNames[i],"_RESULTS.txt",sep=""), append=TRUE)
print(form)
print(dd)
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, rug=TRUE, cex.main=1, font.sub=2, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1]),"; ",ss[2],ss[1],ss[3], sep=""), sub=paste("X1=", round(dd$coef[[2]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ MRCA_mcc + DR_harm
form=(richness ~ MRCA_mcc + DR_harm)
fit<-glm(form, data=cladeData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink(file=paste("cladeLevelFBD_GLMs_bySLICE_", allCladeSetNames[i],"_RESULTS.txt",sep=""), append=TRUE)
print(form)
print(dd)
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, rug=TRUE, cex.main=1, font.sub=2, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1]),"; ",ss[2],ss[1],ss[3], sep=""), sub=paste("X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))


# richness ~ MRCA_mcc + DR_skew
form=(richness ~ MRCA_mcc + DR_skew)
fit<-glm(form, data=cladeData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink(file=paste("cladeLevelFBD_GLMs_bySLICE_", allCladeSetNames[i],"_RESULTS.txt",sep=""), append=TRUE)
print(form)
print(dd)
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, rug=TRUE, cex.main=1, font.sub=2, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1]),"; ",ss[2],ss[1],ss[3], sep=""), sub=paste("X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ MRCA_mcc + DR_harm + DR_skew
form=(richness ~ MRCA_mcc + DR_harm + DR_skew)
fit<-glm(form, data=cladeData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink(file=paste("cladeLevelFBD_GLMs_bySLICE_", allCladeSetNames[i],"_RESULTS.txt",sep=""), append=TRUE)
print(form)
print(dd)
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, rug=TRUE, cex.main=1, font.sub=2, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1]),"; ",ss[2],ss[1],ss[3], sep=""), sub=paste("X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; X3=", round(dd$coef[[4]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))


# richness ~ MRCA_mcc + DR_harm + DR_skew + DR_cv
form=(richness ~ MRCA_mcc + DR_harm + DR_skew + DR_cv)
fit<-glm(form, data=cladeData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink(file=paste("cladeLevelFBD_GLMs_bySLICE_", allCladeSetNames[i],"_RESULTS.txt",sep=""), append=TRUE)
print(form)
print(dd)
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, rug=TRUE, cex.main=0.8, font.sub=2, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1]),"; ",ss[2],ss[1],ss[3], sep=""), sub=paste("X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; X3=", round(dd$coef[[4]],digits=3), "; X4=", round(dd$coef[[5]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

dev.off()

}


###
##
# PHYLOGENETIC --- PGLS !!
###
# get the PHYs for each time slice set.

for (k in 1:length(allCladeSets)){
cladeReps<-vector()
for (i in 1:length(allCladeSets[[k]])){
	cladeSp<-allCladeSets[[k]][[i]]$tip.label
	cladeReps[i]<-cladeSp[1]
	}
toDrop<-setdiff(mamPhy$tip.label,cladeReps)
assign(paste("slice",allCladeSetNames[[k]],"Phy",sep=""),drop.tip(mamPhy,toDrop))
}

slicePhys<-list(slice5MaPhy,slice10MaPhy,slice20MaPhy,slice30MaPhy,slice40MaPhy,slice50MaPhy,slice60MaPhy)

for(i in 1:length(slicePhys)){
	write.tree(slicePhys[[i]], file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_slicePhy-5-10-20-30-40-50-60.trees", append=TRUE)
}
# 5Ma 10Ma 20Ma 30Ma 40Ma 50Ma 60Ma 
#1060  545  196   99   65   43   30 
##
# Load in data

res60Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees60Ma.txt", header=TRUE)
res50Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees50Ma.txt", header=TRUE)
res40Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees40Ma.txt", header=TRUE)
res30Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees30Ma.txt", header=TRUE)
res20Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees20Ma.txt", header=TRUE)
res10Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees10Ma.txt", header=TRUE)
res5Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees5Ma.txt", header=TRUE)

results<-list(res5Ma, res10Ma, res20Ma, res30Ma, res40Ma, res50Ma, res60Ma)

allCladeSetNames<-c("5Ma","10Ma","20Ma","30Ma","40Ma","50Ma","60Ma")

slicePhys<-read.tree("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_slicePhy-5-10-20-30-40-50-60.trees")

for (i in 1:length(allCladeSetNames))
    {
    assign(paste("trees",allCladeSetNames[[i]],sep=""), read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_slice",allCladeSetNames[[i]],"_newick.trees",sep=""))) 
    }

allCladeSets<-list(trees5Ma,trees10Ma,trees20Ma,trees30Ma,trees40Ma,trees50Ma,trees60Ma)
####

# richness ~ MRCA

pdf(file="cladeLevelFBD_slices_by10Ma_PGLS_rich-by-MRCA.pdf", onefile=TRUE, width=5, heigh=5)

for (i in 1:length(results)){

rownames(results[[i]])<-slicePhys[[i]]$tip.label
cladeData<-treedata(slicePhys[[i]],na.omit(results[[i]]))

form<-(log(richness) ~ MRCA_mcc)
fit<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=as.data.frame(cladeData$data), method="ML")
sum<-summary(fit)

plot(form,data=cladeData$data, cex.lab=1.5, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1])," clades", sep=""))
if (sum$tTable[8] < 0.05){
	abline(fit$coef[[1]], fit$coef[[2]], col="red", lwd=2, lty=2)
}
mtext(paste("PGLS p-val = ",round(sum$tTable[8],digits=4),"; X = ",round(fit$coef[[2]],digits=2),"; AIC = ",round(sum$AIC,digits=1),sep=""))

}
dev.off()


# richness ~ DR_harm

pdf(file="cladeLevelFBD_slices_by10Ma_PGLS_rich-by-DR_harm.pdf", onefile=TRUE, width=5, heigh=5)

for (i in 1:length(results)){

rownames(results[[i]])<-slicePhys[[i]]$tip.label
cladeData<-treedata(slicePhys[[i]],na.omit(results[[i]]))

form<-(log(richness) ~ DR_harm)
fit<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=as.data.frame(cladeData$data), method="ML")
sum<-summary(fit)

plot(form,data=cladeData$data, cex.lab=1.5, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1])," clades", sep=""))
if (sum$tTable[8] < 0.05){
	abline(fit$coef[[1]], fit$coef[[2]], col="red", lwd=2, lty=2)
}
mtext(paste("PGLS p-val = ",round(sum$tTable[8],digits=4),"; X = ",round(fit$coef[[2]],digits=2),"; AIC = ",round(sum$AIC,digits=1),sep=""))

}
dev.off()


# richness ~ DR_skew

pdf(file="cladeLevelFBD_slices_by10Ma_PGLS_rich-by-DR_skew.pdf", onefile=TRUE, width=5, heigh=5)

for (i in 1:length(results)){

rownames(results[[i]])<-slicePhys[[i]]$tip.label
cladeData<-treedata(slicePhys[[i]],na.omit(results[[i]]))

form<-(log(richness) ~ DR_skew)
fit<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=as.data.frame(cladeData$data), method="ML")
sum<-summary(fit)

plot(form,data=cladeData$data, cex.lab=1.5, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1])," clades", sep=""))
if (sum$tTable[8] < 0.05){
	abline(fit$coef[[1]], fit$coef[[2]], col="red", lwd=2, lty=2)
}
mtext(paste("PGLS p-val = ",round(sum$tTable[8],digits=4),"; X = ",round(fit$coef[[2]],digits=2),"; AIC = ",round(sum$AIC,digits=1),sep=""))

}
dev.off()


# richness ~ DR_kurt

pdf(file="cladeLevelFBD_slices_by10Ma_PGLS_rich-by-DR_kurt.pdf", onefile=TRUE, width=5, heigh=5)

for (i in 1:length(results)){

rownames(results[[i]])<-slicePhys[[i]]$tip.label
cladeData<-treedata(slicePhys[[i]],na.omit(results[[i]]))

form<-(log(richness) ~ DR_kurt)
fit<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=as.data.frame(cladeData$data), method="ML")
sum<-summary(fit)

plot(form,data=cladeData$data, cex.lab=1.5, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1])," clades", sep=""))
if (sum$tTable[8] < 0.05){
	abline(fit$coef[[1]], fit$coef[[2]], col="red", lwd=2, lty=2)
}
mtext(paste("PGLS p-val = ",round(sum$tTable[8],digits=4),"; X = ",round(fit$coef[[2]],digits=2),"; AIC = ",round(sum$AIC,digits=1),sep=""))

}
dev.off()


# richness ~ DR_cv

pdf(file="cladeLevelFBD_slices_by10Ma_PGLS_rich-by-DR_cv.pdf", onefile=TRUE, width=5, heigh=5)

for (i in 1:length(results)){

rownames(results[[i]])<-slicePhys[[i]]$tip.label
cladeData<-treedata(slicePhys[[i]],na.omit(results[[i]]))

form<-(log(richness) ~ DR_cv)
fit<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=as.data.frame(cladeData$data), method="ML")
sum<-summary(fit)

plot(form,data=cladeData$data, cex.lab=1.5, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1])," clades", sep=""))
if (sum$tTable[8] < 0.05){
	abline(fit$coef[[1]], fit$coef[[2]], col="red", lwd=2, lty=2)
}
mtext(paste("PGLS p-val = ",round(sum$tTable[8],digits=4),"; X = ",round(fit$coef[[2]],digits=2),"; AIC = ",round(sum$AIC,digits=1),sep=""))

}
dev.off()

######
# Now fit ML models to each TIMESLICE clades...
# Yule (empirical), BD (sp + ext, net div), DD 


ymle = function(tree) (.subset2(tree,2)-1L)/sum(.subset2(tree,4)) # this take the # of number of nodes in a tree / sum of branch lengths.














###
# MCMCglmm()
##
# ORDS - start HERE
ordData<-na.omit(ordRes)	
treedata<-treedata(ordPhy,ordData)

# richness ~ MRCA + DR_harm
form=(richness ~ MRCA_mean + DR_harm)
fit<-MCMCglmm(fixed=form, data=ordData, pedigree=treedata$phy, family = "poisson", nodes="ALL", scale=TRUE, nitt=13000, thin=10, burnin=3000, pr=FALSE, pl=FALSE, verbose=TRUE, DIC=TRUE, singular.ok=FALSE, saveX=TRUE,saveZ=TRUE, saveXL=TRUE, slice=FALSE, ginverse=NULL)
dd<-summary(fit)
dd


MCMCglmm(fixed, random=NULL, rcov=~units, family="gaussian", mev=NULL, 
         data,start=NULL, prior=NULL, tune=NULL, pedigree=treedata$phy, nodes="ALL",
         scale=TRUE, nitt=13000, thin=10, burnin=3000, pr=FALSE,
         pl=FALSE, verbose=TRUE, DIC=TRUE, singular.ok=FALSE, saveX=TRUE,
         saveZ=TRUE, saveXL=TRUE, slice=FALSE, ginverse=NULL)


# Example 2: univariate Gaussian model with phylogenetically correlated

    Ainv<-inverseA(ordPhy)$Ainv # inverse matrix of shared phylogenetic history
     
    prior<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))
     
    form=(richness ~ MRCA_mean + DR_harm)

	taxon=rownames(ordData)
    
	ordData2<-cbind(ordData,taxon)
	colnames<-c("DR_harm","DR_cv", "DR_skew", "DR_kurt", "richness","MRCA_mean","MRCA_mcc","taxon")

    model2<-MCMCglmm(form, random=~taxon, ginverse=list(taxon=Ainv), data=ordData2, prior=prior, verbose=FALSE, family="poisson")
     
     plot(model2$VCV)


###
# PHYLOGLM()
##
# ORDS - start HERE
ordData<-na.omit(ordRes)	
treedata<-treedata(ordPhy,ordData)

# richness ~ MRCA + DR_harm
form=(richness ~ MRCA_mean + DR_harm + DR_cv)
fit<-phyloglm(formula=form, data=ordData, phy=treedata$phy, method = "poisson_GEE", start.beta=NULL, boot = 10, full.matrix = TRUE)
dd<-summary(fit)
dd
# FAILS to converge with three vars.


##
# compar.gee approach !!
##
library(ape)
library(geiger)
library(phytools)
library(gee)

#setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors")
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")

ordRes<-read.table("MamPhy_5911sp_FBD_cladeLevel_ORDS_skewKurt_withMRCAs.txt")
famRes<-read.table("MamPhy_5911sp_FBD_cladeLevel_FAMS_skewKurt_withMRCAs.txt")

ordNames<-rownames(ordRes[with(ordRes,order(ordRes$richness, decreasing=TRUE)),]) # is the ORDER names by category (skewness)
famNames<-rownames(famRes[with(famRes,order(famRes$richness, decreasing=TRUE)),]) # is the ORDER names by category (skewness)

ordPhy<-read.tree("MamPhy_BDvr_pcsFIXED_FBD_MCC_target_27ORDERS.tre")
famPhy<-read.tree("MamPhy_BDvr_pcsFIXED_FBD_MCC_target_162FAMILIES.tre")


##
# ORDS - start HERE
ordData<-na.omit(ordRes)	
treedata<-treedata(ordPhy,ordData)

# richness ~ MRCA
form=(richness ~ MRCA)
fit<-compar.gee(form, data=ordData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
sink("cladeLevelFBD_GLMs_byORD_RESULTS.txt", append=TRUE)
fit
dd
sink()
ss<-as.character(form)

pdf(file="cladeLevelFBD_GLMs_byORD.pdf", onefile=TRUE, width=5, heigh=5)

plot(log(richness) ~ MRCA, data=ordData, cex.lab=1.5, main=paste("27 orders; ",ss[2],ss[1],ss[3], "; P = ", round(dd[3],digits=3),sep=""))
if (dd[3] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}
plot(residuals(fit)~MRCA,data=ordData)


# richness ~ DR_harm
form=(richness ~ DR_harm)
fit<-compar.gee(form, data=ordData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
sink("cladeLevelFBD_GLMs_byORD_RESULTS.txt", append=TRUE)
fit
dd
sink()
ss<-as.character(form)

plot(log(richness) ~ DR_harm, data=ordData, cex.lab=1.5, main=paste("27 orders; ",ss[2],ss[1],ss[3], "; P = ", round(dd[3],digits=3),sep=""))
if (dd[3] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}
plot(residuals(fit)~DR_harm,data=ordData)

# richness ~ DR_skew
form=(richness ~ DR_skew)
fit<-compar.gee(form, data=ordData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
sink("cladeLevelFBD_GLMs_byORD_RESULTS.txt", append=TRUE)
fit
dd
sink()
ss<-as.character(form)

plot(log(richness) ~ DR_skew, data=ordData, cex.lab=1.5, main=paste("27 orders; ",ss[2],ss[1],ss[3], "; P = ", round(dd[3],digits=3),sep=""))
if (dd[3] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}
plot(residuals(fit)~DR_skew,data=ordData)


# richness ~ MRCA + DR_harm
form=(richness ~ MRCA + DR_harm)
fit<-compar.gee(form, data=ordData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
sink("cladeLevelFBD_GLMs_byORD_RESULTS.txt", append=TRUE)
fit
dd
sink()
ss<-as.character(form)

plot(log(richness) ~ MRCA, data=ordData, cex.lab=1.5, main=paste("27 orders; ",ss[2],ss[1],ss[3], "; P = ", round(dd[5],digits=3),sep=""))
if (dd[5] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}
plot(log(richness) ~ DR_harm, data=ordData, cex.lab=1.5, main=paste("27 orders; ",ss[2],ss[1],ss[3], "; P = ", round(dd[6],digits=3),sep=""))
if (dd[6] < 0.05){ abline(fit$coef[[1]],fit$coef[[3]])
}
plot(residuals(fit)~MRCA,data=ordData)
plot(residuals(fit)~DR_harm,data=ordData)


# richness ~ MRCA + DR_skew
form=(richness ~ MRCA + DR_skew)
fit<-compar.gee(form, data=ordData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
sink("cladeLevelFBD_GLMs_byORD_RESULTS.txt", append=TRUE)
fit
dd
sink()
ss<-as.character(form)

plot(log(richness) ~ MRCA, data=ordData, cex.lab=1.5, main=paste("27 orders; ",ss[2],ss[1],ss[3], "; P = ", round(dd[5],digits=3),sep=""))
if (dd[5] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}
plot(log(richness) ~ DR_skew, data=ordData, cex.lab=1.5, main=paste("27 orders; ",ss[2],ss[1],ss[3], "; P = ", round(dd[6],digits=3),sep=""))
if (dd[6] < 0.05){ abline(fit$coef[[1]],fit$coef[[3]])
}
plot(residuals(fit)~MRCA,data=ordData)
plot(residuals(fit)~DR_skew,data=ordData)

dev.off()


# richness ~ MRCA + DR_harm + DR_skew
form=(richness ~ MRCA + DR_harm + DR_skew)
fit<-glm(form, data=ordData, family="poisson")#, phy=treedata$phy)

dd<-drop1(fit)
ss<-as.character(form)

plot(log(richness) ~ MRCA, data=ordData, cex.lab=1.5, main=paste("27 orders; ",ss[2],ss[1],ss[3], "; P = ", round(dd[7],digits=3),sep=""))
if (dd[7] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}
plot(log(richness) ~ DR_harm, data=ordData, cex.lab=1.5, main=paste("27 orders; ",ss[2],ss[1],ss[3], "; P = ", round(dd[8],digits=3),sep=""))
if (dd[8] < 0.05){ abline(fit$coef[[1]],fit$coef[[3]])
}
plot(log(richness) ~ DR_skew, data=ordData, cex.lab=1.5, main=paste("27 orders; ",ss[2],ss[1],ss[3], "; P = ", round(dd[9],digits=3),sep=""))
if (dd[9] < 0.05){ abline(fit$coef[[1]],fit$coef[[3]])
}

# richness ~ MRCA + DR_harm + DR_skew + DR_cv
form=(richness ~ MRCA + DR_harm + DR_skew + DR_cv)
fit<-compar.gee(form, data=ordData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
ss<-as.character(form)

plot(log(richness) ~ MRCA, data=ordData, cex.lab=1.5, main=paste("27 orders; ",ss[2],ss[1],ss[3], "; P = ", round(dd[9],digits=3),sep=""))
if (dd[9] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}
plot(log(richness) ~ DR_harm, data=ordData, cex.lab=1.5, main=paste("27 orders; ",ss[2],ss[1],ss[3], "; P = ", round(dd[10],digits=3),sep=""))
if (dd[10] < 0.05){ abline(fit$coef[[1]],fit$coef[[3]])
}
plot(log(richness) ~ DR_skew, data=ordData, cex.lab=1.5, main=paste("27 orders; ",ss[2],ss[1],ss[3], "; P = ", round(dd[11],digits=3),sep=""))
if (dd[11] < 0.05){ abline(fit$coef[[1]],fit$coef[[3]])
}
plot(log(richness) ~ DR_kurt, data=ordData, cex.lab=1.5, main=paste("27 orders; ",ss[2],ss[1],ss[3], "; P = ", round(dd[12],digits=3),sep=""))
if (dd[12] < 0.05){ abline(fit$coef[[1]],fit$coef[[3]])
}

dev.off()

#########




# PHYLOGENETIC
# richness ~ MRCA
form=(richness ~ MRCA)
fit<-compar.gee(form, data=famData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
sink("cladeLevelFBD_GLMs_byFAM_RESULTS.txt", append=TRUE)
fit
dd
sink()
ss<-as.character(form)

pdf(file="cladeLevelFBD_GLMs_byFAM.pdf", onefile=TRUE, width=5, heigh=5)

plot(log(richness) ~ MRCA, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[3],digits=3),sep=""))
if (dd[3] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}
plot(residuals(fit)~MRCA,data=famData)


# richness ~ DR_harm
form=(richness ~ DR_harm)
fit<-compar.gee(form, data=famData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
sink("cladeLevelFBD_GLMs_byFAM_RESULTS.txt", append=TRUE)
fit
dd
sink()
ss<-as.character(form)

plot(log(richness) ~ DR_harm, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[3],digits=3),sep=""))
if (dd[3] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}
plot(residuals(fit)~DR_harm,data=famData)

# richness ~ DR_skew
form=(richness ~ DR_skew)
fit<-compar.gee(form, data=famData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
sink("cladeLevelFBD_GLMs_byFAM_RESULTS.txt", append=TRUE)
fit
dd
sink()
ss<-as.character(form)

plot(log(richness) ~ DR_skew, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[3],digits=3),sep=""))
if (dd[3] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}
plot(residuals(fit)~DR_skew,data=famData)


# richness ~ MRCA + DR_harm
form=(richness ~ MRCA + DR_harm)
fit<-compar.gee(form, data=famData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
sink("cladeLevelFBD_GLMs_byFAM_RESULTS.txt", append=TRUE)
fit
dd
sink()
ss<-as.character(form)

plot(log(richness) ~ MRCA, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[5],digits=3),sep=""))
if (dd[5] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}
plot(log(richness) ~ DR_harm, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[6],digits=3),sep=""))
if (dd[6] < 0.05){ abline(fit$coef[[1]],fit$coef[[3]])
}
plot(residuals(fit)~MRCA,data=famData)
plot(residuals(fit)~DR_harm,data=famData)


# richness ~ MRCA + DR_skew
form=(richness ~ MRCA + DR_skew)
fit<-compar.gee(form, data=famData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
sink("cladeLevelFBD_GLMs_byFAM_RESULTS.txt", append=TRUE)
fit
dd
sink()
ss<-as.character(form)

plot(log(richness) ~ MRCA, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[5],digits=3),sep=""))
if (dd[5] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}
plot(log(richness) ~ DR_skew, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[6],digits=3),sep=""))
if (dd[6] < 0.05){ abline(fit$coef[[1]],fit$coef[[3]])
}
plot(residuals(fit)~MRCA,data=famData)
plot(residuals(fit)~DR_skew,data=famData)

dev.off()









form=(richness ~ MRCA)
fit<-compar.gee(form, data=famData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
sink("cladeLevelFBD_GLMs_byFAM_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

pdf(file="cladeLevelFBD_GLMs_byFAM.pdf", onefile=TRUE, width=5, heigh=5)

plot(log(richness) ~ MRCA, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[3],digits=3),sep=""))
if (dd[3] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}

# richness ~ DR_harm
form=(richness ~ DR_harm)
fit<-compar.gee(form, data=famData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
sink("cladeLevelFBD_GLMs_byFAM_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

plot(log(richness) ~ DR_harm, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[3],digits=3),sep=""))
if (dd[3] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}

# richness ~ DR_skew
form=(richness ~ DR_skew)
fit<-compar.gee(form, data=famData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
sink("cladeLevelFBD_GLMs_byFAM_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

plot(log(richness) ~ DR_skew, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[3],digits=3),sep=""))
if (dd[3] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}


# richness ~ MRCA + DR_harm
form=(richness ~ MRCA + DR_harm)
fit<-compar.gee(form, data=famData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
sink("cladeLevelFBD_GLMs_byFAM_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

plot(log(richness) ~ MRCA, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[5],digits=3),sep=""))
if (dd[5] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}
plot(log(richness) ~ DR_harm, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[6],digits=3),sep=""))
if (dd[6] < 0.05){ abline(fit$coef[[1]],fit$coef[[3]])
}

# richness ~ MRCA + DR_skew
form=(richness ~ MRCA + DR_skew)
fit<-compar.gee(form, data=famData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
sink("cladeLevelFBD_GLMs_byFAM_RESULTS.txt", append=TRUE)
dd
sink()
ss<-as.character(form)

plot(log(richness) ~ MRCA, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[5],digits=3),sep=""))
if (dd[5] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}
plot(log(richness) ~ DR_skew, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[6],digits=3),sep=""))
if (dd[6] < 0.05){ abline(fit$coef[[1]],fit$coef[[3]])
}

dev.off()





# richness ~ MRCA
form=(richness ~ MRCA)
fit<-compar.gee(form, data=famData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
ss<-as.character(form)


plot(log(richness) ~ MRCA, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[3],digits=3),sep=""))
if (dd[3] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}


# richness ~ MRCA + DR_harm
form=(richness ~ MRCA + DR_harm)
fit<-compar.gee(form, data=famData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
ss<-as.character(form)

plot(log(richness) ~ MRCA, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[5],digits=3),sep=""))
if (dd[5] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}
plot(log(richness) ~ DR_harm, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[6],digits=3),sep=""))
if (dd[6] < 0.05){ abline(fit$coef[[1]],fit$coef[[3]])
}


# richness ~ MRCA + DR_harm + DR_skew
form=(richness ~ MRCA + DR_harm + DR_skew)
fit<-compar.gee(form, data=famData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
ss<-as.character(form)

plot(log(richness) ~ MRCA, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[7],digits=3),sep=""))
if (dd[7] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}
plot(log(richness) ~ DR_harm, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[8],digits=3),sep=""))
if (dd[8] < 0.05){ abline(fit$coef[[1]],fit$coef[[3]])
}
plot(log(richness) ~ DR_skew, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[9],digits=3),sep=""))
if (dd[9] < 0.05){ abline(fit$coef[[1]],fit$coef[[3]])
}

# richness ~ MRCA + DR_harm + DR_skew + DR_cv
form=(richness ~ MRCA + DR_harm + DR_skew + DR_cv)
fit<-compar.gee(form, data=famData, family="poisson", phy=treedata$phy)

dd<-drop1(fit)
ss<-as.character(form)

plot(log(richness) ~ MRCA, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[9],digits=3),sep=""))
if (dd[9] < 0.05){ abline(fit$coef[[1]],fit$coef[[2]])
}
plot(log(richness) ~ DR_harm, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[10],digits=3),sep=""))
if (dd[10] < 0.05){ abline(fit$coef[[1]],fit$coef[[3]])
}
plot(log(richness) ~ DR_skew, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[11],digits=3),sep=""))
if (dd[11] < 0.05){ abline(fit$coef[[1]],fit$coef[[3]])
}
plot(log(richness) ~ DR_kurt, data=famData, cex.lab=1.5, main=paste("114 fams; ",ss[2],ss[1],ss[3], "; P = ", round(dd[12],digits=3),sep=""))
if (dd[12] < 0.05){ abline(fit$coef[[1]],fit$coef[[3]])
}

dev.off()






#########
###
# cool-- now plot with PGLS these comparisons...
library(ape)
library(geiger)
library(phytools)
library(gee)
library(nlme)

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")

ordRes<-read.table("MamPhy_5911sp_FBD_cladeLevel_ORDS_skewKurt_withMRCAs.txt")
famRes<-read.table("MamPhy_5911sp_FBD_cladeLevel_FAMS_skewKurt_withMRCAs.txt")

ordNames<-rownames(ordRes[with(ordRes,order(ordRes$richness, decreasing=TRUE)),]) # is the ORDER names by category (skewness)
famNames<-rownames(famRes[with(famRes,order(famRes$richness, decreasing=TRUE)),]) # is the ORDER names by category (skewness)

ordPhy<-read.tree("MamPhy_BDvr_pcsFIXED_FBD_MCC_target_27ORDERS.tre")
famPhy<-read.tree("MamPhy_BDvr_pcsFIXED_FBD_MCC_target_162FAMILIES.tre")

res2<-ordRes

attach(res2)
main="27 orders of mammals"
pdf(file="cladeLevelFBD_DRcompare_ORDS_skewKurt_PGLS.pdf", onefile=TRUE, width=5, heigh=5)
# DR_cv by DR_harm
plot(DR_harm, DR_cv, cex.lab=1.5, main=main)

data2<-data.frame(DR_harm,DR_cv)
rownames(data2)<-ordNames
data3<-na.omit(data2)	
treedata<-treedata(ordPhy,data3)

pglsModel<-gls(DR_cv ~ DR_harm,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_skew by DR_harm
plot(DR_harm, DR_skew, cex.lab=1.5, main=main)

data2<-data.frame(DR_harm,DR_skew)
rownames(data2)<-ordNames
data3<-na.omit(data2)	
treedata<-treedata(ordPhy,data3)

pglsModel<-gls(DR_skew ~ DR_harm,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_kurt by DR_harm
plot(DR_harm, DR_kurt, cex.lab=1.5, main=main)

data2<-data.frame(DR_harm,DR_kurt)
rownames(data2)<-ordNames
data3<-na.omit(data2)	
treedata<-treedata(ordPhy,data3)

pglsModel<-gls(DR_kurt ~ DR_harm,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_harm by ln(richness)
plot(log(richness), DR_harm, cex.lab=1.5, main=main)

data2<-data.frame(log(richness),DR_harm)
rownames(data2)<-ordNames
data3<-na.omit(data2)	
treedata<-treedata(ordPhy,data3)

pglsModel<-gls(DR_harm ~ log.richness.,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# ln(richness) by DR_harm
plot(DR_harm, log(richness), cex.lab=1.5, main=main)

data2<-data.frame(log(richness),DR_harm)
rownames(data2)<-ordNames
data3<-na.omit(data2)	
treedata<-treedata(ordPhy,data3)

pglsModel<-gls(log.richness. ~ DR_harm,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))


# DR_harm by MRCA
plot(MRCA, DR_harm, cex.lab=1.5, main=main)

data2<-data.frame(MRCA,DR_harm)
rownames(data2)<-ordNames
data3<-na.omit(data2)	
treedata<-treedata(ordPhy,data3)

pglsModel<-gls(DR_harm ~ MRCA,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_cv by ln(richness)
plot(log(richness), DR_cv, cex.lab=1.5, main=main)

data2<-data.frame(log(richness),DR_cv)
rownames(data2)<-ordNames
data3<-na.omit(data2)	
treedata<-treedata(ordPhy,data3)

pglsModel<-gls(DR_cv ~ log.richness.,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_cv by MRCA
plot(MRCA, DR_cv, cex.lab=1.5, main=main)

data2<-data.frame(MRCA,DR_cv)
rownames(data2)<-ordNames
data3<-na.omit(data2)	
treedata<-treedata(ordPhy,data3)

pglsModel<-gls(DR_cv ~ MRCA,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_skew by ln richness
plot(log(richness), DR_skew, cex.lab=1.5, main=main)

data2<-data.frame(log(richness),DR_skew)
rownames(data2)<-ordNames
data3<-na.omit(data2)	
treedata<-treedata(ordPhy,data3)

pglsModel<-gls(DR_skew ~ log.richness.,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_skew by MRCA
plot(MRCA, DR_skew, cex.lab=1.5, main=main)

data2<-data.frame(MRCA,DR_skew)
rownames(data2)<-ordNames
data3<-na.omit(data2)	
treedata<-treedata(ordPhy,data3)

pglsModel<-gls(DR_skew ~ MRCA,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_kurt by ln richness
plot(log(richness), DR_kurt, cex.lab=1.5, main=main)

data2<-data.frame(log(richness),DR_kurt)
rownames(data2)<-ordNames
data3<-na.omit(data2)	
treedata<-treedata(ordPhy,data3)

pglsModel<-gls(DR_kurt ~ log.richness.,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_kurt by MRCA
plot(MRCA, DR_kurt, cex.lab=1.5, main=main)

data2<-data.frame(MRCA,DR_kurt)
rownames(data2)<-ordNames
data3<-na.omit(data2)	
treedata<-treedata(ordPhy,data3)

pglsModel<-gls(DR_kurt ~ MRCA,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_kurt by DR_skew
plot(DR_skew, DR_kurt, cex.lab=1.5, main=main)

data2<-data.frame(DR_skew,DR_kurt)
rownames(data2)<-ordNames
data3<-na.omit(data2)	
treedata<-treedata(ordPhy,data3)

pglsModel<-gls(DR_kurt ~ DR_skew,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# ln richness by MRCA
plot(MRCA, log(richness), cex.lab=1.5, main=main)

data2<-data.frame(MRCA,log(richness))
rownames(data2)<-ordNames
data3<-na.omit(data2)	
treedata<-treedata(ordPhy,data3)

pglsModel<-gls(log.richness. ~ MRCA,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

dev.off()

######
## And do the same thing for FAMILIES:
## with PGLS:

famData<-na.omit(famRes)	

res5<-famData

attach(res5)
main="162 families of mammals"
pdf(file="cladeLevelFBD_DRcompare_FAMS_skewKurt_PGLS.pdf", onefile=TRUE, width=5, height=5)
# DR_cv by DR_harm
plot(DR_harm, DR_cv, cex.lab=1.5, main=main)

data2<-data.frame(DR_harm,DR_cv)
rownames(data2)<-rownames(famData)
data3<-na.omit(data2)	
treedata<-treedata(famPhy,data3)

pglsModel<-gls(DR_cv ~ DR_harm,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_skew by DR_harm
plot(DR_harm, DR_skew, cex.lab=1.5, main=main)

data2<-data.frame(DR_harm,DR_skew)
rownames(data2)<-rownames(famData)
data3<-na.omit(data2)	
treedata<-treedata(famPhy,data3)

pglsModel<-gls(DR_skew ~ DR_harm,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_kurt by DR_harm
plot(DR_harm, DR_kurt, cex.lab=1.5, main=main)

data2<-data.frame(DR_harm,DR_kurt)
rownames(data2)<-rownames(famData)
data3<-na.omit(data2)	
treedata<-treedata(famPhy,data3)

pglsModel<-gls(DR_kurt ~ DR_harm,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_harm by ln(richness)
plot(log(richness), DR_harm, cex.lab=1.5, main=main)

data2<-data.frame(log(richness),DR_harm)
rownames(data2)<-rownames(famData)
data3<-na.omit(data2)	
treedata<-treedata(famPhy,data3)

pglsModel<-gls(DR_harm ~ log.richness.,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# ln(richness) by DR_harm
plot(DR_harm, log(richness), cex.lab=1.5, main=main)

data2<-data.frame(log(richness),DR_harm)
rownames(data2)<-rownames(famData)
data3<-na.omit(data2)	
treedata<-treedata(famPhy,data3)

pglsModel<-gls(log.richness. ~ DR_harm,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))


# DR_harm by MRCA
plot(MRCA, DR_harm, cex.lab=1.5, main=main)

data2<-data.frame(MRCA,DR_harm)
rownames(data2)<-rownames(famData)
data3<-na.omit(data2)	
treedata<-treedata(famPhy,data3)

pglsModel<-gls(DR_harm ~ MRCA,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_cv by ln(richness)
plot(log(richness), DR_cv, cex.lab=1.5, main=main)

data2<-data.frame(log(richness),DR_cv)
rownames(data2)<-rownames(famData)
data3<-na.omit(data2)	
treedata<-treedata(famPhy,data3)

pglsModel<-gls(DR_cv ~ log.richness.,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_cv by MRCA
plot(MRCA, DR_cv, cex.lab=1.5, main=main)

data2<-data.frame(MRCA,DR_cv)
rownames(data2)<-rownames(famData)
data3<-na.omit(data2)	
treedata<-treedata(famPhy,data3)

pglsModel<-gls(DR_cv ~ MRCA,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_skew by ln richness
plot(log(richness), DR_skew, cex.lab=1.5, main=main)

data2<-data.frame(log(richness),DR_skew)
rownames(data2)<-rownames(famData)
data3<-na.omit(data2)	
treedata<-treedata(famPhy,data3)

pglsModel<-gls(DR_skew ~ log.richness.,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_skew by MRCA
plot(MRCA, DR_skew, cex.lab=1.5, main=main)

data2<-data.frame(MRCA,DR_skew)
rownames(data2)<-rownames(famData)
data3<-na.omit(data2)	
treedata<-treedata(famPhy,data3)

pglsModel<-gls(DR_skew ~ MRCA,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_kurt by ln richness
plot(log(richness), DR_kurt, cex.lab=1.5, main=main)

data2<-data.frame(log(richness),DR_kurt)
rownames(data2)<-rownames(famData)
data3<-na.omit(data2)	
treedata<-treedata(famPhy,data3)

pglsModel<-gls(DR_kurt ~ log.richness.,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_kurt by MRCA
plot(MRCA, DR_kurt, cex.lab=1.5, main=main)

data2<-data.frame(MRCA,DR_kurt)
rownames(data2)<-rownames(famData)
data3<-na.omit(data2)	
treedata<-treedata(famPhy,data3)

pglsModel<-gls(DR_kurt ~ MRCA,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# DR_kurt by DR_skew
plot(DR_skew, DR_kurt, cex.lab=1.5, main=main)

data2<-data.frame(DR_skew,DR_kurt)
rownames(data2)<-rownames(famData)
data3<-na.omit(data2)	
treedata<-treedata(famPhy,data3)

pglsModel<-gls(DR_kurt ~ DR_skew,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

# ln richness by MRCA
plot(MRCA, log(richness), cex.lab=1.5, main=main)

data2<-data.frame(MRCA,log(richness))
rownames(data2)<-rownames(famData)
data3<-na.omit(data2)	
treedata<-treedata(famPhy,data3)

pglsModel<-gls(log.richness. ~ MRCA,correlation=corBrownian(phy=treedata$phy), data=data3, method="ML")
sum<-summary(pglsModel)
if (sum$tTable[8] < 0.05){
	abline(pglsModel$coefficients[[1]], pglsModel$coefficients[[2]], lty=2)
}
mtext(paste("PGLS p-value = ",round(sum$tTable[8],4)))

dev.off()




###
## with SPEARMAN:

res5<-read.table("MamPhy_5910sp_Exp_cladeLevel_FAMS_skewKurt_withMRCAs.txt", header=TRUE)

attach(res5)
main="162 families of mammals"
pdf(file="cladeLevel_DRcompare_FAMS.pdf", onefile=TRUE)
# DR_cv by DR_harm
plot(DR_harm, DR_cv, cex.lab=1.5, main=main)
xvar=DR_harm
yvar=DR_cv
corr<-cor.test(xvar, yvar, method="spearman")
lm<-aov(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_skew by DR_harm
plot(DR_harm, DR_skew, cex.lab=1.5, main=main)
xvar=DR_harm
yvar=DR_skew
corr<-cor.test(xvar, yvar, method="spearman")
lm<-aov(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_kurt by DR_harm
plot(DR_harm, DR_kurt, cex.lab=1.5, main=main)
xvar=DR_harm
yvar=DR_kurt
corr<-cor.test(xvar, yvar, method="spearman")
lm<-aov(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))


# DR_harm by richness
plot(richness, DR_harm, cex.lab=1.5, main=main)
xvar= richness
yvar= DR_harm
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_harm by MRCA
plot(MRCA, DR_harm, cex.lab=1.5, main=main)
xvar= MRCA
yvar= DR_harm
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_cv by richness
plot(richness, DR_cv, cex.lab=1.5, main=main)
xvar= richness
yvar= DR_cv
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_cv by MRCA
plot(MRCA, DR_cv, cex.lab=1.5, main=main)
xvar=MRCA
yvar=DR_cv
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_skew by richness
plot(richness, DR_skew, cex.lab=1.5, main=main)
xvar= richness
yvar= DR_skew
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_skew by MRCA
plot(MRCA, DR_skew, cex.lab=1.5, main=main)
xvar=MRCA
yvar=DR_skew
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_kurt by richness
plot(richness, DR_kurt, cex.lab=1.5, main=main)
xvar= richness
yvar= DR_kurt
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_kurt by MRCA
plot(MRCA, DR_kurt, cex.lab=1.5, main=main)
xvar=MRCA
yvar=DR_kurt
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_kurt by DR_skew
plot(DR_skew, DR_kurt, cex.lab=1.5, main=main)
xvar=DR_skew
yvar=DR_kurt
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# richness by MRCA
plot(MRCA, richness, cex.lab=1.5, main=main)
xvar=MRCA
yvar=richness
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

dev.off()












########
# BELOW is with SPEARMAN, not PGLS...

###
# cool-- now PLOT these comparisons...

res2<-read.table("MamPhy_5910sp_Exp_cladeLevel_ORDS_skewKurt_withMRCAs.txt", header=TRUE)

attach(res2)
main="27 orders of mammals"
pdf(file="cladeLevel_DRcompare_ORDS_skewKurt.pdf", onefile=TRUE)
# DR_cv by DR_harm
plot(DR_harm, DR_cv, cex.lab=1.5, main=main)
xvar=DR_harm
yvar=DR_cv
corr<-cor.test(xvar, yvar, method="spearman")
lm<-aov(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_skew by DR_harm
plot(DR_harm, DR_skew, cex.lab=1.5, main=main)
xvar=DR_harm
yvar=DR_skew
corr<-cor.test(xvar, yvar, method="spearman")
lm<-aov(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_kurt by DR_harm
plot(DR_harm, DR_kurt, cex.lab=1.5, main=main)
xvar=DR_harm
yvar=DR_kurt
corr<-cor.test(xvar, yvar, method="spearman")
lm<-aov(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))


# DR_harm by richness
plot(richness, DR_harm, cex.lab=1.5, main=main)
xvar= richness
yvar= DR_harm
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_harm by MRCA
plot(MRCA, DR_harm, cex.lab=1.5, main=main)
xvar= MRCA
yvar= DR_harm
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_cv by richness
plot(richness, DR_cv, cex.lab=1.5, main=main)
xvar= richness
yvar= DR_cv
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_cv by MRCA
plot(MRCA, DR_cv, cex.lab=1.5, main=main)
xvar=MRCA
yvar=DR_cv
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_skew by richness
plot(richness, DR_skew, cex.lab=1.5, main=main)
xvar= richness
yvar= DR_skew
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_skew by MRCA
plot(MRCA, DR_skew, cex.lab=1.5, main=main)
xvar=MRCA
yvar=DR_skew
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_kurt by richness
plot(richness, DR_kurt, cex.lab=1.5, main=main)
xvar= richness
yvar= DR_kurt
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_kurt by MRCA
plot(MRCA, DR_kurt, cex.lab=1.5, main=main)
xvar=MRCA
yvar=DR_kurt
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_kurt by DR_skew
plot(DR_skew, DR_kurt, cex.lab=1.5, main=main)
xvar=DR_skew
yvar=DR_kurt
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))


# richness by MRCA
plot(MRCA, richness, cex.lab=1.5, main=main)
xvar=MRCA
yvar=richness
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

dev.off()

######
## And do the same thing for FAMILIES:

##
famsGrThan1<-names(table(cladesDR$fam))

DR_harm<-data.frame(matrix(NA, nrow = length(famsGrThan1), ncol = 1), row.names=famsGrThan1)
DR_cv<-data.frame(matrix(NA, nrow = length(famsGrThan1), ncol = 1), row.names=famsGrThan1)
DR_skew<-data.frame(matrix(NA, nrow = length(famsGrThan1), ncol = 1), row.names=famsGrThan1)
DR_kurt<-data.frame(matrix(NA, nrow = length(famsGrThan1), ncol = 1), row.names=famsGrThan1)
richness<-data.frame(matrix(NA, nrow = length(famsGrThan1), ncol = 1), row.names=famsGrThan1)

for(i in 1:length(famsGrThan1)){
	DR_harm[i,] <- mean(cladesDR[which(cladesDR$fam==famsGrThan1[i]),"harmMeans"])
	DR_cv[i,] <- mean(cladesDR[which(cladesDR$fam==famsGrThan1[i]),"cv"]*100)
	DR_skew[i,] <- skewness(cladesDR[which(cladesDR$fam==famsGrThan1[i]),"harmMeans"])
	DR_kurt[i,] <- kurtosis(cladesDR[which(cladesDR$fam==famsGrThan1[i]),"harmMeans"])
	richness[i,] <- length(cladesDR[which(cladesDR$fam==famsGrThan1[i]),"cv"])
}

res3<-cbind(DR_harm, DR_cv, DR_skew, DR_kurt, richness)

colnames(res3)<-c("DR_harm","DR_cv", "DR_skew", "DR_kurt", "richness")


MRCA<-data.frame(matrix(NA, nrow = length(famsGrThan1), ncol = 1), row.names=famsGrThan1)

famsGrThan1 <- rownames(res3[which(res3$richness!=1),])
famsEq1 <- rownames(res3[which(res3$richness==1),])

for(i in 1:length(famsGrThan1)){
	node <- getMRCA(mamPhy, as.vector(cladesDR[which(cladesDR$fam==famsGrThan1[i]),"tiplabel"]))
	MRCA[match(famsGrThan1[i], rownames(MRCA)),] <- mamPhy$height[node-5911]
}

res4<-cbind(DR_harm, DR_cv, DR_skew, DR_kurt, richness, MRCA)
colnames(res4)<-c("DR_harm","DR_cv", "DR_skew", "DR_kurt", "richness","MRCA")

write.table(res4,"MamPhy_5910sp_Exp_cladeLevel_FAMS_skewKurt_withMRCAs.txt")

### awesome.
## now PLOT these comparisons too...

res5<-read.table("MamPhy_5910sp_Exp_cladeLevel_FAMS_skewKurt_withMRCAs.txt", header=TRUE)

attach(res5)
main="162 families of mammals"
pdf(file="cladeLevel_DRcompare_FAMS.pdf", onefile=TRUE)
# DR_cv by DR_harm
plot(DR_harm, DR_cv, cex.lab=1.5, main=main)
xvar=DR_harm
yvar=DR_cv
corr<-cor.test(xvar, yvar, method="spearman")
lm<-aov(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_skew by DR_harm
plot(DR_harm, DR_skew, cex.lab=1.5, main=main)
xvar=DR_harm
yvar=DR_skew
corr<-cor.test(xvar, yvar, method="spearman")
lm<-aov(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_kurt by DR_harm
plot(DR_harm, DR_kurt, cex.lab=1.5, main=main)
xvar=DR_harm
yvar=DR_kurt
corr<-cor.test(xvar, yvar, method="spearman")
lm<-aov(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))


# DR_harm by richness
plot(richness, DR_harm, cex.lab=1.5, main=main)
xvar= richness
yvar= DR_harm
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_harm by MRCA
plot(MRCA, DR_harm, cex.lab=1.5, main=main)
xvar= MRCA
yvar= DR_harm
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_cv by richness
plot(richness, DR_cv, cex.lab=1.5, main=main)
xvar= richness
yvar= DR_cv
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_cv by MRCA
plot(MRCA, DR_cv, cex.lab=1.5, main=main)
xvar=MRCA
yvar=DR_cv
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_skew by richness
plot(richness, DR_skew, cex.lab=1.5, main=main)
xvar= richness
yvar= DR_skew
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_skew by MRCA
plot(MRCA, DR_skew, cex.lab=1.5, main=main)
xvar=MRCA
yvar=DR_skew
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_kurt by richness
plot(richness, DR_kurt, cex.lab=1.5, main=main)
xvar= richness
yvar= DR_kurt
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_kurt by MRCA
plot(MRCA, DR_kurt, cex.lab=1.5, main=main)
xvar=MRCA
yvar=DR_kurt
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# DR_kurt by DR_skew
plot(DR_skew, DR_kurt, cex.lab=1.5, main=main)
xvar=DR_skew
yvar=DR_kurt
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

# richness by MRCA
plot(MRCA, richness, cex.lab=1.5, main=main)
xvar=MRCA
yvar=richness
corr<-cor.test(xvar, yvar, method="spearman")
lm<-lm(yvar ~ xvar)
abline(lm$coefficients[[1]],lm$coefficients[[2]], lty=2)
mtext(paste("rho p-value = ",round(corr$p.value,4)))

dev.off()





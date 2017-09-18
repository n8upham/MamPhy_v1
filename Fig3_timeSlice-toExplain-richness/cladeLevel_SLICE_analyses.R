
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

bbone<- "NDexp" #"FBD" #

mamPhy<-drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target.tre",sep="")), "_Anolis_carolinensis") ## Use the MCC target node heights one... 

root=max(nodeHeights(mamPhy)[,2])

trees5Ma<-treeSlice(mamPhy, slice=root-5, trivial=FALSE, prompt=FALSE)

trees10Ma<-treeSlice(mamPhy, slice=root-10, trivial=FALSE, prompt=FALSE)

trees15Ma<-treeSlice(mamPhy, slice=root-15, trivial=FALSE, prompt=FALSE)

trees20Ma<-treeSlice(mamPhy, slice=root-20, trivial=FALSE, prompt=FALSE)

trees25Ma<-treeSlice(mamPhy, slice=root-25, trivial=FALSE, prompt=FALSE)

trees30Ma<-treeSlice(mamPhy, slice=root-30, trivial=FALSE, prompt=FALSE)

trees35Ma<-treeSlice(mamPhy, slice=root-35, trivial=FALSE, prompt=FALSE)

trees40Ma<-treeSlice(mamPhy, slice=root-40, trivial=FALSE, prompt=FALSE)

trees45Ma<-treeSlice(mamPhy, slice=root-45, trivial=FALSE, prompt=FALSE)

trees50Ma<-treeSlice(mamPhy, slice=root-50, trivial=FALSE, prompt=FALSE)

trees55Ma<-treeSlice(mamPhy, slice=root-55, trivial=FALSE, prompt=FALSE)

trees60Ma<-treeSlice(mamPhy, slice=root-60, trivial=FALSE, prompt=FALSE)

trees65Ma<-treeSlice(mamPhy, slice=root-65, trivial=FALSE, prompt=FALSE)

trees70Ma<-treeSlice(mamPhy, slice=root-70, trivial=FALSE, prompt=FALSE)

allCladeSets<-list(trees5Ma,trees10Ma,trees15Ma,trees20Ma,trees25Ma,trees30Ma,trees35Ma,trees40Ma,trees45Ma,trees50Ma,trees55Ma,trees60Ma,trees65Ma,trees70Ma)
allCladeSetNames<-c("5Ma","10Ma","15Ma","20Ma","25Ma","30Ma","35Ma","40Ma","45Ma","50Ma","55Ma","60Ma","65Ma","70Ma")

lengths<-vector()
for (i in 1:length(allCladeSets)){
	lengths[i]<-length(allCladeSets[[i]])
}
lengths
names(lengths)<-allCladeSetNames 

# FBD
#  5Ma 10Ma 15Ma 20Ma 25Ma 30Ma 35Ma 40Ma 45Ma 50Ma 55Ma 60Ma 65Ma 70Ma 
# 1060  545  306  196  136   99   77   65   50   43   36   30   26   16 

# ND
# 5Ma 10Ma 15Ma 20Ma 25Ma 30Ma 35Ma 40Ma 45Ma 50Ma 55Ma 60Ma 65Ma 70Ma 
# 1166  656  377  241  159  123   95   75   61   45   38   31   24   17 

for (i in 1:length(allCladeSets)){
	trees<-allCladeSets[[i]]
	write.tree(trees,file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target_slice",allCladeSetNames[[i]],"_newick.trees",sep=""))
	}

# load back in
allCladeSets<-vector("list",length(allCladeSetNames))
for (i in 1:length(allCladeSetNames)){
	allCladeSets[[i]]<-read.tree(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target_slice",allCladeSetNames[[i]],"_newick.trees",sep=""))
	}

#cc<-timeSliceTree(famPhy, sliceTime=50, drop.extinct = FALSE, plot = TRUE)
#	# Does the REVERSE of what I'm interested in ! Still interesting later maybe...
#cutFams10<-timeslice.phy(famPhy,slice.time=10,plot=TRUE)
#	#old BABST function... actually I think this does the CORRECT thing-- same as PHYTOOLS-- just distance above the root.

# Now want to summarize the PER CLADE values for each of these clades...
#[1] "DR_harm"  "DR_cv"    "DR_skew"  "DR_kurt"  "richness" "MRCA"    

cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdGenesSampPC_DR-SUMMARY-expanded_",bbone,".txt",sep=""))
colnames(cladesDR)<-c("tiplabel","gen","fam","ord","clade","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")
head(cladesDR)

#######
## Make those per-clade summaries:

# Next need to expand this to a POSTERIOR sample of trees ~ 100 ...

# Have to get the TIPLABELS for each of the taxa in EACH tree of the time slice...
# use mamPhy as above
btimes<-branching.times(mamPhy)

##
# LOOPING across all these SLICES
##
for(j in 1:length(allCladeSets)){
cladeSet<-allCladeSets[[j]]

	DR_harm<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	DR_cv<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	DR_skew<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	DR_kurt<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	percentSamp<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	richness<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	MRCA_mcc<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))

	for (i in 1:length(cladeSet)){
	cladeSp<-cladeSet[[i]]$tip.label
		DR_harm[i,] <- mean(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"])
		DR_cv[i,] <- mean(cladesDR[match(cladeSp,cladesDR$tiplabel),"cv"]*100)
		DR_skew[i,] <- skewness(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"])
		DR_kurt[i,] <- kurtosis(cladesDR[match(cladeSp,cladesDR$tiplabel),"harmMeans"])
		percentSamp[i,] <- length(which(cladesDR[match(cladeSp,cladesDR$tiplabel),"samp"]=="sampled"))/length(cladeSp)
		richness[i,] <- length(cladeSp)
		node <- getMRCA(mamPhy, cladeSp)
		MRCA_mcc[i,] <- btimes[node-5911] #taking the MCC height
	}

	res<-cbind(DR_harm, DR_cv, DR_skew, DR_kurt, percentSamp, richness, MRCA_mcc)

	colnames(res)<-c("DR_harm","DR_cv", "DR_skew", "DR_kurt", "percentSamp", "richness", "MRCA_mcc")

	write.table(res,paste("MamPhy_5911sp_",bbone,"_cladeLevel_DRstats_trees", allCladeSetNames[j], ".txt", sep=""))

}

## Now LOAD them back in !

results<-vector("list",length(allCladeSets))
for (i in 1: length(allCladeSets)){
	res<-read.table(paste("MamPhy_5911sp_",bbone,"_cladeLevel_DRstats_trees",allCladeSetNames[i],".txt",sep=""), header=TRUE)
	results[[i]]<-res
}


# Per-clade time slice histograms
##

results[[1]][match(order(results[[1]]$DR_skew),rownames(results[[1]])),] ## arranges results as ascending in SKEWNESS


for (k in 1:length(results)){
cladeNames<-order(results[[k]]$DR_skew)

pdf(file=paste("cladeLevel",bbone,"_DRcompare_DRhist_bySkew_trees", allCladeSetNames[[k]], ".pdf",sep=""), onefile=TRUE, width=5, height=5)

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

pdf(file=paste("cladeLevel",bbone,"_DRcompare_DRhist_bySkew_Xsame_trees", allCladeSetNames[[k]], ".pdf",sep=""), onefile=TRUE, width=5, height=5)

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

pdf(file=paste("cladeLevel",bbone,"_DRcompare_DRhist_bySkew_Xsame-DensMamMean_trees", allCladeSetNames[[k]], ".pdf",sep=""), onefile=TRUE, width=5, height=5)

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

pdf(file=paste("cladeLevel",bbone,"_DRcompare_DRhist_bySkew_Xsame-DensMamDist_trees", allCladeSetNames[[k]], ".pdf",sep=""), onefile=TRUE, width=5, height=5)

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

# FreqDists of RICHNESS for each TIMESLICE 

pdf(file=paste("cladeLevel",bbone,"_timeSlice_5Ma-to-70Ma_richnessHIST_combo.pdf",sep=""),onefile=TRUE,width=20,height=20)
#quartz(width=20,height=20)
layout(matrix(c(1:16), 4, 4, byrow = TRUE))

for(i in 1:length(results)){
	hist(results[[i]]$richness, breaks=50, ylab="Frequency",xlab="Clade richness",cex.main=1.2,main=paste("slice ",allCladeSetNames[[i]],sep=""))
}
dev.off()

pdf(file=paste("cladeLevel",bbone,"_timeSlice_5Ma-to-70Ma_richnessDENSITY_combo.pdf",sep=""),onefile=TRUE,width=20,height=20)
#quartz(width=20,height=20)
layout(matrix(c(1:16), 4, 4, byrow = TRUE))

for(i in 1:length(results)){
	dens.rate <- density(results[[i]]$richness)$y
	plot(density(results[[i]]$richness), col="light grey", yaxt="n", xlim=range(density(results[[14]]$richness)$x), ylim=range(dens.rate), ylab="Frequency",xlab="Clade richness",cex.main=1.2,main=paste("slice ",allCladeSetNames[[i]],sep=""))
	polygon(density(results[[i]]$richness), col="light grey", border="light grey", bty="n", xlab="")
	axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=1, las=1, lwd=1, cex.axis=1, tck=-0.05, mgp=c(1,1,0))
	abline(v=mean(results[[i]]$richness),lty=2,lwd=2)
}
dev.off()

pdf(file=paste("cladeLevel",bbone,"_timeSlice_5Ma-to-70Ma_meanDR-DENSITY_combo.pdf",sep=""),onefile=TRUE,width=20,height=20)
#quartz(width=20,height=20)
layout(matrix(c(1:16), 4, 4, byrow = TRUE))

for(i in 1:length(results)){
	dens.rate <- density(results[[i]]$DR_harm)$y
	plot(density(results[[i]]$DR_harm), col="light grey", yaxt="n", xlim=range(cladesDR$harmMeans), ylim=range(dens.rate), ylab="Frequency",xlab="Clade mean DR",cex.main=1.2,main=paste("slice ",allCladeSetNames[[i]],sep=""))
	polygon(density(results[[i]]$DR_harm), col="light grey", border="light grey", bty="n", xlab="")
	axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=1, las=1, lwd=1, cex.axis=1, tck=-0.05, mgp=c(1,1,0))
	abline(v=mean(results[[i]]$DR_harm),lty=2,lwd=2)
}
dev.off()


xlim=range(cladesDR$harmMeans), 

pdf(file=paste("cladeLevel",bbone,"_timeSlice_5Ma-to-70Ma_meanDR-HIST_combo.pdf",sep=""),onefile=TRUE,width=20,height=20)
#quartz(width=20,height=20)
layout(matrix(c(1:16), 4, 4, byrow = TRUE))

for(i in 1:length(results)){
	hist(results[[i]]$DR_harm, breaks=50, ylab="Frequency",xlab="Clade mean DR",cex.main=1.2,main=paste("slice ",allCladeSetNames[[i]],sep=""))
}
dev.off()


x.tick <- as.vector(quantile(cladesDR$harmMeans, c(0.01,0.5,0.99,1)))




#########

# GLM analyses::
##########
#####
# 
# load back in tables
results<-vector("list",length(allCladeSets))
for (i in 1: length(allCladeSets)){
	res<-read.table(paste("MamPhy_5911sp_",bbone,"_cladeLevel_DRstats_trees",allCladeSetNames[i],".txt",sep=""), header=TRUE)
	results[[i]]<-res
}

# load back in trees
allCladeSets<-vector("list",length(allCladeSetNames))
for (i in 1:length(allCladeSetNames)){
	allCladeSets[[i]]<-read.tree(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target_slice",allCladeSetNames[[i]],"_newick.trees",sep=""))
	}

##
# GLMS
# Loop through each of these for the TIME-SLICES
##

# TIMESLICES - start HERE
for (i in 1:length(results)){

cladeData<-na.omit(results[[i]])	

# richness ~ MRCA
form=(richness ~ MRCA_mcc)
fit<-glm(form, data=cladeData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink(file=paste("cladeLevel",bbone,"_GLMs_bySLICE_", allCladeSetNames[i],"_RESULTS.txt",sep=""), append=TRUE)
print(form)
print(dd)
sink()
ss<-as.character(form)

pdf(file=paste("cladeLevel",bbone,"_GLMs_bySLICE_", allCladeSetNames[i],"_withPartialResiduals.pdf",sep=""), onefile=TRUE, width=5, heigh=5)

termplot(fit, partial.resid=TRUE, rug=TRUE, cex.main=1, font.sub=2, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1]),"; ",ss[2],ss[1],ss[3], sep=""), sub=paste("X1=", round(dd$coef[[2]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ DR_harm
form=(richness ~ DR_harm)
fit<-glm(form, data=cladeData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink(file=paste("cladeLevel",bbone,"_GLMs_bySLICE_", allCladeSetNames[i],"_RESULTS.txt",sep=""), append=TRUE)
print(form)
print(dd)
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, rug=TRUE, cex.main=1, font.sub=2, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1]),"; ",ss[2],ss[1],ss[3], sep=""), sub=paste("X1=", round(dd$coef[[2]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ DR_skew
form=(richness ~ DR_skew)
fit<-glm(form, data=cladeData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink(file=paste("cladeLevel",bbone,"_GLMs_bySLICE_", allCladeSetNames[i],"_RESULTS.txt",sep=""), append=TRUE)
print(form)
print(dd)
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, rug=TRUE, cex.main=1, font.sub=2, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1]),"; ",ss[2],ss[1],ss[3], sep=""), sub=paste("X1=", round(dd$coef[[2]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ MRCA_mcc + DR_harm
form=(richness ~ MRCA_mcc + DR_harm)
fit<-glm(form, data=cladeData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink(file=paste("cladeLevel",bbone,"_GLMs_bySLICE_", allCladeSetNames[i],"_RESULTS.txt",sep=""), append=TRUE)
print(form)
print(dd)
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, rug=TRUE, cex.main=1, font.sub=2, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1]),"; ",ss[2],ss[1],ss[3], sep=""), sub=paste("X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))


# richness ~ MRCA_mcc + DR_skew
form=(richness ~ MRCA_mcc + DR_skew)
fit<-glm(form, data=cladeData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink(file=paste("cladeLevel",bbone,"_GLMs_bySLICE_", allCladeSetNames[i],"_RESULTS.txt",sep=""), append=TRUE)
print(form)
print(dd)
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, rug=TRUE, cex.main=1, font.sub=2, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1]),"; ",ss[2],ss[1],ss[3], sep=""), sub=paste("X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# richness ~ MRCA_mcc + DR_harm + DR_skew
form=(richness ~ MRCA_mcc + DR_harm + DR_skew)
fit<-glm(form, data=cladeData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink(file=paste("cladeLevel",bbone,"_GLMs_bySLICE_", allCladeSetNames[i],"_RESULTS.txt",sep=""), append=TRUE)
print(form)
print(dd)
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, rug=TRUE, cex.main=1, font.sub=2, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1]),"; ",ss[2],ss[1],ss[3], sep=""), sub=paste("X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; X3=", round(dd$coef[[4]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))


# richness ~ MRCA_mcc + DR_harm + DR_skew + DR_cv
form=(richness ~ MRCA_mcc + DR_harm + DR_skew + DR_cv)
fit<-glm(form, data=cladeData, family="poisson")#, phy=treedata$phy)

dd<-summary(fit)
sink(file=paste("cladeLevel",bbone,"_GLMs_bySLICE_", allCladeSetNames[i],"_RESULTS.txt",sep=""), append=TRUE)
print(form)
print(dd)
sink()
ss<-as.character(form)

termplot(fit, partial.resid=TRUE, rug=TRUE, cex.main=0.8, font.sub=2, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1]),"; ",ss[2],ss[1],ss[3], sep=""), sub=paste("X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; X3=", round(dd$coef[[4]],digits=3), "; X4=", round(dd$coef[[5]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

dev.off()

}

###
# with GLMS -- get the COEF from each GLMM and the PartialResid = effect size of some kind...
# example
i=14
	cladeData<-na.omit(results[[i]])	

	form=(richness ~ MRCA_mcc + DR_harm + DR_skew)
	fit<-glm(form, data=cladeData, family="poisson")#, phy=treedata$phy)
	dd<-summary(fit)
	ss<-as.character(form)

	parts<-residuals(fit,"partial")

	aa<-termplot(fit, partial.resid=TRUE, se=TRUE, plot=FALSE, rug=TRUE, cex.main=1, font.sub=2, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1]),"; ",ss[2],ss[1],ss[3], sep=""), sub=paste("X1=", round(dd$coef[[2]],digits=3), "; X2=", round(dd$coef[[3]],digits=3), "; X3=", round(dd$coef[[4]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

partSlopes<-vector()
for (i in 1:3){
	partEq<-lm(aa[[i]][[2]] ~ aa[[i]][[1]])
	partSlopes[i]<-partEq$coef[[2]]
}
# ^ Ah! these slopes are EQUIVALENT to the slopes out of the GLM fit... got it.  cool
# partial residuals are just a way to VISUALIZE the deviation from the slope of each term.

# OK, so just grab out all the dd coefficients.
# in the dd object.

partSlopes<-data.frame(matrix(NA, nrow = length(results), ncol = 4))
sliceTimes<-seq(-5,-70,-5)
partSlopesAll<-cbind(sliceTimes,partSlopes)
colnames(partSlopesAll)<-c("time","MRCA_mcc","DR_harm","DR_skew","AIC")

for (i in 1:length(results)){
	cladeData<-na.omit(results[[i]])	
	form=(richness ~ MRCA_mcc + DR_harm + DR_skew)
	fit<-glm(form, data=cladeData, family="poisson")#, phy=treedata$phy)
	dd<-summary(fit)
	partSlopesAll[i,2]<-round(dd$coef[[2]],digits=3)
	partSlopesAll[i,3]<-round(dd$coef[[3]],digits=3)
	partSlopesAll[i,4]<-round(dd$coef[[4]],digits=3)
	partSlopesAll[i,5]<-round(dd$aic,digits=0)
}
write.table(partSlopesAll,paste("cladeLevel",bbone,"_GLMs_timeSlice_5Ma-to-70Ma_partSlopes.txt",sep=""))

# plot through time
pdf(file=paste("cladeLevel",bbone,"_GLMs_timeSlice_5Ma-to-70Ma_partSlopes_PLOTTED_3part.pdf",sep=""),onefile=TRUE, width=10,height=4)

#quartz(width=10,height=4)
#layout(matrix(c(1:3), 3, 1, byrow = TRUE))
par(mfrow = c(1,3), oma = c(5,4,5,0) + 0.1, mar = c(1,1,1,1) + 0.1)
#par(op)

form=(MRCA_mcc ~ time)
y.lab=""
x.lab=""
plot(form, data=partSlopesAll, xlab=x.lab,ylab=y.lab)#, xaxt="n")
#axis(side=1,labels=FALSE)
mtext(text="(a) MRCA", adj=0, font=2)#at=c(-72,0.21)

fit2<-lm(form,data=partSlopesAll)
dd2<-summary(fit2)
if (dd2$coef[8] < 0.05){
	abline(fit2$coef[1],fit2$coef[2], lty=1)
	text(x=-72,y=0.215,pos=4,font=1, labels=bquote( plain("R") ^ plain("2") * plain(" = ") * .(round(dd2$r.squared, digits=3)) ))
	#text(x=-72,y=0.215,adj=c(-0.1,2),font=1, labels=bquote( plain("X = ") * .(round(fit2$coef[2],digits=3)) ))
}

form=(DR_harm ~ time)
plot(form, data=partSlopesAll, xlab=x.lab,ylab=y.lab)#, xaxt="n")
#axis(side=1,labels=FALSE)
mtext(text="(b) DR mean", adj=0, font=2)#at=c(-72,0.21)

fit2<-lm(form,data=partSlopesAll)
dd2<-summary(fit2)
if (dd2$coef[8] < 0.05){
	abline(fit2$coef[1],fit2$coef[2], lty=1)
	text(x=-18,y=11.3,pos=4,font=1, labels=bquote( plain("R") ^ plain("2") * plain(" = ") * .(round(dd2$r.squared, digits=3)) ))
	#text(x=-18,y=11.3,adj=c(-0.1,2),font=1, labels=bquote( plain("X = ") * .(round(fit2$coef[2],digits=3)) ))
}

form=(DR_skew ~ time)
plot(form, data=partSlopesAll, xlab=x.lab,ylab=y.lab)
mtext(text="(c) DR skew", adj=0, font=2)#at=c(-72,0.21)

fit2<-lm(form,data=partSlopesAll)
dd2<-summary(fit2)
if (dd2$coef[8] < 0.05){
	abline(fit2$coef[1],fit2$coef[2], lty=1)
	text(x=-18,y=0.9,pos=4,font=1, labels=bquote( plain("R") ^ plain("2") * plain(" = ") * .(round(dd2$r.squared, digits=3)) ))
	#text(x=-18,y=0.9,adj=c(-0.1,2),font=1, labels=bquote( plain("X = ") * .(round(fit2$coef[2],digits=3)) ))
}
title(main="GLM: poisson(richness) ~ MRCA_mcc + DR_harm + DR_skew", xlab = "Time slices before present (Ma)",
      ylab = "Slope of partial residual",
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)

dev.off()



###
##
# PHYLOGENETIC --- PGLS !!
###
# get the PHYs for each time slice set.
bbone<- "FBD"#"NDexp" # # 
bbone<- "NDexp" # # 

slicePhys<-vector("list",length(allCladeSets))
for (k in 1:length(allCladeSets)){

cladeReps<-vector()
for (i in 1:length(allCladeSets[[k]])){
	cladeSp<-allCladeSets[[k]][[i]]$tip.label
	cladeReps[i]<-cladeSp[1]
	}
toDrop<-setdiff(mamPhy$tip.label,cladeReps)
slicePhys[[k]]<-drop.tip(mamPhy,toDrop)
}

for(i in 1:length(slicePhys)){
	write.tree(slicePhys[[i]], file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target_slicePhy-5to70Ma_by5.trees",sep=""), append=TRUE)
}

# Load in slice tree REP data
slicePhys<-read.tree(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target_slicePhy-5to70Ma_by5.trees",sep=""))

# load back in tables
results<-vector("list",length(allCladeSets))
for (i in 1: length(allCladeSets)){
	res<-read.table(paste("MamPhy_5911sp_",bbone,"_cladeLevel_DRstats_trees",allCladeSetNames[i],".txt",sep=""), header=TRUE)
	results[[i]]<-res
}

# load back in slice tree CLADE data
allCladeSets<-vector("list",length(allCladeSetNames))
for (i in 1:length(allCladeSetNames)){
	allCladeSets[[i]]<-read.tree(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target_slice",allCladeSetNames[[i]],"_newick.trees",sep=""))
	}

allCladeSetNames<-c("5Ma","10Ma","15Ma","20Ma","25Ma","30Ma","35Ma","40Ma","45Ma","50Ma","55Ma","60Ma","65Ma","70Ma")
mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 
root=max(nodeHeights(mamPhy)[,2])

####
# MULTI-VARIATE PGLS...

# example: 1060 clades, richness ~ MRCA + DR_harm
i=1
rownames(results[[i]])<-slicePhys[[i]]$tip.label
cladeData<-treedata(slicePhys[[i]],na.omit(results[[i]]))

form<-(log(richness) ~ MRCA_mcc + DR_harm + DR_skew)
fit<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=as.data.frame(cladeData$data), method="ML")
sum<-summary(fit) #ss<-residuals(fit)

fit2<-gls(form, correlation=corPagel(value=1,phy=cladeData$phy), data=as.data.frame(cladeData$data), method="ML")
sum2<-summary(fit2)

plot(form,data=cladeData$data, cex.lab=1.5, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1])," clades", sep=""))
if (sum$tTable[8] < 0.05){
	abline(fit$coef[[1]], fit$coef[[2]], col="red", lwd=2, lty=2)
}
mtext(paste("PGLS p-val = ",round(sum$tTable[8],digits=4),"; X = ",round(fit$coef[[2]],digits=2),"; AIC = ",round(sum$AIC,digits=1),sep=""))

# example: using phyl.resid()...
i=1
rownames(results[[i]])<-slicePhys[[i]]$tip.label
cladeData<-treedata(slicePhys[[i]],na.omit(results[[i]]))

Y1<-log(cladeData$data[,6]) # log(richness)
x123<-cbind(cladeData$data[,7],cladeData$data[,1],cladeData$data[,3]) # MRCA_mcc + DR_harm + DR_skew
colnames(x123)<-c("MRCA_mcc","DR_harm","DR_skew")
fitP<-phyl.resid(tree=cladeData$phy, x=x123, Y=Y1, method="BM")

fitP$beta # a vector or matrix of regression coefficients.
               [,1]
         0.14554264
MRCA_mcc 0.17692702
DR_harm  5.43026651
DR_skew  0.02434391

sum$coef
(Intercept)    MRCA_mcc     DR_harm     DR_skew 
 0.14554247  0.17692702  5.43026651  0.02434391 



# TIMESLICES - start HERE
for (i in 1:length(results)){

# richness ~ MRCA
form=(log(richness) ~ MRCA_mcc)
	rownames(results[[i]])<-slicePhys[[i]]$tip.label
	cladeData<-treedata(slicePhys[[i]],na.omit(results[[i]]))
	fit<-pgls(form, correlation=corBrownian(phy=cladeData$phy), data=as.data.frame(cladeData$data), method="ML")

dd<-summary(fit)
sink(file=paste("cladeLevel",bbone,"_PGLS_bySLICE_", allCladeSetNames[i],"_RESULTS.txt",sep=""), append=TRUE)
print(form)
print(dd)
sink()
ss<-as.character(form)

plot(fit$residuals[1:length(fit$residuals)]~cladeData$data[,7])

pdf(file=paste("cladeLevel",bbone,"_PGLS_bySLICE_", allCladeSetNames[i],"_withPartialResiduals.pdf",sep=""), onefile=TRUE, width=5, heigh=5)

termplot(fit, partial.resid=TRUE, rug=TRUE, cex.main=1, font.sub=2)#, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1]),"; ",ss[2],ss[1],ss[3], sep=""), sub=paste("X1=", round(dd$coef[[2]],digits=3), "; AIC=",round(dd$aic,digits=0),sep=""))

# ^^ Ended up not being able to use TERMPLOT with GLS object-- no big deal really... the full model was the plan anyway.


# FULL MODEL - across all time slices:
##
# Now Multi-var PGLS across all time slices::

partSlopesPGLS<-data.frame(matrix(NA, nrow = length(results), ncol = 7))
sliceTimes<-seq(-5,-70,-5)
partSlopesPGLS_All<-cbind(sliceTimes,partSlopesPGLS)
colnames(partSlopesPGLS_All)<-c("time","MRCA_mcc","DR_harm","DR_skew","AIC","Pval1","Pval2","Pval3")

for (i in 1:length(results)){
	rownames(results[[i]])<-slicePhys[[i]]$tip.label
	cladeData<-treedata(slicePhys[[i]],na.omit(results[[i]]))

	form<-(log(richness) ~ MRCA_mcc + DR_harm + DR_skew)
	fit<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=as.data.frame(cladeData$data), method="ML")
	sum<-summary(fit)

	partSlopesPGLS_All[i,2]<-round(sum$coef[[2]],digits=3)
	partSlopesPGLS_All[i,3]<-round(sum$coef[[3]],digits=3)
	partSlopesPGLS_All[i,4]<-round(sum$coef[[4]],digits=3)
	partSlopesPGLS_All[i,5]<-round(sum$AIC,digits=0)
	partSlopesPGLS_All[i,6]<-round(sum$tTable[14], digits=3)
	partSlopesPGLS_All[i,7]<-round(sum$tTable[15], digits=3)
	partSlopesPGLS_All[i,8]<-round(sum$tTable[16], digits=3)

}
write.table(partSlopesPGLS_All,paste("cladeLevel",bbone,"_PGLSmulti_timeSlice_5Ma-to-70Ma_partSlopes.txt",sep=""))

# plot through time
pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_timeSlice_5Ma-to-70Ma_partSlopes_PLOTTED_3part.pdf",sep=""),onefile=TRUE, width=10,height=4)

#quartz(width=10,height=4)
#layout(matrix(c(1:3), 3, 1, byrow = TRUE))
par(mfrow = c(1,3), oma = c(5,4,5,0) + 0.1, mar = c(1,1,1,1) + 0.1)
#par(op)

form=(MRCA_mcc ~ time)
y.lab=""
x.lab=""
plot(form, data=partSlopesPGLS_All, xlab=x.lab,ylab=y.lab)
#axis(side=1,labels=FALSE)
mtext(text="(a) MRCA", adj=0, font=2)#at=c(-72,0.21)

fit2<-lm(form,data=partSlopesPGLS_All)
dd2<-summary(fit2)
if (dd2$coef[8] < 0.05){
	abline(fit2$coef[1],fit2$coef[2], lty=2, lwd=2, col="blue")
	text(x=-72,y=0.17,pos=4,font=1, labels=bquote( plain("R") ^ plain("2") * plain(" = ") * .(round(dd2$r.squared, digits=3)) ))
	#text(x=-72,y=0.215,adj=c(-0.1,2),font=1, labels=bquote( plain("X = ") * .(round(fit2$coef[2],digits=3)) ))
}

form=(DR_harm ~ time)
plot(form, data=partSlopesPGLS_All, xlab=x.lab,ylab=y.lab)
#axis(side=1,labels=FALSE)
mtext(text="(b) DR mean", adj=0, font=2)#at=c(-72,0.21)

fit2<-lm(form,data=partSlopesPGLS_All)
dd2<-summary(fit2)
if (dd2$coef[8] < 0.05){
	abline(fit2$coef[1],fit2$coef[2], lty=2, lwd=2, col="blue")
	text(x=-18,y=13.65,pos=4,font=1, labels=bquote( plain("R") ^ plain("2") * plain(" = ") * .(round(dd2$r.squared, digits=3)) ))
	#text(x=-18,y=11.3,adj=c(-0.1,2),font=1, labels=bquote( plain("X = ") * .(round(fit2$coef[2],digits=3)) ))
}

form=(DR_skew ~ time)
plot(form, col="white", data=partSlopesPGLS_All, xlab=x.lab,ylab=y.lab)
points(form, data=partSlopesPGLS_All[which(partSlopesPGLS_All$Pval3 <= 0.05),], xlab=x.lab,ylab=y.lab)
points(form, col="red",pch=2, data=partSlopesPGLS_All[which(partSlopesPGLS_All$Pval3 > 0.05),], xlab=x.lab,ylab=y.lab)
mtext(text="(c) DR skew", adj=0, font=2)#at=c(-72,0.21)

fit2<-lm(form,data=partSlopesPGLS_All)
dd2<-summary(fit2)
if (dd2$coef[8] < 0.05){
	abline(fit2$coef[1],fit2$coef[2], lty=2, lwd=2, col="blue")
	text(x=-18,y=0.65,pos=4,font=1, labels=bquote( plain("R") ^ plain("2") * plain(" = ") * .(round(dd2$r.squared, digits=3)) ))
	#text(x=-18,y=0.9,adj=c(-0.1,2),font=1, labels=bquote( plain("X = ") * .(round(fit2$coef[2],digits=3)) ))
}
title(main="PGLS: log(richness) ~ MRCA_mcc + DR_harm + DR_skew", xlab = "Time slices before present (Ma)",
      ylab = "Slope of partial residual",
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)

dev.off()

# VERTICAL plot through time
pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_timeSlice_5Ma-to-70Ma_partSlopes_PLOTTED_4partVERTICAL.pdf",sep=""),onefile=TRUE, width=4,height=13.33)

#quartz(width=10,height=4)
#layout(matrix(c(1:3), 3, 1, byrow = TRUE))
par(mfrow = c(4,1), oma = c(5,4,5,0) + 0.1, mar = c(1,1,1,1) + 0.1)
#par(op)

form=(MRCA_mcc ~ time)
y.lab=""
x.lab=""
plot(form, data=partSlopesPGLS_All, xlab=x.lab,ylab=y.lab,xaxt="n")
axis(side=1,labels=FALSE)
mtext(text="(a) MRCA", adj=0, font=2)#at=c(-72,0.21)

fit2<-lm(form,data=partSlopesPGLS_All)
dd2<-summary(fit2)
if (dd2$coef[8] < 0.05){
	abline(fit2$coef[1],fit2$coef[2], lty=2, lwd=2, col="blue")
	text(x=-72,y=0.17,pos=4,font=1, labels=bquote( plain("R") ^ plain("2") * plain(" = ") * .(round(dd2$r.squared, digits=3)) ))
	#text(x=-72,y=0.215,adj=c(-0.1,2),font=1, labels=bquote( plain("X = ") * .(round(fit2$coef[2],digits=3)) ))
}
for(i in 1:length(results)){
	abline(v=(-5*i),lty=2, lwd=1.5, col=grey(0.6, alpha=0.5))
}

form=(DR_harm ~ time)
plot(form, data=partSlopesPGLS_All, xlab=x.lab,ylab=y.lab,xaxt="n")
axis(side=1,labels=FALSE)
mtext(text="(b) DR mean", adj=0, font=2)#at=c(-72,0.21)

fit2<-lm(form,data=partSlopesPGLS_All)
dd2<-summary(fit2)
if (dd2$coef[8] < 0.05){
	abline(fit2$coef[1],fit2$coef[2], lty=2, lwd=2, col="blue")
	text(x=-18,y=13.65,pos=4,font=1, labels=bquote( plain("R") ^ plain("2") * plain(" = ") * .(round(dd2$r.squared, digits=3)) ))
	#text(x=-18,y=11.3,adj=c(-0.1,2),font=1, labels=bquote( plain("X = ") * .(round(fit2$coef[2],digits=3)) ))
}
for(i in 1:length(results)){
	abline(v=(-5*i),lty=2, lwd=1.5, col=grey(0.6, alpha=0.5))
}

form=(DR_skew ~ time)
plot(form, col="white", data=partSlopesPGLS_All, xlab=x.lab,ylab=y.lab)
points(form, data=partSlopesPGLS_All[which(partSlopesPGLS_All$Pval3 <= 0.05),], xlab=x.lab,ylab=y.lab)
points(form, col="red",pch=2, data=partSlopesPGLS_All[which(partSlopesPGLS_All$Pval3 > 0.05),], xlab=x.lab,ylab=y.lab)
mtext(text="(c) DR skew", adj=0, font=2)#at=c(-72,0.21)

fit2<-lm(form,data=partSlopesPGLS_All)
dd2<-summary(fit2)
if (dd2$coef[8] < 0.05){
	abline(fit2$coef[1],fit2$coef[2], lty=2, lwd=2, col="blue")
	text(x=-18,y=0.65,pos=4,font=1, labels=bquote( plain("R") ^ plain("2") * plain(" = ") * .(round(dd2$r.squared, digits=3)) ))
	#text(x=-18,y=0.9,adj=c(-0.1,2),font=1, labels=bquote( plain("X = ") * .(round(fit2$coef[2],digits=3)) ))
}
for(i in 1:length(results)){
	abline(v=(-5*i),lty=2, lwd=1.5, col=grey(0.6, alpha=0.5))
}

# for FBD1
#plot(mamPhy, edge.width=0.7, cex=0.1, show.tip.label=FALSE, edge.color=grey(0), no.margin=FALSE, x.lim=c(194,259)) #, tip.color = tcol )
#axisPhylo(cex.axis=1, pos=-70, mgp=c(0,0.2,0))
#for(i in 1:length(results)){
#	abline(v=root-(5*i),lty=2, lwd=1.5, col=grey(0.6))
#}

# for NDexp
plot(mamPhy, edge.width=0.3, cex=0.1, show.tip.label=FALSE, no.margin=FALSE, x.lim=c(111,176)) #, tip.color = tcol )
axisPhylo(cex.axis=1, pos=-70, mgp=c(0,0.2,0))
for(i in 1:length(results)){
	abline(v=root-(5*i),lty=2, lwd=2, col="grey")
}


title(main="PGLS: log(richness) ~ MRCA_mcc + DR_harm + DR_skew", xlab = "Time slices before present (Ma)",
      ylab = "Slope of PGLS term (partial residual)",
      outer = TRUE, line = 3,cex.main=1,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)

dev.off()


##
# Plot MAMPHY MCC without tips, and vertical, but attractive... 
quartz(width=4, height=10)
plot(ladderize(mamPhy), show.tip.label=FALSE)

# Plot with TIMESLICES shown...
jpeg(file=paste("plot_MamPhy_MCC_",bbone,"_showingSlices.jpg",sep=""),width=10,height=10,units="in", res=300, quality=100)
	par(oma=c(0,0,0,0),mar=c(0,0,0,0))
	plot(mamPhy, edge.width=0.7, cex=0.1, show.tip.label=FALSE) #, tip.color = tcol )
	axisPhylo(cex.axis=1, pos=-13, mgp=c(0,0.2,0))
	for(i in 1:length(results)){
		abline(v=root-(5*i),lty=2, lwd=2, col="grey")
	}
dev.off()

## 
# ORDERS and FAMS- now plot the 3-part PGLS for these groups
ordRes<-read.table("MamPhy_5911sp_FBD_cladeLevel_ORDS_skewKurt_withMRCAs.txt")
famRes<-read.table("MamPhy_5911sp_FBD_cladeLevel_FAMS_skewKurt_withMRCAs.txt")

ordNames<-rownames(ordRes[with(ordRes,order(ordRes$richness, decreasing=TRUE)),]) # is the ORDER names by category (richness)
famNames<-rownames(famRes[with(famRes,order(famRes$richness, decreasing=TRUE)),]) # is the ORDER names by category (richness)

ordPhy<-read.tree("MamPhy_BDvr_pcsFIXED_FBD_MCC_target_27ORDERS.tre")
famPhy<-read.tree("MamPhy_BDvr_pcsFIXED_FBD_MCC_target_162FAMILIES.tre")




# colored ORDS

# Can do this in PHYLOCH
# example::
data(bird.orders)
clade1 <- c("Struthioniformes", "Tinamiformes", "Craciformes", "Galliformes", "Anseriformes")
clade2 <- c("Galbuliformes", "Bucerotiformes", "Upupiformes", "Trogoniformes", "Coraciiformes")
clade3 <- c("Coliiformes", "Cuculiformes", "Psittaciformes", "Apodiformes", "Trochiliformes","Musophagiformes", "Strigiformes", "Columbiformes", "Gruiformes", "Ciconiiformes", "Passeriformes")
clade.list <- list(clade1, clade2, clade3)
color <- c("blue", "red", "green")
tcol <- tip.color(phy=bird.orders, tips=clade.list, col=color)
ecol <- edge.color(phy=bird.orders, groups=clade.list, col = color)
par(mfrow=c(1,2)) # a second plot will follow
plot(bird.orders, tip.color = tcol, edge.color = ecol)

# with ORDS
ordNames2<-names(which(sort(table(cladesDR$ord),decreasing=TRUE) > 2))
cols<-palette(rainbow(length(ordNames2)))

ordTipNames<-vector("list",length(ordNames2))
for (i in 1:length(ordTipNames)){
	ordTipNames[[i]]<-as.character(cladesDR[which(cladesDR$ord==ordNames2[i]),"tiplabel"])
}

mamPhy<-drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target.tre",sep="")), "_Anolis_carolinensis") ## Use the MCC target node heights one... 

#tcols<-tip.color(phy=mamPhy, tips=ordTipNames, col=cols)
ecols <- edge.color(phy=mamPhy, groups=ordTipNames, col = cols, bgcol="black", what="crown")
table(ecols)

mamPhy2<-ladderize(mamPhy)
root<-max(nodeHeights(mamPhy)[,2])

pdf(file=paste("plot_MamPhy_MCC_",bbone,"_colOrds_0.7.pdf",sep=""),onefile=TRUE,width=4,height=10)
	par(oma=c(0,0,0,0),mar=c(0,0,0,0))
	plot(mamPhy2, edge.color = ecols[match(mamPhy2$edge[,2],mamPhy$edge[,2])], edge.width=0.7, cex=0.1, show.tip.label=FALSE) #, tip.color = tcol )
	axisPhylo(cex.axis=0.5, pos=-13, mgp=c(0,0.2,0))
	abline(v=root-65.5,lty=2)
dev.off()


pdf(file=paste("plot_MamPhy_MCC_",bbone,"_colOrds_fullSize.pdf",sep=""),onefile=TRUE,width=22,height=150)
	plot(mamPhy2, edge.color = ecols[match(mamPhy2$edge[,2],mamPhy$edge[,2])], cex=0.2, show.tip.label=TRUE) #, tip.color = tcol )
dev.off()




for (i in 1:length(ordNames2)){
	ordSp<-cladesDR[which(cladesDR$ord==ordNames2[i]),"tiplabel"]
	ordNode <- getMRCA(mamPhy, ordSp)
	mamPhy<-paintSubTree(mamPhy,node=ordNode,state=i,stem=FALSE)
}

pdf(file=paste("plot_MamPhy_MCC_",bbone,"_colOrds.pdf",sep=""),onefile=TRUE,width=4,height=10)
plotSimmap(mamPhy, cols, ftype="off")
dev.off()





colTips<-c(rep("black",length(mamPhy$edge[,1])))

	colTips[1:5911][mamPhy$tip.label %in% ordSp]<-cols[i]

   # Assign rates to tips
    tip.idx <- sort(mamPhy$tip.label, index.return=TRUE)$ix

    colTips[match(tip.idx, mamPhy$edge[,2])] <- DR_median[sort(names(DR_median))]





tree<-paintSubTree(tree,node=26,state="3")
> # now let's plot using plotSimmap to ensure
> # that the correct branches were painted
> cols<-c("black","blue","red"); names(cols)<-1:3
> plotSimmap(tree,cols,pts=F,lwd=3,node.numbers=T)



# ORDS - start HERE
ordData<-na.omit(ordRes)	
cladeData<-treedata(ordPhy,ordData)

form<-(log(richness) ~ MRCA_mcc + DR_harm + DR_skew)
fit<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=as.data.frame(cladeData$data), method="ML")
sum<-summary(fit)



plot(form,data=cladeData$data, cex.lab=1.5, main=paste("slice ", allCladeSetNames[i], ", ",length(results[[i]][,1])," clades", sep=""))
if (sum$tTable[14:16] < 0.05){
	abline(fit$coef[[1]], fit$coef[[2]], col="red", lwd=2, lty=2)
}
mtext(paste("PGLS p-val = ",round(sum$tTable[8],digits=4),"; X = ",round(fit$coef[[2]],digits=2),"; AIC = ",round(sum$AIC,digits=1),sep=""))


# FAMS - start HERE
famData<-na.omit(famRes)	
cladeData<-treedata(famPhy,famData)

form<-(log(richness) ~ MRCA_mcc + DR_harm + DR_skew)
fit<-gls(form, correlation=corBrownian(phy=cladeData$phy), data=as.data.frame(cladeData$data), method="ML")
sum<-summary(fit)



##



## Now plot each term separately with PGLS against log(richness)

pdf(file=paste("cladeLevel",bbone,"_PGLSmulti_timeSlice_5Ma-to-70Ma_partSlopes_PLOTTED_3part.pdf",sep=""),onefile=TRUE, width=10,height=4)

#quartz(width=10,height=4)
#layout(matrix(c(1:3), 3, 1, byrow = TRUE))






title(main="PGLS: log(richness) ~ MRCA_mcc + DR_harm + DR_skew", xlab = "Time slices before present (Ma)",
      ylab = "Slope of partial residual",
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)

dev.off()



####
# UNIVARIATE PGLS...

# richness ~ MRCA

pdf(file=paste("cladeLevel",bbone,"_slices_by5Ma_PGLS_rich-by-MRCA.pdf",sep=""), onefile=TRUE, width=5, heigh=5)

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

pdf(file=paste("cladeLevel",bbone,"slices_by5Ma_PGLS_rich-by-DR_harm.pdf",sep=""), onefile=TRUE, width=5, heigh=5)

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

pdf(file=paste("cladeLevel",bbone,"_slices_by5Ma_PGLS_rich-by-DR_skew.pdf",sep=""), onefile=TRUE, width=5, heigh=5)

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

pdf(file=paste("cladeLevel",bbone,"_slices_by5Ma_PGLS_rich-by-DR_kurt.pdf",sep=""), onefile=TRUE, width=5, heigh=5)

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

pdf(file=paste("cladeLevel",bbone,"_slices_by5Ma_PGLS_rich-by-DR_cv.pdf",sep=""), onefile=TRUE, width=5, heigh=5)

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

# original, for newick tree: ymle = function(tree) (.subset2(tree,2)-1L)/sum(.subset2(tree,4)) # this take the # of number of nodes in a tree / sum of branch lengths.













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





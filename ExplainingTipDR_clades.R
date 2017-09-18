# =================================
# Explaining tip DR mean per clade
###
# Time slices
###
# =================================
######
# Load data
library(ape); library(phytools)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"

mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")

tipDataAll_wComments_DETAILS<-read.table("MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass.txt", header=TRUE)
tipDataAll_1<-tipDataAll_wComments_DETAILS[,c("GENUS","FAMILY","ORDER","perVert","perInvert","perPlant","trophicCat","CarnOrNot","HerbOrNot","OmniOrNot","lifemodeCat6","AquaOrNot","ArboOrNot","FlysOrNot","SubtOrNot","TerrOrNot","activityCat3","NoctOrNot","CathOrNot","DiurOrNot","tipDR_mean10k","geoRange_numGrid","geoRange_km2","BM_final_kg","BM_final_g")]
tipDataAll_2<-tipDataAll_wComments_DETAILS[,c(52:81,83:101)]
tipDataAll<-cbind.data.frame(tipDataAll_1,tipDataAll_2)
rownames(tipDataAll)<-tipDataAll_wComments_DETAILS$tiplabel

# Re-calculating the home range to mass relationship in Jetz et al. 2005

datHR_all<-as.data.frame(na.omit(tipDataAll[,c("trophicCat","X22.2_HomeRange_Indiv_km2","BM_final_kg")])) # exclude "X22.1_HomeRange_km2" not corrected for group size
colnames(datHR_all)<-c("trophicCat","HomeRange_Indiv_km2","BodyMass_kg")
	# gives 606 entries

pdf(file="reCalc_of_log-log_HR-vs-BM_mammalsTrophic_4way.pdf", onefile=TRUE, height=8, width=8)

cats<-c("All mammals:", "Herbivores:", "Omnivores:", "Carnivores:")
dats<-list(datHR_all,datHR_all[which(datHR_all$trophicCat=="Herb"),],datHR_all[which(datHR_all$trophicCat=="Omni"),],datHR_all[which(datHR_all$trophicCat=="Carn"),])

HR_equations<-data.frame(matrix(NA, nrow = length(cats), ncol = 5))
colnames(HR_equations)<-c("log10_Int_km2","log10_Slope_kg","nonLog_Int_ha","nonLog_Slope_Kg","R2")
rownames(HR_equations)<-cats

layout(matrix(c(1:4), 2, 2, byrow = TRUE))#, widths=4, heights=4)
par(oma = c(5,4,5,3) + 0.1, mar = rep(1.7,4) + 0.1)

for(i in 1:length(cats)){
	cat<-cats[i]
	dat<-dats[[i]]

	form<-log10(HomeRange_Indiv_km2) ~ log10(BodyMass_kg)
	plot(form, data=dat)
	fit<-lm(form,data=dat)
	sum<-summary(fit)

	n=length(dat[,1])
	a=round(sum$coef[1],2)
	b=round(sum$coef[2],2)
	R2=round(sum$r.squared,2)
	X=seq(min(log10(dat$BodyMass_kg)),max(log10(dat$BodyMass_kg)),length.out=20)
	lines(x=X,y=b*X+a,col="black", lwd=4,lty=1) 

	mtext(side=3, text=bquote(bold(.(cat) ~ n ~ "=" ~ .(n) ~~~~ R^2 ~ "=" ~ .(R2))),adj=0, cex=0.7)

	Y=log10(dat$HomeRange_Indiv_km2)
	text(min(X),max(Y),label=bquote(log10 ~ km^2 ~ bold(y ~ "=" ~ .(a)+.(b)*x)~log10 ~ kg),adj=0, cex=0.8)

	text(min(X),max(Y)-(max(Y)-min(Y))/12,label=bquote(ha ~ bold(y ~ "=" ~ .(round(100*10^a,2))+.(round(10^b,2))*x)~kg),adj=0, cex=0.8)

	HR_equations[i,1]<-a
	HR_equations[i,2]<-b
	HR_equations[i,3]<-round(100*10^a,2)
	HR_equations[i,4]<-round(10^b,2)
	HR_equations[i,5]<-R2

}
title(main="", xlab = "log10(Body Mass in kg)",
      ylab = "log10(Home Range per indiv in km2)",
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)


dev.off()

write.table(HR_equations,file="reCalc_of_log-log_HR-vs-BM_mammalsTrophic_4way_equations.txt")

## Now do those HR calcs
datHR_all<-as.data.frame(na.omit(tipDataAll[,c("trophicCat","BM_final_kg")])) # exclude "X22.1_HomeRange_km2" not corrected for group size
colnames(datHR_all)<-c("trophicCat","BodyMass_kg")
	# gives 4988 entries

cats<-c("All mammals:", "Herbivores:", "Omnivores:", "Carnivores:")
dats<-list(datHR_all,datHR_all[which(datHR_all$trophicCat=="Herb"),],datHR_all[which(datHR_all$trophicCat=="Omni"),],datHR_all[which(datHR_all$trophicCat=="Carn"),])

HR_ext_vals<-data.frame(matrix(NA, nrow = length(tipDataAll[,1]), ncol = 1))
rownames(HR_ext_vals)<-tipDataAll_wComments_DETAILS$tiplabel

for(i in 2:length(cats)){
	cat<-cats[i]
	dat<-dats[[i]]

	a=HR_equations[i,1]
	b=HR_equations[i,2]

	for(j in 1:length(dat[,1])){
		HR_j <- a+(log(dat[j,"BodyMass_kg"])*b)
		names(HR_j)<-rownames(dat[j,])
		HR_ext_vals[which(rownames(HR_ext_vals)==names(HR_j)),]<-HR_j
	}
}
colnames(HR_ext_vals)<-"log10homeRange_km2_ext"

homeRange_km2_ext<-10^HR_ext_vals
colnames(homeRange_km2_ext)<-"homeRange_km2_ext"

HR_only<-cbind(HR_ext_vals,homeRange_km2_ext)
write.table(HR_only,file="MamPhy_5911sp_homeRangeExt_only.txt")

tipDataAllHR<-cbind(tipDataAll,HR_ext_vals,homeRange_km2_ext)
write.table(tipDataAllHR,file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_homeRangeExt.txt")

####
# Now do the same extrapolation for the DISPERSAL DISTANCE
library(ape); library(phytools)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"

tipDataAllHR<-read.table(file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_homeRangeExt.txt")

# BIN the home range data also...
lnHomeRange_km2_ext<-as.vector(log(tipDataAllHR$homeRange_km2_ext))
names(lnHomeRange_km2_ext)<-rownames(tipDataAllHR)
#binsEven
	bins<-quantile(lnHomeRange_km2_ext,probs=seq(0,1,by=0.1),na.rm=TRUE)
	lnHomeRange_rank1to10_ext<-data.frame(matrix(NA, nrow = length(tipDataAllHR[,1]), ncol = 1))
	rownames(lnHomeRange_rank1to10_ext)<-rownames(tipDataAllHR)
	colnames(lnHomeRange_rank1to10_ext)<-"lnHomeRange_rank1to10"

	toBin<-as.matrix(na.omit(lnHomeRange_km2_ext))

	for(i in 1:length(toBin)){
		for(j in 1:(length(bins)-1)){
			if(toBin[[i]] >= bins[[j]] && toBin[[i]] <= bins[[j+1]]){
			lnHomeRange_rank1to10_ext[which(rownames(lnHomeRange_rank1to10_ext)==names(toBin[i,])),]<-j} else {next}
		}
	}

#binsEqual
	range<-range(na.omit(lnHomeRange_km2_ext))
	binsEqual<-seq(from=range[1],to=range[2],length.out=11)
	lnHomeRange_rank1to10_binsEqual<-data.frame(matrix(NA, nrow = length(tipDataAllHR[,1]), ncol = 1))
	rownames(lnHomeRange_rank1to10_binsEqual)<-rownames(tipDataAllHR)
	colnames(lnHomeRange_rank1to10_binsEqual)<-"lnHomeRange_rank1to10_binsEqual"

	toBin<-as.matrix(na.omit(lnHomeRange_km2_ext))

	for(i in 1:length(toBin)){
		for(j in 1:(length(binsEqual)-1)){
			if(toBin[[i]] >= binsEqual[j] && toBin[[i]] <= binsEqual[j+1]){
			lnHomeRange_rank1to10_binsEqual[which(rownames(lnHomeRange_rank1to10_binsEqual)==names(toBin[i,])),]<-j} else {next}
		}
	}

#hist(lnDispDistAll_rank1to10[,1],breaks=30)
pdf(file="lnHomeRange_km2_ext_histoAsBinned.pdf")
hist(toBin,breaks=30, col="grey")
for(j in 1:(length(bins)-1)){
	abline(v=bins[[j+1]], col="blue",lwd=2)
}
dev.off()

pdf(file="lnHomeRange_rank1to10_binsEqual_histoAsBinned.pdf")
hist(toBin,breaks=30, col="grey")
for(j in 1:(length(binsEqual)-1)){
	abline(v=binsEqual[j+1], col="blue",lwd=2)
}
dev.off()



##########
# calc dispersal distance...
datForDispersal<-as.data.frame(tipDataAllHR[,c("BM_final_g","homeRange_km2_ext","geoRange_km2")])
	# gives 4961 values

# Dispersal distance = -1.153 + 0.315*bodyMass(g) + 0.220*homeRange(km2) + 0.252*geoRange(km2)

lnBM<-log(datForDispersal[,"BM_final_g"])
lnHR<-log(datForDispersal[,"homeRange_km2_ext"])
lnGR<-log(datForDispersal[,"geoRange_km2"])
lnDispDistAll_km_ext = -1.153 + 0.315*lnBM + 0.220*lnHR + 0.252*lnGR
names(lnDispDistAll_km_ext)<-rownames(datForDispersal)

DispDistAll_km_ext <- exp(lnDispDistAll_km_ext)

#binsEven
	bins<-quantile(lnDispDistAll_km_ext,probs=seq(0,1,by=0.1),na.rm=TRUE)
	lnDispDistAll_rank1to10_ext<-data.frame(matrix(NA, nrow = length(lnDispDistAll_km_ext), ncol = 1))
	rownames(lnDispDistAll_rank1to10_ext)<-names(lnDispDistAll_km_ext)
	colnames(lnDispDistAll_rank1to10_ext)<-"lnDispDistAll_rank1to10"

	toBin<-na.omit(lnDispDistAll_km_ext)

	for(i in 1:length(toBin)){
		for(j in 1:(length(bins)-1)){
			if(toBin[[i]] >= bins[[j]] && toBin[[i]] <= bins[[j+1]]){
			lnDispDistAll_rank1to10_ext[which(rownames(lnDispDistAll_rank1to10_ext)==names(toBin[i])),]<-j} else {next}
		}
	}

lnDispDistAll_km_ext<-as.vector(log(tipDataAll$DispDistAll_km_ext))
names(lnDispDistAll_km_ext)<-tipDataAll$tiplabel
#binsEqual
	range<-range(na.omit(lnDispDistAll_km_ext))
	binsEqual<-seq(from=range[1],to=range[2],length.out=11)
	lnDispDistAll_rank1to10_binsEqual<-data.frame(matrix(NA, nrow = length(tipDataAll[,1]), ncol = 1))
	rownames(lnDispDistAll_rank1to10_binsEqual)<-tipDataAll$tiplabel
	colnames(lnDispDistAll_rank1to10_binsEqual)<-"lnDispDistAll_rank1to10_binsEqual"

	toBin<-as.matrix(na.omit(lnDispDistAll_km_ext))

	for(i in 1:length(toBin)){
		for(j in 1:(length(binsEqual)-1)){
			if(toBin[[i]] >= binsEqual[j] && toBin[[i]] <= binsEqual[j+1]){
			lnDispDistAll_rank1to10_binsEqual[which(rownames(lnDispDistAll_rank1to10_binsEqual)==names(toBin[i,])),]<-j} else {next}
		}
	}

#hist(lnDispDistAll_rank1to10[,1],breaks=30)
pdf(file="lnDispDistAll_km_ext_histoAsBinned.pdf")
hist(toBin,breaks=30, col="grey")
for(j in 1:(length(bins)-1)){
	abline(v=bins[[j+1]], col="blue",lwd=2)
}
dev.off()

pdf(file="lnDispDistAll_rank1to10_binsEqual_histoAsBinned.pdf")
hist(toBin,breaks=30, col="grey")
for(j in 1:(length(binsEqual)-1)){
	abline(v=binsEqual[j+1], col="blue",lwd=2)
}
dev.off()



DispDist_only<-cbind(lnBM,lnHR,lnGR,lnDispDistAll_km_ext,DispDistAll_km_ext,lnDispDistAll_rank1to10_ext,lnDispDistAll_rank1to10_binsEqual,lnHomeRange_rank1to10_ext,lnHomeRange_rank1to10_binsEqual)
write.table(DispDist_only,file="MamPhy_5911sp_dispersalDistance_Ext_only.txt")

tipDataAll_ext<-cbind(tipDataAll,HR_ext_vals,homeRange_km2_ext,lnDispDistAll_km_ext,DispDistAll_km_ext,lnDispDistAll_rank1to10_ext)
write.table(tipDataAll_ext,file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp.txt")

#####
# validate the home range calcs...
library(ape); library(phytools)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"

tipDataAll_ext<-read.table(file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp.txt")

dat<-na.omit(tipDataAll_ext[,c("homeRange_km2_ext","X22.2_HomeRange_Indiv_km2")])
	# gives 606 species, as expected.
form<-log(homeRange_km2_ext) ~ log(X22.2_HomeRange_Indiv_km2)
plot(form, data=dat)
fit<-lm(form,data=dat)
sum<-summary(fit)
	# Adjusted R-squared:  0.734 >> that makes sense, follows from the 3 prediction equations used.  Cool.

###
# Look at PAIRWISE with tip DR relationships with these new variables
##
library(nlme); library(geiger)
tipDataAll<-read.table(file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)
mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")

# as continuous...
dat1<-tipDataAll[,c("tipDR_mean10k","GenerationLength_d")]
rownames(dat1)<-tipDataAll$tiplabel
cladeData<-treedata(mamMCC,na.omit(dat1)) 
dat<-as.data.frame(cladeData$data)

		#dat<-na.omit(tipDataAll[,c("tipDR_mean10k","GenerationLength_d")])
form<-(tipDR_mean10k ~ log(GenerationLength_d))
plot(form, data=tipDataAll)

fit<-gls(form, data=dat, correlation=corPagel(value=1,phy=cladeData$phy))
sum<-summary(fit)


dat<-na.omit(tipDataAll[,c("tipDR_mean10k","lnHomeRange_rank1to10")])
form<-(tipDR_mean10k ~ (lnHomeRange_rank1to10))
plot(form, data=tipDataAll)

fit<-lm(form, data=dat)
sum<-summary(fit)


dat<-na.omit(tipDataAll[,c("tipDR_mean10k","trophicCat")])
form<-(tipDR_mean10k ~ trophicCat)
plot(form, data=tipDataAll)

fit<-gls(form, data=dat)
sum<-summary(fit)



	a=round(sum$coef[1],2)
	b=round(sum$coef[2],2)
	R2=round(sum$r.squared,2)
	X=seq(min(log10(dat$BodyMass_kg)),max(log10(dat$BodyMass_kg)),length.out=20)
	lines(x=X,y=b*X+a,col="black", lwd=4,lty=1) 

	mtext(side=3, text=bquote(bold(.(cat) ~ n ~ "=" ~ .(n) ~~~~ R^2 ~ "=" ~ .(R2))),adj=0, cex=0.7)

	Y=log10(dat$HomeRange_Indiv_km2)
	text(min(X),max(Y),label=bquote(log10 ~ km^2 ~ bold(y ~ "=" ~ .(a)+.(b)*x)~log10 ~ kg),adj=0, cex=0.8)



######
# =========
# PERFORM those TIMESLICE analyses, with summary stats for the tipDR...
# ======

setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")
#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut")
source("DR_functions.R")
library(ape); library(phytools); library(picante); library(geiger); library(moments); library(nlme)

library(foreach);library(doSNOW)
cl = makeCluster(40, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees=100

foreach(i=1:ntrees, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

# which backbone?
bbone<- "NDexp" #"FBD" # 

# read in 1 of 100 full trees
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")
#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_nexus-and-newickTrees/")
mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_nexus.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 
	#write.tree(mamPhy,file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""))
tree1=scan(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""), what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

# change to NEW directory
###
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/_ExplainingTipDR_timeSliceRuns")
#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_ExplainingTipDR_timeSliceRuns")

# load in stats about TIP SAMPLING
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
head(cladesDR)

# load in stats about TIP TRAITS
###
#tipDataAll_wComments_DETAILS<-read.table("MamPhy_5911sp_tiplabel_DR-range-mass-herb-lifemode-DETAILS.txt", header=TRUE)
#tipDataAll<-tipDataAll_wComments_DETAILS[,c(5:20,24:27)]
#rownames(tipDataAll)<-tipDataAll_wComments_DETAILS$tiplabel

tipDataAll<-read.table(file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)


# ES and DR on per-tip basis
#==================
# gives pairwise clade matrix from CAIC function
clade_matrix = readCAIC(tree1)

ES = ES_v2(clade_matrix)
DR = 1/ES
res = cbind.data.frame(DR,ES)
res1 = res[order(rownames(res)),]

write.table(res1, file=paste(bbone,"_sample100_",i,"_DRtips.txt",sep=""))

# Make time slices to create the clades
#=======================================
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery

root=max(node.age(mamPhy)$ages)

allCladeSets<-vector("list",length=numSlices)
allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSets[[j]]<-treeSlice(mamPhy, slice=root-(sliceEvery*j), trivial=FALSE)
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

# record clade sizes per slice
lengths<-vector()
for (j in 1:length(allCladeSets)){
	lengths[j]<-length(allCladeSets[[j]])
}
names(lengths)<-allCladeSetNames 
lengths
write.table(lengths, file=paste(bbone,"_sample100_",i,"_timeslice_Lengths.txt",sep=""))


#===================================
# Calculate per-SLICE, per-clade summary values 
#===================================

# get node times for tree
btimes<-branching.times(mamPhy)

# yule function
ymle = function(tree){ (.subset2(tree,3)-1L)/sum(.subset2(tree,2)) } # this take the # of number of nodes in a tree (minus 1) / sum of branch lengths.

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
	
	logBodyMass_kg<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))	
	logGeoRange_km2<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	logHomeRange_km2<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	rankHomeRange_1to10<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	rankHomeRange_1to10_binsEq<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	logDispDist_km<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	rankDispDist_1to10<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	rankDispDist_1to10_binsEq<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	logGenLength_d<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))

	trophic123<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	CarnOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	HerbOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	OmniOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	perVert<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	perInvert<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	perPlant<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	AquaOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	ArboOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	FlysOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	SubtOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	TerrOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	NoctOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	CathOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))
	DiurOrNot<-data.frame(matrix(NA, nrow = length(cladeSet), ncol = 1))

	for (k in 1:length(cladeSet)){
	cladeSp<-cladeSet[[k]]$tip.label
		x<-res1[match(cladeSp,rownames(res1)),"DR"]
		DR_harm[k,] <- 1/(mean(1/x))
		DR_cv[k,] <- (sd(x)/mean(x))*100
		DR_skew[k,] <- skewness(x)
		DR_kurt[k,] <- kurtosis(x)
		percentSamp[k,] <- length(which(cladesDR[match(cladeSp,cladesDR$tiplabel),"samp"]=="sampled"))/length(cladeSp)
		richness[k,] <- length(cladeSp)
		node <- getMRCA(mamPhy, cladeSp)
		MRCA[k,] <- btimes[node-5911] #taking the height of SAMPLED tree
	# TRAIT DATA
	# Continuous
		logBodyMass_kg[k,]<-mean(na.omit(log(tipDataAll[match(cladeSp,cladesDR$tiplabel),"BM_final_kg"])))
		logGeoRange_km2[k,]<-mean(na.omit(log(tipDataAll[match(cladeSp,cladesDR$tiplabel),"geoRange_km2"])))
		logHomeRange_km2[k,]<-mean(na.omit(log(tipDataAll[match(cladeSp,cladesDR$tiplabel),"homeRange_km2_ext"])))
		rankHomeRange_1to10[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"lnHomeRange_rank1to10"]))
		rankHomeRange_1to10_binsEq[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"lnHomeRange_rank1to10_binsEqual"]))
		logDispDist_km[k,]<-mean(na.omit(log(tipDataAll[match(cladeSp,cladesDR$tiplabel),"DispDistAll_km_ext"])))
		rankDispDist_1to10[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"lnDispDistAll_rank1to10"]))
		rankDispDist_1to10_binsEq[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"lnDispDistAll_rank1to10_binsEqual"]))
		logGenLength_d[k,]<-mean(na.omit(log(tipDataAll[match(cladeSp,cladesDR$tiplabel),"GenerationLength_d"])))
	# Categorical
		trophic123[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"trophic123"]))
		CarnOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"CarnOrNot"]))
		HerbOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"HerbOrNot"]))
		OmniOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"OmniOrNot"]))
		perVert[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"perVert"]))
		perInvert[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"perInvert"]))
		perPlant[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"perPlant"]))
		AquaOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"AquaOrNot"]))
		ArboOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"ArboOrNot"]))
		FlysOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"FlysOrNot"]))
		SubtOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"SubtOrNot"]))
		TerrOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"TerrOrNot"]))
		NoctOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"NoctOrNot"]))
		CathOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"CathOrNot"]))
		DiurOrNot[k,]<-mean(na.omit(tipDataAll[match(cladeSp,cladesDR$tiplabel),"DiurOrNot"]))

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
		} else NULL
	}

	res2<-cbind(DR_harm, DR_cv, DR_skew, DR_kurt, percentSamp, richness, MRCA, PB_Div, BD_Lam, BD_Mu, BD_Div, BD_Turn, BD.ms0, BD.ms0p5, BD.ms0p9, logBodyMass_kg,logGeoRange_km2,logHomeRange_km2,rankHomeRange_1to10,rankHomeRange_1to10_binsEq,logDispDist_km,rankDispDist_1to10,rankDispDist_1to10_binsEq,logGenLength_d,trophic123,CarnOrNot,HerbOrNot,OmniOrNot,perVert,perInvert,perPlant,AquaOrNot,ArboOrNot,FlysOrNot,SubtOrNot,TerrOrNot,NoctOrNot,CathOrNot,DiurOrNot,i, j*5)

	colnames(res2)<-c("DR_harm","DR_cv", "DR_skew", "DR_kurt", "percentSamp", "richness", "MRCA", "PB_Div", "BD_Lam", "BD_Mu", "BD_Div", "BD_Turn", "BD.ms0", "BD.ms0p5", "BD.ms0p9", "logBodyMass_kg", "logGeoRange_km2", "logHomeRange_km2", "rankHomeRange_1to10", "rankHomeRange_1to10_binsEq", "logDispDist_km", "rankDispDist_1to10", "rankDispDist_1to10_binsEq", "logGenLength_d", "trophic123", "CarnOrNot", "HerbOrNot", "OmniOrNot", "perVert","perInvert","perPlant","AquaOrNot", "ArboOrNot", "FlysOrNot", "SubtOrNot", "TerrOrNot", "NoctOrNot", "CathOrNot", "DiurOrNot", "tree", "slice")

	rownames(res2)<-paste(i,"_",j,"_",c(1:length(res2[,1])),sep="")
	
	write.table(res2,paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[j]],"cladeSTATS_withTraits.txt",sep=""))
} # end 1-tree loop across all 14 slices (every 5 Ma)


# create slice phys
slicePhys<-vector("list",length(allCladeSets))
for (j in 1:length(allCladeSets)){
	cladeSet<-allCladeSets[[j]]

	cladeReps<-vector()
	for (k in 1:length(cladeSet)){
		cladeSp<-cladeSet[[k]]$tip.label
		cladeReps[k]<-cladeSp[1]
		}
	toDrop<-setdiff(mamPhy$tip.label,cladeReps)
	slicePhys[[j]]<-drop.tip(mamPhy,toDrop)
	slicePhys[[j]]$tip.label<-paste(i,"_",j,"_",c(1:length(slicePhys[[j]]$tip.label)),sep="")
}

# write slice phys
for(j in 1:length(slicePhys)){
	write.tree(slicePhys[[j]], file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_slicePhy-5to70Ma_by5.trees",sep=""), append=TRUE)
}
} # end 100-tree loop


#####
# How many slices did we end up with?

lengths<-vector("list", length=100)
for (i in 1:length(lengths)){
	lengths[[i]]<-read.table(file=paste(bbone,"_sample100_",i,"_timeslice_Lengths.txt",sep=""))
}
lengthsALL1<-do.call(cbind,lengths)
colnames(lengthsALL1)<-paste("tree_",1:100,sep="")
lengthsALL<-as.matrix(lengthsALL1)

means<-c()
for(i in 1:length(lengthsALL[,1])){
	means[i]<-mean(lengthsALL[i,])
}
names(means)<-allCladeSetNames 
#   10Ma   15Ma   20Ma   25Ma   30Ma   35Ma   40Ma   45Ma   50Ma   55Ma   60Ma 
# 655.98 379.38 238.00 158.98 115.71  89.60  71.33  57.84  44.30  35.92  28.88 


# ==============
# do the PGLS
#######
# intialize
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/_ExplainingTipDR_timeSliceRuns")
#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_ExplainingTipDR_timeSliceRuns")
library(ape); library(phytools); library(picante); library(geiger); library(moments); library(nlme)

# which backbone?
bbone<- "NDexp" #"FBD" # 

cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
head(cladesDR)

tipDataAll<-read.table(file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)

sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery
allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

# start cluster
library(foreach);library(doSNOW)
cl2 = makeCluster(5, type = 'SOCK', outfile="")
registerDoSNOW(cl2)


#i=1; q=1
#res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[q]],"cladeSTATS_withTraits.txt",sep=""), header=TRUE)
#predictors<-colnames(res[,c(5,15:38)])#c("MRCA", "logBodyMass", "trophicIndex1", "trophicIndex2", "perVert", "perInvert", "perPlant", "numLifemode", "numActivity", "ArboOrNot","SubtOrNot","DiurOrNot")#"majorActivity") #"majorLifemode"

setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/_ExplainingTipDR_timeSliceRuns")
i=1
# Z-score STANDARDIZED from RAW
results<-vector("list",length(allCladeSetNames))
for (q in 1: length(allCladeSetNames)){
	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[q]],"cladeSTATS_withTraits.txt",sep=""), header=TRUE)
	vars<-cbind.data.frame(res[,c(5,16:39)]) # 
	varsScale<-scale(vars, center=TRUE,scale=TRUE)
	results[[q]]<-cbind(res$DR_harm,log(res$DR_harm),res$BD_Div,log(res$BD_Div),res$BD_Lam,res$BD_Mu,(res$BD_Mu/res$BD_Lam),log(res$richness),varsScale)
	colnames(results[[q]])<-c("tipDR_mean","logTipDR","BD_Div","logBD_Div","BD_Lam","BD_Mu","BD_Turn","logRichness",colnames(res[,c(5,16:39)]))
}

sliceTimes<-seq(-5,-70,-5)
predictors<-colnames(res[,c(5,16:39)])#c("MRCA", "logBodyMass", "trophicIndex1", "trophicIndex2", "perVert", "perInvert", "perPlant", "numLifemode", "numActivity", "ArboOrNot","SubtOrNot","DiurOrNot")#"majorActivity") #"majorLifemode"
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



#ntrees=100
#foreach(i=1:ntrees, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

ranTrees<-c(2, 6, 10, 11, 12, 21, 23, 24, 25, 27, 29, 31, 33, 40, 41, 48, 53, 54, 59, 60, 63, 65, 67, 68, 72, 73, 74, 77, 82, 83, 84, 85, 86, 91, 94, 96, 98)
toRun<-setdiff(1:100,ranTrees)

foreach(i=toRun, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

respVar<-c("tipDR_mean","logTipDR","BD_Div","logBD_Div","BD_Lam","BD_Mu","BD_Turn","logRichness") #c("tipDR_mean","logTipDR","BD_Div","logBD_Div","BD_Lam","BD_Mu","BD_Turn","logRichness") #c("tipDR_mean","BD_Div","BD_Lam","BD_Mu","BD_Turn") #,"richness") #
#resp<-c("TipDR", "BDdiv", "BDlam", "BDmu","BDturn") #"TraitsAsResponse"
#respVar<-c("logTipDR","logRichness","BD_Turn","BD_Lam") #c("tipDR_mean","logTipDR","BD_Div","logBD_Div","BD_Lam","BD_Mu","BD_Turn","logRichness") #c("tipDR_mean","BD_Div","BD_Lam","BD_Mu","BD_Turn") #,"richness") #
#resp<-c("logTipDR", "logRichness", "BDturn", "BDlam") #c("TipDR", "BDdiv", "BDlam", "BDmu","BDturn") #"TraitsAsResponse"

predictors<-read.table(file="predictors_names.txt")
nPreds<-length(predictors$x)
uniPGLS<-data.frame(matrix(NA, nrow = length(respVar), ncol = nPreds*5),row.names=respVar)
colnames(uniPGLS)<-c(colnames(uniPGLS_allInts),colnames(uniPGLS_allSlopes),colnames(uniPGLS_allPs),colnames(uniPGLS_allSEs),colnames(uniPGLS_allLams))

for(z in 1:length(respVar)){

# read back in slice phys
slicePhys<-read.tree(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_tree_",i,"_slicePhy-5to70Ma_by5.trees",sep=""))

## multiply edges by 100 for corPagel
#modTrees<-vector("list",length(slicePhys))
#for (j in 1:length(slicePhys)){
#	modTree<-slicePhys[[j]]
#	modTree$edge.length<-modTree$edge.length*100
#	modTrees[[j]]<-modTree
#}

# Z-score STANDARDIZED from RAW
results<-vector("list",length(allCladeSetNames))
for (q in 1: length(allCladeSetNames)){
	res<-read.table(paste(bbone,"_sample100_",i,"_slice",allCladeSetNames[[q]],"cladeSTATS_withTraits.txt",sep=""), header=TRUE)
	vars<-cbind.data.frame(res[,c(5,16:39)]) # 
	varsScale<-scale(vars, center=TRUE,scale=TRUE)
	results[[q]]<-cbind(res$DR_harm,log(res$DR_harm),res$BD_Div,log(res$BD_Div),res$BD_Lam,res$BD_Mu,(res$BD_Mu/res$BD_Lam),log(res$richness),varsScale)
	colnames(results[[q]])<-c("tipDR_mean","logTipDR","BD_Div","logBD_Div","BD_Lam","BD_Mu","BD_Turn","logRichness",colnames(res[,c(5,16:39)]))
}

###
# Predictors of ** tip DR mean -- or diversification rate in general **
##
# UNIVARIATE first
sliceTimes<-seq(-5,-70,-5)

predictors<-colnames(res[,c(5,16:39)])#c("MRCA", "logBodyMass", "trophicIndex1", "trophicIndex2", "perVert", "perInvert", "perPlant", "numLifemode", "numActivity", "ArboOrNot","SubtOrNot","DiurOrNot")#"majorActivity") #"majorLifemode"
nCols<-length(predictors)
#predictors<-c("tipDR_mean","BD_Div","BD_Lam","BD_Mu")
#nCols<-length(predictors)

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

# MULTIVARIATE
#sliceTimes<-seq(-5,-70,-5)
#
#Cols<-c("int","CarnOrNot", "HerbOrNot", "SE1","SE2","SE3","Lam","AIC","Pval1","Pval2")
#partSlopesPGLS_trophic2<-data.frame(matrix(NA, nrow = length(results), ncol = length(Cols)),row.names=sliceTimes)
#colnames(partSlopesPGLS_trophic2)<-Cols
#
#Cols<-c("int","DiurOrNot", "NoctOrNot","SE1","SE2","SE3","Lam","AIC","Pval1","Pval2")
#partSlopesPGLS_activity2<-data.frame(matrix(NA, nrow = length(results), ncol = length(Cols)),row.names=sliceTimes)
#colnames(partSlopesPGLS_activity2)<-Cols
#
#Cols<-c("int","SubtOrNot", "AquaOrNot", "TerrOrNot", "ArboOrNot", "SE1","SE2","SE3","SE4","SE5","Lam","AIC","Pval1","Pval2","Pval3","Pval4")
#partSlopesPGLS_lifemode4<-data.frame(matrix(NA, nrow = length(results), ncol = length(Cols)),row.names=sliceTimes)
#colnames(partSlopesPGLS_lifemode4)<-Cols

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

#} # cycles each time slice

#corr<-"BROWNIAN"
#corr<-"PAGEL"
corr<-"PAGEL-traitsAsResponse"
#corr<-"NO_TREE"
sliceN<-"45Ma-only" #"10to60Ma" # # #"35Ma-only" ##"all" #"60Ma"#"all" #
#resp<-"TipDR" 
#resp<-"BDdiv"
#resp<-"BDlam"
#resp<-"BDmu"

#uniPGLS<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs,uniPGLS_allSEs,uniPGLS_allLams)
#write.table(uniPGLS,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining",resp[z],"_timeSlices_",sliceN,"_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))

#write.table(partSlopesPGLS_trophic2,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining",resp,"_timeSlices_",sliceN,"_SCALED_trophic2.txt",sep=""))
#write.table(partSlopesPGLS_activity2,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining",resp,"_timeSlices_",sliceN,"_SCALED_activity2.txt",sep=""))
#write.table(partSlopesPGLS_lifemode4,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining",resp,"_timeSlices_",sliceN,"_SCALED_lifemode4.txt",sep=""))

uniPGLS[z,]<-cbind(uniPGLS_allInts,uniPGLS_allSlopes,uniPGLS_allPs,uniPGLS_allSEs,uniPGLS_allLams)[j,]

} # 4 response vars

write.table(uniPGLS,paste(bbone,"_sample100_",i,"_PGLS_",corr,"_explaining8_timeSlices_",sliceN,"_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))

} # cycles the 100 trees




#####
# Now want to be able to SUMMARIZE that one time slice 10 Ma... 

#===================================
# SUMMARIZE
# Load data back in and COMBINE all slice clades per-tree, then across all trees
# Post doing PGLS
#===================================
#setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/_ExplainingTipDR_timeSliceRuns")
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_ExplainingTipDR_timeSliceRuns")

bbone<- "NDexp" #"FBD" # 
library(nlme); library(MASS); library(AICcmodavg); library(ape); library(phytools)

numTrees<-100
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery
sliceTimes<-seq(-5,-70,-5)

allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

##
# MULTIVAR PENDANTS
#####
# 10 Ma or 35 Ma slice
######
library(plotrix)
# load data

model<-"trophic2"
#model<-"activity2"
#model<-"lifemode4"

slice<-c(NA,"10Ma",NA,NA,NA,NA,"35Ma",NA,NA,NA,NA,"60Ma",NA,NA)
#resp<-"TipDR"
resp<-"BDdiv" 

multiDat_All<-vector("list",length(slice))

#for(j in c(2,7)){

j=7
	allDat_MULTI_i<-vector("list",length=numTrees)
	for(i in 1:numTrees){
		tree<-i
		dat<-read.table(paste(bbone,"_sample100_",i,"_PGLS_PAGEL_explaining",resp,"_timeSlices_",slice[j],"_SCALED_",model,".txt",sep=""))
		allDat_MULTI_i[[i]]<-cbind(dat[j,],tree)
	}
	multiDat<-do.call(rbind,allDat_MULTI_i)

multiDat_All[[j]]<-multiDat
}

# SEs to calc 95 CIs...
###
#j=2 # 10 Ma time slice
j=7 # 35 Ma time slice
#j=12 # 60 Ma time slice

dat<-multiDat_All[[j]]
#slopesToPlot<-dat[,2:5]
#pvalsToPlot<-dat[,13:16]
#SEsToPlot<-dat[,7:10]
slopesToPlot<-dat[,2:3]
pvalsToPlot<-dat[,9:10]
SEsToPlot<-dat[,5:6]

empirPointCol<-grey(0.3,alpha=0.3)
redCol<-rgb(1,0,0,alpha=0.2)

# The 95% CI for each estimate of each var, plotted together
ALLdata_ToPlot_ALL<-vector("list",length(slopesToPlot))

for (i in 1:length(slopesToPlot[,1])){
	slopes_i<-slopesToPlot[i,]
	SEs_i<-SEsToPlot[i,]
	pVals_i<-pvalsToPlot[i,]
	ALLdata_ToPlot_i<-vector("list",length(slopesToPlot))
	for (j in 1:length(slopesToPlot)){
		color<-if(pVals_i[,j] < 0.05){ empirPointCol } else { redCol } 
		lowHigh95_j<-cbind.data.frame(slopes_i[,j]-(1.96*SEs_i[,j]),slopes_i[,j],slopes_i[,j]+(1.96*SEs_i[,j]),color,j)
		colnames(lowHigh95_j)<-c("low","mean","high","color","xVal")
		ALLdata_ToPlot_i[[j]]<-lowHigh95_j
	}
	ALLdata_ToPlot_ALL[[i]]<-ALLdata_ToPlot_i
}
varNames<-colnames(slopesToPlot)


####
# MULTIVAR PLOT with 2 variables and ERROR BARS-- at the 10 Ma and 35 Ma time slice...
###

sliceLabs<-c("10Ma", "35Ma", "60Ma")
q<-2
height<-0.1

#layout(matrix(c(1:3), 1, 3, byrow = TRUE), widths=c(2,2,2), heights=2)
pdf(file=paste("cladeLevel",bbone,"_PGLS_MULTI_PAGEL_explaining",resp,"_timeSlices_",sliceLabs[q],"_SCALED_",model,"_full95CI-fromSE.pdf",sep=""),onefile=TRUE, width=length(ALLdata_ToPlot_ALL[[1]])+3,height=5)
#par(oma = rep(5,4) + 0.1, mar = rep(0.5,4) + 0.1)

LwdCI<-1.5
pchs<-1 #c(rep(1,4),rep(2,4),rep(0,4)) # circle, triangle, square

par(oma = c(6,4,2,2) + 0.1, mar = c(2.5,0,0,0) + 0.1) # c(bottom, left, top, right)

#quartz(width=10,height=5)
dat<-ALLdata_ToPlot_ALL[[1]][[1]]
plotCI(x=dat$xVal, y=dat$mean,ui=dat$high,li=dat$low, cex=1.1,xlim=c(0.5,length(ALLdata_ToPlot_ALL[[i]])+0.4),ylim=c(-height,height), sfrac=0,xaxt="n", xlab="", ylab="Standardized effect on tip DR mean", err="y", lwd=LwdCI,col=as.character(dat$color),scol=as.character(dat$color),pch=pchs[1],font.lab=2,cex.axis=0.9,cex.lab=1.1)
abline(h=0,lty=1, lwd=2, col=grey(0.6, alpha=0.5)) #rgb(1,0,0,alpha=0.5))#

for(i in 1:length(ALLdata_ToPlot_ALL)){
	for(j in 1:length(ALLdata_ToPlot_ALL[[i]])){
		dat<-ALLdata_ToPlot_ALL[[i]][[j]]
		plotCI(x=dat$xVal, add=TRUE, y=dat$mean,ui=dat$high,li=dat$low, cex=1.1,xlim=c(1,15.4),ylim=c(-0.5,2), sfrac=0,xaxt="n", xlab="", ylab="", err="y", lwd=LwdCI,col=as.character(dat$color),scol=as.character(dat$color),pch=pchs[1],font.lab=2,cex.axis=1.1,cex.lab=1.1)
	}
}
axis(1, at=c(1:length(ALLdata_ToPlot_ALL[[i]]))+0.1, labels = FALSE)
text(x= c(1:length(ALLdata_ToPlot_ALL[[i]]))+0.1, y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), cex=1.1,font=2,srt = 45, labels = varNames, xpd = NA, adj=c(0.95,0.05))
#text(x= c(1:length(ALLdata_ToPlot_ALL[[i]]))+0.1, y = -0.05, cex=1.1,font=1, labels = round(lamMeans, 2), xpd = TRUE)#, adj=c(0.95,0.05))

mtext(side=2, font=2, text=paste("Standardized effect on ",resp,sep=""), line=2)

mtext(side=3, font=2, text=paste("PGLS Pagel: ",sliceLabs[q]," slice",sep=""), adj=0)
#mtext(side=3, font=2, text="PGLS Pagel: 35 Ma slice", adj=0)
#mtext(side=3, font=2, text="PGLS Pagel: 60 Ma slice", adj=0)

dev.off()





##
# UNIVAR PENDANTS
#####
# 10 Ma or 35 Ma slice
######
library(plotrix)

# load data
#q=2 # 10 Ma time slice
#q=7 # 35 Ma time slice
q=9 # 45 Ma slice
#q=12 # 60 Ma time slice
slice<-c(NA,"10Ma",NA,NA,NA,NA,"35Ma-only",NA,"45Ma-only",NA,NA,"60Ma",NA,NA)
#resp<-"BDdiv" #"TipDR"
#resp<-"TipDR"
#resp<-"BDlam"
#resp<-"BDmu"

respVar<-c("tipDR_mean","logTipDR","BD_Div","logBD_Div","BD_Lam","BD_Mu","BD_Turn","logRichness") #c("tipDR_mean","logTipDR","BD_Div","logBD_Div","BD_Lam","BD_Mu","BD_Turn","logRichness") #c("tipDR_mean","BD_Div","BD_Lam","BD_Mu","BD_Turn") #,"richness") #

numTrees<-c(2, 6, 10, 11, 12, 21, 23, 24, 25, 27, 29, 31, 33, 40, 41, 48, 53, 54, 59, 60, 63, 65, 67, 68, 72, 73, 74, 77, 82, 83, 84, 85, 86, 91, 94, 96, 98)
#numTrees<-c(2:100)
#numTrees<-c(1, 2, 3, 6, 7, 8, 10, 11, 12, 21, 23, 24, 25, 27, 29, 30, 31, 33, 35, 40, 41, 42, 45, 48, 53, 54, 59, 60, 61, 63, 65, 67, 68, 72, 73, 74, 75, 77, 78, 79, 82, 83, 84, 85, 86, 90, 91, 94, 96, 98)
#numTrees<-c(12, 13, 15, 16, 17, 19, 20, 22, 24, 25, 26, 28, 29, 30, 31, 32, 33, 36, 39, 41, 1, 2, 4, 5, 6, 7, 8, 9, 11, 14, 21, 27, 34, 35, 37, 42, 43, 44, 46, 48, 58, 3, 10, 23, 47, 49, 50, 51, 52, 55, 56, 61, 62, 63, 67, 71, 73, 54, 57, 59, 60, 65, 66, 69, 72, 74, 76, 77, 78, 81, 86, 68, 75, 80, 83, 84, 85, 87, 88, 89, 90, 91, 95, 97, 100, 96, 98, 99)

allDat_UNI_varALL<-vector("list",length(respVar))
for(j in 1:length(respVar)){
	
	allDat_UNI_varTrees<-vector("list",length(numTrees))
	for(i in numTrees){
		tree<-i
		dat<-read.table(paste(bbone,"_sample100_",i,"_PGLS_PAGEL-traitsAsResponse_explaining8_timeSlices_",slice[q],"_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))
		allDat_UNI_varTrees[[i]]<-(cbind(dat[j,],tree))
	}
	res_j<-do.call(rbind,allDat_UNI_varTrees)
	allDat_UNI_varALL[[j]]<-res_j
}


#uniDat_means<-as.data.frame(colMeans(uniDat))
#lamMeans<-uniDat_means[c(paste("lam_",1:20,sep="")),]

#predictors<-c("percentSamp", "logBodyMass_kg", "logGeoRange_km2", "logHomeRange_km2", "rankHomeRange_1to10", "rankHomeRange_1to10_binsEq", "logDispDist_km", "rankDispDist_1to10", "rankDispDist_1to10_binsEq", "logGenLength_d", "trophic123", "CarnOrNot", "HerbOrNot", "OmniOrNot", "perVert", "perInvert", "perPlant", "AquaOrNot", "ArboOrNot", "FlysOrNot", "SubtOrNot", "TerrOrNot", "NoctOrNot", "CathOrNot", "DiurOrNot")
#nCols<-length(predictors)
#write.table(predictors,file="predictors_names.txt")
predictors<-read.table(file="predictors_names.txt")

ALLdata_ToPlot_ALL_perResp<-vector("list",length(respVar))
for (z in 1:length(respVar)){
uniDat<-allDat_UNI_varALL[[z]]

# SEs to calc 95 CIs...
slopesToPlot<-uniDat[,c(as.vector(predictors$x))]
pvalsToPlot<-uniDat[,c(paste("p_",1:length(predictors[,1]),sep=""))]
SEsToPlot<-uniDat[,c(paste("SE_",1:length(predictors[,1]),sep=""))]

empirPointCol<-grey(0.3,alpha=0.3)
redCol<-rgb(1,0,0,alpha=0.05)

# The 95% CI for each estimate of each var, plotted together
ALLdata_ToPlot_ALL<-vector("list",length(slopesToPlot))

for (i in 1:length(slopesToPlot[,1])){
	slopes_i<-slopesToPlot[i,]
	SEs_i<-SEsToPlot[i,]
	pVals_i<-pvalsToPlot[i,]
	ALLdata_ToPlot_i<-vector("list",length(slopesToPlot))
	for (j in 1:length(slopesToPlot)){
		if(j==11){ next } else {
		color<-if(pVals_i[,j] <= 0.05){ empirPointCol } else { redCol } 
		lowHigh95_j<-cbind.data.frame(slopes_i[,j]-(1.96*SEs_i[,j]),slopes_i[,j],slopes_i[,j]+(1.96*SEs_i[,j]),color,j)
		colnames(lowHigh95_j)<-c("low","mean","high","color","xVal")
		ALLdata_ToPlot_i[[j]]<-lowHigh95_j
		}
	}
	ALLdata_ToPlot_ALL[[i]]<-ALLdata_ToPlot_i
}
varNames<-colnames(slopesToPlot)

ALLdata_ToPlot_ALL_perResp[[z]]<-ALLdata_ToPlot_ALL

}

# get the num per var signif-- for each of the 8 respVars

sumOfSignif_ALL_respVars<-vector("list",length(respVar))
for(z in 1:length(respVar)){
ALLdata_ToPlot_ALL<-ALLdata_ToPlot_ALL_perResp[[z]]

sumOfSignif_ALL<-vector("list",length(predictors$x))
sumOfSignif_ALL[[11]]<-NA
for(j in 1:length(predictors$x)){
	if(j==11){ next } else {
	sumOfSignif<-c()
	for(k in 1:(length(numTrees))){
		dat<-ALLdata_ToPlot_ALL[[k]][[j]]
		if(dat$color==redCol){sumOfSignif[k]<-0} else {sumOfSignif[k]<-1}
		}
	SUM<-sum(sumOfSignif)
	}	
	sumOfSignif_ALL[[j]]<-round(((SUM/length(numTrees))*100),0)
	}
sumOfSignif_ALL_respVars[[z]]<-sumOfSignif_ALL
}



####
# UNIVAR PLOT with 25 variables and ERROR BARS-- for 8 different RespVars... at SPECIFIED time slices...
###

height<-2

respVar<-c("tipDR_mean","logTipDR","BD_Div","logBD_Div","BD_Lam","BD_Mu","BD_Turn","logRichness") #c("tipDR_mean","logTipDR","BD_Div","logBD_Div","BD_Lam","BD_Mu","BD_Turn","logRichness") #c("tipDR_mean","BD_Div","BD_Lam","BD_Mu","BD_Turn") #,"richness") #

for(z in 1:length(respVar)){

ALLdata_ToPlot_ALL<-ALLdata_ToPlot_ALL_perResp[[z]]
sumOfSignif_ALL<-unlist(sumOfSignif_ALL_respVars[[z]])

#layout(matrix(c(1:3), 1, 3, byrow = TRUE), widths=c(2,2,2), heights=2)
pdf(file=paste("cladeLevel",bbone,"_PGLS_UNI_PAGEL_explaining",respVar[z],"_timeSlices_",slice[q],"_SCALED_full95CI-fromSE_25vars.pdf",sep=""),onefile=TRUE, width=15,height=6)
#pdf(file=paste("cladeLevel",bbone,"_PGLS_UNI_PAGEL_explaining",resp,"_timeSlices_",slice[q],"_SCALED_full95CI-fromSE_20vars.pdf",sep=""),onefile=TRUE, width=15,height=6)
#par(oma = rep(5,4) + 0.1, mar = rep(0.5,4) + 0.1)

LwdCI<-1.5
pchs<-1 #c(rep(1,4),rep(2,4),rep(0,4)) # circle, triangle, square

par(oma = c(6,4,2,2) + 0.1, mar = c(2.5,0,0,0) + 0.1) # c(bottom, left, top, right)

#quartz(width=10,height=5)
dat<-ALLdata_ToPlot_ALL[[1]][[1]]
plotCI(x=dat$xVal, y=dat$mean,ui=dat$high,li=dat$low, cex=1.1,xlim=c(1,length(ALLdata_ToPlot_ALL[[i]])+0.4),ylim=c(-height,height), sfrac=0,xaxt="n", xlab="", ylab=paste("Standardized effect on ",respVar[z],sep=""), err="y", lwd=LwdCI,col=as.character(dat$color),scol=as.character(dat$color),pch=pchs[1],font.lab=2,cex.axis=0.9,cex.lab=1.1)
abline(h=0,lty=1, lwd=2, col=grey(0.6, alpha=0.5)) #rgb(1,0,0,alpha=0.5))#

for(i in 1:length(ALLdata_ToPlot_ALL)){
	for(j in 1:length(ALLdata_ToPlot_ALL[[i]])){
		if(j==11){ next } else {
		dat<-ALLdata_ToPlot_ALL[[i]][[j]]
		plotCI(x=dat$xVal, add=TRUE, y=dat$mean,ui=dat$high,li=dat$low, cex=1.1,sfrac=0,xaxt="n", xlab="", ylab="", err="y", lwd=LwdCI,col=as.character(dat$color),scol=as.character(dat$color),pch=pchs[1],font.lab=2,cex.axis=1.1,cex.lab=1.1)
		}
	}
}
axis(1, at=c(1:length(ALLdata_ToPlot_ALL[[i]]))+0.1, labels = FALSE)
text(x= c(1:length(ALLdata_ToPlot_ALL[[i]]))+0.1, y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), cex=1.1,font=2,srt = 45, labels = varNames, xpd = NA, adj=c(0.95,0.05))
text(x= c(1:length(ALLdata_ToPlot_ALL[[i]]))+0.1, y = -height, cex=1.3,font=2, labels = sumOfSignif_ALL, xpd = TRUE)#, adj=c(0.95,0.05))

mtext(side=2, font=2, text=paste("Standardized effect on ",respVar[z],sep=""), line=2)

mtext(side=3, font=2, text=paste("PGLS Pagel: predictors ~ ",respVar[z]," for ",slice[q]," slice",sep=""), adj=0)

dev.off()

}



####
# UNIVAR from 10-60 Ma across each of 20 predictors of DR
##
####
# co-opting code from before...
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_ExplainingTipDR_timeSliceRuns")

bbone<- "NDexp" #"FBD" # 
library(nlme); library(MASS); library(AICcmodavg); library(ape); library(phytools)

numTrees<-100
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery
sliceTimes<-seq(-5,-70,-5)

allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

# load data
slice<-"10to60Ma"
numTrees<-c(1:100)

respVar<-c("logTipDR_mean","logRichness","BD_Turn","BD_Lam") #c("tipDR_mean","BD_Div","BD_Lam","BD_Mu")
resp<-c("logTipDR","logRichness","BDturn", "BDlam")
#respVar<-c("tipDR_mean","logTipDR_mean","BD_Div","logBD_Div","logRichness","BD_Turn","BD_Lam","BD_Mu") #c("tipDR_mean","BD_Div","BD_Lam","BD_Mu")
#resp<-c("TipDR", "LogTipDR","BDdiv", "LogBDdiv","LogRichness","BDturn", "BDlam", "BDmu")
#respVar<-c("logTipDR_mean","logBD_Div","logRichness") #c("tipDR_mean","BD_Div","BD_Lam","BD_Mu")
#resp<-c("LogTipDR", "LogBDdiv", "LogRichness") #"BDlam", "BDmu")

model<-"PAGEL"
#model<-"NO_TREE"

pdf(file=paste("cladeLevel",bbone,"_PGLSuni_",model,"_explaining_4resp-logDR-richness-BD_timeSlices_",slice,"_values_andSE-BAR_25vars.pdf",sep=""),onefile=TRUE, width=(4*4),height=(4*25))
#pdf(file=paste("cladeLevel",bbone,"_PGLSuni_",model,"_explaining_3resp-wLogVars_timeSlices_",slice,"_values_andSE-BAR_25vars.pdf",sep=""),onefile=TRUE, width=(4*4),height=(4*25))
#pdf(file=paste("cladeLevel",bbone,"_PGLSuni_",model,"_explaining_8resp-wLogVarsTurnover_timeSlices_",slice,"_values_andSE-BAR_25vars.pdf",sep=""),onefile=TRUE, width=(4*8),height=(4*25))

#quartz(width=(6*4),height=(4*25))
layout(matrix(c(1:100), 25, 4, byrow = FALSE), widths=rep(4,100), heights=rep(3,100))
par(oma = c(5,4,5,3) + 0.1, mar = c(4,1,1,1) + 0.1)

for(z in 1:length(respVar)){

#numTrees<-100
sliceNums<-c(1:14)
allDat_UNI_i<-vector("list",length(numTrees))
#for(i in 1:numTrees){
#ranTrees<-c(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27, 28, 30, 31, 32, 33) # BD_mu
#ranTrees<-c(1:49,51:100) # BD_div
ranTrees<-c(1:100) # BD_lam
#ranTrees<-c(2, 7, 8, 16, 18, 5, 11, 15, 22, 23, 30, 35, 36, 37, 40, 41, 42, 45, 1, 3, 4, 6, 9, 12, 13, 14, 19, 20, 26, 27, 31, 46, 48, 17, 21, 25, 28, 33, 34, 38, 39, 43, 44, 47, 53, 56, 64, 66, 54, 59, 68, 71, 73, 75, 79, 90, 24, 52, 62, 65, 74, 76, 77, 80, 83, 84, 85, 97, 51, 58, 61, 70, 78, 81, 82, 89, 91, 93, 98, 57, 60, 63, 88, 92, 94, 95, 96, 100)
numTrees<-length(ranTrees)
for(i in ranTrees){
	tree<-i
	dat<-read.table(paste(bbone,"_sample100_",i,"_PGLS_",model,"_explaining",resp[z],"_timeSlices_",slice,"_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))
	allDat_UNI_i[[i]]<-cbind(sliceTimes,dat,sliceTimes,tree)
}
uniDat_all<-do.call(rbind,allDat_UNI_i)
#uniDat_all[,c(paste("lam_",c(1:25),sep=""))]<-list(NULL)
uniDat_all[,c("i_11","trophic123","p_11","SE_11","lam_11")]<-list(NULL)
uniDat<-na.omit(uniDat_all)
write.table(uniDat,file=paste(bbone,"_sample100_ALL100-TABLE_",model,"_explaining",resp[z],"_timeSlices_",slice,"_SCALED_uniINTS_uniSLO_uniP.txt",sep=""))

predictorNames<-read.table(file="predictors_names.txt")
predictors<-as.vector(predictorNames$x)[-11]

# SEs to calc 95 CIs...
slopesToPlot<-uniDat[,predictors]
pvalsToPlot<-uniDat[,c(paste("p_",(1:25)[-11],sep=""))]
SEsToPlot<-uniDat[,c(paste("SE_",(1:25)[-11],sep=""))]
sliceNumsToPlot<-uniDat[,"sliceTimes"]

# The 95% CI for each estimate of each var, plotted together
ALLdata_ToPlot_ALL<-vector("list",length(slopesToPlot))

for (i in 1:length(slopesToPlot[,1])){ # so for a given slice / tree... taking each of 20 vars and calc 95%CI
	slopes_i<-slopesToPlot[i,]
	SEs_i<-SEsToPlot[i,]
	pVals_i<-pvalsToPlot[i,]
	ALLdata_ToPlot_i<-vector("list",length(slopesToPlot))
	for (j in 1:length(slopesToPlot)){
		color<-if(pVals_i[,j] <= 0.05){ "grey" } else { "red" } 
		lowHigh95_j<-cbind.data.frame(slopes_i[,j]-(1.96*SEs_i[,j]),slopes_i[,j],slopes_i[,j]+(1.96*SEs_i[,j]),color,j,sliceNumsToPlot[i])
		colnames(lowHigh95_j)<-c("low","mean","high","color","xVal","slice")
		ALLdata_ToPlot_i[[j]]<-lowHigh95_j
	}
	ALLdata_ToPlot_ALL[[i]]<-ALLdata_ToPlot_i
}
varNames<-colnames(slopesToPlot)

# get the number of significant runs per time slice per variable across 100 trees
sumOfSignif_ALL<-vector("list",length(predictors))
for(i in 1:length(predictors)){
	sumOfSignif_slicesPerVar<-vector(length=(numSlices-3))
	for(k in 1:11){
		sumOfSignif<-c()
		for(j in 0:(numTrees-1)){
		q<-k+(j*11)
		dat<-ALLdata_ToPlot_ALL[[q]][[i]]
		if(dat$color=="red"){sumOfSignif[j]<-0} else {sumOfSignif[j]<-1}
		}
		sumOfSignif_slicesPerVar[k]<-sum(sumOfSignif)
	}
	sumOfSignif_ALL[[i]]<-round(((sumOfSignif_slicesPerVar/length(ranTrees))*100),0)
}


vertLwd<-0.5
lwdCI<-1.5
horizLine<-0
empirPointCol<-grey(0.3,alpha=0.8)
redCol<-rgb(1,0,0,alpha=0.3)
yLims1<-c(-0.5,0.5)

#respVar<-"Tip DR mean"
#respVar<-"BD div"
#respVar<-"BD lam"
#respVar<-"BD mu"


#pdf(file=paste("cladeLevel",bbone,"_PGLSuni_",model,"_explaining",resp[z],"_timeSlices_",slice,"_values_andSE-BAR_22vars.pdf",sep=""),onefile=TRUE, width=6,height=4)
#quartz(width=6,height=4)

#quartz(width=(4*4),height=(4*25))
#layout(matrix(c(1:100), 25, 4, byrow = FALSE), widths=rep(4,100), heights=rep(3,100))
#par(oma = c(5,4,5,3) + 0.1, mar = c(4,1,1,1) + 0.1)

for(i in 1:length(predictors)){

plot(as.formula(paste(predictors[i]," ~ sliceTimes")), data=uniDat, type="n",ylim=yLims1, ylab="",xlab="", yaxt="n",xaxt="n")
axis(side=2,at=NULL,labels=TRUE)
axis(side=1,at=NULL,labels=TRUE)

for(j in 1:numTrees){
	dat<-ALLdata_ToPlot_ALL[[j]][[i]]
	if(dat$color=="red"){col<-redCol} else {col<-empirPointCol}
	plotCI(x=dat$slice, add=TRUE,y=dat$mean,ui=dat$high,li=dat$low, cex=1.3,sfrac=0, err="y", lwd=lwdCI,col=col,scol=col,pch=1,font.lab=2,cex.axis=1.1,cex.lab=1.1)
}
mtext(side=3,text=paste(respVar[z]," ~ "),font=2,adj=0,line=1)
mtext(side=3,text=predictors[i],font=2,adj=0)
text(x=sliceTimes[2:12], y = yLims1[2], cex=1.4,font=2, labels = sumOfSignif_ALL[[i]], col="dark grey")

for(i in 2:(numSlices-2)){
	abline(v=(-5*i),lty=2, lwd=vertLwd, col=grey(0.6, alpha=0.5))
}
	abline(h=horizLine,lty=1, lwd=1, col=grey(0.6, alpha=0.5))

}

} # end 4 response var loop

dev.off()




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
























#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy v1 -- Upham et al. 2017
###
# Figure X - distributionns of clade richnness per time slice
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Visualize the EMPIRICAL time slice distributions
######
# - Load in richness counts per clade for each timeSlice (pre-calculated) 
# - 
# - 
# - 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# packages
library(moments); library(nlme); library(ape); library(picante); library(phytools); library(geiger)
library(plotrix)

# ~~~~~~~~~~~~
# EMPIRICAL vs SIMULATIONS 
# ===================

# directory and source
dirname="/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/cladeLevel_SLICES_explainingDR/slicesBy5Ma_EMPIRICAL"
setwd(dirname)

bbone<-"NDexp"
ntrees<-100

# Make time slices NAMES
#=======================================
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery

allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
}

# read in richness counts (saved as Rdata)
#=======================================
allCladeN_mamPhy<-vector("list",length=ntrees)
for(i in 1:ntrees){
	load(paste(bbone,"_sample100_",i,"_timeslice_cladeRichnesses_wSingletons.Rda",sep=""))
	allCladeN_mamPhy[[i]]<-allLengths_i
}

# Now on PER TREE basis, calculate expected N (mean) and variance of N (var) PER SLICE
#=======================================
sliceStatsALL<-vector("list",length=ntrees)
for(i in 1:ntrees){
	
	sliceStats_i<-vector("list",length(allLengths_i))
	for(j in 1:length(allLengths_i)){
	x<-allCladeN_mamPhy[[i]][[j]]
	stats<-c(mean(x), var(x), sd(x), std.error(x), quantile(x,0.5)[[1]], quantile(x,0.025)[[1]],quantile(x,0.975)[[1]],sum(x))
	names(stats)<-c("mean", "var", "sd", "stdErr", "median", "low95", "up95", "sum")
	sliceStats_i[[j]]<-stats
	} # end 14 slices
	names(sliceStats_i)<-allCladeSetNames
	sliceStatsALL[[i]]<-do.call(rbind, sliceStats_i)

} # end 100 trees
write.table(sliceStatsALL, file=paste("SUMMARY_",bbone,"_sample100_ALL100_timeslice_cladeRichnesses.txt",sep=""))


# for SIMS-- do the same read-in richness counts and summary (MamPhy-based extinction fraction)
# ======================================
# set dir
setwd("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/cladeLevel_SLICES_explainingDR/slicesBy5Ma_SIMS")
sims<-c("mamPhyE","lowE_0p2","highE_0p8")

sliceStatsALL_SIM_all3sims<-vector("list",length(sims))
for(z in 1:length(sims)){
#z=1
# read in
allCladeN_SIM<-vector("list",length=ntrees)
for(i in 1:ntrees){
	load(paste("MamPhy_SIMS_",sims[z],"_",bbone,"_sample100_",i,"_timeslice_cladeRichnesses_wSingletons.Rda",sep=""))
	allCladeN_SIM[[i]]<-allLengths_i
}

# summary
sliceStatsALL_SIM<-vector("list",length=ntrees)
for(i in 1:ntrees){
	
	sliceStats_i<-vector("list",length(allLengths_i))
	for(j in 1:length(allLengths_i)){
	x<-allCladeN_SIM[[i]][[j]]
	stats<-c(mean(x), var(x), sd(x), std.error(x), quantile(x,0.5)[[1]], quantile(x,0.025)[[1]],quantile(x,0.975)[[1]],sum(x))
	names(stats)<-c("mean", "var", "sd", "stdErr", "median", "low95", "up95", "sum")
	sliceStats_i[[j]]<-stats
	} # end 14 slices
	names(sliceStats_i)<-allCladeSetNames
	sliceStatsALL_SIM[[i]]<-do.call(rbind, sliceStats_i)

} # end 100 trees
write.table(sliceStatsALL_SIM, file=paste("SUMMARY_","MamPhy_SIMS_",sims[z],"_",bbone,"_sample100_ALL100_timeslice_cladeRichnesses.txt",sep=""))
sliceStatsALL_SIM_all3sims[[z]]<-sliceStatsALL_SIM

} # end 3 sim loop
#########

# graph the clade sizes per slice, EMPIRICAL vs SIM
# =============
library(plotrix)
#sliceStatsALL<-read.table(file=paste("SUMMARY_",bbone,"_sample100_ALL100_timeslice_cladeRichnesses.txt",sep=""))
#sliceStatsALL_SIM<-read.table(file=paste("SUMMARY_","MamPhy_SIMS_",sims[z],"_",bbone,"_sample100_ALL100_timeslice_cladeRichnesses.txt",sep=""))

LwdCI<-1.5
greyAlpha<-0.2
empirPointCol<-grey(0.3,alpha=0.3)

simAlpha<-0.2
Col<-"darkgoldenrod2"
Col1<-rgb((col2rgb(Col)/255)[1],(col2rgb(Col)/255)[2],(col2rgb(Col)/255)[3],alpha=simAlpha)
Col<-"deepskyblue2"
Col2<-rgb((col2rgb(Col)/255)[1],(col2rgb(Col)/255)[2],(col2rgb(Col)/255)[3],alpha=simAlpha)
Col<-"darkorchid3"
Col3<-rgb((col2rgb(Col)/255)[1],(col2rgb(Col)/255)[2],(col2rgb(Col)/255)[3],alpha=simAlpha)
ColorsSIM<-c(Col1,Col2,Col3)
ColorsSIM_full<-c("darkgoldenrod2","deepskyblue2","darkorchid3")

#pdf(file="Expected_CladeN_acrossTimeSlices_mamPhy-and-1SIM_nonLOG.pdf",height=5, width=5)
#pdf(file="Expected_CladeN_acrossTimeSlices_mamPhy-and-1SIM_LOG.pdf",height=5, width=5)
#pdf(file="Expected_CladeN_acrossTimeSlices_mamPhy-vs-3SIM_LOG.pdf",height=5, width=5)
#pdf(file="Expected_CladeN_acrossTimeSlices_mamPhy-vs-3SIM_nonLOG.pdf",height=5, width=5)
#pdf(file="Expected_CladeN_acrossTimeSlices_mamPhy-vs-3SIM_LOG_stDevs.pdf",height=5, width=5)
pdf(file="Expected_CladeN_acrossTimeSlices_mamPhy-vs-3SIM_LOG_variances.pdf",height=5, width=5)

dat<-sliceStatsALL[[1]]
yMax<-100000
plot(x=rev(1:14), y=dat[,"mean"],col=empirPointCol, ylim=c(1,yMax),xaxt="n", xlab="Time before present (Ma)", ylab="Variance of clade richness", log="y")

# EMPIRICAL
#plotCI(x=rev(1:14), y=dat[,"mean"], liw=(dat[,"mean"]-dat[,"stdErr"]), uiw=(dat[,"mean"]+dat[,"stdErr"]), sfrac=0, ylim=c(0,1500), lwd=LwdCI,col=empirPointCol,scol=grey(0.3,alpha=greyAlpha),pch=1)
for(i in 1:ntrees){
	dat<-sliceStatsALL[[i]]
	#plotCI(add=TRUE, x=rev(1:14), y=dat[,"mean"], liw=(dat[,"mean"]-dat[,"sd"]), uiw=(dat[,"mean"]+dat[,"sd"]), sfrac=0, lwd=LwdCI,col=empirPointCol,scol=grey(0.3,alpha=greyAlpha),pch=1)
	points(x=rev(1:14), y=dat[,"var"],col=empirPointCol,pch=1)
}

# SIMS
for(z in 1:length(sims)){
sliceStatsALL_SIM<-sliceStatsALL_SIM_all3sims[[z]]

	for(i in 1:ntrees){
		dat<-sliceStatsALL_SIM[[i]]
		#plotCI(add=TRUE, x=rev(1:14), y=dat[,"mean"], liw=(dat[,"mean"]-dat[,"sd"]), uiw=(dat[,"mean"]+dat[,"sd"]), sfrac=0, lwd=LwdCI,col=ColorsSIM[z],scol=ColorsSIM[z],pch=1)
		points(x=rev(1:14), y=dat[,"var"],col=ColorsSIM[z],pch=1)
	}
}

axis(1, at=rev(1:14), labels = seq(from=5,to=70,by=5))

legend(x=10, y=yMax, legend=c("empirical","RC 0.65", "RC 0.2","RC 0.8"),fill=c(grey(0.3),ColorsSIM_full))

dev.off()









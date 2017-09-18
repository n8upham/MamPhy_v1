#######
# STRAPP analyses... on the 10 trees
###
library(BAMMtools); library(coda); library(phytools); library(ape); library(phangorn)

folders<-c(1:10)

#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/grace_2.5_1week_NDexp_1")

## BODY MASS
bbone = "NDexp" # "FBD" 
res<-vector("list",length(folders))
resAll<-data.frame(matrix(NA,nrow=length(folders),ncol=2))
obsCorrs<-vector("list",length(folders))
nullDists<-vector("list",length(folders))

for (i in 1:length(folders)){

#setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/grace_2.5_1week_",dir,folders[i],sep=""))
setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/grace_2.5_1week_NDexp_",folders[i],sep=""))

# load edata
tree <- ladderize(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample10_",folders[i],".tre",sep="")))
edata <- getEventData(tree, eventdata = paste("mamPhy_",bbone,"_event_data.txt",sep=""), burnin=0.33, nsamples=1000)

# load traits
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/_NDexp_STRAPP_analyses")

tipDataAll_wComments<-read.table(file="MamPhy_5911sp_tiplabel_DR-range-mass-herb-lifemode.txt", header=TRUE)
tipDataAll<-tipDataAll_wComments[,5:9]
rownames(tipDataAll)<-tipDataAll_wComments$tiplabel

BMall<-log(tipDataAll[,"bodyMass"])
names(BMall)<-rownames(tipDataAll)

# run STRAPP
res[[i]]<-traitDependentBAMM(ephy=edata, traits=BMall, reps=1000, rate = "speciation", return.full = TRUE, method = "spearman", logrates = TRUE, two.tailed = TRUE, traitorder = NA, nthreads = 1)
resAll[i,1]<-res[[i]][1]
resAll[i,2]<-res[[i]][2]
obsCorrs[[i]]<-res[[i]][6]
nullDists[[i]]<-res[[i]][8]
}

colnames(resAll)<-c("estimate","pVal")
write.table(resAll,file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-bodyMass_ALL.txt")
write.table(obsCorrs,file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-bodyMass_ALL_obsCorr.txt")
write.table(nullDists,file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-bodyMass_ALL_nullDists.txt")

## GEO RANGE
bbone = "NDexp" # "FBD" 
res<-vector("list",length(folders))
resAll<-data.frame(matrix(NA,nrow=length(folders),ncol=2))
obsCorrs<-vector("list",length(folders))
nullDists<-vector("list",length(folders))

for (i in 1:length(folders)){

#setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/grace_2.5_1week_",dir,folders[i],sep=""))
setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/grace_2.5_1week_NDexp_",folders[i],sep=""))

# load edata
tree <- ladderize(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample10_",folders[i],".tre",sep="")))
edata <- getEventData(tree, eventdata = paste("mamPhy_",bbone,"_event_data.txt",sep=""), burnin=0.33, nsamples=1000)
#cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
#subtreeBAMM(edata, tips=)
#>>>

# load traits
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/_NDexp_STRAPP_analyses")

tipDataAll_wComments<-read.table(file="MamPhy_5911sp_tiplabel_DR-range-mass-herb-lifemode.txt", header=TRUE)
tipDataAll<-tipDataAll_wComments[,5:9]
rownames(tipDataAll)<-tipDataAll_wComments$tiplabel

GEOall<-log(tipDataAll[,"geoRange"])
names(GEOall)<-rownames(tipDataAll)

# run STRAPP
res[[i]]<-traitDependentBAMM(ephy=edata, traits=GEOall, reps=1000, rate = "speciation", return.full = TRUE, method = "spearman", logrates = TRUE, two.tailed = TRUE, traitorder = NA, nthreads = 1)
resAll[i,1]<-res[[i]][1]
resAll[i,2]<-res[[i]][2]
obsCorrs[[i]]<-res[[i]][6]
nullDists[[i]]<-res[[i]][8]
}

colnames(resAll)<-c("estimate","pVal")
#write.table(resAll,file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-geoRange_ALL.txt")
#write.table(obsCorrs,file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-geoRange_ALL_obsCorr.txt")
#write.table(nullDists,file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-geoRange_ALL_nullDists.txt")


## LIFEMODE
library(BAMMtools); library(coda); library(phytools); library(ape); library(phangorn)
folders<-c(1:10)

bbone = "NDexp" # "FBD" 
res<-vector("list",length(folders))
resAll<-data.frame(matrix(NA,nrow=length(folders),ncol=2))
obsCorrs<-vector("list",length(folders))
nullDists<-vector("list",length(folders))

for (i in 1:length(folders)){

#setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/grace_2.5_1week_",dir,folders[i],sep=""))
setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/grace_2.5_1week_NDexp_",folders[i],sep=""))

# load edata
tree <- ladderize(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample10_",folders[i],".tre",sep="")))
edata <- getEventData(tree, eventdata = paste("mamPhy_",bbone,"_event_data.txt",sep=""), burnin=0.33, nsamples=1000)
#cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
#subtreeBAMM(edata, tips=)
#>>>

# load traits
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/_NDexp_STRAPP_analyses")

tipDataAll_wComments<-read.table(file="MamPhy_5911sp_tiplabel_DR-range-mass-herb-lifemode.txt", header=TRUE)
tipDataAll<-tipDataAll_wComments[,5:9]
rownames(tipDataAll)<-tipDataAll_wComments$tiplabel

LMall<-tipDataAll[,"LIFEMODE"]
names(LMall)<-rownames(tipDataAll)

# run STRAPP
res[[i]]<-traitDependentBAMM(ephy=edata, traits=LMall, reps=1000, rate = "speciation", return.full = TRUE, method = "kruskal", logrates = TRUE, two.tailed = TRUE, traitorder = NA, nthreads = 1)
resAll[i,1]<-res[[i]][1]
resAll[i,2]<-res[[i]][2]
obsCorrs[[i]]<-res[[i]][6]
nullDists[[i]]<-res[[i]][8]
}

colnames(resAll)<-c("estimate","pVal")
write.table(resAll,file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-LifeMode_ALL.txt")
write.table(obsCorrs,file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-LifeMode_ALL_obsCorr.txt")
write.table(nullDists,file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-LifeMode_ALL_nullDists.txt")


## TROPHIC LEVEL
library(BAMMtools); library(coda); library(phytools); library(ape); library(phangorn)
folders<-c(1:10)

bbone = "NDexp" # "FBD" 
res<-vector("list",length(folders))
resAll<-data.frame(matrix(NA,nrow=length(folders),ncol=2))
obsCorrs<-vector("list",length(folders))
nullDists<-vector("list",length(folders))

for (i in 1:length(folders)){

#library(foreach);library(doSNOW)
#cl = makeCluster(10, type = 'SOCK', outfile="")
#registerDoSNOW(cl)
#ntrees = length(folders)
#foreach(i=1:ntrees, .packages=c('BAMMtools', 'coda', 'phytools', 'ape', 'phangorn'), .combine=cbind, .verbose=TRUE) %dopar% {

setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/grace_2.5_1week_NDexp_",folders[i],sep=""))
#setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/grace_2.5_1week_NDexp_",folders[i],sep=""))

# load edata
tree <- ladderize(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample10_",folders[i],".tre",sep="")))
edata <- getEventData(tree, eventdata = paste("mamPhy_",bbone,"_event_data.txt",sep=""), burnin=0.33, nsamples=1000)
#cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
#subtreeBAMM(edata, tips=)
#>>>

# load traits
setwd("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/")
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/_NDexp_STRAPP_analyses")

tipDataAll_wComments<-read.table(file="MamPhy_5911sp_tiplabel_DR-range-mass-herb-lifemode.txt", header=TRUE)
tipDataAll<-tipDataAll_wComments[,5:9]
rownames(tipDataAll)<-tipDataAll_wComments$tiplabel

trophic_cats<-c("Herbivore","Omnivore","Carnivore")

trophicSpecies<-data.frame(matrix(NA,nrow=5911,ncol=1))
trophicSpecies[which(tipDataAll$percentHerb==100),]<-"Herbivore"
trophicSpecies[which(tipDataAll$percentHerb <= 90 & tipDataAll$percentHerb >= 10 ),]<-"Omnivore"
trophicSpecies[which(tipDataAll$percentHerb==0),]<-"Carnivore"

trophicALL<-cbind.data.frame(tipDataAll$percentHerb,trophicSpecies)
rownames(trophicALL)<-rownames(tipDataAll)
colnames(trophicALL)<-c("percentHerb","trophicCat")

trophicOnly<-as.factor(trophicALL[,2])
names(trophicOnly)<-rownames(tipDataAll)

# run STRAPP
res[[i]]<-traitDependentBAMM(ephy=edata, traits=trophicOnly, reps=1000, rate = "speciation", return.full = TRUE, method = "kruskal", logrates = TRUE, two.tailed = TRUE, traitorder = NA, nthreads = 1)
resAll[i,1]<-res[[i]][1]
resAll[i,2]<-res[[i]][2]
obsCorrs[[i]]<-res[[i]][6]
nullDists[[i]]<-res[[i]][8]
}

colnames(resAll)<-c("estimate","pVal")
write.table(as.matrix(resAll),file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-Trophic_ALL.txt")
write.table(obsCorrs,file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-Trophic_ALL_obsCorr.txt")
write.table(nullDists,file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-Trophic_ALL_nullDists.txt")






##
#####
# Repeat STRAPP for the SUBSET categories of Lifemode and Trophic level... same as did for PGLS.
##
## >> Crap, so this is actually a FLAWED PREMISE.
# > Can't subset edata but to monophyletic clades-- and these LM and trophic cats are by no means Mono
# So this is innappropriate.
library(BAMMtools); library(coda); library(phytools); library(ape); library(phangorn)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/_NDexp_STRAPP_analyses")




##
# Load back in the OBS and NULL dists...
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/_NDexp_STRAPP_analyses")

#obsCorrs<-read.table(file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-bodyMass_ALL_obsCorr.txt")
#nullDists<-read.table(file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-bodyMass_ALL_nullDists.txt")
#obsCorrs<-read.table(file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-geoRange_ALL_obsCorr.txt")
#nullDists<-read.table(file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-geoRange_ALL_nullDists.txt")
#obsCorrs<-read.table(file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-Lifemode_ALL_obsCorr.txt")
#nullDists<-read.table(file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-Lifemode_ALL_nullDists.txt")
obsCorrs<-read.table(file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-Trophic_ALL_obsCorr.txt")
nullDists<-read.table(file="res-STRAPP_10trees_NDexp_1000perms_speciation-vs-Trophic_ALL_nullDists.txt")

obsCorrsALL<-as.vector(unlist(obsCorrs))
nullDistsALL<-as.vector(unlist(nullDists))

obsStats<-c(mean(obsCorrsALL),quantile(obsCorrsALL,c(0.025,0.975)))
nullStats<-c(mean(nullDistsALL),quantile(nullDistsALL,c(0.025,0.975)))

# this is how STRAPP function calculates pVals-- in traitDependentBAMM()
# two-tailed
pval <- sum(abs(obsCorrsALL) <= abs(nullDistsALL))/length(nullDistsALL)
    # bodyMass-- [1] 0.0854
	# geoRange-- [1] 0.2739
	# lifeMode-- [1] 0.4031
	# trophicLevel-- [1] 0.5459

# one-tail positive
pval <- sum(obsCorrsALL <= nullDistsALL)/length(nullDistsALL)
	# bodyMass-- [1] 0.0421

# one-tail negative
pval <- sum(obsCorrsALL >= nullDistsALL)/length(nullDistsALL)
	# geoRange-- [1] 0.1353



obsCol<-rgb(1,0,0,alpha=0.5)
nullCol<-grey(0.7,alpha=0.5)

#pdf(file="plot-STRAPP_10trees_NDexp_1000perms_speciation-vs-Lifemode_ALLmams.pdf")
#pdf(file="plot-STRAPP_10trees_NDexp_1000perms_speciation-vs-bodyMass_ALLmams.pdf")
#pdf(file="plot-STRAPP_10trees_NDexp_1000perms_speciation-vs-geoRange_ALLmams.pdf")
pdf(file="plot-STRAPP_10trees_NDexp_1000perms_speciation-vs-Trophic_ALLmams.pdf")
yMax<- 0.0055 #20 #0.002 7
xMin<--100 #-0.4 #100 0.8
xMax<-2000 #0.4 #4000 0.8
#plot(density(nullDistsALL), col=nullCol, main="", xlab="",ylab="",bty="n", axes=F, xlim=c(-xDim,xDim),ylim=c(0,yMax))#7))
plot(density(nullDistsALL), col=nullCol, main="", xlab="",ylab="",bty="n", axes=F, xlim=c(xMin,xMax),ylim=c(0,yMax))#7))
polygon(density(nullDistsALL), col=nullCol, border=nullCol, bty="n")

polygon(density(obsCorrsALL), col=obsCol, border=obsCol, bty="n")

x.tick <- seq(from=xMin,to=xMax,length.out=5)
axis(at=x.tick, labels=TRUE, side=1, line=0.3, cex.axis=1, lwd=1)
#mtext(side=1,text="Spearman correlation of BAMM speciation rate ~ body size", line =3, font=2)
#mtext(side=1,text="Spearman correlation", line =3, font=2)
mtext(side=1,text="Kuskall-Wallis test statistic", line =3, font=2)
dens.rate1 <- density(obsCorrsALL)$y
dens.rate2 <- density(nullDistsALL)$y
axis(at=c(min(dens.rate1),0.45*max(dens.rate1),0.9*max(dens.rate1)), labels=c(0,0.45,0.9), side=2, las=1, lwd=1, cex.axis=1)

lines(c(nullStats[[1]],nullStats[[1]]),c(0,max(dens.rate2)),col=grey(0.7),lwd=4)
lines(c(obsStats[[1]],obsStats[[1]]),c(0,max(dens.rate1)),col="red",lwd=4)

mtext(text="BAMM speciation rate ~ trophic category", font=2,side=3, adj=0,padj=4.5)
#mtext(text="BAMM speciation rate ~ lifemode category", font=2,side=3, adj=0,padj=4.5)
#mtext(text="BAMM speciation rate ~ body mass", font=2,side=3, adj=0,padj=4.5)
#mtext(text="BAMM speciation rate ~ geo range", font=2,side=3, adj=0,padj=4.5)

tTest<-t.test(obsCorrsALL,nullDistsALL,alternative="two.sided",conf.level=0.95)
if(tTest$p.val[[1]] < 0.001){ pVal = "< 0.001"} else { pVal = round(tTest$p.val[[1]],3)}
mtext(text=paste("t = ",round(tTest$stat[[1]],1),", df = ",round(tTest$par[[1]],0),", P = ", pVal,sep=""),side=3, adj=0,padj=6.5)

twoTail <- sum(abs(obsCorrsALL) <= abs(nullDistsALL))/length(nullDistsALL)
mtext(text=paste("STRAPP two-tailed P = ",round(twoTail,3),sep=""),side=3, adj=0,padj=8.5)

legend(x=xMin,y=yMax-(yMax/4.5),fill=c(obsCol,nullCol),legend=c("All mammals","Null"))

dev.off()




# summarize

    # estimate A numeric value for continous trait data: the average observed correlation between tip rates and the trait across the posterior samples. For categorical traits, it is a list showing the median species-specific rates for each trait state.
    # p.value A numeric value. The probability that the observed correlation is less than or equal to a sample from the null distribution.
    # obs.corr A numeric vector, the observed correlation coefficents for each posterior sample. Only present when return.full is TRUE. For binary traits, centered U statistics (U - n1* n2/2; where n1 and n2 are the number of species in each state of the binary trait) is reported.
    # null A numeric vector. The null distribution of correlation coefficients (or centered U statistics for binary traits) from permutation. Only present when return.full is TRUE.








########
# BAMMtools analyses...
###
library(BAMMtools)
library(coda)
library(phytools)
library(ape)
library(phangorn)

folders<-c(1:10)

#folders<-1
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/grace_2.5_1week_1")

library(foreach);library(doSNOW)
cl = makeCluster(10, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees = length(folders)

bbone = "NDexp" # "FBD" 
dir = "NDexp_" #"" #

foreach(i=1:ntrees, .packages=c('BAMMtools', 'coda', 'phytools', 'ape', 'phangorn'), .combine=cbind, .verbose=TRUE) %dopar% {

#for (i in 1:length(folders)){

setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/grace_2.5_1week_",dir,folders[i],sep=""))
#setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/grace_2.5_1week_NDexp_",folders[i],sep=""))

tree <- ladderize(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample10_",folders[i],".tre",sep="")))
edata <- getEventData(tree, eventdata = paste("mamPhy_",bbone,"_event_data.txt",sep=""), burnin=0.33, nsamples=1000)
#tree <- ladderize(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample10_",folders[i],".tre",sep="")))
#edata <- getEventData(tree, eventdata = "mamPhy_NDexp_event_data.txt", burnin=0.33, nsamples=1000)

# convergence
mcmcout <- read.csv(paste("mamPhy_",bbone,"_mcmc_out.txt",sep=""), header=T)

# COMMON directory
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens")
setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/",bbone,sep=""))

pdf(file=paste("mamPhy_",bbone,"_mcmc_out.",folders[i],".pdf",sep=""))
plot(mcmcout$logLik ~ mcmcout$generation)
dev.off()

burnstart <- floor(0.33 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

pdf(file=paste("prior-posterior_compare",folders[i],".pdf",sep=""))
XX<-plotPrior(postburn, expectedNumberOfShifts = 1, burnin = 0.0)
dev.off()

post_probs <- table(postburn$N_shifts) / nrow(postburn)
allrates <- getCladeRates(edata)

treeNum<-folders[i]
ESS_nShifts<-effectiveSize(postburn$N_shifts)[[1]]
ESS_logLik<-effectiveSize(postburn$logLik)[[1]]
postShifts_mean<-mean(postburn$N_shifts)
postShifts_low95<-quantile(postburn$N_shifts,probs=c(0.025,0.975))[[1]]
postShifts_up95<-quantile(postburn$N_shifts,probs=c(0.025,0.975))[[2]]
likNumShifts<-as.numeric(names(post_probs[order(post_probs, decreasing=TRUE)][1]))
freqLikShifts<-post_probs[order(post_probs, decreasing=TRUE)][[1]]
allLam_mean<-mean(allrates$lambda)
allLam_low95<-quantile(allrates$lambda, c(0.025, 0.975))[[1]]
allLam_up95<-quantile(allrates$lambda, c(0.025, 0.975))[[2]]
allMu_mean<-mean(allrates$mu)
allMu_low95<-quantile(allrates$mu, c(0.025, 0.975))[[1]]
allMu_up95<-quantile(allrates$mu, c(0.025, 0.975))[[2]]

RES_SUM<-cbind.data.frame(treeNum, ESS_nShifts, ESS_logLik, postShifts_mean, postShifts_low95, postShifts_up95, likNumShifts, freqLikShifts, allLam_mean, allLam_low95, allLam_up95, allMu_mean, allMu_low95, allMu_up95)
colnames(RES_SUM)<-c("treeNum", "ESS_nShifts", "ESS_logLik", "postShifts_mean", "postShifts_low95", "postShifts_up95", "likNumShifts", "freqLikShifts", "allLam_mean", "allLam_low95", "allLam_up95", "allMu_mean", "allMu_low95", "allMu_up95")
write.table(RES_SUM, file=paste("results.treeSUMMARY.",folders[i],".txt",sep=""), append=FALSE)

shift_probs <- summary(edata)
write.table(shift_probs, file=paste("results.shiftProbs.",folders[i],".txt",sep=""), append=TRUE)

#bfmat <- computeBayesFactors(postburn, expectedNumberOfShifts=1, burnin=0.0) #burnin already done
#write.table(bfmat, file=paste("results.",folders[i],".txt",sep=""), append=TRUE)

# phylorate

pdf(file=paste("mean-phylorate_LIN_small_spec.",folders[i],".pdf",sep=""), width=8.5, height=11)
plot.bammdata(edata, spex="s", lwd=2, color.interval=c(0,1), legend=TRUE, labels=FALSE,cex=0.2)
title(main=paste("BAMM tree ",folders[i],", ",bbone," backbone",sep=""))
mtext("Mean speciation rate")
axisPhylo()
dev.off()

pdf(file=paste("mean-phylorate_LIN_small_netDiv.",folders[i],".pdf",sep=""), width=8.5, height=11)
plot.bammdata(edata, spex="netdiv", lwd=2, legend=TRUE, color.interval=c(0,1), labels=FALSE,cex=0.2)
title(main=paste("BAMM tree ",folders[i],", ",bbone," backbone",sep=""))
mtext("Mean net diversification rate")
axisPhylo()
dev.off()

pdf(file=paste("mean-phylorate_LIN_small_ext.",folders[i],".pdf",sep=""), width=8.5, height=11)
plot.bammdata(edata, spex="e", lwd=2, legend=TRUE, color.interval=c(0,1), labels=FALSE,cex=0.2)
title(main=paste("BAMM tree ",folders[i],", ",bbone," backbone",sep=""))
mtext("MSC extinction rate")
axisPhylo()
dev.off()

pdf(file=paste("mean-phylorate_LIN_large_netdiv.",folders[i],".pdf",sep=""), width=8.5, height=200)
plot.bammdata(edata, spex="netdiv", lwd=2, legend=TRUE, color.interval=c(0,1), labels=TRUE,cex=0.2)
title(main=paste("BAMM tree ",folders[i],", ",bbone," backbone",sep=""))
mtext("MSC net diversification  rate")
axisPhylo()
dev.off()

pdf(file=paste("mean-phylorate_POL_large_netdiv.",folders[i],".pdf",sep=""), width=35, height=35)
plot.bammdata(edata, spex="netdiv", lwd=2, legend=TRUE, color.interval=c(0,1), method="polar", labels=TRUE,cex=0.11)
dev.off()

#pdf(file="mean-phylorate_25th_wShifts.pdf")
#index <- 25
#e2 <- subsetEventData(edata, index = index)
#plot.bammdata(e2, lwd=2, legend=TRUE)
#addBAMMshifts(e2, cex=2)
#dev.off()

# # Credible shifts
# css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
# #write.table(summary(css), file=paste("results.",folders[i],".txt",sep=""), append=TRUE)
# 
# pdf(file=paste("credibleShiftSet",folders[i],".pdf",sep=""))
# plot.credibleshiftset(css)
# dev.off()
# 
# # RATES through time...
# pdf(file=paste("rateThroughTime_spec.",folders[i],".pdf",sep=""))
# plotRateThroughTime(edata, ratetype="speciation")
# dev.off()
# 
# pdf(file=paste("rateThroughTime_ext.",folders[i],".pdf",sep=""))
# plotRateThroughTime(edata, ratetype="extinction")
# dev.off()
# 
# pdf(file=paste("rateThroughTime_div.",folders[i],".pdf",sep=""))
# plotRateThroughTime(edata, ratetype="netdiv")
# dev.off()
# 
# cohort matrix...
cmat <- getCohortMatrix(edata)
pdf(file=paste("cohortMatrixPlot.",folders[i],".pdf",sep=""))
cohorts(cmat, edata)
title(main=paste("BAMM tree ",folders[i],", ",bbone," backbone",sep=""))
dev.off()


# maximum shift sets...

msc.set <- maximumShiftCredibility(edata, maximize='product')
msc.config <- subsetEventData(edata, index = msc.set$sampleindex)

#write.table("MSC - likeihood", file=paste("results.MSC.",folders[i],".txt",sep=""), append=FALSE)
#write.table(msc.set$bestconfigs, file=paste("results.MSC.",folders[i],".txt",sep=""), append=TRUE)
#write.table(msc.set$scores[msc.set$sampleindex], file=paste("results.MSC.",folders[i],".txt",sep=""), append=TRUE)
#write.table("MSC - shifts", file=paste("results.MSC.",folders[i],".txt",sep=""), append=TRUE)
#write.table(shiftNodes, file=paste("results.MSC.",folders[i],".txt",sep=""), append=TRUE)
write.table(msc.config$eventData, file=paste("results.MSC.nodeSHIFTS.",folders[i],".txt",sep=""), append=FALSE)


pdf(file=paste("MSC_ShiftSet_LIN_small_spec.",folders[i],".pdf",sep=""),width=8.5,height=11) 
X<-plot.bammdata(msc.config, legend=FALSE, color.interval=c(0,1), lwd=2,spex="s")
addBAMMshifts(msc.config, cex = 2)
title(main=paste("BAMM tree ",folders[i],", ",bbone," backbone",sep=""))
mtext("MSC speciation rate")
axisPhylo()
#addBAMMlegend(X)
#nodelabels(cex=0.3)
dev.off()

pdf(file=paste("MSC_ShiftSet_LIN_small_netDiv.",folders[i],".pdf",sep=""),width=8.5,height=11) 
X<-plot.bammdata(msc.config, legend=FALSE,color.interval=c(0,1), lwd=2,spex="netdiv")
addBAMMshifts(msc.config, cex = 2)
title(main=paste("BAMM tree ",folders[i],", ",bbone," backbone",sep=""))
mtext("MSC net diversification rate")
axisPhylo()
#nodelabels(cex=0.3)
dev.off()

pdf(file=paste("MSC_ShiftSet_LIN_small_ext.",folders[i],".pdf",sep=""),width=8.5,height=11) 
X<-plot.bammdata(msc.config, legend=FALSE,color.interval=c(0,1), lwd=2,spex="e")
addBAMMshifts(msc.config, cex = 2)
title(main=paste("BAMM tree ",folders[i],", ",bbone," backbone",sep=""))
mtext("MSC extinction rate")
axisPhylo()
#nodelabels(cex=0.3)
dev.off()

pdf(file=paste("MSC_ShiftSet_LIN_large_netdiv.",folders[i],".pdf",sep=""),width=8.5,height=200) 
X<-plot.bammdata(msc.config, legend=FALSE, color.interval=c(0,1), lwd=2,spex="e")
addBAMMshifts(msc.config, cex = 2)
title(main=paste("BAMM tree ",folders[i],", ",bbone," backbone",sep=""))
mtext("MSC net diversification rate")
axisPhylo()
#nodelabels(cex=0.3)
dev.off()

pdf(file=paste("MSC_ShiftSet_POL_large_netdiv.",folders[i],".pdf",sep=""), width=35, height=35)
X<-plot.bammdata(msc.config, legend=FALSE,lwd=2, color.interval=c(0,1), method="polar", labels=TRUE,cex=0.11)
title(main=paste("BAMM tree ",folders[i],", ",bbone," backbone",sep=""))
mtext("MSC net diversification rate")
addBAMMshifts(msc.config, cex = 2, method="polar")
dev.off()


 #   For more information about these color palettes visit <URL:
 #    http://colorbrewer2.org> and <URL:
 #    http://geography.uoregon.edu/datagraphics/color_scales.htm> or use
 #    the help files of the R packages ‘RColorBrewer’ and ‘dichromat’.

#color.interval=c(-0.2,1) # netdiv
#color.interval=c(0,1) # speciation
#color.interval=c(0,0.5) # netdiv

######
# BREAK out the MSC nodes to test.

MSC_shiftNodes<-getShiftNodesFromIndex(edata, index= msc.set$sampleindex)

cladeTaxNums<-vector("list",length(MSC_shiftNodes))
cladeTipNames<-vector("list",length(MSC_shiftNodes))
for (j in 1:length(MSC_shiftNodes)){
	xx<-Descendants(tree,node=MSC_shiftNodes[j], type="tips") 
	cladeTaxNums[[j]]<-length(xx[[1]])
	cladeTipNames[[j]]<-tree$tip.label[xx[[1]]]
}

nodeShifted<-data.frame(matrix(NA, nrow = length(MSC_shiftNodes), ncol = 1))
numTaxa<-data.frame(matrix(NA, nrow = length(MSC_shiftNodes), ncol = 1))
from<-data.frame(matrix(NA, nrow = length(MSC_shiftNodes), ncol = 1))
to<-data.frame(matrix(NA, nrow = length(MSC_shiftNodes), ncol = 1))
numNontaxa<-data.frame(matrix(NA, nrow = length(MSC_shiftNodes), ncol = 1))

Lam_mean<-data.frame(matrix(NA, nrow = length(MSC_shiftNodes), ncol = 1))
Lam_low95<-data.frame(matrix(NA, nrow = length(MSC_shiftNodes), ncol = 1))
Lam_up95<-data.frame(matrix(NA, nrow = length(MSC_shiftNodes), ncol = 1))
Lam_nonmean<-data.frame(matrix(NA, nrow = length(MSC_shiftNodes), ncol = 1))
Lam_nonlow95<-data.frame(matrix(NA, nrow = length(MSC_shiftNodes), ncol = 1))
Lam_nonup95<-data.frame(matrix(NA, nrow = length(MSC_shiftNodes), ncol = 1))
Mu_mean<-data.frame(matrix(NA, nrow = length(MSC_shiftNodes), ncol = 1))
Mu_low95<-data.frame(matrix(NA, nrow = length(MSC_shiftNodes), ncol = 1))
Mu_up95<-data.frame(matrix(NA, nrow = length(MSC_shiftNodes), ncol = 1))
Mu_nonmean<-data.frame(matrix(NA, nrow = length(MSC_shiftNodes), ncol = 1))
Mu_nonlow95<-data.frame(matrix(NA, nrow = length(MSC_shiftNodes), ncol = 1))
Mu_nonup95<-data.frame(matrix(NA, nrow = length(MSC_shiftNodes), ncol = 1))

for (j in 1:length(MSC_shiftNodes)){
	nodeShifted[j,]<-MSC_shiftNodes[j]
	numTaxa[j,]<-length(cladeTipNames[[j]])
	from[j,]<-cladeTipNames[[j]][1]
	to[j,]<-cladeTipNames[[j]][length(cladeTipNames[[j]])]
	numNontaxa[j,]<-5911-length(cladeTipNames[[j]])

	claderate <- getCladeRates(edata, node = MSC_shiftNodes[j], nodetype = "include") 
	NONcladerate <- getCladeRates(edata, node = MSC_shiftNodes[j], nodetype = "exclude") 
	
	Lam_mean[j,]<-mean(claderate$lambda)
	Lam_low95[j,]<-quantile(claderate$lambda, c(0.025, 0.975))[[1]]
	Lam_up95[j,]<-quantile(claderate$lambda, c(0.025, 0.975))[[2]]
	Lam_nonmean[j,]<-mean(NONcladerate$lambda)
	Lam_nonlow95[j,]<-quantile(NONcladerate$lambda, c(0.025, 0.975))[[1]]
	Lam_nonup95[j,]<-quantile(NONcladerate$lambda, c(0.025, 0.975))[[2]]
	Mu_mean[j,]<-mean(claderate$mu)
	Mu_low95[j,]<-quantile(claderate$mu, c(0.025, 0.975))[[1]]
	Mu_up95[j,]<-quantile(claderate$mu, c(0.025, 0.975))[[2]]
	Mu_nonmean[j,]<-mean(NONcladerate$mu)
	Mu_nonlow95[j,]<-quantile(NONcladerate$mu, c(0.025, 0.975))[[1]]
	Mu_nonup95[j,]<-quantile(NONcladerate$mu, c(0.025, 0.975))[[2]]

}
MSC_RES<-cbind(nodeShifted, numTaxa, from, to, numNontaxa, Lam_mean, Lam_low95, Lam_up95, Lam_nonmean, Lam_nonlow95, Lam_nonup95, Mu_mean, Mu_low95, Mu_up95, Mu_nonmean, Mu_nonlow95, Mu_nonup95)

colnames(MSC_RES)<-c("nodeShifted", "numTaxa", "from", "to", "numNontaxa", "Lam_mean", "Lam_low95", "Lam_up95", "Lam_nonmean", "Lam_nonlow95", "Lam_nonup95", "Mu_mean", "Mu_low95", "Mu_up95", "Mu_nonmean", "Mu_nonlow95", "Mu_nonup95")

write.table(MSC_RES,file=paste("results.MSC.nodeRATES.",folders[i],".txt",sep=""))



#for (j in 1:length(cladeEdata)){
#	if(length(cladeEdata[[j]]$tip.label) > 2 )
#
#	pdf(file=paste("mean-phylorate_LIN_small_netDiv.",folders[i],".node",MSC_shiftNodes[j],".pdf",sep=""), width=8.5, height=8.5)
#	plot.bammdata(cladeEdata[[j]], spex="netdiv", lwd=2, legend=TRUE, labels=FALSE,cex=0.2)
#	dev.off()
#
#	pdf(file=paste("mean-phylorate_LIN_large_netDiv.",folders[i],".node",MSC_shiftNodes[j],".pdf",sep=""), width=8.5, height=8.5*(length(cladeEdata[[j]]$tip.label)/300))
#	plot.bammdata(cladeEdata[[j]], spex="netdiv", lwd=2, legend=TRUE, labels=TRUE,cex=0.2)
#	dev.off()
#
#	css_node <- credibleShiftSet(cladeEdata[[j]], expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
#	write.table(summary(css_node), file=paste("results.MSC.perNodeCSS.",folders[i],".node",MSC_shiftNodes[j],".txt",sep=""), append=FALSE)
#
#	pdf(file=paste("credibleShiftSet.",folders[i],".node",MSC_shiftNodes[j],".pdf",sep=""))
#	plot.credibleshiftset(css_node)
#	dev.off()
#
#	# RATES through time...
#	pdf(file=paste("rateThroughTime_spec.",folders[i],".node",MSC_shiftNodes[j],".pdf",sep=""))
#	plotRateThroughTime(cladeEdata[[j]], ratetype="speciation")
#	dev.off()
#	
#	pdf(file=paste("rateThroughTime_ext.",folders[i],".node",MSC_shiftNodes[j],".pdf",sep=""))
#	plotRateThroughTime(cladeEdata[[j]], ratetype="extinction")
#	dev.off()
#	
#	pdf(file=paste("rateThroughTime_div.",folders[i],".node",MSC_shiftNodes[j],".pdf",sep=""))
#	plotRateThroughTime(cladeEdata[[j]], ratetype="netdiv")
#	dev.off()
#	
#	# cohort matrix...
#	cmat <- getCohortMatrix(cladeEdata[[j]])
#	pdf(file=paste("cohortMatrixPlot.",folders[i],".node",MSC_shiftNodes[j],".pdf",sep=""))
#	cohorts(cmat, cladeEdata[[j]])
#	dev.off()
#	
#}


# DOLPHIN compare
##
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdGenesSampPC_DR-SUMMARY-expanded_",bbone,".txt",sep=""))
colnames(cladesDR)<-c("tiplabel","gen","fam","ord","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")
head(cladesDR)

whaleSp<-cladesDR[which(cladesDR$fam=="BALAENIDAE" | cladesDR$fam=="BALAENOPTERIDAE" | cladesDR$fam=="DELPHINIDAE" | cladesDR$fam=="ESCHRICHTIIDAE" | cladesDR$fam=="INIIDAE" | cladesDR$fam=="MONODONTIDAE" | cladesDR$fam=="NEOBALAENIDAE" | cladesDR$fam=="PHOCOENIDAE" | cladesDR$fam=="PHYSETERIDAE" | cladesDR$fam=="PLATANISTIDAE" | cladesDR$fam=="ZIPHIIDAE"),"tiplabel"]
NODE_whales<-getMRCA(tree, as.vector(whaleSp)) 
#Descendants(tree,node=(NODE_whales), type="tips") # 91 species

dolphinSp<-cladesDR[which(cladesDR$fam=="DELPHINIDAE" & cladesDR$tiplabel!="Lissodelphis_borealis_DELPHINIDAE_CETARTIODACTYLA" & cladesDR$tiplabel!="Lissodelphis_peronii_DELPHINIDAE_CETARTIODACTYLA"),"tiplabel"]
	# excludes 2 sp that are NON MONO with the rest. >> Lissodelphis = right whales...
NODE_dolph<-getMRCA(tree, tip=as.vector(dolphinSp)) 
#Descendants(tree,node=(NODE_dolph), type="tips") # 36 species

edata_whale<-subtreeBAMM(edata, node= NODE_whales)
NODE_nondolph<-getMRCA(as.phylo(edata_whale), tip=as.vector(dolphinSp)) 

whalerates <- getCladeRates(edata, node= NODE_whales) #141)
dolphinrates <- getCladeRates(edata, node= NODE_dolph) #141)
nondolphinwhalerate <- getCladeRates(edata_whale, node = NODE_nondolph, nodetype = "exclude") #141, nodetype = "exclude")

whaleLam_mean<-mean(whalerates$lambda)
whaleLam_low95<-quantile(whalerates$lambda, c(0.025, 0.975))[[1]]
whaleLam_up95<-quantile(whalerates$lambda, c(0.025, 0.975))[[2]]
whaleMu_mean<-mean(whalerates$mu)
whaleMu_low95<-quantile(whalerates$mu, c(0.025, 0.975))[[1]]
whaleMu_up95<-quantile(whalerates$mu, c(0.025, 0.975))[[2]]

dolphLam_mean<-mean(dolphinrates$lambda)
dolphLam_low95<-quantile(dolphinrates$lambda, c(0.025, 0.975))[[1]]
dolphLam_up95<-quantile(dolphinrates$lambda, c(0.025, 0.975))[[2]]
dolphMu_mean<-mean(dolphinrates$mu)
dolphMu_low95<-quantile(dolphinrates$mu, c(0.025, 0.975))[[1]]
dolphMu_up95<-quantile(dolphinrates$mu, c(0.025, 0.975))[[2]]

nondolphLam_mean<-mean(nondolphinwhalerate$lambda)
nondolphLam_low95<-quantile(nondolphinwhalerate$lambda, c(0.025, 0.975))[[1]]
nondolphLam_up95<-quantile(nondolphinwhalerate$lambda, c(0.025, 0.975))[[2]]
nondolphMu_mean<-mean(nondolphinwhalerate$mu)
nondolphMu_low95<-quantile(nondolphinwhalerate$mu, c(0.025, 0.975))[[1]]
nondolphMu_up95<-quantile(nondolphinwhalerate$mu, c(0.025, 0.975))[[2]]

RES_WHA<-cbind.data.frame(NODE_whales, whaleLam_mean, whaleLam_low95, whaleLam_up95, whaleMu_mean, whaleMu_low95, whaleMu_up95, NODE_dolph, dolphLam_mean, dolphLam_low95, dolphLam_up95, dolphMu_mean, dolphMu_low95, dolphMu_up95, nondolphLam_mean, nondolphLam_low95, nondolphLam_up95, nondolphMu_mean, nondolphMu_low95, nondolphMu_up95)
names(RES_WHA)<-c("NODE_whales","whaleLam_mean", "whaleLam_low95", "whaleLam_up95", "whaleMu_mean", "whaleMu_low95", "whaleMu_up95", "NODE_dolph", "dolphLam_mean", "dolphLam_low95", "dolphLam_up95", "dolphMu_mean", "dolphMu_low95", "dolphMu_up95", "nondolphLam_mean", "nondolphLam_low95", "nondolphLam_up95", "nondolphMu_mean", "nondolphMu_low95", "nondolphMu_up95")
write.table(RES_WHA, file=paste("results.dolphins.RATES.",folders[i],".txt",sep=""), append=FALSE)

pdf(file=paste("mean-phylorate_LIN_small_netDiv_whale.",folders[i],".pdf",sep=""), width=8.5, height=8.5)
plot.bammdata(edata_whale, spex="netdiv", color.interval=c(0,0.3) ,lwd=2, legend=TRUE, labels=TRUE,cex=0.2)
dev.off()

css_whale <- credibleShiftSet(edata_whale, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
write.table(summary(css_whale), file=paste("results.dolphins.cssSHIFTS.",folders[i],".txt",sep=""), append=FALSE)

pdf(file=paste("credibleShiftSet_whale",folders[i],".pdf",sep=""))
plot.credibleshiftset(css_whale)
dev.off()


folders<-c(1:10)

ntrees = length(folders)

foreach(i=1:ntrees, .packages=c('BAMMtools', 'coda', 'phytools', 'ape', 'phangorn'), .combine=cbind, .verbose=TRUE) %dopar% {

#for (i in 1:length(folders)){

setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/grace_2.5_1week_",dir,folders[i],sep=""))
#setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/grace_2.5_1week_NDexp_",folders[i],sep=""))

tree <- ladderize(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample10_",folders[i],".tre",sep="")))
edata <- getEventData(tree, eventdata = paste("mamPhy_",bbone,"_event_data.txt",sep=""), burnin=0.33, nsamples=1000)
#tree <- ladderize(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample10_",folders[i],".tre",sep="")))
#edata <- getEventData(tree, eventdata = "mamPhy_NDexp_event_data.txt", burnin=0.33, nsamples=1000)

# convergence
mcmcout <- read.csv(paste("mamPhy_",bbone,"_mcmc_out.txt",sep=""), header=T)

# COMMON directory
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens")
setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/",bbone,sep=""))

burnstart <- floor(0.33 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdGenesSampPC_DR-SUMMARY-expanded_",bbone,".txt",sep=""))
colnames(cladesDR)<-c("tiplabel","gen","fam","ord","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")
head(cladesDR)


######
# BREAK out the ORDERS to test.

ordNames<-names(which(sort(table(cladesDR$ord),decreasing=TRUE) > 2))

ordTipNames<-vector("list",length(ordNames))
for (j in 1:length(ordTipNames)){
	ordTipNames[[j]]<-as.character(cladesDR[which(cladesDR$ord==ordNames[j]),"tiplabel"])
}

ordNODES<-vector("list",length(ordNames))
edata_ORDS<-vector("list",length(ordNames))
for (j in 1:length(ordTipNames)){
	ordNODES[[j]]<-getMRCA(tree, as.vector(ordTipNames[[j]]))
	edata_ORDS[[j]]<-subtreeBAMM(edata, node= ordNODES[[j]])
}
ord<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1))

Lam_mean<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1))
Lam_low95<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1))
Lam_up95<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1))
Lam_nonmean<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1))
Lam_nonlow95<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1))
Lam_nonup95<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1))
Mu_mean<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1))
Mu_low95<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1))
Mu_up95<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1))
Mu_nonmean<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1))
Mu_nonlow95<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1))
Mu_nonup95<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 1))

for (j in 1:length(ordNames)){
	ordrate <- getCladeRates(edata, node= ordNODES[[j]])
	NONordrate <- getCladeRates(edata, node= ordNODES[[j]], nodetype = "exclude") 
	
	ord[j,]<-ordNames[j]

	Lam_mean[j,]<-mean(ordrate$lambda)
	Lam_low95[j,]<-quantile(ordrate$lambda, c(0.025, 0.975))[[1]]
	Lam_up95[j,]<-quantile(ordrate$lambda, c(0.025, 0.975))[[2]]
	Lam_nonmean[j,]<-mean(NONordrate$lambda)
	Lam_nonlow95[j,]<-quantile(NONordrate$lambda, c(0.025, 0.975))[[1]]
	Lam_nonup95[j,]<-quantile(NONordrate$lambda, c(0.025, 0.975))[[2]]
	Mu_mean[j,]<-mean(ordrate$mu)
	Mu_low95[j,]<-quantile(ordrate$mu, c(0.025, 0.975))[[1]]
	Mu_up95[j,]<-quantile(ordrate$mu, c(0.025, 0.975))[[2]]
	Mu_nonmean[j,]<-mean(NONordrate$mu)
	Mu_nonlow95[j,]<-quantile(NONordrate$mu, c(0.025, 0.975))[[1]]
	Mu_nonup95[j,]<-quantile(NONordrate$mu, c(0.025, 0.975))[[2]]

}
ORD_RES<-cbind(ord, Lam_mean, Lam_low95, Lam_up95, Lam_nonmean, Lam_nonlow95, Lam_nonup95, Mu_mean, Mu_low95, Mu_up95, Mu_nonmean, Mu_nonlow95, Mu_nonup95)

colnames(ORD_RES)<-c("ord", "Lam_mean", "Lam_low95", "Lam_up95", "Lam_nonmean", "Lam_nonlow95", "Lam_nonup95", "Mu_mean", "Mu_low95", "Mu_up95", "Mu_nonmean", "Mu_nonlow95", "Mu_nonup95")

write.table(ORD_RES,file=paste("results.ordRATES.",folders[i],".txt",sep=""), append=FALSE)


for (j in 1:length(edata_ORDS)){
	pdf(file=paste("mean-phylorate_LIN_small_netDiv_",ordNames[j],".",folders[i],".pdf",sep=""), width=8.5, height=11)
	plot.bammdata(edata_ORDS[[j]], spex="netdiv", color.interval=c(0,1) ,lwd=2, legend=TRUE, labels=FALSE,cex=0.2)
	dev.off()

	css_ord <- credibleShiftSet(edata_ORDS[[j]], expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
	write.table(summary(css_ord), file=paste("results.",ordNames[j],".cssSHIFTS.",folders[i],".txt",sep=""), append=TRUE)

	pdf(file=paste("credibleShiftSet_",ordNames[j],".",folders[i],".pdf",sep=""))
	plot.credibleshiftset(css_ord)
	dev.off()

	# RATES through time...
	pdf(file=paste("rateThroughTime_spec.",ordNames[j],".",folders[i],".pdf",sep=""))
	plotRateThroughTime(edata_ORDS[[j]], ratetype="speciation")
	dev.off()
	
	pdf(file=paste("rateThroughTime_ext.",ordNames[j],".",folders[i],".pdf",sep=""))
	plotRateThroughTime(edata_ORDS[[j]], ratetype="extinction")
	dev.off()
	
	pdf(file=paste("rateThroughTime_div.",ordNames[j],".",folders[i],".pdf",sep=""))
	plotRateThroughTime(edata_ORDS[[j]], ratetype="netdiv")
	dev.off()
	
	# cohort matrix...
	cmat <- getCohortMatrix(edata_ORDS[[j]])
	pdf(file=paste("cohortMatrixPlot.",ordNames[j],".",folders[i],".pdf",sep=""))
	cohorts(cmat, edata_ORDS[[j]])
	dev.off()
	

# Breakout MSC **within** the ORDS:: 
##
clade_msc.set <- maximumShiftCredibility(edata_ORDS[[j]], maximize='product')
clade_msc.config <- subsetEventData(edata_ORDS[[j]], index = clade_msc.set$sampleindex)
clade_shiftNodes<-getShiftNodesFromIndex(edata_ORDS[[j]], index= clade_msc.set$sampleindex)

if (length(clade_shiftNodes) > 0){

cladeTaxNums<-vector("list",length(clade_shiftNodes))
cladeTipNames<-vector("list",length(clade_shiftNodes))
for (k in 1:length(clade_shiftNodes)){
	xx<-Descendants(x=as.phylo(edata_ORDS[[j]]),node=clade_shiftNodes[k], type="tips") 
	cladeTaxNums[[k]]<-length(xx[[1]])
	cladeTipNames[[k]]<-edata_ORDS[[j]]$tip.label[xx[[1]]]
}

nodeShifted<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
numTaxa<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
from<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
to<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
numNontaxa<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))

Lam_mean<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
Lam_low95<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
Lam_up95<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
Lam_nonmean<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
Lam_nonlow95<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
Lam_nonup95<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
Mu_mean<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
Mu_low95<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
Mu_up95<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
Mu_nonmean<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
Mu_nonlow95<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
Mu_nonup95<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
mean_Div<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
mean_nonDiv<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
mean_rateChange<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))

for (k in 1:length(clade_shiftNodes)){
	nodeShifted[k,]<-clade_shiftNodes[k]
	numTaxa[k,]<-length(cladeTipNames[[k]])
	from[k,]<-cladeTipNames[[k]][1]
	to[k,]<-cladeTipNames[[k]][length(cladeTipNames[[k]])]
	numNontaxa[k,]<-length(edata_ORDS[[j]]$tip.label)-length(cladeTipNames[[k]])

	claderate <- getCladeRates(edata_ORDS[[j]], node = clade_shiftNodes[k], nodetype = "include") 
	NONcladerate <- getCladeRates(edata_ORDS[[j]], node = clade_shiftNodes[k], nodetype = "exclude") 
	
	Lam_mean[k,]<-mean(claderate$lambda)
	Lam_low95[k,]<-quantile(claderate$lambda, c(0.025, 0.975))[[1]]
	Lam_up95[k,]<-quantile(claderate$lambda, c(0.025, 0.975))[[2]]
	Lam_nonmean[k,]<-mean(NONcladerate$lambda)
	Lam_nonlow95[k,]<-quantile(NONcladerate$lambda, c(0.025, 0.975))[[1]]
	Lam_nonup95[k,]<-quantile(NONcladerate$lambda, c(0.025, 0.975))[[2]]
	Mu_mean[k,]<-mean(claderate$mu)
	Mu_low95[k,]<-quantile(claderate$mu, c(0.025, 0.975))[[1]]
	Mu_up95[k,]<-quantile(claderate$mu, c(0.025, 0.975))[[2]]
	Mu_nonmean[k,]<-mean(NONcladerate$mu)
	Mu_nonlow95[k,]<-quantile(NONcladerate$mu, c(0.025, 0.975))[[1]]
	Mu_nonup95[k,]<-quantile(NONcladerate$mu, c(0.025, 0.975))[[2]]
	mean_Div[k,]<-(Lam_mean[k,]-Mu_mean[k,])
	mean_nonDiv[k,]<-(Lam_nonmean[k,]-Mu_nonmean[k,])
	mean_rateChange[k,]<-round(mean_Div[k,]/mean_nonDiv[k,], digits=4)

}
ORD_MSC_RES<-cbind(nodeShifted, numTaxa, from, to, numNontaxa, Lam_mean, Lam_low95, Lam_up95, Lam_nonmean, Lam_nonlow95, Lam_nonup95, Mu_mean, Mu_low95, Mu_up95, Mu_nonmean, Mu_nonlow95, Mu_nonup95, mean_Div, mean_nonDiv, mean_rateChange)

colnames(ORD_MSC_RES)<-c("nodeShifted", "numTaxa", "from", "to", "numNontaxa", "Lam_mean", "Lam_low95", "Lam_up95", "Lam_nonmean", "Lam_nonlow95", "Lam_nonup95", "Mu_mean", "Mu_low95", "Mu_up95", "Mu_nonmean", "Mu_nonlow95", "Mu_nonup95", "mean_Div", "mean_nonDiv", "mean_rateChange")

write.table(ORD_MSC_RES,file=paste("results.MSC.nodeRATES.",ordNames[j],".",folders[i],".txt",sep=""))

}

}

}



#######
## BREAK out the PATCH CLADES to test.
#
#patchNames<-names(which(sort(table(cladesDR$PC),decreasing=TRUE) > 2))
#
#patchTipNames<-vector("list",length(patchNames))
#for (j in 1:length(patchTipNames)){
#	patchTipNames[[j]]<-as.character(cladesDR[which(cladesDR$PC==patchNames[j]),"tiplabel"])
#}
#
#patchNODES<-vector("list",length(patchNames))
#edata_PCS<-vector("list",length(patchNames))
#for (j in 1:length(patchTipNames)){
#	patchNODES[[j]]<-getMRCA(tree, as.vector(patchTipNames[[j]]))
#	edata_PCS[[j]]<-subtreeBAMM(edata, node= patchNODES[[j]])
#}
#
#pc<-data.frame(matrix(NA, nrow = length(patchNames), ncol = 1))
#mean<-data.frame(matrix(NA, nrow = length(patchNames), ncol = 1))
#low95<-data.frame(matrix(NA, nrow = length(patchNames), ncol = 1))
#up95<-data.frame(matrix(NA, nrow = length(patchNames), ncol = 1))
#nonmean<-data.frame(matrix(NA, nrow = length(patchNames), ncol = 1))
#nonlow95<-data.frame(matrix(NA, nrow = length(patchNames), ncol = 1))
#nonup95<-data.frame(matrix(NA, nrow = length(patchNames), ncol = 1))
#
#for (j in 1:length(patchNames)){
#	pcrate <- getCladeRates(edata, node= patchNODES[[j]])
#	NONpcrate <- getCladeRates(edata, node= patchNODES[[j]], nodetype = "exclude") 
#	
#	pc[j,]<-patchNames[j]
#	mean[j,]<-mean(pcrate$lambda)
#	low95[j,]<-quantile(pcrate$lambda, c(0.025, 0.975))[[1]]
#	up95[j,]<-quantile(pcrate$lambda, c(0.025, 0.975))[[2]]
#	nonmean[j,]<--mean(NONpcrate$lambda)
#	nonlow95[j,]<-quantile(NONpcrate$lambda, c(0.025, 0.975))[[1]]
#	nonup95[j,]<-quantile(NONpcrate$lambda, c(0.025, 0.975))[[2]]
#}
#PC_RES<-cbind(pc, mean, low95, up95, nonmean, nonlow95, nonup95)
#colnames(PC_RES)<-c("pc", "mean", "low95", "up95", "nonmean", "nonlow95", "nonup95")
#write.table(PC_RES,file=paste("results.patchRATES.",folders[i],".txt",sep=""), append=FALSE)
#
# 
#for (j in 1:length(edata_PCS)){
#	pdf(file=paste("mean-phylorate_LIN_small_netDiv_",patchNames[j],".",folders[i],".pdf",sep=""), width=8.5, height=8.5)
#	plot.bammdata(edata_PCS[[j]], spex="netdiv", lwd=2, legend=TRUE, labels=FALSE,cex=0.2)
#	dev.off()
#
#	pdf(file=paste("mean-phylorate_LIN_large_netDiv_",patchNames[j],".",folders[i],".pdf",sep=""), width=8.5, height=8.5*(length(edata_PCS[[j]]$tip.label)/300))
#	plot.bammdata(edata_PCS[[j]], spex="netdiv", lwd=2, legend=TRUE, labels=TRUE,cex=0.2)
#	dev.off()
#
#	css_pc <- credibleShiftSet(edata_PCS[[j]], expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
#	write.table(summary(css_pc), file=paste("results.",patchNames[j],".cssSHIFTS.",folders[i],".txt",sep=""), append=TRUE)
#
#	pdf(file=paste("credibleShiftSet_",patchNames[j],".",folders[i],".pdf",sep=""))
#	plot.credibleshiftset(css_pc)
#	dev.off()
#
#	# RATES through time...
#	pdf(file=paste("rateThroughTime_spec.",patchNames[j],".",folders[i],".pdf",sep=""))
#	plotRateThroughTime(edata_PCS[[j]], ratetype="speciation")
#	dev.off()
#	
#	pdf(file=paste("rateThroughTime_ext.",patchNames[j],".",folders[i],".pdf",sep=""))
#	plotRateThroughTime(edata_PCS[[j]], ratetype="extinction")
#	dev.off()
#	
#	pdf(file=paste("rateThroughTime_div.",patchNames[j],".",folders[i],".pdf",sep=""))
#	plotRateThroughTime(edata_PCS[[j]], ratetype="netdiv")
#	dev.off()
#	
#	# cohort matrix...
#	cmat <- getCohortMatrix(edata_PCS[[j]])
#	pdf(file=paste("cohortMatrixPlot.",patchNames[j],".",folders[i],".pdf",sep=""))
#	cohorts(cmat, edata_PCS[[j]])
#	dev.off()
#
## Breakout MSC **within** the PCS:: 
###
#clade_msc.set <- maximumShiftCredibility(edata_PCS[[j]], maximize='product')
#clade_msc.config <- subsetEventData(edata_PCS[[j]], index = clade_msc.set$sampleindex)
#clade_shiftNodes<-getShiftNodesFromIndex(edata_PCS[[j]], index= clade_msc.set$sampleindex)
#
#cladeTipNames<-vector("list",length(clade_shiftNodes))
#cladeEdata<-vector("list",length(clade_shiftNodes))
#for (k in 1:length(cladeEdata)){
#	cladeEdata[[k]]<-subtreeBAMM(edata_PCS[[j]], node= clade_shiftNodes[k])
#	cladeTipNames[[k]]<-cladeEdata[[k]]$tip.label
#}
#
#nodes<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
#taxa<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
#from<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
#to<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
#
#mean<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
#low95<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
#up95<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
#nontaxa<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
#nonmean<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
#nonlow95<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
#nonup95<-data.frame(matrix(NA, nrow = length(clade_shiftNodes), ncol = 1))
#
#for (k in 1:length(clade_shiftNodes)){
#	nodes[k,]<-clade_shiftNodes[k]
#	taxa[k,]<-length(cladeTipNames[[k]])
#	from[k,]<-cladeTipNames[[k]][1]
#	to[k,]<-cladeTipNames[[k]][length(cladeTipNames[[k]])]
#
#	claderate <- getCladeRates(edata_PCS[[j]], node = clade_shiftNodes[k], nodetype = "include") 
#	NONcladerate <- getCladeRates(edata_PCS[[j]], node = clade_shiftNodes[k], nodetype = "exclude") 
#	
#	mean[k,]<-mean(claderate$lambda)
#	low95[k,]<-quantile(claderate$lambda, c(0.025, 0.975))[[1]]
#	up95[k,]<-quantile(claderate$lambda, c(0.025, 0.975))[[2]]
#	nontaxa[k,]<-length(edata_PCS[[j]]$tip.label)-length(cladeTipNames[[k]])
#	nonmean[k,]<--mean(NONcladerate$lambda)
#	nonlow95[k,]<-quantile(NONcladerate$lambda, c(0.025, 0.975))[[1]]
#	nonup95[k,]<-quantile(NONcladerate$lambda, c(0.025, 0.975))[[2]]
#}
#PC_MSC_RES<-cbind(nodes, taxa, from, to, mean, low95, up95, nontaxa, nonmean, nonlow95, nonup95)
#colnames(PC_MSC_RES)<-c("nodes", "taxa", "from", "to", "mean", "low95", "up95", "nontaxa", "nonmean", "nonlow95", "nonup95")
#
#write.table(PC_MSC_RES,file=paste("results.MSC.nodeRATES.",patchNames[j],".",folders[i],".txt",sep=""))
#
#}
#
#
#}
#
#
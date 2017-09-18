##
# COMET in TESS
#####
library(ape); library(phytools); library(picante); library(plotrix); library(phangorn); library(phyloch); library(TESS)
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")

#bbone<- "NDexp" # "FBD"
#mamMCC_orig<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")
#
#is.ultrametric(mamMCC_orig) ## fails
#
#mamMCC<-nnls.tree(cophenetic(mamMCC_orig),mamMCC_orig,rooted=TRUE)
#is.ultrametric(mamMCC) ## should pass

# DO the 10-tree loop...

library(foreach);library(doSNOW)
cl = makeCluster(10, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees=10

foreach(i=1:ntrees, .packages=c('TESS','phangorn','geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

# which backbone?
bbone<- "NDexp" #"FBD" # 

# read in 1 of 100 full trees
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")
#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_nexus-and-newickTrees/")
#mamPhy_orig<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_nexus.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 
mamPhy_orig<-ladderize(drop.tip(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample10_",i,".tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 

is.ultrametric(mamPhy_orig) ## fails

mamPhy<-nnls.tree(cophenetic(mamPhy_orig),mamPhy_orig,rooted=TRUE)
is.ultrametric(mamPhy) ## should pass

#tips<-mamMCC$tip.label
#cor(as.vector(cophenetic(mamMCC_orig)[tips,tips]), as.vector(cophenetic(mamMCC)[tips,tips])) # gives value of 1 ! Okay awesome, it was just a rounding error GIVEN the linux machine vs. my mac machine in their numeric treatments.
times <- as.numeric( branching.times(mamPhy) )

## Run CoMet
# CPP on Mass-Extinction Times (CoMET) model (May et al., 2015)
 
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/CoMET_analyses_MamPhy_7shifts-initME66-MajorMinor")
#####
# set parameters
samplingFraction <- 1.0 # all tips are eincluded in the full MamPhy

numExpectedMassExtinctions <- 1 # K-Pg boundary
#numExpectedRateChanges <- 1 # at the K-Pg, otherwise assuming not
#numExpectedMassExtinctions <- 0 # not including reference here-- make analogous to TREEPAR analyses...
numExpectedRateChanges <- 7 # seven possible shift points.
 
#expectedSurvivalProbability <- 0.05 # MAJOR-- gives a beta prior with 95% of density between 0.01 and 0.10 
expectedSurvivalProbability <- 0.30 # MINOR-- gives a beta prior with 95% of density between 0.23 and 0.38
pMassExtinctionPriorShape2 <- 100
pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 * expectedSurvivalProbability / (expectedSurvivalProbability - 1)


# Plot the density function of our beta distribution.
#curve(dbeta(x,shape1=pMassExtinctionPriorShape1, shape2=pMassExtinctionPriorShape2),n=1001, xlab= 'survival probability' ,ylab= 'density' ,las=1)
# Plot the 95% prior interval on the survival probability.
#abline(v = qbeta(c(0.025,0.975),shape1=pMassExtinctionPriorShape1, shape2=pMassExtinctionPriorShape2),lty=2)
 
### Do CoMet analysis with these parameters

## CoMet
## Use Empirical Hyperpriors, and auto-stopping...

#param<-"noShift-initME66-Major0p05"
param<-"noShift-initME66-Minor0p30"
  #param<-"7Shift-initME66-Major0p05"
  #param<-"7Shift-initME66-Minor0p30"

#posterior_directories <- paste("MamPhy_CoMET_estimHyper_10trees_tree",i,"_minESS500_run",1:4,sep="")
posterior_directories <- paste("MamPhy_CoMET_estimHyper_",param,"_10trees_tree",i,"_minESS500_run",1:4,sep="")

for(j in 1:length(posterior_directories)){

tess.analysis(mamPhy,
  empiricalHyperPriors = TRUE,
  #empiricalHyperPriors = FALSE, empiricalHyperPriorForm = "normal", speciationRatePriorMean = 0.2857622, speciationRatePriorStDev = 0.02137434, extinctionRatePriorMean = 0.1897766, extinctionRatePriorStDev = 0.02876494, initialSpeciationRate = 0.2857622, initialExtinctionRate = 0.1897766,
  samplingProbability = 1,
  #numExpectedRateChanges = 7, # seven possible shift points.
  numExpectedRateChanges = 0, # constant rates
  #estimateNumberRateChanges = TRUE,
  estimateNumberRateChanges = FALSE,
  numExpectedMassExtinctions = 1, # K-Pg boundary
  estimateNumberMassExtinctions = TRUE,
  estimateMassExtinctionTimes = TRUE,
  tInitialMassExtinction = (max(times)-66),
  #pInitialMassExtinction = 0.05, # MAJOR event
  pInitialMassExtinction = 0.30, # MINOR event
  pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
  pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
  MAX_ITERATIONS = 100000000,
  MAX_TIME = 24*60*60,
  MIN_ESS = 500,
  dir = posterior_directories[j])

}

} # close multiple trees loop



######
# Repeat on CLADES...
##########
# =================================
# GLOBAL TREE - 10 trees, ND
# CoMET model RATE SHIFTS
###
# ready the trees and subsets first.
library(ape); library(picante); library(phytools); library(geiger); library(DDD); library(pspline);library(TreePar); library(TESS); library(phangorn)
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")

#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_10trees")

library(foreach);library(doSNOW)
cl = makeCluster(10, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees=10

foreach(i=1:ntrees, .packages=c('TESS','phangorn','geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

# set tree dir
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")

# which backbone?
bbone<- "NDexp" #"FBD" # 

# load tree
#mamPhy<-ladderize(drop.tip(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample10_",i,".tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 
mamPhy_orig<-ladderize(drop.tip(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample10_",i,".tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 

is.ultrametric(mamPhy_orig) ## fails

mamPhy<-nnls.tree(cophenetic(mamPhy_orig),mamPhy_orig,rooted=TRUE)
is.ultrametric(mamPhy) ## should pass

# tip categories
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
head(cladesDR)

ordNames<-names(which(table(cladesDR$ord) > 25))[c(1:4,8:11)] # [8:11] # excluding "DASYUROMORPHIA"  "DIDELPHIMORPHIA" "DIPROTODONTIA" bc crowns are ~ 20, 21, 30.1 Ma < 30 Ma needed for 7 splits of 5 Ma...
ordNames<-ordNames[-c(2:4,6)] # removing.. c("CARNIVORA","CETARTIODACTYLA","CHIROPTERA","LAGOMORPHA")] << all w crown times < 66 Ma
higherNames<-names(which(table(cladesDR$higher) > 1000)) #25))
setNames<-c(higherNames, ordNames)

# BREAKOUTS
mamPhy_cladeSets<-vector("list",length(setNames))
for (j in 1:length(setNames)){
    if(j < 3){ setCat<-"higher"} else {setCat<-"ord"}
    setTipNames<-cladesDR[which(as.vector(cladesDR[,setCat])==setNames[j]),"tiplabel"]
    toDrop<-setdiff(mamPhy$tip.label,as.vector(setTipNames))
    mamPhy_cladeSets[[j]]<-drop.tip(mamPhy,toDrop)
    }

## Run CoMet
# CPP on Mass-Extinction Times (CoMET) model (May et al., 2015)
 
# set results dir
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/CoMET_analyses_MamPhy_7shifts-initME66-MajorMinor_CLADES")
#####
# set parameters
samplingFraction <- 1.0 # all tips are eincluded in the full MamPhy

numExpectedMassExtinctions <- 1 # K-Pg boundary
##numExpectedRateChanges <- 1 # at the K-Pg, otherwise assuming not
##numExpectedMassExtinctions <- 0 # not including reference here-- make analogous to TREEPAR analyses...
#numExpectedRateChanges <- 7 # seven possible shift points.
# 
#expectedSurvivalProbability <- 0.5 # gives a beta prior with 95% of density between 0.45 and 0.55-- 0.01 and 0.10 
#pMassExtinctionPriorShape2 <- 100
#pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 * expectedSurvivalProbability / (expectedSurvivalProbability - 1)


####
param_all<-c("7Shift-initME66-Major0p05","7Shift-initME66-Minor0p30", "noShift-initME66-Major0p05", "noShift-initME66-Minor0p30")

numExpectedRateChanges_all <- c(7,7,0,0)
estimateNumberRateChanges_all<-c(TRUE,TRUE,FALSE,FALSE)
expectedSurvivalProbability_all <- c(0.05,0.30,0.05,0.30)

for(z in c(2,4)){

numExpectedRateChanges<-numExpectedRateChanges_all[z]
estimateNumberRateChanges<-estimateNumberRateChanges_all[z]

expectedSurvivalProbability<-expectedSurvivalProbability_all[z]
pMassExtinctionPriorShape2 <- 100
pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 * expectedSurvivalProbability / (expectedSurvivalProbability - 1)

param<-param_all[z]


for(q in 1:length(setNames)){
#foreach(q=1:length(setNames), .packages=c('ape', 'picante', 'phytools','geiger','DDD','pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar% { 

tree_file<-mamPhy_cladeSets[[q]]
posterior_directories <- paste("MamPhy_CoMET_estimHyper_",param,"_",setNames[q],"_10trees_tree",i,"_minESS500_run",1:4,sep="")

times <- as.numeric( branching.times(tree_file) )

for(j in 1:length(posterior_directories)){

tess.analysis(tree_file,
  empiricalHyperPriors = TRUE,
  #empiricalHyperPriors = FALSE,
  #empiricalHyperPriorForm = "normal",
  #speciationRatePriorMean = 0.2857622,
  #speciationRatePriorStDev = 0.02137434,
  #extinctionRatePriorMean = 0.1897766,
  #extinctionRatePriorStDev = 0.02876494,
  #initialSpeciationRate = 0.2857622,
  #initialExtinctionRate = 0.1897766,
  samplingProbability = 1,
  numExpectedRateChanges = numExpectedRateChanges,
  estimateNumberRateChanges = estimateNumberRateChanges,
  numExpectedMassExtinctions = 1, # K-Pg boundary
  estimateNumberMassExtinctions = TRUE,
  estimateMassExtinctionTimes = TRUE,
  tInitialMassExtinction = (max(times)-66),
  pInitialMassExtinction = expectedSurvivalProbability,
  pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
  pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
  MAX_ITERATIONS = 100000000,
  MAX_TIME = 24*60*60,
  MIN_ESS = 500,
  dir = posterior_directories[j])

}

} # close loop on 6 clades

} # close loop on multiple models (4)

} # close multiple trees loop (parallelized)




########
# REPEAT on other MAM PHYS...
#setwd("/mnt/data/personal/nateu/Nate_Backup/DR_calcOnOther_mamPhys/")

# subset out 10 trees...
setwd("/Users/Nate/Desktop/JETZ-pdfs/Faurby_Svenning2015_SUPP/AppendixD_producedPhylogenies_Rcode")
FS2015_phy1k<-scan("FS2015_Fully_resolved_phylogeny.nex", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") # is newick format

FS2015_10trees<-sample(FS2015_phy1k,size=10)
for(i in 1:length(FS2015_10trees)){
  phy<-read.tree(text=FS2015_10trees[[i]])
  write.tree(phy,file=paste("FS2015_Fully_resolved_phylogeny_10trees_tree",i,".tre",sep=""),append=FALSE)
  #write.tree(phy,file="FS2015_Fully_resolved_phylogeny_10trees.trees",append=TRUE)
}

setwd("/Users/Nate/Desktop/JETZ-pdfs/KUHN_et_al_2011_Polytomy_Resolver")
Kuhn_phy1k<-scan("Fritz.Resolved.Normal.cr1r5_v2b.nwk", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

Kuhn_10trees<-sample(Kuhn_phy1k,size=10)
for(i in 1:length(Kuhn_10trees)){
  phy<-read.tree(text=Kuhn_10trees[[i]])
  #write.tree(phy,file=paste("Fritz.Resolved.Normal.cr1r5_v2b_10trees_tree",i,".tre",sep=""),append=FALSE)
  write.tree(phy,file="Fritz.Resolved.Normal.cr1r5_v2b_10trees.trees",append=TRUE)
}

# NOW, run the same on both OTHER MAMPHYS... 

library(ape); library(phytools); library(picante); library(plotrix); library(phangorn); library(phyloch)
library(TESS)
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")

# KUHN...
# DO the 10-tree loop...
library(foreach);library(doSNOW)
cl = makeCluster(10, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees=10

foreach(i=1:ntrees, .packages=c('TESS','phangorn','geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

# which backbone?
bbone<- "NDexp" #"FBD" # 

# read in 1 of 10 full trees
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")
mamPhy_orig<-ladderize(read.tree(paste("Fritz.Resolved.Normal.cr1r5_v2b_10trees_tree",i,".tre",sep="")))
previousTree<-"KuhnEtAl2011"
#mamPhy_orig<-ladderize(read.tree(paste("FS2015_Fully_resolved_phylogeny_10trees_tree",i,".tre",sep="")))
#previousTree<-"FS2015"

if(is.ultrametric(mamPhy_orig)==FALSE) {
  mamPhy<-nnls.tree(cophenetic(mamPhy_orig),mamPhy_orig,rooted=TRUE)
  is.ultrametric(mamPhy) ## should pass
} else {
  mamPhy<-mamPhy_orig
}

#tips<-mamMCC$tip.label
#cor(as.vector(cophenetic(mamMCC_orig)[tips,tips]), as.vector(cophenetic(mamMCC)[tips,tips])) # gives value of 1 ! Okay awesome, it was just a rounding error GIVEN the linux machine vs. my mac machine in their numeric treatments.
times <- as.numeric( branching.times(mamPhy) )

## Run CoMet
# CPP on Mass-Extinction Times (CoMET) model (May et al., 2015)
 
#####
# set parameters
 samplingFraction <- 1.0 # all tips are eincluded in the full MamPhy

numExpectedMassExtinctions <- 1 # K-Pg boundary
numExpectedRateChanges <- 1 # at the K-Pg, otherwise assuming not
 
expectedSurvivalProbability <- 0.05 # gives a beta prior with 95% of density between 0.01 and 0.10 
pMassExtinctionPriorShape2 <- 100
pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 * expectedSurvivalProbability / (expectedSurvivalProbability - 1)
 
# Plot the density function of our beta distribution.
#curve(dbeta(x,shape1=pMassExtinctionPriorShape1, shape2=pMassExtinctionPriorShape2),n=1001, xlab= 'survival probability' ,ylab= 'density' ,las=1)
# Plot the 95% prior interval on the survival probability.
#abline(v = qbeta(c(0.025,0.975),shape1=pMassExtinctionPriorShape1, shape2=pMassExtinctionPriorShape2),lty=2)
 
### Do CoMet analysis with these parameters

## CoMet
## Use Empirical Hyperpriors, and auto-stopping...

posterior_directories <- paste(previousTree,"_CoMET_estimHyper_10trees_tree",i,"_minESS500_run",1:4,sep="")

for(j in 1:length(posterior_directories)){

tess.analysis(mamPhy,
  empiricalHyperPriors = TRUE,
  samplingProbability = samplingFraction,
  numExpectedRateChanges = numExpectedRateChanges,
  numExpectedMassExtinctions = numExpectedMassExtinctions,
  pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
  pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
  MAX_ITERATIONS = 100000000,
  MAX_TIME = 24*60*60,
  MIN_ESS = 500,
  dir = posterior_directories[j])

}

} # close multiple trees loop


####################
# SUMMARIZING
################################
# Summarize each run SEPARATE
##
library(ape); library(phytools); library(picante); library(plotrix); library(phangorn); library(phyloch)
library(TESS); library(matrixStats)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/")
source("tess.plot.output2.R")

#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_MamPhy")
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_MamPhy_subclades")
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_MamPhy_7shift-noME")
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_MamPhy_7shifts-initME66-MajorMinor")
previousTree<-"MamPhy"
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_KuhnEtAl2011")
#previousTree<-"KuhnEtAl2011"
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_FS2015")
#previousTree<-"FS2015"

# set parameters
samplingFraction <- 1.0 # all tips are eincluded in the full MamPhy

numExpectedMassExtinctions <- 1 # K-Pg boundary
#numExpectedRateChanges <- 1 # at the K-Pg, otherwise assuming not
#numExpectedMassExtinctions <- 0 # not including reference here-- make analogous to TREEPAR analyses...
numExpectedRateChanges <- 7 # seven possible shift points.
 
# expectedSurvivalProbability <- 0.05 # MAJOR-- gives a beta prior with 95% of density between 0.01 and 0.10 
# #expectedSurvivalProbability <- 0.30 # MINOR-- gives a beta prior with 95% of density between 0.23 and 0.38
# pMassExtinctionPriorShape2 <- 100
# pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 * expectedSurvivalProbability / (expectedSurvivalProbability - 1)

param<-"7shift-noME"
#param<-c("7Shift-initME66-Major0p05","7Shift-initME66-Minor0p30", "noShift-initME66-Major0p05", "noShift-initME66-Minor0p30")

numExpectedRateChanges_all <- c(7,7,0,0)
#expectedSurvivalProbability_all <- c(0.05,0.30,0.05,0.30)

binSize<-1 #5

for(z in 2:length(param)){

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_MamPhy_subclades")
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_MamPhy_7shifts-initME66-MajorMinor")
#clades<-c("_","_Euarchontoglires_","_Laurasiatheria_","_AFROSORICIDA_", "_EULIPOTYPHLA_","_PRIMATES_","_RODENTIA_")

clades<-c("_","_Euarchontoglires_","_Laurasiatheria_","_AFROSORICIDA_", "_EULIPOTYPHLA_","_PRIMATES_","_RODENTIA_","_CARNIVORA_","_CETARTIODACTYLA_","_CHIROPTERA_","_LAGOMORPHA_")

#numExpectedRateChanges<-numExpectedRateChanges_all[z]
#numExpectedMassExtinctions<-1
z=1

for(q in 2:length(clades)){
#q<-1
spShifts_per10trees<-vector("list",length=10)
exShifts_per10trees<-vector("list",length=10)
for (i in c(1:10)) {  #start multiple trees loop
#for (i in c(2,4,5,8)) {  #start multiple trees loop
#for (i in c(1:3,6:10)) {  #start multiple trees loop
posterior_directories <- paste(previousTree,"_CoMET_estimHyper_",param[z],clades[q],"10trees_tree",i,"_minESS500_run",1:4,sep="")

folder<-"subclades"

setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_",previousTree,"_",folder,"/",posterior_directories[1],sep=""))
tre<-read.nexus("input_tree.tre")
root<-max(branching.times(tre))
setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_",previousTree,"_",folder,sep=""))

for(j in 1:length(posterior_directories)){

assign(paste("output_",j,sep=""), tess.process.output(posterior_directories[j],
  numExpectedRateChanges = numExpectedRateChanges,
  numExpectedMassExtinctions = numExpectedMassExtinctions, 
  burnin=0.25,
  numIntervals=round((root/binSize),0),
  criticalBayesFactors=c(2,6,10)))
 
}

output_list <- list(output_1,output_2,output_3,output_4)

for(j in 1:length(output_list)){

pdf(file=paste(previousTree,"_CoMET_estimHyper_",param[z],clades[q],"10trees_tree",i,"_minESS500_run",j,"_results_",binSize,"Ma.pdf",sep=""),width=8,height=8)

layout.mat <- matrix(1:8,nrow=4,ncol=2,byrow=TRUE)
layout(layout.mat)
  tess.plot.output(output_list[[j]],
  fig.types = c("speciation rates","speciation shift times","extinction rates","extinction shift times","net-diversification rates","relative-extinction rates","mass extinction Bayes factors","mass extinction times"),
  las=2)
 
dev.off()

}

pdf(file=paste(previousTree,"_CoMET_estimHyper_",param[z],clades[q],"10trees_tree",i,"_minESS500_all4runsCombined_chainAssessment.pdf",sep=""), width=8,height=16)

layout.mat <- matrix(1:7,nrow=7,ncol=1,byrow=TRUE)
layout(layout.mat)
tess.plot.multichain.diagnostics(output_list, parameters = c("speciation rates"
    ,"speciation shift times","extinction rates", "extinction shift times", "net-diversification rates","relative-extinction rates","mass extinction times"), 
    las=2)
 
dev.off()



# Combine runs TOGETHER
##
# 4 runs of each trees...
###

#clades<-c("_","_Euarchontoglires_","_Laurasiatheria_","_AFROSORICIDA_", "_EULIPOTYPHLA_","_PRIMATES_","_RODENTIA_")
#for(q in 2:length(clades)){

#q=1
#for (i in c(1:10)) {  #start multiple trees loop
#for (i in c(2,4,5,8)) {  #start multiple trees loop

#for (i in c(8:10)) {  #start multiple trees loop
#  posterior_directories <- paste(previousTree,"_CoMET_estimHyper_",param,clades[q],"10trees_tree",i,"_minESS500_run",1:4,sep="")
#  
#  setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_",previousTree,"_all_ME66-7shifts/",posterior_directories[1],sep=""))
#  tre<-read.nexus("input_tree.tre")
#  root<-max(branching.times(tre))
#  setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_",previousTree,"_all_ME66-7shifts",sep=""))
#  
#  for(j in 1:length(posterior_directories)){
#  
#  assign(paste("output_",j,sep=""), tess.process.output(posterior_directories[j],
#    numExpectedRateChanges = numExpectedRateChanges,
#    numExpectedMassExtinctions = numExpectedMassExtinctions, 
#    burnin=0.25,
#    numIntervals=round((root/binSize),0),
#    criticalBayesFactors=c(2,6,10)))
#   
#  }

if(i==8){
  output_list <- list(output_1,output_2,output_3,output_4)
#  output_list[[1]][["extinction shift times"]]<-output_list[[1]][["extinction shift times"]][,-43] #[,-101]
#  output_list[[1]][["extinction Bayes factors"]]<-output_list[[1]][["extinction Bayes factors"]][1:42]
} else {
  output_list <- list(output_1,output_2,output_3,output_4)
}

keysAll <- unique(unlist(lapply(output_list, names)))

join4runs<-vector("list",(length(keysAll)))
names(join4runs)<-keysAll
for(j in c(1:18)){
  key<-keysAll[j]
  if(j %in% c(1:4)){
    join4runs[[j]]<-c(output_list[[1]][[key]],output_list[[2]][[key]],output_list[[3]][[key]],output_list[[4]][[key]])
  } else {
  join4runs[[j]]<-rbind(output_list[[1]][[key]],output_list[[2]][[key]],output_list[[3]][[key]],output_list[[4]][[key]])
  }
}
join4runs[[19]]<-output_list[[1]]$tree
join4runs[[20]]<-output_list[[1]][[20]]

for(j in c(7,8,11,12,16,17,18)){
join4runs[[j]]<-colMeans(join4runs[[j]])
}

assign(paste("join4runs_tree",i,sep=""),join4runs)

#pdf(file=paste(previousTree,"_CoMET_estimHyper_10trees_minESS500_results_join4runs_ALL10trees.pdf",sep=""))
pdf(file=paste(previousTree,"_CoMET_estimHyper_",param[z],clades[q],"10trees_tree",i,"_minESS500_results_join4runs_",binSize,"Ma.pdf",sep=""),width=8,height=8)

layout.mat <- matrix(1:8,nrow=4,ncol=2,byrow=TRUE)
layout(layout.mat)
  tess.plot.output(join4runs, #join4runs_ALL, #join4runs
  fig.types = c("speciation rates","speciation shift times","extinction rates","extinction shift times","net-diversification rates","relative-extinction rates","mass extinction Bayes factors","mass extinction times"),
  las=2)
 
dev.off()

## GET the RATE shifts > 2 BF per tree, as averaged across 4 runs each
# Speciation
type<-"speciation shift times"
criticalPP <- join4runs[[grep(strsplit(type, " ")[[1]][1], grep("CriticalPosteriorProbabilities", names(join4runs), value = TRUE), value = TRUE)]]

thisOutput <- join4runs[[type]]
meanThisOutput <- colMeans(thisOutput)
names(meanThisOutput)<-c((length(meanThisOutput)-1):0)
#meanThisOutput[(length(meanThisOutput)-1):length(meanThisOutput)]<-NA # excludes last 2 Ma from the plot
BF_greaterThan2<-as.numeric(names(meanThisOutput[which(meanThisOutput >= criticalPP[1])])) # gets shifts >= 2 BF
if(length(BF_greaterThan2)>0){ shiftSubstantial_sp<- BF_greaterThan2 } else { shiftSubstantial_sp<- NA}


# Extinction
type<-"extinction shift times"

thisOutput <- output[[type]]
meanThisOutput <- colMeans(thisOutput)
names(meanThisOutput)<-c((length(meanThisOutput)-1):0)
#meanThisOutput[(length(meanThisOutput)-1):length(meanThisOutput)]<-NA # excludes last 2 Ma from the plot
BF_greaterThan2<-as.numeric(names(meanThisOutput[which(meanThisOutput >= criticalPP[1])])) # gets shifts >= 2 BF
if(length(BF_greaterThan2)>0){ shiftSubstantial_ex<- BF_greaterThan2 } else { shiftSubstantial_ex<- NA}

spShifts_per10trees[[i]]<-shiftSubstantial_sp

exShifts_per10trees[[i]]<-shiftSubstantial_ex

} # close multiple trees loop

assign(paste("spShifts_",clades[q],sep=""),unlist(spShifts_per10trees))
assign(paste("exShifts_",clades[q],sep=""),unlist(exShifts_per10trees))

write.table(get(paste("exShifts_",clades[q],sep="")), file=paste("exShifts",clades[q],"CoMET_10trees_",param,".txt",sep=""))
write.table(get(paste("spShifts_",clades[q],sep="")), file=paste("spShifts",clades[q],"CoMET_10trees_",param,".txt",sep=""))

} #close multiple clades loop... 

# >> write each of THOSE to a list too...


###
# Combine 10 trees together

#allFileNames<-paste("join4runs_tree",1:10,sep="")
#RUNS<-c(2,4,5,8)
#RUNS<-c(1:10)
RUNS<-c(1:3,6:10)
allFileNames<-paste("join4runs_tree",RUNS,sep="")

output_list_ALL <- list(join4runs_tree1, join4runs_tree2, join4runs_tree3, 
                    #join4runs_tree4, join4runs_tree5, 
                    join4runs_tree6, join4runs_tree7, join4runs_tree8, join4runs_tree9, join4runs_tree10)
#output_list_ALL <- list(join4runs_tree2, join4runs_tree4, join4runs_tree5, join4runs_tree8)

treeAge_ALL<-c()
for(j in 1:length(output_list_ALL)){
  treeAge_ALL[j] <- max(branching.times(output_list_ALL[[j]]$tree))
}
treeAgeMed<-median(treeAge_ALL)
treeAgeMax<-max(treeAge_ALL)

#allCols_ALL<-vector("list",length(output_list_ALL)); allRows_ALL<-vector("list",length(output_list_ALL))
#for(j in 1:length(output_list_ALL)){
#    allCols<-c(); allRows<-c()
#    for(k in c(5:6,9:10,13:15)){ # the variables based on bins where NEED to even them out!
#    allCols[k-4]<-dim(output_list_ALL[[j]][[k]])[2]
#    allRows[k-4]<-dim(output_list_ALL[[j]][[k]])[1]
#    }
#allCols_ALL[[j]]<-allCols
#allRows_ALL[[j]]<-allRows
#}

allCols<-c(); 
for(j in 1:length(RUNS)){
    allCols[j]<-dim(output_list_ALL[[j]][[5]])[2]
}
allCols # this is the number of 5 Ma bins per tree
maxBin<-max(na.omit(allCols))

intToUse<-which(allCols==maxBin)[1]
treeAge<-treeAge_ALL[intToUse] # using tree 4, 41 intervals
output_list_ALL[[intToUse]][["tree"]]$root.age<-treeAge

#output_list_ALL <- list(join4runs_tree1, join4runs_tree2, join4runs_tree3, join4runs_tree4, join4runs_tree5, join4runs_tree6, join4runs_tree7, join4runs_tree8, join4runs_tree9, join4runs_tree10)

for(j in 1:length(RUNS)){
  colToAdd<-maxBin-allCols[j] #dim(output_list_ALL[[j]][[k]])[2]

  for(k in c(5:7,9:11,13:16)){ # the variables based on bins where NEED to even them out!
    if(k %in% c(7,11,16)){
      output_list_ALL[[j]][[k]]<-c(rep(output_list_ALL[[j]][[k]][1],colToAdd),output_list_ALL[[j]][[k]])
    } else {

      if(colToAdd==0){ addCol<-NULL } else if(colToAdd==1) { 
      rowToAdd<-dim(output_list_ALL[[j]][[k]])[1]
      addCol<-rep(output_list_ALL[[j]][[k]][1,1],rowToAdd)
  #    names(addCol)<-paste("col",dim(output_list_ALL[[j]][[k]])[2]+1,sep="")
    } else {
      rowToAdd<-dim(output_list_ALL[[j]][[k]])[1]
      addCol<-rep(output_list_ALL[[j]][[k]][1,1],rowToAdd)
      for(r in 1:(colToAdd-1)){
        addCol<-cbind(addCol,rep(output_list_ALL[[j]][[k]][1,1],rowToAdd))
      }
      colnames(addCol)<-paste("col",(allCols[j]+(1:colToAdd)),sep="")
    }
      output_list_ALL[[j]][[k]]<-cbind(addCol,output_list_ALL[[j]][[k]])
    }
  }
}

allCols2<-c(); 
for(j in 1:length(RUNS)){
#    allCols2[j]<-dim(output_list_ALL[[j]][[5]])[2]
    allCols2[j]<-length(output_list_ALL[[j]][[7]])
}
allCols2 # this is the number of 5 Ma bins per tree
maxBin2<-max(na.omit(allCols2))



Q<-18
join4runs_ALL<-vector("list",(length(keysAll)))
names(join4runs_ALL)<-keysAll
for(j in c(1:Q)){
  key<-keysAll[j]
  if(j %in% c(1:4)){
    join4runs_ALL[[j]]<-c(output_list_ALL[[1]][[key]], output_list_ALL[[2]][[key]], output_list_ALL[[3]][[key]], 
        output_list_ALL[[4]][[key]], output_list_ALL[[5]][[key]], 
        output_list_ALL[[6]][[key]], output_list_ALL[[7]][[key]], output_list_ALL[[8]][[key]])
        #, output_list_ALL[[9]][[key]], output_list_ALL[[10]][[key]])
    #(join4runs_tree1[[key]], join4runs_tree2[[key]], join4runs_tree3[[key]], join4runs_tree4[[key]], join4runs_tree5[[key]], join4runs_tree6[[key]], join4runs_tree7[[key]], join4runs_tree8[[key]], join4runs_tree9[[key]], join4runs_tree10[[key]])
  } else {
  join4runs_ALL[[j]]<-rbind(output_list_ALL[[1]][[key]], output_list_ALL[[2]][[key]], output_list_ALL[[3]][[key]], 
        output_list_ALL[[4]][[key]], output_list_ALL[[5]][[key]], 
        output_list_ALL[[6]][[key]], output_list_ALL[[7]][[key]], output_list_ALL[[8]][[key]])
        #, output_list_ALL[[9]][[key]], output_list_ALL[[10]][[key]])
}
}

join4runs_ALL[[19]]<-output_list_ALL[[intToUse]][["tree"]]
join4runs_ALL[[20]]<-output_list_ALL[[intToUse]][[20]] #using the deepest rooted tree to set the intervals

for(j in c(7,8,11,12,16,17,18)){
  m<-join4runs_ALL[[j]]
  #m[!is.finite(m)] <- 0
  join4runs_ALL[[j]]<-colMeans(m)
}

previousTree<-"MamPhy"
#previousTree<-"KuhnEtAl2011"
#join4runs_ALL<-KUHN_10tree_CoMET_5Ma
#previousTree<-"FS2015"
#join4runs_ALL<-FS_10tree_CoMET_5Ma

pdf(file=paste(previousTree,"_CoMET_estimHyper_",param[z],clades[q],"10trees_minESS500_results_join4runs_ALL10trees_",binSize,"Ma_treeAgeUsing.pdf",sep=""),width=8,height=8)
#pdf(file=paste(previousTree,"_CoMET_estimHyper_10trees_tree",i,"_minESS500_results_join4runs.pdf",sep=""))

op <- par(oma = c(1,1,5,1) + 0.1,   #c(bottom, left, top, right)
          mar = rep(4,4) + 0.1)

layout.mat <- matrix(1:8,nrow=4,ncol=2,byrow=TRUE)
layout(layout.mat)
  tess.plot.output2(join4runs_ALL, #join4runs
  fig.types = c("speciation rates","speciation shift times","extinction rates","extinction shift times","net-diversification rates","relative-extinction rates","mass extinction Bayes factors","mass extinction times"),
  las=2)

title(main=paste(previousTree,clades[q],param[z]," (10 trees, 4 runs each)",sep=""), xlab = "",
      ylab = "",
      outer = TRUE, line = 1.5,cex.main=1.5,font.main=2)

dev.off()

clades[1]<-"globalTree_"
assign(paste(clades[q],param[z],"_10tree_1Ma",sep=""),join4runs_ALL)
do.call(save,list(paste(clades[q],param[z],"_10tree_1Ma",sep=""),file=paste("MAMPHY_10tree_CoMET_summary_1Ma_",param[z],clades[q],".Rda",sep="")))



} # close multiple models loop (4)


#} # close multi- clade loop (10 trees each)



#KUHN_10tree_CoMET<-join4runs_ALL
KUHN_10tree_CoMET_5Ma<-join4runs_ALL
#FS_10tree_CoMET<-join4runs_ALL
#FS_10tree_CoMET_5Ma<-join4runs_ALL
MAMPHY_10tree_CoMET_5Ma<-join4runs_ALL
MAMPHY_10tree_CoMET_5Ma_7shiftnoME<-join4runs_ALL
MAMPHY_10tree_CoMET_1Ma_7shiftnoME<-join4runs_ALL

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/")
save(MAMPHY_10tree_CoMET_1Ma_7shiftnoME,file="MAMPHY_10tree_CoMET_summary_1Ma_7shiftnoME.Rda")

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/")
save(MAMPHY_10tree_CoMET_5Ma_7shiftnoME,file="MAMPHY_10tree_CoMET_summary_5Ma_7shiftnoME.Rda")

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/")
load(file="MAMPHY_10tree_CoMET_summary_5Ma.Rda")

save(MAMPHY_10tree_CoMET_5Ma,file="MAMPHY_10tree_CoMET_summary_5Ma.Rda")

save(KUHN_10tree_CoMET,file="KUHN_10tree_CoMET_summary.Rda")
save(KUHN_10tree_CoMET_5Ma,file="KUHN_10tree_CoMET_summary_5Ma.Rda")
save(FS_10tree_CoMET,file="FS_10tree_CoMET_summary.Rda")
save(FS_10tree_CoMET_5Ma,file="FS_10tree_CoMET_summary_5Ma.Rda")

#save.image(file="workspace_CoMET_v10tree_summaries.Rda")

#} # close multiple trees loop

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/")
#load(file="MAMPHY_10tree_CoMET_summary_5Ma.Rda")
previousTree<-"MamPhy"
join4runs_ALL<-MAMPHY_10tree_CoMET_5Ma

#load(file="KUHN_10tree_CoMET_summary_5Ma.Rda")
previousTree<-"KuhnEtAl2011"
join4runs_ALL<-KUHN_10tree_CoMET_5Ma

#load(file="FS_10tree_CoMET_summary_5Ma.Rda")
previousTree<-"FS2015"
join4runs_ALL<-FS_10tree_CoMET_5Ma


pdf(file=paste(previousTree,"_CoMET_estimHyper_10trees_minESS500_results_join4runs_ALL10trees_5Ma_justMassExt_wTree.pdf",sep=""),width=8,height=10)

 tess.plot.output2(join4runs_ALL, plot.tree=TRUE,
  fig.types = c("mass extinction times"), las=2)

mtext(text=paste(previousTree," (10 trees, 4 runs each)",sep=""),side=3)
dev.off()

# all three
pdf(file=paste("MamPhy-vs-others_CoMET_estimHyper_10trees_minESS500_results_join4runs_ALL10trees_5Ma_justMassExt_wTree.pdf",sep=""),width=7,height=10)

layout.mat <- matrix(1:3,nrow=3,ncol=1,byrow=TRUE)
layout(layout.mat)

previousTree<-"MamPhy"
join4runs_ALL<-MAMPHY_10tree_CoMET_5Ma
tess.plot.output2(join4runs_ALL, plot.tree=TRUE,
  fig.types = c("mass extinction times"), las=2)
mtext(text=paste(previousTree," (10 trees, 4 runs each, 5 Ma bins)",sep=""),side=3)

previousTree<-"KuhnEtAl2011"
join4runs_ALL<-KUHN_10tree_CoMET_5Ma
tess.plot.output2(join4runs_ALL, plot.tree=TRUE,
  fig.types = c("mass extinction times"), las=2)
mtext(text=paste(previousTree," (10 trees, 4 runs each, 5 Ma bins)",sep=""),side=3)

previousTree<-"FS2015"
join4runs_ALL<-FS_10tree_CoMET_5Ma
tess.plot.output2(join4runs_ALL, plot.tree=TRUE,
  fig.types = c("mass extinction times"), las=2)
mtext(text=paste(previousTree," (10 trees, 4 runs each, 5 Ma bins)",sep=""),side=3)

dev.off()






####
## Calc the MEAN and 95% CI for the inferred extinction divTimes...
#####

output_list_ALL <- list(join4runs_tree1, join4runs_tree2, join4runs_tree3, join4runs_tree4, join4runs_tree5, join4runs_tree6, join4runs_tree7, join4runs_tree8, join4runs_tree9, join4runs_tree10)
type<-"mass extinction times"

res<-vector("list",length=10)
for (i in 1:10){
  output<-output_list_ALL[[i]]
  numIntervals<-length(output[[5]][1,])
  intervalTimes<-seq(0,(numIntervals*5),by=5)
  thisOutput <- output[[type]]
  
  colnames(thisOutput)<-rev(intervalTimes)[2:length(intervalTimes)]
  
  sumThisOutput <- rev(colSums(thisOutput))[1:31] # takes only from 0-150 Mya

  allVals<-vector("list",length(sumThisOutput))
    for(j in 1:length(sumThisOutput)){
      val<-as.numeric(names(sumThisOutput[j]))
      allVals[[j]]<-rep(val,sumThisOutput[[j]])
    }
  freqDist<-do.call(c,allVals)

  Mean<-mean(freqDist)
  quants<-quantile(freqDist, c(0.025,0.5,0.975))
  Median<-quants[[2]]
  low95<-quants[[1]]
  up95<-quants[[3]]

  res[[i]]<-c(Mean,Median,low95,up95)
}
resAll<-do.call(rbind,res)
colnames(resAll)<-c("Mean","Median","low95","up95")

toWrite<-rbind(resAll,colMeans(resAll))
rownames(toWrite)<-c(1:10,"Means")

write.table(toWrite,file=paste(previousTree,"_CoMET_sumAcrossTrees_all10.txt",sep=""))

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/")
previousTree<-"MamPhy"
toWriteMAM<-read.table(file=paste(previousTree,"_CoMET_sumAcrossTrees_all10.txt",sep=""))

previousTree<-"KuhnEtAl2011"
toWriteKUHN<-read.table(file=paste(previousTree,"_CoMET_sumAcrossTrees_all10.txt",sep=""))







######
## Diagnostics
# Compute the effective sample size (ESS) and Geweke diagnostic
# Higher ESS values should provide more precise inferences from the posterior sample. As
# a rule of thumb, the ESS should be larger than 200.
 
# for the number of speciation-rate shifts.
effectiveSize(output$numSpeciationCategories)
geweke.diag(output$numSpeciationCategories)
geweke.diag(output$numExtinctionCategories)
# for the number of mass-extinctionevents.
effectiveSize(output$numMassExtinctions)
geweke.diag(output$numMassExtinctions)
 
 
# plot
 
pdf(paste("C:\\Woidl\\AlexPyron\\Amphibians\\MS\\tess_analysis\\CoMET",Taxon,i,"emphyper_100k_assess.pdf",sep="_"))
 

 
 









# non-empirical priors w CoMET
######
 #Specify the mean and standard deviation of the lognormal
# prior on the speciation rate in real space
speciationPriorMu <- 0.2
speciationPriorSigma <- 0.5
 
# Specify the mean and standard deviation of the lognormal
# prior on the extinction rate in real space
extinctionPriorMu <- 0.15
extinctionPriorSigma <- 0.5

#Transform the priors on the speciation rate into log space.
speciationRatePriorMean <- log((speciationPriorMu^2)/sqrt(speciationPriorSigma^2+ speciationPriorMu^2))
speciationRatePriorStDev <- sqrt( log(1+speciationPriorSigma^2 /(speciationPriorMu^2)))
# Transform the priors on the extinction rate into log space.
extinctionRatePriorMean <- log((extinctionPriorMu^2) /sqrt(extinctionPriorSigma^2+extinctionPriorMu^2))
extinctionRatePriorStDev <- sqrt( log(1+extinctionPriorSigma^2 /(extinctionPriorMu^2)))
 
tess.analysis(mamMCC,
  empiricalHyperPriors = FALSE,
  initialSpeciationRate = speciationPriorMu,
  speciationRatePriorMean = speciationRatePriorMean,
  speciationRatePriorStDev = speciationRatePriorStDev,
 initialExtinctionRate = extinctionPriorMu,
  extinctionRatePriorMean = extinctionRatePriorMean,
  extinctionRatePriorStDev = extinctionRatePriorStDev,
  samplingProbability = samplingFraction,
  numExpectedRateChanges = numExpectedRateChanges,
  numExpectedMassExtinctions = numExpectedMassExtinctions,
  pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
  pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
  MAX_ITERATIONS = 1000,
dir = "tess_analysis")
 
output <- tess.process.output("tess_analysis",
       numExpectedRateChanges = numExpectedRateChanges,
       numExpectedMassExtinctions = numExpectedMassExtinctions)
 
layout.mat <- matrix(1:6,nrow=3,ncol=2,byrow=TRUE)
layout(layout.mat)
tess.plot.output(output,
  fig.types = c("speciation rates",
  "speciation shift times",
  "extinction rates",
  "extinction shift times",
  "mass extinction Bayes factors",
  "mass extinction times"),
las=2)


 



## Constant rates:
prior_delta <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_tau <- function(x) { dexp(x,rate=10.0,log=TRUE) }
priorsConstBD <- c("diversification"=prior_delta, "turnover"=prior_tau)
 
likelihoodConstBD <- function(params) {
  speciation <- params[1] + params[2]
  extinction <- params[2]
  lnl <- tess.likelihood(times,
  lambda = speciation,
  mu = extinction,
  samplingProbability = 1.0,
  log = TRUE)
  return (lnl)
}
 
samplesConstBD <- tess.mcmc(likelihoodFunction = likelihoodConstBD,
  priors = priorsConstBD,
  parameters = runif(2,0,1),
  logTransforms = c(TRUE,TRUE),
  delta = c(1,1),
  iterations = 10000,
  burnin = 1000,
  thinning = 10,
  adaptive = TRUE,
  verbose = TRUE)
 
summary(samplesConstBD)
plot(samplesConstBD)
  
## Varying rates: exponentially decreasing speciation rate
# 2.4.2 Birth-death processes with continuously varying rates
prior_delta <- function(x) { dexp(x,rate=0.1,log=TRUE) }
prior_lambda <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_alpha <- function(x) { dexp(x,rate=0.1,log=TRUE) }
priorsDecrBD <- c("turnover"=prior_delta,
  "initial speciation"=prior_lambda,
  "speciation decay"=prior_alpha)
 
likelihoodDecrBD <- function(params) {
  speciation <- function(t) params[1] + params[2] * exp(-params[3]*t)
  extinction <- function(t) params[1]
  lnl <- tess.likelihood(times,
  lambda = speciation,
  mu = extinction,
  samplingProbability = 1.0,
  log = TRUE)
  return (lnl)
}
 
samplesDecrBD <- tess.mcmc(likelihoodFunction = likelihoodDecrBD,
  priors = priorsDecrBD,
  parameters = runif(3,0,1),
  logTransforms = c(TRUE,TRUE,TRUE),
  delta = c(1,1,1),
  iterations = 10000,
  burnin = 1000,
  thinning = 10,
  adaptive = TRUE,
  verbose = TRUE)

summary(samplesConstBD) # slow in running
plot(samplesConstBD)
  

# 2.4.3 Birth-death processes with episodically varying rates
rateChangeTime <- max( times ) / 2

prior_delta_before <- function(x) { dexp(x,rate=10.0,log=TRUE) } prior_tau_before <- function(x) { dexp(x,rate=10.0,log=TRUE) } prior_delta_after <- function(x) { dexp(x,rate=10.0,log=TRUE) } prior_tau_after <- function(x) { dexp(x,rate=10.0,log=TRUE) } priorsEpisodicBD <- c("diversification before"=prior_delta_before,
                      "turnover before"=prior_tau_before,
                      "diversification after"=prior_delta_after,
                      "turnover after"=prior_tau_after)

likelihoodEpisodicBD <- function(params) {
speciation <- c(params[1]+params[2],params[3]+params[4])
extinction <- c(params[2],params[4])
lnl <- tess.likelihood.rateshift(times,
                                   lambda = speciation,
                                   mu = extinction,
                                   rateChangeTimesLambda = rateChangeTime,
                                   rateChangeTimesMu = rateChangeTime,
                                   samplingProbability = 1.0,
                                   log=TRUE)
    return (lnl) 
    } 

samplesEpisodicBD <- tess.mcmc(likelihoodFunction = likelihoodEpisodicBD,
                               priors = priorsEpisodicBD,
                               parameters = runif(4,0,1),
                               logTransforms = c(TRUE,TRUE,TRUE,TRUE),
                               delta = c(1,1,1,1),
                               iterations = 10000,
                               burnin = 1000,
                               thinning = 10,
                               adaptive = TRUE,
                               verbose = TRUE)

summary(samplesEpisodicBD)
plot(samplesEpisodicBD)



##
## Explicit Mass extinction
# 2.4.4 Birth-death processes with explicit mass-extinction events 

survivalProbability <- 0.1
 
prior_delta <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_tau <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_time <- function(x) { dunif(x,min=min(times)/4,max=max(times),log=TRUE)}
priorsMassExtinctionBD <- c("diversification"=prior_delta,
  "turnover"=prior_tau,
  "mass-extinction time"=prior_time)
 
likelihoodMassExtinctionBD <- function(params) {
  speciation <- params[1]+params[2]
  extinction <- params[2]
  time <- params[3]
  lnl <- tess.likelihood(times,
  lambda = speciation,
  mu = extinction,
  massExtinctionTimes = time,
  massExtinctionSurvivalProbabilities =
  survivalProbability,
  samplingProbability = 1.0,
  log = TRUE)
  return (lnl)
}
 
samplesMassExtinctionBD <- tess.mcmc(likelihoodFunction =
  likelihoodMassExtinctionBD,
  priors = priorsMassExtinctionBD,
  parameters = c(runif(2,0,1),max(times)*3/4),
  logTransforms = c(TRUE,TRUE,FALSE),
  delta = c(1,1,1),
  iterations = 1000000,
  burnin = 10000,
  thinning = 10,
  adaptive = TRUE,
  verbose = TRUE)


 
summary(samplesMassExtinctionBD)
plot(samplesMassExtinctionBD)

# 100,000 gens...
# Parameter | delta | Acceptance Probability
# ==========================================
# diversification         | 0.059         | 0.449
# turnover                | 0.087         | 0.439
# mass-extinction time            | 0.475         | 0.442

# Iterations = 1:100001
# Thinning interval = 1 
# Number of chains = 1 
# Sample size per chain = 100001 
# 
# 1. Empirical mean and standard deviation for each variable,
#    plus standard error of the mean:
# 
#                           Mean       SD  Naive SE Time-series SE
# diversification        0.09886 0.003646 1.153e-05      1.701e-05
# turnover               0.15571 0.008151 2.577e-05      3.513e-05
# mass-extinction time 130.05399 1.223627 3.869e-03      2.942e-02
# 
# 2. Quantiles for each variable:
# 
#                           2.5%      25%       50%      75%    97.5%
# diversification        0.09176   0.0964   0.09885   0.1013   0.1060
# turnover               0.13990   0.1502   0.15564   0.1612   0.1718
# mass-extinction time 128.35168 129.2870 129.46594 131.4310 132.4992






# WALTERS CODE::
#####

####################
##### Div analysis
 library(TESS)
 
setwd("C:\\Woidl\\AlexPyron\\Amphibians\\MS")
AmpED <- read.table("Amphibian_Data_20May.csv", header = TRUE, sep = ",")
 
 
#### Select trees:
# species taxonomic assignments, arranged in a table
#speciesClades<-read.table("C:\\Woidl\\AlexPyron\\Amphibians\\2014\\amphibian_species.csv", header = TRUE, sep = ",")
#head(speciesClades)
 
all <- AmpED$Scientific_Name
anu <- AmpED[which(AmpED$Taxon=="Anura"),"Scientific_Name"]
caud <- AmpED[which(AmpED$Taxon=="Caudata"),"Scientific_Name"]
gym <- AmpED[which(AmpED$Taxon=="Gymnophiona"),"Scientific_Name"]
 
notanu <- setdiff(all,anu)
notcaud <- setdiff(all,caud)
notgym <- setdiff(all,gym)
 
gym.tree <- drop.tip(tree,notgym)
caud.tree <- drop.tip(tree,notcaud)
anu.tree <- drop.tip(tree,notanu)
 
 
times <- as.numeric( branching.times(conifers) )
times <- as.numeric( branching.times(gym.tree))
times <- as.numeric( branching.times(caud.tree))
times <- as.numeric( branching.times(anu.tree))
 
 
##
## Constant rates:
prior_delta <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_tau <- function(x) { dexp(x,rate=10.0,log=TRUE) }
priorsConstBD <- c("diversification"=prior_delta, "turnover"=prior_tau)
 
likelihoodConstBD <- function(params) {
  speciation <- params[1] + params[2]
  extinction <- params[2]
  lnl <- tess.likelihood(times,
  lambda = speciation,
  mu = extinction,
  samplingProbability = 1.0,
  log = TRUE)
  return (lnl)
}
 
samplesConstBD <- tess.mcmc(likelihoodFunction = likelihoodConstBD,
  priors = priorsConstBD,
  parameters = runif(2,0,1),
  logTransforms = c(TRUE,TRUE),
  delta = c(1,1),
  iterations = 10000,
  burnin = 1000,
  thinning = 10,
  adaptive = TRUE,
  verbose = TRUE)
 
summary(samplesConstBD)
plot(samplesConstBD)
 
samplesConstBD.gym <- samplesConstBD
samplesConstBD.caud <- samplesConstBD
samplesConstBD.anu <- samplesConstBD
 
##
## Varying rates: exponentially decreasing speciation rate
 
prior_delta <- function(x) { dexp(x,rate=0.1,log=TRUE) }
prior_lambda <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_alpha <- function(x) { dexp(x,rate=0.1,log=TRUE) }
priorsDecrBD <- c("turnover"=prior_delta,
  "initial speciation"=prior_lambda,
  "speciation decay"=prior_alpha)
 
likelihoodDecrBD <- function(params) {
  speciation <- function(t) params[1] + params[2] * exp(-params[3]*t)
  extinction <- function(t) params[1]
  lnl <- tess.likelihood(times,
  lambda = speciation,
  mu = extinction,
  samplingProbability = 1.0,
  log = TRUE)
  return (lnl)
}
 
 
samplesDecrBD <- tess.mcmc(likelihoodFunction = likelihoodDecrBD,
  priors = priorsDecrBD,
  parameters = runif(3,0,1),
  logTransforms = c(TRUE,TRUE,TRUE),
  delta = c(1,1,1),
  iterations = 10000,
  burnin = 1000,
  thinning = 10,
  adaptive = TRUE,
  verbose = TRUE)
 
##
## Explicit Mass extinction
 
survivalProbability <- 0.1
 
prior_delta <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_tau <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_time <- function(x) { dunif(x,min=min(times)/4,max=max(times),log=TRUE)}
priorsMassExtinctionBD <- c("diversification"=prior_delta,
  "turnover"=prior_tau,
  "mass-extinction time"=prior_time)
 
likelihoodMassExtinctionBD <- function(params) {
  speciation <- params[1]+params[2]
  extinction <- params[2]
  time <- params[3]
  lnl <- tess.likelihood(times,
  lambda = speciation,
  mu = extinction,
  massExtinctionTimes = time,
  massExtinctionSurvivalProbabilities =
  survivalProbability,
  samplingProbability = 1.0,
  log = TRUE)
  return (lnl)
}
 
 
 
samplesMassExtinctionBD <- tess.mcmc(likelihoodFunction =
  likelihoodMassExtinctionBD,
  priors = priorsMassExtinctionBD,
  parameters = c(runif(2,0,1),max(times)*3/4),
  logTransforms = c(TRUE,TRUE,FALSE),
  delta = c(1,1,1),
  iterations = 10000,
  burnin = 1000,
  thinning = 10,
  adaptive = TRUE,
  verbose = TRUE)
 
summary(samplesMassExtinctionBD)
plot(samplesMassExtinctionBD)
 
#################################
 
#####
##CoMet
# CPP on Mass-Extinction Times (CoMET) model (May et al., 2015)
 
#####
# set parameters
 
test.tree <- conifers
samplingFraction <- (conifers$Nnode + 1) / 630  # proportion of extant species in tree
 
test.tree <- caud.tree
test.tree <- gym.tree
test.tree <- anu.tree
 
samplingFraction <- 1
numExpectedMassExtinctions <- 2
numExpectedRateChanges <- 2
 
 
#Specify the mean and standard deviation of the lognormal
# prior on the speciation rate in real space
speciationPriorMu <- 0.2
speciationPriorSigma <- 0.5
 
# Specify the mean and standard deviation of the lognormal
# prior on the extinction rate in real space
extinctionPriorMu <- 0.15
extinctionPriorSigma <- 0.5
#Transform the priors on the speciation rate into log space.
speciationRatePriorMean <- log((speciationPriorMu^2)/sqrt(speciationPriorSigma^2+ speciationPriorMu^2))
speciationRatePriorStDev <- sqrt( log(1+speciationPriorSigma^2 /(speciationPriorMu^2)))
# Transform the priors on the extinction rate into log space.
extinctionRatePriorMean <- log((extinctionPriorMu^2) /sqrt(extinctionPriorSigma^2+extinctionPriorMu^2))
extinctionRatePriorStDev <- sqrt( log(1+extinctionPriorSigma^2 /(extinctionPriorMu^2)))
 
expectedSurvivalProbability <- 0.05
 
pMassExtinctionPriorShape2 <- 100
pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 * expectedSurvivalProbability / (expectedSurvivalProbability - 1)
 
# Plot the density function of our beta distribution.
curve(dbeta(x,shape1=pMassExtinctionPriorShape1, shape2=pMassExtinctionPriorShape2),n=1001, xlab= 'survival probability' ,ylab= 'density' ,las=1)
# Plot the 95% prior interval on the survival probability.
abline(v = qbeta(c(0.025,0.975),shape1=pMassExtinctionPriorShape1, shape2=pMassExtinctionPriorShape2),lty=2)
 
 
########
## CoMet
### Do CoMet analysis with these parameters
 
tess.analysis(test.tree,
  empiricalHyperPriors = FALSE,
  initialSpeciationRate = speciationPriorMu,
  speciationRatePriorMean = speciationRatePriorMean,
  speciationRatePriorStDev = speciationRatePriorStDev,
 initialExtinctionRate = extinctionPriorMu,
  extinctionRatePriorMean = extinctionRatePriorMean,
  extinctionRatePriorStDev = extinctionRatePriorStDev,
  samplingProbability = samplingFraction,
  numExpectedRateChanges = numExpectedRateChanges,
  numExpectedMassExtinctions = numExpectedMassExtinctions,
  pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
  pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
  MAX_ITERATIONS = 1000,
dir = "tess_analysis")
 
output <- tess.process.output("tess_analysis",
       numExpectedRateChanges = numExpectedRateChanges,
       numExpectedMassExtinctions = numExpectedMassExtinctions)
 
layout.mat <- matrix(1:6,nrow=3,ncol=2,byrow=TRUE)
layout(layout.mat)
tess.plot.output(output,
  fig.types = c("speciation rates",
  "speciation shift times",
  "extinction rates",
  "extinction shift times",
  "mass extinction Bayes factors",
  "mass extinction times"),
las=2)
 
output.anu1 <- output
 
 
#####
## CoMet
## Use Empirical Hyperpriors
 
test.tree <- caud.tree; Taxon <- "caud"
test.tree <- gym.tree; Taxon <- "gym"
test.tree <- anu.tree; Taxon <- "anu"
 
for (i in 1:10) {  #start multiple trees loop
print(i)
tree <- trees10[[i]]
 
tess.analysis(test.tree,
  empiricalHyperPriors = TRUE,
  samplingProbability = samplingFraction,
  numExpectedRateChanges = numExpectedRateChanges,
  numExpectedMassExtinctions = numExpectedMassExtinctions,
  pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
  pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
  MAX_ITERATIONS = 100000,
  dir = "comet_hyperpriors")
 
output <- tess.process.output("comet_hyperpriors",
  numExpectedRateChanges = numExpectedRateChanges,
  numExpectedMassExtinctions = numExpectedMassExtinctions)
 
assign(paste("output.",Taxon,i,".emphyp.100k",sep="."),output)
 
#output.caud1.emphyp.200k <- output
 
pdf(paste("C:\\Woidl\\AlexPyron\\Amphibians\\MS\\tess_analysis\\CoMET",Taxon,i,"emphyper_100k.pdf",sep="_"))
 
layout.mat <- matrix(1:6,nrow=3,ncol=2,byrow=TRUE)
layout(layout.mat)
  tess.plot.output(output,
  fig.types = c("speciation rates",
  "speciation shift times",
  "extinction rates",
  "extinction shift times",
  "mass extinction Bayes factors",
  "mass extinction times"),
  las=2)
 
dev.off()
 
 
 
######
## Diagnostics
# Compute the effective sample size (ESS) and Geweke diagnostic
# Higher ESS values should provide more precise inferences from the posterior sample. As
# a rule of thumb, the ESS should be larger than 200.
 
# for the number of speciation-rate shifts.
effectiveSize(output$numSpeciationCategories)
geweke.diag(output$numSpeciationCategories)
geweke.diag(output$numExtinctionCategories)
# for the number of mass-extinctionevents.
effectiveSize(output$numMassExtinctions)
geweke.diag(output$numMassExtinctions)
 
 
# plot
 
pdf(paste("C:\\Woidl\\AlexPyron\\Amphibians\\MS\\tess_analysis\\CoMET",Taxon,i,"emphyper_100k_assess.pdf",sep="_"))
 
layout.mat <- matrix(1:6,nrow=3,ncol=2,byrow=TRUE)
layout(layout.mat)
  tess.plot.singlechain.diagnostics(output,
  parameters = c("speciation rates",
  "extinction rates",
  "mass extinction times"),
  las=2)
 
dev.off()
 
 
} # close multiple trees loop
 
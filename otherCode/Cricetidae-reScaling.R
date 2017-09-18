library(ape)
library(phytools)
library(phyloch)

#####
# BD, with rateVar, and ONLY estimating the Sigmos + Tylos + Rattus -- as FREE of constraints too.
#####
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/9_patch-PASTIS-ing/_BIRTH-DEATH_varRates_fixedTopos/__TO_upload/PC12a_Sigmodontinae_BDfree_33M")

rescaleTree<-function(tree,scale){
	tree$edge.length<-tree$edge.length/max(nodeHeights(tree)[,2])*scale
	return(tree)
	}

sigmo10k <- read.nexus("postburn10k_PC12a_Sigmodontinae_BDfree_33M.trees")


sigmo10k_Mer<-list()
for(i in 1:length(sigmo10k)){
	sigmo10k_Mer[[i]] <- rescaleTree(sigmo10k[[i]],scale=25.88) # this is the MEAN GLOBAL div time between Cricetus and Mus
}

write.nexus(sigmo10k_Mer, file="postburn10k_PC12a_Sigmodontinae_BDfree_33M_rescaleMer11.trees")



#####
# BD, with rateVar, and thomasomyines fixed
#####
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/9_patch-PASTIS-ing/_BIRTH-DEATH_varRates/PC12_Cricetidae_BD_HPC_re2_33M-X")

rescaleTree<-function(tree,scale){
	tree$edge.length<-tree$edge.length/max(nodeHeights(tree)[,2])*scale
	return(tree)
	}

cricetid10k <- read.nexus("postburn10k_PC12_Cricetidae_BD.trees")


cricetid10k_Mer<-list()
for(i in 1:length(cricetid10k)){
	cricetid10k_Mer[[i]] <- rescaleTree(cricetid10k[[i]],scale=25.88) # this is the MEAN GLOBAL div time between Cricetus and Mus
}

#write.nexus(cricetid10k_Mer, file="postburn10k_PC12_Cricetidae_BD_rescaleMer11-static.trees")


FBDnofossils<-read.nexus("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/10_tipDating/BbonePASTIS_FIN4_FBD_156taxa_re4_50M-X/postburn10k_Bbone_FBD_156taxa_ALL_inMYR_pruneBP_noFossils_53taxa_nex.trees")

FBDnofossils[[1]]$edge.length[getMRCA(FBDnofossils[[1]], c("Cricetulus_barabensis__CRICETIDAE__RODENTIA","Rattus_norvegicus__MURIDAE__RODENTIA"))]

divTimes<-vector(length=10000)
for(i in 1:length(FBDnofossils)){
	node<-as.character(getMRCA(FBDnofossils[[i]], c("Cricetulus_barabensis__CRICETIDAE__RODENTIA","Rattus_norvegicus__MURIDAE__RODENTIA")))
	divTimes[i]<-round(branching.times(FBDnofossils[[i]]),2)[node]
}

plot(density(divTimes))

abline(median(divTimes))
[1] 31.9

abline(quantile(divTimes, c(0.025,0.975)))
    2.5%    97.5% 
24.12975 48.62050 

## COOL.

cricetid10k_FBD<-list()
for(i in 1:length(cricetid10k)){
	cricetid10k_FBD[[i]] <- rescaleTree(cricetid10k[[i]],scale=divTimes[i]) # this is EACH of the 10k FBD trees as mrca of Rattus and Cricetidae
}

write.nexus(cricetid10k_FBD, file="postburn10k_PC12_Cricetidae_BD_rescaleFBD.trees")

## Prune out all the unsampled Cricetidae -- 530 taxa

genCricetids<-read.table("PC12_Cricetidae_geneticSp-530.txt")

toDrop<-setdiff(cricetid10k_FBD[[1]]$tip.label,genCricetids$V1) # these are all the non-Sigmos

cricetid10k_FBD_530<-list()
for(i in 1:length(cricetid10k_FBD)){
	cricetid10k_FBD_530[[i]] <- drop.tip(cricetid10k_FBD[[i]],toDrop)
}

write.nexus(cricetid10k_FBD_530, file="postburn10k_PC12_Cricetidae_BD_rescaleFBD_530taxa.trees")



## Prune out all but the Sigmos -- all 413 taxa

allSigs<-read.table("PC12_Cricetidae_Sigmos_allSp-413.txt")

toDrop2<-setdiff(cricetid10k_FBD[[1]]$tip.label,allSigs$V1) # these are all the non-Sigmos

sigmo10k_FBD_413<-list()
for(i in 1:length(cricetid10k_FBD)){
	sigmo10k_FBD_413[[i]] <- drop.tip(cricetid10k_FBD[[i]],toDrop2)
}

write.nexus(sigmo10k_FBD_413, file="postburn10k_PC12_Cricetidae_BD_rescaleFBD_Sigmos_413taxa.trees")

## Prune out all but the Sigmos -- 297 taxa

sampledSigs<-read.table("PC12_Cricetidae_Sigmos_geneticSp-279.txt")

toDrop3<-setdiff(cricetid10k_FBD[[1]]$tip.label,sampledSigs$V1) # these are all the non-Sigmos + unsampled Sigmos

sigmo10k_FBD_279<-list()
for(i in 1:length(cricetid10k_FBD)){
	sigmo10k_FBD_279[[i]] <- drop.tip(cricetid10k_FBD[[i]],toDrop3)
}

write.nexus(sigmo10k_FBD_279, file="postburn10k_PC12_Cricetidae_BD_rescaleFBD_Sigmos_279taxa.trees")

#########
# FINAL work on the sigmodtinae tree....
######
# OK wait, go back to these FBD rescaled posteriors and see whats up...

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/9_patch-PASTIS-ing/_BIRTH-DEATH_varRates/PC12_Cricetidae")

SigmoFBD<-read.nexus("postburn10k_PC12_Cricetidae_BD_rescaleFBD_Sigmos_413taxa.trees")

#>> These trees all VARY in their ROOT age... because I used the distribution of root ages from the FBD to rescale the distribution of patch trees... non-bueno

# Instead-- rescale all to 25.88, then prune out the non-sigmos.

cricetid10k <- read.nexus("postburn10k_PC12_Cricetidae_BD.trees")

cricetid10k_Mer<-list()
for(i in 1:length(cricetid10k)){
	cricetid10k_Mer[[i]] <- rescaleTree(cricetid10k[[i]],scale=25.88) # this is the MEAN GLOBAL div time between Cricetus and Mus
}

write.nexus(cricetid10k_Mer, file="postburn10k_PC12_Cricetidae_BD_rescaleMer11-static.trees")

# now prune to just the Sigmos.

allSigs<-read.table("PC12_Cricetidae_Sigmos_allSp-413.txt")

toDrop2<-setdiff(cricetid10k_Mer[[1]]$tip.label,allSigs$V1) # these are all the non-Sigmos

sigmo10k_Mer_413<-list()
for(i in 1:length(cricetid10k_Mer)){
	sigmo10k_Mer_413[[i]] <- drop.tip(cricetid10k_Mer[[i]],toDrop2)
}

write.nexus(sigmo10k_Mer_413, file="postburn10k_PC12_Cricetidae_BD_rescaleMer11-static_Sigmos_413taxa.trees")

# Now grab out 100 of these

randomSig413_100<-sample(sigmo10k_Mer_413,size=100)

# get a distribution of those root ages:

rootAges<-vector()
for (i in 1:length(randomSig413_100)){
	rootAges[i]<-max(randomSig413_100[[i]][[2]])
}

plot(density(rootAges))

quantile(rootAges, c(0.025,0.975))
     2.5%     97.5% 
 8.336342 10.188662 

rootAges<-vector()
for (i in 1:length(sigmo10k_Mer_413)){
	rootAges[i]<-max(sigmo10k_Mer_413[[i]][[2]])
}

plot(density(rootAges))

quantile(rootAges, c(0.025,0.975))
     2.5%     97.5% 
 8.178472 10.719992 
> mean(rootAges)
[1] 9.346108
> median(rootAges)
[1] 9.301606

# >> That is all TOO YOUNG though, compared to the 12.7 Ma crown age that I derived from the MCC of the free analyses PC12a...

# So then rescale to THAT:
# 12.7 (11.35,14.18)

# First prune the UNSCALED Cricetidae tree:

allSigs<-read.table("PC12_Cricetidae_Sigmos_allSp-413.txt")

toDrop2<-setdiff(cricetid10k_Mer[[1]]$tip.label,allSigs$V1) # these are all the non-Sigmos

sigmo10k_unscaled<-list()
for(i in 1:length(cricetid10k)){
	sigmo10k_unscaled[[i]] <- drop.tip(cricetid10k[[i]],toDrop2)
	}

# Now get just the root Sigmodontine times from the PC12a posterior of 10k trees

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/9_patch-PASTIS-ing/_BIRTH-DEATH_varRates_FREE-topos/PC12a_Sigmodontinae_BDfree_33M")

PC12a_MerAges<-read.nexus("postburn10k_PC12a_Sigmodontinae_BDfree_33M_rescaleMer11.trees") # the freely estimated DNA only 10k, as rescaled to Mer 11

divTimes<-vector(length=10000)
for(i in 1:length(PC12a_MerAges)){
	node<-as.character(getMRCA(PC12a_MerAges[[i]], c("Sigmodon_alstoni_CRICETIDAE_RODENTIA","Akodon_sylvanus_CRICETIDAE_RODENTIA")))
	divTimes[i]<-round(branching.times(PC12a_MerAges[[i]]),2)[node]
}

plot(density(divTimes))

quantile(divTimes, c(0.025,0.975))
2.5% 97.5% 
11.38 14.22 

mean(divTimes)
[1] 12.69634 # BINGO.

median(divTimes)
[1] 12.66

# Now use those DIVTIMES to rescale each root in the sigmo10k_unscaled...

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/9_patch-PASTIS-ing/_BIRTH-DEATH_varRates/PC12_Cricetidae")


sigmo10k_12p7rootDist<-list()
for(i in 1:length(sigmo10k_unscaled)){
	sigmo10k_12p7rootDist[[i]] <- rescaleTree(sigmo10k_unscaled[[i]],scale=divTimes[i]) # this is EACH of the 10k FBD trees as mrca of Rattus and Cricetidae
}

write.nexus(sigmo10k_12p7rootDist, file="postburn10k_PC12_Cricetidae_BD_Sigmos413_rescale-12p7rootDist.trees")


# NOW sample 100 of those trees

sigmo10k_12p7rootDist_rand100<-sample(sigmo10k_12p7rootDist,size=100)

###
#
# USING this... !! 6 Oct 2016

phy12p7<-read.nexus("postburn10k_PC12_Cricetidae_BD_Sigmos413_rescale-12p7rootDist.trees")

phy12p7_100<-sample(phy12p7, size=100)

ltt.plot(phy12p7_100[[1]], log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:length(phy12p7_100)) {
	ltt.lines(phy12p7_100[[i]], col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}

write.nexus(phy12p7_100, file="postburn10k_PC12_Cricetidae_BD_Sigmos413_rescale-12p7rootDist_sample100.trees")

#
#
#####



# get a distribution of those root ages:

rootAges<-vector()
for (i in 1:length(sigmo10k_12p7rootDist)){
	rootAges[i]<-max(sigmo10k_12p7rootDist[[i]][[2]])
}

plot(density(rootAges))

quantile(rootAges, c(0.025,0.975))
     2.5%     97.5% 
 7.662385 10.783707 

mean(rootAges)
[1] 9.114627
 
median(rootAges)
[1] 9.062463

# >>> those ROOTS are actually INACCURATE !! That is just the max edge length! Which is sorta arbitary and not indicative of a node height.








# OR, could just use a STATIC estimate of 12.7 for every root tree-- but that would be inconsistent with the other MCC tree
# Try it... might actually be MORE consistent...
# ** 16.5 (arbitrary) **

sigmo10k_16p5rootPoint<-list()
for(i in 1:length(sigmo10k_unscaled)){
	sigmo10k_16p5rootPoint[[i]] <- rescaleTree(sigmo10k_unscaled[[i]],scale=16.5)
	}

write.nexus(sigmo10k_16p5rootPoint, file="postburn10k_PC12_Cricetidae_BD_Sigmos413_rescale-16p5rootPoint.trees")

#max(sigmo10k_unscaled[[1]]$edge.length/max(nodeHeights(sigmo10k_unscaled[[1]])[,2])*16.5)

# Now again get 100 trees

sigmo10k_16p5rootPoint_rand100<-sample(sigmo10k_16p5rootPoint,size=100)

# get a distribution of those root ages:

rootAges<-vector()
for (i in 1:length(sigmo10k_16p5rootPoint)){
	rootAges[i]<-max(sigmo10k_16p5rootPoint[[i]][[2]])
}

plot(density(rootAges))

quantile(rootAges, c(0.025,0.975))
    2.5%    97.5% 
10.42229 13.47395 

mean(rootAges)
[1] 11.84606

median(rootAges)
[1] 11.80596

###
# ** 18.3 (divTime from Tylomyinae) **

sigmo10k_18p3rootPoint<-list()
for(i in 1:length(sigmo10k_unscaled)){
	sigmo10k_18p3rootPoint[[i]] <- rescaleTree(sigmo10k_unscaled[[i]],scale=18.3)
	}

write.nexus(sigmo10k_18p3rootPoint, file="postburn10k_PC12_Cricetidae_BD_Sigmos413_rescale-18p3rootPoint.trees")



# Now again get 100 trees

sigmo10k_18p3rootPoint_rand100<-sample(sigmo10k_18p3rootPoint,size=100)

# get a distribution of those root ages:

rootAges<-vector()
for (i in 1:length(sigmo10k_18p3rootPoint)){
	rootAges[i]<-max(sigmo10k_18p3rootPoint[[i]][[2]])
}

plot(density(rootAges))

quantile(rootAges, c(0.025,0.975))
    2.5%    97.5% 
11.55927 14.94384 
 
mean(rootAges)
[1] 13.13836
 
median(rootAges)
[1] 13.09389


###
# Get 100 random UNSCALED

sigmo10k_unscaled_rand100<-sample(sigmo10k_unscaled,size=100)


# ** 17.65 (arbitary...) **

sigmo10k_17p65rootPoint<-list()
for(i in 1:100){
	sigmo10k_17p65rootPoint[[i]] <- rescaleTree(sigmo10k_unscaled_rand100[[i]],scale=17.65)
	}

rootAges<-vector()
for (i in 1:length(sigmo10k_17p65rootPoint)){
	rootAges[i]<-max(sigmo10k_17p65rootPoint[[i]][[2]])
}

plot(density(rootAges))

quantile(rootAges, c(0.025,0.975))
    2.5%    97.5% 
11.21458 15.56419 
  
mean(rootAges)
[1] 12.69572
  
median(rootAges)
[1] 12.53915


write.nexus(sigmo10k_17p65rootPoint, file="postburn10k_PC12_Cricetidae_BD_Sigmos413_rescale-17p65rootPoint_100trees.trees")
# == FINAL.

###
# NOW, plot the LTT.

jpeg(file="LTT_Sigmos413_rescale-17p65rootPoint_100trees.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(sigmo10k_17p65rootPoint[[1]], log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:length(sigmo10k_17p65rootPoint)) {
	ltt.lines(sigmo10k_17p65rootPoint[[i]], col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}
title(main="Sigmodontinae, completed to 413 spp, 100 trees")

dev.off()

# Now get the GAMMA statistics

gammaStats<-vector()
Pvals<-vector()
for (i in 1: length(sigmo10k_17p65rootPoint)){
	phy<-sigmo10k_17p65rootPoint[[i]]
	gammaStats[i]<-gammaStat(phy)
	Pvals[i]<-2*(1 - pnorm(abs(gammaStat(phy)))) # two-tailed test
}

mean(gammaStats)
[1] -1.86826

quantile(gammaStats, c(0.025,0.975))
      2.5%      97.5% 
-2.8555698 -0.8022379 

mean(Pvals)
[1] 0.1066223

quantile(Pvals, c(0.025,0.975))
       2.5%       97.5% 
0.004297275 0.423316492 

# use gammatest()


phy<-read.nexus("postburn10k_PC12_Cricetidae_BD_Sigmos413_rescale-17p65rootPoint_100trees.trees")

stats<-matrix(NA,nrow=length(phy), ncol=2)
colnames(stats)<-c("gamma","p")

dd<-ltt(phy)
for (i in 1:length(phy)){
	stats[i,1]<-dd[[i]]$gamma
	stats[i,2]<-dd[[i]]$p
}

mean(stats[,"gamma"])
[1] -1.86826

qq<-quantile(stats[,"gamma"], c(0.025,0.975))
      2.5%      97.5% 
-2.8555697 -0.8022378 

mean(stats[,"p"])
[1] 0.1066224

quantile(stats[,"p"], c(0.025,0.975))
       2.5%       97.5% 
0.004297277 0.423316568 


plot(density(stats[,"gamma"]), main="gamma statistics on 100 completed trees")
abline(v=-1.645, lty=2, col="red")
abline(v=mean(stats[,"gamma"]), lty=1, lwd=2, col="black")
abline(v=qq[[1]], lty=2, col="black")
abline(v=qq[[2]], lty=2, col="black")



obj<-ltt95(phy, alpha=0.05, log=TRUE, method="lineages", mode="mean")

plot(obj)







######
######
# Now get the SAME but just with the MER11 re-scaling (same mean one) across all nodes...
##

sampledSigs<-read.table("PC12_Cricetidae_Sigmos_geneticSp-279.txt")

toDrop4<-setdiff(cricetid10k_Mer[[1]]$tip.label,sampledSigs$V1) # these are all the non-Sigmos + unsampled Sigmos

sigmo10k_Mer_279<-list()
for(i in 1:length(cricetid10k_Mer)){
	sigmo10k_Mer_279[[i]] <- drop.tip(cricetid10k_Mer[[i]],toDrop3)
}

write.nexus(sigmo10k_Mer_279, file="postburn10k_PC12_Cricetidae_BD_rescaleMer11_Sigmos_279taxa.trees")









######
# YULE
##########
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/9_patch-PASTIS-ing/PC12_Cricetidae_MB_1re3_33M-fin/")

rescaleTree<-function(tree,scale){
  tree$edge.length<-tree$edge.length/max(nodeHeights(tree)[,2])*scale
  return(tree)
  }

### ON ONE TREE ###
tree<-read.nexus("infile.nex.runs1-4_postburn_lastTree.tre")

tree2<-rescaleTree(tree,scale=25.88) # this is the MEAN GLOBAL div time between Cricetus and Mus

quartz(width=8.5, height=22) 
plot(tree2, cex=0.15, label.offset=0.5, x.lim=c(0,30), y.lim=c(0,730))
node.support(branching.times(tree2), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.4)
data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-12.5, mgp=c(0,0.2,0))

## NOW, do this just for the SIGMO part of the tree...
tree3<-extract.clade(tree2,1044) 

quartz(width=8.5, height=22) 
plot(tree3, cex=0.3, label.offset=0.5, x.lim=c(0,30), y.lim=c(0,400))
node.support(branching.times(tree3), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.4)
data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-12.5, mgp=c(0,0.2,0))

###
# Cool... but the GOAL of this is the what, 
1) have the MCC of the full tree posterior -- post-burnin, scaled to time -- with all species incl PASTIS (and lots of polytomies) 
2) have the MCC of the full tree posterior -- post-burnin, scaled to time -- as based ONLY on the species with genetic data 
3) the full posterior of 10k trees -- post-burnin, scaled to time -- with all species incl PASTIS
4) the full posterior of 10k trees -- post-burnin, scaled to time -- pruned to just the species with genetic data

sampled<-read.table("PC12_Cricetidae_sampledSp-530.txt")

toDrop<-setdiff(tree2$tip.label,sampled$V1)

tree2_genetic<-drop.tip(tree2,toDrop)

quartz(width=8.5, height=22) 
plot(tree2_genetic, cex=0.25, label.offset=0.5, x.lim=c(0,30), y.lim=c(0,530))
nodelabels(cex=0.3)
node.support(branching.times(tree2_genetic), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.4)
data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-12.5, mgp=c(0,0.2,0))

tree3_genetic<-extract.clade(tree2_genetic,782) #appropriate node for Sigmodontinae in the genetic only sampling tree 

quartz(width=8.5, height=22) 
plot(tree3_genetic, cex=0.3, label.offset=0.5, x.lim=c(0,30), y.lim=c(0,300))
node.support(branching.times(tree3_genetic), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.4)
data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-8.5, mgp=c(0,0.2,0))

###
# 1) MCC of the full tree posterior -- post-burnin, scaled to time -- with all species incl PASTIS (and lots of polytomies) 
###
MCC<-read.beast("infile.nex.con_MCC_25p_treeAnOut_corr.tre")
MCC_table<-read.beast.table("infile.nex.con_MCC_25p_treeAnOut_corr.tre")

MCCX<-read.beast("infile.nex.con_MCC_delUnsamp-530taxa_corr.tre")

MCC2<-rescaleTree(MCC,scale=25.88) # this is the MEAN GLOBAL div time between Cricetus and Mus

quartz(width=8.5, height=22) 
plot(MCC2, cex=0.15, label.offset=0.5, x.lim=c(0,30), y.lim=c(0,730))
node.support(branching.times(MCC2), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.4)
node.support(MCC2$posterior, mode="numbers", digits=2,pos="below", col = "darkgray", font=2, cex=0.4)
data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-12.5, mgp=c(0,0.2,0))

quartz(width=8.5, height=22) 
plot(MCC2, cex=0.15, label.offset=0.5, x.lim=c(0,30), y.lim=c(0,730))
nodelabels(cex=0.3)

# Node 1044 is Sigmodontinae
MCC3<-extract.clade(MCC2,1044) #dont use extract.clade2 as is causes segfault... though may be useful for just extracting the values and then plotting the node values in a sep OBJ

MCC3_data<-extract.clade2(MCC2,1044)
MCC3_data$posterior
#MCC3_posteriors2<-MCC_table[316:727,20] # for nodes 1044 - 1455 (7283 tips in original tree)

# OR, extract all Sigmos from MCC2 using drop.tip2...
allSigs<-read.table("PC12_Cricetidae_Sigmos_allSp-413.txt")

toDrop<-setdiff(MCC2$tip.label,allSigs$V1) # these are all the non-Sigmos

MCC3_Sigs<-drop.tip2(MCC2,toDrop) # sweet, and no segFault

# AND, extract all the sampled Sigmos from MCC2...
sampledSigs<-read.table("PC12_Cricetidae_Sigmos_geneticSp-279.txt")

toDrop<-setdiff(MCC2$tip.label,sampledSigs$V1) # these are all the non-Sigmos + unsampled Sigmos

write.table(toDrop,"PC12_Cricetidae_nonSigmos_unsampledSigmos-449.txt")

MCC3_Sigs_gen<-drop.tip2(MCC2,toDrop) # sweet, and no segFault -- awesome.

quartz(width=8.5, height=22) 
plot(MCC3_Sigs, cex=0.3, label.offset=0.5, x.lim=c(0,30), y.lim=c(0,400))
node.support(branching.times(MCC3_Sigs), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.4)
node.support(MCC3_Sigs$posterior, mode="numbers", digits=2,pos="below", col = "darkgrey", font=2, cex=0.4)
node.support(MCC3_Sigs$posterior, mode="numbers", digits=2,cutoff=0.90,pos="below", col = "blue", font=2, cex=0.4)
data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-10.5, mgp=c(0,0.2,0))

quartz(width=8.5, height=22) 
plot(MCC3_Sigs_gen, cex=0.4, label.offset=0.5, x.lim=c(0,30), y.lim=c(0,270))
node.support(branching.times(MCC3_Sigs_gen), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.4)
node.support(MCC3_Sigs_gen$posterior, mode="numbers", digits=2,pos="below", col = "darkgrey", font=2, cex=0.4)
node.support(MCC3_Sigs_gen$posterior, mode="numbers", digits=2,cutoff=0.90,pos="below", col = "blue", font=2, cex=0.4)
data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-8.5, mgp=c(0,0.2,0))


# save the MCC with all species - 413 taxa
write.tree(MCC3, "infile.nex.con_MCC_25p_treeAnOut_corr_sigmos_413taxa.newick")

new<-read.tree("infile.nex.con_MCC_25p_treeAnOut_corr_sigmos_413taxa.newick")
plot(new, cex=0.2, label.offset=0.5)
node.support(branching.times(new), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.2)
data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-10.5, mgp=c(0,0.2,0))

# Prune to just the species sampled with genetic data - 279 taxa I think
sampled<-read.table("PC12_Cricetidae_sampledSp.txt")

toDrop<-setdiff(MCC3$tip.label,sampled$V1)

MCC3_gen<-drop.tip(MCC3,toDrop)

MCC3_gen_data<-drop.tip(MCC3_data,toDrop)


#####
# Get the RESCALED MCC tree for SIGMOS
###

MCCsig<-read.nexus("infile.nex.con_MCC_delNonUnsampSigs-279taxa_SIGMO.tre")

MCCsig2<-rescaleTree(MCCsig,scale=25.88) # this is the MEAN GLOBAL div time between Cricetus and Mus

plot(MCCsig2, cex=0.3, label.offset=0.5)
node.support(branching.times(MCCsig2), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.2)
data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-10.5, mgp=c(0,0.2,0))

quartz(width=8.5, height=22) 
plot(MCCsig2, cex=0.4, label.offset=0.5, x.lim=c(0,30), y.lim=c(0,270))
node.support(branching.times(MCCsig2), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.4)


##
###########
# DO for the FULL POSTERIORS now.
#####
# Read in the posterior sample of trees for CRICETIDAE (10,000)

sigUn10k<-read.nexus("infile.nex.runs1-4_postburn.trees")

cricetid10k<-vector("list",10000)
for(i in 1:length(cricetid10k)){
	cricetid10k[[i]]<-rescaleTree(sigUn10k[[i]],scale=25.88) # this is the MEAN GLOBAL div time between Cricetus and Mus
	}

for(i in 1:length(cricetid10k)){
	write.tree(cricetid10k[[i]], "infile.nex.runs1-4_postburn_reScaled-CRICETIDAE-all.trees", append=TRUE)
	}
	
# cricetid10k_c<-.compressTipLabel(cricetid10k)
# TRICK, then, was to just use write.tree rather than nexus format... but the tip.labels will take up more space :(

#####
# Prune those re-scaled trees to just the genetic sampled taxa

sampled<-read.table("PC12_Cricetidae_sampledSp-530.txt")

toDrop<-setdiff(cricetid10k[[1]]$tip.label,sampled$V1)

cricetid10k_gen<-vector("list",10000)
for(i in 1:length(cricetid10k_gen)){
	cricetid10k_gen[[i]]<-drop.tip(cricetid10k[[i]],toDrop)
	}

for(i in 1:length(cricetid10k_gen)){
	write.tree(cricetid10k_gen[[i]], "infile.nex.runs1-4_postburn_reScaled-CRICETIDAE-sampled.trees", append=TRUE)
	}

#####
# Prune to the sigmodontinae clade from the ALL (w pastis) posterior

sampled<-read.table("PC12_Cricetidae_Sigmos_allSp-413.txt")
toDrop<-setdiff(cricetid10k[[1]]$tip.label,sampled$V1)

sigmo10k<-vector("list",10000)
for(i in 1:length(sigmo10k)){
	sigmo10k[[i]]<-drop.tip(cricetid10k[[i]],toDrop)
	}

for(i in 1:length(sigmo10k)){
	write.tree(sigmo10k[[i]], "infile.nex.runs1-4_postburn_reScaled-SIGMO-all.trees", append=TRUE)
	}

#####
# Prune to the sigmodontinae clade from the SAMPLED (w genetics) posterior

sampled<-read.table("PC12_Cricetidae_Sigmos_geneticSp-279.txt")
toDrop<-setdiff(cricetid10k_gen[[1]]$tip.label,sampled$V1)

sigmo10k_gen<-vector("list",10000)
for(i in 1:length(sigmo10k_gen)){
	sigmo10k_gen[[i]]<-drop.tip(cricetid10k_gen[[i]],toDrop)
	}

for(i in 1:length(sigmo10k_gen)){
	write.tree(sigmo10k_gen[[i]], "infile.nex.runs1-4_postburn_reScaled-SIGMO-sampled.trees", append=TRUE)
	}

#########

for(i in 1:length(cricetid10k_gen)){
	write.tree(cricetid10k_gen[[i]], "infile.nex.runs1-4_postburn_reScaled-CRICETIDAE-sampled.trees", append=TRUE)
	}
for(i in 1:length(sigmo10k)){
	write.tree(sigmo10k[[i]], "infile.nex.runs1-4_postburn_reScaled-SIGMO-all.trees", append=TRUE)
	}
for(i in 1:length(sigmo10k_gen)){
	write.tree(sigmo10k_gen[[i]], "infile.nex.runs1-4_postburn_reScaled-SIGMO-sampled.trees", append=TRUE)
	}

########
# Still need to RESCALE the SIMPLE output of MCC sigmodontinae tree:

MCCsig<-read.nexus("infile.nex.con_MCC_delNonUnsampSigs-279taxa_SIGMO-simple.tre")

MCCsig2<-rescaleTree(MCCsig,scale=17.6) # this is the inferred age for Sigmos from the MEAN GLOBAL div time between Cricetus and Mus

write.tree(MCCsig2, "infile.nex.con_MCC_delNonUnsampSigs-279taxa_SIGMO-simple_reScaled.tre")

quartz(width=8.5, height=22) 
plot(MCCsig2, cex=0.4, label.offset=0.5, x.lim=c(0,30), y.lim=c(0,270))
node.support(branching.times(MCCsig2), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.4)
node.support(MCCsig2$node.label, mode="dots", digits=2,cutoff=0.90,pos="below", col = c("black","red"), font=2, cex=0.4)
node.support(MCCsig2$node.label, mode="numbers", digits=2,pos="below", col = "grey30", font=2, cex=0.4)
data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-6.5, mgp=c(0,0.2,0))

#######
# AND rescale the simple output of MCC for CRICETIDAE:

MCCcric<-read.nexus("infile.nex.con_MCC_delUnsamp-530taxa_CRICETIDAE-simple.tre")

MCCcric2<-rescaleTree(MCCcric,scale=25.88) # the MEAN GLOBAL div time between Cricetus and Mus

write.tree(MCCcric2, "nfile.nex.con_MCC_delUnsamp-530taxa_CRICETIDAE-simple_reScaled.tre")

quartz(width=8.5, height=22) 
plot(MCCcric2, cex=0.2, label.offset=0.5, x.lim=c(0,30), y.lim=c(0,530))
node.support(branching.times(MCCcric2), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.4)
node.support(MCCcric2$node.label, mode="dots", digits=2,cutoff=0.90,pos="below", col = c("black","red"), font=2, cex=0.4)
node.support(MCCcric2$node.label, mode="numbers", digits=2,pos="below", col = "grey30", font=2, cex=0.4)
data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-12.5, mgp=c(0,0.2,0))




library(ape)

setwd("/mnt/data2/scratch/Nate_Backup/")

#setwd("/MamPhy_BDvarRatesALL_fullPosterior")

#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/11_makingFullPosteriors")

# scan in trees
#trees1=scan("MamPhy_fullPosterior_BDvarRates_17Exp_Last100.trees", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

trees1=scan("MamPhy_fullPosterior_BDvarRates_17Uni_all10k.trees", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

#trees=strsplit(trees1,"[;]")  #trees now indexed #tree1=trees[[1]]  #first tree

#==================
# newick to nexus
#==================

#######
# global -- with imputed-- 5910 species
treesPhylo<-vector("list",length(trees1))
for(i in 1:length(trees1)) {
	treesPhylo[[i]]<-read.tree(text=trees1[[i]])
}
write.nexus(treesPhylo, file="MamPhy_fullPosterior_BDvarRates_17Uni_all10k_nexus.trees")


###
# DROPPING TAXA

ss<-as.vector(read.table("MamPhy_FIN4_1812sp_missing_LIST.txt", header=FALSE))

# For each tree, read in as PHYLO and then prune and write to file (append)
ntrees = length(trees1)
phy<-treesPhylo[[1]]

trees_drop<-vector("list",ntrees)
for(i in 1:ntrees) {
	phy_drop<-drop.tip(treesPhylo[[i]], phy$tip.label[match(ss[,1], phy$tip.label)])
	write.tree(phy_drop, "MamPhy_fullPosterior_BDvarRates_17Uni_all10k_pruned_4098spp.trees", append = TRUE)
	trees_drop[[i]] <- phy_drop
}
write.nexus(trees_drop, file="MamPhy_fullPosterior_BDvarRates_17Uni_all10k_pruned_4098spp_nexus.trees")


q()

n

library(ape)

setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors")

# scan in trees
trees1=scan("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_topoAsZhou_all10k.trees", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

ntrees = length(trees1)

trees_toWrite<-vector("list",ntrees)
for(i in 1:ntrees) {
	phy<-read.tree(text=trees1[[i]])
	trees_toWrite[[i]] <- phy
}

write.nexus(trees_toWrite, file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_topoAsZhou_all10k_nexus.trees")


##
# scan in trees
trees1=scan("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_all10k.trees", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

ntrees = length(trees1)

trees_toWrite<-vector("list",ntrees)
for(i in 1:ntrees) {
	phy<-read.tree(text=trees1[[i]])
	trees_toWrite[[i]] <- phy
}

write.nexus(trees_toWrite, file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_all10k_nexus.trees")



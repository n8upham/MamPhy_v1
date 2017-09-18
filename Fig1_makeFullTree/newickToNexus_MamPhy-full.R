library(ape)

#setwd("~/MamPhy_BDvarRatesALL_fullPosterior")

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/11_makingFullPosteriors")

BDall10k <- read.tree("MamPhy_fullPosterior_BDvarRates_17Exp_Last100.trees")

write.nexus(BDall10k, file="MamPhy_fullPosterior_BDvarRates_17Exp_Last100_nexus.trees")

q()

n
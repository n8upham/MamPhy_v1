library(ape)
library(RPANDA)

####

# TUTORIAL:

data(Phyllostomidae)
data(Phyllostomidae_genera)
res1<-spectR(Phyllostomidae)
pdf(file="spectR_example.pdf")
plot_spectR(Phyllostomidae)
dev.off()
res1$eigengap # = 3 modalities

res2<-BICompare(Phyllostomidae,3)
pdf(file="spectR_example_BIC.pdf")
plot_BICompare(Phyllostomidae,res2) # shows branches colored by given modalities, default to FAN
dev.off()

plot_BICompare <- function (phylo, BICompare) 
{
    if (!inherits(BICompare, "BICompare")) 
        stop("object \"BICompare\" is not of class \"BICompare\"")
    t <- max(BICompare[[2]])
    col_edge <- rainbow(t)[BICompare[[2]][phylo$edge[, 2]]]
    col_tip <- rainbow(t)[BICompare[[2]][1:length(phylo$tip.label)]]
    plot(phylo, edge.color = col_edge, tip.color = col_tip, type = "phylogram", 
        cex = 0.4)
} # changes to normal phylogram

res3<-JSDtree(Phyllostomidae_genera)
pdf(file="spectR_example_JSD-25trees.pdf")
JSDtree_cluster(res3)
dev.off()

######
# with MAMPHY:

setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors")

#mamPhy10k<-read.tree("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_all10k.trees")

trees1=scan("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_all10k.trees", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

phy<-read.tree(text=trees1[[1]])

res1<-spectR(phy)

pdf(file="spectR_mamPhy_FBD_1tree.pdf")
plot_spectR(phy)
dev.off()

res1$eigengap # = how many modalities?

write.table(res1$eigengap,"res_mamPhy_FBD_1tree.txt")

res2<-BICompare(phy,res1$eigengap)

pdf(file="spectR_BIC_mamPhy_FBD_1tree.pdf")
plot_BICompare(phy,res2) # shows branches colored by given modalities, default to FAN
dev.off()

plot_BICompare <- function (phylo, BICompare) 
{
    if (!inherits(BICompare, "BICompare")) 
        stop("object \"BICompare\" is not of class \"BICompare\"")
    t <- max(BICompare[[2]])
    col_edge <- rainbow(t)[BICompare[[2]][phylo$edge[, 2]]]
    col_tip <- rainbow(t)[BICompare[[2]][1:length(phylo$tip.label)]]
    plot(phylo, edge.color = col_edge, tip.color = col_tip, type = "phylogram", 
        cex = 0.4)
} # changes to normal phylogram





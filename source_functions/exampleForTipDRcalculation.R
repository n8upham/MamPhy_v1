# Code for running tip DR calculations
####
# NS Upham, 8 Feb 2018
#####
# load packages
library(ape)

# directory and source
dirname = "/set/dir/name/"
setwd(dirname)

source("DR_functions.R")


# set tree number (do this in loop and them summarize across ~100 trees!)
i<-1


#==================
# Load in 1 tree of 100
#==================
mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_nexus.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 
write.tree(mamPhy,file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""))
tree1=scan(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""), what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 


#==================
# Calculate ES and DR on per-tip basis
#==================
# gives pairwise clade matrix from CAIC function
clade_matrix = readCAIC(tree1)
#
## calculate and write to file
DR = 1/ES_v2(clade_matrix)
ES = ES_v2(clade_matrix)
res = cbind.data.frame(DR,ES)
res1 = res[order(rownames(res)),]
#
write.table(res1, file=paste(bbone,"_sample100_",i,"_DRtips.txt",sep=""))
res1<-read.table(file=paste(bbone,"_sample100_",i,"_DRtips.txt",sep=""), header=TRUE)

# tip DR values to use for this tree's comparisons
res2<-res1[order(rownames(res1)),]

tipDR_i<-res2$DR
names(tipDR_i)<-rownames(res2)

# etc... !





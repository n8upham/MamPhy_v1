library(ape)
library(phytools)
library(phyloch)

#Static code
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/9_patch-PASTIS-ing/PC12_Cricetidae_MB_1re3_33M-fin/")

rescaleTree<-function(tree,scale)
{
tree$edge.length<-tree$edge.length/max(nodeHeights(tree)[,2])*scale
return(tree)
}

tree<-read.nexus("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/9_patch-PASTIS-ing/PC12_Cricetidae_MB_1re3_33M-fin/infile.nex.runs1-4_postburn_lastTree.tre")

tree2<-rescaleTree(tree,scale=25.88) # this is the MEAN GLOBAL div time between Cricetus and Mus

quartz(width=8.5,height=11)
plot(tree, cex=0.1)
axisPhylo()

quartz(width=8.5,height=22)
plot(tree2, cex=0.2, label.offset=0.1)
axisPhylo()

## >>May want to prune out RATTUS first... nope, use that Cricetidae-Muridae split





#Branching times of clades, unscaled
#clade.times <- lapply(clade.trees, branching.times)
clade.times <- branching.times(tree)

#Scale of outgroup to ingroup
#clade.scales <- lapply(clade.times, function(x) {x[2]/x[1]})



#Clade trees unscaled without outgroup
#clade.trees.add <- lapply(clade.trees, drop.tip, 1)
clade.trees.add <- drop.tip(tree, 1)


#Root ages of trees
clade.root.ages <- ages[as.character(root)]

#Tip.ages of trees
clade.tip.ages <- setNames(backbone$edge.length[unlist(lapply(tips,function(x){which.edge(backbone,x)}))],tips)

#Scaled root age of tree
clade.root.scale <- unlist(clade.scales) * clade.root.ages

#Clade trees rescaled, ready to paste
clade.trees.scaled <- list()
for (i in 1:length(clade.trees))
{
clade.trees.scaled[[i]] <- rescaleTree(clade.trees.add[[i]], clade.root.scale[i])
}

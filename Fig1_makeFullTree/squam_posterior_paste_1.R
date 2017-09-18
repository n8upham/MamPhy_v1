library(ape)
library(phytools)

#Static code
setwd("C://Users//Alex/Desktop//squam//Posterior")

rescaleTree<-function(tree,scale)
{
tree$edge.length<-tree$edge.length/max(nodeHeights(tree)[,2])*scale
return(tree)
}

clade <- read.csv("clades.csv")
clades <- as.vector(clade[,1])
ingroup <- as.vector(clade[,2])
outgroup <- as.vector(clade[,3])

#Changing code

#Backbone run, 2500 trees
backbones <- read.nexus("squam_shl_new_Backbone_Myr.run1.t")

#Clade trees
for (i in 1:length(clades))
    {
    assign(paste(clades[i], ".trees", sep=""), read.nexus(paste(clades[i], ".nex.run1.t", sep="")))
    }
    

#Loop to generate 2500 trees 

for (z in 0:2499)
{   

#Backbone
backbone <- backbones[[2500 - z]]
ages <- branching.times(backbone) 
nodes <- mrca(backbone)

tips <- vector()
root <- vector()
for(i in 1:length(clades)) 
{
tips[i] <- match(ingroup[i],backbone$tip.label)
root[i] <- nodes[ingroup[i],outgroup[i]]
}

#Clades
clade.trees <- list()
for (i in 1:length(clades))
{
#Clade trees
clade.trees[[i]] <- get(paste(clades[i],".trees",sep=""))[[length(get(paste(clades[i],".trees",sep=""))) - z]]
}

#Branching times of clades, unscaled
clade.times <- lapply(clade.trees, branching.times)

#Scale of outgroup to ingroup
clade.scales <- lapply(clade.times, function(x) {x[2]/x[1]})

#Clade trees unscaled without outgroup
clade.trees.add <- lapply(clade.trees, drop.tip, 1)

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


#Loop to bind trees
paste_tree <- backbone 

for (i in 1:length(clades))
{

#Tip edge
clade.tip <- match(ingroup[i],paste_tree$tip.label)
tip.edge <- which.edge(paste_tree, clade.tip)

if (clade.root.scale[i] < clade.tip.ages[i])
{
#Shorten branch
paste_tree$edge.length[tip.edge] <- as.vector(clade.tip.ages[i] - clade.root.scale[i])

#Bind tree
paste_tree <- bind.tree(paste_tree, clade.trees.scaled[[i]], clade.tip, 0)

} else {

#Make tree fit tip age
re.tree <- rescaleTree(clade.trees.scaled[[i]], as.vector(paste_tree$edge.length[tip.edge] * 0.999999))

#Shorten branch
paste_tree$edge.length[tip.edge] <- as.vector(paste_tree$edge.length[tip.edge] * 0.000001)

#Bind tree
paste_tree <- bind.tree(paste_tree, re.tree, clade.tip, 0)
}

}

write.tree(paste_tree, "squam_shl_new_Posterior_9755.run1.trees", append = TRUE)

}

q()

n

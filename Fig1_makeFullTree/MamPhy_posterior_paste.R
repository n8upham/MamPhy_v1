library(ape)
library(phytools)
library(phyloch)

#Static code
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/11_makingFullPosteriors")
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")

rescaleTree<-function(tree,scale)
{
tree$edge.length<-tree$edge.length/max(nodeHeights(tree)[,2])*scale
return(tree)
}

#clade <- read.csv("clades.csv")
#clades <- as.vector(clade[,1])
#ingroup <- as.vector(clade[,2])
#outgroup <- as.vector(clade[,3])

clade<-read.table("clade_inRep_outgroup_BDvr_Fixed-n-Free.txt")
clades <- as.vector(clade[,1])
ingroup <- as.vector(clade[,2])
outgroup <- as.vector(clade[,3])


#Changing code

#Backbone run, FIRST with 100 trees
backbones <- read.nexus("postburnLast100_Bbone_ND_60taxa_17calExp_BDallDNA.trees")

## PLOT one of the backbones:
tree<-ladderize(backbones[[1]])
#tree<-ladderize(bboneDrop[[1]])

quartz(width=8.5, height=11) 
plot(tree, cex=0.5, label.offset=0.2)#, x.lim=c(0,15))#, y.lim=c(0,530))
node.support(branching.times(tree), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.6)
data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-12.5, mgp=c(0,0.2,0))

nodelabels(cex=0.4)
tiplabels(cex=0.4)

# Drop taxa to just the INGROUP-OUTGROUP reps ...
toKeep <- c("_Anolis_carolinensis", "Anomalurus_beecrofti__ANOMALURIDAE__RODENTIA", "Bison_bison__BOVIDAE__CETARTIODACTYLA", "Callithrix_jacchus__CALLITRICHIDAE__PRIMATES", "Calomyscus_baluchi__CALOMYSCIDAE__RODENTIA", "Castor_canadensis__CASTORIDAE__RODENTIA", "Cricetomys_gambianus__NESOMYIDAE__RODENTIA", "Cricetulus_barabensis__CRICETIDAE__RODENTIA", "Dasypus_novemcinctus__DASYPODIDAE__CINGULATA", "Eptesicus_fuscus__VESPERTILIONIDAE__CHIROPTERA", "Equus_caballus__EQUIDAE__PERISSODACTYLA", "Erethizon_dorsatum__ERETHIZONTIDAE__RODENTIA", "Erinaceus_europaeus__ERINACEIDAE__EULIPOTYPHLA", "Felis_catus__FELIDAE__CARNIVORA", "Galeopterus_variegatus__CYNOCEPHALIDAE__DERMOPTERA", "Ictidomys_tridecemlineatus__SCIURIDAE__RODENTIA", "Jaculus_jaculus__DIPODIDAE__RODENTIA", "Manis_pentadactyla__MANIDAE__PHOLIDOTA", "Ornithorhynchus_anatinus__ORNITHORHYNCHIDAE__MONOTREMATA", "Oryctolagus_cuniculus__LEPORIDAE__LAGOMORPHA", "Pteronotus_parnellii__MORMOOPIDAE__CHIROPTERA", "Pteropus_alecto__PTEROPODIDAE__CHIROPTERA", "Rattus_norvegicus__MURIDAE__RODENTIA", "Saccopteryx_bilineata__EMBALLONURIDAE__CHIROPTERA", "Spalax_ehrenbergi__SPALACIDAE__RODENTIA", "Trichechus_manatus__TRICHECHIDAE__SIRENIA", "Tupaia_belangeri__TUPAIIDAE__SCANDENTIA", "Typhlomys_cinereus__PLATACANTHOMYIDAE__RODENTIA", "Vombatus_ursinus__VOMBATIDAE__DIPROTODONTIA")
toDrop <- setdiff(backbones[[1]]$tip.label,toKeep)

# drop non-IN-OUT reps, first for 100 trees
bboneDrop<-vector("list",100)
for (i in 1:length(backbones)){
	bboneDrop[[i]]<-drop.tip(backbones[[i]],toDrop)
	}

#Clade trees
for (i in 1:length(clades))
    {
    assign(paste(clades[i], ".trees", sep=""), read.nexus(paste("postburnLast100_BDvr_FIXED_",clades[i], ".trees", sep=""))) 
    }


#Prune out extra taxa from PC25_Dermoptera and PC28_Platacanthomyidae
PC25_DermopteraOK<-vector("list",length(PC25_Dermoptera.trees))
PC28_PlatacanthomyidaeOK<-vector("list",length(PC28_Platacanthomyidae.trees))

for (i in 1:length(PC25_Dermoptera.trees)){
	PC25_DermopteraOK[[i]]<-drop.tip(PC25_Dermoptera.trees[[i]],c("Callithrix_jacchus_CALLITRICHIDAE_PRIMATES","Gorilla_gorilla_HOMINIDAE_PRIMATES"))
	PC28_PlatacanthomyidaeOK[[i]]<-drop.tip(PC28_Platacanthomyidae.trees[[i]],c("Rattus_norvegicus_MURIDAE_RODENTIA","Spalax_ehrenbergi_SPALACIDAE_RODENTIA"))
}
PC25_Dermoptera.trees<-vector("list",length(PC25_DermopteraOK))
PC28_Platacanthomyidae.trees<-vector("list",length(PC28_PlatacanthomyidaeOK))
for (i in 1:length(PC25_Dermoptera.trees)){
	PC25_Dermoptera.trees[[i]]<-PC25_DermopteraOK[[i]]
	PC28_Platacanthomyidae.trees[[i]]<-PC28_PlatacanthomyidaeOK[[i]]
}



######
#Loop to generate 100 trees 

for (z in 0:99)
{   

#Backbone
backbone <- bboneDrop[[100 - z]]
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
pasted_fulltree <- backbone 

for (i in 1:length(clades))
{

#Tip edge
clade.tip <- match(ingroup[i],pasted_fulltree$tip.label)
tip.edge <- which.edge(pasted_fulltree, clade.tip)

if (clade.root.scale[i] < clade.tip.ages[i])
{
#Shorten branch
pasted_fulltree$edge.length[tip.edge] <- as.vector(clade.tip.ages[i] - clade.root.scale[i])

#Bind tree
pasted_fulltree <- bind.tree(pasted_fulltree, clade.trees.scaled[[i]], clade.tip, 0)

} else {

#Make tree fit tip age
re.tree <- rescaleTree(clade.trees.scaled[[i]], as.vector(pasted_fulltree$edge.length[tip.edge] * 0.999999))

#Shorten branch
pasted_fulltree$edge.length[tip.edge] <- as.vector(pasted_fulltree$edge.length[tip.edge] * 0.000001)

#Bind tree
pasted_fulltree <- bind.tree(pasted_fulltree, re.tree, clade.tip, 0)
}

}

write.tree(pasted_fulltree, "MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_Last100.trees", append = TRUE)

}


######
# Plot a FULL tree, so that is read-able.

mamPhy<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_commonAnc2.tre")
tree2 <- drop.tip(mamPhy,"_Anolis_carolinensis")
tree <- ladderize(tree2, right=TRUE)

pdf(file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_commonAnc_5911taxa_linear.pdf", width=8.5, height=130, onefile=TRUE)

#quartz(width=8.5, height=130)
plot(tree, cex=0.15, label.offset=0.4, y.lim=c(0,5690))
node.support(branching.times(tree), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.15)
data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-18.5, mgp=c(0,0.2,0))

dev.off()


pdf(file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_commonAnc_5911taxa_radial.pdf", width=25, height=25, onefile=TRUE)

plot(tree, cex=0.05, label.offset=0.4, type="fan", edge.width=0.5, show.tip.label=TRUE)#, x.lim=c(0,15))#, y.lim=c(0,530))
node.support(branching.times(tree), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.2)

dev.off()

##########
###
# Plot small phylo of MCC
###
# For MULTILEVEL figure of my semnar talk at OU
#########
library(ape)
library(phyloch)

tree <- read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target.tre")
tree2 <- drop.tip(tree,"_Anolis_carolinensis")
tree3 <- ladderize(tree2, right=TRUE)
#quartz(height=10,width=5)

pdf(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911taxa_linear_noTips_wide80_greenTipLabs.pdf", width=5, height=80, onefile=TRUE)
plot(tree3, edge.width=0.5, cex=0.05, edge.color="dark green", show.tip.label=TRUE)
dev.off()

png(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911taxa_linear_noTips_wide10_greenTipLabs.png", width=10, height=7.5, units="in", res=200, bg="transparent")
plot(tree3, edge.width=0.5, cex=0.05, edge.color="dark green", show.tip.label=FALSE, direction="downwards")
axisPhylo(side=2, cex=2)
dev.off()

cladesDR<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_DR-SUMMARY-expanded_NDexp.txt",sep="")
colnames(cladesDR)<-c("tiplabel","gen","fam","ord","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")
head(cladesDR)

# prune to OCTODONTIDAE
octos<-cladesDR[which(cladesDR$fam=="OCTODONTIDAE"),"tiplabel"]
nonOctos<-setdiff(tree3$tip.label,octos)

octoTree<-drop.tip(tree3,nonOctos)
rootnode<-getMRCA(octoTree,as.vector(octos))
octoTree<-paintSubTree(octoTree,node=rootnode,state="1")

col<-"red"; names(col)<-"1"

png(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_subtoOctodontidae_Red_left.png", width=3, height=3, units="in", res=200, bg="transparent")
plotSimmap(octoTree, col, lwd=3, ftype="off",mar=c(5.1,1.1,1.1,1.1),direction="leftwards")
dev.off()


png(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_subtoOctodontidae_Red.png", width=7.5, height=10, units="in", res=200, bg="transparent")
plotSimmap(octoTree, col, lwd=2, ftype="off",mar=c(5.1,1.1,1.1,1.1))
axisPhylo(side=1, cex=2)
dev.off()

png(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_subtoOctodontidae_Red_sm.png", width=3, height=10, units="in", res=200, bg="transparent")
plotSimmap(octoTree, col, lwd=2, ftype="off",mar=c(5.1,1.1,1.1,1.1))
axisPhylo(side=1, cex=2)
dev.off()

png(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_subtoOctodontidae_Red_fan.png", width=5, height=5, units="in", res=200, bg="transparent")
plotSimmap(octoTree, col, type="fan", lwd=2, ftype="off",mar=c(5.1,1.1,1.1,1.1))
dev.off()

# prune to CTENOHYSTRICA
ctenos<-cladesDR[which(cladesDR$PC=="PC14_GuineaPigRelated"),"tiplabel"]
nonCtenos<-setdiff(tree3$tip.label,ctenos)

ctenoTree<-drop.tip(tree3,nonCtenos)
nonOctoCtenos<-setdiff(ctenoTree$tip.label,octos)

octoNode_inCtenos<-getMRCA(ctenoTree, as.vector(octos))

ctenoTree2<-paintSubTree(ctenoTree,node=octoNode_inCtenos,state="2",stem=FALSE)

cols<-c("dodgerblue4","red"); names(cols)<-1:2

png(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_subtoCtenoBlue_octoRed_tips_fan.png", width=5, height=5, units="in", res=200, bg="transparent")
plotSimmap(ctenoTree2, cols, type="fan",lwd=2, ftype="off",mar=c(5.1,1.1,1.1,1.1))
dev.off()

png(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_subtoCtenoBlue_octoRed_tips_sm.png", width=3, height=10, units="in", res=200, bg="transparent")
plotSimmap(ctenoTree2, cols, lwd=2, ftype="off",mar=c(5.1,1.1,1.1,1.1))
axisPhylo(side=1, cex=2)
dev.off()

png(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_subtoCtenoBlue_octoRed_tips.png", width=7.5, height=10, units="in", res=200, bg="transparent")
plotSimmap(ctenoTree2, cols, lwd=2, ftype="off",mar=c(5.1,1.1,1.1,1.1))
axisPhylo(side=1, cex=2)
dev.off()


# all of MAMMALIA, with colors
ctenoNode <- getMRCA(tree3, as.vector(ctenos))
octoNode <- getMRCA(tree3, as.vector(octos))

tree4<-paintSubTree(tree3,node=ctenoNode,state="2",stem=FALSE)
tree4<-paintSubTree(tree4,node=octoNode,state="3",stem=FALSE)

cols<-c("black","dodgerblue4","red"); names(cols)<-1:3

png(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_allMammBlk_octoRedCtenoBlue_fan.png", width=5, height=5, units="in", res=200, bg="transparent")
plotSimmap(tree4, cols, lwd=0.5, type="fan", ftype="off",mar=c(5.1,1.1,1.1,1.1))
dev.off()

png(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_allMammBlk_octoRedCtenoBlue_sm.png", width=3, height=10, units="in", res=200, bg="transparent")
plotSimmap(tree4, cols, lwd=0.5, ftype="off",mar=c(5.1,1.1,1.1,1.1))
axisPhylo(side=1, cex=2)
dev.off()

png(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_allMammBlk_octoRedCtenoBlue_placental.png", width=7.5, height=10, units="in", res=200, bg="transparent")
plotSimmap(tree4, cols, lwd=0.5, ftype="off",mar=c(5.1,1.1,1.1,1.1))
axisPhylo(side=1, cex=2)
dev.off()


# to PLACENTALIA
marMono<-cladesDR[which(cladesDR$PC=="PC1_Marsupials" | cladesDR$PC=="PC23_Monotremata"),"tiplabel"]

placentalTree<-drop.tip(tree3,as.vector(marMono))
ctenoNode2 <- getMRCA(placentalTree, as.vector(ctenos))
octoNode2 <- getMRCA(placentalTree, as.vector(octos))

placentalTree2<-paintSubTree(placentalTree,node=ctenoNode2,state="2",stem=FALSE)
placentalTree2<-paintSubTree(placentalTree2,node=octoNode2,state="3",stem=FALSE)

cols<-c("black","dodgerblue4","red"); names(cols)<-1:3

png(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_allPlacentalsBlk_octoRedCtenoBlue_sm.png", width=3, height=10, units="in", res=200, bg="transparent")
plotSimmap(placentalTree2, cols, lwd=0.5, ftype="off",mar=c(5.1,1.1,1.1,1.1))
axisPhylo(side=1, cex=2)
dev.off()

png(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_allPlacentalsBlk_octoRedCtenoBlue_placental.png", width=7.5, height=10, units="in", res=200, bg="transparent")
plotSimmap(placentalTree2, cols, lwd=0.5, ftype="off",mar=c(5.1,1.1,1.1,1.1))
axisPhylo(side=1, cex=2)
dev.off()




#> tree<-paintSubTree(tree,node=22,state="2")
#> tree<-paintSubTree(tree,node=26,state="3")
#> # now let's plot using plotSimmap to ensure
#> # that the correct branches were painted
#> cols<-c("black","blue","red"); names(cols)<-1:3
#> plotSimmap(tree,cols,pts=F,lwd=3,node.numbers=T)





#####
# Plot the MCC of the FBD and NDexp TARGET=MCC topology and node height trees... in a nice way:
# this MCC made with TreeAnnotator.
tree <- read.beast("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target.tre")
tree2 <- drop.tip2(tree,"_Anolis_carolinensis")
tree3 <- ladderize(tree2, right=TRUE)


missing<-read.table("MamPhy_FIN4_1813sp_missing_LIST.txt", header=FALSE)

tipColors <- rep("black",5912)
tipColors[match(missing$V1,tree$tip.label)] <- "red" 

pdf(file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_5911taxa_linear_ladderPP_dotsTargetAges.pdf", width=11, height=260, onefile=TRUE)
pdf(file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_5912taxa_linear_nodelabels.pdf", width=11, height=260, onefile=TRUE)

#quartz(width=8.5, height=130)
plot(mamPhy, cex=0.25, tip.color=tipColors, label.offset=0.4, y.lim=c(0,5690), x.lim=c(-28.47793, 337.00159))

nodelabels(cex=0.2)
#HPDbars(tree3, label="height_95%_HPD", broken=T, lwd=1.5, col=hsv(0.65,1,1,alpha=0.7))
#node.support(tree3$posterior, mode="dots", col = "red", cex=0.18)
#node.support(tree3$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.2)
#node.support(branching.times(tree3), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.15)
#node.support(tree3$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.15)

data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-18.5, mgp=c(0,0.2,0))

dev.off()

## Now with MEAN ages annotated:
pdf(file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_5911taxa_linear_ladderPP_dotsMeanAges.pdf", width=11, height=260, onefile=TRUE)

#quartz(width=8.5, height=130)
plot(tree3, cex=0.15, tip.color=tipColors, label.offset=0.4, y.lim=c(0,5690), x.lim=c(-28.47793, 337.00159))

HPDbars(tree3, label="height_95%_HPD", broken=T, lwd=1.5, col=hsv(0.65,1,1,alpha=0.7))
node.support(tree3$posterior, mode="dots", col = "red", cex=0.18)
node.support(tree3$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.2)
node.support(tree3$height, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.15)
node.support(tree3$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.15)

data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-18.5, mgp=c(0,0.2,0))

dev.off()



#####
# Plot the MCC of the FBD and NDexp COMMONANC trees... in a nice way:
# this MCC made with TreeAnnotator.

tree <- read.beast("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_commonAnc2.tre")

tree2 <- drop.tip2(tree,"_Anolis_carolinensis")
tree3 <- ladderize(tree2, right=TRUE)


missing<-read.table("MamPhy_FIN4_1813sp_missing_LIST.txt", header=FALSE)

tipColors <- rep("black",5912)
tipColors[match(missing$V1,tree$tip.label)] <- "red" 

pdf(file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_commonAnc2_5911taxa_linear_NOladderPP_dots.pdf", width=11, height=260, onefile=TRUE)

#quartz(width=8.5, height=130)
plot(tree, cex=0.15, tip.color=tipColors, label.offset=0.4, y.lim=c(0,5690), x.lim=c(-28.47793, 337.00159))

HPDbars(tree, label="height_95%_HPD", broken=T, lwd=1.5, col=hsv(0.65,1,1,alpha=0.7))
node.support(tree$posterior, mode="dots", col = "red", cex=0.18)
node.support(tree$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.2)
node.support(branching.times(tree), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.15)
node.support(tree$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.15)

data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-18.5, mgp=c(0,0.2,0))

dev.off()

# >> BUT something is wrong with the CommonAnc node heights... so don't use them.




#######
# Same MCC plotting of last100, but with PRUNED 4098 taxon tree
##
tree <- read.beast("MamPhy_fullPosterior_BDvarRates_17Exp_last100_pruned_4098spp_nexus_MCC_treeAnOut.tre")

tree2 <- drop.tip2(tree,"_Anolis_carolinensis")
tree3 <- ladderize(tree2, right=TRUE)

#node.trans(tree2, tree3, index = FALSE) # did not need this

pdf(file="MamPhy_fullPosterior_BDvarRates_Last100_MCC-tree_4098taxa_linear.pdf", width=8.5, height=100, onefile=TRUE)

#quartz(width=8.5, height=100)
plot(tree2, cex=0.15, label.offset=0.4, y.lim=c(0,3950))

node.support(branching.times(tree2), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.15)
#node.support(tree3$posterior, mode="numbers", digits=2,pos="above", col = "black", font=2, cex=0.15)
HPDbars(tree2, label="height_95%_HPD", broken=T, lwd=1.5, col=hsv(0.65,1,1,alpha=0.7))

data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-18.5, mgp=c(0,0.2,0))

dev.off()

####
# NOW: do the MCC of the all10k EXP pruned 4099taxa tree.
##
#tree <- read.beast("MamPhy_fullPosterior_BDvarRates_17Exp_all10k_pruned_4098spp_nexus_MCC_treeAnOut.tre")

tree <- read.beast("MamPhy_fullPosterior_BDvarRates_17Exp_all10k_pruned_4098spp_nexus_MCC-targetHeights.tre")
tree2 <- drop.tip2(tree,"_Anolis_carolinensis")
tree3 <- ladderize(tree2, right=TRUE)

pdf(file="MamPhy_fullPosterior_BDvarRates_17Exp_all10k_MCC-targetHeights_4098taxa_linear_agesPP.pdf", width=8.5, height=122, onefile=TRUE)

#quartz(width=8.5, height=122)
plot(tree3, cex=0.25, label.offset=0.4, y.lim=c(130,3950))

HPDbars(tree3, label="height_95%_HPD", broken=T, lwd=1.5, col=hsv(0.65,1,1,alpha=0.7))
#node.support(tree3$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.15)
#node.support(tree3$posterior, mode="numbers", digits=2,pos="below", cutoff=0.95, col = "black", font=2, cex=0.19)
node.support(tree3$posterior, mode="dots", col = "red", cex=0.2)
node.support(tree3$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.26)
node.support(tree3$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.25)
node.support(tree3$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.15)

data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-18.5, mgp=c(0,0.2,0))

dev.off()

####
#  MCC of the all10k EXP fully sampled 5910 taxa tree. -- HAS NEGATIVE BRANCH LENGTHS
##

#tree_prune <- read.beast("MamPhy_fullPosterior_BDvarRates_17Exp_all10k_pruned_4098spp_nexus_MCC_treeAnOut.tre")
tree <- read.beast("MamPhy_fullPosterior_BDvarRates_17Exp_all10k_nexus_MCC_treeAnOut.tre")

tree2 <- drop.tip2(tree,"_Anolis_carolinensis")
tree3 <- ladderize(tree2, right=TRUE)
#tree4 <- untangle(tree3, method="reorder")

missing<-read.table("MamPhy_FIN4_1812sp_missing_LIST.txt", header=FALSE)

tipColors <- rep("black",5910)
tipColors[match(missing$V1,tree3$tip.label)] <- "red" 


#node.trans(tree2, tree3, index = FALSE) # did not need this

pdf(file="MamPhy_fullPosterior_BDvarRates_17Exp_all10k_MCC-tree_5910taxa_linear_agesPP_unsampRed.pdf", width=8.5, height=144, onefile=TRUE)

#quartz(width=8.5, height=144)
plot(tree3, cex=0.2, tip.color=tipColors, label.offset=0.4, y.lim=c(130,5700))

HPDbars(tree3, label="height_95%_HPD", broken=T, lwd=1.5, col=hsv(0.65,1,1,alpha=0.7))
#node.support(tree3$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.15)
#node.support(tree3$posterior, mode="numbers", digits=2,pos="below", cutoff=0.95, col = "black", font=2, cex=0.19)
node.support(tree3$posterior, mode="dots", col = "red", cex=0.2)
node.support(tree3$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.26)
node.support(branching.times(tree3), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.25)
node.support(tree3$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.15)


data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-18.5, mgp=c(0,0.2,0))

dev.off()

####
#  MCC--target node heights-- of the all10k EXP fully sampled 5910 taxa tree. -- OK
##
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")
#tree_prune <- read.beast("MamPhy_fullPosterior_BDvarRates_17Exp_all10k_pruned_4098spp_nexus_MCC_treeAnOut.tre")
#tree <- read.beast("MamPhy_fullPosterior_BDvarRates_17Exp_all10k_nexus_MCC-targetHeights_treeAnOut.tre")
tree <- read.beast("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target.tre")
tree <- read.beast("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target.tre")


tree2 <- drop.tip2(tree,"_Anolis_carolinensis")
tree3 <- ladderize(tree2, right=TRUE)
#tree4 <- untangle(tree3, method="reorder")

missing<-read.table("MamPhy_FIN4_1813sp_missing_LIST.txt", header=FALSE)

tipColors <- rep("black",5911)
tipColors[match(missing$V1,tree3$tip.label)] <- "red" 


#node.trans(tree2, tree3, index = FALSE) # did not need this

pdf(file="MamPhy_fullPosterior_BDvrFIXED_FBD_MCC_target_5911taxa_linear_agesPP.pdf", width=8.5, height=144, onefile=TRUE)
pdf(file="MamPhy_fullPosterior_BDvrFIXED_NDexp_MCC_target_5911taxa_linear_agesPP.pdf", width=8.5, height=144, onefile=TRUE)

#quartz(width=8.5, height=144)
plot(tree3, cex=0.2, tip.color=tipColors, label.offset=0.4, y.lim=c(130,5700))

HPDbars(tree3, label="height_95%_HPD", broken=T, lwd=1.5, col=hsv(0.65,1,1,alpha=0.7))
#node.support(tree3$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.15)
#node.support(tree3$posterior, mode="numbers", digits=2,pos="below", cutoff=0.95, col = "black", font=2, cex=0.19)
node.support(tree3$posterior, mode="dots", col = "red", cex=0.2)
node.support(tree3$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.26)
node.support(tree3$height, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.2)
node.support(tree3$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.15)


data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-18.5, mgp=c(0,0.2,0))

dev.off()





####
#######
# LTT plots of the entire trees...

# load up the NEXUS posteriors. >> first using 100 trees...
library(ape)
library(phytools)
library(phyloch)

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")

mamExp <- read.nexus("MamPhy_fullPosterior_BDvarRates_17Exp_first100_nexus.trees")

mamUni <- read.nexus("MamPhy_fullPosterior_BDvarRates_17Uni_first100_nexus.trees")

pdf(file="LTT_MamPhy_BDvarRates_Exp-v-Uni_first100.pdf", width=8, height=8, onefile=TRUE)

#quartz(width=8,height=8)
ltt.plot(mamExp[[1]], log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:100) {
	ltt.lines(mamExp[[i]], col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}
for (i in 1:100) {
	ltt.lines(mamUni[[i]], col=hsv(0.9708995,0.2470588,1,alpha=0.5), lty = 2, lwd=1)
}

title(main="MamPhy, BD varRates, 5910 spp, first100 - Exp (blue) vs Uni (pink)")

dev.off()


#col2rgb("pink")
#rgb2hsv(255,192,203)
       [,1]
h 0.9708995
s 0.2470588
v 1.0000000

#########
#####
# DensityTree plots of the PCs pre-re-scaling -- so that can visualize the amount of variation in the posteriors...

# possible code from here: http://blog.phytools.org/2016/01/alternative-density-tree-plotting.html

#Clade trees
clades2<-data.frame(rep(NA,length(clades)))
for (i in 1:length(clades)){
	clades2[i,]<-(paste(clades[i], ".trees", sep=""))
}
names(clades2)<-"PC"

clades3<-vector("list",length(clades))
clades4<-vector("list",length(clades))
for (i in 1:length(clades))
    {
    clades3[[i]]<-assign(paste(clades[i], ".trees", sep=""), read.nexus(paste("postburnLast100_BDvr_FIXED_",clades[i], ".trees", sep=""))) 
    for (j in 1:length(clades3[[i]])){
    	clades4[[i]][[j]]<-ladderize(clades3[[i]][[j]])
    	}
    }
 

#Plot density trees of last 100 -- FIXED

for(i in 1:length(clades2[,1])){
    pdf(file=paste("densityTrees_FIXED_last100_inner_", clades2[i,1], ".pdf", sep=""), width=8.5, height=(length(clades4[[i]][[1]]$tip.label)/10), onefile=FALSE)
    #quartz(width=8.5, height=(length(clades4[[i]][[1]]$tip.label)/10))
    densityTree(clades4[[i]], colors="blue", alpha=NULL, method="plotTree", fix.depth=FALSE, use.edge.length=FALSE, compute.consensus=FALSE, fsize=0.5, ftype="i", nodes="inner")
    add.scale.bar(0,10,cex=0.7,pos=4)
    text(x=0,y=2,labels=paste(clades[i],": BD, DNA spp fixed + unsampled spp imputed",sep=""), pos=4, cex=0.8, font=2)
    dev.off()
	}

for(i in 1:length(clades2[,1])){
    pdf(file=paste("densityTrees_FIXED_last100_nodeHeights_inner_", clades2[i,1], ".pdf", sep=""), width=8.5, height=(length(clades4[[i]][[1]]$tip.label)/10), onefile=FALSE)
    #quartz(width=8.5, height=(length(clades4[[i]][[1]]$tip.label)/10))
    densityTree(clades4[[i]], colors="blue", alpha=NULL, method="plotTree", fix.depth=FALSE, use.edge.length=TRUE, compute.consensus=FALSE, fsize=0.5, ftype="i", nodes="inner")
    add.scale.bar(0,10,cex=0.7,pos=4)
    text(x=0,y=2,labels=paste(clades[i],": BD, DNA spp fixed + unsampled spp imputed",sep=""), pos=4, cex=0.8, font=2)
    dev.off()
	}

###
# And for the FREE patches:

clades5<-vector("list",length(clades))
clades6<-vector("list",length(clades))
for (i in 1:length(clades))
    {
    clades5[[i]]<-assign(paste(clades[i], ".trees", sep=""), read.nexus(paste("postburnLast100_BDvr_FREE_",clades[i], ".trees", sep=""))) 
    for (j in 1:length(clades5[[i]])){
    	clades6[[i]][[j]]<-ladderize(clades5[[i]][[j]])
    	}
    }
 
for(i in 1:length(clades2[,1])){
    pdf(file=paste("densityTrees_FREE_last100_", clades2[i,1], ".pdf", sep=""), width=8.5, height=(length(clades6[[i]][[1]]$tip.label)/10), onefile=FALSE)
    #quartz(width=8.5, height=(length(clades4[[i]][[1]]$tip.label)/10))
    densityTree(clades6[[i]], colors="blue", alpha=NULL, method="plotTree", fix.depth=FALSE, use.edge.length=FALSE, compute.consensus=FALSE, fsize=0.5, ftype="i") #nodes="inner")
    add.scale.bar(0,10,cex=0.7,pos=4)
    text(x=0,y=2,labels=paste(clades[i],": BD, DNA spp only (unconstrained)",sep=""), pos=4, cex=0.8, font=2)
    dev.off()
	}

for(i in 1:length(clades2[,1])){
    pdf(file=paste("densityTrees_FREE_last100_nodeHeights_inner_", clades2[i,1], ".pdf", sep=""), width=8.5, height=(length(clades6[[i]][[1]]$tip.label)/10), onefile=FALSE)
    #quartz(width=8.5, height=(length(clades4[[i]][[1]]$tip.label)/10))
    densityTree(clades6[[i]], colors="blue", alpha=NULL, method="plotTree", fix.depth=FALSE, use.edge.length=TRUE, compute.consensus=FALSE, fsize=0.5, ftype="i", nodes="inner")
    add.scale.bar(0,10,cex=0.7,pos=4)
    text(x=0,y=2,labels=paste(clades[i],": BD, DNA spp only (unconstrained)",sep=""), pos=4, cex=0.8, font=2)
    dev.off()
	}

for(i in 1:length(clades2[,1])){
    pdf(file=paste("densityTrees_FREE_last100_inner_", clades2[i,1], ".pdf", sep=""), width=8.5, height=(length(clades6[[i]][[1]]$tip.label)/10), onefile=FALSE)
    #quartz(width=8.5, height=(length(clades4[[i]][[1]]$tip.label)/10))
    densityTree(clades6[[i]], colors="blue", alpha=NULL, method="plotTree", fix.depth=FALSE, use.edge.length=FALSE, compute.consensus=FALSE, fsize=0.5, ftype="i", nodes="inner")
    add.scale.bar(0,10,cex=0.7,pos=4)
    text(x=0,y=2,labels=paste(clades[i],": BD, DNA spp only (unconstrained)",sep=""), pos=4, cex=0.8, font=2)
    dev.off()
	}

######
## NOW plot just all the MCC's of the PC's and EXAMINE THEM for influence of outgroup coming in-- wuite possible some of these don't need to be redone! Like MARSUPIALS for example.

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/11_makingFullPosteriors/PC_MCCs_FREE_mean")

cladesMCC1<-vector("list",length(clades))
cladesMCC<-vector("list",length(clades))
for (i in 1:length(clades)){
    cladesMCC1[[i]]<-assign(paste(clades[i], ".MCC", sep=""), read.beast(paste("postburn10k_BDvr_FREE_MCC_",clades[i], ".ok.tre", sep=""))) 
    cladesMCC[[i]]<-ladderize(cladesMCC1[[i]])
    }


for(i in 1:length(cladesMCC)){
    pdf(file=paste("postburn10k_BDvr_FREE_MCC_", clades[i], "_plot",".pdf", sep=""), width=8.5, height=((length(cladesMCC[[i]]$tip.label)/10)+2), onefile=FALSE)
    
    #quartz(width=8.5, height=(length(cladesMCC[[i]]$tip.label)/10+2))
	plot(cladesMCC[[i]], cex=0.5, label.offset=0.01)
	node.support(cladesMCC[[i]]$posterior, mode="dots", col = "red", cex=0.4)
	node.support(cladesMCC[[i]]$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.41)
	node.support(cladesMCC[[i]]$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.5)
    add.scale.bar(0,10,cex=0.7,pos=4)
    text(x=0,y=2,labels=paste(clades[i],": BD, DNA spp only (unconstrained)",sep=""), pos=4, cex=0.8, font=2)
    
    dev.off()
	}

####
## AND the SAME for the FIXED MCCs... why not. Examine these to be sure they fit your expectations of the constraints, for example.

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/11_makingFullPosteriors")

clade<-read.table("clade_inRep_outgroup_BDvr_Fixed-n-Free.txt")
clades <- as.vector(clade[,1])

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/11_makingFullPosteriors/PC_MCCs_fixed_mean")

cladesMCC1<-vector("list",length(clades))
cladesMCC<-vector("list",length(clades))
for (i in 1:length(clades)){
    cladesMCC1[[i]]<-assign(paste(clades[i], ".MCC", sep=""), read.beast(paste("postburn10k_BDvr_FIXED_MCC_",clades[i], ".ok.tre", sep=""))) 
    cladesMCC[[i]]<-ladderize(cladesMCC1[[i]])
    }


for(i in 1:length(cladesMCC)){
    pdf(file=paste("postburn10k_BDvr_FIXED_MCC_", clades[i], "_plot",".pdf", sep=""), width=8.5, height=((length(cladesMCC[[i]]$tip.label)/10)+2), onefile=FALSE)
    
    #quartz(width=8.5, height=(length(cladesMCC[[i]]$tip.label)/10+2))
	plot(cladesMCC[[i]], cex=0.5, label.offset=0.01)
	node.support(cladesMCC[[i]]$posterior, mode="dots", col = "red", cex=0.4)
	node.support(cladesMCC[[i]]$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.41)
	node.support(cladesMCC[[i]]$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.5)
    add.scale.bar(0,10,cex=0.7,pos=4)
    text(x=0,y=2,labels=paste(clades[i],": BD, ALL spp, incl imputed (no DNA)",sep=""), pos=4, cex=0.8, font=2)
    
    dev.off()
	}







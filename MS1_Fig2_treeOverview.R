#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy MS1 -- Upham et al. 2019 -- PLOS Biology
###
# Figure 2 - overview of methods for new phylogeny of Mammalia
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# part a == by hand.
# ===
# ===
# part b: plotting the simplified backbone phylogenies for comparison.
# ===
	#install.packages(repos=NULL, pkgs="/Users/nate/Desktop/Phylogen-programs/phyloch_1.5-4.tar.gz")
	library(ape); library(phyloch); library(phytools)
	library(paleotree); library(strap); #library(OutbreakTools)

setwd("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/11_makingFullPosteriors")

#FBD0<-read.beast("postburn10k_Bbone_FBD_60taxa_BDallDNA.inMYR_MCC_treeAnOut.tre")
FBD0<-read.beast("postburn10k_Bbone_FBD_136taxa_topoAsZhou_FIN.inMYR_MCC_treeAnOut.tre")
FBD1<-drop.tip2(FBD0, "_Anolis_carolinensis")
FBD<-ladderize(FBD1)

NDexp0<-read.beast("postburn10k_Bbone_ND_60taxa_17calExp_BDallDNA.inMYR_MCC_treeAnOut.tre")
NDexp1<-drop.tip2(NDexp0, "_Anolis_carolinensis")
NDexp<-ladderize(NDexp1)

####
# Then drop all the non-PatchClade reps-- giving just the REAL backbones...
# 28 TAXA
####

toKeep <- c("Anomalurus_beecrofti__ANOMALURIDAE__RODENTIA", "Bison_bison__BOVIDAE__CETARTIODACTYLA", "Callithrix_jacchus__CALLITRICHIDAE__PRIMATES", "Calomyscus_baluchi__CALOMYSCIDAE__RODENTIA", "Castor_canadensis__CASTORIDAE__RODENTIA", "Cricetomys_gambianus__NESOMYIDAE__RODENTIA", "Cricetulus_barabensis__CRICETIDAE__RODENTIA", "Dasypus_novemcinctus__DASYPODIDAE__CINGULATA", "Eptesicus_fuscus__VESPERTILIONIDAE__CHIROPTERA", "Equus_caballus__EQUIDAE__PERISSODACTYLA", "Erethizon_dorsatum__ERETHIZONTIDAE__RODENTIA", "Erinaceus_europaeus__ERINACEIDAE__EULIPOTYPHLA", "Felis_catus__FELIDAE__CARNIVORA", "Galeopterus_variegatus__CYNOCEPHALIDAE__DERMOPTERA", "Ictidomys_tridecemlineatus__SCIURIDAE__RODENTIA", "Jaculus_jaculus__DIPODIDAE__RODENTIA", "Manis_pentadactyla__MANIDAE__PHOLIDOTA", "Ornithorhynchus_anatinus__ORNITHORHYNCHIDAE__MONOTREMATA", "Oryctolagus_cuniculus__LEPORIDAE__LAGOMORPHA", "Pteronotus_parnellii__MORMOOPIDAE__CHIROPTERA", "Pteropus_alecto__PTEROPODIDAE__CHIROPTERA", "Rattus_norvegicus__MURIDAE__RODENTIA", "Saccopteryx_bilineata__EMBALLONURIDAE__CHIROPTERA", "Spalax_ehrenbergi__SPALACIDAE__RODENTIA", "Trichechus_manatus__TRICHECHIDAE__SIRENIA", "Tupaia_belangeri__TUPAIIDAE__SCANDENTIA", "Typhlomys_cinereus__PLATACANTHOMYIDAE__RODENTIA", "Vombatus_ursinus__VOMBATIDAE__DIPROTODONTIA")
toDropFBD <- setdiff(FBD0$tip.label,toKeep)
toDropND <- setdiff(NDexp0$tip.label,toKeep)

repFamOrdSp<-read.table("PCreps_toFamOrd_toSp_toSpFam.txt")
names(repFamOrdSp)<-c("tip","famOrd","sp","spFam")

FBD_full0<-read.beast("postburn10k_Bbone_FBD_136taxa_topoAsZhou_FIN.inMYR_MCC_treeAnOut.tre")
FBD_full1<-drop.tip2(FBD_full0, "_Anolis_carolinensis")
FBD_full<-ladderize(FBD_full1)

#FBD0<-read.beast("postburn10k_Bbone_FBD_60taxa_BDallDNA.inMYR_MCC_treeAnOut.tre")
FBD0<-read.beast("postburn10k_Bbone_FBD_136taxa_topoAsZhou_FIN.inMYR_MCC_treeAnOut.tre")
FBD1<-drop.tip2(FBD0, toDropFBD)
FBD<-ladderize(FBD1)
FBD$tip.label<-as.vector(repFamOrdSp[match(FBD$tip.label, as.vector(repFamOrdSp$tip)),"sp"])

NDexp0<-read.beast("postburn10k_Bbone_ND_60taxa_17calExp_BDallDNA.inMYR_MCC_treeAnOut.tre")
NDexp1<-drop.tip2(NDexp0,toDropND)
NDexp<-ladderize(NDexp1)
NDexp$tip.label<-as.vector(repFamOrdSp[match(NDexp$tip.label, as.vector(repFamOrdSp$tip)),"spFam"])

NDuni0<-read.beast("postburn10k_Bbone_ND_60taxa_17calUni_BDallDNA.inMYR_MCC_treeAnOut.tre")
NDuni1<-drop.tip2(NDuni0,toDropND)
NDuni<-ladderize(NDuni1)
NDuni$tip.label<-as.vector(repFamOrdSp[match(NDuni$tip.label, as.vector(repFamOrdSp$tip)),"sp"])

bbonesPC <- c(FBD, NDexp, NDuni)

# GEOSCALE, do each separate, then combine...

pdf(file="Bbone_geoscaleR_FBD-135taxa_topoAsZhou_agesPP_noTipLabels_upwards_new.pdf", width=4, height=8.5, onefile=TRUE)
	FBD_full$root.time <- max(FBD_full$height)
	geoscalePhylo(tree=ladderize(FBD_full,right=FALSE), units=c("Era","Period", "Epoch"), direction="upwards", boxes="Epoch", x.lim=c(-50,FBD_full$root.time), quat.rm=TRUE, cex.tip=0.7, cex.age=0.7, cex.ts=0.7, tick.scale=50, label.offset=0.4, lwd=3, width=1,show.tip.label=FALSE)
	#geoscalePhylo(tree=ladderize(t,right=FALSE), units=c("Era","Period", "Epoch"), direction= "upwards", boxes="Epoch", cex.tip=0.9,cex.age=0.7, cex.ts=0.7,label.offset=0,tick.scale=50,x.lim=c(-0,325),lwd=3, width=2,quat.rm=TRUE, show.tip.label=FALSE)#, ages=ages,ranges=TRUE)
dev.off()

pdf(file="Bbone_geoscaleR_FBD-28taxa_topoAsZhou_agesPP_sp_3colorSupport_fullFBDroot_new.pdf", width=7, height=8.5, onefile=TRUE)
	FBD$root.time <- max(FBD$height)
	geoscalePhylo(tree=ladderize(FBD,right=FALSE), units=c("Era","Period", "Epoch"), direction="upwards", boxes="Epoch", x.lim=c(-50,FBD_full$root.time), quat.rm=TRUE, cex.tip=0.7, cex.age=0.7, cex.ts=0.7, tick.scale=50, label.offset=0.4, lwd=3, width=2)
	HPDbars(FBD, label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))

	nodeLabel<-FBD$posterior
	p <- character(length(FBD$node.label))
	p[nodeLabel >= 0.95] <- "black"
	p[nodeLabel < 0.95 & nodeLabel >= 0.75] <- "gray"
	p[nodeLabel < 0.75] <- "white"
	nodelabels(frame="circle",pch=21, cex=1, bg=p)
	#node.support(FBD$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.7)
	#node.support(FBD$posterior, mode="dots", col = "red", cex=0.7)
	#node.support(FBD$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.71)
	#node.support(FBD$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.7)
	node.support(FBD$posterior, mode="numbers", digits=2,pos="pretty", col = "black", font=3, cex=0.7)#, adj=0.5)
dev.off()

pdf(file="Bbone_geoscaleR_NDexp-28taxa_agesPP_sp_3colorSupport_fullFBDroot_new.pdf", width=7, height=8.5, onefile=TRUE)
#pdf(file="Bbone_geoscaleR_NDexp-28taxa_agesPP_sp_3colorSupport_fullFBDroot_wAges_new.pdf", width=9, height=6, onefile=TRUE)
	NDexp$root.time <- max(NDexp$height)
	geoscalePhylo(tree=ladderize(NDexp,right=FALSE), units=c("Era","Period", "Epoch"), direction="upwards", boxes="Epoch", x.lim=c(-50,NDexp$root.time+25), quat.rm=TRUE, cex.tip=0.7, cex.age=0.7, cex.ts=0.7, tick.scale=50, label.offset=0.4, lwd=3, width=2)
	HPDbars(NDexp, label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))

	nodeLabel<-NDexp$posterior
	p <- character(length(NDexp$node.label))
	p[nodeLabel >= 0.95] <- "black"
	p[nodeLabel < 0.95 & nodeLabel >= 0.75] <- "gray"
	p[nodeLabel < 0.75] <- "white"
	nodelabels(frame="circle",pch=21, cex=1, bg=p)
	#node.support(NDexp$height_median, mode="numbers", digits=1,pos="below", col = "black", font=2, cex=0.7)
	node.support(NDexp$posterior, mode="numbers", digits=2,pos="pretty", col = "black", font=3, cex=0.7)
dev.off()

pdf(file="Bbone_geoscaleR_NDuni-28taxa_agesPP_sp_3colorSupport_fullFBDroot_new.pdf", width=7, height=8.5, onefile=TRUE)
	NDuni$root.time <- max(NDuni$height)
	geoscalePhylo(tree=ladderize(NDuni,right=FALSE), units=c("Era","Period", "Epoch"), direction="upwards", boxes="Epoch", x.lim=c(-50,FBD_full$root.time), quat.rm=TRUE, cex.tip=0.7, cex.age=0.7, cex.ts=0.7, tick.scale=50, label.offset=0.4, lwd=3, width=2)
	HPDbars(NDuni, label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))

	nodeLabel<-NDuni$posterior
	p <- character(length(NDuni$node.label))
	p[nodeLabel >= 0.95] <- "black"
	p[nodeLabel < 0.95 & nodeLabel >= 0.75] <- "gray"
	p[nodeLabel < 0.75] <- "white"
	nodelabels(frame="circle",pch=21, cex=1, bg=p)
	node.support(NDuni$posterior, mode="numbers", digits=2,pos="pretty", col = "black", font=3, cex=0.7)
dev.off()




# ===
# part c: plotting the NDexp backbone as densitree uncertainty
# ===
setwd("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")
library(ape); library(phangorn); library(picante)

# Plot as densitree, but manual...
###
# load mamPhy 100 tree sample
mamPhy_100<-read.nexus("MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_v2_sample100_nexus.trees")
mamPhy100<-lapply(mamPhy_100,ladderize)
class(mamPhy100)<-"multiPhylo"

# load MCC tree
bbone <- "NDexp" #"FBDasZhouEtAl"
mamMCC<-drop.tip(read.nexus(file=paste0("MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_",bbone,"_MCC_v2_target.tre")),"_Anolis_carolinensis")
plottree <- ladderize(mamMCC, right=TRUE)

# define PATCHES
patchSp<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_1spMostGenes_perPC.txt")
newNames<-c("Nesomyidae", "Muridae", "Cricetidae", "Squirrel-related", "Guinea_pig-related", "Eulipotyphla", "Phyllostomid-related", "Vespertilionid-related", "Emballonurid-related", "Yinpterochiroptera", "Marsupialia", "Cetartiodactyla", "Perissodactyla", "Carnivora", "Monotremata", "Pholidota", "Dermoptera", "Anomaluromorpha", "Calomyscidae", "Platacanthomyidae", "Afrotheria", "Xenarthra", "Scandentia", "Primates", "Lagomorpha", "Castorimorpha", "Dipodidae", "Spalacidae")
patchSpPatch<-cbind(patchSp[1],newNames)
colnames(patchSpPatch)<-c("tip","PC")

patchPhy_100<-vector("list",length=100)
for (i in 1:length(mamPhy100)){
	patchPhy_100[[i]]<-drop.tip(mamPhy100[[i]],setdiff(mamPhy100[[i]]$tip.label,patchSpPatch$tip))
	patchPhy_100[[i]]$tip.label<-as.vector(patchSpPatch[match(patchPhy_100[[i]]$tip.label, as.vector(patchSpPatch$tip)),"PC"])
}
class(patchPhy_100)<-"multiPhylo"
write.nexus(patchPhy_100, file="MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_sample100_28PATCHES_patchNames.tre")

patchPhy<-drop.tip(mamMCC,setdiff(mamMCC$tip.label,patchSpPatch$tip))
write.tree(patchPhy, "MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target_28PATCHES_fullNames.tre")
patchPhy$tip.label<-as.vector(patchSpPatch[match(patchPhy$tip.label, as.vector(patchSpPatch$tip)),"PC"])
write.tree(patchPhy, "MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target_28PATCHES_patchNames.tre")

# READ back in the 100 patch trees
patchPhy_100<-read.nexus(file="MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_sample100_28PATCHES_patchNames.tre")

# 100trees: find which trees are A(X+B) vs (A+X)B
nodeSizes<-c()
for(i in 1:length(patchPhy_100)){
  tree<-patchPhy_100[[i]]
  node<-getMRCA(tree,tip=c("Afrotheria","Xenarthra"))
  nodeTips<-tree$tip.label[Descendants(tree,node)[[1]]]  
  nodeSizes[i]<-length(nodeTips)
}
# nodeSizes
# 2 26 
# 51 49 

treeNum<-c(1:length(patchPhy_100))
ID<-cbind.data.frame(treeNum,nodeSizes)
treesAX<-ID[which(ID$nodeSizes==2),]
treesXB<-ID[which(ID$nodeSizes==26),]
tipOrder<-rev(c("Vespertilionid−related", "Emballonurid−related", "Phyllostomid−related", "Yinpterochiroptera", "Pholidota", "Carnivora", "Perissodactyla", "Cetartiodactyla", "Eulipotyphla", "Guinea_pig−related", "Squirrel−related", "Castorimorpha", "Anomaluromorpha", "Dipodidae", "Platacanthomyidae", "Spalacidae", "Calomyscidae", "Nesomyidae", "Cricetidae", "Muridae", "Lagomorpha", "Primates", "Dermoptera", "Scandentia", "Xenarthra", "Afrotheria", "Marsupialia", "Monotremata"))
	#"Monotremata", "Marsupialia", "Afrotheria", "Xenarthra", "Eulipotyphla", "Pholidota", "Carnivora", "Perissodactyla", "Cetartiodactyla", "Yinpterochiroptera", "Phyllostomid−related", "Emballonurid−related", "Vespertilionid−related", "Primates", "Dermoptera", "Scandentia", "Lagomorpha", "Guinea_pig−related", "Squirrel−related", "Anomaluromorpha", "Castorimorpha", "Dipodidae", "Platacanthomyidae", "Spalacidae", "Calomyscidae", "Nesomyidae", "Cricetidae", "Muridae")
	
patch_XB<-patchPhy_100[treesXB[,1]]
class(patch_XB)<-"multiPhylo"
patchXB_l<-lapply(patch_XB,ladderize)
patchXB<-lapply(patchXB_l,rotateConstr,constraint=tipOrder)

patch_AX<-patchPhy_100[treesAX[,1]]
class(patch_AX)<-"multiPhylo"
patchAX_l<-lapply(patch_AX,ladderize)
patchAX<-lapply(patchAX_l,rotateConstr,constraint=tipOrder)

# load the MCC for starting the plot 
#patchPhy<-read.tree("MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_v2_sample100_28PATCHES_patchNames.tre")
rootTime<-max(node.age(patchPhy)$ages)

# with densiTree (phangorn)::
# red-blue together
# pdf(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_28PATCHES_densiTree.pdf",sep=""),width=7, height=10)
# densiTree(patchXB, type="phylogram",col="red",alpha=0.1,  consensus=patchPhy, no.margin=FALSE, cex=0.9, label.offset=0.4, direction="rightwards", show.tip.label=FALSE)
# par(new=TRUE)
# densiTree(patchAX, type="phylogram",col="blue",alpha=0.1,  consensus=patchPhy, no.margin=FALSE, cex=0.9, label.offset=0.4, direction="rightwards", show.tip.label=FALSE)
# dev.off()

# MANUALLY AS LINEAR...
edgeCol1<-rgb(1,0,0,alpha=0.1)
edgeCol2<-rgb(0,0,1,alpha=0.1)

pdf(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_28PATCHES_LINEAR_redBlue_withTips_clade_OverlapNew.pdf",sep=""),width=10, height=7)
#jpeg(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_28PATCHES_LINEAR_redBlue_withTips_clade_OverlapNew.jpg",sep=""),width=10, height=7, units="in",quality=100, res=600)

#layout(matrix(1:3,1,3),widths=c(0.42,0.16,0.42))
layout(matrix(1:2,1,2),widths=c(0.42,0.10))
par(mar=c(0,0,0,0))
plot(patchXB[[1]], type="cladogram", cex=0.7, font=2, direction="rightwards",label.offset=0.5,tip.color="black",edge.col=edgeCol1,show.tip.label=FALSE,no.margin=TRUE)
for (j in 2:length(patchXB)){
	par(new=TRUE)
	plot(patchXB[[j]], type="cladogram", cex=0.7, font=2, direction="rightwards", label.offset=0.5,tip.color=grey(0.5, alpha=0.01),edge.col=edgeCol1,show.tip.label=FALSE,no.margin=TRUE)
}
#axisPhylo(cex=0.6, line=-1)
#plot.new()
#plot.window(xlim=c(-1,1),ylim=c(1, length(patchXB[[1]]$tip.label)))
#par(cex=0.6,mar=c(0,0,0,0))
#text(rep(0,length(patchXB[[1]]$tip.label)), 1:length(patchXB[[1]]$tip.label),(tipOrder))

par(new=TRUE)
plot(patchAX[[1]], type="cladogram", cex=0.7, font=2, direction="rightwards", label.offset=0.5,tip.color="white",edge.col=edgeCol2,show.tip.label=FALSE,no.margin=TRUE)
for (j in 2:length(patchAX)){
	par(new=TRUE)
	plot(patchAX[[j]], type="cladogram", cex=0.7, font=2, direction="rightwards", label.offset=0.5,tip.color=grey(0.5, alpha=0.01),edge.col=edgeCol2,show.tip.label=FALSE,no.margin=TRUE)
}
axisPhylo(cex=0.6, line=-1)

plot.new(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot.window(xlim=c(-2.1,0.1),ylim=c(1, length(patchXB[[1]]$tip.label)))
par(cex=0.6)
text(rep(0,length(patchXB[[1]]$tip.label)), 1:length(patchXB[[1]]$tip.label),(tipOrder), pos=4)

dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy MS1 -- Upham et al. 2019 -- PLOS Biology
###
# Supplemental Figures S1-S9 - plotting code
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SUPP FIGURES PLOTTING
###

# ============
# FIG S1 and S2
# ============
# By hand

# ============
# FIG S3 == per GENUS / FAMILY fractions plotted as barplots
# ============

setwd("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")
mamData<-read.csv(file="traits_datReady_mamPhy_5911species.csv",header=TRUE)[,c(1:11)]

# GENERA
genNames<-names(table(mamData$gen))
numSp<-c();numSampled<-c()
for(j in 1:length(genNames)){
	genBySp<-mamData[which(mamData$gen==genNames[j]),c("Species_Name","genes")]
	numSp[j]<-length(genBySp[,1])
	numSampled[j]<-length(genBySp[which(genBySp$genes!=0),1])
}
allRes_unSort<-cbind.data.frame(species=numSp,sampled=numSampled)
allRes<-allRes_unSort[order(allRes_unSort$species, decreasing=TRUE),]
rownames(allRes)<-genNames[order(allRes_unSort$species, decreasing=TRUE)]

pdf(file="DNAsampling_perGenus_mamPhy.pdf", onefile=TRUE, width=6, height=4)
	barplot(height=allRes$species, width=0.1,cex.names=0.3, names.arg=allRes$genNames, angle=45, col=grey(0.15), border=NA)
	barplot(add=TRUE,height=allRes$sampled, width=0.1,col="cyan3", border=NA)

	x<-barplot(height=allRes$species[c(1:20)], width=0.1,cex.names=0.3, names.arg=rownames(allRes)[c(1:20)], angle=45, col=grey(0.15), border=NA, xaxt="n")
	barplot(add=TRUE,height=allRes$sampled[c(1:20)], width=0.1,col="cyan3", border=NA)
	legend(x=0.8, y=150, legend=c("total species","DNA sampled"), fill=c(grey(0.15),"cyan3"))
	text(cex=1, x=x-0, y=-8, rownames(allRes)[c(1:20)], xpd=TRUE, srt=45, font=3, adj=1)
dev.off()

# FAMILIES
genNames<-names(table(mamData$fam))
numSp<-c();numSampled<-c()
for(j in 1:length(genNames)){
	genBySp<-mamData[which(mamData$fam==genNames[j]),c("Species_Name","genes")]
	numSp[j]<-length(genBySp[,1])
	numSampled[j]<-length(genBySp[which(genBySp$genes!=0),1])
}
allRes_unSort<-cbind.data.frame(species=numSp,sampled=numSampled)
allRes<-allRes_unSort[order(allRes_unSort$species, decreasing=TRUE),]
rownames(allRes)<-genNames[order(allRes_unSort$species, decreasing=TRUE)]

pdf(file="DNAsampling_perFamily_mamPhy.pdf", onefile=TRUE, width=6, height=4)
	barplot(height=allRes$species, width=0.1,cex.names=0.3, names.arg=allRes$genNames, angle=45, col=grey(0.15), border=NA)
	barplot(add=TRUE,height=allRes$sampled, width=0.1,col="cyan3", border=NA)

	x<-barplot(height=allRes$species[c(1:20)], width=0.1,cex.names=0.3, names.arg=rownames(allRes)[c(1:20)], angle=45, col=grey(0.15), border=NA, xaxt="n")
	barplot(add=TRUE,height=allRes$sampled[c(1:20)], width=0.1,col="cyan3", border=NA)
	#legend(x=0.8, y=150, legend=c("total species","DNA sampled"), fill=c(grey(0.15),"cyan3"))
	text(cex=0.7, x=x-0, y=-8, rownames(allRes)[c(1:20)], xpd=TRUE, srt=45, font=3, adj=1)
dev.off()


# ORDERS
genNames<-names(table(mamData$ord))
numSp<-c();numSampled<-c()
for(j in 1:length(genNames)){
	genBySp<-mamData[which(mamData$ord==genNames[j]),c("Species_Name","genes")]
	numSp[j]<-length(genBySp[,1])
	numSampled[j]<-length(genBySp[which(genBySp$genes!=0),1])
}
allRes_unSort<-cbind.data.frame(species=numSp,sampled=numSampled)
allRes<-allRes_unSort[order(allRes_unSort$species, decreasing=TRUE),]
rownames(allRes)<-genNames[order(allRes_unSort$species, decreasing=TRUE)]


pdf(file="DNAsampling_perOrder_mamPhy.pdf", onefile=TRUE, width=6, height=4)
	barplot(height=allRes$species, width=0.1,cex.names=0.3, names.arg=allRes$genNames, angle=45, col=grey(0.15), border=NA)
	barplot(add=TRUE,height=allRes$sampled, width=0.1,col="cyan3", border=NA)

	x<-barplot(height=allRes$species[c(1:10)], width=0.1,cex.names=0.3, angle=45, col=grey(0.15), border=NA, xaxt="n")
	barplot(add=TRUE,height=allRes$sampled[c(1:10)], width=0.1,col="cyan3", border=NA)
	#legend(x=0.8, y=150, legend=c("total species","DNA sampled"), fill=c(grey(0.15),"cyan3"))
	text(cex=0.7, x=x-0, y=-8, rownames(allRes)[c(1:10)], xpd=TRUE, srt=45, font=3, adj=1)
dev.off()



# ============
# FIG S4 == plotting node-dated backbone trees
# ============
	library(ape); library(phyloch); library(phytools)
	library(paleotree); library(strap); #library(OutbreakTools)

setwd("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/11_makingFullPosteriors")

NDexp0<-read.beast("postburn10k_Bbone_ND_60taxa_17calExp_BDallDNA.inMYR_MCC_treeAnOut.tre")
NDexp1<-drop.tip2(NDexp0,"_Anolis_carolinensis")
NDexp<-ladderize(NDexp1)

NDuni0<-read.beast("postburn10k_Bbone_ND_60taxa_17calUni_BDallDNA.inMYR_MCC_treeAnOut.tre")
NDuni1<-drop.tip2(NDuni0,"_Anolis_carolinensis")
NDuni<-ladderize(NDuni1)

# GEOSCALE for the 59 TAXA
###
# NDexp
#####
pdf(file="Bbone_geoscaleR_NDexp-59taxa_agesPP_3colorSupport_median_new.pdf", width=8.5, height=11, onefile=TRUE)
	NDexp$root.time <- max(NDexp$height)
	geoscalePhylo(tree=NDexp, units=c("Era","Period", "Epoch"), direction="rightwards", boxes="Epoch", x.lim=c(-150,NDexp$root.time+30), quat.rm=TRUE, cex.tip=0.7, cex.age=0.7, cex.ts=0.7, tick.scale=50, label.offset=0.4, lwd=3, width=2)
	HPDbars(NDexp, label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))

	NDexp$node.label<-NDexp$posterior
	p <- character(length(NDexp$node.label))
	p[NDexp$node.label >= 0.95] <- "black"
	p[NDexp$node.label < 0.95 & NDexp$node.label >= 0.75] <- "gray"
	p[NDexp$node.label < 0.75] <- "white"
	nodelabels(pch=21, cex=1, bg=p)
	node.support(NDexp$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.7)
	#node.support(NDexp$posterior, mode="dots", col = "red", cex=0.7)
	#node.support(NDexp$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.71)
	#node.support(NDexp$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.7)
dev.off()

###
# NDuni
#####
pdf(file="Bbone_geoscaleR_NDuni-59taxa_agesPP_3colorSupport_median_new.pdf", width=8.5, height=11, onefile=TRUE)
	NDuni$root.time <- max(NDuni$height)
	geoscalePhylo(tree=NDuni, units=c("Era","Period", "Epoch"), direction="rightwards", boxes="Epoch", x.lim=c(-150,NDexp$root.time+30), quat.rm=TRUE, cex.tip=0.7, cex.age=0.7, cex.ts=0.7, tick.scale=50, label.offset=0.4, lwd=3, width=2)
	HPDbars(NDuni, label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))

	NDuni$node.label<-NDuni$posterior
	p <- character(length(NDuni$node.label))
	p[NDuni$node.label >= 0.95] <- "black"
	p[NDuni$node.label < 0.95 & NDuni$node.label >= 0.75] <- "gray"
	p[NDuni$node.label < 0.75] <- "white"
	nodelabels(pch=21, cex=1, bg=p)
	node.support(NDuni$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.7)
	#node.support(NDuni$posterior, mode="dots", col = "red", cex=0.7)
	#node.support(NDuni$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.71)
	#node.support(NDuni$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.7)
dev.off()



# ============
# FIG S5 == plotting tip-dated backbone trees
# ============

	library(ape); library(phyloch); library(phytools)
	library(paleotree); library(strap); #library(OutbreakTools)

setwd("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/11_makingFullPosteriors")

#FBD0<-read.beast("postburn10k_Bbone_FBD_60taxa_BDallDNA.inMYR_MCC_treeAnOut.tre")
FBD0<-read.beast("postburn10k_Bbone_FBD_136taxa_topoAsZhou_FIN.inMYR_MCC_treeAnOut.tre")
FBD1<-drop.tip2(FBD0, "_Anolis_carolinensis")
FBD<-ladderize(FBD1)

stratRanges<-read.table("tipStratiRanges.txt")
colnames(stratRanges)<-c("tip","min","max")

ages<-cbind(stratRanges$min,stratRanges$max)
rownames(ages)<-stratRanges$tip
colnames(ages)<-c("FAD","LAD")

FBD$root.time <- max(FBD$height)
#geoscalePhylo(tree=FBD,ages=ages,ranges=TRUE,cex.tip=0.8,cex.ts=0.55,cex.age=0.5,width=2)

###
t<-FBD
num_taxa <- length(t$tip.label)
t$root.time <- max(t$height)

stem_length <- 0
origin_HPD<-c(t$"height_95%_HPD_MIN"[1],t$"height_95%_HPD_MAX"[1])
names(origin_HPD)<-c("lower","upper")

display_all_node_bars <- TRUE

names_list <-vector()
for (name in t$tip){
  v <- strsplit(name,"_")[[1]]
  if(display_all_node_bars){
  	names_list = c(names_list,name)
  }
  else if(v[length(v)]=="0"){
  	names_list = c(names_list,name)
  }
}

nids <- vector()
pos <- 1
len_nl <- length(names_list)
for(n in names_list){
  for(nn in names_list[pos:len_nl]){
    if(n != nn){
      m <- getMRCA(t,c(n,nn))
      if(m %in% nids == FALSE){
        nids <- c(nids,m)
      }
    }
  }
  pos<-pos+1
}

root_max <- t$"height_95%_HPD_MAX"[1]
x_max <- origin_HPD[2] * 0.02 + origin_HPD[2]

stratRanges<-read.table("tipStratiRanges.txt")
#stratRanges<-read.table("tipStratiRanges_all4098Species.txt")
colnames(stratRanges)<-c("tip","min","max")
ages<-cbind(stratRanges$min,stratRanges$max)
rownames(ages)<-stratRanges$tip
colnames(ages)<-c("FAD","LAD")


pdf("geoscaled_Mammalia_136taxa_topoAsZhou_FIN_wFossilAges_new.pdf", width=20, height=30)
geoscalePhylo(tree=ladderize(t,right=FALSE), units=c("Era","Period","Epoch"), boxes="Period", y.lim=c(0,140),cex.tip=1,cex.age=1.5, cex.ts=1.5,label.offset=0,tick.scale=50,x.lim=c(-150,x_max),lwd=1, width=1,quat.rm=TRUE)#, ages=ages,ranges=TRUE)

lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)

for(nv in nids){
  bar_xx_a <- c(lastPP$xx[nv]+t$height[nv-num_taxa]-t$"height_95%_HPD_MIN"[nv-num_taxa], lastPP$xx[nv]-(t$"height_95%_HPD_MAX"[nv-num_taxa]-t$height[nv-num_taxa]))
  lines(bar_xx_a,c(lastPP$yy[nv],lastPP$yy[nv]),col=rgb(0,0,1,alpha=0.6),lwd=2)
}

t$node.label<-t$posterior
p <- character(length(t$node.label))
p[t$node.label >= 0.95] <- "black"
p[t$node.label < 0.95 & t$node.label >= 0.75] <- "gray"
p[t$node.label < 0.75] <- "white"
nodelabels(pch=21, cex=1.2, bg=p, lwd=0.2)
node.support(t$height, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.9)


dev.off()



# ============
# FIG S6 == plotting PRUNED node- and tip-dated backbone trees for comparison
# ============


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


#########
# THEN-- For the 28-taxon trees... Get the TABLE values for this.

pdf(file="Bbone_FBD-28taxa_nodelabels.pdf", width=6, height=8, onefile=TRUE)
plot(FBD, cex=0.6)
nodelabels(cex=0.6)
dev.off()

pdf(file="Bbone_NDexp-28taxa_nodelabels_shorter.pdf", width=6, height=6, onefile=TRUE)
plot(NDexp, cex=0.6)
nodelabels(cex=0.6)
dev.off()

pdf(file="Bbone_NDuni-28taxa_nodelabels.pdf", width=6, height=8, onefile=TRUE)
plot(NDuni, cex=0.6)
nodelabels(cex=0.6)
dev.off()

#>> nodes are not exactly identical...
#all.equal.phylo(FBD,NDexp, use.edge.length=FALSE, use.tip.label=FALSE, index.return=TRUE)

NDexp_NDuni <- matchNodes(NDexp,NDuni, "descendants")

NDexp_FBD <- matchNodes(NDexp,FBD, "descendants")

NDuni_FBD <- matchNodes(NDuni,FBD, "descendants")


# Get node data:
FBD_table<-cbind(FBD$"height_95%_HPD_MIN", FBD$"height_95%_HPD_MAX",FBD$height,FBD$height_median)
colnames(FBD_table)<-c("height_95%_HPD_MIN","height_95%_HPD_MAX","height_mean","height_median")
rownames(FBD_table)<-NDexp_NDuni[,1]

NDexp_table<-cbind(NDexp$"height_95%_HPD_MIN", NDexp$"height_95%_HPD_MAX",NDexp$height,NDexp$height_median)
colnames(NDexp_table)<-c("height_95%_HPD_MIN","height_95%_HPD_MAX","height_mean","height_median")
rownames(NDexp_table)<-NDexp_NDuni[,1] # this is just 60:117 continuous

NDuni_table<-cbind(NDuni$"height_95%_HPD_MIN", NDuni$"height_95%_HPD_MAX",NDuni$height,NDuni$height_median)
colnames(NDuni_table)<-c("height_95%_HPD_MIN","height_95%_HPD_MAX","height_mean","height_median")
rownames(NDuni_table)<-NDexp_NDuni[,1] # this is just 60:117 continuous

# BASE this on the NDexp topology...
mins <- cbind(NDexp_table[,"height_95%_HPD_MIN"],NDuni_table[match(NDexp_NDuni[,2],rownames(NDuni_table)),"height_95%_HPD_MIN"],FBD_table[match(NDexp_FBD[,2],rownames(FBD_table)),"height_95%_HPD_MIN"])
colnames(mins)<-c("NDexp_min","NDuni_min","FBD_min")

maxs <- cbind(NDexp_table[,"height_95%_HPD_MAX"],NDuni_table[match(NDexp_NDuni[,2],rownames(NDuni_table)),"height_95%_HPD_MAX"],FBD_table[match(NDexp_FBD[,2],rownames(FBD_table)),"height_95%_HPD_MAX"])
colnames(maxs)<-c("NDexp_max","NDuni_max","FBD_max")

means <- cbind(NDexp_table[,"height_mean"],NDuni_table[match(NDexp_NDuni[,2],rownames(NDuni_table)),"height_mean"],FBD_table[match(NDexp_FBD[,2],rownames(FBD_table)),"height_mean"])
colnames(means)<-c("NDexp_mean","NDuni_mean","FBD_mean")

medians <- cbind(NDexp_table[,"height_median"],NDuni_table[match(NDexp_NDuni[,2],rownames(NDuni_table)),"height_median"],FBD_table[match(NDexp_FBD[,2],rownames(FBD_table)),"height_median"])
colnames(medians)<-c("NDexp_median","NDuni_median","FBD_median")

summary<-cbind(mins,maxs,means,medians)

sum_Ordered<-summary[,c("NDexp_median","NDexp_min","NDexp_max","NDuni_median","NDuni_min","NDuni_max","FBD_median","FBD_min","FBD_max")]

write.table(sum_Ordered,"bbone_nodeCompare_NDexp-NDuni-FBD_withFBD-topoAsZhou_ordered_28taxa_medians.txt")



# ============
# FIG S7 == plotting effect of gene sampling upon tip DR estimates
# ============

# look at GENE-SAMPLING at the TIP LEVEL...
###

# directory and source
dirname = "/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/"
setwd(dirname)

# which backbone?
bbone<- "NDexp" #"FBD" # 

# other packages
 library(ape); library(geiger); library(picante); library(phangorn); library(ggplot2); library(phylolm); library(phytools)

#==================
# Load in stats about tip TAXONOMY and GENE SAMPLING
#==================
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
head(cladesDR)
sorted<-cladesDR[order(cladesDR$tiplabel),]

taxonomyVars<-c("tiplabel","gen", "fam", "ord", "higher", "genes","extinct.", "harmMeans", "medians","variance", "cv")
taxFinal<-sorted[,taxonomyVars]

#==================
# load in stats about TIP TRAITS
#==================
tipDataAll_orig<-read.table(file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)
tipDataAll<-tipDataAll_orig[order(tipDataAll_orig$tiplabel),]

allDatFinal_all<-cbind.data.frame(taxFinal,tipDataAll) #[,varsSelected],tipDataAll[,catsSelected])
allDatFinal_extant<-allDatFinal_all[which(allDatFinal_all$extinct.==0),] # keeping EXTANT species only.

allDatFinal<-allDatFinal_extant[which(allDatFinal_extant$MarineOrNot==0),] #allDatFinal_extant$ord!="CHIROPTERA" & 
	# gives 5675 species


toPlot<-allDatFinal_all[,c("harmMeans","medians","variance","cv","genes")]
toPlot_DNAonly<-allDatFinal_all[which(allDatFinal_all$genes!=0),c("harmMeans","medians","variance","cv","genes")]

pdf(file=paste("test-for-tipDR-bias_10kTrees_withMedian.pdf",sep=""), onefile=TRUE, width=9, height=10)	

#layout(matrix(c(1:8), 2, 4, byrow = TRUE), widths=rep(3,8), heights=rep(1,8))
layout( matrix( c(1:4), ncol=2, nrow=2, byrow=TRUE) )
#par(oma = c(2,2,2,2) + 0.1, mar = c(0,1,0,1) + 0.1) #c(bottom, left, top, right)

# VARIANCE
plot(variance ~ genes, data=toPlot, pch=".", ylab="Variance in tip DR (10,000 trees)", xlab="Gene sampling in supermatrix", main="With unsampled")
corr1<-cor.test(y=toPlot$variance, x=toPlot$genes, method="spearman", exact=FALSE)
mtext(text=paste("n = ",length(toPlot[,1]),"; r = ", round(corr1$estimate[[1]],3),"; P = ", round(corr1$p.value[[1]],3),sep=""), side=3, )

plot(variance ~ genes, data=toPlot_DNAonly, pch=".", ylab="Variance in tip DR (10,000 trees)", xlab="Gene sampling in supermatrix", main="Excluding unsampled")
corr2<-cor.test(y=toPlot_DNAonly$variance, x=toPlot_DNAonly$genes, method="spearman", exact=FALSE)
mtext(text=paste("n = ",length(toPlot_DNAonly[,1]),"; r = ", round(corr2$estimate[[1]],3),"; P = ", round(corr2$p.value[[1]],3),sep=""), side=3, )

#lines(lowess(x= toPlot_DNAonly$genes, y = toPlot_DNAonly$variance, f = 1/3, iter = 100))#, delta = 0.01 * diff(range(datToPlot[,j]))) )
#fit2<-lm(variance ~ genes, data=toPlot_DNAonly)

# MEDIAN
plot(medians ~ genes, data=toPlot, pch=".", ylab="Median tip DR (10,000 trees)", xlab="Gene sampling in supermatrix", main="With unsampled")
corr1<-cor.test(y=toPlot$medians, x=toPlot$genes, method="spearman", exact=FALSE)
mtext(text=paste("n = ",length(toPlot[,1]),"; r = ", round(corr1$estimate[[1]],3),"; P = ", round(corr1$p.value[[1]],3),sep=""), side=3, )

plot(medians ~ genes, data=toPlot_DNAonly, pch=".", ylab="Median tip DR (10,000 trees)", xlab="Gene sampling in supermatrix", main="Excluding unsampled")
corr2<-cor.test(y=toPlot_DNAonly$medians, x=toPlot_DNAonly$genes, method="spearman", exact=FALSE)
mtext(text=paste("n = ",length(toPlot_DNAonly[,1]),"; r = ", round(corr2$estimate[[1]],3),"; P = ", round(corr2$p.value[[1]],3),sep=""), side=3, )

## HARM MEAN
#plot(harmMeans ~ genes, data=toPlot, pch=".", ylab="Mean tip DR (10,000 trees)", xlab="Gene sampling in supermatrix", main="With unsampled")
#corr1<-cor.test(y=toPlot$harmMeans, x=toPlot$genes, method="spearman", exact=FALSE)
#mtext(text=paste("n = ",length(toPlot[,1]),"; r = ", round(corr1$estimate[[1]],3),"; P = ", round(corr1$p.value[[1]],3),sep=""), side=3, )
#
#plot(harmMeans ~ genes, data=toPlot_DNAonly, pch=".", ylab="Mean tip DR (10,000 trees)", xlab="Gene sampling in supermatrix", main="Excluding unsampled")
#corr2<-cor.test(y=toPlot_DNAonly$harmMeans, x=toPlot_DNAonly$genes, method="spearman", exact=FALSE)
#mtext(text=paste("n = ",length(toPlot_DNAonly[,1]),"; r = ", round(corr2$estimate[[1]],3),"; P = ", round(corr2$p.value[[1]],3),sep=""), side=3, )


dev.off()


# ============
# FIG S8 == plotting the per-author contributions to this supermatrix study
# ============

# Code to extract just the NCBI accession numbers from the 31 genes in our final list
######

library(dplyr)

# set wd
dirname<-"/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses"
setwd(dirname)
	
	# load in the full results
	authorDatAll2<-read.csv(file="/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/NCBI_authorFetching/ALL31genes_authorRes_firstBlock_Ready_reParsed_Mer11Mod.csv", header=FALSE)
	colnames(authorDatAll2)<-c("gene","NCBI","authors")

	# try just the authors alone
	authorsSplit2<-strsplit( as.character(authorDatAll2$authors), split="\\|", perl=FALSE) # double escape!
	authorsSplit_ALL2<-unlist(authorsSplit2)
		# dim(as.data.frame(table(authorDatAll2$authors)))
		# [1] 1934    2 <<< 1,934 individual studies contributed sequences
		# head(as.data.frame(table(authorDatAll2$authors)[order(table(authorDatAll2$authors), decreasing=TRUE)]))
			# TOP 5 studies::
			###
			# Meredith,R.W.|Janecka,J.E.|Gatesy,J.|Springer,M.S.|Murphy,W.J. 
			# 1963 
			# Perelman,P.|Johnson,W.E.|Roos,C.|Seuanez,H.N.|Horvath,J.E.|Moreira,M.A.|Kessing,B.|Pontius,J.|Roelke,M.|Rumpler,Y.|Schneider,M.P.|Silva,A.|O'Brien,S.J.|Pecon-Slattery,J.
			# 1507
			# Murphy,W.J.|Eizirik,E.|Johnson,W.E.|Zhang,Y.P.|Ryder,O.A.|O'Brien,S.J.
			# 625
			# Schenk,J.J.|Rowe,K.C.|Steppan,S.J.
			# 392
			# Hassanin,A.|Delsuc,F.|Ropiquet,A.|Hammer,C.|JansenvanVuuren,B.|Matthee,C.|Ruiz-Garcia,M.|Catzeflis,F.|Areskoug,V.|Nguyen,T.T.|Couloux,A.
			# 288
			# Meredith,R.W.|Westerman,M.|Springer,M.S.
			# 233

	tableSort2<-table(authorsSplit_ALL2)[order(table(authorsSplit_ALL2), decreasing=TRUE)]
		# dim(as.data.frame(tableSort2))
		# [1] 6069    2 <<< 6,069 individual authors contributed to the studies

	tableSort2[1:30]
		#       Murphy,W.J.     Springer,M.S.      O'Brien,S.J.     Meredith,R.W. 
		#              3516              3222              2640              2512 
		#      Johnson,W.E.         Gatesy,J.      Janecka,J.E.           Roos,C. 
		#              2465              2231              2128              1730 
		# Pecon-Slattery,J.      Seuanez,H.N.      Horvath,J.E.        Rumpler,Y. 
		#              1607              1542              1538              1529 
		#      Moreira,M.A.    Schneider,M.P.        Kessing,B.       Perelman,P. 
		#              1510              1509              1507              1507 
		#        Pontius,J.         Roelke,M.          Silva,A.        Eizirik,E. 
		#              1507              1507              1507              1147 
		#        Ryder,O.A.        Zhang,Y.P.      Steppan,S.J.        Wayne,R.K. 
		#               936               820               775               664 
		#         Rowe,K.C.      Teeling,E.C.         Madsen,O.          Lim,B.K. 
		#               647               641               586               580 
		#      Westerman,M.      Koepfli,K.P. 
		#               544               520 

	# Just FIRST AUTHORS?
	firstAuthors<-c()
	for(j in 1:length(authorsSplit2)){
		firstAuthors[j]<-authorsSplit2[[j]][1]
	}
	as.data.frame(table(firstAuthors))

	tableSort_firstAuthors<-table(firstAuthors)[order(table(firstAuthors), decreasing=TRUE)]
	tableSort_firstAuthors[1:30]
		# Meredith,R.W.    Perelman,P.    Murphy,W.J.    Schenk,J.J.    Hassanin,A. 
		#          2318           1507            684            392            363 
		#  Teeling,E.C.      Sato,J.J.     Eizirik,E.   Koepfli,K.P.     Foley,N.M. 
		#           342            279            264            249            231 
		#   Dumont,E.R.      Rowe,K.C.    Fulton,T.L.       Dubey,S.     Arnason,U. 
		#           222            215            204            201            182 
		#  Steppan,S.J.   Almeida,F.C.  Mitchell,K.J.     Jansa,S.A.          He,K. 
		#           181            177            177            173            149 
		#         Wu,S.  Springer,M.S. Esselstyn,J.A.     Clare,E.L.     Upham,N.S. 
		#           141            140            124            119            117 
		#     Voss,R.S.   Velazco,P.M.      Eger,J.L.  Giannini,N.P.   Johnson,W.E. 
		#           117            111            108            101            100 

	# Just LAST AUTHORS?
	lastAuthors<-c()
	for(j in 1:length(authorsSplit2)){
		lastAuthors[j]<-tail(authorsSplit2[[j]],n=1)
	}
	as.data.frame(table(lastAuthors))
	
	tableSort_lastAuthors<-table(lastAuthors)[order(table(lastAuthors), decreasing=TRUE)]
	tableSort_lastAuthors[1:30]
    	#     Murphy,W.J.    Pecon-Slattery,J.         O'Brien,S.J. 
    	#            2398                 1507                 1009 
    	#   Springer,M.S.         Steppan,S.J.           Wayne,R.K. 
    	#             668                  595                  309 
    	#      Couloux,A. Van Den Bussche,R.A.           Jansa,S.A. 
    	#             288                  280                  273 
    	#       Suzuki,H.         Teeling,E.C.          Strobeck,C. 
    	#             268                  259                  254 
    	#    Davalos,L.M.       Patterson,B.D.         Douzery,E.J. 
    	#             251                  224                  218 
    	#    Simmons,N.B.          Hebert,P.D.             Vogel,P. 
    	#             216                  203                  200 
    	#  Borisenko,A.V.            Cooper,A.         Bradley,R.D. 
    	#             198                  180                  157 
    	#      Organ,C.L.            Cook,J.A.             Janke,A. 
    	#             137                  135                  129 
    	#      Baker,R.J.              Roos,C.             Mayer,F. 
    	#             128                  120                  115 
    	#       Voss,R.S.            Gatesy,J.           Weksler,M. 
    	#             110                  109                  105 

#######
# VISUALIZE
pdf(file="authorRes_barplot_firstAuthors_top30.pdf", width=6, height=6)
	bp<-barplot(tableSort_firstAuthors[1:30], names.arg="", ylim=c(0,2750))
	par(new=TRUE, xpd=NA)
	text(x=bp,y=(as.vector(tableSort_firstAuthors[1:30])+50),names(tableSort_firstAuthors[1:30]), adj=0, offset=1, srt=45, cex=0.75)
	mtext(side=2, text="Number of DNA sequences", line=3)
	mtext(side=1, text="Author frequency", line=1)
	mtext(side=3, text="First authors (top 30 in supermatrix)", line=-3, font=2)
dev.off()


pdf(file="authorRes_barplot_first-and-lastAuthors_top30.pdf", width=12, height=6)
	par(mar = c(3,4,0.5,0.5),oma=rep(0.5,4), xpd=NA)
	layout(matrix(c(1, 2), 1, 2, byrow = TRUE), heights = c(1,1,1,1))

	bp<-barplot(tableSort_firstAuthors[1:30], names.arg="", ylim=c(0,2900))
	#par(new=TRUE, xpd=NA,mar = c(3,4,0,0),oma=rep(0.5,4))
	text(x=bp,y=(as.vector(tableSort_firstAuthors[1:30])+50),names(tableSort_firstAuthors[1:30]), adj=0, offset=1, srt=45, cex=0.75)
	mtext(side=2, text="Number of DNA sequences", line=3)
	mtext(side=1, text="Author frequency", line=1)
	mtext(side=3, text="First authors", line=-2, font=2)
	mtext(side=3, text="(top 30 contributors in supermatrix)", line=-3, font=1, cex=0.9)

	bp<-barplot(tableSort_lastAuthors[1:30], names.arg="", ylim=c(0,2900))#, yaxt="n")
	text(x=bp,y=(as.vector(tableSort_lastAuthors[1:30])+50),names(tableSort_lastAuthors[1:30]), adj=0, offset=1, srt=45, cex=0.75)
	#mtext(side=2, text="Number of DNA sequences", line=3)
	mtext(side=1, text="Author frequency", line=1)
	mtext(side=3, text="Last authors", line=-2, font=2)
	mtext(side=3, text="(top 30 contributors in supermatrix)", line=-3, font=1, cex=0.9)
dev.off()


# TREEMAPS
# library
library(treemap)
library(d3treeR)
	#install.packages("remotes")
	#remotes::install_github("d3treeR/d3treeR")
	#install.packages("/Users/nate/Downloads/gridSVG_1.6-0.tar.gz", repos = NULL, type="source")

pdf(file="authorRes_TREEMAP_first-and-lastAuthors_top30.pdf", width=12, height=12, onefile=TRUE)
	#par(mar = c(3,4,0.5,0.5),oma=rep(0.5,4), xpd=NA)
	#layout(matrix(c(1, 2), 1, 2, byrow = TRUE), heights = c(1,1,1,1))

plotdata1 <- as.data.frame(tableSort_firstAuthors[1:50])
plotdata <- cbind(plotdata1,LABEL=paste(plotdata1$firstAuthors, plotdata1$Freq, sep = "\n"))
treemap(plotdata, algorithm="pivotSize",#drop.unused.levels = FALSE,
       index=c("LABEL"),type="index", inflate.labels=FALSE,#align.labels=c("center", "left"),
       vSize="Freq", fontsize.labels = 16,bg.labels=220,fontsize.title = 22,palette="YlOrRd",
       title="First authors")

plotdata1 <- as.data.frame(tableSort_lastAuthors[1:50])
plotdata <- cbind(plotdata1,LABEL=paste(plotdata1$lastAuthors, plotdata1$Freq, sep = "\n"))
treemap(plotdata, algorithm="pivotSize",#drop.unused.levels = FALSE,
       index=c("LABEL"),type="index",inflate.labels=FALSE,#align.labels=c("center", "left"),
       vSize="Freq",fontsize.labels = 16,bg.labels=220,fontsize.title = 22,palette="YlOrRd",
       title="Last authors")

dev.off()



# ============
# FIG S9 and S10 == plotting the tip-level big trees (MCC) with node support and tip labels readable.
# ============

#####
# Plot MCC trees made with TreeAnnotator.
setwd("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/11_makingFullPosteriors/10k_tree_andDR_distributions_FINAL")
library(ape); library(phyloch)

tree <- read.beast("MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre")
#tree <- read.beast("MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre")

tree2 <- drop.tip2(tree,"_Anolis_carolinensis")
mamPhy <- ladderize(tree2, right=TRUE)

#missing<-read.table("MamPhy_FIN4_1813sp_missing_LIST.txt", header=FALSE)
#tipColors <- rep("black",5912)
#tipColors[match(missing$V1,tree$tip.label)] <- "red" 

setEPS()
postscript(file="MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_PLOTTED.eps", width=8.5, height=150)#, onefile=TRUE)
#pdf(file="MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_PLOTTED.pdf", width=8.5, height=150, onefile=TRUE)

#quartz(width=8.5, height=130)
#plot(mamPhy, cex=0.25, label.offset=0.4, y.lim=c(0,5690), x.lim=c(-28.47793, 337.00159))
plot(mamPhy, cex=0.25, label.offset=0.4, y.lim=c(0,3960),x.lim=c(-35, 220))#, x.lim=c(-28.47793, 337.00159))

mtext(text="Fig. S9", side=3, line=2, font=2, cex=0.75, adj=0)
mtext(text="Node-dating exponential backbone -- Upham et al. 2019", side=3, line=0, font=2, cex=0.75, adj=0)
mtext(text="MCC tree of DNA-only (4098 species) TopoFree distribution of 10,000 trees", side=3, line=-1, font=1, cex=0.75, adj=0)
mtext(text="Mean node ages; Bayesian PP >= 0.95 (black), < 0.95 (red)", side=3, line=-2, font=1, cex=0.75, adj=0)

#nodelabels(cex=0.2)
#HPDbars(mamPhy, label="height_95%_HPD", broken=T, lwd=1.5, col=hsv(0.65,1,1,alpha=0.7))
HPDbars(mamPhy, label="height_95%_HPD", broken=T, lwd=1.5, col=hsv(0.65,1,1))
node.support(mamPhy$posterior, mode="dots", col = "red", cex=0.18)
node.support(mamPhy$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.2)
#node.support(branching.times(mamPhy), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.15)
node.support(mamPhy$height, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.17)
node.support(mamPhy$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.17)

data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-18.5, mgp=c(0,0.2,0))

dev.off()


# For plotting the FBD full tree 
# =================
library(ape)
library(phyloch)
library(phytools)
library(paleotree)
library(strap)
#library(OutbreakTools)

setwd("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/11_makingFullPosteriors/10k_tree_andDR_distributions_FINAL")

FBD0<-read.beast("MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_MCC_v2_target.tre")
FBD1<-drop.tip2(FBD0, "_Anolis_carolinensis")
FBD<-ladderize(FBD1)
#FBD<-mamPhy

###
t<-FBD
num_taxa <- length(t$tip.label)
t$root.time <- max(t$height)

stem_length <- 0
origin_HPD<-c(t$"height_95%_HPD_MIN"[1],t$"height_95%_HPD_MAX"[1])
names(origin_HPD)<-c("lower","upper")

display_all_node_bars <- TRUE

names_list <-vector()
for (name in t$tip){
  v <- strsplit(name,"_")[[1]]
  if(display_all_node_bars){
  	names_list = c(names_list,name)
  }
  else if(v[length(v)]=="0"){
  	names_list = c(names_list,name)
  }
}

nids <- vector()
pos <- 1
len_nl <- length(names_list)
for(n in names_list){
  for(nn in names_list[pos:len_nl]){
    if(n != nn){
      m <- getMRCA(t,c(n,nn))
      if(m %in% nids == FALSE){
        nids <- c(nids,m)
      }
    }
  }
  pos<-pos+1
}

root_max <- t$"height_95%_HPD_MAX"[1]
x_max <- origin_HPD[2] * 0.02 + origin_HPD[2]

#stratRanges<-read.table("tipStratiRanges.txt")
stratRanges<-read.table("tipStratiRanges_all4098Species.txt")
colnames(stratRanges)<-c("tip","min","max")
ages<-cbind(stratRanges$min,stratRanges$max)
rownames(ages)<-stratRanges$tip
colnames(ages)<-c("FAD","LAD")

#
setEPS()
postscript(file="MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_MCC_v2_PLOTTED.eps", width=8.5, height=150)#, onefile=TRUE)
#pdf(file="MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_MCC_v2_PLOTTED.pdf", width=8.5, height=150, onefile=TRUE)

plot(t, cex=0.25, label.offset=0.4, y.lim=c(0,4030),x.lim=c(0, 400.04995))#, x.lim=c(-28.47793, 337.00159))

mtext(text="Fig. S10", side=3, line=2, font=2, cex=0.75, adj=0)
mtext(text="Fossilized birth-death backbone (as Zhou et al. 2013) -- Upham et al. 2019", side=3, line=0, font=2, cex=0.75, adj=0)
mtext(text="MCC tree of DNA-only (4098 species + 76 fossils) TopoFree distribution of 10,000 trees", side=3, line=-1, font=1, cex=0.75, adj=0)
mtext(text="Mean node ages; Bayesian PP >= 0.95 (black), < 0.95 (red)", side=3, line=-2, font=1, cex=0.75, adj=0)
#dev.off()

lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)

for(nv in nids){
  bar_xx_a <- c(lastPP$xx[nv]+t$height[nv-num_taxa]-t$"height_95%_HPD_MIN"[nv-num_taxa], lastPP$xx[nv]-(t$"height_95%_HPD_MAX"[nv-num_taxa]-t$height[nv-num_taxa]))
#  lines(bar_xx_a,c(lastPP$yy[nv],lastPP$yy[nv]),col=rgb(0,0,1,alpha=0.6),lwd=2)
  lines(bar_xx_a,c(lastPP$yy[nv],lastPP$yy[nv]),col=rgb(0,0,1),lwd=2)
}

node.support(t$posterior, mode="dots", col = "red", cex=0.18)
node.support(t$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.2)
node.support(t$height, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.17)
node.support(t$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.17)

#t$node.label<-t$posterior
#p <- character(length(t$node.label))
#p[t$node.label >= 0.95] <- "black"
#p[t$node.label < 0.95 & t$node.label >= 0.75] <- "gray"
#p[t$node.label < 0.75] <- "white"
#nodelabels(pch=21, cex=0.3, bg=p, lwd=0.2)

data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.5, pos=-18.5, mgp=c(0,0.2,0))


dev.off()









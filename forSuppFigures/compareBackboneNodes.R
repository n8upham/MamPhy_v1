library(ape)
library(phyloch)
library(phytools)
library(paleotree)
library(strap)
library(OutbreakTools)

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/11_makingFullPosteriors")

FBD0<-read.beast("postburn10k_Bbone_FBD_60taxa_BDallDNA.inMYR_MCC_treeAnOut.tre")
FBD1<-drop.tip2(FBD0, "_Anolis_carolinensis")
FBD<-ladderize(FBD1)

stratRanges<-read.table("tipStratiRanges.txt")
colnames(stratRanges)<-c("tip","min","max")

ages<-cbind(stratRanges$min,stratRanges$max)
rownames(ages)<-stratRanges$tip
colnames(ages)<-c("FAD","LAD")

FBD$root.time <- max(FBD$height)
geoscalePhylo(tree=FBD,ages=ages,ranges=TRUE,cex.tip=0.8,cex.ts=0.55,cex.age=0.5,width=2)

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
x_max <- origin_HPD[2] * 0.1 + origin_HPD[2]

stratRanges<-read.table("tipStratiRanges.txt")
colnames(stratRanges)<-c("tip","min","max")
ages<-cbind(stratRanges$min,stratRanges$max)
rownames(ages)<-stratRanges$tip
colnames(ages)<-c("FAD","LAD")

#
pdf("geoscaled_Mammalia.pdf", width=20, height=30)
geoscalePhylo(tree=ladderize(t,right=FALSE), units=c("Era","Period", "Epoch"), boxes="Epoch", cex.tip=0.9,cex.age=2, cex.ts=2,label.offset=0,tick.scale=50,x.lim=c(-150,x_max),lwd=4, width=2,quat.rm=TRUE)#, ages=ages,ranges=TRUE)

lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)

for(nv in nids){
  bar_xx_a <- c(lastPP$xx[nv]+t$height[nv-num_taxa]-t$"height_95%_HPD_MIN"[nv-num_taxa], lastPP$xx[nv]-(t$"height_95%_HPD_MAX"[nv-num_taxa]-t$height[nv-num_taxa]))
  lines(bar_xx_a,c(lastPP$yy[nv],lastPP$yy[nv]),col=rgb(0,0,1,alpha=0.6),lwd=4)
}

t$node.label<-t$posterior
p <- character(length(t$node.label))
p[t$node.label >= 0.95] <- "black"
p[t$node.label < 0.95 & t$node.label >= 0.75] <- "gray"
p[t$node.label < 0.75] <- "white"
nodelabels(pch=21, cex=1.5, bg=p)
node.support(t$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.9)

dev.off()

##
pdf("geoscaled_Mammalia_smallNoTip.pdf", width=7.5, height=5)
geoscalePhylo(tree=ladderize(t,right=FALSE), units=c("Era","Period", "Epoch"), direction= "upwards", boxes="Epoch", cex.tip=0.9,cex.age=0.7, cex.ts=0.7,label.offset=0,tick.scale=50,x.lim=c(-0,325),lwd=3, width=2,quat.rm=TRUE, show.tip.label=FALSE)#, ages=ages,ranges=TRUE)

lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)

#for(nv in nids){
#  bar_xx_a <- c(lastPP$xx[nv]+t$height[nv-num_taxa]-t$"height_95%_HPD_MIN"[nv-num_taxa], lastPP$xx[nv]-(t$"height_95%_HPD_MAX"[nv-#num_taxa]-t$height[nv-num_taxa]))
#  lines(bar_xx_a,c(lastPP$yy[nv],lastPP$yy[nv]),col=rgb(0,0,1,alpha=0.6),lwd=4)
#}
dev.off()


#####
quartz(width=8.5, height=22)
geoscalePhylo(tree=FBD, ages=ages, units=c("Era","Period", "Epoch"), direction="rightwards", boxes="Epoch", x.lim=c(-150,FBD$root.time+30), quat.rm=TRUE, cex.tip=0.5, cex.age=0.7, cex.ts=0.7, tick.scale=50, label.offset=0, lwd=3, width=2)
#HPDbars(FBD, label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))
node.support(FBD$posterior, mode="dots", col = "red", cex=0.5)
node.support(FBD$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.51)
node.support(FBD$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.5)




#########
FBD0<-read.beast("postburn10k_Bbone_FBD_60taxa_BDallDNA.inMYR_MCC_treeAnOut.tre")
FBD_toDrop<-c("_Anolis_carolinensis","X_Adelobasileus", "X_Aegialodon", "X_Akidolestes", "X_Albertatherium", "X_Ambondro", "X_Amphilestes", "X_Amphitherium", "X_Anchistodelphys", "X_Andinodelphys", "X_Asfaltomylos", "X_Asiatherium", "X_Asioryctes", "X_Atokatheridium", "X_Ausktribosphenos", "X_Bishops", "X_Castorocauda", "X_Cimolodontidae", "X_Daulestes", "X_Deltatheridium", "X_Didelphodon", "X_Dryolestes", "X_Eomaia", "X_Fruitafossor", "X_Gobiconodon", "X_Hadrocodium", "X_Haldanodon", "X_Haramiyavia", "X_Henkelotherium", "X_Holoclemensia", "X_Jeholodens", "X_Juramaia", "X_Kennalestes", "X_Kielantherium", "X_Kokopellia", "X_Kuehneodon", "X_Maotherium", "X_Massetognathus", "X_Mayulestes", "X_Megaconus", "X_Megazostrodon", "X_Montanalestes", "X_Morganucodon", "X_Murtoilestes", "X_Nanolestes", "X_Obdurodon", "X_Pachygenelus", "X_Pediomys", "X_Peramus", "X_Plagiaulacidae", "X_Priacodon", "X_Probainognathus", "X_Prokennalestes", "X_Pseudotribos", "X_Pucadelphys", "X_Repenomamus", "X_Rugosodon", "X_Shuotherium", "X_Sineleutherus", "X_Sinobaatar", "X_Sinoconodon", "X_Sinodelphys", "X_Spalacotherium", "X_Steropodon", "X_Sulestes", "X_Teinolophos", "X_Thomasia", "X_Thrinaxodon", "X_Tinodon", "X_Trioracodon", "X_Tritylodontidae", "X_Turgidodon", "X_Ukhaatherium", "X_Vincelestes", "X_Yanoconodon", "X_Zalambdalestes", "X_Zhangheotherium")
FBD1<-drop.tip2(FBD0, FBD_toDrop)
FBD<-ladderize(FBD1)

pdf(file="postburn10k_Bbone_FBD_59taxa_BDallDNA_agesPP.pdf", width=11, height=17, onefile=TRUE)

plot(FBD, cex=0.7, label.offset=0.4, x.lim=c(-27.21775, 455.10395))

HPDbars(FBD, label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))
node.support(FBD$posterior, mode="dots", col = "red", cex=0.7)
node.support(FBD$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.71)
node.support(FBD$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.7)
node.support(FBD$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.7)

data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period", "era"), cex = 1, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.8, pos=-2.5, mgp=c(0,0.3,0))

mtext("FBD, 60 taxa", side=3, outer=TRUE, cex=1)

dev.off()

#####

NDexp0<-read.beast("postburn10k_Bbone_ND_60taxa_17calExp_BDallDNA.inMYR_MCC_treeAnOut.tre")
NDexp1<-drop.tip2(NDexp0,"_Anolis_carolinensis")
NDexp<-ladderize(NDexp1)

pdf(file="postburn10k_Bbone_NDexp_59taxa_BDallDNA_agesPP.pdf", width=11, height=17, onefile=TRUE)

plot(NDexp, cex=0.7, label.offset=0.4, x.lim=c(-23.58487, 334.2138))

HPDbars(NDexp, label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))
node.support(NDexp$posterior, mode="dots", col = "red", cex=0.7)
node.support(NDexp$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.71)
node.support(NDexp$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.7)
node.support(NDexp$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.7)

data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period", "era"), cex = 1, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.8, pos=-2.5, mgp=c(0,0.3,0))

dev.off()


#####

NDuni0<-read.beast("postburn10k_Bbone_ND_60taxa_17calUni_BDallDNA.inMYR_MCC_treeAnOut.tre")
NDuni1<-drop.tip2(NDuni0,"_Anolis_carolinensis")
NDuni<-ladderize(NDuni1)

pdf(file="postburn10k_Bbone_NDuni_59taxa_BDallDNA_agesPP.pdf", width=11, height=17, onefile=TRUE)

plot(NDuni, cex=0.7, label.offset=0.4, x.lim=c(-23.58487, 334.2138))

HPDbars(NDuni, label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))
node.support(NDuni$posterior, mode="dots", col = "red", cex=0.7)
node.support(NDuni$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.71)
node.support(NDuni$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.7)
node.support(NDuni$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.7)

data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period", "era"), cex = 1, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.8, pos=-2.5, mgp=c(0,0.3,0))

dev.off()

#######
# Just get the trees
# 59 TAXA - comparable
#####

FBD0<-read.beast("postburn10k_Bbone_FBD_60taxa_BDallDNA.inMYR_MCC_treeAnOut.tre")
FBD_toDrop<-c("_Anolis_carolinensis","X_Adelobasileus", "X_Aegialodon", "X_Akidolestes", "X_Albertatherium", "X_Ambondro", "X_Amphilestes", "X_Amphitherium", "X_Anchistodelphys", "X_Andinodelphys", "X_Asfaltomylos", "X_Asiatherium", "X_Asioryctes", "X_Atokatheridium", "X_Ausktribosphenos", "X_Bishops", "X_Castorocauda", "X_Cimolodontidae", "X_Daulestes", "X_Deltatheridium", "X_Didelphodon", "X_Dryolestes", "X_Eomaia", "X_Fruitafossor", "X_Gobiconodon", "X_Hadrocodium", "X_Haldanodon", "X_Haramiyavia", "X_Henkelotherium", "X_Holoclemensia", "X_Jeholodens", "X_Juramaia", "X_Kennalestes", "X_Kielantherium", "X_Kokopellia", "X_Kuehneodon", "X_Maotherium", "X_Massetognathus", "X_Mayulestes", "X_Megaconus", "X_Megazostrodon", "X_Montanalestes", "X_Morganucodon", "X_Murtoilestes", "X_Nanolestes", "X_Obdurodon", "X_Pachygenelus", "X_Pediomys", "X_Peramus", "X_Plagiaulacidae", "X_Priacodon", "X_Probainognathus", "X_Prokennalestes", "X_Pseudotribos", "X_Pucadelphys", "X_Repenomamus", "X_Rugosodon", "X_Shuotherium", "X_Sineleutherus", "X_Sinobaatar", "X_Sinoconodon", "X_Sinodelphys", "X_Spalacotherium", "X_Steropodon", "X_Sulestes", "X_Teinolophos", "X_Thomasia", "X_Thrinaxodon", "X_Tinodon", "X_Trioracodon", "X_Tritylodontidae", "X_Turgidodon", "X_Ukhaatherium", "X_Vincelestes", "X_Yanoconodon", "X_Zalambdalestes", "X_Zhangheotherium")
FBD1<-drop.tip2(FBD0, FBD_toDrop)
FBD<-ladderize(FBD1)

taxNames<-read.table("bbone_taxNames.txt")
colnames(taxNames)<-c("tip","sp","spFam")

FBD$tip.label<-as.vector(taxNames[match(FBD$tip.label, as.vector(taxNames$tip)),"spFam"])

NDexp0<-read.beast("postburn10k_Bbone_ND_60taxa_17calExp_BDallDNA.inMYR_MCC_treeAnOut.tre")
NDexp1<-drop.tip2(NDexp0,"_Anolis_carolinensis")
NDexp<-ladderize(NDexp1)
NDexp$tip.label<-as.vector(taxNames[match(NDexp$tip.label, as.vector(taxNames$tip)),"spFam"])

NDuni0<-read.beast("postburn10k_Bbone_ND_60taxa_17calUni_BDallDNA.inMYR_MCC_treeAnOut.tre")
NDuni1<-drop.tip2(NDuni0,"_Anolis_carolinensis")
NDuni<-ladderize(NDuni1)
NDuni$tip.label<-as.vector(taxNames[match(NDuni$tip.label, as.vector(taxNames$tip)),"spFam"])

bbones <- c(FBD, NDexp, NDuni)

##
# Get the TABLE values for this.

plot(FBD, cex=0.6)
nodelabels(cex=0.6)
plot(NDexp, cex=0.6)
nodelabels(cex=0.6)
plot(NDuni, cex=0.6)
nodelabels(cex=0.6)

#>> nodes are not exactly identical...
#all.equal.phylo(FBD,NDexp, use.edge.length=FALSE, use.tip.label=FALSE, index.return=TRUE)

FBD_NDexp <- matchNodes(FBD,NDexp, "descendants")

FBD_NDuni <- matchNodes(FBD,NDuni, "descendants")

NDexp_NDuni <- matchNodes(NDexp,NDuni, "descendants")

# Get node data:
FBD_table<-cbind(FBD$"height_95%_HPD_MIN", FBD$"height_95%_HPD_MAX",FBD$height,FBD$height_median)
colnames(FBD_table)<-c("height_95%_HPD_MIN","height_95%_HPD_MAX","height_mean","height_median")
rownames(FBD_table)<-FBD_NDexp[,1]

NDexp_table<-cbind(NDexp$"height_95%_HPD_MIN", NDexp$"height_95%_HPD_MAX",NDexp$height,NDexp$height_median)
colnames(NDexp_table)<-c("height_95%_HPD_MIN","height_95%_HPD_MAX","height_mean","height_median")
rownames(NDexp_table)<-FBD_NDexp[,1] # this is just 60:117 continuous

NDuni_table<-cbind(NDuni$"height_95%_HPD_MIN", NDuni$"height_95%_HPD_MAX",NDuni$height,NDuni$height_median)
colnames(NDuni_table)<-c("height_95%_HPD_MIN","height_95%_HPD_MAX","height_mean","height_median")
rownames(NDuni_table)<-FBD_NDexp[,1] # this is just 60:117 continuous

mins <- cbind(FBD_table[,"height_95%_HPD_MIN"],NDexp_table[match(FBD_NDexp[,2],rownames(NDexp_table)),"height_95%_HPD_MIN"],NDuni_table[match(FBD_NDuni[,2],rownames(NDuni_table)),"height_95%_HPD_MIN"])
colnames(mins)<-c("FBD_min","NDexp_min","NDuni_min")

maxs <- cbind(FBD_table[,"height_95%_HPD_MAX"],postburn10k_Bbone_NDexp_59taxa_BDallDNA_agesPPle[match(FBD_NDexp[,2],rownames(NDexp_table)),"height_95%_HPD_MAX"],NDuni_table[match(FBD_NDuni[,2],rownames(NDuni_table)),"height_95%_HPD_MAX"])
colnames(maxs)<-c("FBD2_max","NDexp_max","NDuni_max")

means <- cbind(FBD_table[,"height_mean"],NDexp_table[match(FBD_NDexp[,2],rownames(NDexp_table)),"height_mean"],NDuni_table[match(FBD_NDuni[,2],rownames(NDuni_table)),"height_mean"])
colnames(means)<-c("FBD_mean","NDexp_mean","NDuni_mean")

medians <- cbind(FBD_table[,"height_median"],NDexp_table[match(FBD_NDexp[,2],rownames(NDexp_table)),"height_median"],NDuni_table[match(FBD_NDuni[,2],rownames(NDuni_table)),"height_median"])
colnames(medians)<-c("FBD_med","NDexp_med","NDuni_med")

summary<-cbind(mins,maxs,means,medians)

write.table(summary,"bbone_nodeCompare_FBD-NDexp-NDuni.txt")


# Now plot those 3 trees above each other...
#m <- matrix(c(1, 2, 3), nrow = 3, ncol = 1)
#
#pdf(file="Bbone_normal_FBD-NDexp-NDuni-59taxa_agesPP.pdf", width=8.5, height=22, onefile=TRUE)
#layout(m)

# Plot full 59 TAXA
# NORMAL

pdf(file="Bbone_normal_FBD-59taxa_agesPP.pdf", width=8.5, height=11, onefile=TRUE)
plot(bbones[[1]], cex=0.9, label.offset=0.4, no.margin=FALSE, x.lim=c(-27.21775, 449.1694))
HPDbars(bbones[[1]], label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))
node.support(bbones[[1]]$posterior, mode="dots", col = "red", cex=0.7)
node.support(bbones[[1]]$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.71)
node.support(bbones[[1]]$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.8)
node.support(bbones[[1]]$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.7)
data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period", "era"), cex = 1, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.8, pos=-4.5, mgp=c(0,0.3,0))
dev.off()

pdf(file="Bbone_normal_NDexp-59taxa_agesPP.pdf", width=8.5, height=11, onefile=TRUE)
plot(bbones[[2]], cex=0.9, label.offset=0.4, no.margin=FALSE, x.lim=c(-23.58487, 329.90706))
HPDbars(bbones[[2]], label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))
node.support(bbones[[2]]$posterior, mode="dots", col = "red", cex=0.7)
node.support(bbones[[2]]$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.71)
node.support(bbones[[2]]$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.8)
node.support(bbones[[2]]$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.7)
data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period", "era"), cex = 1, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.8, pos=-4.5, mgp=c(0,0.3,0))
dev.off()

pdf(file="Bbone_normal_NDuni-59taxa_agesPP.pdf", width=8.5, height=11, onefile=TRUE)
plot(bbones[[3]], cex=0.9, label.offset=0.4, no.margin=FALSE, x.lim=c(-8.78774, 323.87127))
HPDbars(bbones[[3]], label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))
node.support(bbones[[3]]$posterior, mode="dots", col = "red", cex=0.7)
node.support(bbones[[3]]$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.71)
node.support(bbones[[3]]$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.8)
node.support(bbones[[3]]$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.7)
data(gradstein04)
axisGeo(GTS = gradstein04, unit = c("epoch","period", "era"), cex = 1, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.8, pos=-4.5, mgp=c(0,0.3,0))
dev.off()

# GEOSCALE for the 59 TAXA

pdf(file="Bbone_geoscaleR_FBD-59taxa_agesPP.pdf", width=8.5, height=11, onefile=TRUE)
FBD$root.time <- max(FBD$height)
geoscalePhylo(tree=FBD, units=c("Era","Period", "Epoch"), direction="rightwards", boxes="Epoch", x.lim=c(-150,FBD$root.time+30), quat.rm=TRUE, cex.tip=0.7, cex.age=0.7, cex.ts=0.7, tick.scale=50, label.offset=0.4, lwd=3, width=2)
HPDbars(FBD, label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))
node.support(FBD$posterior, mode="dots", col = "red", cex=0.7)
node.support(FBD$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.71)
node.support(FBD$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.7)
dev.off()

pdf(file="Bbone_geoscaleR_NDexp-59taxa_agesPP.pdf", width=8.5, height=11, onefile=TRUE)
NDexp$root.time <- max(NDexp$height)
geoscalePhylo(tree=NDexp, units=c("Era","Period", "Epoch"), direction="rightwards", boxes="Epoch", x.lim=c(-150,NDexp$root.time+30), quat.rm=TRUE, cex.tip=0.7, cex.age=0.7, cex.ts=0.7, tick.scale=50, label.offset=0.4, lwd=3, width=2)
HPDbars(NDexp, label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))
node.support(NDexp$posterior, mode="dots", col = "red", cex=0.7)
node.support(NDexp$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.71)
node.support(NDexp$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.7)
dev.off()

pdf(file="Bbone_geoscaleR_NDuni-59taxa_agesPP.pdf", width=8.5, height=11, onefile=TRUE)
NDuni$root.time <- max(NDuni$height)
geoscalePhylo(tree=NDuni, units=c("Era","Period", "Epoch"), direction="rightwards", boxes="Epoch", x.lim=c(-150,NDuni$root.time+30), quat.rm=TRUE, cex.tip=0.7, cex.age=0.7, cex.ts=0.7, tick.scale=50, label.offset=0.4, lwd=3, width=2)
HPDbars(NDuni, label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))
node.support(NDuni$posterior, mode="dots", col = "red", cex=0.7)
node.support(NDuni$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.71)
node.support(NDuni$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.7)
dev.off()



####
# Then drop all the non-PatchClade reps-- giving just the REAL backbones...
# 28 TAXA
####

toKeep <- c("Anomalurus_beecrofti__ANOMALURIDAE__RODENTIA", "Bison_bison__BOVIDAE__CETARTIODACTYLA", "Callithrix_jacchus__CALLITRICHIDAE__PRIMATES", "Calomyscus_baluchi__CALOMYSCIDAE__RODENTIA", "Castor_canadensis__CASTORIDAE__RODENTIA", "Cricetomys_gambianus__NESOMYIDAE__RODENTIA", "Cricetulus_barabensis__CRICETIDAE__RODENTIA", "Dasypus_novemcinctus__DASYPODIDAE__CINGULATA", "Eptesicus_fuscus__VESPERTILIONIDAE__CHIROPTERA", "Equus_caballus__EQUIDAE__PERISSODACTYLA", "Erethizon_dorsatum__ERETHIZONTIDAE__RODENTIA", "Erinaceus_europaeus__ERINACEIDAE__EULIPOTYPHLA", "Felis_catus__FELIDAE__CARNIVORA", "Galeopterus_variegatus__CYNOCEPHALIDAE__DERMOPTERA", "Ictidomys_tridecemlineatus__SCIURIDAE__RODENTIA", "Jaculus_jaculus__DIPODIDAE__RODENTIA", "Manis_pentadactyla__MANIDAE__PHOLIDOTA", "Ornithorhynchus_anatinus__ORNITHORHYNCHIDAE__MONOTREMATA", "Oryctolagus_cuniculus__LEPORIDAE__LAGOMORPHA", "Pteronotus_parnellii__MORMOOPIDAE__CHIROPTERA", "Pteropus_alecto__PTEROPODIDAE__CHIROPTERA", "Rattus_norvegicus__MURIDAE__RODENTIA", "Saccopteryx_bilineata__EMBALLONURIDAE__CHIROPTERA", "Spalax_ehrenbergi__SPALACIDAE__RODENTIA", "Trichechus_manatus__TRICHECHIDAE__SIRENIA", "Tupaia_belangeri__TUPAIIDAE__SCANDENTIA", "Typhlomys_cinereus__PLATACANTHOMYIDAE__RODENTIA", "Vombatus_ursinus__VOMBATIDAE__DIPROTODONTIA")
toDropFBD <- setdiff(FBD0$tip.label,toKeep)
toDropND <- setdiff(NDexp0$tip.label,toKeep)

repFamOrdSp<-read.table("PCreps_toFamOrd_toSp_toSpFam.txt")
names(repFamOrdSp)<-c("tip","famOrd","sp","spFam")

FBD0<-read.beast("postburn10k_Bbone_FBD_60taxa_BDallDNA.inMYR_MCC_treeAnOut.tre")
FBD1<-drop.tip2(FBD0, toDropFBD)
FBD<-ladderize(FBD1)
FBD$tip.label<-as.vector(repFamOrdSp[match(FBD$tip.label, as.vector(repFamOrdSp$tip)),"sp"])

NDexp0<-read.beast("postburn10k_Bbone_ND_60taxa_17calExp_BDallDNA.inMYR_MCC_treeAnOut.tre")
NDexp1<-drop.tip2(NDexp0,toDropND)
NDexp<-ladderize(NDexp1)
NDexp$tip.label<-as.vector(repFamOrdSp[match(NDexp$tip.label, as.vector(repFamOrdSp$tip)),"sp"])

NDuni0<-read.beast("postburn10k_Bbone_ND_60taxa_17calUni_BDallDNA.inMYR_MCC_treeAnOut.tre")
NDuni1<-drop.tip2(NDuni0,toDropND)
NDuni<-ladderize(NDuni1)
NDuni$tip.label<-as.vector(repFamOrdSp[match(NDuni$tip.label, as.vector(repFamOrdSp$tip)),"sp"])

bbonesPC <- c(FBD, NDexp, NDuni)

# GEOSCALE, do each separate, then combine...

pdf(file="Bbone_geoscaleR_FBD-28taxa_agesPP_sp.pdf", width=8.5, height=7, onefile=TRUE)
FBD$root.time <- max(FBD$height)
geoscalePhylo(tree=FBD, units=c("Era","Period", "Epoch"), direction="rightwards", boxes="Epoch", x.lim=c(-50,FBD$root.time+30), quat.rm=TRUE, cex.tip=0.7, cex.age=0.7, cex.ts=0.7, tick.scale=50, label.offset=0.4, lwd=3, width=2)
HPDbars(FBD, label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))
node.support(FBD$posterior, mode="dots", col = "red", cex=0.7)
node.support(FBD$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.71)
node.support(FBD$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.7)
node.support(FBD$posterior, mode="numbers", digits=2,pos="right", col = "grey20", font=2, cex=0.7)
dev.off()

pdf(file="Bbone_geoscaleR_NDexp-28taxa_agesPP_sp.pdf", width=8.5, height=7, onefile=TRUE)
NDexp$root.time <- max(NDexp$height)
geoscalePhylo(tree=NDexp, units=c("Era","Period", "Epoch"), direction="rightwards", boxes="Epoch", x.lim=c(-50,NDexp$root.time+30), quat.rm=TRUE, cex.tip=0.7, cex.age=0.7, cex.ts=0.7, tick.scale=50, label.offset=0.4, lwd=3, width=2)
HPDbars(NDexp, label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))
node.support(NDexp$posterior, mode="dots", col = "red", cex=0.7)
node.support(NDexp$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.71)
node.support(NDexp$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.7)
node.support(NDexp$posterior, mode="numbers", digits=2,pos="right", col = "grey20", font=2, cex=0.7)
dev.off()

pdf(file="Bbone_geoscaleR_NDuni-28taxa_agesPP_sp.pdf", width=8.5, height=7, onefile=TRUE)
NDuni$root.time <- max(NDuni$height)
geoscalePhylo(tree=NDuni, units=c("Era","Period", "Epoch"), direction="rightwards", boxes="Epoch", x.lim=c(-50,NDuni$root.time+30), quat.rm=TRUE, cex.tip=0.7, cex.age=0.7, cex.ts=0.7, tick.scale=50, label.offset=0.4, lwd=3, width=2)
HPDbars(NDuni, label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))
node.support(NDuni$posterior, mode="dots", col = "red", cex=0.7)
node.support(NDuni$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.71)
node.support(NDuni$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.7)
node.support(NDuni$posterior, mode="numbers", digits=2,pos="right", col = "grey20", font=2, cex=0.7)
dev.off()

#########

#####
# DENSITREES, backbone levels
library(ape)
library(phangorn)
library(phytools)
library(phyloch)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
##
# Load in the 10k trees of the BACKBONE
NDexp10k0<-read.nexus("postburn10k_Bbone_ND_60taxa_17calExp_BDallDNA.inMYR.trees")
# prune to 28 taxa, ladderize
toKeep <- c("Anomalurus_beecrofti__ANOMALURIDAE__RODENTIA", "Bison_bison__BOVIDAE__CETARTIODACTYLA", "Callithrix_jacchus__CALLITRICHIDAE__PRIMATES", "Calomyscus_baluchi__CALOMYSCIDAE__RODENTIA", "Castor_canadensis__CASTORIDAE__RODENTIA", "Cricetomys_gambianus__NESOMYIDAE__RODENTIA", "Cricetulus_barabensis__CRICETIDAE__RODENTIA", "Dasypus_novemcinctus__DASYPODIDAE__CINGULATA", "Eptesicus_fuscus__VESPERTILIONIDAE__CHIROPTERA", "Equus_caballus__EQUIDAE__PERISSODACTYLA", "Erethizon_dorsatum__ERETHIZONTIDAE__RODENTIA", "Erinaceus_europaeus__ERINACEIDAE__EULIPOTYPHLA", "Felis_catus__FELIDAE__CARNIVORA", "Galeopterus_variegatus__CYNOCEPHALIDAE__DERMOPTERA", "Ictidomys_tridecemlineatus__SCIURIDAE__RODENTIA", "Jaculus_jaculus__DIPODIDAE__RODENTIA", "Manis_pentadactyla__MANIDAE__PHOLIDOTA", "Ornithorhynchus_anatinus__ORNITHORHYNCHIDAE__MONOTREMATA", "Oryctolagus_cuniculus__LEPORIDAE__LAGOMORPHA", "Pteronotus_parnellii__MORMOOPIDAE__CHIROPTERA", "Pteropus_alecto__PTEROPODIDAE__CHIROPTERA", "Rattus_norvegicus__MURIDAE__RODENTIA", "Saccopteryx_bilineata__EMBALLONURIDAE__CHIROPTERA", "Spalax_ehrenbergi__SPALACIDAE__RODENTIA", "Trichechus_manatus__TRICHECHIDAE__SIRENIA", "Tupaia_belangeri__TUPAIIDAE__SCANDENTIA", "Typhlomys_cinereus__PLATACANTHOMYIDAE__RODENTIA", "Vombatus_ursinus__VOMBATIDAE__DIPROTODONTIA")
toDropND <- setdiff(NDexp10k0[[1]]$tip.label,toKeep)
NDexp10k1<-lapply(NDexp10k0,drop.tip,tip=toDropND)
NDexp10k<-lapply(NDexp10k1,ladderize)
class(NDexp10k)<-"multiPhylo"

# MCC tree
NDexp0<-read.beast("postburn10k_Bbone_ND_60taxa_17calExp_BDallDNA.inMYR_MCC_treeAnOut.tre")
NDexp1<-drop.tip2(NDexp0,toDropND)
NDexp<-ladderize(NDexp1)

repPC<-read.table("PCreps_toPC.txt")
names(repPC)<-c("tip","PC")

NDexp100<-sample(NDexp10k,size=100)
NDexp500<-sample(NDexp10k,size=500)
NDexp1k<-sample(NDexp10k,size=1000)

for (i in 1:length(NDexp500)){
  NDexp500[[i]]$tip.label<-as.vector(repPC[match(NDexp500[[i]]$tip.label, as.vector(repPC$tip)),"PC"])
}
NDexp$tip.label<-as.vector(repPC[match(NDexp$tip.label, as.vector(repPC$tip)),"PC"])
for (i in 1:length(NDexp1k)){
  NDexp1k[[i]]$tip.label<-as.vector(repPC[match(NDexp1k[[i]]$tip.label, as.vector(repPC$tip)),"PC"])
}
for (i in 1:length(NDexp10k)){
  NDexp10k[[i]]$tip.label<-as.vector(repPC[match(NDexp10k[[i]]$tip.label, as.vector(repPC$tip)),"PC"])
}
# 500trees: find which trees are A(X+B) vs (A+X)B
nodeSizes<-c()
for(i in 1:length(NDexp500)){
  tree<-NDexp500[[i]]
  node<-getMRCA(tree,tip=c("Afrotheria","Xenarthra"))
  nodeTips<-tree$tip.label[Descendants(tree,node)[[1]]]  
  nodeSizes[i]<-length(nodeTips)
}
# nodeSizes
#  2  26 
# 270 230
treeNum<-c(1:length(NDexp500))
ID<-cbind.data.frame(treeNum,nodeSizes)
treesAX<-ID[which(ID$nodeSizes==2),]
treesXB<-ID[which(ID$nodeSizes==26),]

# 1000 trees: find which trees are A(X+B) vs (A+X)B
nodeSizes<-c()
for(i in 1:length(NDexp1k)){
  tree<-NDexp1k[[i]]
  node<-getMRCA(tree,tip=c("Afrotheria","Xenarthra"))
  nodeTips<-tree$tip.label[Descendants(tree,node)[[1]]]  
  nodeSizes[i]<-length(nodeTips)
}
# nodeSizes
#  2  26 
# 541 459 
treeNum<-c(1:length(NDexp1k))
ID<-cbind.data.frame(treeNum,nodeSizes)
treesAX<-ID[which(ID$nodeSizes==2),]
treesXB<-ID[which(ID$nodeSizes==26),]

# 10k trees: find which trees are A(X+B) vs (A+X)B
nodeSizes<-c()
for(i in 1:length(NDexp10k)){
  tree<-NDexp10k[[i]]
  node<-getMRCA(tree,tip=c("Afrotheria","Xenarthra"))
  nodeTips<-tree$tip.label[Descendants(tree,node)[[1]]]  
  nodeSizes[i]<-length(nodeTips)
} 
# nodeSizes
#   2   26 
# 5314 4686
treeNum<-c(1:length(NDexp10k))
ID<-cbind.data.frame(treeNum,nodeSizes)
treesAX<-ID[which(ID$nodeSizes==2),]
treesXB<-ID[which(ID$nodeSizes==26),]


# with densiTree (phangorn)::
# red-blue together
pdf(file="postburn10k_Bbone_ND_60taxa_17calExp_BDallDNA.inMYR_DensityTree_1k_a01_PCs_axis_blueAX-redXB-Together_densi.pdf",width=7, height=10)
#quartz(width=8,height=8)
#par(mfcol=c(1,2))
#densiTree(NDexpM, type="phylogram",col="black",alpha=0.5, consensus=NDexp, no.margin=FALSE, cex=0.7, label.offset=0.4, direction="rightwards", show.tip.label=TRUE)
densiTree(NDexp1k[treesXB[,1]], type="phylogram",col="red",alpha=0.01, consensus=NDexp, no.margin=FALSE, cex=0.9, label.offset=0.4, direction="rightwards", show.tip.label=TRUE)
par(new=TRUE)
densiTree(NDexp1k[treesAX[,1]], type="phylogram",col="blue",alpha=0.01, consensus=NDexp, no.margin=FALSE, cex=0.9, label.offset=0.4, direction="rightwards", show.tip.label=TRUE)
par(new=TRUE)
plot(NDexp, type="phylogram",ftype="reg", fsize=0.5,direction="rightwards",colors="black")
dev.off()

# blueAX only 
pdf(file="postburn10k_Bbone_ND_60taxa_17calExp_BDallDNA.inMYR_DensityTree_1k_a01_PCs_axis_blueAX-Only_densi.pdf",width=7, height=10)
#quartz(width=8,height=8)
#par(mfcol=c(1,2))
#densiTree(NDexpM, type="phylogram",col="black",alpha=0.5, consensus=NDexp, no.margin=FALSE, cex=0.7, label.offset=0.4, direction="rightwards", show.tip.label=TRUE)
densiTree(NDexp1k[treesAX[,1]], type="phylogram",col="blue",alpha=0.01, consensus=NDexp, no.margin=FALSE, cex=0.9, label.offset=0.4, direction="rightwards", show.tip.label=TRUE)
dev.off()

# redXB only 
pdf(file="postburn10k_Bbone_ND_60taxa_17calExp_BDallDNA.inMYR_DensityTree_1k_a01_PCs_axis_redXB-Only_densi.pdf",width=7, height=10)
#quartz(width=8,height=8)
#par(mfcol=c(1,2))
#densiTree(NDexpM, type="phylogram",col="black",alpha=0.5, consensus=NDexp, no.margin=FALSE, cex=0.7, label.offset=0.4, direction="rightwards", show.tip.label=TRUE)
densiTree(NDexp1k[treesXB[,1]], type="phylogram",col="red",alpha=0.01, consensus=NDexp, no.margin=FALSE, cex=0.9, label.offset=0.4, direction="rightwards", show.tip.label=TRUE)
dev.off()

# MCC only -- doesnt work to plot NODE info in the densiTree
NDexpM<-c(NDexp,NDexp)
class(NDexpM)<-"multiPhylo"

pdf(file="postburn10k_Bbone_ND_60taxa_17calExp_BDallDNA.inMYR_DensityTree_1k_a01_PCs_axis_MCC-Only_densi.pdf",width=7, height=10)
#quartz(width=8,height=8)
#par(mfcol=c(1,2))
#densiTree(NDexpM, type="phylogram",col="black",alpha=0.5, consensus=NDexp, no.margin=FALSE, cex=0.7, label.offset=0.4, direction="rightwards", show.tip.label=TRUE)
densiTree(NDexpM, type="phylogram",col="black",alpha=0.1, consensus=NDexp, no.margin=FALSE, cex=0.9, label.offset=0.4, direction="rightwards", show.tip.label=TRUE)

plot(NDexp, cex=0.9, label.offset=0.4, direction="rightwards", edge.color="white",show.tip.label=FALSE)#, x.lim=c(-23.58487, 218.42446))#, x.lim=c(-23.58487, 218.42446))

par(new=TRUE)
densiTree(NDexp1k[treesXB[,1]], type="phylogram",col="red",alpha=0.01, consensus=NDexp, no.margin=FALSE, cex=0.9, label.offset=0.4, direction="rightwards", show.tip.label=TRUE)
par(new=TRUE)
densiTree(NDexp1k[treesAX[,1]], type="phylogram",col="blue",alpha=0.01, consensus=NDexp, no.margin=FALSE, cex=0.9, label.offset=0.4, direction="rightwards", show.tip.label=TRUE)
node.support(NDexp$posterior, mode="dots", col = "red", cex=0.7)

node.support(NDexp$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.71)
#node.support(NDexp$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.7)
node.support(NDexp$posterior, mode="numbers", digits=2,pos="above", col = "black", font=2, cex=0.7)
dev.off()

#HPDbars(NDexp, label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))




# with DensityTree (phytools)::
pdf(file="postburn10k_Bbone_ND_60taxa_17calExp_BDallDNA.inMYR_DensityTree_500_a01_PCs_axis_redBlue_3col_noAgeUncert.pdf",width=10, height=5)
par(mfcol=c(1,3))
densityTree(NDexp500[treesAX[,1]], type="phylogram",ftype="reg", fsize=0.8,alpha=0.01, use.edge.length=FALSE, colors="blue")
axisPhylo(side=1)
#par(new=TRUE)
densityTree(NDexp500[treesXB[,1]], type="phylogram",ftype="reg", offset=2,fsize=0.8,alpha=0.01, use.edge.length=FALSE,direction="rightwards",colors="red")
#plot(NDexp, type="phylogram",ftype="reg", fsize=0.5,direction="rightwards",colors="black")

plot(NDexp, cex=0.7, label.offset=0.4, direction="rightwards", show.tip.label=TRUE, x.lim=c(-23.58487, 218.42446))#, x.lim=c(-23.58487, 218.42446))
HPDbars(NDexp, label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))
node.support(NDexp$posterior, mode="dots", col = "red", cex=0.7)
node.support(NDexp$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.71)
#node.support(NDexp$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.7)
node.support(NDexp$posterior, mode="numbers", digits=2,pos="above", col = "black", font=2, cex=0.7)

dev.off()




densiTree(NDexp500, type="phylogram",show.tip.label=FALSE, alpha=0.01)
densityTree(NDexp500, type="phylogram",ftype="off", alpha=0.01)


pdf(file="postburn10k_Bbone_ND_60taxa_17calExp_BDallDNA.inMYR_DensityTree_500_a01_PCs_axis.pdf")
densityTree(NDexp500, type="phylogram",ftype="reg", fsize=0.8,alpha=0.01)
axisPhylo(side=1)
dev.off()


pdf(file="postburn10k_Bbone_ND_60taxa_17calExp_BDallDNA.inMYR_density.tree_500.pdf")
density.tree(NDexp500, type="phylogram", ftype="off",use.edge.length=TRUE) # alpha=0.01, 
dev.off()


density.tree(NDexp500, type="phylogram", ftype="i",fsize=0.2,use.edge.length=TRUE) # alpha=0.01, 


make.transparent<-function(color,alpha){
    RGB<-col2rgb(color)[,1]/255
    rgb(RGB[1],RGB[2],RGB[3],alpha)
}

density.tree<-function(trees,colors="blue",alpha=NULL,method="plotTree",fix.depth=FALSE,...){
    N<-length(trees)
    if(hasArg(use.edge.length)) use.edge.length<-list(...)$use.edge.length
    else use.edge.length<-TRUE
    if(!use.edge.length) trees<-lapply(trees,compute.brlen)
    if(!fix.depth){
        h<-sapply(trees,function(x) max(nodeHeights(x)))
        ii<-order(h,decreasing=TRUE)
        trees<-trees[ii]
        h<-h[ii]
    }
    if(is.null(alpha)) alpha<-1/N
    colors<-setNames(sapply(colors,make.transparent,alpha),names(colors))
    if(method=="plotTree") foo<-plotTree else foo<-plot
    foo(trees[[1]],color=colors,...)
    xlim<-get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim
    xlim[1]<-xlim[1]+0.03703704*diff(xlim)
    xlim[2]<-xlim[2]-0.03703704*diff(xlim)
    par(fg="transparent")
    for(i in 2:length(trees))
        foo(trees[[i]],
            tips=setNames(1:Ntip(trees[[1]]),trees[[1]]$tip.label),
            color=colors,add=TRUE,
            xlim=if(fix.depth) NULL else xlim-(h[1]-h[i]),...)
    par(fg="black")
}






##
# Load in the BACKBONE MCC - for plotting with error bars...
library(phyloch)

NDexp0<-read.beast("postburn10k_Bbone_ND_60taxa_17calExp_BDallDNA.inMYR_MCC_treeAnOut.tre")
NDexp1<-drop.tip2(NDexp0,"_Anolis_carolinensis")
NDexp<-ladderize(NDexp1)

# If I wanted to add the backbone tree inside...
par(fig=c(0.53, 0.70, 0.52, 0.72),mar=c(0,0,0,0), mar=c(0,0,0,0), new=T) #c(x1, x2, y1, y2)

plot(NDexp, cex=0.7, label.offset=0.4, direction="upwards", show.tip.label=FALSE, y.lim=c(-23.58487, 218.42446))#, x.lim=c(-23.58487, 218.42446))

HPDbars(NDexp, label="height_95%_HPD", broken=T, lwd=4, col=hsv(0.65,1,1,alpha=0.7))
node.support(NDexp$posterior, mox`de="dots", col = "red", cex=0.7)
node.support(NDexp$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.71)
#node.support(NDexp$height_median, mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.7)
node.support(NDexp$posterior, mode="numbers", digits=2,pos="above", col = "black", font=2, cex=0.7)

#data(gradstein04)
#axisGeo(GTS = gradstein04, unit = c("epoch","period", "era"), cex = 1, ages=FALSE, gridty=3, gridcol="grey50")
axisPhylo(cex.axis=0.8, pos=-2.5, mgp=c(0,0.3,0))








# comparing these 53 tip and 49 tip trees

compare1<-

matchNodes(FBD,NDexp, "descendants")

all.equal.phylo(FBD,NDexp, use.edge.length=FALSE, use.tip.label=FALSE, index.return=TRUE)






#=========
# PLOT the TRAIT distributions for RATE-SHIFTED clades as compared to NON-CLADE (background) rates
#=========
library(ape); library(phytools); library(picante); library(plotrix); library(phangorn); library(phyloch)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")
#setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")

bbone<- "NDexp" # "FBD"
mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")

allTraits<-read.table("MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)
head(allTraits)

# get rateShifts
agesAllClades_orig<-read.table(file=paste("divTime_24BAMMnodes_",bbone,"_MCC_target.txt", sep=""))
nodeNums<-agesAllClades_orig$nodes[c(1,3,5,7,8,10,11,13,15,17:29,31,32)]
cladeNames<-agesAllClades_orig$ID[c(1,3,5,7,8,10,11,13,15,17:29,31,32)]

cladeTips_all<-vector("list",length(nodeNums))
for(i in 1:length(nodeNums)){
	tipNums<-getDescendants(mamMCC,nodeNums[i])
	cladeTips<-as.vector(na.omit(mamMCC$tip.label[tipNums]))
	cladeTips_all[[i]]<-cladeTips
}
cladeTips_all[[6]]<-NA # Because no species in CETACEA have geographic ranges-- AH, because focused on LAND areas...
cladeTips_all[[3]]<-NA # Is all Placentals
cladeTips_all[[1]]<-NA # is all Marsupials
cladeTips_all[[4]]<-NA # is all Carnivora, and they have BIGGER ranges... 

cladeTips_all_Combo<-unlist(cladeTips_all)

cladeTips_all_uniq<-unique(cladeTips_all_Combo) # 4205 species involved in some rate shift... 

shiftedDat<-tipDataAll[match(cladeTips_all_uniq,rownames(tipDataAll)),]


library(phytools)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"
mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")

agesAllClades_orig<-read.table(file=paste("divTime_24BAMMnodes_",bbone,"_MCC_target.txt", sep=""))
nodeNums<-agesAllClades_orig$nodes[c(1,3,5,7,8,10,11,13,15,17:29,31,32)]
IDs<-agesAllClades_orig$ID[c(1,3,5,7,8,10,11,13,15,17:29,31,32)]
CLADEs<-substr(agesAllClades_orig$CLADE[c(1,3,5,7,8,10,11,13,15,17:29,31,32)],1,10)

###
##
# for BODY MASSES
###
BM_table<-read.table(file="MamPhy_5911sp_GenSp_bodyMasses_FIN_5327species.txt", header=TRUE)

cladeBM_all<-vector("list",length(nodeNums))
for(i in 1:length(nodeNums)){
	tipNums<-getDescendants(mamMCC,nodeNums[i])
	cladeTips<-as.vector(na.omit(mamMCC$tip.label[tipNums]))
	cladeBM_i<-as.numeric(na.omit(BM_table[match(cladeTips,as.vector(BM_table$mamPhy_tiplabel)),"BM_final"]))
	cladeBM_all[[i]]<-cladeBM_i
}












#=========
# OUTPUT DIVERGENCE TIMES for this study and then compile for others
#=========
library(ape); library(phytools); library(picante); library(plotrix); library(phangorn); library(phyloch)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")
#setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")

bbone<- "NDexp" # "FBD"

#cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
#head(cladesDR)
allTraits<-read.table("MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)


# Load mamPhy MCC
NDexp0<-read.beast(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target.tre",sep=""))
#NDexp1<-drop.tip2(NDexp0,"_Anolis_carolinensis")
#mamPhy<-ladderize(NDexp1)
mamPhy<-NDexp0
btimes<-branching.times(mamPhy)

pdf(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target_nodeNumbers_ages.pdf",sep=""), onefile=TRUE, width=22, height=122)
plot(mamPhy, show.tip.label=TRUE, cex=0.2)
node.support(mamPhy$height,digits=2,cex=0.4,font=2, mode="numbers", pos="above" )
nodelabels(cex=0.2, adj=c(1,1))
dev.off()

#NDexpNEX<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target.tre",sep=""))
#NDexp1<-drop.tip(NDexpNEX,"_Anolis_carolinensis")
#mamPhyNEX<-ladderize(NDexp1)
#btimesNEX<-branching.times(mamPhyNEX)

means<-mamPhy$height
hpd95_mins<-mamPhy$"height_95%_HPD_MIN"
hpd95_maxs<-mamPhy$"height_95%_HPD_MAX"
ages<-cbind.data.frame(means,hpd95_mins,hpd95_maxs)

write.table(ages, file=paste("divTime_mean95pHPD_MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target.txt",sep=""))

# get the TABLE from the MCC tree.
mamPhyMCC_Table<-read.beast.table(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target.tre",sep=""))

write.table(mamPhyMCC_Table, file=paste("divTime_allTable_MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target.txt",sep=""))


##
# OUTPUT DIVERGENCE TIMES for this study and then compile for others
# Really you should do this for every taxonomic unit...
# ORDERS first
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")


#save.image("divTime_workspace.Rda")
load("divTime_workspace.Rda")

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")

bbone<- "NDexp" # "FBD"

cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""))
colnames(cladesDR)<-c("tiplabel","gen","fam","famLabel","famNumLabel","famNumAll","ord","ordNums","ordLabel1","ordLabel2","clade","cladeCommonAll","cladeCommonSubs","cladeCommonSubs2","cladeCommonSubs2Nums","cladeCombo","higher","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")
head(cladesDR)

mamPhyMCC_Table<-read.table(file=paste("divTime_allTable_MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target.txt",sep=""))

mamPhyMCC_Table_sub<-read.table(file=paste("divTime_mean95pHPD_MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target.txt",sep=""))

ord<-order(mamPhyMCC_Table_sub[,1], decreasing=TRUE)
sorted<-mamPhyMCC_Table_sub[ord,]
write.table(sorted,file=paste("divTime_mean95pHPD_MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target_SORTED-mean.txt",sep=""))

ord2<-order(mamPhyMCC_Table_sub[,3], decreasing=TRUE)
sorted2<-mamPhyMCC_Table_sub[ord,]
write.table(sorted2,file=paste("divTime_mean95pHPD_MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target_SORTED-max.txt",sep=""))



toPlotFIN<-read.table(file=paste(bbone,"_BAMMshifts_24shifts_ReadyToPlot.txt",sep=""), header=TRUE)

ages<-vector("list",length(toPlotFIN$nodes))
for (i in 1:length(toPlotFIN$nodes)){
	ages[[i]]<-mamPhyMCC_Table[which(mamPhyMCC_Table$node==(toPlotFIN$nodes[i]+2)),c("height","height_95._HPD_MIN","height_95._HPD_MAX")]
}
agesAll<-do.call(rbind,ages)
agesAllClades<-cbind(agesAll,toPlotFIN)
colnames(agesAllClades)<-c("mean", "lower","upper","nodes","ID","factor","shift","CLADE")

write.table(agesAllClades,file=paste("divTime_24BAMMnodes_",bbone,"_MCC_target.txt", sep=""))

##
# GET the tip names of the taxa in each rate-shift...
library(phytools)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"
mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")

agesAllClades_orig<-read.table(file=paste("divTime_24BAMMnodes_",bbone,"_MCC_target.txt", sep=""))
nodeNums<-agesAllClades_orig$nodes[c(1,3,5,7,8,10,11,13,15,17:29,31,32)]
IDs<-agesAllClades_orig$ID[c(1,3,5,7,8,10,11,13,15,17:29,31,32)]
CLADEs<-substr(agesAllClades_orig$CLADE[c(1,3,5,7,8,10,11,13,15,17:29,31,32)],1,10)

###
##
# for BODY MASSES
###
BM_table<-read.table(file="MamPhy_5911sp_GenSp_bodyMasses_FIN_5327species.txt", header=TRUE)

cladeBM_all<-vector("list",length(nodeNums))
for(i in 1:length(nodeNums)){
	tipNums<-getDescendants(mamMCC,nodeNums[i])
	cladeTips<-as.vector(na.omit(mamMCC$tip.label[tipNums]))
	cladeBM_i<-as.numeric(na.omit(BM_table[match(cladeTips,as.vector(BM_table$mamPhy_tiplabel)),"BM_final"]))
	cladeBM_all[[i]]<-cladeBM_i
}

# DENS plot...
allMams<-log(na.omit(BM_table[,"BM_final"]))
pdf(file="bodyMasses_24BAMMnodes_vs_allMam_Densitys.pdf", onefile=TRUE, width=8,height=8)
#quartz(width=8,height=8)
op <- par(oma = c(0,0,3,0) + 0.1,		#c(bottom, left, top, right)
          mar = rep(1,4) + 0.1)
layout(matrix(1:25,5,5, byrow=TRUE))

for (i in 1:length(cladeGEO_all)){
	dat<-log(cladeBM_all[[i]])

	plot(density(dat), col="dark grey", main="", bty="n", axes=F, xlim=range(allMams),xlab="",ylab="", zero.line=FALSE)
	polygon(density(dat), col="dark grey", border="dark grey", bty="n")
	x.tick <- quantile(allMams, c(0.01,0.5,0.99))
	axis(at=c(x.tick), labels=round(x.tick,digits=1), side=1, line=0, cex=2.5, lwd=1, tck=0.01, cex.axis=0.8, mgp=c(1,0,0))

	dens.rate <- density(dat)$y
	#axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=1, las=1, lwd=1, cex.axis=1, tck=-0.05, mgp=c(1,2,0))
	seg.tick <- quantile(allMams, c(0.01,0.5,0.99))
	segments(seg.tick[[1]],0,seg.tick[[1]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")
	segments(seg.tick[[2]],0,seg.tick[[2]],max(dens.rate)*1, lty=2, lwd=2,col="black")
	segments(seg.tick[[3]],0,seg.tick[[3]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")

	abline(v=median(dat), lty=1,col="red",lwd=2)
	
	test<-wilcox.test(x=allMams, y=dat, alternative="two.sided", conf.int=FALSE)
	if(test$p.value < 0.001){ 
		mtext(side=3,text=paste(IDs[i]," ***",sep=""), line=0, cex=1, font=2)
	} else if(test$p.value < 0.01) { 
			mtext(side=3,text=paste(IDs[i]," **",sep=""), line=0, cex=1, font=2)
	} else if(test$p.value < 0.05) {
			mtext(side=3,text=paste(IDs[i]," *",sep=""), line=0, cex=1, font=2)
	} else {
		mtext(side=3,text=paste(IDs[i],"; ns",sep=""), line=0, cex=1, font=2)
		}
}
title(main="Body masses (F&S 2016, EltonTraits 2014); 5327 spp (90.1%)", xlab = "",
      ylab = "",
      outer = TRUE, line = 1,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)

dev.off()



###
##
# for GEO RANGES
###
GEO_table<-read.table(file="MamPhy_5911sp_GenSp_RangesIUCN2012_FIN_5097species.txt", header=TRUE)

cladeGEO_all<-vector("list",length(nodeNums))
for(i in 1:length(nodeNums)){
	tipNums<-getDescendants(mamMCC,nodeNums[i])
	cladeTips<-as.vector(na.omit(mamMCC$tip.label[tipNums]))
	cladeGEO_i<-as.numeric(na.omit(GEO_table[match(cladeTips,as.vector(GEO_table$mamPhy_tiplabel)),"RangeIUCN"]))
	cladeGEO_all[[i]]<-cladeGEO_i
}
cladeGEO_all[[6]]<-NA # Because no species in CETACEA have geographic ranges-- AH, because focused on LAND areas...


allMams<-log(as.numeric(na.omit(GEO_table[,"RangeIUCN"])))
allMams_lines<-quantile(allMams,c(0.025,0.50,0.975))

# DENS plot...
pdf(file="geoRanges_24BAMMnodes_vs_allMam_Densitys.pdf", onefile=TRUE, width=8,height=8)
#quartz(width=8,height=8)
op <- par(oma = c(0,0,3,0) + 0.1,		#c(bottom, left, top, right)
          mar = rep(1.5,4) + 0.1)
layout(matrix(1:25,5,5, byrow=TRUE))

for (i in 1:length(cladeGEO_all)){
	dat<-log(cladeGEO_all[[i]])
	if(i==6) {	plot(density(log(cladeGEO_all[[5]])), col="dark grey", type="n", main="", bty="n", axes=F, xlim=range(allMams),xlab="",ylab="", zero.line=FALSE)
		mtext(side=3,text=IDs[i], line=-3, cex=1, font=2)
	} else {

	plot(density(dat), col="dark grey", main="", bty="n", axes=F, xlim=range(allMams),xlab="",ylab="", zero.line=FALSE)
	polygon(density(dat), col="dark grey", border="dark grey", bty="n")
	x.tick <- quantile(allMams, c(0.01,0.5,0.99))
	axis(at=c(x.tick), labels=round(x.tick,digits=1), side=1, line=0, cex=2.5, lwd=1, tck=0.01, cex.axis=1, mgp=c(1,0,0))

	dens.rate <- density(dat)$y
	#axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=1, las=1, lwd=1, cex.axis=1, tck=-0.05, mgp=c(1,2,0))
	seg.tick <- quantile(allMams, c(0.01,0.5,0.99))
	segments(seg.tick[[1]],0,seg.tick[[1]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")
	segments(seg.tick[[2]],0,seg.tick[[2]],max(dens.rate)*1, lty=2, lwd=2,col="black")
	segments(seg.tick[[3]],0,seg.tick[[3]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")

	abline(v=median(dat), lty=1,col="red",lwd=2)

	test<-wilcox.test(x=allMams, y=dat, alternative="two.sided", conf.int=FALSE)
	if(test$p.value < 0.001){ 
		mtext(side=3,text=paste(IDs[i]," ***",sep=""), line=0, cex=1, font=2)
	} else if(test$p.value < 0.01) { 
			mtext(side=3,text=paste(IDs[i]," **",sep=""), line=0, cex=1, font=2)
	} else if(test$p.value < 0.05) {
			mtext(side=3,text=paste(IDs[i]," *",sep=""), line=0, cex=1, font=2)
	} else {
		mtext(side=3,text=paste(IDs[i],"; ns",sep=""), line=0, cex=1, font=2)
		}
	}
}
title(main="Geographic ranges (IUCN 2012, 360 x 114 grids); 5097 spp (86.2%)", xlab = "",
      ylab = "",
      outer = TRUE, line = 1,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)
dev.off()




###
##
# for TIP DR... 
###
bbone<- "NDexp" # "FBD"
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
head(cladesDR)

cladesDR[which(cladesDR$gen=="Ctenomys"),"harmMeans"]

cladesDR[which(cladesDR$gen=="Homo"),"harmMeans"]

cladesDR[which(cladesDR$gen=="Trachypithecus"),"harmMeans"]


cerco<-cladesDR[which(cladesDR$fam=="CERCOPITHECIDAE"),"harmMeans"]



cladeDR_all<-vector("list",length(nodeNums))
for(i in 1:length(nodeNums)){
	tipNums<-getDescendants(mamMCC,nodeNums[i])
	cladeTips<-as.vector(na.omit(mamMCC$tip.label[tipNums]))
	cladeDR_i<-cladesDR[match(cladeTips,as.vector(cladesDR$tiplabel)),"harmMeans"]
	names(cladeDR_i)<-cladeTips
	cladeDR_all[[i]]<-cladeDR_i
}

allMams<-cladesDR[,"harmMeans"]
allMams_lines<-quantile(allMams,c(0.025,0.50,0.975))

# DENS plot...
pdf(file="tipDR_24BAMMnodes_vs_allMam_Densitys.pdf", onefile=TRUE, width=8,height=8)
#quartz(width=8,height=8)
op <- par(oma = c(0,0,3,0) + 0.1,		#c(bottom, left, top, right)
          mar = rep(1,4) + 0.1)
layout(matrix(1:25,5,5, byrow=TRUE))

for (i in 1:length(cladeGEO_all)){
	dat<-(cladeDR_all[[i]])

	plot(density(dat), col="dark grey", main="", bty="n", axes=F, xlim=range(allMams),xlab="",ylab="", zero.line=FALSE)
	polygon(density(dat), col="dark grey", border="dark grey", bty="n")
	x.tick <- quantile(allMams, c(0.01,0.5,0.99))
	axis(at=c(x.tick), labels=round(x.tick,digits=1), side=1, line=0, cex=2.5, lwd=1, tck=0.01, cex.axis=1, mgp=c(1,0,0))

	dens.rate <- density(dat)$y
	#axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=1, las=1, lwd=1, cex.axis=1, tck=-0.05, mgp=c(1,2,0))
	seg.tick <- quantile(allMams, c(0.01,0.5,0.99))
	segments(seg.tick[[1]],0,seg.tick[[1]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")
	segments(seg.tick[[2]],0,seg.tick[[2]],max(dens.rate)*1, lty=2, lwd=2,col="black")
	segments(seg.tick[[3]],0,seg.tick[[3]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")

	abline(v=median(dat), lty=1,col="red",lwd=2)

	test<-wilcox.test(x=allMams, y=dat, alternative="two.sided", conf.int=FALSE)
	if(test$p.value < 0.001){ 
		mtext(side=3,text=paste(IDs[i]," ***",sep=""), line=0, cex=1, font=2)
	} else if(test$p.value < 0.01) { 
			mtext(side=3,text=paste(IDs[i]," **",sep=""), line=0, cex=1, font=2)
	} else if(test$p.value < 0.05) {
			mtext(side=3,text=paste(IDs[i]," *",sep=""), line=0, cex=1, font=2)
	} else {
		mtext(side=3,text=paste(IDs[i],"; ns",sep=""), line=0, cex=1, font=2)
		}
}
title(main="Tip DR (this study) for all 5911 spp", xlab = "",
      ylab = "",
      outer = TRUE, line = 1,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)

dev.off()


##
# with TIP DR-- make plotCI BARS rather than density plots to show DIFF in medians from all Mammalia.
setwd()

bbone<- "NDexp" # "FBD"
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
head(cladesDR)

#orderToPlot<-c("RODENTIA", "MouseRelated", "SquirrelRelated", "GuineaPigRelated", "CHIROPTERA", "Yinpterochiroptera", "Yangochiroptera", "EULIPOTYPHLA", "PRIMATES", "OldWorldMonkey2", "NewWorldMonkey2", "Strepsirrhini", "CETARTIODACTYLA", "Ruminantia", "Cetaceans", "AFROSORICIDA", "SCANDENTIA", "MACROSCELIDEA", "PILOSA", "CARNIVORA","CatRelated","DogRelated", "LAGOMORPHA", "PERISSODACTYLA", "PHOLIDOTA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA")
orderToPlot_ords<-c("RODENTIA", "CHIROPTERA", "EULIPOTYPHLA", "PRIMATES", "CETARTIODACTYLA", "AFROSORICIDA", "SCANDENTIA", "MACROSCELIDEA", "PILOSA", "CARNIVORA","LAGOMORPHA", "PERISSODACTYLA", "PHOLIDOTA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA")
orderToPlot_clades<-c("Mouse-related", "Squirrel-related", "Guinea_pig-related", "Yinpterochiroptera", "Yangochiroptera", "Catarrhini", "Playtrrhini", "Strepsirrhini", "Ruminantia", "Whippomorpha", "Feliformes","Caniformes")

orderToPlot<-c("RODENTIA", "Mouse-related", "Squirrel-related", "Guinea_pig-related", "CHIROPTERA", "Yinpterochiroptera", "Yangochiroptera", "EULIPOTYPHLA", "PRIMATES", "Catarrhini", "Playtrrhini", "Strepsirrhini", "CETARTIODACTYLA", "Ruminantia", "Whippomorpha", "AFROSORICIDA", "SCANDENTIA", "MACROSCELIDEA", "PILOSA", "CARNIVORA","Feliformes","Caniformes", "LAGOMORPHA", "PERISSODACTYLA", "PHOLIDOTA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA")

cladeDR_all<-vector("list",length(orderToPlot))
for(i in 1:length(orderToPlot)){
	if (orderToPlot[i] %in% cladesDR$ord) {
		cladeDR_i<-cladesDR[which(cladesDR$ord==orderToPlot[i]),"harmMeans"]
		cladeDR_all[[i]]<-cladeDR_i		
	} else {
		cladeDR_i<-cladesDR[which(cladesDR$clade==orderToPlot[i]),"harmMeans"]
		cladeDR_all[[i]]<-cladeDR_i		
	}
}
medianLowHigh<-cbind.data.frame(unlist(lapply(cladeDR_all,length)),unlist(lapply(cladeDR_all,median)),as.vector(unlist(lapply(cladeDR_all,quantile,c(0.025)))),as.vector(unlist(lapply(cladeDR_all,quantile,c(0.975)))))
rownames(medianLowHigh)<-orderToPlot
colnames(medianLowHigh)<-c("richness","median","low","high")
write.table(medianLowHigh,file="medianLowHigh_DRtable_28clades_inFig2.txt")

marsupials<-c("DIPROTODONTIA", "DIDELPHIMORPHIA", "DASYUROMORPHIA", "PERAMELEMORPHIA", "NOTORYCTEMORPHIA", "PAUCITUBERCULATA", "MICROBIOTHERIA")
cladeDR_all_marsup<-vector("list",length(marsupials))
for(i in 1:length(marsupials)){
	cladeDR_i<-cladesDR[which(cladesDR$ord==marsupials[i]),"harmMeans"]
	cladeDR_all_marsup[[i]]<-cladeDR_i		
}
medianLowHigh_marsup<-cbind.data.frame(unlist(lapply(cladeDR_all_marsup,length)),unlist(lapply(cladeDR_all_marsup,median)),as.vector(unlist(lapply(cladeDR_all_marsup,quantile,c(0.025)))),as.vector(unlist(lapply(cladeDR_all_marsup,quantile,c(0.975)))))
rownames(medianLowHigh_marsup)<-marsupials
colnames(medianLowHigh_marsup)<-c("richness","median","low","high")
write.table(medianLowHigh_marsup,file="medianLowHigh=marsup_DRtable_28clades_inFig2.txt")






sorted<-medianLowHigh[order(medianLowHigh[,1]),]

                       median        low       high
PILOSA             0.04263498 0.02424200 0.05781400
PHOLIDOTA          0.05737702 0.05375663 0.07193926
MACROSCELIDEA      0.06070134 0.04431158 0.06752350
HYRACOIDEA         0.06254040 0.05627257 0.06262370
SIRENIA            0.06585778 0.04241215 0.06755565
AFROSORICIDA       0.08249114 0.03240265 0.14896045
PROBOSCIDEA        0.08925329 0.05876137 0.09305938
SCANDENTIA         0.11041779 0.02228780 0.13666994
PERISSODACTYLA     0.11556987 0.06332399 0.17091039
Yangochiroptera    0.16187852 0.05914970 0.35054232
CHIROPTERA         0.18068497 0.06158593 0.44167188
Mouse-related      0.20318177 0.06057697 0.42950652
RODENTIA           0.20845656 0.06070805 0.48100742
Squirrel-related   0.21288570 0.07555553 0.50097712
Yinpterochiroptera 0.22780187 0.07177901 0.64673471
EULIPOTYPHLA       0.22849398 0.04968758 0.48496836
Caniformes         0.23190544 0.08431320 0.42373238
Guinea_pig-related 0.23316745 0.05835549 0.68335047
CARNIVORA          0.23914381 0.08478853 0.45797944
Whippomorpha       0.24061360 0.07541643 0.60184203
CETARTIODACTYLA    0.24698987 0.07873310 0.53007021
Strepsirrhini      0.25098542 0.07687970 0.39405249
Feliformes         0.25640619 0.08968520 0.47219152
Ruminantia         0.25998217 0.09103792 0.49871943
LAGOMORPHA         0.26418723 0.11028174 0.54269508
PRIMATES           0.32456199 0.09101976 0.48014675
Playtrrhini        0.34160779 0.17159784 0.51147736
Catarrhini         0.34531093 0.18465632 0.48566215

## awesome.





###########
#####
# PAIRWISE plots of tipDR -- geoRange -- bodyMass -- LIFEMODE
bbone<- "NDexp" # "FBD"
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)

#priDR<-cladesDR[which(cladesDR$ord=="PRIMATES"),c("tiplabel","harmMeans")]
#cercDR<-cladesDR[which(cladesDR$fam=="CERCOPITHECIDAE"),"harmMeans"]
#traDR<-cladesDR[which(cladesDR$gen=="Trachypithecus"),]

DRrange <- quantile(cladesDR[,"harmMeans"], seq(0,1, 0.01))[2:101]

GEO_table<-read.table(file="MamPhy_5911sp_GenSp_RangesIUCN2012_FIN_5097species.txt", header=TRUE)

BM_table<-read.table(file="MamPhy_5911sp_GenSp_bodyMasses_FIN_5327species.txt", header=TRUE)

tipDR1<-cladesDR[,c("tiplabel","harmMeans")]
ordered<-order(tipDR1[,1])
tipDR<-tipDR1[ordered,]

tipRange1<-GEO_table[,c("mamPhy_tiplabel","RangeIUCN")]
ordered<-order(tipRange1[,1])
tipRange<-tipRange1[ordered,]

tipMass1<-BM_table[,c("mamPhy_tiplabel","BM_final")]
ordered<-order(tipMass1[,1])
tipMass<-tipMass1[ordered,]

tipDataAll<-cbind.data.frame(tipDR[,2],tipRange[,2],tipMass[,2])
rownames(tipDataAll)<-tipDR[,1]
colnames(tipDataAll)<-c("tipDR","geoRange","bodyMass")
write.table(tipDataAll,file="MamPhy_5911sp_tiplabel_DR-range-mass.txt")

tipDataAllHerb<-read.table(file="MamPhy_5911sp_tiplabel_DR-range-mass-herb.txt")

##
# Get JUST the rate-shifted data subset...
agesAllClades_orig<-read.table(file=paste("divTime_24BAMMnodes_",bbone,"_MCC_target.txt", sep=""))
nodeNums<-agesAllClades_orig$nodes[c(1,3,5,7,8,10,11,13,15,17:29,31,32)]
cladeNames<-agesAllClades_orig$ID[c(1,3,5,7,8,10,11,13,15,17:29,31,32)]

cladeTips_all<-vector("list",length(nodeNums))
for(i in 1:length(nodeNums)){
	tipNums<-getDescendants(mamMCC,nodeNums[i])
	cladeTips<-as.vector(na.omit(mamMCC$tip.label[tipNums]))
	cladeTips_all[[i]]<-cladeTips
}
cladeTips_all[[6]]<-NA # Because no species in CETACEA have geographic ranges-- AH, because focused on LAND areas...
cladeTips_all[[3]]<-NA # Is all Placentals
cladeTips_all[[1]]<-NA # is all Marsupials
cladeTips_all[[4]]<-NA # is all Carnivora, and they have BIGGER ranges... 

cladeTips_all_Combo<-unlist(cladeTips_all)

cladeTips_all_uniq<-unique(cladeTips_all_Combo) # 4205 species involved in some rate shift... 

shiftedDat<-tipDataAll[match(cladeTips_all_uniq,rownames(tipDataAll)),]
 

# ==================
#######
# PGLS start here...
library(ape); library(phytools); library(picante); library(geiger); library(moments); library(nlme)
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
source("DR_functions.R")

library(foreach);library(doSNOW)
cl = makeCluster(45, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees=100

foreach(i=1:ntrees, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")

# which backbone?
bbone<- "NDexp" #"FBD" # 

#mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)


# read in 1 of 100 full trees
#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_nexus-and-newickTrees/")
mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_nexus.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 
	#write.tree(mamPhy,file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""))
tree1=scan(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""), what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

# change to NEW directory
###
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/_ExplainingTipDR_speciesLevelRuns")
#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_ExplainingTipDR_speciesLevelRuns")


# ES and DR on per-tip basis
#==================
# gives pairwise clade matrix from CAIC function
#clade_matrix = readCAIC(tree1)
#
#ES = ES_v2(clade_matrix)
#DR = 1/ES
#res = cbind.data.frame(DR,ES)
#res1 = res[order(rownames(res)),]
#
#write.table(res1, file=paste(bbone,"_sample100_",i,"_DRtips.txt",sep=""))
res1<-read.table(file=paste(bbone,"_sample100_",i,"_DRtips.txt",sep=""))


# -ALL- traits
tipDataAll<-read.table("MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)
head(tipDataAll)

## Traits -TO LOOK AT- 
## >> ADDING in the per-tree tipDR action::
traitVars1<-c("BM_final_kg","homeRange_km2_ext","geoRange_km2","DispDistAll_km_ext","GenerationLength_d","trophicCat","activityCat3","lifemode1234","lifemodeCat6")
traitDat1<-tipDataAll[,c("tiplabel",traitVars1)]
traitDat<-cbind.data.frame(traitDat1$tiplabel,log(res1$DR),log(traitDat1$BM_final_kg),log(traitDat1$homeRange_km2_ext),log(traitDat1$geoRange_km2),log(traitDat1$DispDistAll_km_ext),(traitDat1$lifemode1234*log(traitDat1$DispDistAll_km_ext)),log(traitDat1$GenerationLength_d), traitDat1$trophicCat,traitDat1$activityCat3,traitDat1$lifemode1234,traitDat1$lifemodeCat6)
traitNames<-c("Tip DR (log species/Ma)","Body mass (log kg)","Home range (log km2)","Geographic range (log km2)","Dispersal distance (log km)","Dispersal index (log km x lifemode)","Generation Length (log d)","trophicCat","activityCat3","lifemode1234","lifemodeCat6")
traitVars<-c("logTipDR","logBodyMass","logHomeRange","logGeoRange","logDispDist","logDispIndex","logGenLength","trophicCat","activityCat3","lifemode1234","lifemodeCat6")
colnames(traitDat)<-c("tiplabel",traitVars)

##
# Subset out the tips per CATEGORY
categoricalVars<-c("trophicCat","activityCat3","lifemode1234","lifemodeCat6")

ALLvars_perCat_Data<-vector("list",length(categoricalVars))
ALLvars_perCat_sums<-vector("list",length(categoricalVars))

for(j in 1:length(categoricalVars)){
	varCats<-names(table(traitDat[,categoricalVars[j]]))
	varDat<-traitDat[,categoricalVars[j]]

	if(j==3){
	varCatNames<-c("TerrArboScan")
	perCat_Data<-traitDat[which(varDat==2),]
	rownames(perCat_Data)<-as.vector(perCat_Data$tiplabel)
	ALLvars_perCat_Data[[j]]<-perCat_Data
	sum<-length(perCat_Data$tiplabel)
	names(sum)<-varCatNames
	ALLvars_perCat_sums[[j]]<-sum
	} else {
 
	varCatNames<-names(table(traitDat[,categoricalVars[j]]))
	perCat_Data<-vector("list",length(varCats))
	sum<-vector()
	for (k in 1:length(varCats)){
		perCat_Data[[k]]<-traitDat[which(varDat==varCats[k]),]
		rownames(perCat_Data[[k]])<-as.vector(perCat_Data[[k]]$tiplabel)
		sum[k]<-length(perCat_Data[[k]]$tiplabel)
	}
	names(sum)<-varCatNames
	ALLvars_perCat_Data[[j]]<-perCat_Data
	ALLvars_perCat_sums[[j]]<-sum
	}
}
# got it.

predictorVars<-c("logBodyMass","logHomeRange","logGeoRange","logDispDist","logDispIndex","logGenLength")

## do the lambda calcs-- on the RESIDUALS of the GLS
ALL_perPred_lams<-vector("list",length(categoricalVars))
for (j in 1:length(categoricalVars)){
	if(j==3){varCatNames<-c("Subt", "TerrArboScan", "Aquatic", "Flying")} else {
		varCatNames<-names(table(traitDat[,categoricalVars[j]]))}
	varCats<-names(table(traitDat[,categoricalVars[j]]))

	perPred_lams<-vector("list",length(varCats))
	for (k in 1:length(varCats)){
		lam<-c()
		for (m in 1:length(predictorVars)){
		dat<-na.omit(ALLvars_perCat_Data[[j]][[k]][,c("logTipDR",predictorVars[m])])

		treeDat<-treedata(mamPhy,dat)
		datPlot<-as.data.frame(treeDat$dat)

		form<-as.formula(paste(predictorVars[m], " ~ ", "logTipDR",sep=""))
		fitGLS<-gls(form, data=dat, method="ML")

		resids<-residuals(fitGLS,type="n")[1:length(datPlot[,1])]

		#sigTest<-phylosig(tree=treeDat$phy,x=resids, method="K")#,test=TRUE)#,start=0)
		sigTest<-phylosig(tree=treeDat$phy,x=resids, method="lambda",test=TRUE)#,start=0)
		lam[m]<-sigTest$lambda
		}
	perPred_lams[[k]]<-lam
	}
ALL_perPred_lams[[j]]<-perPred_lams
}
write.table(ALL_perPred_lams,file=paste("ALL_perPred_lams_tree",i,".txt",sep=""))

###
# Do the PGLS then...
ALL_fitDat_perCat<-vector("list",length(categoricalVars))
for (j in 1:length(categoricalVars)){
	if(j==3){
	varCatNames<-c("TerrArboScan")
	} else {
	varCatNames<-names(table(traitDat[,categoricalVars[j]]))
	}

	fitDat_perCat<-vector("list",length(varCatNames))
	for (k in 1:length(varCatNames)){
		fitDat_perPred<-vector("list",length(predictorVars))
		for (m in 1:length(predictorVars)){

		if(j==4 && k==1){
			fitDat_perPred[[m]]<-rep(NA,5)
		} 

		else if (j==3) {
		dat<-na.omit(ALLvars_perCat_Data[[j]][,c("logTipDR",predictorVars[m])])
		
		treeDat<-treedata(mamPhy,dat)
		datPlot<-as.data.frame(treeDat$dat)

		form<-as.formula(paste(predictorVars[m], " ~ ", "logTipDR",sep=""))
		
#		fit1<-gls(form, data=dat, method="ML")
		for (p in seq(0,1,by=0.01)) {possibleError <- tryCatch(
		      gls(form, correlation=corPagel(value=p,phy=treeDat$phy), data=dat, method="ML"),
		      error=function(e) e)
		if(inherits(possibleError, "gls")) break		
		if(inherits(possibleError, "error")) next}
		fit1<-possibleError

		sum1<-summary(fit1)

		a=round(sum1$tTable[1], digits=3)
		b=round(sum1$tTable[2], digits=3)
		SE<-round(sum1$tTable[4], digits=3)
		pVal=round(sum1$tTable[8], digits=3)
		lambda=round(sum1$modelStruct[[1]][[1]], digits=3)
#		lambda=NA

		fitDat_perPred[[m]]<-c(a,b,SE,pVal,lambda)

		} else {

		dat<-na.omit(ALLvars_perCat_Data[[j]][[k]][,c("logTipDR",predictorVars[m])])

		treeDat<-treedata(mamPhy,dat)
		datPlot<-as.data.frame(treeDat$dat)

		form<-as.formula(paste(predictorVars[m], " ~ ", "logTipDR",sep=""))
		
#		fit1<-gls(form, data=dat, method="ML")
		for (p in seq(0,1,by=0.01)) {possibleError <- tryCatch(
		      gls(form, correlation=corPagel(value=p,phy=treeDat$phy), data=dat, method="ML"),
		      error=function(e) e)
		if(inherits(possibleError, "gls")) break		
		if(inherits(possibleError, "error")) next}
		fit1<-possibleError

		sum1<-summary(fit1)

		a=round(sum1$tTable[1], digits=3)
		b=round(sum1$tTable[2], digits=3)
		SE<-round(sum1$tTable[4], digits=3)
		pVal=round(sum1$tTable[8], digits=3)
		lambda=round(sum1$modelStruct[[1]][[1]], digits=3)
#		lambda=NA

		fitDat_perPred[[m]]<-c(a,b,SE,pVal,lambda)
		}
		}

	res<-do.call(rbind,fitDat_perPred)
	rownames(res)<-paste(predictorVars,"_",varCatNames[k],sep="")
	colnames(res)<-c("a","b","SE","pVal","lam")
	fitDat_perCat[[k]]<-res
	}
ALL_fitDat_perCat[[j]]<-do.call(rbind,fitDat_perCat)
}
Final_ALL_fitDat_perCat<-do.call(rbind,ALL_fitDat_perCat)

#model<-"GLS"
model<-"PGLS-Pagel"
write.table(Final_ALL_fitDat_perCat,file=paste("ExplainingDR_spLevel_",model,"predAsResp-vs-tipDR_NDexp_sample100_tree",i,".txt",sep=""))

}






# by TROPHIC LEVEL
# load data
trophic_cats<-c("Herbivore","Omnivore","Carnivore")
trophicData1<-vector("list",length(trophic_cats))
	trophicData1[[1]]<-tipDataAll[which(tipDataAll$percentHerb==100),c("tipDR",var)]
	trophicData1[[2]]<-tipDataAll[which(tipDataAll$percentHerb <= 90 & tipDataAll$percentHerb >= 10 ),c("tipDR",var)]
	trophicData1[[3]]<-tipDataAll[which(tipDataAll$percentHerb==0),c("tipDR",var)]
	# this already removes all NAs

trophicData<-vector("list",length(trophic_cats))
for (j in 1:length(trophic_cats)){
	tips<-rownames(trophicData1[[j]])
	tipDR<-res1[match(tips,rownames(res1)),][,1]
	geoRange<-trophicData1[[j]][,2]
	trophicData[[j]]<-cbind.data.frame(tipDR,geoRange)
	rownames(trophicData[[j]])<-tips
}

# do calcs
abP_trophic_i<-vector("list",length(trophic_cats))

for (j in 1:length(trophic_cats)){
	dat<-na.omit(trophicData[[j]])
	treeDat_trophic<-treedata(mamPhy,dat)
	datPlot<-as.data.frame(treeDat_trophic$dat)

	fitPGLS_trophic<-gls(form, data=datPlot,correlation=corPagel(value=0.5,phy=treeDat_trophic$phy))
	#fitPGLS_trophic<-gls(form, data=datPlot,correlation=corBrownian(phy=treeDat_trophic$phy))
	#fitPGLS_trophic<-gls(form, data=datPlot)
	sumPGLS_trophic<-summary(fitPGLS_trophic)

	sum<-sumPGLS_trophic
	a=round(sum$coef[[1]],3)
	b=round(sum$coef[[2]],3)
	pVal=round(sum$tTable[[8]],3)
	lambda=sum$modelStruct[[1]][[1]]
#	lambda=NA

	abP_trophic_i[[j]]<-c(a,b,pVal,lambda)
}
res_trophic_i<-do.call(rbind,abP_trophic_i)
rownames(res_trophic_i)<-trophic_cats
colnames(res_trophic_i)<-c("a","b","pVal","lam")

resCats<-rbind(res_LM_i,res_trophic_i,ALL_mamPhy)

write.table(resCats,file=paste("ExplainingDR_logBM-vs-tipDR_NDexp_sample100_tree",i,"_PGLS-Pagel_estimDR_allRes.txt",sep=""))
#write.table(resCats,file=paste("ExplainingDR_logGeoRangeKM2-vs-tipDR_NDexp_sample100_tree",i,"_PGLS-Pagel_estimDR_allRes.txt",sep=""))
#write.table(resCats,file=paste("ExplainingDR_logGeoRangeKM2-vs-tipDR_NDexp_sample100_tree",i,"_PGLS-Brownian_estimDR_allRes.txt",sep=""))
#write.table(resCats,file=paste("ExplainingDR_logDR-vs-logGeoRange_NDexp_sample100_tree",i,"_PGLS-Brownian_estimDR_allRes.txt",sep=""))
#write.table(resCats,file=paste("ExplainingDR_logDR-vs-logGeoRangeKM2_NDexp_sample100_tree",i,"_GLS_estimDR_allRes.txt",sep=""))

#write.table(resCats,file=paste("ExplainingDR_tipDR-vs-logBM_NDexp_sample100_tree",i,"_PGLS-Brownian_estimDR_allRes.txt",sep=""))
#write.table(resCats,file=paste("ExplainingDR_tipDR-vs-logBM_NDexp_sample100_tree",i,"_GLS_estimDR_allRes.txt",sep=""))

#write.table(resCats,file=paste("ExplainingDR_logDR-vs-logBM_NDexp_sample100_tree",i,"_GLS_estimDR_allRes.txt",sep=""))
#write.table(resCats,file=paste("ExplainingDR_logDR-vs-logBM_NDexp_sample100_tree",i,"_PGLS-Brownian_estimDR_allRes.txt",sep=""))
#write.table(resCats,file=paste("ExplainingDR_logDR-vs-logBM_NDexp_sample100_tree",i,"_PGLS-Pagel_allRes.txt",sep=""))

}







##
# Subset out the tips of TROPHIC LEVEL
trophic_cats<-c("Herbivore","Omnivore","Carnivore")

trophicData<-vector("list",length(trophic_cats))
trophicTips<-vector("list",length(trophic_cats))
	trophicData[[1]]<-tipDataAll[which(tipDataAll$percentHerb==100),]
	trophicTips[[1]]<-rownames(trophicData[[1]])
	trophicData[[2]]<-tipDataAll[which(tipDataAll$percentHerb < 100 & tipDataAll$percentHerb > 0 ),]
	trophicTips[[2]]<-rownames(trophicData[[2]])
	trophicData[[3]]<-tipDataAll[which(tipDataAll$percentHerb==0),]
	trophicTips[[3]]<-rownames(trophicData[[3]])

# Herb = 1620
# Omni = 1805
# Carn = 1561



###
# RUN PGLS on 100 trees::
library(geiger); library(nlme); library(ape)
#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_nexus-and-newickTrees")
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/_ExplainingTipDR_timeSliceRuns")

# ready nodes
#library(foreach);library(doSNOW)
#cl = makeCluster(45, type = 'SOCK', outfile="")
#registerDoSNOW(cl)
#

ntrees=10
foreach(i=1:ntrees, .packages=c('ape', 'geiger', 'nlme'), .combine=cbind, .verbose=TRUE) %dopar% {

setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")

# which backbone?
bbone<- "NDexp" #"FBD" # 

# load tree
mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_nexus.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 
#write.tree(mamPhy,file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""))
#tree1=scan(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_newick.tre",sep=""), what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_RESULTS_tipDR-vs-bodyMass_PGLS/DR_per_100trees")
#==================
# ES and DR on per-tip basis
#==================
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/_ExplainingTipDR_timeSliceRuns")

tipDR<-res1$DR
names(tipDR)<-rownames(res1)

#cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
#head(cladesDR)

tipDataAll<-read.table(file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)

datSelected<-tipDataAll[,c(7:14,28:32,37:39,47:48,54:56,58:60,61)] # 25 vars

datLogged<-log(datSelected[,c("BM_final_kg","homeRange_km2_ext","DispDistAll_km_ext","GenerationLength_d")]) #"geoRange_km2",
colnames(datLogged)<-c("logBM_final_kg","logHomeRange_km2_ext","logDispDistAll_km_ext","logGenerationLength_d") #"logGeoRange_km2",

#datFinal<-cbind.data.frame(datSelected[,c(1:16,20:21,23:24)],datLogged)
datFinal<-datLogged
predictors<-colnames(datFinal)

perPredictor<-vector("list",length(predictors))
for(j in 1:length(predictors)){

var<-predictors[j]

dat1<-cbind.data.frame(tipDR,datFinal[,var])
colnames(dat1)<-c("tipDR",var)
rownames(dat1)<-names(tipDR)

dat<-na.omit(dat1)
treeDat<-treedata(mamPhy,dat)
datPlot<-as.data.frame(treeDat$dat)

#	form<-as.formula(paste(respVar[z], " ~ ", predictors[k], sep=""))
	form<-as.formula(paste(predictors[j], " ~ ","log(tipDR)",sep=""))
#	fit1<-gls(form, data=dat, method="ML")
	for (p in seq(0,1,by=0.01)) {possibleError <- tryCatch(
	      gls(form, correlation=corPagel(value=p,phy=treeDat$phy), data=dat, method="ML"),
	      error=function(e) e)
	if(inherits(possibleError, "gls")) break		
	if(inherits(possibleError, "error")) next}
	fit1<-possibleError

	sum<-summary(fit1)

	a=round(sum$tTable[1], digits=3)
	b=round(sum$tTable[2], digits=3)
	SE=round(sum$tTable[4], digits=3)
	pVal=round(sum$tTable[8], digits=3)
	lambda=round(sum1$modelStruct[[1]][[1]], digits=3)
#	lambda=NA

	perPredictor[[j]]<-c(a,b,SE,pVal,lambda)
}
allRes_i<-cbind(do.call(rbind,perPredictor),rep(i,length(predictors)))
colnames(allRes_i)<-c("a","b","SE","pVal","lam","tree")
rownames(allRes_i)<-predictors

setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")

write.table(allRes_i,file=paste("speciesLevel_PGLS-PAGEL_logTipDR_conVars-as-response_tree",i,".txt",sep=""))

} # cycle 100 trees



allRes_ALL[[i]]<-



#####
# Script to SUMMARIZE those 100 tree outputs
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_RESULTS_tipDR-vs-bodyMass_PGLS")
ntrees<-100
resALL<-vector("list",length(ntrees))
Terrestrial <-data.frame(matrix(NA,ncol=5,nrow=ntrees))
Flying      <-data.frame(matrix(NA,ncol=5,nrow=ntrees))
Arboreal    <-data.frame(matrix(NA,ncol=5,nrow=ntrees))
Scansorial  <-data.frame(matrix(NA,ncol=5,nrow=ntrees))
Subterranean<-data.frame(matrix(NA,ncol=5,nrow=ntrees))
Aquatic     <-data.frame(matrix(NA,ncol=5,nrow=ntrees))
Herbivore   <-data.frame(matrix(NA,ncol=5,nrow=ntrees))
Omnivore    <-data.frame(matrix(NA,ncol=5,nrow=ntrees))
Carnivore   <-data.frame(matrix(NA,ncol=5,nrow=ntrees))
ALL_mamPhy  <-data.frame(matrix(NA,ncol=5,nrow=ntrees))

for(i in 1:ntrees){
	#resCats_i<-read.table(file=paste("ExplainingDR_logDR-vs-logBM_NDexp_sample100_tree",i,"_PGLS-Brownian_estimDR_allRes.txt",sep=""))
	#resCats_i<-read.table(file=paste("ExplainingDR_logDR-vs-logBM_NDexp_sample100_tree",i,"_PGLS-Brownian_allRes.txt",sep=""))
	#resCats_i<-read.table(file=paste("ExplainingDR_tipDR-vs-logBM_NDexp_sample100_tree",i,"_PGLS-Brownian_estimDR_allRes.txt",sep=""))
	#resCats_i<-read.table(file=paste("ExplainingDR_tipDR-vs-logBM_NDexp_sample100_tree",i,"_GLS_estimDR_allRes.txt",sep=""))
	#resCats_i<-read.table(file=paste("ExplainingDR_logDR-vs-logBM_NDexp_sample100_tree",i,"_GLS_estimDR_allRes.txt",sep=""))
	#resCats_i<-read.table(file=paste("ExplainingDR_logDR-vs-logBM_NDexp_sample100_tree",i,"_GLS_allRes.txt",sep=""))
	#resCats_i<-read.table(file=paste("ExplainingDR_logDR-vs-logGeoRangeKM2_NDexp_sample100_tree",i,"_GLS_estimDR_allRes.txt",sep=""))
	#resCats_i<-read.table(file=paste("ExplainingDR_logDR-vs-logGeoRange_NDexp_sample100_tree",i,"_PGLS-Brownian_estimDR_allRes.txt",sep=""))
	#resCats_i<-read.table(file=paste("ExplainingDR_logGeoRangeKM2-vs-tipDR_NDexp_sample100_tree",i,"_PGLS-Brownian_estimDR_allRes.txt",sep=""))
	resCats_i<-read.table(file=paste("ExplainingDR_logGeoRangeKM2-vs-tipDR_NDexp_sample100_tree",i,"_PGLS-Pagel_estimDR_allRes.txt",sep=""))

#	tree<-rep(i,length=10)
	tree<-rep(i,length=9)
	resALL[[i]]<-cbind.data.frame(resCats_i,tree)

	Terrestrial[i,] <-resALL[[i]][1,]
	Flying[i,]      <-resALL[[i]][2,]
	Arboreal[i,]    <-resALL[[i]][3,]
	Scansorial[i,]  <-resALL[[i]][4,]
	Subterranean[i,]<-resALL[[i]][5,]
	Herbivore[i,]   <-resALL[[i]][6,]
	Omnivore[i,]    <-resALL[[i]][7,]
	Carnivore[i,]   <-resALL[[i]][8,]
	ALL_mamPhy[i,]  <-resALL[[i]][9,]
}


#	Aquatic[i,]     <-resALL[[i]][6,]
#	Herbivore[i,]   <-resALL[[i]][7,]
#	Omnivore[i,]    <-resALL[[i]][8,]
#	Carnivore[i,]   <-resALL[[i]][9,]
#	ALL_mamPhy[i,]  <-resALL[[i]][10,]
#
#}

allCats<-c("Terrestrial","Flying","Arboreal","Scansorial","Subterranean"#,"Aquatic"
	,"Herbivore","Omnivore","Carnivore","ALL_mamPhy")
res<-vector("list",length(allCats))
for(i in 1:length(allCats)){
	cat<-get(allCats[i])
	mean<-apply(cat[,1:4],MARGIN=2,FUN=mean)
	res[[i]]<-rbind(mean,apply(cat[,1:4],MARGIN=2,FUN=quantile,c(0.025,0.975), na.rm=TRUE))
	rownames(res[[i]])<-c(allCats[i],"low","up")
}
resCats_sum100<-do.call(rbind,res)
colnames(resCats_sum100)<-c("a","b","pVal","lam")

#write.table(resCats_sum100,file=paste("sum100_ExplainingDR_logDR-vs-logBM_NDexp_sample100_PGLS-Brownian_estimDR_allRes.txt",sep=""))
#write.table(resCats_sum100,file=paste("sum100_ExplainingDR_logDR-vs-logBM_NDexp_sample100_PGLS-Brownian_allRes.txt",sep=""))
#write.table(resCats_sum100,file=paste("sum100_ExplainingDR_tipDR-vs-logBM_NDexp_sample100_PGLS-Brownian_estimDR_allRes.txt",sep=""))
#write.table(resCats_sum100,file=paste("sum100_ExplainingDR_tipDR-vs-logBM_NDexp_sample100_GLS_estimDR_allRes.txt",sep=""))
#write.table(resCats_sum100,file=paste("sum100_ExplainingDR_logDR-vs-logBM_NDexp_sample100_GLS_estimDR_allRes.txt",sep=""))
#write.table(resCats_sum100,file=paste("sum100_ExplainingDR_logDR-vs-logBM_NDexp_sample100_GLS_allRes.txt",sep=""))
#write.table(resCats_sum100,file=paste("sum100_ExplainingDR_tipDR-vs-logGeoRangeKM2_NDexp_sample100_GLS_estimDR_allRes.txt",sep=""))
#write.table(resCats_sum100,file=paste("sum100_ExplainingDR_tipDR-vs-logGeoRange_NDexp_sample100_PGLS-Brownian_estimDR_allRes.txt",sep=""))
#write.table(resCats_sum100,file=paste("sum100_ExplainingDR_logGeoRangeKM2-vs-tipDR_NDexp_sample100_PGLS-Brownian_estimDR_allRes.txt",sep=""))
write.table(resCats_sum100,file=paste("sum100_ExplainingDR_logGeoRangeKM2-vs-tipDR_NDexp_sample100_PGLS-Pagel_estimDR_allRes.txt",sep=""))


resCats_sum100[seq(1,27,3),]

resCats_sum100[seq(1,30,3),]



# GLS_estimDR
	# LOG tipDR
                    a        b    pVal lam
Terrestrial  -1.58651  0.00053 0.29494  NA
Flying       -2.07608  0.13094 0.00000  NA
Arboreal     -2.21023  0.10812 0.00000  NA
Scansorial   -1.85435  0.05108 0.00542  NA 
Subterranean -3.14626  0.32111 0.00000  NA
Aquatic       0.02294 -0.11000 0.00004  NA
Herbivore    -1.59178  0.01647 0.11489  NA <<
Omnivore     -1.84530  0.04851 0.00000  NA
Carnivore    -1.74752  0.01357 0.08447  NA
ALL_mamPhy   -1.72660  0.02591 0.00000  NA
	# non-log tipDR
                    a        b    pVal lam
Terrestrial   0.23406  0.00117 0.22330  NA<<
Flying        0.10724  0.03869 0.00000  NA
Arboreal      0.10200  0.02415 0.00000  NA
Scansorial    0.14467  0.01681 0.00004  NA
Subterranean -0.07171  0.07312 0.00000  NA
Aquatic       0.68672 -0.03075 0.00000  NA
Herbivore     0.23849  0.00481 0.08682  NA<<
Omnivore      0.18053  0.01044 0.00000  NA
Carnivore     0.20562  0.00386 0.03868  NA
ALL_mamPhy    0.20583  0.00668 0.00000  NA	
	# log tip DR vs log geoRange
                    a        b    pVal lam
Terrestrial  -1.48517 -0.02773 0.00017  NA
Flying       -1.48230 -0.05732 0.00000  NA
Arboreal     -1.31382 -0.06615 0.00000  NA
Scansorial   -1.44494 -0.03832 0.06112  NA
Subterranean -0.86335 -0.26965 0.00000  NA
Herbivore    -1.24689 -0.07180 0.00000  NA
Omnivore     -1.44962 -0.04453 0.00000  NA
Carnivore    -1.61433 -0.03011 0.00803  NA
ALL_mamPhy   -1.41630 -0.05181 0.00000  NA


# GLS_sameDR
                  a      b  pVal lam
Terrestrial  -1.641  0.003 0.368  NA
Flying       -2.126  0.132 0.000  NA
Arboreal     -2.254  0.109 0.000  NA
Scansorial   -1.903  0.053 0.000  NA
Subterranean -3.188  0.323 0.000  NA
Aquatic      -0.008 -0.110 0.000  NA
Herbivore    -1.633  0.018 0.001  NA
Omnivore     -1.894  0.051 0.000  NA
Carnivore    -1.805  0.016 0.000  NA
ALL_mamPhy   -1.779  0.028 0.000  NA

# Brownian_sameDR
                    a        b    pVal lam
Terrestrial  -3.13846  0.02115 0.10588   1
Flying       -2.53003  0.00641 0.17932   1
Arboreal     -2.50680 -0.00392 0.25240   1
Scansorial   -2.26920 -0.00551 0.31489   1
Subterranean -2.41125 -0.04421 0.21290   1
Aquatic      -1.95374 -0.03085 0.15180   1
Herbivore    -2.36355  0.00226 0.20137   1
Omnivore     -2.40716  0.00605 0.20289   1
Carnivore    -3.15554  0.02194 0.12257   1
ALL_mamPhy   -3.13042  0.00645 0.15424   1

# Brownian_estimDR
	# log(tip DR) ~ log(bodyMass)
                    a        b    pVal lam
Terrestrial  -3.01082  0.00086 0.73826   1
Flying       -2.57407  0.00139 0.72025   1
Arboreal     -2.54876  0.00206 0.71944   1
Scansorial   -2.28830 -0.00174 0.75941   1
Subterranean -2.58021 -0.00409 0.74341   1
Aquatic      -1.91459 -0.03330 0.06372   1 <<
Herbivore    -2.35438  0.00040 0.77719   1
Omnivore     -2.38849  0.00132 0.75376   1
Carnivore    -3.01342 -0.00250 0.63337   1
ALL_mamPhy   -3.10312 -0.00072 0.74412   1 <<<
	# tip DR ~ log(bodyMass)
                   a        b    pVal lam
Terrestrial  0.07145  0.00043 0.64190   1
Flying       0.09242 -0.00003 0.68169   1
Arboreal     0.10498  0.00010 0.72646   1
Scansorial   0.13101 -0.00041 0.70459   1
Subterranean 0.11123 -0.00483 0.53163   1
Aquatic      0.25070 -0.00948 0.09392   1 <<
Herbivore    0.12355 -0.00048 0.67823   1
Omnivore     0.11842  0.00017 0.74021   1
Carnivore    0.08147 -0.00069 0.59387   1
ALL_mamPhy   0.06735 -0.00034 0.58950   1
	# log(geoRange*(110^2)) ~ tipDR
                    a        b    pVal lam
Terrestrial  13.51566 -1.14518 0.33332   1
Flying       13.95102 -2.39151 0.27706   1
Arboreal     13.68036 -1.80449 0.38443   1
Scansorial   13.76858 -1.44563 0.59629   1
Subterranean 13.11614 -1.36663 0.50373   1
Herbivore    13.64889 -1.48749 0.24515   1
Omnivore     13.48334 -1.26807 0.37602   1
Carnivore    13.61649 -1.84568 0.33422   1
ALL_mamPhy   13.54113 -1.53555 0.08032   1
	# PAGEL lambda
	# log(geoRange*(110^2)) ~ tipDR
                    a        b    pVal       lam
Terrestrial  13.32056 -1.48557 0.00009 0.6343496
Flying       13.98207 -1.93083 0.00442 0.1836828
Arboreal     13.50570 -2.26974 0.00048 0.5110007
Scansorial   13.63550 -1.46342 0.21656 0.5997943 <<
Subterranean 13.07232 -2.30132 0.02101 0.5797682 
Herbivore    13.54746 -1.87173 0.00000 0.6150002
Omnivore     13.40481 -1.28463 0.01101 0.5636225
Carnivore    13.35321 -1.50744 0.01988 0.5221902
ALL_mamPhy   13.31217 -1.70110 0.00000 0.5870455


##
# VISUALIZE those results other than via lines... 95% CIs on the SLOPE (i.e., the EFFECT)
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_RESULTS_tipDR-vs-bodyMass_PGLS")
library(plotrix)


# (non-log / log) tip DR ~ log body mass
#========

#name<-"sum100_ExplainingDR_tipDR-vs-logBM_NDexp_sample100_GLS"
#name<-"sum100_ExplainingDR_logDR-vs-logBM_NDexp_sample100_GLS"
#name<-"sum100_ExplainingDR_tipDR-vs-logBM_NDexp_sample100_PGLS-Brownian"
#name<-"sum100_ExplainingDR_logDR-vs-logBM_NDexp_sample100_PGLS-Brownian"
#name<-"sum100_ExplainingDR_tipDR-vs-logGeoRange_NDexp_sample100_PGLS-Brownian"
#name<-"sum100_ExplainingDR_tipDR-vs-logGeoRangeKM2_NDexp_sample100_GLS"
#name<-"sum100_ExplainingDR_logGeoRangeKM2-vs-tipDR_NDexp_sample100_PGLS-Brownian"
name<-"sum100_ExplainingDR_logGeoRangeKM2-vs-tipDR_NDexp_sample100_PGLS-Pagel"

perCatRes<-read.table(file=paste(name,"_estimDR_allRes.txt",sep=""), row.names=NULL)

allCats<-c("All mammals","Terrestrial","Flying","Arboreal","Scansorial","Subterranean"#,"Aquatic"
	,"Herbivore","Omnivore","Carnivore")

means<-c(); low95<-c(); up95<-c(); pVals<-c(); lams<-c()
for(i in 1:length(allCats)){
	means[i]<-perCatRes[i*3-2,"b"]
	low95[i]<-perCatRes[i*3-1,"b"]
	up95[i]<-perCatRes[i*3,"b"]
	pVals[i]<-perCatRes[i*3-2,"pVal"]
	lams[i]<-perCatRes[i*3-2,"lam"]
}
cols<-rep("black",length(allCats))
cols[which(pVals > 0.05)]<-"red"

#reOrd_index<-c(10,1:9)
reOrd_index<-c(9,1:8)

yMax<- 1.5 #0.05 # 0.4
yMin<- -4.5 #-0.10 # -0.15 # 
yLab<-"Effect on tip Geo Range" #"Effect on tip DR"

pdf(file=paste("plotCI_",name,".pdf",sep=""),onefile=TRUE, width=6,height=4)
#jpeg(file=paste("plotCI_",name,".jpg",sep=""),units="in", res=200,width=6,height=4)

plotCI(x=1:length(allCats), y=means[reOrd_index],ui=up95[reOrd_index],li=low95[reOrd_index], ylim=c(yMin,yMax),xaxt="n", xlab="", ylab=yLab, sfrac=0,err="y", lwd=2,col="white",scol="white",pch=20,font.lab=2,cex.axis=0.95,cex.lab=1)
abline(h=0,lty=1, lwd=2, col=grey(0.6, alpha=0.5)) #rgb(1,0,0,alpha=0.5))#
plotCI(add=TRUE,x=1:length(allCats),y=means[reOrd_index],ui=up95[reOrd_index],li=low95[reOrd_index], sfrac=0, err="y", lwd=2,col=cols[reOrd_index],scol="black",pch=20,font.lab=2,cex.axis=0.95,cex.lab=1)

axis(1, at=c(1:length(allCats)), labels = FALSE)
text(x= c(1:length(allCats)), y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), cex=0.9,font=2,srt = 45, labels = allCats, xpd = TRUE, adj=c(0.95,0.05))
#mtext(text="GLS: tip DR ~ log(Body Mass)", side=3, font=2, adj=0)
#mtext(text="GLS: log(tip DR) ~ log(Body Mass)", side=3, font=2, adj=0)
#mtext(text="PGLS Brownian: tip DR ~ log(Body Mass)", side=3, font=2, adj=0)
#mtext(text="PGLS Brownian: log(tip DR) ~ log(Body Mass)", side=3, font=2, adj=0)
#mtext(text="PGLS Brownian: tip DR ~ log(Geo Range)", side=3, font=2, adj=0)
#mtext(text="GLS: tip DR ~ log(Geo Range km2)", side=3, font=2, adj=0)
#mtext(text="PGLS Brownian: log(tip Geo Range km2) ~ tip DR", side=3, font=2, adj=0)
mtext(text="PGLS Pagel: log(tip Geo Range km2) ~ tip DR", side=3, font=2, adj=0)
text(x= c(1:length(allCats)), y = yMin, cex=0.9,font=2, labels = round(lams[reOrd_index],2), xpd = TRUE)#, adj=c(0.95,0.05))

dev.off()







####
# FOR THE MCC TREE ONLY... (want on a posterior of trees, will do above)
library(phytools); library(geiger); library(nlme); library(ape)
#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_nexus-and-newickTrees")
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")

tipDataAll_wComments<-read.table(file="MamPhy_5911sp_tiplabel_DR-range-mass-herb-lifemode.txt", header=TRUE)
tipDataAll<-tipDataAll_wComments[,5:9]
rownames(tipDataAll)<-tipDataAll_wComments$tiplabel

mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")
bbone<- "NDexp" # "FBD"

##
# ALL - log(geoRange) ~ tipDR
###
dat<-na.omit(tipDataAll[,c("tipDR","geoRange")])
treeDat<-treedata(mamMCC,dat)
datPlot<-as.data.frame(treeDat$dat)

form<-log(geoRange*(110^2)) ~ tipDR

fitGLS<-gls(form, data=datPlot)
sumGLS<-summary(fitGLS)
sumGLS

fitPGLS<-gls(form, data=datPlot,correlation=corBrownian(phy=treeDat$phy))
sumPGLS<-summary(fitPGLS)
sumPGLS
	sum<-sumPGLS
	a=round(sum$coef[[1]],3)
	b=round(sum$coef[[2]],3)
	pVal=round(sum$tTable[[8]],3)
	lambda=sum$modelStruct[[1]][[1]]
	abP_All_brn<-c(a,b,pVal,lambda)

fitPGLS_lam<-gls(form, data=datPlot,correlation=corPagel(value=1, phy=treeDat$phy))
sumPGLS_lam<-summary(fitPGLS_lam)
sumPGLS_lam
	sum<-sumPGLS_lam
	a=round(sum$coef[[1]],3)
	b=round(sum$coef[[2]],3)
	pVal=round(sum$tTable[[8]],3)
	lambda=sum$modelStruct[[1]][[1]]
	abP_All_lam<-c(a,b,pVal,lambda)

# Subset out the tips per LIFEMODE
LM_cats<-c("Terrestrial","Flying","Arboreal","Scansorial","Subterranean","Aquatic")

LM_catData<-vector("list",length(LM_cats))
LM_catTips<-vector("list",length(LM_cats))
sum<-vector()
for (i in 1:length(LM_cats)){
	LM_catData[[i]]<-tipDataAll[which(tipDataAll$LIFEMODE==LM_cats[i]),]
	LM_catTips[[i]]<-rownames(LM_catData[[i]])
}

# Make the PER-LIFEMODE PGLS lines
## LIFEMODE
datPlot_LM_all<-vector("list",length(LM_catTips))
abP_LM_all<-vector("list",length(LM_catTips))

for (i in 1:length(LM_catTips)){
	if(i==6){ next } else {
	tips<-LM_catTips[[i]]
	LM_dat<-tipDataAll[match(tips,rownames(tipDataAll)),c("tipDR","geoRange")]
	dat<-na.omit(LM_dat)
	treeDat_LM<-treedata(mamMCC,dat)
	datPlot_LM<-as.data.frame(treeDat_LM$dat)
	datPlot_LM_all[[i]]<-datPlot_LM

	form<-log(geoRange*(110^2)) ~ tipDR
	#form<-log(bodyMass) ~ tipDR
	
	fitPGLS_LM<-gls(form, data=datPlot_LM,correlation=corPagel(value=1,phy=treeDat_LM$phy))
	sumPGLS_LM<-summary(fitPGLS_LM)

	sum<-sumPGLS_LM
	a=round(sum$coef[[1]],2)
	b=round(sum$coef[[2]],2)
	pVal=round(sum$tTable[[8]],2)
	lambda=sumPGLS_LM$modelStruct[[1]][[1]]

	abP_LM_all[[i]]<-c(a,b,pVal,lambda)
	}
}

# Subset out the tips of TROPHIC LEVEL
trophic_cats<-c("Herbivore","Omnivore","Carnivore")

trophicData<-vector("list",length(trophic_cats))
trophicTips<-vector("list",length(trophic_cats))
	trophicData[[1]]<-tipDataAll[which(tipDataAll$percentHerb==100),]
	trophicTips[[1]]<-rownames(trophicData[[1]])
	trophicData[[2]]<-tipDataAll[which(tipDataAll$percentHerb < 100 & tipDataAll$percentHerb > 0 ),]
	trophicTips[[2]]<-rownames(trophicData[[2]])
	trophicData[[3]]<-tipDataAll[which(tipDataAll$percentHerb==0),]
	trophicTips[[3]]<-rownames(trophicData[[3]])

## TROPHIC
datPlot_trophic_all<-vector("list",length(trophicTips))
abP_trophic_all<-vector("list",length(trophicTips))

for (i in 1:length(trophicTips)){
	tips<-trophicTips[[i]]
	dat_all<-tipDataAll[match(tips,rownames(tipDataAll)),c("tipDR","geoRange")]
	dat<-na.omit(dat_all)
	treeDat_trophic<-treedata(mamMCC,dat)
	datPlot_trophic<-as.data.frame(treeDat_trophic$dat)
	datPlot_trophic_all[[i]]<-datPlot_trophic

	form<-log(geoRange*(110^2)) ~ tipDR
	#form<-log(bodyMass) ~ tipDR
	
	fitPGLS_trophic<-gls(form, data=datPlot_trophic,correlation=corPagel(value=1,phy=treeDat_trophic$phy))
	sumPGLS_trophic<-summary(fitPGLS_trophic)

	sum<-sumPGLS_trophic
	a=round(sum$coef[[1]],2)
	b=round(sum$coef[[2]],2)
	pVal=round(sum$tTable[[8]],2)
	lambda=sum$modelStruct[[1]][[1]]

	abP_trophic_all[[i]]<-c(a,b,pVal,lambda)
}


lifemodes<-do.call(rbind,abP_LM_all)
trophic<-do.call(rbind,abP_trophic_all)

resMCC_lambda<-rbind(abP_All_brn,abP_All_lam,lifemodes,trophic)
rownames(resMCC_lambda)<-c("All_brn","All_lam",LM_cats[1:5],trophic_cats)
colnames(resMCC_lambda)<-c("a","b","pVal","lam")
	#                  a      b pVal       lam
	#All_brn      13.574 -2.694 0.00 1.0000000
	#All_lam      13.288 -1.878 0.00 0.5887112
	#Terrestrial  13.250 -1.600 0.00 0.6330396
	#Flying       13.940 -1.810 0.01 0.2096184
	#Arboreal     13.590 -3.010 0.00 0.5062276
	#Scansorial   13.640 -1.810 0.18 0.5783366
	#Subterranean 13.160 -3.030 0.00 0.5050454
	#Herbivore    13.580 -2.300 0.00 0.6185685
	#Omnivore     13.380 -1.190 0.04 0.5745476
	#Carnivore    13.320 -1.200 0.09 0.5017407

write.table(resMCC_lambda, file="MCCtree_logGeoRangeKM2-vs-tipDR_allRes.txt")

save.image(file="range_vs_tipDR_LIFEMODE_workspace.Rdata")


# test of the Subterranean niche... not actually a positive slope!! 
datSub<-datPlot_LM_all[[5]][,c("tipDR","geoRange")]
treeDat_sub<-treedata(mamMCC,datSub)

form<-log(geoRange*(110^2)) ~ tipDR
plot(form,data=datSub)
	fitGLS_sub<-gls(form, data=datSub)
	sumGLS_sub<-summary(fitGLS_sub)
	fitPGLS_sub<-gls(form, data=datSub,correlation=corPagel(value=1,phy=treeDat_sub$phy))
	sumPGLS_sub<-summary(fitPGLS_sub)

	sum<-sumPGLS_sub
	#sum<-sumGLS_sub
	a=round(sum$coef[[1]],2)
	b=round(sum$coef[[2]],2)
X=seq(min((datSub$tipDR)),max((datSub$tipDR)),length.out=20)
lines(x=X,y=b*X+a,col="black", lwd=4,lty=1) 



# PLOT
lineCol<-"black"
cexEq<-0.75
adjEq<-0.02
cexPts<-0.5
#LM_cols<-rainbow(length(LM_catTips))
LM_cols<-c("firebrick2","deepskyblue2","chartreuse3","darkorchid3","darkgoldenrod3","black")
trophic_cols<-c("sienna4","sienna3","sienna1")

pointCol<-rep(grey(0.3, alpha=0.5),length(datPlot[,1]))
pchs<-rep(16,length(datPlot[,1]))

#pdf(file="tipDR-vs-geoRange_Overall_withLMcats_andTrophic.pdf", onefile=TRUE, width=4.5,height=4) #5.333333,height=5) #
#pdf(file="tipDR-vs-geoRange_Overall_withLMcats.pdf", onefile=TRUE, width=4.5,height=4) #5.333333,height=5) #
pdf(file="tipDR-vs-geoRange_Overall_trophicOnly.pdf", onefile=TRUE, width=4.5,height=4) #5.333333,height=5) #

#quartz(width=4.5,height=4) #5.333333,height=5) #
op <- par(oma = c(2,4,1,1) + 0.1,		#c(bottom, left, top, right)
          mar = c(1,0,0,0) + 0.1, cex=0.85)
#layout(matrix(1:2,1,2), heights=c(1,1), widths=c(2/3,1/3))
#layout(matrix(1:2,2,1), heights=c(0.5,0.5), widths=c(1,1))

plot(log(geoRange*(110^2)) ~ tipDR, data=datPlot, type="n",col=pointCol, yaxt="n",ylab="" ,xlab="", cex=cexPts, pch=pchs)
mtext(side=2, text=bquote(bold(paste("Geographic range (10"^4," km"^"2",")",sep=""))),cex=1,font=2,line=2.5)
axis(side=2,at=log(c(10^4,10^5,10^6,10^7,10^8)), labels=FALSE)
text(xpd=NA, x=rep(-0.08,5), y=log(c(10^4,10^5,10^6,10^7,10^8)), labels=c(1,10,100,expression(10^3),expression(10^4))) 

mtext(side=1, text="Tip DR (species/Ma)",cex=1,font=2,line=2.5)

sum<-sumPGLS_lam
a=round(sum$coef[[1]],2)
b=round(sum$coef[[2]],2)
pVal=round(sum$tTable[[8]],3)
lam=round(sum$modelStruct[[1]][[1]],2)
X=seq(min((datPlot$tipDR)),max((datPlot$tipDR)),length.out=20)
lines(x=X,y=b*X+a,col=lineCol, lwd=4,lty=1) 
#mtext(text=bquote("y = "*.(a)*" + "*.(b)*"x; P < 0.001"), side=3, col="black", cex=cexEq, adj=adjEq, line=-1.5)
mtext(text=bquote("All: y = "*.(a)*" + "*.(b)*"x; P = "*.(pVal)*"; lam = "*.(lam)), side=3, col="black", cex=cexEq, adj=adjEq, line=-1.5)
	Xlab<-max((datPlot$tipDR))
	Ylab<-b*Xlab+a
	#text(x=Xlab, y=Ylab, labels="all",font=2, cex=1.5, col=grey(0))#shiftCols[i]) 

#for (i in 1:length(LM_catTips)){
#	if(i==6){ next } else {
#	
#	dat<-datPlot_LM_all[[i]]
#	#points(log(geoRange*(110^2)) ~ tipDR, data=dat, col=LM_cols[i], cex=cexPts, pch=pchs)
#
#	X=seq(min((dat$tipDR)),max((dat$tipDR)),length.out=20)
#	mtext(text=bquote(.(LM_cats[i])*": y = "*.(abP_LM_all[[i]][1])*" + "*.(abP_LM_all[[i]][2])*"x; P = "*.(abP_LM_all[[i]][3])*"; lam = "*.(round(abP_LM_all[[i]][4],2))), side=3, col=LM_cols[i], cex=cexEq, adj=adjEq, line=-1.5-i)
#
#	if(abP_LM_all[[i]][3] > 0.05){ 
#	lines(x=X,y=abP_LM_all[[i]][2]*X+abP_LM_all[[i]][1],col=LM_cols[i], lwd=3,lty=2) 
#	} else {
#	lines(x=X,y=abP_LM_all[[i]][2]*X+abP_LM_all[[i]][1],col=LM_cols[i], lwd=3,lty=1)
#	}	
#	}
#}

for (i in 1:length(datPlot_trophic_all)){
	dat<-datPlot_trophic_all[[i]]

	X=seq(min((dat$tipDR)),max((dat$tipDR)),length.out=20)
	#mtext(text=bquote(.(trophic_cats[i])*": y = "*.(abP_trophic_all[[i]][1])*" + "*.(abP_trophic_all[[i]][2])*"x; P = "*.(abP_trophic_all[[i]][3])*"; lam = "*.(round(abP_LM_all[[i]][4],2))), side=3, col=trophic_cols[i], cex=cexEq, adj=adjEq, line=-6.5-i)
	mtext(text=bquote(.(trophic_cats[i])*": y = "*.(abP_trophic_all[[i]][1])*" + "*.(abP_trophic_all[[i]][2])*"x; P = "*.(abP_trophic_all[[i]][3])*"; lam = "*.(round(abP_LM_all[[i]][4],2))), side=3, col=trophic_cols[i], cex=cexEq, adj=adjEq, line=-1.5-i)

	if(abP_trophic_all[[i]][3] > 0.05){  
 	lines(x=X,y=abP_trophic_all[[i]][2]*X+abP_trophic_all[[i]][1],col=trophic_cols[i], lwd=3,lty=2)
	} else {
	lines(x=X,y=abP_trophic_all[[i]][2]*X+abP_trophic_all[[i]][1],col=trophic_cols[i], lwd=3,lty=1)
	}	
}


dev.off()






#############

# percentCarn ~ tipDR >>> Never actually completed this one...
dat<-na.omit(tipDataAllHerb[,c(1,4)])
#dat<-na.omit(shiftedDat[,1:2])
treeDat<-treedata(mamMCC,dat)
datPlot<-as.data.frame(treeDat$dat)

percentCarn<-100-datPlot$percentHerb

form<- tipDR ~ percentCarn

fitGLS<-gls(form, data=datPlot)
sumGLS<-summary(fitGLS)
sumGLS

plot(form, data=datPlot)

fitPGLS<-gls(form, data=datPlot,correlation=corBrownian(phy=treeDat$phy))
sumPGLS<-summary(fitPGLS)
sumPGLS


# log(geoRange) ~ tipDR
dat<-na.omit(tipDataAll[,1:2])
#dat<-na.omit(shiftedDat[,1:2])
treeDat<-treedata(mamMCC,dat)
datPlot<-as.data.frame(treeDat$dat)

form<-log(geoRange*(110^2)) ~ tipDR

fitGLS<-gls(form, data=datPlot)
sumGLS<-summary(fitGLS)
sumGLS

fitPGLS<-gls(form, data=datPlot,correlation=corBrownian(phy=treeDat$phy))
sumPGLS<-summary(fitPGLS)
sumPGLS

#fitPGLM<-phyloglm(geoRange ~ tipDR, data=datPlot,phy=treeDat$phy, method="poisson_GEE")#, start.beta=c(0.1,0.2))

# log(bodyMass) ~ tipDR
dat2<-na.omit(tipDataAll[,c(1,3)])
treeDat2<-treedata(mamMCC,dat2)
datPlot2<-as.data.frame(treeDat2$dat)

form2<-log(bodyMass) ~ tipDR

fitGLS2<-gls(form2, data=datPlot2)
sumGLS2<-summary(fitGLS2)
sumGLS2

fitPGLS2<-gls(form2, data=datPlot2,correlation=corBrownian(phy=treeDat2$phy))
sumPGLS2<-summary(fitPGLS2)
sumPGLS2

# log(geoRange) ~ log(bodyMass)
dat3<-na.omit(tipDataAll[,2:3])
treeDat3<-treedata(mamMCC,dat3)
datPlot3<-as.data.frame(treeDat3$dat) # 5051 values

form3<-log(geoRange*(110^2)) ~ log(bodyMass)

fitGLS3<-gls(form3, data=datPlot3)
sumGLS3<-summary(fitGLS3)
sumGLS3

fitPGLS3<-gls(form3, data=datPlot3,correlation=corBrownian(phy=treeDat3$phy))
sumPGLS3<-summary(fitPGLS3)
sumPGLS3

save.image(file="rangeMass_vs_tipDR_workspace.Rdata")

# load
library(geiger); library(nlme)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
load(file="rangeMass_vs_tipDR_workspace.Rdata")
# Now PLOT-- dual

## get COLORS
# rate shifts greater than 2x
B 2.3, 2.7 - Macropod.-Potoridae
J 3.2 - Pteropus
Q 4.0 - Ctenomyidae
W 2.3 - Rattus-Srilankamys

subList<-vector("list",length=4)
greaterThan2x<-c(2,10,17,23)
for (i in greaterThan2x){
	subList[[i]]<-cladeTips_all[[i]]
}
cladeTips_gr2x<-unlist(subList)

# make LINES just for those 2x rate-shifted taxa...
shifted2x_dat<-tipDataAll[match(cladeTips_gr2x,rownames(tipDataAll)),1:2]
dat<-na.omit(shifted2x_dat)
treeDat_2x<-treedata(mamMCC,dat)
datPlot_2x<-as.data.frame(treeDat_2x$dat)
form<-log(geoRange*(110^2)) ~ tipDR

fitGLS_2x<-gls(form, data=datPlot_2x)
sumGLS_2x<-summary(fitGLS_2x)
sumGLS_2x

fitPGLS_2x<-gls(form, data=datPlot_2x,correlation=corBrownian(phy=treeDat_2x$phy))
sumPGLS_2x<-summary(fitPGLS_2x)
sumPGLS_2x

shifted2x_dat<-tipDataAll[match(cladeTips_gr2x,rownames(tipDataAll)),c(1,3)]
dat<-na.omit(shifted2x_dat)
treeDat2_2x<-treedata(mamMCC,dat)
datPlot2_2x<-as.data.frame(treeDat2_2x$dat)
form<-log(bodyMass) ~ tipDR

fitGLS2_2x<-gls(form, data=datPlot2_2x)
sumGLS2_2x<-summary(fitGLS2_2x)
sumGLS2_2x

fitPGLS2_2x<-gls(form, data=datPlot2_2x,correlation=corBrownian(phy=treeDat2_2x$phy))
sumPGLS2_2x<-summary(fitPGLS2_2x)
sumPGLS2_2x

## Just the GeoRanges...
library(viridis)

lineCol<-"black"
cexEq<-0.75
adjEq<-0.98
cexPts<-0.65
shiftCols<-rainbow(length(cladeTips_all))

#datPlot_shift_all<-vector("list",length(cladeTips_all))
#abP<-vector("list",length(cladeTips_all))
datPlot_shift_all_2<-vector("list",length(cladeTips_all))
abP_2<-vector("list",length(cladeTips_all))

for (i in 1:length(cladeTips_all)){
#	if(i==6 || i==3){ next } else {
	if(i==3){ next } else {
	tips<-cladeTips_all[[i]]
	shifted_dat<-tipDataAll[match(tips,rownames(tipDataAll)),c(1,3)]
	dat<-na.omit(shifted_dat)
	treeDat_shift<-treedata(mamMCC,dat)
	datPlot_shift<-as.data.frame(treeDat_shift$dat)
	datPlot_shift_all_2[[i]]<-datPlot_shift

	#form<-log(geoRange*(110^2)) ~ tipDR
	form<-log(bodyMass) ~ tipDR
	
	fitPGLS_shift<-gls(form, data=datPlot_shift,correlation=corPagel(value=1,phy=treeDat_shift$phy))
	sumPGLS_shift<-summary(fitPGLS_shift)

	sum<-sumPGLS_shift
	a=round(sum$coef[[1]],2)
	b=round(sum$coef[[2]],2)
	pVal=round(sum$tTable[[8]],2)
	abP_2[[i]]<-c(a,b,pVal)
	}
}
save.image(file="rangeMass_vs_tipDR_workspace.Rdata")

load(file="rangeMass_vs_tipDR_workspace.Rdata")

#pdf(file="tipDR-vs-Traits_geoRange_ONLY_sepRateShiftSlopes.pdf", onefile=TRUE, width=4.5,height=4) #5.333333,height=5) #
pdf(file="tipDR-vs-Traits_bodyMass_ONLY_sepRateShiftSlopes.pdf", onefile=TRUE, width=4.5,height=4) #5.333333,height=5) #

#quartz(width=4.5,height=4) #5.333333,height=5) #
op <- par(oma = c(2,4,1,1) + 0.1,		#c(bottom, left, top, right)
          mar = c(1,0,0,0) + 0.1, cex=0.85)
#layout(matrix(1:2,1,2), heights=c(1,1), widths=c(2/3,1/3))
#layout(matrix(1:2,2,1), heights=c(0.5,0.5), widths=c(1,1))

pointCol<-rep(grey(0.3, alpha=0.3),length(datPlot2[,1]))
#pointCol[match(cladeTips_gr2x,rownames(datPlot))]<-rgb(1,0,0,alpha=1)
pchs<-rep(1,length(datPlot2[,1]))
#pchs[match(cladeTips_gr2x,rownames(datPlot))]<-16

#plot(log(geoRange*(110^2)) ~ tipDR, data=datPlot, col=pointCol, yaxt="n",ylab="" ,xlab="", cex=cexPts, pch=pchs)
#mtext(side=2, text=bquote(bold(paste("Geographic range (10"^4," km"^"2",")",sep=""))),cex=1,font=2,line=2.5)
#axis(side=2,at=log(c(10^4,10^5,10^6,10^7,10^8)), labels=FALSE)
#text(xpd=NA, x=rep(-0.08,5), y=log(c(10^4,10^5,10^6,10^7,10^8)), labels=c(1,10,100,expression(10^3),expression(10^4))) 
plot(log(bodyMass) ~ tipDR, data=datPlot2, col=pointCol, yaxt="n",ylab="" ,xlab="", cex=cexPts, pch=pchs)
mtext(side=2, text="Body mass (g)",cex=1,font=2,line=2.5)
axis(side=2,at=log(c(1,100,10000,10^6,10^8)), labels=FALSE)
text(xpd=NA, x=rep(-0.08,5), y=log(c(1,100,10000,10^6,10^8)), labels=c(1,100,expression(10^4),expression(10^6),expression(10^8)))
#axis(side=1,labels=FALSE)
mtext(side=1, text="Tip DR (species/Ma)",cex=1,font=2,line=2.5)

sum<-sumPGLS2
a=round(sum$coef[[1]],2)
b=round(sum$coef[[2]],2)
pVal=round(sum$tTable[[8]],3)
X=seq(min((datPlot2$tipDR)),max((datPlot2$tipDR)),length.out=20)
lines(x=X,y=b*X+a,col=lineCol, lwd=4,lty=1) 
#mtext(text=bquote("y = "*.(a)*" + "*.(b)*"x; P < 0.001"), side=3, col="black", cex=cexEq, adj=adjEq, line=-1.5)
mtext(text=bquote("y = "*.(a)*" + "*.(b)*"x; P = "*.(pVal)), side=3, col="black", cex=cexEq, adj=adjEq, line=-1.5)
	Xlab<-max((datPlot2$tipDR))
	Ylab<-b*Xlab+a
	text(x=Xlab, y=Ylab, labels="all",font=2, cex=1.5, col=grey(0))#shiftCols[i]) 


for (i in 1:length(cladeTips_all)){
#	if(i==6 || i==3){ next } else {
	if(i==3){ next } else {

	if(abP_2[[i]][3] > 0.05){ next } else {
	dat<-datPlot_shift_all_2[[i]]
	#points(log(geoRange*(110^2)) ~ tipDR, data=dat, col=shiftCols[i], cex=cexPts, pch=pchs)

	X=seq(min((dat$tipDR)),max((dat$tipDR)),length.out=20)
	lines(x=X,y=abP_2[[i]][2]*X+abP_2[[i]][1],col=shiftCols[i], lwd=2,lty=1)
	Xlab<-max((dat$tipDR))
	Ylab<-abP_2[[i]][2]*Xlab+(abP_2[[i]][1])
	text(x=Xlab, y=Ylab, labels=cladeNames[i],font=2, cex=1.5, col=grey(0))#shiftCols[i]) 
		}
	}
}

dev.off()




## DO THE CO-PLOT...
lineCol<-"black"
cexEq<-0.75
adjEq<-0.98
cexPts<-0.65

pdf(file="tipDR-vs-Traits_geoRange_bodyMass_coPlot_gr2x_fixedAxis.pdf", onefile=TRUE, width=5,height=7) #5.333333,height=5) #

quartz(width=4,height=7) #5.333333,height=5) #
op <- par(oma = c(2,4,1,1) + 0.1,		#c(bottom, left, top, right)
          mar = c(1,0,0,0) + 0.1, cex=0.85)
#layout(matrix(1:2,1,2), heights=c(1,1), widths=c(2/3,1/3))
layout(matrix(1:2,2,1), heights=c(0.5,0.5), widths=c(1,1))

pointCol<-rep(grey(0.3, alpha=0.3),length(datPlot[,1]))
pointCol[match(cladeTips_gr2x,rownames(datPlot))]<-rgb(1,0,0,alpha=1)
pchs<-rep(1,length(datPlot[,1]))
pchs[match(cladeTips_gr2x,rownames(datPlot))]<-16

plot(log(geoRange*(110^2)) ~ tipDR, data=datPlot, col=pointCol, xaxt="n", yaxt="n",ylab="" ,xlab="", cex=cexPts, pch=pchs)
axis(side=1,labels=FALSE)
axis(side=2,at=log(c(10^4,10^5,10^6,10^7,10^8)), labels=FALSE)
text(xpd=NA, x=rep(-0.08,5), y=log(c(10^4,10^5,10^6,10^7,10^8)), labels=c(1,10,100,expression(10^3),expression(10^4))) 
mtext(side=2, text=bquote(bold(paste("Geographic range (10"^4," km"^"2",")",sep=""))),cex=1,font=2,line=2.5)

sum<-sumGLS
a=round(sum$coef[[1]],2)
b=round(sum$coef[[2]],2)
pVal=round(sum$tTable[[8]],6)
X=seq(min((datPlot$tipDR)),max((datPlot$tipDR)),length.out=20)
lines(x=X,y=b*X+a,col=lineCol, lwd=2,lty=2) 
mtext(text=bquote("GLS: y = "*.(a)*" + "*.(b)*"x; P < 0.001"), side=3, col="black", cex=cexEq, adj=adjEq, line=-2)

sum<-sumPGLS
a=round(sum$coef[[1]],2)
b=round(sum$coef[[2]],2)
pVal=round(sum$tTable[[8]],6)
X=seq(min((datPlot$tipDR)),max((datPlot$tipDR)),length.out=20)
lines(x=X,y=b*X+a,col=lineCol, lwd=2,lty=1) 
mtext(text=bquote("PGLS: y = "*.(a)*" + "*.(b)*"x; P < 0.001"), side=3, col="black", cex=cexEq, adj=adjEq, line=-1)

#sum<-sumGLS_2x
#a=round(sum$coef[[1]],2)
#b=round(sum$coef[[2]],2)
#pVal=round(sum$tTable[[8]],2)
#X=seq(min((datPlot_2x$tipDR)),max((datPlot_2x$tipDR)),length.out=20)
#lines(x=X,y=b*X+a,col="red", lwd=2,lty=2) 
#mtext(text=bquote("GLS: y = "*.(a)*" + "*.(b)*"x; P = "*.(pVal)), side=3, col="red", cex=cexEq, adj=adjEq, line=-4)
#
#sum<-sumPGLS_2x
#a=round(sum$coef[[1]],2)
#b=round(sum$coef[[2]],2)
#pVal=round(sum$tTable[[8]],2)
#X=seq(min((datPlot_2x$tipDR)),max((datPlot_2x$tipDR)),length.out=20)
#lines(x=X,y=b*X+a,col="red", lwd=2,lty=1) 
#mtext(text=bquote("PGLS: y = "*.(a)*" + "*.(b)*"x; P = "*.(pVal)), side=3, col="red", cex=cexEq, adj=adjEq, line=-3)



pointCol<-rep(grey(0.3, alpha=0.3),length(datPlot2[,1]))
pointCol[match(cladeTips_gr2x,rownames(datPlot2))]<-rgb(1,0,0,alpha=1)
pchs<-rep(1,length(datPlot2[,1]))
pchs[match(cladeTips_gr2x,rownames(datPlot2))]<-16

plot(log(bodyMass) ~ tipDR, data=datPlot2, col=pointCol, yaxt="n",ylab="" ,xlab="", cex=cexPts, pch=pchs)
axis(side=2,at=log(c(1,100,10000,10^6,10^8)), labels=FALSE)
text(xpd=NA, x=rep(-0.08,5), y=log(c(1,100,10000,10^6,10^8)), labels=c(1,100,expression(10^4),expression(10^6),expression(10^8)))
mtext(side=2,text="Body mass (g)",cex=1, font=2,line=2.5)

sum<-sumGLS2
a=round(sum$coef[[1]],2)
b=round(sum$coef[[2]],2)
pVal=round(sum$tTable[[8]],6)
X=seq(min((datPlot2$tipDR)),max((datPlot2$tipDR)),length.out=20)
lines(x=X,y=b*X+a,col=lineCol, lwd=2,lty=2) 
mtext(text=bquote("GLS: y = "*.(a)*" + "*.(b)*"x; P < 0.001"), side=3, col="black", cex=cexEq, adj=adjEq, line=-2)

sum<-sumPGLS2
a=round(sum$coef[[1]],2)
b=round(sum$coef[[2]],2)
pVal=round(sum$tTable[[8]],2)
X=seq(min((datPlot2$tipDR)),max((datPlot2$tipDR)),length.out=20)
lines(x=X,y=b*X+a,col=lineCol, lwd=2,lty=1) 
mtext(text=bquote("PGLS: y = "*.(a)*" + "*.(b)*"x; P = "*.(pVal)), side=3, col="black", cex=cexEq, adj=adjEq, line=-1)

#sum<-sumGLS2_2x
#a=round(sum$coef[[1]],2)
#b=round(sum$coef[[2]],2)
#pVal=round(sum$tTable[[8]],2)
#X=seq(min((datPlot2_2x$tipDR)),max((datPlot2_2x$tipDR)),length.out=20)
#lines(x=X,y=b*X+a,col="red", lwd=2,lty=2) 
#mtext(text=bquote("GLS: y = "*.(a)*" + "*.(b)*"x; P = "*.(pVal)), side=3, col="red", cex=cexEq, adj=adjEq, line=-4)
#
#sum<-sumPGLS2_2x
#a=round(sum$coef[[1]],2)
#b=round(sum$coef[[2]],2)
#pVal=round(sum$tTable[[8]],2)
#X=seq(min((datPlot2_2x$tipDR)),max((datPlot2_2x$tipDR)),length.out=20)
#lines(x=X,y=b*X+a,col="red", lwd=2,lty=1) 
#mtext(text=bquote("PGLS: y = "*.(a)*" + "*.(b)*"x; P = "*.(pVal)), side=3, col="red", cex=cexEq, adj=adjEq, line=-3)

title(main="", xlab = "Tip DR (species/Ma)",
      ylab = "",
      outer = TRUE, line = 1,cex.main=1.5,font.main=1,cex.axis=1,cex.lab=1,font.axis=2, font.lab=2)

dev.off()


########
# log(geoRange) vs log(bodyMass)... hmmm...
pdf(file="Traits_geoRange_bodyMass_coPlot_gr2x_fixedAxis.pdf", onefile=TRUE, width=5,height=5) #5.333333,height=5) #

pointCol<-rep(grey(0.3, alpha=0.3),length(datPlot3[,1]))
pointCol[match(cladeTips_gr2x,rownames(datPlot3))]<-rgb(1,0,0,alpha=1)
pchs<-rep(1,length(datPlot3[,1]))
pchs[match(cladeTips_gr2x,rownames(datPlot3))]<-16

plot(log(geoRange*(110^2)) ~ log(bodyMass), data=datPlot3, col=pointCol, xaxt="n",yaxt="n",ylab="" ,xlab="", cex=cexPts, pch=pchs)

axis(side=2,at=log(c(10^4,10^5,10^6,10^7,10^8)), labels=c(1,10,100,expression(10^3),expression(10^4)))
mtext(side=2, text=bquote(bold(paste("Geographic range (10"^4," km"^"2",")",sep=""))),cex=1,font=2,line=2.5)

axis(side=1,at=log(c(1,100,10000,10^6,10^8)), labels=c(1,100,expression(10^4),expression(10^6),expression(10^8)))
mtext(side=1,text="Body mass (g)",cex=1, font=2,line=2.5)

sum<-sumGLS3
a=round(sum$coef[[1]],2)
b=round(sum$coef[[2]],2)
pVal=round(sum$tTable[[8]],6)
X=seq(min(log(datPlot3$bodyMass)),max(log(datPlot3$bodyMass)),length.out=20)
lines(x=X,y=b*X+a,col=lineCol, lwd=2,lty=2) 
mtext(text=bquote("GLS: y = "*.(a)*" + "*.(b)*"x; P < 0.001"), side=3, col="black", cex=cexEq, adj=adjEq, line=-2)

sum<-sumPGLS3
a=round(sum$coef[[1]],2)
b=round(sum$coef[[2]],2)
pVal=round(sum$tTable[[8]],2)
X=seq(min(log(datPlot3$bodyMass)),max(log(datPlot3$bodyMass)),length.out=20)
lines(x=X,y=b*X+a,col=lineCol, lwd=2,lty=1) 
mtext(text=bquote("PGLS: y = "*.(a)*" + "*.(b)*"x; P = "*.(pVal)), side=3, col="black", cex=cexEq, adj=adjEq, line=-1)

dev.off()

# tipDR ~ log(geoRange) 
fitPagel<-gls(tipDR ~ log(geoRange), data=datPlot, correlation=corPagel(value=0.5,phy=treeDat$phy))
sumPagel<-summary(fitPagel)
sumPagel$coef
# corPagel # P = 0.0003237138 # a = 0.07576442 # b = -0.00136883 # lam = 0.9895197

fitBrownian2<-gls(tipDR ~ log(bodyMass), data=datPlot, correlation=corBrownian(phy=treeDat$phy))
sumBrownian2<-summary(fitBrownian)
# corBrownian # P = 0.6393165 # a = 0.0655446449 # b = 0.0004441946 





#########
# ==============
# RATE-SHIFTED clades --24 of them--and identify what's going on with them in terms of COVARIATES
#####
# Clade DivRate
library(ape); library(phytools)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"

tipDataAll<-read.table("MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)
head(tipDataAll)

mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")

agesAllClades_origALL<-read.table(file=paste("divTime_24BAMMnodes_",bbone,"_MCC_target_DETAILS.txt", sep=""), header=TRUE)
firstShift<-c(1,3,5,7,8,10,11,13,15,17:29,31,32) #5 == rm Placentalia
agesAllClades_ALL<-agesAllClades_origALL[firstShift,]
#	agesAllClades_ALL<-agesAllClades_origALL[which(agesAllClades_origALL$Indep==1),]

nodeNums<-agesAllClades_ALL$nodeMCC
#nodeNumsAll<-agesAllClades_ALL$nodeMCC
#firstShift<-c(1,3,5,7,8,10,11,13,15,17:29,31,32) #5 == rm Placentalia
#nodeNums<-nodeNumsAll[firstShift]
#shiftDivDat1<-agesAllClades_ALL[firstShift,c("CLADE_label","ID_label","nodeMCC","numTaxa","avgFactor","avgIncDec","avgCladeDiv","cladeDiv_low95","cladeDiv_high95","mean","lower","upper")]
#colnames(shiftDivDat1)<-c("clade","ID","nodeMCC","numTaxa","factor","shift","cladeDiv","cladeDiv_low","cladeDiv_high","mean","lower","upper")
	
# gather the BAMM rate-shift data... 
#	shiftDivDat<-agesAllClades_ALL[order(agesAllClades_ALL$avgCladeDiv),]
#	nodeNums<-shiftDivDat$nodeMCC # re-ordered
shiftDivDat<-agesAllClades_ALL # keeping original phylo order

shiftedNames<-shiftDivDat$CLADE_label
shiftedIDs<-shiftDivDat$ID

# Use the MCC tree to identify tips that SUBTEND each clade:
tipNames<-vector("list",length(nodeNums))
lengths<-c()
for(i in 1:length(nodeNums)){
	tipNums<-getDescendants(mamMCC,node=nodeNums[i])
	tipNames[[i]]<-as.vector(na.omit(mamMCC$tip.label[tipNums]))
	lengths[i]<-length(tipNames[[i]])
}
#checkLengths<-cbind.data.frame(lengths,agesAllClades_ALL$numTaxa[firstShift]) # okay! >> but not after you REORDER things
checkLengths<-cbind.data.frame(lengths,shiftDivDat$numTaxa) # okay!

# aaand the non-clade tips!!
nonClade_tipNames<-vector("list",length(nodeNums))
lengths<-c()
for(i in 1:length(nodeNums)){
	nonClade_tipNames[[i]]<-as.vector(setdiff(mamMCC$tip.label,tipNames[[i]]))
	lengths[i]<-length(nonClade_tipNames[[i]])
}
checkLengths_nonClade<-cbind.data.frame(lengths,checkLengths$lengths) # okay!
rowSums(checkLengths_nonClade) # groovy.


## TRAITS to look at
#traitVars1<-c("tipDR_mean10k","BM_final_kg","geoRange_km2","homeRange_km2_ext","DispDistAll_km_ext","GenerationLength_d","trophic123", "activity123","lifemode1234") # "perVert","perInvert","perPlant",
traitVars1<-c("BM_final_kg","homeRange_km2_ext","geoRange_km2","DispDistAll_km_ext","GenerationLength_d","lifemode1234") # "perVert","perInvert","perPlant",
traitDat1<-tipDataAll[,c("tiplabel",traitVars1)]

#traitDat<-cbind.data.frame(traitDat1$tiplabel,traitDat1$tipDR_mean10k,log(traitDat1$BM_final_kg),log(traitDat1$geoRange_km2),log(traitDat1$homeRange_km2_ext),log(traitDat1$DispDistAll_km_ext),log(traitDat1$GenerationLength_d),traitDat1$trophic123, traitDat1$activity123,traitDat1$lifemode1234) # traitDat1$perVert,traitDat1$perInvert,traitDat1$perPlant, 
#traitVars<-c("tipDR_mean10k","logBM_kg","logGR_km2","logHR_km2","logDispDin_km","logGenLen_d","trophic123", "activity123","lifemode1234") #"perVert","perInvert","perPlant",
#traitDat<-cbind.data.frame(traitDat1$tiplabel,log(traitDat1$BM_final_kg),log(traitDat1$homeRange_km2_ext),log(traitDat1$geoRange_km2),log(traitDat1$DispDistAll_km_ext),(traitDat1$lifemode1234*log(traitDat1$DispDistAll_km_ext)),log(traitDat1$GenerationLength_d)) # traitDat1$perVert,traitDat1$perInvert,traitDat1$perPlant, 
#traitVars<-c("Body mass (log kg)","Home range (log km2)","Geographic range (log km2)","Dispersal distance (log km)","Dispersal index (log km x lifemode)","Generation Length (log d)") #"perVert","perInvert","perPlant",
unscaled<-cbind.data.frame(log(traitDat1$BM_final_kg),log(traitDat1$homeRange_km2_ext),log(traitDat1$geoRange_km2),(traitDat1$lifemode1234*log(traitDat1$DispDistAll_km_ext))) # traitDat1$perVert,traitDat1$perInvert,traitDat1$perPlant, 
traitDat<-cbind.data.frame(traitDat1$tiplabel,scale(unscaled, center=TRUE,scale=TRUE))

traitVars<-c("Body mass (log kg)","Home range (log km2)","Geographic range (log km2)","Dispersal index (log km x lifemode)") #"perVert","perInvert","perPlant",
colnames(traitDat)<-c("tiplabel",traitVars)

## standard errors
sterror <- function(x){
			se <- sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)]))
			return(se)
			}

# shifted clade summary
medianLowHigh_CLADE<-vector("list",length(traitVars))
for(i in 1:length(traitVars)){	# across 9 vars
	tipDat_ALL<-vector("list",length(tipNames))
	for(j in 1:length(tipNames)){ # for 18 indep rate-shifted clades
		tipDat<-traitDat[match(tipNames[[j]],traitDat$tiplabel),traitVars[i]]
		names(tipDat)<-tipNames[[j]]
		tipDat_ALL[[j]]<-tipDat
		}
#	medianLowHigh<-cbind.data.frame(unlist(lapply(tipDat_ALL,length)),unlist(lapply(tipDat_ALL,median,na.rm=TRUE)),as.vector(unlist(lapply(tipDat_ALL,quantile,c(0.025),na.rm=TRUE))),as.vector(unlist(lapply(tipDat_ALL,quantile,c(0.975),na.rm=TRUE))))
	medianLowHigh<-cbind.data.frame(unlist(lapply(tipDat_ALL,length)),unlist(lapply(tipDat_ALL,median,na.rm=TRUE)),(unlist(lapply(tipDat_ALL,median,na.rm=TRUE))-unlist(lapply(tipDat_ALL,sterror))),(unlist(lapply(tipDat_ALL,median,na.rm=TRUE))+unlist(lapply(tipDat_ALL,sterror))))
	rownames(medianLowHigh)<-shiftedIDs
	colnames(medianLowHigh)<-c(traitVars[i],"median","low","high")
medianLowHigh_CLADE[[i]]<-medianLowHigh
}

# non-shifted (background) clade summary
medianLowHigh_nonCLADE<-vector("list",length(traitVars))
for(i in 1:length(traitVars)){	# across 9 vars
	tipDat_ALL<-vector("list",length(nonClade_tipNames))
	for(j in 1:length(nonClade_tipNames)){ # for 18 indep rate-shifted clades
		tipDat<-traitDat[match(nonClade_tipNames[[j]],traitDat$tiplabel),traitVars[i]]
		names(tipDat)<-nonClade_tipNames[[j]]
		tipDat_ALL[[j]]<-tipDat
		}
#	medianLowHigh<-cbind.data.frame(unlist(lapply(tipDat_ALL,length)),unlist(lapply(tipDat_ALL,median,na.rm=TRUE)),as.vector(unlist(lapply(tipDat_ALL,quantile,c(0.025),na.rm=TRUE))),as.vector(unlist(lapply(tipDat_ALL,quantile,c(0.975),na.rm=TRUE))))
	medianLowHigh<-cbind.data.frame(unlist(lapply(tipDat_ALL,length)),unlist(lapply(tipDat_ALL,median,na.rm=TRUE)),(unlist(lapply(tipDat_ALL,median,na.rm=TRUE))-unlist(lapply(tipDat_ALL,sterror))),(unlist(lapply(tipDat_ALL,median,na.rm=TRUE))+unlist(lapply(tipDat_ALL,sterror))))
	rownames(medianLowHigh)<-shiftedIDs
	colnames(medianLowHigh)<-c(traitVars[i],"median","low","high")
medianLowHigh_nonCLADE[[i]]<-medianLowHigh
}

# >> BUT-- Standard errors end up being more a function of SAMPLE SIZE...

# LOOK AT DIFF IN MEDIANS-- WILCOXON / MANN-WHITNEY U
###
diffInMedians_ALL<-vector("list",length(traitVars))
for(i in 1:length(traitVars)){	# across 4 vars
	meanDat<-rep(NA,length(tipNames)); lowDat<-rep(NA,length(tipNames)); upDat<-rep(NA,length(tipNames)); pVals<-rep(NA,length(tipNames)); overlapZero<-rep(NA,length(tipNames))
	for(j in 1:length(tipNames)){ # for 18 indep rate-shifted clades
		if(i==3 & j==7 | i==4 & j==7 ) { next } else { # excludes Cetacea for geoRange and dispersal (need IUCN shape data)

		cladeDat<-traitDat[match(tipNames[[j]],traitDat$tiplabel),traitVars[i]]
		names(cladeDat)<-tipNames[[j]]

		nonCladeDat<-traitDat[match(nonClade_tipNames[[j]],traitDat$tiplabel),traitVars[i]]
		names(nonCladeDat)<-nonClade_tipNames[[j]]

		test<-wilcox.test(x=cladeDat, y=nonCladeDat, alternative="two.sided", conf.int=TRUE,conf.level=0.95)

		meanDat[j]<-test$estimate[[1]]
		lowDat[j]<-test$conf.int[1]
		upDat[j]<-test$conf.int[2]
		pVals[j]<-test$p.value
		overlapZero[j]<-if (test$conf.int[1] < 0 && test$conf.int[2] < 0 ) {"bigger"} else if (test$conf.int[1] > 0 && test$conf.int[2] > 0) {"smaller"} else {"NS"}
		}
	}
	diffInMedians_ALL[[i]]<-cbind.data.frame(meanDat,lowDat,upDat,pVals,overlapZero)
}

#shiftDivDat1<-agesAllClades_ALL[firstShift,c("CLADE_label","ID_label","nodeMCC","avgFactor","avgIncDec","avgCladeDiv","cladeDiv_low95","cladeDiv_high95","mean","lower","upper")]
#colnames(shiftDivDat1)<-c("clade","ID","nodeMCC","factor","shift","cladeDiv","cladeDiv_low","cladeDiv_high","mean","lower","upper")
#shiftDivDat<-shiftDivDat1[order(shiftDivDat1$factor),]

## First do just the rate-shifts plot
means<-shiftDivDat[,"mean"]
low95s<-shiftDivDat[,"lower"]
up95s<-shiftDivDat[,"upper"]
numTaxa<-shiftDivDat$numTaxa

alphaCol<-1
up<-rgb(1,0,0,alphaCol) #red
down<-rgb(0,0,1,alphaCol) #blue
upDown<-grey(0.5, alpha = alphaCol)
bgcols<-rep(NA,length(shiftDivDat$avgIncDec))

bgcols[which(shiftDivDat$avgIncDec=="upShift")]<-up
bgcols[which(shiftDivDat$avgIncDec=="downShift")]<-down
bgcols[which(shiftDivDat$avgIncDec=="up-or-down")]<-upDown

bgsizes<-shiftDivDat$avgFactor*1.2 # but you need to INVERSE it for the down-shifts...
#val<-bgsizes[which(dat$shift=="downShift")]
#bgsizes[which(dat$shift=="downShift")]<-1/val

cladeNames<-as.vector(shiftDivDat$CLADE_label)
#cladeNames[6]<-"Rhinoloph.-Hipposid."
cladeNames[15]<-"Macropod.-Potoridae" 
nameAtY<-rev(1:18)
#nameAtY<-c(1.5, 3.5, 5, 6.5, 8, 9.5, 11.5, 13.5, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27.5, 29, 30.5)
#nameAtY<-c(1.5, 2, 3.5, 4, 5.5, 6, 7, 8.5, 9, 10, 11.5, 12, 13.5, 14, 15.5, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29.5, 30, 31, 32.5, 33)

pdf(file="divTime_18BAMMnodes_NDexp_plot95pCIs_byDivTime.pdf", onefile=TRUE, width=4,height=5)
#pdf(file="divTime_23BAMMnodes_NDexp_plot95pCIs_byDivTime.pdf", onefile=TRUE, width=4,height=5)
#pdf(file="divTime_23BAMMnodes_NDexp_plot95pCIs_byFactor.pdf", onefile=TRUE, width=4,height=5)
op <- par(oma = c(2,1,1,1) + 0.1,		#c(bottom, left, top, right)
          mar = c(0,0,0,0) + 0.1, lwd=2)
plotCI(x=rev(-means),y=1:length(shiftDivDat[,1]), ui=rev(-low95s),li=rev(-up95s), cex=rev(bgsizes), err="x",sfrac=0, col=NA,xlim=c(-70,10),bty="n",xaxt="n",yaxt="n",xlab="Time before present (Ma)", ylab="",lwd=2, pch=NA)
for(i in nameAtY){
    segments(x0=-70,y0=i,x1=0,y1=i,lty=1, lwd=0.5, col=grey(0.9, alpha=1))
}
	segments(x0=-66.043,y0=0,x1=-66.043,y1=32,lty=2, lwd=2, col=grey(0.5, alpha=1))

plotCI(add=TRUE,x=rev(-means),y=1:length(shiftDivDat[,1]), ui=rev(-low95s),li=rev(-up95s), cex=rev(bgsizes), err="x",sfrac=0, col=rev(bgcols),xlim=c(-110,80),bty="n",xaxt="n",yaxt="n",xlab="Time before present (Ma)", ylab="",lwd=2, pch=20)
axis(side=1,at=seq(-70,0,by=10),labels=c(NA,-60,NA,-40,NA,-20,NA,0), cex.axis=0.85)
text(x=rev(-means)+1,y=(1:length(shiftDivDat[,1]))+0.2,labels=rev(shiftDivDat$ID), cex=0.7, font=2, adj=c(0,0))
text(x=-62,y=(1:length(shiftDivDat[,1]))+0.2,labels=rev(numTaxa), cex=0.7, font=2, adj=c(0,0))
#text(x=rep(0,length(agesAllClades[,1])),y=(1:length(agesAllClades[,1]))+0.2,labels=rev(agesAllClades$ID_label), cex=0.7, font=2, adj=c(0,0))
par(xpd=TRUE)
text(x=rep(0,length(cladeNames)),y=nameAtY,labels=cladeNames,cex=0.6, font=1, adj=0)#, adj=1)

#text(x=rep(0,length(agesAllClades[,1])),y=33.7-(nameAtY),labels=cladeNames,cex=0.6, font=1, adj=c(0,0))
dev.off()


#########
## CO-plot -- 24 rate shifts in PHYLO ORDER
cladeNames<-as.vector(shiftDivDat$CLADE_label)
cladeNames[2]<-"Macropod.-Potoridae" 
nameAtY<-rev(1:24)

# Load in CoMET results...
library(TESS); library(plotrix)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMet_inTESS/")
load(file="MAMPHY_10tree_CoMET_summary_5Ma.Rda")
source("tess.plot.output2.R")
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")

##
xMin<- -155
xMax<- 30
labelCex<-0.7

pdf(file="divTime_24BAMMnodes_NDexp_plot95pCIs_byDivTime_Phylo_wExtinctions.pdf", onefile=TRUE, width=6,height=8)
#pdf(file="divTime_24BAMMnodes_NDexp_plot95pCIs_byDivTime_Phylo.pdf", onefile=TRUE, width=6,height=5)
op <- par(oma = c(1,1,1,1) + 0.1,		#c(bottom, left, top, right)
          mar = c(4,4,1,1) + 0.1, lwd=2)
layout(matrix(1:2,2,1), heights=c(0.7,0.3), widths=c(1,1))

# BAMM shifts
plotCI(x=rev(-means),y=1:length(shiftDivDat[,1]), ui=rev(-low95s),li=rev(-up95s), cex=rev(bgsizes), err="x",sfrac=0, col=NA,xlim=c(xMin,xMax),bty="n",xaxt="n",yaxt="n",xlab="Time before present (Ma)", ylab="",lwd=2, pch=NA)
for(i in nameAtY){
    segments(x0=xMin,y0=i,x1=0,y1=i,lty=1, lwd=0.5, col=grey(0.9, alpha=1))
}
	segments(x0=-66.043,y0=0,x1=-66.043,y1=32,lty=2, lwd=2, col=grey(0.5, alpha=1))

plotCI(add=TRUE,x=rev(-means),y=1:length(shiftDivDat[,1]), ui=rev(-low95s),li=rev(-up95s), cex=rev(bgsizes), err="x",sfrac=0, col=rev(bgcols),xlim=c(-110,80),bty="n",xaxt="n",yaxt="n",xlab="Time before present (Ma)", ylab="",lwd=2, pch=20)
axis(side=1,at=seq(-150,0,by=50),labels=c(150,100,50,0), cex.axis=0.85)
text(x=rev(-means)+1,y=(1:length(shiftDivDat[,1]))+0.2,labels=rev(shiftDivDat$ID), cex=0.7, font=2, adj=c(0,0))
text(x=0,y=(1:length(shiftDivDat[,1]))-0.4,labels=rev(numTaxa), cex=labelCex, font=2, adj=c(0,0))
#text(x=rep(0,length(agesAllClades[,1])),y=(1:length(agesAllClades[,1]))+0.2,labels=rev(agesAllClades$ID_label), cex=0.7, font=2, adj=c(0,0))
#par(xpd=TRUE)
text(x=rep(0,length(cladeNames)),y=nameAtY+0.30,labels=cladeNames,cex=labelCex, font=1, adj=0)#, adj=1)

# CoMET shifts
output<-MAMPHY_10tree_CoMET_5Ma
	
	treeAge <- output[["tree"]]$root.age
	treeAge
    numIntervals <- length(output$intervals) - 1
    numIntervals
    plotAt <- 0:numIntervals
    plotAt
    intervalTimes<-seq(0,205,by=5)
    intervalTimes
    intervalSize <- max(intervalTimes)/numIntervals
    intervalSize
 
    type<-"mass extinction times"
    thisOutput <- output[[type]]
    colnames(thisOutput)<-rev(intervalTimes[1:41])

    meanThisOutput <- colMeans(thisOutput)
    toPlot<-round(meanThisOutput,3)[11:41] # plots from 150-0 Mya

    /length(thisOutput[,1]),3)
   	
   	labels <- pretty(c(0, names(toPlot[1])))
    labels
    labelsAt <- c(41,31,21,11) -10
    labelsAt

    criticalPP<- output[[grep(strsplit(type, " ")[[1]][1], grep("CriticalPosteriorProbabilities", names(output), value = TRUE), value = TRUE)]]

barplot(toPlot, space = 0, xaxt = "n", col = "#4DAF4A", 
        border = "#4DAF4A", main = "", ylab = "posterior probability", 
        xlab = "million years ago", yaxt = "s", pch = 19, ylim=c(0,0.4), xlim=c(0, (185/5)))
axis(1, at = labelsAt, labels = labels)
axis(4, at = c(0,criticalPP[1:2]), labels = c(0,(2 * log(output$criticalBayesFactors[1:2]))), las = 1, tick = TRUE, line = -4)
lines(x=c(max(labelsAt),min(labelsAt)),y=c(criticalPP[1],criticalPP[1]),lty = 2,col="grey")
lines(x=c(max(labelsAt),min(labelsAt)),y=c(criticalPP[2],criticalPP[2]),lty = 2,col="grey")
mtext("Bayes factors",side=4, line=-2)
abline(v = (155/5) - (66.02/5), lty = 2, lwd=3, col="dark grey")

#tess.plot.output2(MAMPHY_10tree_CoMET_5Ma, fig.types = c("mass extinction times"), las=2, plot.tree=FALSE,xlim=c(-xMax,-xMin))

dev.off()


## Calc the MEAN and 95% CI for the inferred extinction divTimes...
#####

    type<-"mass extinction times"
    thisOutput <- output[[type]]
    colnames(thisOutput)<-rev(intervalTimes[1:41])

    meanThisOutput <- colMeans(thisOutput)

    sumThisOutput <- colSums(thisOutput)[11:41] # plots from 150-0 Mya
   
    allVals<-vector("list",length(sumThisOutput))
    for(j in 1:length(sumThisOutput)){
    	val<-as.numeric(names(sumThisOutput[j]))
    	allVals[[j]]<-rep(val,sumThisOutput[[j]])
    }
    freqDist<-do.call(c,allVals)

    mean(freqDist)
    #[1] 99.83792
	quantile(freqDist, c(0.025,0.5,0.975))
	# 2.5%  50%  97.5% 
   	# 80    95   135  

    toPlot<-round(meanThisOutput,3)[11:41] # plots from 150-0 Mya

######


type<-"speciation shift times"
    thisOutput <- output[[type]]
    colnames(thisOutput)<-rev(intervalTimes[1:41])
	meanThisOutput <- colSums(thisOutput)
    toPlot<-round(meanThisOutput/length(thisOutput[,1]),3)[11:41] # plots from 150-0 Mya
   	
barplot(add=TRUE, toPlot, space = 0, xaxt = "n", col = "#4DAF4A", 
        border = "#4DAF4A", main = type, ylab = "posterior probability", 
        xlab = "million years ago", yaxt = "s", pch = 19, ylim=c(0,0.5), xlim=c(0, (180/5)))

            





#########
## CO-plot -- 4 vars TOGETHER as SCALED-- 18 rate shifts in TIME ORDER

library(viridis)
## 
numVars<-length(traitVars)
pchVar<-rep(23,4)
colVar<-rep("black",4)
bgFills<-magma(4)
labelCex<-1
xMin<- -3.5 #min(do.call(rbind,diffInMedians_ALL)$lowDat,na.rm=TRUE)
xMax<-  3.5 #max(do.call(rbind,diffInMedians_ALL)$lowDat,na.rm=TRUE)

pdf(file="divTime_18BAMMnodes_NDexp_byDivTime_co-plot_diffInMedians_w4vars_scaled.pdf", onefile=TRUE, width=7,height=6) #5.333333,height=5) #

#quartz(width=7,height=6) #5.333333,height=5) #
op <- par(oma = c(4,0,2,0) + 0.1,		#c(bottom, left, top, right)
          mar = c(0,0,0,0) + 0.1)
#layout(matrix(1:2,1,2), heights=c(1,1), widths=c(2/3,1/3))
layout(matrix(1:3,1,3), heights=rep(1,3), widths=c(0.7,0.2,0.2))


plotCI(x=rev(-means),y=1:length(shiftDivDat[,1]), ui=rev(-low95s),li=rev(-up95s), cex=rev(bgsizes), err="x",sfrac=0, col=NA,xlim=c(-60,30),bty="n",xaxt="n",yaxt="n",xlab="Time before present (Ma)", ylab="",lwd=2, pch=NA)
for(i in nameAtY){
    segments(x0=-60,y0=i,x1=0,y1=i,lty=1, lwd=0.5, col=grey(0.9, alpha=1))
}
#	segments(x0=-66.043,y0=0,x1=-66.043,y1=32,lty=2, lwd=2, col=grey(0.5, alpha=1))

plotCI(add=TRUE,x=rev(-means),y=1:length(shiftDivDat[,1]), ui=rev(-low95s),li=rev(-up95s), cex=rev(bgsizes), err="x",sfrac=0, col=rev(bgcols),xlim=c(-110,80),bty="n",xaxt="n",yaxt="n",xlab="Time before present (Ma)", ylab="",lwd=2, pch=20)
axis(side=1,at=seq(-60,0,by=10),labels=c(-60,NA,-40,NA,-20,NA,0), cex.axis=0.85)
text(x=rev(-means)+1,y=(1:length(shiftDivDat[,1]))+0.2,labels=rev(shiftDivDat$ID), cex=0.7, font=2, adj=c(0,0))
text(x=0,y=(1:length(shiftDivDat[,1]))-0.35,labels=rev(numTaxa), cex=labelCex, font=2, adj=c(0,0))
#text(x=rep(0,length(agesAllClades[,1])),y=(1:length(agesAllClades[,1]))+0.2,labels=rev(agesAllClades$ID_label), cex=0.7, font=2, adj=c(0,0))
par(xpd=TRUE)
text(x=rep(0,length(cladeNames)),y=nameAtY+0.25,labels=cladeNames,cex=labelCex, font=1, adj=0)#, adj=1)

mtext(side=3,text="(a) Rate-shift timing and magnitude", adj=0,cex=labelCex, line=0, font=2)
#mtext(side=2,text="Rate-shifted clades", adj=0.4,cex=1, line=0, font=2)
mtext(side=1,text="Time (millions of years)", adj=0.25,cex=labelCex, font=2, line=3)

	dat<-diffInMedians_ALL[[1]]
	pchs<-pchVar[1]

	#bgFills<-rep(grey(0.5),length(dat$overlapZero))
	#bgFills[which(dat$overlapZero=="bigger")]<-grey(0,alphaCol)
	#bgFills[which(dat$overlapZero=="smaller")]<-grey(1,alphaCol)
	#lineCols<-rep(grey(0.5),length(dat$overlapZero))
	#lineCols[which(dat$overlapZero=="bigger")]<-grey(0,alphaCol)
	#lineCols[which(dat$overlapZero=="smaller")]<-grey(0,alphaCol)

	plotCI(x=dat$meanDat,y=nameAtY,ui=dat$upDat,li=dat$lowDat, cex=0.5, err="x",sfrac=0, col=NA,ylim=c(1,length(nameAtY)),xlim=c(xMin,xMax),bty="o",xlab="", xaxt="n", yaxt="n",ylab="",lwd=2.5, pch=NA) #
	axis(side=1,at=seq(-3,3,by=1),label=c(NA,-2,NA,0,NA,2,NA), cex.axis=0.85)
	for(j in nameAtY){
	    segments(x0=xMin,y0=j,x1=xMax,y1=j,lty=1, lwd=0.5, col=grey(0.9, alpha=alphaCol))
	}
		segments(x0=0,y0=0,x1=0,y1=length(nameAtY)+1,lty=2, lwd=2, col=grey(0.5, alpha=alphaCol))

	plotCI(add=TRUE,x=dat$meanDat,y=nameAtY+0.4,ui=dat$upDat,li=dat$lowDat, cex=1, err="x",sfrac=0, col=colVar[i],pt.bg=bgFills[i],bty="n",lwd=1.5, pch=pchs)#21)#15) #xaxt="n",
	mtext(side=3,text="(b) Clade traits", adj=0,cex=labelCex, line=0, font=2)
	mtext(side=1,text="Trait departure from non-clade", adj=0.8,cex=labelCex, font=2, line=3)

for(i in 2:length(traitVars)){
	dat<-diffInMedians_ALL[[i]]
	pchs<-pchVar[i]
	
	#bgFills<-rep(grey(0.5),length(dat$overlapZero))
	#bgFills[which(dat$overlapZero=="bigger")]<-grey(0,alphaCol)
	#bgFills[which(dat$overlapZero=="smaller")]<-grey(1,alphaCol)
	#lineCols<-rep(grey(0.5),length(dat$overlapZero))
	#lineCols[which(dat$overlapZero=="bigger")]<-grey(0,alphaCol)
	#lineCols[which(dat$overlapZero=="smaller")]<-grey(0,alphaCol)

	plotCI(add=TRUE,x=dat$meanDat,y=(nameAtY+0.4)-(0.15*i),ui=dat$upDat,li=dat$lowDat, cex=1, err="x",sfrac=0, col=colVar[i],pt.bg=bgFills[i],bty="n",lwd=1, pch=pchs)#21)#15) #xaxt="n",
}
plotCI(x=dat$meanDat,y=nameAtY,ui=dat$upDat,li=dat$lowDat, cex=0.5, err="x",sfrac=0, col=NA,ylim=c(1,length(nameAtY)),xlim=c(xMin,xMax),bty="n",xlab="", xaxt="n", yaxt="n",ylab="",lwd=2.5, pch=NA) #
legend("left", legend=traitVars, fill=bgFills, border=colVar, cex=0.8)

dev.off()





#########
## CO-plot -- 6 vars, diff in medians, separate
numVars<-length(traitVars)
pdf(file="divTime_18BAMMnodes_NDexp_plot95pCIs_byDivTime_co-plot_diffInMedians_w6vars.pdf", onefile=TRUE, width=4+(1.333333*numVars),height=5) #5.333333,height=5) #
#pdf(file="divTime_18BAMMnodes_NDexp_plot95pCIs_byDivTime_co-plot_w9vars_wSEs.pdf", onefile=TRUE, width=4+(1.333333*numVars),height=5) #5.333333,height=5) #
#pdf(file="divTime_23BAMMnodes_NDexp_plot95pCIs_byDivTime_co-plot_w9vars_refine.pdf", onefile=TRUE, width=4+(1.333333*numVars),height=5) #5.333333,height=5) #
#pdf(file="divTime_23BAMMnodes_NDexp_plot95pCIs_byDivTime_co-plot_w2vars.pdf", onefile=TRUE, width=6.6666666,height=5) #5.333333,height=5) #
#pdf(file="divTime_24BAMMnodes_NDexp_MCC_plot95pCIs_noPlacentals_COplot_wRangesMasses.pdf", onefile=TRUE, width=6.6666666,height=5) #5.333333,height=5) #

quartz(width=4+(1.333333*numVars),height=5) #5.333333,height=5) #
op <- par(oma = c(4,1,2,1) + 0.1,		#c(bottom, left, top, right)
          mar = c(0,1,0,1) + 0.1)
#layout(matrix(1:2,1,2), heights=c(1,1), widths=c(2/3,1/3))
layout(matrix(1:(1+numVars),1,(1+numVars)), heights=rep(1,(1+numVars)), widths=c(0.7,rep(0.2,numVars)))

plotCI(x=rev(-means),y=1:length(dat[,1]), ui=rev(-low95s),li=rev(-up95s), cex=rev(bgsizes), err="x",sfrac=0, col=NA,xlim=c(-60,20),bty="n",xaxt="n",yaxt="n",xlab="Time before present (Ma)", ylab="",lwd=2, pch=NA)
for(i in nameAtY){
    segments(x0=-60,y0=i,x1=0,y1=i,lty=1, lwd=0.5, col=grey(0.9, alpha=1))
}
#	segments(x0=-66.043,y0=0,x1=-66.043,y1=32,lty=2, lwd=2, col=grey(0.5, alpha=1))

plotCI(add=TRUE,x=rev(-means),y=1:length(dat[,1]), ui=rev(-low95s),li=rev(-up95s), cex=rev(bgsizes), err="x",sfrac=0, col=rev(bgcols),xlim=c(-110,80),bty="n",xaxt="n",yaxt="n",xlab="Time before present (Ma)", ylab="",lwd=2, pch=20)
axis(side=1,at=seq(-60,0,by=10),labels=c(-60,NA,-40,NA,-20,NA,0), cex.axis=0.85)
text(x=rev(-means)+1,y=(1:length(dat[,1]))+0.2,labels=rev(dat$ID), cex=0.7, font=2, adj=c(0,0))
text(x=0,y=(1:length(dat[,1]))-0.35,labels=rev(numTaxa), cex=1, font=2, adj=c(0,0))
#text(x=rep(0,length(agesAllClades[,1])),y=(1:length(agesAllClades[,1]))+0.2,labels=rev(agesAllClades$ID_label), cex=0.7, font=2, adj=c(0,0))
par(xpd=TRUE)
text(x=rep(0,length(cladeNames)),y=nameAtY+0.25,labels=cladeNames,cex=1, font=1, adj=0)#, adj=1)

mtext(side=3,text="(a) Rate-shift timing and magnitude", adj=0,cex=1, line=0, font=2)
#mtext(side=2,text="Rate-shifted clades", adj=0.4,cex=1, line=0, font=2)
mtext(side=1,text="Time (millions of years)", adj=0.4,cex=1, font=2, line=3)

for(i in 1:length(traitVars)){
	dat<-diffInMedians_ALL[[i]]
	
	bgFills<-rep(grey(0.5),length(dat$overlapZero))
	bgFills[which(dat$overlapZero=="bigger")]<-grey(0,alphaCol)
	bgFills[which(dat$overlapZero=="smaller")]<-grey(1,alphaCol)
	lineCols<-rep(grey(0.5),length(dat$overlapZero))
	lineCols[which(dat$overlapZero=="bigger")]<-grey(0,alphaCol)
	lineCols[which(dat$overlapZero=="smaller")]<-grey(0,alphaCol)

	plotCI(x=dat$meanDat,y=nameAtY,ui=dat$upDat,li=dat$lowDat, cex=0.5, err="x",sfrac=0, col=NA,ylim=c(1,length(nameAtY)),xlim=c(min(dat$lowDat,na.rm=TRUE),max(dat$upDat,na.rm=TRUE)),bty="o",xlab="", xaxt="n", yaxt="n",ylab="",lwd=2.5, pch=NA) #
	axis(side=1,at=c(min(round(dat$upDat,0),na.rm=TRUE),0,round((max(round(dat$lowDat,0),na.rm=TRUE)/2),0),max(round(dat$lowDat,0),na.rm=TRUE)), cex.axis=0.85)
	for(j in nameAtY){
	    segments(x0=min(dat$lowDat,na.rm=TRUE),y0=j,x1=max(dat$upDat,na.rm=TRUE),y1=j,lty=1, lwd=0.5, col=grey(0.9, alpha=alphaCol))
	}
		segments(x0=0,y0=0,x1=0,y1=length(nameAtY)+1,lty=2, lwd=2, col=grey(0.5, alpha=alphaCol))

	plotCI(add=TRUE,x=dat$meanDat,y=nameAtY,ui=dat$upDat,li=dat$lowDat, cex=1, err="x",sfrac=0, col=lineCols,pt.bg=bgFills,bty="n",lwd=1.5, pch=21)#15) #xaxt="n",
	mtext(side=2,text=traitVars[i], cex=0.85, line=0, font=2)
	if(i==1){
		mtext(side=3,text="(b) Clade traits", adj=0,cex=1, line=0, font=2)
	} 
}
title(main="", xlab = "Difference in medians (clade - rest of mammals)",
      ylab = "",
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2, adj=0.8)

dev.off()


##########
#############
# Now same 6 vars but back at the SPECIES LEVEL in explaining tip DR
#######
library(ape); library(phytools)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"

tipDataAll<-read.table("MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)
head(tipDataAll)

mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")

## ALL TRAITS to look at
traitVars1<-c("tipDR_mean10k","BM_final_kg","homeRange_km2_ext","geoRange_km2","DispDistAll_km_ext","GenerationLength_d","trophic123","activity123","lifemode1234")
traitDat1<-tipDataAll[,c("tiplabel",traitVars1)]
traitDat<-cbind.data.frame(traitDat1$tiplabel,log(traitDat1$tipDR_mean10k),log(traitDat1$BM_final_kg),log(traitDat1$homeRange_km2_ext),log(traitDat1$geoRange_km2),log(traitDat1$DispDistAll_km_ext),(traitDat1$lifemode1234*log(traitDat1$DispDistAll_km_ext)),log(traitDat1$GenerationLength_d), traitDat1$trophic123,traitDat1$activity123,traitDat1$lifemode1234)
traitVars<-c("Tip DR (log species/Ma)","Body mass (log kg)","Home range (log km2)","Geographic range (log km2)","Dispersal distance (log km)","Dispersal index (log km x lifemode)","Generation Length (log d)", "trophic123","activity123","lifemode1234")
colnames(traitDat)<-c("tiplabel",traitVars)


>> SEE ABOVE UNDER "PGLS start here"


















###
# PLOT means with CIs to visualize
##
# edited this table some...
library(ape); library(phytools); library(plotrix); library(geiger); library(nlme)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"

allDat<-read.table(file=paste("divTime_24BAMMnodes_",bbone,"_MCC_target_cladeDR-mass-range-carn.txt", sep=""), header=TRUE, stringsAsFactors=FALSE)

#pdf(file=paste("plotCI_cladeBAMMdiv-vs-cladeBM_withBars.pdf",sep=""),onefile=TRUE, width=6,height=6)
#pdf(file=paste("plotCI_cladeBAMMdiv-vs-cladeGEO_withBars.pdf",sep=""),onefile=TRUE, width=6,height=6)
#pdf(file=paste("plotCI_cladeBAMMdiv-vs-cladePerCarn_withBars.pdf",sep=""),onefile=TRUE, width=6,height=6)
pdf(file=paste("plotCI_cladeBAMMdiv-vs-cladeAge_withBars.pdf",sep=""),onefile=TRUE, width=6,height=6)
pchs<-rep(1,length(allDat[,1]))
pchs[which(allDat$Indep==1)]<-16

#yLab<-"Clade BAMM divRate"
#xLab<-"Clade log(Body mass in g)"
#xLab<-"Clade log(Geo range in km2)"
#xLab<-"Clade percent carnivory"
xLab<-"Clade crown age (Ma)"

yVar<-allDat$cladeDiv
yLow<-allDat$cladeDiv_low
yHigh<-allDat$cladeDiv_high
#xVar<-log(allDat$bodyMass_mean)
#xLow<-log(allDat$bodyMass_low)
#xHigh<-log(allDat$bodyMass_high)
#xVar<-log(allDat$geoRange_mean)
#xLow<-log(allDat$geoRange_low)
#xHigh<-log(allDat$geoRange_high)
#xVar<-(allDat$perCarn_mean)
#xLow<-(allDat$perCarn_low)
#xHigh<-(allDat$perCarn_high)
xVar<-(allDat$mean)
xLow<-(allDat$lower)
xHigh<-(allDat$upper)

pchs<-rep(21,length(allDat[,1]))
pchs[which(allDat$Indep==1)]<-16

plotCI(x=xVar, y=yVar,ui=yHigh,li=yLow, xlim=c(min(na.omit(xLow)),max(na.omit(xHigh))), ylim=c(min(yLow),max(yHigh)), xlab=xLab, ylab=yLab, sfrac=0,err="y", lwd=2,col="black",scol="black",pch=NA,font.lab=2,cex.axis=0.95,cex.lab=1)

plotCI(add=TRUE, x=xVar, y=yVar,ui=xHigh,li=xLow, sfrac=0,err="x", lwd=2,col="black", pt.bg="white", scol="black",pch=pchs,font.lab=2,cex.axis=0.95,cex.lab=1)

dev.off()


# For the PGLS lines...
mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")

indepDat<-allDat[which(allDat$Indep==1),]

forms<-list(cladeDiv ~ log(bodyMass_mean), cladeDiv ~ log(geoRange_mean), cladeDiv ~ perCarn_mean, cladeDiv ~ mean)
abP_clades<-vector("list",length=length(forms))
for(i in 1:length(forms)){
	if(i==3 | i==4){ var<-as.character(forms[[i]][[3]])
	} else { var<-as.character(forms[[i]][[3]][2]) }
	
	dat<-cbind.data.frame(indepDat[,c("cladeDiv", var)])
	rownames(dat)<-indepDat$repTaxon
	dat2<-na.omit(dat)

	treeDat<-treedata(mamMCC, dat2)
	cladePhy<-treeDat$phy
	plotDat<-as.data.frame(treeDat$dat)

		fitPGLS_lam<-gls(forms[i][[1]], data=plotDat,correlation=corPagel(value=0.5,phy=cladePhy))#, control = list(singular.ok = TRUE))
		sumPGLS_lam<-summary(fitPGLS_lam)
		sumPGLS_lam

		sum<-sumPGLS_lam
		a=round(sum$coef[[1]],3)
		b=round(sum$coef[[2]],3)
		pVal=round(sum$tTable[[8]],3)
		lambda=round(sum$modelStruct[[1]][[1]],3)

		abP_clades[[i]]<-c(a,b,pVal,lambda)
}
abP_clades_all<-do.call(rbind,abP_clades)
colnames(abP_clades_all)<-c("a","b","pVal","lam")

# if GLS:
         a      b  pVal
[1,] 0.253 -0.002 0.755
[2,] 1.162 -0.066 0.008
[3,] 0.267 -0.001 0.235
[4,] 0.402 -0.007 0.000
# PGLS:
         a      b  pVal    lam
[1,] 0.225 -0.002 0.740 -0.179
[2,] 1.146 -0.065 0.008 -0.120
[3,] 0.264 -0.001 0.234 -0.113
[4,] 0.402 -0.007 0.000 -0.113

##
# Plot again but with PGLS lines...

names<-c("cladeBM", "cladeGEO", "cladePerCarn", "cladeAge")

pchs<-rep(21,length(allDat[,1]))
pchs[which(allDat$Indep==1)]<-16

yLab<-"Clade BAMM divRate"
xLabs<-c("Clade log(Body mass in g)", "Clade log(Geo range in km2)", "Clade percent carnivory", "Clade crown age (Ma)")

yVar<-allDat$cladeDiv
yLow<-allDat$cladeDiv_low
yHigh<-allDat$cladeDiv_high

xVars<-cbind.data.frame(log(allDat$bodyMass_mean), log(allDat$geoRange_mean), allDat$perCarn_mean, allDat$mean)
xLows<-cbind.data.frame(log(allDat$bodyMass_low), log(allDat$geoRange_low), allDat$perCarn_low, allDat$lower)
xHighs<-cbind.data.frame(log(allDat$bodyMass_high), log(allDat$geoRange_high), allDat$perCarn_high, allDat$upper)

for(i in 1:length(names)){

pdf(file=paste("plotCI_cladeBAMMdiv-vs-",names[i],"_wBars_wPGLS.pdf",sep=""),onefile=TRUE, width=6,height=6)

plotCI(x=xVars[,i], y=yVar,ui=yHigh,li=yLow, xlim=c(min(na.omit(xLows[,i])),max(na.omit(xHighs[,i]))), ylim=c(min(yLow),max(yHigh)), xlab=xLabs[i], ylab=yLab, sfrac=0,err="y", lwd=2,col="black",scol="black",pch=NA,font.lab=2,cex.axis=0.95,cex.lab=1)

plotCI(add=TRUE, x=xVars[,i], y=yVar,ui=xHighs[,i],li=xLows[,i], sfrac=0,err="x", lwd=2,col="black", pt.bg="white", scol="black",pch=pchs,font.lab=2,cex.axis=0.95,cex.lab=1)

	a=abP_clades_all[i,"a"]
	b=abP_clades_all[i,"b"]
	pVal=abP_clades_all[i,"pVal"]
	lam=abP_clades_all[i,"lam"]

if(pVal < 0.05){ 
	X=seq(min(na.omit(xLows[,i])),max(na.omit(xHighs[,i])),length.out=20)
	lines(x=X,y=b*X+a,col="blue", lwd=2,lty=2)
	} 

mtext(text=paste("PGLS: y = ", a," + ", b, "x; P = ", pVal,"; lam = ",lam,sep=""), side=3, adj=0, font=1)

dev.off()

}

######
# Repeat the same with tipDR... except that it shant be the response var... 
library(ape); library(phytools)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"

tipDataAll_wComments<-read.table(file="MamPhy_5911sp_tiplabel_DR-range-mass-herb-lifemode.txt", header=TRUE)
tipDataAll<-tipDataAll_wComments[,5:9]
rownames(tipDataAll)<-tipDataAll_wComments$tiplabel

setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_RESULTS_tipDR-vs-bodyMass_PGLS")
for(i in 1:ntrees){
	

	resCats_i<-read.table(file=paste("ExplainingDR_logGeoRangeKM2-vs-tipDR_NDexp_sample100_tree",i,"_PGLS-Pagel_estimDR_allRes.txt",sep=""))

#	tree<-rep(i,length=10)
	tree<-rep(i,length=9)
	resALL[[i]]<-cbind.data.frame(resCats_i,tree)

}
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")

allCats<-c("Terrestrial","Flying","Arboreal","Scansorial","Subterranean"#,"Aquatic"
	,"Herbivore","Omnivore","Carnivore", "ALL_mamPhy")


#pdf(file=paste("plotLines_logGeoRangeKM2-vs-tipDR_PGLS-Pagel.pdf",sep=""),onefile=TRUE, width=6,height=6)
pdf(file=paste("plotLines_tipDR-vs-logBodyMass_PGLS-Brownian.pdf",sep=""),onefile=TRUE, width=6,height=6)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")

#datPlot<-tipDataAll[,c("tipDR", "geoRange")]
datPlot<-na.omit(tipDataAll[,c("tipDR", "bodyMass")])

Cols<-c(rainbow(n=8,alpha=0.2), grey(0.5,alpha=0.5))
ColNames<-c(rainbow(n=8,alpha=1), grey(0.5,alpha=1))

#plot(log(geoRange*(110^2)) ~ tipDR, data=datPlot, xlim=c(0,1.0), xaxt="n", yaxt="n",ylab="" ,xlab="Tip DR", font.lab=2, cex=0)#, pch=pchs)
#axis(side=1,labels=TRUE)
#axis(side=2,at=log(c(10^4,10^5,10^6,10^7,10^8)), labels=FALSE)
#text(xpd=NA, x=rep(-0.08,5), y=log(c(10^4,10^5,10^6,10^7,10^8)), labels=c(1,10,100,expression(10^3),expression(10^4))) 
#mtext(side=2, text=bquote(bold(paste("Geographic range (10"^4," km"^"2",")",sep=""))),cex=1,font=2,line=2.5)

plot(tipDR ~ log(bodyMass), data=datPlot, xlim=c(0, max(log(bodyMass))+5), xaxt="n", yaxt="n",ylab="Tip DR", font.lab=2,xlab="", cex=0)#, pch=pchs)
axis(side=1,at=log(c(1,100,10000,10^6,10^8)), labels=FALSE)
text(xpd=NA, y=rep(-0.08,5), x=log(c(1,100,10000,10^6,10^8)), labels=c(1,100,expression(10^4),expression(10^6),expression(10^8)))
mtext(side=1,text="Body mass (g)",cex=1, font=2,line=2.5)
axis(side=2,labels=TRUE)

for(j in 1:length(allCats)){
Cat<-allCats[j]
lineCol<-Cols[j]

for(i in 1:ntrees){
#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_RESULTS_tipDR-vs-bodyMass_PGLS")
#resCats_i<-read.table(file=paste("ExplainingDR_logGeoRangeKM2-vs-tipDR_NDexp_sample100_tree",i,"_PGLS-Pagel_estimDR_allRes.txt",sep=""))
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_RESULTS_tipDR-vs-bodyMass_PGLS/Brownian_estimDR_tipDR-bodyMass")
resCats_i<-read.table(file=paste("ExplainingDR_tipDR-vs-logBM_NDexp_sample100_tree",i,"_PGLS-Brownian_estimDR_allRes.txt",sep=""))
dat<-resCats_i

a=dat[which(rownames(dat)==Cat),][[1]]
b=dat[which(rownames(dat)==Cat),][[2]]
pVal=dat[which(rownames(dat)==Cat),][[3]]
lam=dat[which(rownames(dat)==Cat),][[4]]

if(pVal < 0.05){ 
	X=seq(min(log(datPlot$bodyMass)),max(log(datPlot$bodyMass)),length.out=20)
	lines(x=X,y=b*X+a,col=lineCol, lwd=1,lty=1)
	} else { 
	X=seq(min(log(datPlot$bodyMass)),max(log(datPlot$bodyMass)),length.out=20)
	lines(x=X,y=b*X+a,col=rgb(1,0,0,alpha=0.5), lwd=1,lty=1)
	}
if(i==100){
	text(x=max(log(datPlot$bodyMass))+0.03,y=b*max(log(datPlot$bodyMass))+a, labels=Cat, col=ColNames[j], cex=0.4, font=2, adj=0)
	}	
}
}

dev.off()

######
# Get a SIMPLE CATEGORICAL showing of the tipDR per Category...

library(ape); library(phytools)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"

tipDataAll_wComments<-read.table(file="MamPhy_5911sp_tiplabel_DR-range-mass-herb-lifemode.txt", header=TRUE)
tipDataAll<-tipDataAll_wComments[,5:9]
rownames(tipDataAll)<-tipDataAll_wComments$tiplabel

# get LOGS of BM and GR-- and multiply geoRange to get km2
tipDataMod<-cbind.data.frame(tipDataAll[,c("tipDR", "LIFEMODE", "percentHerb")],log(tipDataAll$bodyMass),log(tipDataAll$geoRange*(110^2)))
colnames(tipDataMod)<-c("tipDR", "LIFEMODE", "percentHerb","bodyMass", "geoRange")

# subset the data
allCats<-c("Terrestrial","Flying","Arboreal","Scansorial","Subterranean","Aquatic","Herbivore","Omnivore","Carnivore", "HerbLess50", "CarnMore50")
allCatsTraits<-vector("list", length(allCats))

LMcats<-c("Terrestrial","Flying","Arboreal","Scansorial","Subterranean","Aquatic")

for(i in 1:length(LMcats)){
	allCatsTraits[[i]]<-tipDataMod[which(tipDataMod$LIFEMODE==LMcats[i]), c("tipDR", "bodyMass", "geoRange")]
}

trophic_cats<-c("Herbivore","Omnivore","Carnivore", "HerbLess50", "CarnMore50")
	
	allCatsTraits[[7]]<-tipDataMod[which(tipDataMod$percentHerb==100), c("tipDR", "bodyMass", "geoRange")]
	allCatsTraits[[8]]<-tipDataMod[which(tipDataMod$percentHerb < 100 & tipDataMod$percentHerb > 0 ), c("tipDR", "bodyMass", "geoRange")]
	allCatsTraits[[9]]<-tipDataMod[which(tipDataMod$percentHerb==0), c("tipDR", "bodyMass", "geoRange")]
	allCatsTraits[[10]]<-tipDataMod[which(tipDataMod$percentHerb < 50), c("tipDR", "bodyMass", "geoRange")]
	allCatsTraits[[11]]<-tipDataMod[which(tipDataMod$percentHerb >= 50), c("tipDR", "bodyMass", "geoRange")]

# summarize as means and 95% CIs

means<-c(NA); lows<-c(NA); highs<-c(NA)
#var<-"tipDR"
#var<-"bodyMass"
var<-"geoRange"
for(i in 1:length(allCatsTraits)){
	dat<-na.omit(allCatsTraits[[i]])
	
		means[i]<-mean(dat[,var])
		quantDat<-quantile(dat[,var],c(0.035,0.975))
		lows[i]<-quantDat[[1]]
		highs[i]<-quantDat[[2]]
}
means[6]<-NA
lows[6]<-NA
highs[6]<-NA


pdf(file=paste("plotCI_sumAcrossCats_",var,".pdf",sep=""),onefile=TRUE, width=6,height=4)

plotCI(x=1:length(means), y=means,ui=highs,li=lows, ylab=var,xaxt="n", xlab="", sfrac=0,err="y", lwd=2,col="black",scol="black",pch=20,font.lab=2,cex.axis=0.95,cex.lab=1)

axis(1, at=c(1:length(means)), labels = FALSE)
text(x= c(1:length(means)), y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), cex=0.9,font=2,srt = 45, labels = allCats, xpd = TRUE, adj=c(0.95,0.05))

dev.off()





abline(h=0,lty=1, lwd=2, col=grey(0.6, alpha=0.5)) #rgb(1,0,0,alpha=0.5))#
plotCI(add=TRUE,x=1:length(allCats),y=means[reOrd_index],ui=up95[reOrd_index],li=low95[reOrd_index], sfrac=0, err="y", lwd=2,col=cols[reOrd_index],scol="black",pch=20,font.lab=2,cex.axis=0.95,cex.lab=1)

#mtext(text="GLS: tip DR ~ log(Body Mass)", side=3, font=2, adj=0)



	means<-vector("list", length(dat[,1:3]))
	lows<-vector("list", length(dat[,1:3]))
	highs<-vector("list", length(dat[,1:3]))

	for(j in 1:length(means)){
		means[[j]]<-mean(dat[,j])
		quantDat<-quantile(dat[,j],c(0.035,0.975))
		lows[[j]]<-quantDat[[1]]
		highs[[j]]<-quantDat[[2]]
	}



}





  






##
# VISUALIZE those results other than via lines... 95% CIs on the SLOPE (i.e., the EFFECT)
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_RESULTS_tipDR-vs-bodyMass_PGLS")
library(plotrix)


# (non-log / log) tip DR ~ log body mass
#========

#name<-"sum100_ExplainingDR_tipDR-vs-logBM_NDexp_sample100_GLS"
#name<-"sum100_ExplainingDR_logDR-vs-logBM_NDexp_sample100_GLS"
#name<-"sum100_ExplainingDR_tipDR-vs-logBM_NDexp_sample100_PGLS-Brownian"
#name<-"sum100_ExplainingDR_logDR-vs-logBM_NDexp_sample100_PGLS-Brownian"
#name<-"sum100_ExplainingDR_tipDR-vs-logGeoRange_NDexp_sample100_PGLS-Brownian"
#name<-"sum100_ExplainingDR_tipDR-vs-logGeoRangeKM2_NDexp_sample100_GLS"
#name<-"sum100_ExplainingDR_logGeoRangeKM2-vs-tipDR_NDexp_sample100_PGLS-Brownian"
name<-"sum100_ExplainingDR_logGeoRangeKM2-vs-tipDR_NDexp_sample100_PGLS-Pagel"

perCatRes<-read.table(file=paste(name,"_estimDR_allRes.txt",sep=""), row.names=NULL)

allCats<-c("All mammals","Terrestrial","Flying","Arboreal","Scansorial","Subterranean"#,"Aquatic"
	,"Herbivore","Omnivore","Carnivore")

means<-c(); low95<-c(); up95<-c(); pVals<-c(); lams<-c()
for(i in 1:length(allCats)){
	means[i]<-perCatRes[i*3-2,"b"]
	low95[i]<-perCatRes[i*3-1,"b"]
	up95[i]<-perCatRes[i*3,"b"]
	pVals[i]<-perCatRes[i*3-2,"pVal"]
	lams[i]<-perCatRes[i*3-2,"lam"]
}
cols<-rep("black",length(allCats))
cols[which(pVals > 0.05)]<-"red"

#reOrd_index<-c(10,1:9)
reOrd_index<-c(9,1:8)

yMax<- 1.5 #0.05 # 0.4
yMin<- -4.5 #-0.10 # -0.15 # 
yLab<-"Effect on tip Geo Range" #"Effect on tip DR"

pdf(file=paste("plotCI_",name,".pdf",sep=""),onefile=TRUE, width=6,height=4)
#jpeg(file=paste("plotCI_",name,".jpg",sep=""),units="in", res=200,width=6,height=4)

plotCI(x=1:length(allCats), y=means[reOrd_index],ui=up95[reOrd_index],li=low95[reOrd_index], ylim=c(yMin,yMax),xaxt="n", xlab="", ylab=yLab, sfrac=0,err="y", lwd=2,col="white",scol="white",pch=20,font.lab=2,cex.axis=0.95,cex.lab=1)
abline(h=0,lty=1, lwd=2, col=grey(0.6, alpha=0.5)) #rgb(1,0,0,alpha=0.5))#
plotCI(add=TRUE,x=1:length(allCats),y=means[reOrd_index],ui=up95[reOrd_index],li=low95[reOrd_index], sfrac=0, err="y", lwd=2,col=cols[reOrd_index],scol="black",pch=20,font.lab=2,cex.axis=0.95,cex.lab=1)

axis(1, at=c(1:length(allCats)), labels = FALSE)
text(x= c(1:length(allCats)), y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), cex=0.9,font=2,srt = 45, labels = allCats, xpd = TRUE, adj=c(0.95,0.05))
#mtext(text="GLS: tip DR ~ log(Body Mass)", side=3, font=2, adj=0)
#mtext(text="GLS: log(tip DR) ~ log(Body Mass)", side=3, font=2, adj=0)
#mtext(text="PGLS Brownian: tip DR ~ log(Body Mass)", side=3, font=2, adj=0)
#mtext(text="PGLS Brownian: log(tip DR) ~ log(Body Mass)", side=3, font=2, adj=0)
#mtext(text="PGLS Brownian: tip DR ~ log(Geo Range)", side=3, font=2, adj=0)
#mtext(text="GLS: tip DR ~ log(Geo Range km2)", side=3, font=2, adj=0)
#mtext(text="PGLS Brownian: log(tip Geo Range km2) ~ tip DR", side=3, font=2, adj=0)
mtext(text="PGLS Pagel: log(tip Geo Range km2) ~ tip DR", side=3, font=2, adj=0)
text(x= c(1:length(allCats)), y = yMin, cex=0.9,font=2, labels = round(lams[reOrd_index],2), xpd = TRUE)#, adj=c(0.95,0.05))

dev.off()


















######
# ===============================
##
# PLOT the COs for each of those divTimes...
library(plotrix)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"
agesAllClades_ALL<-read.table(file=paste("divTime_24BAMMnodes_",bbone,"_MCC_target_DETAILS.txt", sep=""), header=TRUE)

agesAllClades<-agesAllClades_ALL[-c(5:6),] # removes Placentalia bs is too deep time

means<-agesAllClades[,"mean"]
low95s<-agesAllClades[,"lower"]
up95s<-agesAllClades[,"upper"]

up<-rgb(1,0,0,0.5) #red
down<-rgb(0,0,1,0.5) #blue
upDown<-grey(0.5, alpha = 0.5)
bgcols<-rep(NA,length(agesAllClades$shift))

bgcols[which(agesAllClades$shift=="upShift")]<-up
bgcols[which(agesAllClades$shift=="downShift")]<-down
bgcols[which(agesAllClades$shift=="up-or-down")]<-upDown

bgsizes<-agesAllClades$factor*1.2 # but you need to INVERSE it for the down-shifts...
val<-bgsizes[which(agesAllClades$shift=="downShift")]
bgsizes[which(agesAllClades$shift=="downShift")]<-1/val


cladeNames<-as.vector(na.omit(agesAllClades$CLADE_label))
cladeNames[10]<-"Rhinoloph.-Hipposid."
cladeNames[2]<-"Macropod.-Potoridae" 
nameAtY<-c(1.5, 3.5, 5, 6.5, 8, 9.5, 11.5, 13.5, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27.5, 29, 30.5)
#nameAtY<-c(1.5, 2, 3.5, 4, 5.5, 6, 7, 8.5, 9, 10, 11.5, 12, 13.5, 14, 15.5, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29.5, 30, 31, 32.5, 33)

pdf(file="divTime_24BAMMnodes_NDexp_MCC_plot95pCIs_noPlacentals.pdf", onefile=TRUE, width=4,height=5)
op <- par(oma = c(2,1,1,1) + 0.1,		#c(bottom, left, top, right)
          mar = c(0,0,0,0) + 0.1, lwd=2)
plotCI(x=rev(-means),y=1:length(agesAllClades[,1]), ui=rev(-low95s),li=rev(-up95s), cex=rev(bgsizes), err="x",sfrac=0, col=NA,xlim=c(-70,10),bty="n",xaxt="n",yaxt="n",xlab="Time before present (Ma)", ylab="",lwd=2, pch=NA)
for(i in 32-(nameAtY)){
    segments(x0=-70,y0=i,x1=0,y1=i,lty=1, lwd=0.5, col=grey(0.9, alpha=1))
}
	segments(x0=-66.043,y0=0,x1=-66.043,y1=32,lty=2, lwd=2, col=grey(0.5, alpha=1))

plotCI(add=TRUE,x=rev(-means),y=1:length(agesAllClades[,1]), ui=rev(-low95s),li=rev(-up95s), cex=rev(bgsizes), err="x",sfrac=0, col=rev(bgcols),xlim=c(-110,80),bty="n",xaxt="n",yaxt="n",xlab="Time before present (Ma)", ylab="",lwd=2, pch=20)
axis(side=1,at=seq(-70,0,by=10),labels=c(NA,-60,NA,-40,NA,-20,NA,0), cex.axis=0.85)
text(x=rev(-means)+1,y=(1:length(agesAllClades[,1]))+0.2,labels=rev(agesAllClades$ID_label), cex=0.7, font=2, adj=c(0,0))
#text(x=rep(0,length(agesAllClades[,1])),y=(1:length(agesAllClades[,1]))+0.2,labels=rev(agesAllClades$ID_label), cex=0.7, font=2, adj=c(0,0))
par(xpd=TRUE)
text(x=rep(0,length(cladeNames)),y=32-(nameAtY[1:23]),labels=cladeNames,cex=0.6, font=1, adj=0)#, adj=1)

#text(x=rep(0,length(agesAllClades[,1])),y=33.7-(nameAtY),labels=cladeNames,cex=0.6, font=1, adj=c(0,0))
dev.off()




########
# AS a CO-PLOT with GEO RANGE and TROPHIC LEVEL-- not BODY SIZE...
##
# Get the DIFF in MEDIANS-- from Mann-Whitney U, plus 95% CI
library(plotrix); library(phytools)

agesAllClades_orig<-read.table(file=paste("divTime_24BAMMnodes_",bbone,"_MCC_target.txt", sep=""))
nodeNums<-agesAllClades_orig$nodes[c(1,3,5,7,8,10,11,13,15,17:29,31,32)]
tipDataAll<-read.table(file="MamPhy_5911sp_tiplabel_DR-range-mass.txt", header=TRUE)

#eltonAllUnmatched<-read.table(file="EltonTraits_MamFuncDat.txt", header=TRUE)
eltonAll<-read.table(file="EltonTraits2014_mamPhyMatch5911_4949spp_wData.txt", header=TRUE)

eltonHerb<-eltonAll[,9:12]
percentHerb<-rowSums(eltonHerb, na.rm=FALSE)

names(percentHerb)<-eltonAll[,"tiplabel"]
tipDataAllHerb<-cbind.data.frame(tipDataAll,percentHerb)
write.table(tipDataAllHerb, file="MamPhy_5911sp_tiplabel_DR-range-mass-herb.txt")
tipDataAllHerb<-read.table(file="MamPhy_5911sp_tiplabel_DR-range-mass-herb.txt")


library(plotrix); library(phytools)

agesAllClades_orig<-read.table(file=paste("divTime_24BAMMnodes_",bbone,"_MCC_target.txt", sep=""))
nodeNums<-agesAllClades_orig$nodes[c(1,3,5,7,8,10,11,13,15,17:29,31,32)]
tipDataAll_wComments<-read.table(file="MamPhy_5911sp_tiplabel_DR-range-mass-herb-lifemode.txt", header=TRUE)
tipDataAll<-tipDataAll_wComments[,1:9]



mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")

cladeHERB_all<-vector("list",length(nodeNums))
for(i in 1:length(nodeNums)){
	tipNums<-getDescendants(mamMCC,nodeNums[i])
	cladeTips<-as.vector(na.omit(mamMCC$tip.label[tipNums]))
	cladeHERB_i<-as.numeric(na.omit(tipDataAllHerb[match(cladeTips,rownames(tipDataAllHerb)),"percentHerb"]))
	cladeHERB_all[[i]]<-cladeHERB_i
}

cladeGEO_all<-vector("list",length(nodeNums))
for(i in 1:length(nodeNums)){
	tipNums<-getDescendants(mamMCC,nodeNums[i])
	cladeTips<-as.vector(na.omit(mamMCC$tip.label[tipNums]))
	cladeGEO_i<-as.numeric(na.omit(tipDataAll[match(cladeTips,rownames(tipDataAll)),"geoRange"]))
	cladeGEO_all[[i]]<-cladeGEO_i
}
cladeGEO_all[[6]]<-NA # Because no species in CETACEA have geographic ranges-- AH, because focused on LAND areas...

cladeBM_all<-vector("list",length(nodeNums))
for(i in 1:length(nodeNums)){
	tipNums<-getDescendants(mamMCC,nodeNums[i])
	cladeTips<-as.vector(na.omit(mamMCC$tip.label[tipNums]))
	cladeBM_i<-as.numeric(na.omit(tipDataAll[match(cladeTips,rownames(tipDataAll)),"bodyMass"]))
	cladeBM_all[[i]]<-cladeBM_i
}

datTrait<-cladeHERB_all
allMams<-(as.numeric(na.omit(tipDataAllHerb[,"percentHerb"])))
#datTrait<-cladeGEO_all
#allMams<-log(as.numeric(na.omit(tipDataAll[,"geoRange"])))
#datTrait<-cladeBM_all
#allMams<-log(as.numeric(na.omit(tipDataAll[,"bodyMass"])))

meanDat<-rep(NA,length(nodeNums)); lowDat<-rep(NA,length(nodeNums)); upDat<-rep(NA,length(nodeNums))
for (i in 1:length(nodeNums)){
	#dat<-(as.numeric(na.omit(datTrait[[i]])))
	dat<-100-(as.numeric(na.omit(datTrait[[i]])))
		quantDat<-quantile(dat,c(0.025,0.975))
		meanDat[i]<-mean(dat)
		lowDat[i]<-quantDat[[1]]
		upDat[i]<-quantDat[[2]]
}
#cladeHERB_mean95CI<-cbind.data.frame(meanDat,lowDat,upDat)
cladeCARN_mean95CI<-cbind.data.frame(meanDat,lowDat,upDat)


meanDat<-rep(NA,length(nodeNums)); lowDat<-rep(NA,length(nodeNums)); upDat<-rep(NA,length(nodeNums)); pVals<-rep(NA,length(nodeNums)); overlapZero<-rep(NA,length(nodeNums))
for (i in 1:length(nodeNums)){
	#dat<-log(as.numeric(na.omit(datTrait[[i]])))
	dat<-(as.numeric(na.omit(datTrait[[i]])))
	#if(i==6) { next	} else {
		test<-wilcox.test(x=allMams, y=dat, alternative="two.sided", conf.int=TRUE)
		meanDat[i]<-test$estimate[[1]]
		lowDat[i]<-test$conf.int[1]
		upDat[i]<-test$conf.int[2]
		pVals[i]<-test$p.value
		overlapZero[i]<-if (test$conf.int[1] < 0 && test$conf.int[2] < 0 ) {"bigger"} else if (test$conf.int[1] > 0 && test$conf.int[2] > 0) {"smaller"} else {"NS"}
	#}
}
#cladeBM_diffMedians<-cbind.data.frame(meanDat,lowDat,upDat,pVals,overlapZero)
#cladeGEO_diffMedians<-cbind.data.frame(meanDat,lowDat,upDat,pVals,overlapZero)
cladeHERB_diffMedians<-cbind.data.frame(meanDat,lowDat,upDat,pVals,overlapZero)
	# >> NOT useful for herbivory though !!
	# I just want to look at SPECTRUM of herbivory within the clade-- i.e., 95% CI on the distribution, per clade
	# So not relative to all Mams, just on its own.

save.image(file="divTime_24BAMMnodes_workspace.Rdata")
load(file="divTime_24BAMMnodes_workspace.Rdata")

rownames(cladeGEO_diffMedians)<-agesAllClades_orig$ID[c(1,3,5,7,8,10,11,13,15,17:29,31,32)]

geoSmaller<-cladeGEO_diffMedians[which(cladeGEO_diffMedians$overlapZero=="smaller"),]
	#> mean(geoSmaller[,1])
	#[1] 0.8874973
	#> range(geoSmaller[,1])
	#[1] 0.3067098 2.1282401
# But need to convert back out of log units and from grid cells to km2	
	#> exp(0.887)*(110*110)
	#[1] 29,376.81
	#> exp(0.307)*(110*110)
	#[1] 16448.03
	#> exp(2.128)*(110*110)
	#[1] 101616.5
# vs ALL mammals:
	#allMams<-log(as.numeric(na.omit(tipDataAll[,"geoRange"])))
	#> median(allMams)
	#[1] 3.583519
	#> exp(3.583519)*(110*110)
	#[1] 435,600
	#allMams<-(as.numeric(na.omit(tipDataAll[,"geoRange"])))
	#> median(allMams)
	#[1] 36
	#> mean(allMams)
	#[1] 159.8527
	#> (159.8527)*(110*110)
	#[1] 1,934,218
# vs ACTUAL rate-shifted clades ranges...

smaller<-which(cladeGEO_diffMedians$overlapZero=="smaller")
cladeGEO_smaller<-vector("list",length(smaller))
for (i in smaller){
	cladeGEO_smaller[[i]]<-cladeGEO_all[[i]]
}
allSmaller<-unlist(cladeGEO_smaller)
	#> mean(allSmaller)
	#[1] 86.09116
	#> quantile(allSmaller,c(0.025,0.975))
	#  2.5%  97.5% 
	#  1.00 581.65 
	#> median(allSmaller)
	#[1] 19
	#> 19*(110*110)
	#[1] 229,900


## CO-plot
pdf(file="divTime_24BAMMnodes_NDexp_MCC_plot95pCIs_noPlacentals_COplot_wRangesTrophic.pdf", onefile=TRUE, width=6.6666666,height=5) #5.333333,height=5) #
#pdf(file="divTime_24BAMMnodes_NDexp_MCC_plot95pCIs_noPlacentals_COplot_wRangesMasses.pdf", onefile=TRUE, width=6.6666666,height=5) #5.333333,height=5) #

#quartz(width=6.6666666,height=5) #5.333333,height=5) #
op <- par(oma = c(2,1,1,1) + 0.1,		#c(bottom, left, top, right)
          mar = c(0,0,0,1) + 0.1)
#layout(matrix(1:2,1,2), heights=c(1,1), widths=c(2/3,1/3))
layout(matrix(1:3,1,3), heights=c(1,1), widths=c(0.6,0.2,0.2))

yPlaces<-32-(nameAtY)
fonts<-c(rep(2,8),4,rep(2,11),4,4,2)

plotCI(x=rev(-means),y=1:length(agesAllClades[,1]), ui=rev(-low95s),li=rev(-up95s), cex=rev(bgsizes), err="x",sfrac=0, col=NA,ylim=c(0,31),xlim=c(-70,21),bty="n",xaxt="n",yaxt="n",xlab="Time before present (Ma)", ylab="",lwd=2, pch=NA)
for(i in yPlaces){
    segments(x0=-70,y0=i,x1=0,y1=i,lty=1, lwd=0.5, col=grey(0.9, alpha=1))
}
	segments(x0=-66.043,y0=0,x1=-66.043,y1=32,lty=2, lwd=2, col=grey(0.5, alpha=1))

plotCI(add=TRUE,x=rev(-means),y=1:length(agesAllClades[,1]), ui=rev(-low95s),li=rev(-up95s), cex=rev(bgsizes), err="x",sfrac=0, col=rev(bgcols),ylim=c(0,31),xlim=c(-70,10),bty="n",xaxt="n",yaxt="n",xlab="Time before present (Ma)", ylab="",lwd=3, pch=20)
axis(side=1,at=seq(-70,0,by=10),labels=c(NA,-60,NA,-40,NA,-20,NA,0), cex.axis=1)
text(x=rev(-means)+1,y=(1:length(agesAllClades[,1]))+0.2,labels=rev(agesAllClades$ID_label), cex=0.95, font=2, adj=c(0,0))
#text(x=rep(0,length(agesAllClades[,1])),y=(1:length(agesAllClades[,1]))+0.2,labels=rev(agesAllClades$ID_label), cex=0.7, font=2, adj=c(0,0))
par(xpd=TRUE)
text(x=rep(15,length(cladeNames)),y=32-(nameAtY),labels=cladeNames,cex=0.95, font=fonts, adj=0.5)#, adj=1)
mtext(side=3,text="Rate shift timing and magnitude", adj=0.4,cex=0.85)

dat<-cladeGEO_diffMedians[-c(3:3),]

plotCI(x=-dat$meanDat,y=yPlaces, ui=-dat$lowDat,li=-dat$upDat, cex=0.5, err="x",sfrac=0, col=NA,ylim=c(0,31),xlim=c(-3,3),bty="n",xlab="", xaxt="n",yaxt="n",ylab="",lwd=2, pch=NA)
for(i in 32-(nameAtY)){
    segments(x0=-3,y0=i,x1=20,y1=i,lty=1, lwd=0.5, col=grey(0.9, alpha=1))
}
	segments(x0=0,y0=0,x1=0,y1=32,lty=1, lwd=2, col=rgb(1,0,0,alpha=0.5))

axis(side=1,at=seq(-3,3,by=1), labels=c(NA,-2,NA,0,NA,2,NA), cex.axis=1)
plotCI(add=TRUE,x=-dat$meanDat,y=yPlaces, ui=-dat$lowDat,li=-dat$upDat, cex=1, err="x",sfrac=0, col="black",bty="n",xaxt="n",yaxt="n",xlab="", ylab="",lwd=2.5, pch=15)
mtext(side=3,text="Geographic range", cex=0.85)

#dat<-cladeBM_diffMedians[-c(3:3),]
#dat<-cladeHERB_diffMedians[-c(3:3),]
#dat<-cladeHERB_mean95CI[-c(3:3),]
dat<-cladeCARN_mean95CI[-c(3:3),]

#plotCI(x=-dat$meanDat,y=yPlaces, ui=-dat$lowDat,li=-dat$upDat, cex=0.5, err="x",sfrac=0, col=NA,ylim=c(0,31),xlim=c(-4,10),bty="n",xlab="", xaxt="n",yaxt="n",ylab="",lwd=3, pch=NA)
#mtext(side=3,text="Body mass", cex=0.85)
plotCI(x=dat$meanDat,y=yPlaces, ui=dat$upDat, li=dat$lowDat , cex=0.5, err="x",sfrac=0, col=NA,ylim=c(0,31),xlim=c(0,100),bty="n",xlab="", xaxt="n",yaxt="n",ylab="",lwd=3, pch=NA)
mtext(side=3,text="% Carnivory", cex=0.85)
for(i in 32-(nameAtY)){
    segments(x0=-10,y0=i,x1=100,y1=i,lty=1, lwd=0.5, col=grey(0.9, alpha=1))
}
	segments(x0=mean(100-allMams),y0=0,x1=mean(100-allMams),y1=32,lty=1, lwd=2, col=rgb(1,0,0,alpha=0.5))

axis(side=1,at=seq(0,100,by=25), labels=TRUE, cex.axis=1)
plotCI(add=TRUE,x=dat$meanDat,y=yPlaces, ui=dat$upDat, li=dat$lowDat , cex=1, err="x",sfrac=0, col="black",bty="n",xaxt="n",yaxt="n",xlab="", ylab="",lwd=2.5, pch=15)

dev.off()


##
# BIVARIATE of FAMILY-LEVEL spRichness ~ geoRange // and tipDR ~ geoRange

library(plotrix); library(phytools)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"

tipDataAll<-read.table(file="MamPhy_5911sp_tiplabel_DR-range-mass.txt", header=TRUE)

mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")

cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)

famNames<-names(which(table(cladesDR$fam) > 10)) # 74 data points -- could ADJUST this threshold -- 114 > 2

cladeGEO_FAMS<-vector("list",length(famNames))
cladeDIV_FAMS<-vector("list",length(famNames))
cladeTipDR_FAMS<-vector("list",length(famNames))
famNames2<-c()
for(i in 1:length(famNames)){
	famTips<-as.vector(cladesDR[which(cladesDR$fam==famNames[i]),"tiplabel"])
	cladeDIV_FAMS[[i]]<-length(famTips)
	famNames2[i]<-famNames[i]
	cladeGEO_FAMS_i<-as.numeric(na.omit(tipDataAll[match(famTips,rownames(tipDataAll)),"geoRange"]))
	cladeGEO_FAMS[[i]]<-cladeGEO_FAMS_i
	cladeTipDR_FAMS_i<-as.numeric(na.omit(tipDataAll[match(famTips,rownames(tipDataAll)),"tipDR"]))
	cladeTipDR_FAMS[[i]]<-cladeTipDR_FAMS_i
}


cond <- sapply(cladeGEO_FAMS, function(x) length(x) > 4)
cladeGEO_FAMS_10sp5pts<-cladeGEO_FAMS[cond] ## now have 70 data points
cladeDIV_FAMS_10sp5pts<-cladeDIV_FAMS[cond]
famNamesAll<-famNames2[cond]
cladeTipDR_FAMS_10sp5pts<-cladeTipDR_FAMS[cond]

cladeGEO_FAMdat<-data.frame(matrix(NA, nrow = length(cladeGEO_FAMS_10sp5pts), ncol = 6))
for (i in 1:length(cladeGEO_FAMS_10sp5pts)){
	geoDat<-((110^2)*cladeGEO_FAMS_10sp5pts[[i]])
	cladeGEO_FAMdat[i,1]<-mean(geoDat)
	cladeGEO_FAMdat[i,2]<-quantile(geoDat,0.025)[[1]]
	cladeGEO_FAMdat[i,3]<-quantile(geoDat,0.975)[[1]]
	tipDRDat<-cladeTipDR_FAMS_10sp5pts[[i]]
	cladeGEO_FAMdat[i,4]<-mean(tipDRDat)
	cladeGEO_FAMdat[i,5]<-quantile(tipDRDat,0.025)[[1]]
	cladeGEO_FAMdat[i,6]<-quantile(tipDRDat,0.975)[[1]]
}
colnames(cladeGEO_FAMdat)<-c("geoMean","geoLower","geoUpper","DRmean","DRlower","DRupper")
rownames(cladeGEO_FAMdat)<-famNamesAll

famRichness<-unlist(cladeDIV_FAMS_10sp5pts)

allDat<-cbind.data.frame(cladeGEO_FAMdat,famRichness)
write.table(allDat,file="mamPhyData_geoRange_tipDR_richness_70famsGr10sp5pts.txt")

##
# Now PLOT 
library(geiger); library(nlme); library(phylolm)
#setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")

famDataAll<-read.table(file="mamPhyData_geoRange_tipDR_richness_70famsGr10sp5pts.txt", header=TRUE)
mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")

# Get fam-level tree of those taxa for the PGLS



# family mean tipDR ~ family mean log(geoRange)

form<-DRmean ~ log(geoMean)
plot(form,data=famDataAll)


# family log(sp richness) ~ family mean log(geoRange)

form<-log(famRichness) ~ log(geoMean)
plot(form,data=famDataAll)



dat<-na.omit(famDataAll[,1:2])
#dat<-na.omit(shiftedDat[,1:2])


treeDat<-treedata(mamMCC,dat)
datPlot<-as.data.frame(treeDat$dat)

form<-log(geoRange*(110^2)) ~ tipDR

fitGLS<-gls(form, data=datPlot)
sumGLS<-summary(fitGLS)
sumGLS

fitPGLS<-gls(form, data=datPlot,correlation=corBrownian(phy=treeDat$phy))
sumPGLS<-summary(fitPGLS)
sumPGLS




save.image(file="rangeMass_vs_tipDR_workspace.Rdata")

# load
library(geiger); library(nlme)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
load(file="rangeMass_vs_tipDR_workspace.Rdata")
# Now PLOT-- dual



## Just the GeoRanges...
library(viridis)

lineCol<-"black"
cexEq<-0.75
adjEq<-0.98
cexPts<-0.65
shiftCols<-rainbow(length(cladeTips_all))



save.image(file="rangeMass_vs_tipDR_workspace.Rdata")

#pdf(file="tipDR-vs-Traits_geoRange_ONLY_sepRateShiftSlopes.pdf", onefile=TRUE, width=4.5,height=4) #5.333333,height=5) #
pdf(file="tipDR-vs-Traits_bodyMass_ONLY_sepRateShiftSlopes.pdf", onefile=TRUE, width=4.5,height=4) #5.333333,height=5) #

#quartz(width=4.5,height=4) #5.333333,height=5) #
op <- par(oma = c(2,4,1,1) + 0.1,		#c(bottom, left, top, right)
          mar = c(1,0,0,0) + 0.1, cex=0.85)
#layout(matrix(1:2,1,2), heights=c(1,1), widths=c(2/3,1/3))
#layout(matrix(1:2,2,1), heights=c(0.5,0.5), widths=c(1,1))

pointCol<-rep(grey(0.3, alpha=0.3),length(datPlot2[,1]))
#pointCol[match(cladeTips_gr2x,rownames(datPlot))]<-rgb(1,0,0,alpha=1)
pchs<-rep(1,length(datPlot2[,1]))
#pchs[match(cladeTips_gr2x,rownames(datPlot))]<-16

#plot(log(geoRange*(110^2)) ~ tipDR, data=datPlot, col=pointCol, yaxt="n",ylab="" ,xlab="", cex=cexPts, pch=pchs)
#mtext(side=2, text=bquote(bold(paste("Geographic range (10"^4," km"^"2",")",sep=""))),cex=1,font=2,line=2.5)
#axis(side=2,at=log(c(10^4,10^5,10^6,10^7,10^8)), labels=FALSE)
#text(xpd=NA, x=rep(-0.08,5), y=log(c(10^4,10^5,10^6,10^7,10^8)), labels=c(1,10,100,expression(10^3),expression(10^4))) 
plot(log(bodyMass) ~ tipDR, data=datPlot2, col=pointCol, yaxt="n",ylab="" ,xlab="", cex=cexPts, pch=pchs)
mtext(side=2, text="Body mass (g)",cex=1,font=2,line=2.5)
axis(side=2,at=log(c(1,100,10000,10^6,10^8)), labels=FALSE)
text(xpd=NA, x=rep(-0.08,5), y=log(c(1,100,10000,10^6,10^8)), labels=c(1,100,expression(10^4),expression(10^6),expression(10^8)))
#axis(side=1,labels=FALSE)
mtext(side=1, text="Tip DR (species/Ma)",cex=1,font=2,line=2.5)

sum<-sumPGLS2
a=round(sum$coef[[1]],2)
b=round(sum$coef[[2]],2)
pVal=round(sum$tTable[[8]],3)
X=seq(min((datPlot2$tipDR)),max((datPlot2$tipDR)),length.out=20)
lines(x=X,y=b*X+a,col=lineCol, lwd=4,lty=1) 
#mtext(text=bquote("y = "*.(a)*" + "*.(b)*"x; P < 0.001"), side=3, col="black", cex=cexEq, adj=adjEq, line=-1.5)
mtext(text=bquote("y = "*.(a)*" + "*.(b)*"x; P = "*.(pVal)), side=3, col="black", cex=cexEq, adj=adjEq, line=-1.5)
	Xlab<-max((datPlot2$tipDR))
	Ylab<-b*Xlab+a
	text(x=Xlab, y=Ylab, labels="all",font=2, cex=1.5, col=grey(0))#shiftCols[i]) 


for (i in 1:length(cladeTips_all)){
#	if(i==6 || i==3){ next } else {
	if(i==3){ next } else {

	if(abP_2[[i]][3] > 0.05){ next } else {
	dat<-datPlot_shift_all_2[[i]]
	#points(log(geoRange*(110^2)) ~ tipDR, data=dat, col=shiftCols[i], cex=cexPts, pch=pchs)

	X=seq(min((dat$tipDR)),max((dat$tipDR)),length.out=20)
	lines(x=X,y=abP_2[[i]][2]*X+abP_2[[i]][1],col=shiftCols[i], lwd=2,lty=1)
	Xlab<-max((dat$tipDR))
	Ylab<-abP_2[[i]][2]*Xlab+(abP_2[[i]][1])
	text(x=Xlab, y=Ylab, labels=cladeNames[i],font=2, cex=1.5, col=grey(0))#shiftCols[i]) 
		}
	}
}

dev.off()







## SINGLE plot
pdf(file="geoRanges_24BAMMnodes_vs_allMam_DiffInMedians_95pCI.pdf", onefile=TRUE, width=4,height=5)
#quartz(width=8,height=8)
op <- par(oma = c(2,1,1,1) + 0.1,		#c(bottom, left, top, right)
          mar = c(0,0,0,0) + 0.1, lwd=2)

plotCI(x=rev(-means),y=1:length(cladeGEO_all), ui=rev(-lower),li=rev(-upper), cex=1, err="x",sfrac=0, col=NA,xlim=c(-3,3),bty="n",xlab="", yaxt="n",ylab="",lwd=2, pch=NA)
	segments(x0=0,y0=0,x1=0,y1=24,lty=1, lwd=2, col=rgb(1,0,0,alpha=0.5))
plotCI(add=TRUE,x=rev(-means),y=1:length(cladeGEO_all), ui=rev(-lower),li=rev(-upper), cex=1, err="x",sfrac=0, col="black",xlim=c(-3,3),bty="n",xaxt="n",yaxt="n",xlab="", ylab="",lwd=2, pch=1)

dev.off()


text(x=rev(-means)+1,y=(1:length(agesAllClades[,1]))+0.2,labels=rev(agesAllClades$ID_label), cex=0.7, font=2, adj=c(0,0))





#####
# Get the CLIMATE data too.

O18TempdataAll<-read.table("Zachos2008_Hansen2013_O18Temp_data_forR.txt", header=TRUE)

CO2dataAll<-read.table("BeerlingRoyer2011_CO2_data_forR.txt", header=TRUE)
dataSets<-c("B/Ca", "Boron", "Paleosols", "Phytoplankton", "Stomata")

# Plot the O18 Temp data (Zachos 2008) first:

Time<-(-O18TempdataAll$Time_runMean)
TimeOrig<-(-O18TempdataAll$Time_orig)
TempSurface<-O18TempdataAll$TempC_surface
TempDeep<-O18TempdataAll$TempC_deepOcean
SeaLevel<-O18TempdataAll$SeaLevelM

tempCol<-grey(0.3,alpha=0.5)

par(mfrow=c(2,3))
plot(TempSurface ~ Time, col=tempCol)
plot(TempDeep ~ Time, col=tempCol)
plot(SeaLevel ~ Time, col=tempCol)
plot(delta_18O_orig ~ TimeOrig, data=O18TempdataAll, col=tempCol)
plot(delta_18O_runMean ~ Time, data=O18TempdataAll, col=tempCol)



# And plot the CO2 data also
CO2dataAll[,"Age_mean"]<-(-CO2dataAll$Age_mean)

Boron<-CO2dataAll[which(CO2dataAll$Method=="Boron"),]
Paleosols<-CO2dataAll[which(CO2dataAll$Method=="Paleosols"),]
Phytoplankton<-CO2dataAll[which(CO2dataAll$Method=="Phytoplankton"),]
Stomata<-CO2dataAll[which(CO2dataAll$Method=="Stomata"),]

c("Boron")

dataSets<-c("B/Ca", "Boron", "Paleosols", "Phytoplankton", "Stomata")

library(viridis); library(RColorBrewer)
source("DR_circleTree_functions.R")

Cols<-viridis(length(dataSets))
#Cols<-rich.colors(length(dataSets))

dataToPlot<-CO2dataAll[which(CO2dataAll$Method==dataSets[1]),]
plot(CO2_mean ~ Age_mean, data=dataToPlot, col=Cols[1], pch=1+20, xlim=c(-65,0), ylim=c(0,1600),xaxt="n")
axis(side=1,at=seq(from=-60,to=0,by=10),labels=TRUE)
for (i in 2:length(dataSets)){
	dataToPlot<-CO2dataAll[which(CO2dataAll$Method==dataSets[i]),]
	par(new=TRUE)
	plot(CO2_mean ~ Age_mean, data=dataToPlot, col=Cols[i], pch=i+20, xlim=c(-65,0), ylim=c(0,1600),xaxt="n",yaxt="n", xlab="",ylab="")
}




##
# COMBINE the data...
##
library(viridis)
O18TempdataAll<-read.table("Zachos2008_Hansen2013_O18Temp_data_forR.txt", header=TRUE)
O18TempdataAll[,"Time_runMean"]<-(-O18TempdataAll$Time_runMean)
O18TempdataAll_pre<-O18TempdataAll[which(O18TempdataAll$Time_runMean < -35),]
O18TempdataAll_post<-O18TempdataAll[which(O18TempdataAll$Time_runMean > -35),]

CO2dataAll<-read.table("BeerlingRoyer2011_CO2_data_forR.txt", header=TRUE)
CO2dataAll[,"Age_mean"]<-(-CO2dataAll$Age_mean)
dataSets<-c("B/Ca", "Boron", "Paleosols", "Phytoplankton", "Stomata")

data(gradstein04)
library(phyloch); library(strap); library(phytools)
mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")
root<-max(nodeHeights(mamMCC))
mamMCC$root.time<-root

#########
pdf(file="divTime_24BAMMnodes_NDexp_MCC_plot95pCIs_TriPlot_wTempCO2.pdf", onefile=TRUE, width=4,height=8)

hotClim<-"darkolivegreen2" #"darkorange3" # "mediumorchid3"
coldClim<-"lightblue2" #"mediumorchid4"
#quartz(width=4,height=8)
op <- par(oma = c(6,6,0,0) + 0.1,		#c(bottom, left, top, right)
          mar = c(0,0,1,0) + 0.1, lwd=2)
layout(matrix(1:3,3,1),heights=c(0.5,0.25,0.25), widths=c(1,0.5,0.5))

# part a
plot(mamMCC, plot=FALSE, show.tip.label=FALSE, x.lim=c(root-70,root))
axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 1, ages=FALSE, gridty=3, gridcol="grey50")

par(new=TRUE)
plotCI(x=rev(-means),y=1:length(agesAllClades[,1]), ui=rev(-low95s),li=rev(-up95s), cex=rev(bgsizes), err="x",sfrac=0, col=NA,xlim=c(-70,10),bty="n",xaxt="n",yaxt="n",xlab="Time before present (Ma)", ylab="",lwd=2, pch=NA)
for(i in 32-(nameAtY)){
    segments(x0=-70,y0=i,x1=0,y1=i,lty=1, lwd=0.5, col=grey(0.9, alpha=1))
}
	#segments(x0=-66.043,y0=-2,x1=-66.043,y1=32,lty=2, lwd=2, col=grey(0.5, alpha=1))
	rect(xleft=-55.5, ybottom=-2, xright=-54.5, ytop=50, col=hotClim, border=NA) # PETM: Paleocene-Eocene Thermal Maximum / EECO - early Eocene Climatic Optimum
	rect(xleft=-42, ybottom=-2, xright=-41, ytop=50, col=hotClim, border=NA) # MECO: Mid-Eocene Climatic Optimum
	rect(xleft=-34, ybottom=-2, xright=-33, ytop=50, col=coldClim, border=NA) # O1-1 Glaciation
	rect(xleft=-16, ybottom=-2, xright=-15, ytop=50, col=hotClim, border=NA) # MMCO - Mid-Miocene Climatic Optimum

plotCI(add=TRUE,x=rev(-means),y=1:length(agesAllClades[,1]), ui=rev(-low95s),li=rev(-up95s), cex=rev(bgsizes), err="x",sfrac=0, col=rev(bgcols),xlim=c(-110,80),bty="n",xaxt="n",yaxt="n",xlab="Time before present (Ma)", ylab="",lwd=2, pch=20)
axis(side=1,at=seq(-70,0,by=10),cex.axis=0.85, labels=FALSE)#labels=c(NA,-60,NA,-40,NA,-20,NA,0), )
text(x=rev(-means)+1,y=(1:length(agesAllClades[,1]))+0.2,labels=rev(agesAllClades$ID_label), cex=0.7, font=2, adj=c(0,0))
#text(x=rep(0,length(agesAllClades[,1])),y=(1:length(agesAllClades[,1]))+0.2,labels=rev(agesAllClades$ID_label), cex=0.7, font=2, adj=c(0,0))
par(xpd=TRUE)
text(x=rep(-60,length(cladeNames)),y=32-(nameAtY[1:23]),labels=cladeNames,cex=0.7, font=2, adj=1)#, adj=1)


geoscalePlot(-O18TempdataAll$Time_runMean, O18TempdataAll$TempC_surface, boxes="Age",ts.col=TRUE,type="p",xlim=c(0,70))



# part b
tempCol<-grey(0.5,alpha=0.1)
plot(TempC_surface ~ Time_runMean, data=O18TempdataAll, ylim=c(7,32), xlim=c(-70,10), type="n",bty="n", col=tempCol, xaxt="n", yaxt="n")
	segments(x0=-66.043,y0=0,x1=-66.043,y1=52,lty=2, lwd=2, col=grey(0.5, alpha=1))
	rect(xleft=-55.5, ybottom=0, xright=-54.5, ytop=50, col=hotClim, border=NA) # PETM: Paleocene-Eocene Thermal Maximum / EECO - early Eocene Climatic Optimum
	rect(xleft=-42, ybottom=0, xright=-41, ytop=50, col=hotClim, border=NA) # MECO: Mid-Eocene Climatic Optimum
	rect(xleft=-34, ybottom=0, xright=-33, ytop=50, col=coldClim, border=NA) # O1-1 Glaciation
	rect(xleft=-16, ybottom=0, xright=-15, ytop=50, col=hotClim, border=NA) # MMCO - Mid-Miocene Climatic Optimum
points(TempC_surface ~ Time_runMean, data=O18TempdataAll, ylim=c(7,32), xlim=c(-70,10), bty="n", col=tempCol, xaxt="n", yaxt="n")
points(type="l",col = "blue", lwd = 3, lty = 1,with(O18TempdataAll_pre, loess.smooth(y=TempC_surface, x=Time_runMean, span=0.01)))
points(type="l",col = "light blue", lwd = 3, lty = 1,with(O18TempdataAll_post, loess.smooth(y=TempC_surface, x=Time_runMean, span=0.01)))

axis(side=1,at=seq(-70,0,by=10),cex.axis=0.85, labels=FALSE)#labels=c(NA,-60,NA,-40,NA,-20,NA,0), )
axis(side=2, at=c(5,10,15,20,25,30), labels=c(NA,10,NA,20,NA,30))
mtext(side=2,line=3, cex=0.8, text="Surface temperature (deg C)")#, adj=0)
#with(O18TempdataAll, scatter.smooth(y=TempC_surface, x=Time_runMean, span=0.001, lpars=list(col = "yellow", lwd = 3, lty = 1)))

# part c
Cols<-viridis(length(dataSets))

dataToPlot<-CO2dataAll[which(CO2dataAll$Method==dataSets[1]),]
plot(CO2_mean ~ Age_mean, data=dataToPlot, col=Cols[1], pch=1+20, xlim=c(-70,10), ylim=c(0,1600),xaxt="n", bty="n", type="n")
	segments(x0=-66.043,y0=0,x1=-66.043,y1=2500,lty=2, lwd=2, col=grey(0.5, alpha=1))
	rect(xleft=-55.5, ybottom=0, xright=-54.5, ytop=2500, col=hotClim, border=NA) # PETM: Paleocene-Eocene Thermal Maximum / EECO - early Eocene Climatic Optimum
	rect(xleft=-42, ybottom=0, xright=-41, ytop=2500, col=hotClim, border=NA) # MECO: Mid-Eocene Climatic Optimum
	rect(xleft=-34, ybottom=0, xright=-33, ytop=2500, col=coldClim, border=NA) # O1-1 Glaciation
	rect(xleft=-16, ybottom=0, xright=-15, ytop=2500, col=hotClim, border=NA) # MMCO - Mid-Miocene Climatic Optimum
points(CO2_mean ~ Age_mean, data=dataToPlot, col=Cols[1], pch=1+20, xlim=c(-70,10), ylim=c(0,1600),xaxt="n", bty="n", cex=0.5)
axis(side=1,at=seq(-70,0,by=10),labels=c(NA,-60,NA,-40,NA,-20,NA,0), cex.axis=1)
for (i in 2:length(dataSets)){
	dataToPlot<-CO2dataAll[which(CO2dataAll$Method==dataSets[i]),]
	par(new=TRUE)
	points(CO2_mean ~ Age_mean, data=dataToPlot, col=Cols[i], pch=i+20, xlim=c(-65,0), ylim=c(0,1600),xaxt="n",yaxt="n", xlab="",ylab="", bty="n")
}
points(type="l",col = "black", lwd = 3, lty = 1,with(CO2dataAll, loess.smooth(y=CO2_mean, x=Age_mean, span=0.1)))
legend(x=-10,y=1600,legend=dataSets,pch=21:25,col=Cols, cex=0.7)
mtext(side=2,line=3, cex=0.8, text="CO2 ppm")#, adj=0)
mtext(side=1,line=3, cex=0.8, text="Time before present (Ma)")#, adj=0)

dev.off()






# plot JUST the 2 climate layers...
pdf(file="Climates_dual_wTempCO2.pdf", onefile=TRUE, width=4,height=5)

op <- par(oma = c(2,2,1.5,1) + 0.1,		#c(bottom, left, top, right)
          mar = c(1,0,1,0) + 0.1, lwd=2)
layout(matrix(1:2,2,1))#,heights=c(0.5,0.25,0.25), widths=c(1,0.5,0.5))

# part b
tempCol<-grey(0.5,alpha=0.1)
plot(TempC_surface ~ Time_runMean, data=O18TempdataAll, ylim=c(7,32), col=tempCol, xaxt="n", yaxt="n")
axis(side=1, labels=FALSE)
axis(side=2, at=c(5,10,15,20,25,30), labels=c(NA,10,NA,20,NA,30))
mtext(side=3,text="(b) Surface temperature (deg C)", adj=0)

# part c
Cols<-viridis(length(dataSets))

dataToPlot<-CO2dataAll[which(CO2dataAll$Method==dataSets[1]),]
plot(CO2_mean ~ Age_mean, data=dataToPlot, col=Cols[1], pch=1+20, xlim=c(-65,0), ylim=c(0,1600),xaxt="n")
axis(side=1,at=seq(from=-60,to=0,by=10),labels=TRUE)
for (i in 2:length(dataSets)){
	dataToPlot<-CO2dataAll[which(CO2dataAll$Method==dataSets[i]),]
	par(new=TRUE)
	points(CO2_mean ~ Age_mean, data=dataToPlot, col=Cols[i], pch=i+20, xlim=c(-65,0), ylim=c(0,1600),xaxt="n",yaxt="n", xlab="",ylab="", bty="n")
}
mtext(side=3,text="(c) CO2 ppm", adj=0)
legend(x=-14,y=1600,legend=dataSets,pch=21:25,col=Cols, cex=0.5)

dev.off()







######
# ====================
##
#orderToPlot<-c("DIPROTODONTIA", "DIDELPHIMORPHIA", "DASYUROMORPHIA", "PERAMELEMORPHIA", "PAUCITUBERCULATA", "RODENTIA", "MouseRelated", "SquirrelRelated", "GuineaPigRelated", "PRIMATES", "OldWorldMonkey", "NewWorldMonkey", "Strepsirrhini", "LAGOMORPHA", "SCANDENTIA", "CHIROPTERA", "EmballonuridRelated", "VespertilionidRelated", "Yinpterochiroptera", "PhyllostomidRelated", "EULIPOTYPHLA", "CETARTIODACTYLA", "TerrestrialCetartios", "Cetaceans", "CARNIVORA","CatRelated","DogRelated", "PERISSODACTYLA", "PHOLIDOTA")
orderToPlot<-redoNames

means<-finALL[orderToPlot,"mean_dAICrc"]
low95s<-finALL[orderToPlot,"low95_dAICrc"]
up95s<-finALL[orderToPlot,"up95_dAICrc"]

pdf(file="LTT_deltaAICrc_95dists.pdf", onefile=TRUE, width=4,height=10)
op <- par(oma = c(2,2,2,2) + 0.1,
          mar = c(0,0,0,0) + 0.1, lwd=2)

plotCI(x=rev(means),y=1:length(orderToPlot), ui=rev(up95s),li=rev(low95s), cex=1.3, err="x",sfrac=0, xlim=c(-55,55),bty="n",yaxt="n",xlab="delta AICrc", ylab="",lwd=2, pch=20)
abline(v=0,lty=2, lwd=2, col=rgb(1,0,0, alpha=0.5)) #rgb(1,0,0,alpha=0.5))#
text(x=rep(-10,26),y=1:29,labels=rev(orderToPlot), pos=2,cex=0.7)
for(i in 1:length(orderToPlot)){
    abline(h=i,lty=1, lwd=0.5, col=grey(0.6, alpha=0.5))
}
dev.off()


pdf(file="LTT_pieCharts_ofMLdivModels_plot29clades_NDexp_blue_2cat.pdf", onefile=TRUE, width=10,height=10)
op <- par(mfrow = c(5,6),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1, lwd=2)

for (j in 1:length(orderToPlot)){
    percent<-finALL[orderToPlot[j],"X.timeCon"]
    pie(x=c(percent,1-percent),labels=NA,col=c("black","deepskyblue2"))
    mtext(orderToPlot[j],cex=0.6)
}
dev.off()




#######
# Script to extract...
## >> this is SLOW >> needs to be switched to just TWO taxa per clade (family, order, etc) to do the node-finding

ordNames<-names(which(table(cladesDR$ord) > 2))

library(foreach);library(doSNOW)
cl2 = makeCluster(30, type = 'SOCK', outfile="")
registerDoSNOW(cl2)

#divDataORD<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 4))
#colnames(divDataORD)<-c("node","height","height_95%_HPD_MIN","height_95%_HPD_MAX")
#rownames(divDataORD)<-ordNames
#for (i in 1:length(ordNames)){

ordDIVS= foreach(i=1:length(ordNames), .packages=c('ape', 'phytools'), .combine=rbind, .verbose=TRUE) %dopar% {
divDataORD<-vector()
	tips<-cladesDR[which(cladesDR$ord==ordNames[i]),"tiplabel"]
	node<-findMRCA(mamPhy,tips,type="node")
	#node<-getMRCA(mamPhy,tips)
	divDataORD[1]<-node
	divDataORD[2]<-mamPhyMCC_Table[which(mamPhyMCC_Table[,"node"]==node),"height"][[1]]
	divDataORD[3]<-mamPhyMCC_Table[which(mamPhyMCC_Table[,"node"]==node),"height_95%_HPD_MIN"][[1]]
	divDataORD[4]<-mamPhyMCC_Table[which(mamPhyMCC_Table[,"node"]==node),"height_95%_HPD_MAX"][[1]]
	return(divDataORD)
}
colnames(ordDIVS)<-c("node","height","height_95%_HPD_MIN","height_95%_HPD_MAX")
rownames(ordDIVS)<-ordNames
write.table(ordDIVS,file="divDataORD_MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target.txt")



g_es_dr= foreach(i=1:ntrees, .combine = cbind) %dopar% {   
   ES = ES_v2(readCAIC(trees[[i]]))
   DR = 1/ES
   res1 = cbind.data.frame(DR,ES)
   res2 = as.matrix(res1[match(cm$tip.label, rownames(res1)),]) 
   return(res2)
}
write.table(g_es_dr, file="FS2015_resolved_tipCalcs_1kTrees_DR-ES.txt")





# FAMILIES next 
famNames<-names(which(table(cladesDR$fam) > 2))

famDIVS= foreach(i=1:length(famNames), .packages=c('ape', 'phytools'), .combine=rbind, .verbose=TRUE) %dopar% {
divDataFAM<-vector()
	tips<-cladesDR[which(cladesDR$fam==famNames[i]),"tiplabel"]
	node<-findMRCA(mamPhy,tips,type="node")
	#node<-getMRCA(mamPhy,tips)
	divDataFAM[1]<-node
	divDataFAM[2]<-mamPhyMCC_Table[which(mamPhyMCC_Table[,"node"]==node),"height"][[1]]
	divDataFAM[3]<-mamPhyMCC_Table[which(mamPhyMCC_Table[,"node"]==node),"height_95%_HPD_MIN"][[1]]
	divDataFAM[4]<-mamPhyMCC_Table[which(mamPhyMCC_Table[,"node"]==node),"height_95%_HPD_MAX"][[1]]
	return(divDataFAM)
}
colnames(famDIVS)<-c("node","height","height_95%_HPD_MIN","height_95%_HPD_MAX")
rownames(famDIVS)<-famNames
write.table(famDIVS,file="divDataFAM_MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target.txt")



divDataFAM<-data.frame(matrix(NA, nrow = length(famNames), ncol = 3))
colnames(divDataFAM)<-colnames(ages)
rownames(divDataFAM)<-famNames
for (i in 1:length(famNames)){
	tips<-cladesDR[which(cladesDR$ord==famNames[i]),"tiplabel"]
	node<-findMRCA(mamPhy,tips,type="node")
	divDataFAM[i,]<-ages[node-5912,]
}



# CLADES now
cladeNames<-names(which(table(cladesDR$clade) > 2))

cladeDIVS= foreach(i=1:length(cladeNames), .packages=c('ape', 'phytools'), .combine=rbind, .verbose=TRUE) %dopar% {
divDataCLADE<-vector()
	tips<-cladesDR[which(cladesDR$clade==cladeNames[i]),"tiplabel"]
	node<-findMRCA(mamPhy,tips,type="node")
	#node<-getMRCA(mamPhy,tips)
	divDataCLADE[1]<-node
	divDataCLADE[2]<-mamPhyMCC_Table[which(mamPhyMCC_Table[,"node"]==node),"height"][[1]]
	divDataCLADE[3]<-mamPhyMCC_Table[which(mamPhyMCC_Table[,"node"]==node),"height_95%_HPD_MIN"][[1]]
	divDataCLADE[4]<-mamPhyMCC_Table[which(mamPhyMCC_Table[,"node"]==node),"height_95%_HPD_MAX"][[1]]
	return(divDataCLADE)
}
colnames(cladeDIVS)<-c("node","height","height_95%_HPD_MIN","height_95%_HPD_MAX")
rownames(cladeDIVS)<-cladeNames
write.table(cladeDIVS,file="divDataCLADE_MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target.txt")




divDataCLADE<-data.frame(matrix(NA, nrow = length(cladeNames), ncol = 3))
colnames(divDataFAM)<-colnames(ages)
rownames(divDataFAM)<-cladeNames
for (i in 1:length(cladeNames)){
	tips<-cladesDR[which(cladesDR$clade==cladeNames[i]),"tiplabel"]
	node<-findMRCA(mamPhy,tips,type="node")
	divDataCLADE[i,]<-ages[node-5912,]
}


NDexpNEWICK<-read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target_5911sp_newick.tre",sep=""))
NDexpNEXUS<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_MCC_target_5911sp.tre",sep=""))

	node<-getMRCA(NDexpNEXUS,tips)
	node2<-getMRCA(NDexpNEXUS,tips)


Marsupialia
DIPROTODONTIA
Macropodidae-Potoridae
Placentalia
CARNIVORA
CETARTIODACTYLA
Pecora
Cetacea
CHIROPTERA
Vespertilionid-related
Phyllostomidae
Pteropodidae
Pteropus
Rhinolophidae-Hipposideridae
Rhinolophidae
EULIPOTYPHLA
Crocidurinae
PRIMATES
Strepsirrhini
LAGOMORPHA
RODENTIA
Ctenomyidae
Geomyoidea
Cricetidae-Muridae
Cricetidae
Muridae
Murinae
Murinae: Apomys-Melomys
Murinae: Rattus-Srilankamys
Gerbillinae


# prune mamPhy to an ordPhy
## Now prune big tree to ORDER reps
#ordSp<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_1spMostGenes_perORD.txt")
#ordSpOrd<-cbind(ordSp[1],ordSp[4])
#colnames(ordSpOrd)<-c("tip","ord")
#
#ordPhy<-drop.tip(mamPhy,setdiff(mamPhy$tip.label,ordSpOrd$tip))
#write.tree(ordPhy, paste("MamPhy_BDvr_pcsFIXED_",bbone,"_MCC_target_27ORDERS_fullName.tre",sep=""))
#
#ordPhy$tip.label<-as.vector(ordSpOrd[match(ordPhy$tip.label, as.vector(ordSpOrd$tip)),"ord"])
#write.tree(ordPhy, paste("MamPhy_BDvr_pcsFIXED_",bbone,"_MCC_target_27ORDERS_ordNames.tre",sep=""))

# PATCHES
patchSp<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_1spMostGenes_perPC.txt")
newNames<-c("Nesomyidae", "Muridae", "Cricetidae", "Squirrel-related", "Guinea_pig-related", "Eulipotyphla", "Phyllostomid-related", "Vespertilionid-related", "Emballonurid-related", "Yinpterochiroptera", "Marsupialia", "Cetartiodactyla", "Perissodactyla", "Carnivora", "Monotremata", "Pholidota", "Dermoptera", "Anomaluromorpha", "Calomyscidae", "Platacanthomyidae", "Afrotheria", "Xenarthra", "Scandentia", "Primates", "Lagomorpha", "Castorimorpha", "Dipodidae", "Spalacidae")
patchSpPatch<-cbind(patchSp[1],newNames)
colnames(patchSpPatch)<-c("tip","PC")

patchPhy<-drop.tip(mamPhy,setdiff(mamPhy$tip.label,patchSpPatch$tip))
write.tree(patchPhy, paste("MamPhy_BDvr_pcsFIXED_",bbone,"_MCC_target_28PATCHES_fullNames.tre",sep=""))
patchPhy$tip.label<-as.vector(patchSpPatch[match(patchPhy$tip.label, as.vector(patchSpPatch$tip)),"PC"])
write.tree(patchPhy, paste("MamPhy_BDvr_pcsFIXED_",bbone,"_MCC_target_28PATCHES_patchNames.tre",sep=""))


###############
# load back in ordPhy
ordPhy<-read.tree(paste("MamPhy_BDvr_pcsFIXED_",bbone,"_MCC_target_27ORDERS_ordNames.tre",sep=""))
rootTime<-max(node.age(ordPhy)$ages)

# PLOT in a CIRCLE...
pdf(file=paste("MamPhy_BDvr_pcsFIXED_",bbone,"_MCC_target_27ORDERS_ordNames_CIRCLEplot.pdf",sep=""),width=5,height=5)
plot(ordPhy, type="fan", cex=0.7, font=2, label.offset=0.5,open.angle=20,rotate.tree=33, tip.color="white",show.tip.label=TRUE,no.margin=TRUE)
dd<-0.95
cc<-0.9
draw.circle(0,0, radius=rootTime-0, lty=2,col=grey(dd), border=grey(dd))
draw.circle(0,0, radius=rootTime-66, lty=2,col=grey(cc), border=grey(cc))
draw.circle(0,0, radius=rootTime-66, lty=2,col=NA, border=grey(0.5))
#draw.circle(0,0, radius=rootTime-100, lty=2,col=grey(dd), border=grey(dd))
#draw.circle(0,0, radius=rootTime-150, lty=2,col=grey(cc), border=grey(cc))
par(new=TRUE)
plot(ordPhy, type="fan", cex=0.7, font=2, label.offset=0.5,open.angle=20,rotate.tree=33, no.margin=TRUE)
axisPhylo(cex=0.6)

dev.off()

# PLOT THE PC PHY THE SAME WAY...
patchPhy<-read.tree(paste("MamPhy_BDvr_pcsFIXED_",bbone,"_MCC_target_28PATCHES_patchNames.tre",sep=""))
rootTime<-max(node.age(patchPhy)$ages)

# radial
pdf(file=paste("MamPhy_BDvr_pcsFIXED_",bbone,"_MCC_target_28PATCHES_patchNames_CIRCLEplot.pdf",sep=""),width=5,height=5)
plot(patchPhy, type="fan", cex=0.7, font=2, label.offset=0.5,open.angle=20,rotate.tree=33, tip.color="white",show.tip.label=TRUE,no.margin=TRUE)
dd<-0.95
cc<-0.9
draw.circle(0,0, radius=rootTime-0, lty=2,col=grey(dd), border=grey(dd))
draw.circle(0,0, radius=rootTime-66, lty=2,col=grey(cc), border=grey(cc))
draw.circle(0,0, radius=rootTime-66, lty=2,col=NA, border=grey(0.5))
#draw.circle(0,0, radius=rootTime-100, lty=2,col=grey(dd), border=grey(dd))
#draw.circle(0,0, radius=rootTime-150, lty=2,col=grey(cc), border=grey(cc))
par(new=TRUE)
plot(patchPhy, type="fan", cex=0.7, font=2, label.offset=0.5,open.angle=20,rotate.tree=33, no.margin=TRUE)
axisPhylo(cex=0.6)

dev.off()


# =============================
############
# Plot all 100 patchPhys....
# like a densitree, but manual...
###
#
# same for 100 trees
# load mamPhy 100 tree sample
mamPhy_100<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_nexus.trees",sep=""))
mamPhy100<-lapply(mamPhy_100,ladderize)
class(mamPhy100)<-"multiPhylo"

patchPhy_100<-vector("list",length=100)
for (i in 1:length(mamPhy100)){
	patchPhy_100[[i]]<-drop.tip(mamPhy100[[i]],setdiff(mamPhy100[[i]]$tip.label,patchSpPatch$tip))
	patchPhy_100[[i]]$tip.label<-as.vector(patchSpPatch[match(patchPhy_100[[i]]$tip.label, as.vector(patchSpPatch$tip)),"PC"])
}
class(patchPhy_100)<-"multiPhylo"
write.nexus(patchPhy_100, file=paste("MamPhy_BDvr_pcsFIXED_",bbone,"_sample100_28PATCHES_patchNames.tre",sep=""))

# READ back in the 100 patch trees
patchPhy_100<-read.nexus(file=paste("MamPhy_BDvr_pcsFIXED_",bbone,"_sample100_28PATCHES_patchNames.tre",sep=""))

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
# 59 41 

treeNum<-c(1:length(patchPhy_100))
ID<-cbind.data.frame(treeNum,nodeSizes)
treesAX<-ID[which(ID$nodeSizes==2),]
treesXB<-ID[which(ID$nodeSizes==26),]
tipOrder<-rev(c("Vespertilionid−related", "Emballonurid−related", "Phyllostomid−related", "Yinpterochiroptera", "Pholidota", "Carnivora", "Perissodactyla", "Cetartiodactyla", "Eulipotyphla", "Guinea_pig−related", "Squirrel−related", "Castorimorpha", "Anomaluromorpha", "Dipodidae", "Platacanthomyidae", "Spalacidae", "Calomyscidae", "Nesomyidae", "Cricetidae", "Muridae", "Lagomorpha", "Primates", "Dermoptera", "Scandentia", "Xenarthra", "Afrotheria", "Marsupialia", "Monotremata"))

	"Monotremata", "Marsupialia", "Afrotheria", "Xenarthra", "Eulipotyphla", "Pholidota", "Carnivora", "Perissodactyla", "Cetartiodactyla", "Yinpterochiroptera", "Phyllostomid−related", "Emballonurid−related", "Vespertilionid−related", "Primates", "Dermoptera", "Scandentia", "Lagomorpha", "Guinea_pig−related", "Squirrel−related", "Anomaluromorpha", "Castorimorpha", "Dipodidae", "Platacanthomyidae", "Spalacidae", "Calomyscidae", "Nesomyidae", "Cricetidae", "Muridae"))

patch_XB<-patchPhy_100[treesXB[,1]]
class(patch_XB)<-"multiPhylo"
patchXB_l<-lapply(patch_XB,ladderize)
patchXB<-lapply(patchXB_l,rotateConstr,constraint=tipOrder)

patch_AX<-patchPhy_100[treesAX[,1]]
class(patch_AX)<-"multiPhylo"
patchAX_l<-lapply(patch_AX,ladderize)
patchAX<-lapply(patchAX_l,rotateConstr,constraint=tipOrder)


# load the MCC for starting the plot 
patchPhy<-ladderize(read.tree(paste("MamPhy_BDvr_pcsFIXED_",bbone,"_MCC_target_28PATCHES_patchNames.tre",sep="")))
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

pdf(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_28PATCHES_LINEAR_redBlue_noTips.pdf",sep=""),width=7, height=10)
pdf(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_28PATCHES_LINEAR_redBlue_withTips_facing_clado.pdf",sep=""),width=10, height=7)
pdf(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_28PATCHES_LINEAR_redBlue_withTips_clade_Overlap.pdf",sep=""),width=10, height=7)
jpeg(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_28PATCHES_LINEAR_redBlue_withTips_clade_Overlap.jpg",sep=""),width=10, height=7, units="in",quality=100, res=600)

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


> plotTree(tree,ftype="off")
> plot.new()
> plot.window(xlim=c(-0.1,0.1),ylim=c(1, length(tree$tip.label)))
> par(cex=1)
> text(rep(0,length(tree$tip.label)), 1:length(tree$tip.label),tree$tip.label)
> plotTree(tree,ftype="off",direction="leftwards")

# linear
pdf(file=paste("MamPhy_BDvr_pcsFIXED_",bbone,"_MCC_target_28PATCHES_patchNames_CIRCLEplot.pdf",sep=""),width=5,height=5)
plot(patchPhy, type="phylogram", cex=0.7, font=2, label.offset=0.5,tip.color="white",show.tip.label=TRUE,no.margin=TRUE)
par(new=TRUE)
plot(patchPhy, type="phylogram", cex=0.7, font=2, label.offset=0.5, no.margin=TRUE)
axisPhylo(cex=0.6)

dev.off()


# MANUALLY AS A CIRCLE...
edgeCol1<-rgb(1,0,0,alpha=0.1)
edgeCol2<-rgb(0,0,1,alpha=0.05)

pdf(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_28PATCHES_Circular_redBlue.pdf",sep=""),width=7, height=10)
pdf(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_28PATCHES_Circular_redBlue_noTips.pdf",sep=""),width=7, height=10)
#jpeg(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_28PATCHES_Circular_redBlue.jpg",sep=""),width=7, height=10, units="in",quality=100, res=600)

plot(patchPhy, type="fan", cex=0.7, font=2, label.offset=0.5,open.angle=20,rotate.tree=33, tip.color="white",edge.col=edgeCol1,show.tip.label=TRUE,no.margin=TRUE)

dd<-0.89
cc<-0.95
draw.circle(0,0, radius=rootTime-0,lty=1,col=grey(cc), border=NA)
#draw.circle(0,0, radius=rootTime-66, lty=2,col=grey(cc), border=grey(cc))
draw.circle(0,0, radius=rootTime-66, lwd=2.5,lty=2,col=NA, border=grey(0.5))

for (j in 1:length(patchXB)){
	par(new=TRUE)
	plot(patchXB[[j]], type="fan", cex=0.7, font=2, label.offset=0.5,open.angle=20,rotate.tree=33, tip.color=grey(0.5, alpha=0.001),edge.col=edgeCol1,show.tip.label=TRUE,no.margin=TRUE)
}
par(new=TRUE)
plot(patchAX[[1]], type="fan", cex=0.7, font=2, label.offset=0.5,open.angle=20,rotate.tree=33, tip.color=grey(0.5, alpha=0.001),edge.col=edgeCol2,show.tip.label=TRUE,no.margin=TRUE)
for (j in 2:length(patchAX)){
	par(new=TRUE)
	plot(patchAX[[j]], type="fan", cex=0.7, font=2, label.offset=0.5,open.angle=20,rotate.tree=33, tip.color=grey(0.5, alpha=0.001),edge.col=edgeCol2,show.tip.label=TRUE,no.margin=TRUE)
}
axisPhylo(cex=0.6)

dev.off()


##
# same one but LIGHT, for working with.
pdf(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_28PATCHES_Circular_redBlue_LIGHT.pdf",sep=""),width=7, height=10)

plot(patchPhy, type="fan", cex=0.7, font=2, label.offset=0.5,open.angle=20,rotate.tree=33, tip.color="black",edge.col=edgeCol1,show.tip.label=TRUE,no.margin=TRUE)

dd<-0.89
cc<-0.95
draw.circle(0,0, radius=rootTime-0,lty=1,col=grey(cc), border=NA)
#draw.circle(0,0, radius=rootTime-66, lty=2,col=grey(cc), border=grey(cc))
draw.circle(0,0, radius=rootTime-66, lwd=2.5,lty=2,col=NA, border=grey(0.5))

for (j in 1:2){
	par(new=TRUE)
	plot(patchXB[[j]], type="fan", cex=0.7, font=2, label.offset=0.5,open.angle=20,rotate.tree=33, tip.color=grey(0.5, alpha=0.001),edge.col=edgeCol1,show.tip.label=TRUE,no.margin=TRUE)
}
par(new=TRUE)
for (j in 1:2){
	par(new=TRUE)
	plot(patchAX[[j]], type="fan", cex=0.7, font=2, label.offset=0.5,open.angle=20,rotate.tree=33, tip.color=grey(0.5, alpha=0.001),edge.col=edgeCol2,show.tip.label=TRUE,no.margin=TRUE)
}
axisPhylo(cex=0.6)

dev.off()






# NOW, plot the density plots all in an ARRAY, so you can swing them into INKSCAPE and move around.
###
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""))
colnames(cladesDR)<-c("tiplabel","gen","fam","famLabel","famNumLabel","famNumAll","ord","ordLabel1","ordLabel2","clade","cladeCommon","cladeCombo","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")
head(cladesDR)

# for REDOS of monkeys...
#######
OWM<-cladesDR[which(cladesDR$fam=="CERCOPITHECIDAE" | cladesDR$fam=="HYLOBATIDAE" |  cladesDR$fam=="HOMINIDAE"),"harmMeans"] 
NWM<-cladesDR[which(cladesDR$fam=="PITHECIIDAE" | cladesDR$fam=="ATELIDAE" |  cladesDR$fam=="AOTIDAE" |  cladesDR$fam=="CEBIDAE" |  cladesDR$fam=="CALLITRICHIDAE"),"harmMeans"] 

DR<-cladesDR$harmMeans
monkeyNames<-c("OWM","NWM")
monkeysDR<-list(OWM,NWM)

pdf(file=paste("DR-SUM_PLOT_harmMean_",bbone, "_REDOS-monkeys.pdf", sep=""), width=10, height=10, onefile=TRUE)
op <- par(mfrow = c(6,5),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)

for(i in 1:length(monkeyNames)){
CladeDR<-monkeysDR[[i]]

plot(density(CladeDR), col="dark grey", main="", bty="n", axes=F, xlim=range(DR),xlab="",ylab="", zero.line=FALSE)
polygon(density(CladeDR), col="dark grey", border="dark grey", bty="n")
x.tick <- quantile(DR, c(0.01,0.5,0.99))
axis(at=c(x.tick), labels=NA, side=1, line=0, cex=2.5, lwd=1, tck=0.01, cex.axis=1.5, mgp=c(1,2,0))
#axis(side=1, labels=FALSE, at=FALSE)
dens.rate <- density(CladeDR)$y
#axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=2.5, las=1, lwd=1, cex.axis=1.5, tck=-0.05, mgp=c(1,2,0))
seg.tick <- quantile(DR, c(0.01,0.5,0.99))
segments(seg.tick[[1]],0,seg.tick[[1]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")
segments(seg.tick[[2]],0,seg.tick[[2]],max(dens.rate)*1, lty=2, lwd=2,col="black")
segments(seg.tick[[3]],0,seg.tick[[3]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")

mtext(side=1,text=paste(monkeyNames[i],", N = ",length(CladeDR),sep=""), cex=0.6, line=0, adj=0)
}

par(op)
dev.off()

# for CLADES
#######
DR<-cladesDR$harmMeans
cladeNames1<-names(which(sort(table(cladesDR$clade),decreasing=TRUE) > 2))
cladeNames<-cladeNames1[c(10, 13, 18, 23, 29, 1 , 5 , 6 , 8 , 14, 11, 17, 25, 4 , 2 , 7 , 15, 24, 28, 9 , 12, 22, 27)]


pdf(file=paste("DR-SUM_PLOT_harmMean_",bbone, "_CLADES-23new_arrayedTogether.pdf", sep=""), width=6, height=6, onefile=TRUE)
op <- par(oma = c(1,1,1,1) + 0.1,
          mar = c(0,0,1,1) + 0.1)

layout(matrix(1:25,5,5, byrow=FALSE))#,widths=rep(1,0.10))
for(i in 1:length(cladeNames)){
CladeDR<-cladesDR[which(cladesDR$clade==cladeNames[i]),"harmMeans"]

plot(density(CladeDR), col="dark grey", main="", bty="n", axes=F, xlim=range(DR),xlab="",ylab="", zero.line=FALSE)
polygon(density(CladeDR), col="dark grey", border="dark grey", bty="n")
x.tick <- quantile(DR, c(0.01,0.5,0.99))
axis(at=c(x.tick), labels=NA, side=1, line=0, cex=2.5, lwd=1, tck=0.01, cex.axis=1.5, mgp=c(1,2,0))
#axis(side=1, labels=FALSE, at=FALSE)
dens.rate <- density(CladeDR)$y
#axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=2.5, las=1, lwd=1, cex.axis=1.5, tck=-0.05, mgp=c(1,2,0))
seg.tick <- quantile(DR, c(0.01,0.5,0.99))
segments(seg.tick[[1]],0,seg.tick[[1]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")
segments(seg.tick[[2]],0,seg.tick[[2]],max(dens.rate)*1, lty=2, lwd=2,col="black")
segments(seg.tick[[3]],0,seg.tick[[3]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")

mtext(side=1,text=paste(cladeNames[i],", n = ",length(CladeDR),sep=""), cex=0.6, line=0, adj=0)
}

par(op)
dev.off()

###
# for ORDERS
#######
DR<-cladesDR$harmMeans
ordNames<-names(which(sort(table(cladesDR$ord),decreasing=TRUE) > 2))

pdf(file=paste("DR-SUM_PLOT_harmMean_",bbone, "_ORDS-all23_arrayedTogether.pdf", sep=""), width=6, height=6, onefile=TRUE)
op <- par(oma = c(1,1,1,1) + 0.1,
          mar = c(0,0,1,1) + 0.1)

layout(matrix(1:25,5,5, byrow=FALSE))#,widths=rep(1,0.10))

for(i in 1:length(ordNames)){
CladeDR<-cladesDR[which(cladesDR$ord==ordNames[i]),"harmMeans"]

plot(density(CladeDR), col="dark grey", main="", bty="n", axes=F, xlim=range(DR),xlab="",ylab="", zero.line=FALSE)
polygon(density(CladeDR), col="dark grey", border="dark grey", bty="n")
x.tick <- quantile(DR, c(0.01,0.5,0.99))
axis(at=c(x.tick), labels=NA, side=1, line=0, cex=2.5, lwd=1, tck=0.01, cex.axis=1.5, mgp=c(1,2,0))
#axis(side=1, labels=FALSE, at=FALSE)
dens.rate <- density(CladeDR)$y
#axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=2.5, las=1, lwd=1, cex.axis=1.5, tck=-0.05, mgp=c(1,2,0))
seg.tick <- quantile(DR, c(0.01,0.5,0.99))
segments(seg.tick[[1]],0,seg.tick[[1]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")
segments(seg.tick[[2]],0,seg.tick[[2]],max(dens.rate)*1, lty=2, lwd=2,col="black")
segments(seg.tick[[3]],0,seg.tick[[3]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")

mtext(side=1,text=paste(ordNames[i],", n = ",length(CladeDR),sep=""), cex=0.6, font=1, line=0, adj=0)
}

par(op)
dev.off()

###
# for PATCHES
#######
DR<-cladesDR$harmMeans
patchNames<-names(which(sort(table(cladesDR$PC),decreasing=TRUE) > 2))

pdf(file=paste("DR-SUM_PLOT_harmMean_",bbone, "_PCS-all26_arrayedTogether.pdf", sep=""), width=10, height=10, onefile=TRUE)
op <- par(mfrow = c(6,5),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)

for(i in 1:length(patchNames)){
CladeDR<-cladesDR[which(cladesDR$PC==patchNames[i]),"harmMeans"]

plot(density(CladeDR), col="dark grey", main="", bty="n", axes=F, xlim=range(DR),xlab="",ylab="", zero.line=FALSE)
polygon(density(CladeDR), col="dark grey", border="dark grey", bty="n")
x.tick <- quantile(DR, c(0.01,0.5,0.99))
axis(at=c(x.tick), labels=NA, side=1, line=0, cex=2.5, lwd=1, tck=0.01, cex.axis=1.5, mgp=c(1,2,0))
#axis(side=1, labels=FALSE, at=FALSE)
dens.rate <- density(CladeDR)$y
#axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=2.5, las=1, lwd=1, cex.axis=1.5, tck=-0.05, mgp=c(1,2,0))
seg.tick <- quantile(DR, c(0.01,0.5,0.99))
segments(seg.tick[[1]],0,seg.tick[[1]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")
segments(seg.tick[[2]],0,seg.tick[[2]],max(dens.rate)*1, lty=2, lwd=2,col="black")
segments(seg.tick[[3]],0,seg.tick[[3]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")

mtext(side=1,text=paste(patchNames[i],", N = ",length(CladeDR),sep=""), cex=0.6, font=2, line=0, adj=0)
}

par(op)
dev.off()





#=========
# Getting the density plots themselves...
#=========

library(ape)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")

#clades<-c("MAMMALIA", "MARSUPIALS", "PLACENTALS", "AFROTHERIA", "XENARTHRA", "CARNIVORA", "CETARTIODACTYLA", "EULIPOTYPHLA",  "CHIROPTERA", "YANGOCHIROPTERA", "YINPTEROCHIROPTERA", "PRIMATES", "STREPSIRRHINI", "HAPLORRHINI", "NWMonkeys", "OWMonkeys", "LAGOMORPHA", "RODENTIA", "GP-RELATED", "SQ-RELATED","MO-RELATED")

bbone<- "NDexp" # "FBD"

cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdGenesSampPC_DR-SUMMARY-expanded_",bbone,".txt",sep=""))
head(cladesDR)
colnames(cladesDR)<-c("tiplabel","gen","fam","ord","clade","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")

DR<-cladesDR$harmMeans

ordNames<-names(which(sort(table(cladesDR$ord),decreasing=TRUE) > 2))
famNames<-names(which(sort(table(cladesDR$fam),decreasing=TRUE) > 2))
patchNames<-names(which(sort(table(cladesDR$PC),decreasing=TRUE) > 2))
cladeNames<-names(which(sort(table(cladesDR$clade),decreasing=TRUE) > 2))

cladeNames<-cladeNames1[c(1,5,6,18,2,4,8,20,3,7,17,12,9,15,10,13,16,21,22,24,25,11,14,19,23,26,27)]


pdf(file=paste("DR-SUM_PLOT_harmMean_",bbone, "_ORDS-all23.pdf", sep=""), width=11, height=8.5, onefile=TRUE)
for(i in 1:length(ordNames)){
ordDR<-cladesDR[which(cladesDR$ord==ordNames[i]),"harmMeans"]

plot(density(ordDR), col="dark grey", main="", bty="n", axes=F, xlim=range(DR),xlab="",ylab="", zero.line=FALSE)
polygon(density(ordDR), col="dark grey", border="dark grey", bty="n")
x.tick <- quantile(DR, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick), labels=NA, side=1, line=0, cex=2.5, lwd=1, tck=-0.05, cex.axis=1.5, mgp=c(1,2,0))
#dens.rate <- density(CladeDR)$y
#axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=2.5, las=1, lwd=1, cex.axis=1.5, tck=-0.05, mgp=c(1,2,0))
segments(mean(DR),0,mean(DR),100, lty=2, lwd=3,col="black")

title(main=paste(ordNames[i],", N = ",length(ordDR),sep=""), cex.main=2)
}
dev.off()


pdf(file=paste("DR-SUM_PLOT_harmMean_",bbone, "_PCS-all26.pdf", sep=""), width=11, height=8.5, onefile=TRUE)
for(i in 1:length(patchNames)){
patchDR<-cladesDR[which(cladesDR$PC==patchNames[i]),"harmMeans"]

plot(density(patchDR), col="dark grey", main="", bty="n", axes=F, xlim=range(DR),xlab="",ylab="", zero.line=FALSE)
polygon(density(patchDR), col="dark grey", border="dark grey", bty="n")
x.tick <- quantile(DR, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick), labels=NA, side=1, line=0, cex=2.5, lwd=1, tck=-0.05, cex.axis=1.5, mgp=c(1,2,0))
#dens.rate <- density(CladeDR)$y
#axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=2.5, las=1, lwd=1, cex.axis=1.5, tck=-0.05, mgp=c(1,2,0))
segments(mean(DR),0,mean(DR),100, lty=2, lwd=3,col="black")

title(main=paste(patchNames[i],", N = ",length(patchDR),sep=""), cex.main=2)
}
dev.off()


pdf(file=paste("DR-SUM_PLOT_harmMean_",bbone, "_CLADES-all30.pdf", sep=""), width=11, height=8.5, onefile=TRUE)
for(i in 1:length(cladeNames)){

CladeDR<-cladesDR[which(cladesDR$clade==cladeNames[i]),"harmMeans"]

plot(density(CladeDR), col="dark grey", main="", bty="n", axes=F, xlim=range(DR),xlab="",ylab="", zero.line=FALSE)
polygon(density(CladeDR), col="dark grey", border="dark grey", bty="n")
x.tick <- quantile(DR, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick), labels=NA, side=1, line=0, cex=2.5, lwd=1, tck=-0.05, cex.axis=1.5, mgp=c(1,2,0))
#dens.rate <- density(CladeDR)$y
#axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=2.5, las=1, lwd=1, cex.axis=1.5, tck=-0.05, mgp=c(1,2,0))
segments(mean(DR),0,mean(DR),100, lty=2, lwd=3,col="black")

title(main=paste(cladeNames[i],", N = ",length(CladeDR),sep=""), cex.main=2)
}
dev.off()


###
# JUST the desired clades

Rodents
- 

ordNames[1], patchNames[9], patchNames[10], 
patchNames[1]+patchNames[2]+patchNames[13]+patchNames[17]+patchNames[18]+patchNames[21]+patchNames[23]+patchNames[25]

"RODENTIA", "Squirrel-related", "Guinea pig-related", "Mouse-related"


MARSUPIALS
- the 5 orders, great






######
# Get densiTree of the 100 mammal trees...
library(ape)
library(phangorn)
library(phytools)
library(phyloch)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
# specify the backbone to use - FBD or NDexp
bbone <- "NDexp" # "FBD"

# sample of 100 trees (different backbones)
mamPhy_100<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_nexus.trees",sep=""))
mamPhy100<-lapply(mamPhy_100,ladderize)
class(mamPhy100)<-"multiPhylo"

# MCC tree
NDexp0<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target.tre")
NDexp1<-drop.tip(NDexp0,"_Anolis_carolinensis")
NDexp_mam<-ladderize(NDexp1)

# 100trees: find which trees are A(X+B) vs (A+X)B
nodeSizes<-c()
for(i in 1:length(mamPhy100)){
  tree<-mamPhy100[[i]]
  node<-getMRCA(tree,tip=c("Amblysomus_corriae_CHRYSOCHLORIDAE_AFROSORICIDA","Cabassous_centralis_DASYPODIDAE_CINGULATA"))
  nodeTips<-tree$tip.label[Descendants(tree,node)[[1]]]  
  nodeSizes[i]<-length(nodeTips)
}
# nodeSizes
# 125 5544 
#  59   41 
treeNum<-c(1:length(mamPhy100))
ID<-cbind.data.frame(treeNum,nodeSizes)
treesAX<-ID[which(ID$nodeSizes==125),]
treesXB<-ID[which(ID$nodeSizes==5544),]

mamPhy100z <- .compressTipLabel(mamPhy100) # moves tiplabels to one element, nice...!!

# with densiTree (phangorn)::
# red-blue together
pdf(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_DensityTree_100_a01_blueAX-redXB-Together_densi.pdf",sep=""),width=7, height=10)
densiTree(mamPhy100[treesXB[,1]], type="phylogram",col="red",alpha=0.01,  consensus=NDexp_mam, no.margin=FALSE, cex=0.9, label.offset=0.4, direction="rightwards", show.tip.label=FALSE)
par(new=TRUE)
densiTree(mamPhy100[treesAX[,1]], type="phylogram",col="blue",alpha=0.01,  consensus=NDexp_mam, no.margin=FALSE, cex=0.9, label.offset=0.4, direction="rightwards", show.tip.label=FALSE)
dev.off()

# with densityTree (phytools):
pdf(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_DensityTree_100_a01_blueAX-redXB-Together_densi.pdf",sep=""),width=7, height=10)
densityTree(mamPhy100[treesXB[,1]], type="phylogram",ftype="off", fsize=0.8,alpha=0.01, use.edge.length=FALSE, colors="red")
axisPhylo(side=1)
#par(new=TRUE)
densityTree(mamPhy100[treesAX[,1]], type="phylogram",ftype="off", offset=2,fsize=0.8,alpha=0.01, use.edge.length=FALSE,direction="rightwards",colors="blue")
dev.off()



densiTree<-function (x, type = "cladogram", alpha = 1/length(x), consensus = NULL, 
    optim = FALSE, scaleX = FALSE, col = 1, width = 1, cex = 0.8, 
    ...) 
{
    if (!inherits(x, "multiPhylo")) 
        stop("x must be of class multiPhylo")
    compressed <- ifelse(is.null(attr(x, "TipLabel")), FALSE, 
        TRUE)
    if (is.null(consensus)) 
        consensus <- superTree(x)
    if (inherits(consensus, "multiPhylo")) 
        consensus = consensus[[1]]
    type <- match.arg(type, c("phylogram", "cladogram"))
    consensus = reorder(consensus, "postorder")
    e2 = reorder(consensus)$edge[, 2]
    nTip = as.integer(length(consensus$tip))
    tiporder = e2[e2 <= nTip]
    maxBT = max(getAges(x))
    if (scaleX) 
        maxBT = 1
    label = rev(pretty(c(maxBT, 0)))
    maxBT = max(label)
    xy = plotPhyloCoor(consensus, ...)
    yy = xy[, 2]
    plot.new()
    tl = which.max(nchar(consensus$tip.label))
    sw <- strwidth(consensus$tip.label[tl], cex = cex) * 1.1
    plot.window(xlim = c(0, 1 + sw), ylim = c(0, nTip + 1))
    axis(side = 1, at = seq(0, 1, length.out = length(label)), 
        labels = label)
    text(x = rep(1, Ntip(consensus)), y = yy[1:nTip], labels = consensus$tip.label, 
        pos = 4, cex = cex)
    tip.order = yy[1:nTip]
    for (treeindex in 1:length(x)) {
        tmp <- reorder(x[[treeindex]], "postorder")
        if (!compressed) 
            tip.order <- match(tmp$tip.label, consensus$tip.label)
        xy <- plotPhyloCoor(tmp, tip.order = tiporder, ...)
        xx = xy[, 1]
        yy = xy[, 2]
        if (scaleX) 
            xx <- xx/max(xx)
        else xx <- xx/maxBT
        xx <- xx + (1 - max(xx))
        e1 = tmp$edge[, 1]
        e2 = tmp$edge[, 2]
        if (type == "cladogram") 
            cladogram.plot(tmp$edge, xx, yy, edge.color = adjustcolor(col, 
                alpha.f = alpha), edge.width = width, edge.lty = 1)
        if (type == "phylogram") {
            Ntip <- min(e1) - 1L
            Nnode <- tmp$Nnode
            phylogram.plot(tmp$edge, Ntip, Nnode, xx, yy, TRUE, 
                edge.color = adjustcolor(col, alpha.f = alpha), 
                edge.width = width, 1)
        }
    }
}







##########
for(i in 1:length(ordNames)){

data<-read.table(paste("DR-SUM_Exp_harmMean_", ordNames[i], ".txt", sep=""), header=FALSE)

subsetDR<-data$V2

pdf(file=paste("DR-SUM_PLOT_Exp_harmMean_", i, "_", clades[i], ".pdf", sep=""), width=11, height=8.5, onefile=TRUE)

plot(density(subsetDR), col="dark grey", main="", bty="n", axes=F, xlim=c(0,0.78))
polygon(density(subsetDR), col="dark grey", border="dark grey", bty="n")

x.tick <- quantile(subsetDR, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick), labels=c(0,round(x.tick,2)), side=1, line=1.3, cex=2.5, lwd=1, tck=-0.05, cex.axis=1.5, mgp=c(1,2,0))

dens.rate <- density(subsetDR)$y
axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=2.5, las=1, lwd=1, cex.axis=1.5, tck=-0.05, mgp=c(1,2,0))

segments(0.2246393,0,0.2246393,5, lty=2, lwd=3,col="red")

title(main=clades[i], cex.main=3)

dev.off()

}


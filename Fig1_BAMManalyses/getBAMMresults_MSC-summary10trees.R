library(BAMMtools)
library(coda)
library(phytools)
library(ape)

folders<-c(1:10)

bbone<-"NDexp" # "FBD"

#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
#tree<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_5911sp.tre")
#write.tree(tree, "MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_5911sp_newick.tre")
#tree2<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp.tre")
#write.tree(tree2, "MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")

#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/_FBD_10trees_MSC-focus")
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/_NDexp_10trees_MSC-focus")

###
# Tree-wide rate and shift data summary

RES_SUM<-vector("list",length(folders))
for (i in 1:length(folders)){
	RES_SUM[[i]]<-read.table(file=paste("results.treeSUMMARY.",folders[i],".txt",sep=""), header=TRUE)
}

treewideRES<-mapply(cbind,RES_SUM)

# take the mean of each column
treewideRES_means

RES<-vector(length=length(treewideRES[,1]))
for (i in 1:length(treewideRES[,1])){
	RES[i]<-mean(unlist(treewideRES[i,]))
}
names(RES)<-rownames(treewideRES)

RES_final<-setNames(as.data.frame(round(RES,digits=4)),"mean")

write.table(RES_final,file=paste("SUM10_results.treeSUMMARY.",bbone,".txt",sep=""))

###
# Tree-wide rate and shift data summary
RES_SUM_MSC<-vector("list",length(folders))
for (i in 1:length(folders)){
	RES_SUM_MSC[[i]]<-read.table(file=paste("results.MSC.nodeRATES.",folders[i],".txt",sep=""), header=TRUE)
}

resMSC<-vector("list",length(folders))
for (i in 1:length(folders)){
	# summarize to get divRate mean of shift clade and non-shift clade= background rate
	res<-RES_SUM_MSC[[i]]

	cladeDiv_mean<-res[,"Lam_mean"]-res[,"Mu_mean"] 
	noncladeDiv_mean<-res[,"Lam_nonmean"]-res[,"Mu_nonmean"]
	factor <- cladeDiv_mean/noncladeDiv_mean
	incDec<-rep(NA,length(res[,1]))
	incDec[factor < 1] <- "downShift" 
	incDec[factor > 1] <- "upShift" 
	tree<-rep(i,length(res[,1]))
	divs<-cbind.data.frame(cladeDiv_mean,noncladeDiv_mean,factor,incDec,tree)

	resMSC[[i]]<-cbind.data.frame(res[,1:5], divs)
}

all10MSC<-do.call(rbind,resMSC)
write.table(all10MSC,file=paste("SUM10_results.MSC.nodeRATES.",bbone,".txt",sep=""))


###
# per Order rate data summary

RES_SUM_ord<-vector("list",length(folders))
for (i in 1:length(folders)){
	RES_SUM_ord[[i]]<-read.table(file=paste("results.ordRATES.",folders[i],".txt",sep=""), header=TRUE)
}

RES_SUM_ord_ALL<-do.call(rbind,RES_SUM_ord) # concatenates the tables...

means<-vector("list",length(table(RES_SUM_ord_ALL$ord)))
sigCounts<-vector("list",length(table(RES_SUM_ord_ALL$ord)))
for (i in 1:length(table(RES_SUM_ord_ALL$ord))){
	ordRes<-RES_SUM_ord_ALL[which(RES_SUM_ord_ALL$ord==names(table(RES_SUM_ord_ALL$ord))[i]),]
	signifs<-data.frame(matrix(NA, nrow = 10, ncol = 2))
	for (j in 1:length(folders)){
		if (ordRes$Lam_mean[j] < ordRes$Lam_nonmean[j]){
		signifs[j,1]<-ordRes$Lam_up95[j] < ordRes$Lam_nonlow95[j]
		} else
		signifs[j,1]<-ordRes$Lam_low95[j] > ordRes$Lam_nonup95[j]
		if (ordRes$Mu_mean[j] < ordRes$Mu_nonmean[j]){
		signifs[j,2]<-ordRes$Mu_up95[j] < ordRes$Mu_nonlow95[j]
		} else
		signifs[j,2]<-ordRes$Mu_low95[j] > ordRes$Mu_nonup95[j]
	}
	sigCounts[[i]]<-colSums(signifs)
	means[[i]]<-colSums(ordRes[2:13])/10
}
allOrds_meanRates<-t(mapply(cbind,means))
colnames(allOrds_meanRates)<-names(means[[1]])
rownames(allOrds_meanRates)<-names(table(RES_SUM_ord_ALL$ord))
allOrds_meanRates

allOrds_sigCounts<-t(mapply(cbind,sigCounts))
colnames(allOrds_sigCounts)<-c("Lam_diff","Mu_diff")
rownames(allOrds_sigCounts)<-names(table(RES_SUM_ord_ALL$ord))
allOrds_sigCounts

finalResOrds<-cbind(allOrds_sigCounts,allOrds_meanRates)
write.table(finalResOrds,file=paste("SUM10_results.ordRATES.",bbone,".txt",sep=""))


###########
#############
# Load back in processed SHIFT data to R
library(phangorn)
library(phytools)
bbone <- "NDexp" #"FBD"

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/")
# load shifts
shiftDat<-read.csv(file=paste("BAMM_summary_",bbone,"_24shifts_R.csv",sep=""), header=TRUE)

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
# load MCC tree
mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")

# plot the mamMCC for easy viewing.
pdf(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick_linear.pdf",width=22,height=150)
plot(mamMCC, cex=0.2)
axisPhylo()
dev.off()

# get node nums
nodes<-c()
lengths<-c()
tips<-vector("list",length(shiftDat$to))
for (i in 1:length(shiftDat$to)){
	taxa<-c(as.vector(shiftDat$from[i]),as.vector(shiftDat$to[i]))
	if (taxa[1]==taxa[2])
		{nodes[i]<-findMRCA(mamMCC,tips=taxa)} 
	else 
		{nodes[i]<-getMRCA(mamMCC,taxa)}
	lengths[i]<-length(Descendants(mamMCC,nodes[i])[[1]])
	tips[[i]]<-mamMCC$tip.label[Descendants(mamMCC,nodes[i])[[1]]]
}

toPlot<-cbind.data.frame(nodes,lengths,shiftDat$numTaxa,as.vector(shiftDat$CLADE),as.vector(shiftDat$from), as.vector(shiftDat$to), stringsAsFactors=FALSE)
colnames(toPlot)<-c("nodes","lengths","numE","CLADE","from","to")

misMatch<-toPlot[which(toPlot[,2] != toPlot[,3]),]

# Second time:
   nodes lengths numE                                      CLADE
5  11337     125 5544                                Placentalia
14  7824      92   93    Phyllostomatidae sub1 (Stenodermatinae)
29 10098     196  206     Muridae sub2 (Murinae: Apomys-Melomys)
31  9759     104  109 Muridae sub4 (Murinae: Rattus-Srilankamys)
32  9542      59   80                 Muridae sub5 (Gerbillinae)
33  9550      51   50      Muridae sub6 (Gerbillinae: Gerbillus)
                                              from
5      Dasypus_septemcinctus_DASYPODIDAE_CINGULATA
14 Uroderma_magnirostrum_PHYLLOSTOMIDAE_CHIROPTERA
29          Chiropodomys_muroides_MURIDAE_RODENTIA
31               Berylmys_bowersi_MURIDAE_RODENTIA
32            Taterillus_congicus_MURIDAE_RODENTIA
33                Gerbillus_dunni_MURIDAE_RODENTIA
                                                      to
5  Chrysochloris_stuhlmanni_CHRYSOCHLORIDAE_AFROSORICIDA
14             Sturnira_bidens_PHYLLOSTOMIDAE_CHIROPTERA
29                Apomys_gracilirostris_MURIDAE_RODENTIA
31                 Srilankamys_ohiensis_MURIDAE_RODENTIA
32                  Ammodillus_imbellis_MURIDAE_RODENTIA
33                    Microdillus_peeli_MURIDAE_RODENTIA

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/")
shiftDatEdit<-as.matrix(read.csv(file=paste("BAMM_summary_",bbone,"_24shifts_R.csv",sep=""), header=TRUE))

shiftDatEdit[which(shiftDatEdit[,12]==misMatch$to[1]),"to"]<- "Apomys_gracilirostris_MURIDAE_RODENTIA" # a boreo
#shiftDat[which(shiftDat$to==misMatch$to[2]),"to"] # is OK
#shiftDat[which(shiftDat$to==misMatch$to[3]),"to"]<- # is OK
#shiftDat[which(shiftDat$to==misMatch$to[4]),"to"]<- # is OK
shiftDatEdit[which(shiftDatEdit[,12]==misMatch$to[5]),"to"]<- "Meriones_rex_MURIDAE_RODENTIA" # get the whole gerbillinae clade
#shiftDat[which(shiftDat$to==misMatch$to[6]),"to"]<- #is OK


# First time:
      nodes   lengths                                                         
 [1,] "11337" "125"   "5544" "Dasypus_septemcinctus_DASYPODIDAE_CINGULATA"   "Chrysochloris_stuhlmanni_CHRYSOCHLORIDAE_AFROSORICIDA"
 [2,] "7824"  "92"    "93"   "Uroderma_magnirostrum_PHYLLOSTOMIDAE_CHIROPTERA"
 [3,] "6267"  "122"   "128"  "Sorex_hoyi_SORICIDAE_EULIPOTYPHLA"              
 [4,] "11090" "50"    "51"   "Paradipus_ctenodactylus_DIPODIDAE_RODENTIA"     
 [5,] "10098" "196"   "206"  "Chiropodomys_muroides_MURIDAE_RODENTIA"         
 [6,] "9759"  "104"   "109"  "Berylmys_bowersi_MURIDAE_RODENTIA"              
 [7,] "9542"  "59"    "80"   "Taterillus_congicus_MURIDAE_RODENTIA"           
 [8,] "9550"  "51"    "50"   "Gerbillus_dunni_MURIDAE_RODENTIA"               
 [9,] "9646"  "16"    "23"   "Lophuromys_medicaudatus_MURIDAE_RODENTIA"       
            
shiftDatEdit<-as.matrix(read.csv(file=paste("BAMM_summary_",bbone,"_35shifts_R.csv",sep=""), header=TRUE))

shiftDatEdit[which(shiftDatEdit[,12]==misMatch$to[1]),"to"]<- "Apomys_gracilirostris_MURIDAE_RODENTIA" # a boreo
#shiftDat[which(shiftDat$to==misMatch$to[2]),"to"] # is OK
#shiftDat[which(shiftDat$to==misMatch$to[3]),"to"]<- # is OK
shiftDatEdit[which(shiftDatEdit[,12]==misMatch$to[4]),"to"]<- "Salpingotulus_michaelis_DIPODIDAE_RODENTIA" # full dipodidae
#shiftDat[which(shiftDat$to==misMatch$to[5]),"to"]<- # is OK
#shiftDat[which(shiftDat$to==misMatch$to[6]),"to"]<- # is OK
shiftDatEdit[which(shiftDatEdit[,12]==misMatch$to[7]),"to"]<- "Meriones_rex_MURIDAE_RODENTIA" # get the whole gerbillinae clade
#shiftDat[which(shiftDat$to==misMatch$to[8]),"to"]<- #is OK
shiftDatEdit[which(shiftDatEdit[,12]==misMatch$to[9]),"to"]<- "Lophuromys_luteogaster_MURIDAE_RODENTIA" # get the whole genus

shiftDatEdit2<-as.data.frame(shiftDatEdit)

# now CHECK it.
# get node nums
nodes<-c()
lengths<-c()
tips<-vector("list",length(shiftDatEdit2$to))
for (i in 1:length(shiftDatEdit2$to)){
	taxa<-c(as.vector(shiftDatEdit2$from[i]),as.vector(shiftDatEdit2$to[i]))
	if (taxa[1]==taxa[2])
		{nodes[i]<-findMRCA(mamMCC,tips=taxa)} 
	else 
		{nodes[i]<-getMRCA(mamMCC,taxa)}
	lengths[i]<-length(Descendants(mamMCC,nodes[i])[[1]])
	tips[[i]]<-mamMCC$tip.label[Descendants(mamMCC,nodes[i])[[1]]]
}

toPlot2<-cbind.data.frame(nodes,lengths,shiftDat$numTaxa,as.vector(shiftDat$CLADE),as.vector(shiftDat$from), as.vector(shiftDat$to), stringsAsFactors=FALSE)
colnames(toPlot2)<-c("nodes","lengths","numE","CLADE","from","to")

toPlotFIN<-cbind.data.frame(nodes,shiftDat$ID,shiftDat$avgFactor,shiftDat$avgIncDec,as.vector(shiftDat$CLADE))
colnames(toPlotFIN)<-c("nodes","ID","factor","shift","CLADE")
write.table(toPlotFIN,file=paste(bbone,"_BAMMshifts_24shifts_ReadyToPlot.txt",sep=""))

#David BapstJanuary 29, 2012 at 12:39 PM
#When one wants to deal with the same 'node' on different trees, one define nodes as being the same based on their tip descendants. prop.part() in ape becomes indispensible for this, as one can use the output to find the nodes which match on multiple trees. Obviously, we may want to define nodes somewhat differently in other cases. (For example, consider two phylogenies which do not include the same plant species but still each have a node which includes all angiosperms.)

###
# LOAD in the FAMILY names as NUMBER equivalencies...
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")
famNums<-read.table("mamPhy_tipGenFam_FamNum_Ord_5911sp_NDexpMCC.txt")
colnames(famNums)<-c("tiplabel","gen","fam","num","ord")

famNumOnly<-cbind.data.frame(famNums$fam,famNums$num)
xx<-unique(famNumOnly)
colnames(xx)<-c("fam","num")

famCount<-as.data.frame(table(famNums$fam))
colnames(famCount)<-c("fam","Freq")

zz<-merge(famCount,xx,by="fam",all.x=TRUE)


numCount<-as.data.frame(table(famNums$num))

numCount$var

famNums[which(famNums$num==numCount$var),"fam"]


dd<-merge(numCount,famCount,by="Freq")

###
# READY to plot on the main tree
###
# Lets do a LINEAR one first.
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")

bbone<-"NDexp"
toPlotFIN<-read.table(file=paste(bbone,"_BAMMshifts_24shifts_ReadyToPlot.txt",sep=""), header=TRUE)

up<-rgb(1,0,0,0.5) #red
down<-rgb(0,0,1,0.5) #blue
upDown<-grey(0.5, alpha = 0.5)
bgcols<-rep(NA,length(toPlotFIN$shift))

bgcols[which(toPlotFIN$shift=="upShift")]<-up
bgcols[which(toPlotFIN$shift=="downShift")]<-down
bgcols[which(toPlotFIN$shift=="up-or-down")]<-upDown

bgsizes<-toPlotFIN$factor

#pdf(file="NDexp_MCC_target_5911sp_linear_with35shifts.pdf",width=22,height=150)
pdf(file="NDexp_MCC_target_5911sp_linear_with24shifts_label.pdf",width=22,height=150)
plot(ladderize(mamMCC), cex=0.2)
nodelabels(pch=21, node=toPlotFIN$node, frame = "circle", bg = bgcols, cex= bgsizes)
nodelabels(text=toPlotFIN$ID, node=toPlotFIN$node, adj=c(1.5,1.5), font=2, frame = "none", col = "black", cex= 1.5)
nodelabels(text=toPlotFIN$CLADE, node=toPlotFIN$node, adj=c(1.5,2.5), font=1, frame = "none", col = "black", cex= 1.5)
axisPhylo()
dev.off()

# now as a RADIAL plot
#pdf(file="NDexp_MCC_target_5911sp_FAN_with35shifts.pdf",width=30,height=30)
pdf(file="NDexp_MCC_target_5911sp_FAN_with24shifts_label.pdf",width=30,height=30)
plot(ladderize(mamMCC), cex=0.1,type="fan")
nodelabels(pch=21, node=toPlotFIN$node, frame = "circle", bg = bgcols, cex= bgsizes)
nodelabels(text=toPlotFIN$ID, node=toPlotFIN$node, adj=c(1.5,1.5), font=2, frame = "none", col = "black", cex= 1.5)
nodelabels(text=toPlotFIN$CLADE, node=toPlotFIN$node, adj=c(1.5,2.5), font=1, frame = "none", col = "black", cex= 1.5)
dev.off()

##########
# now as a DR plot !!!
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
source("DR_circleTree_functions.R")
library(RColorBrewer)

# plotting tree
bbone <- "NDexp" #"FBD"
mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")
plottree <- ladderize(mamMCC)

# for the BAMM shifts
toPlotFIN<-read.table(file=paste(bbone,"_BAMMshifts_24shifts_ReadyToPlot.txt",sep=""), header=TRUE)

up<-rgb(1,0,0,0.5) #red
down<-rgb(0,0,1,0.5) #blue
upDown<-grey(0.5, alpha = 0.5)
bgcols<-rep(NA,length(toPlotFIN$shift))

bgcols[which(toPlotFIN$shift=="upShift")]<-up
bgcols[which(toPlotFIN$shift=="downShift")]<-down
bgcols[which(toPlotFIN$shift=="up-or-down")]<-upDown

bgsizes<-toPlotFIN$factor

# read in the data for colouring branches
#cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdGenesSampPC_DR-SUMMARY-expanded_",bbone,".txt",sep=""))
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
#colnames(cladesDR)<-c("tiplabel","gen","fam","famLabel","famNumLabel","famNumAll","ord","ordNums","ordNumsSelect","cladeOrdNums","ordLabel1","ordLabel2","clade","cladeCommonAll","cladeCommonSubs","cladeCommonSubs2","cladeCommonSubs2Nums","cladeCombo","higher","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")
head(cladesDR)

#===
# FAM labels
# ===
# if truncate all family names...
#fams <- cladesDR$fam[match(plottree$tip.label,cladesDR$tiplabel)]
#famLow<-tolower(fams)
firstup <- function(x) {
	substr(x, 1, 1) <- toupper(substr(x, 1, 1))
	x}
#famFirst<-firstup(famLow)
#fam3<-substr(famFirst,1,3)
#length(names(table(fam3))) # only gives 135 fams-- so some redundancy...
#fam8<-substr(famFirst,1,8)
#length(names(table(fam8))) # now gives 162. Myrmecob and Myrmecop are the last most signif.

# if all full fam names
#fams <- cladesDR$fam[match(plottree$tip.label,cladesDR$tiplabel)]
# if some fams as Nums (< 10 species)
famUp <- cladesDR$famLabel[match(plottree$tip.label,cladesDR$tiplabel)]
famLow<-tolower(famUp)
fams <- firstup(famLow)
famNums <- cladesDR$famNumLabel[match(plottree$tip.label,cladesDR$tiplabel)]

# if all fams as nums
famNumAll <- cladesDR$famNumAll[match(plottree$tip.label,cladesDR$tiplabel)]

#===
# ORD labels
# ===
# if full ord
#ords <- cladesDR$ord[match(plottree$tip.label,cladesDR$tiplabel)]
# if some abbreviated ords (2 chars)
ords1 <- cladesDR$ordLabel1[match(plottree$tip.label,cladesDR$tiplabel)]
ords2 <- cladesDR$ordLabel2[match(plottree$tip.label,cladesDR$tiplabel)]
ords <- cladesDR$ord[match(plottree$tip.label,cladesDR$tiplabel)]
ordNums<- cladesDR$ordNums[match(plottree$tip.label,cladesDR$tiplabel)]
ordNumsSelect<- cladesDR$ordNumsSelect[match(plottree$tip.label,cladesDR$tiplabel)]

ordsAbrev<-substr(ords,1,4)


# =====
# CLADE labels
# =====
#PCS <- cladesDR$PC[match(plottree$tip.label,cladesDR$tiplabel)]
#xx<-strsplit(as.vector(PCS),"[_]")
#PC_names<-c()
#for (i in 1:length(xx)){
#	PC_names[i]<-xx[[i]][2]
#}
clades <- cladesDR$clade[match(plottree$tip.label,cladesDR$tiplabel)]
cladesCommon <- cladesDR$cladeCommonAll[match(plottree$tip.label,cladesDR$tiplabel)]
cladesComSubs <- cladesDR$cladeCommonSubs[match(plottree$tip.label,cladesDR$tiplabel)]
cladesComSubs2 <- cladesDR$cladeCommonSubs2[match(plottree$tip.label,cladesDR$tiplabel)]
cladesComSubs2Nums <- cladesDR$cladeCommonSubs2Nums[match(plottree$tip.label,cladesDR$tiplabel)]
cladeOrdNums<- cladesDR$cladeOrdNums[match(plottree$tip.label,cladesDR$tiplabel)]

higherTax<- cladesDR$higher[match(plottree$tip.label,cladesDR$tiplabel)]

# load DR
DR <- cladesDR[,"harmMeans"]
names(DR) <- cladesDR$tiplabel

DR_ordered <- DR[match(plottree$tip.label,names(DR))]
DR_anc <- ace(DR_ordered, plottree, method="pic", scaled=TRUE)$ace

# Match ancestors totree edges
    match.anc <- match(as.numeric(names(DR_anc)), plottree$edge[,2])[-1]

    # Assign rates to each internal node
    reconRate <- vector(mode="numeric", length=length(plottree$edge[,2]))
    reconRate[match.anc] <- DR_anc[2:5910] #[2:7237]

    # Assign rates to tips
    tip.idx <- sort(plottree$tip.label, index.return=TRUE)$ix

    reconRate[match(tip.idx, plottree$edge[,2])] <- DR[sort(names(DR))]

    # Create colour palette
    reconColors <- reconRate

	#cols<-rev(colorRampPalette(c('red','orange','yellow','green','deepskyblue1','blue','black'))(100))
	
	#library(viridis)
	#cols<-viridis(100)
	cols<-rich.colors(100)
	range <- quantile(reconRate, seq(0,1, 0.01))[2:101]
    for (i in 1:100) {
        if (i==100) {range[i] <- range[i]+0.1}
        if (i==1) {reconColors[which(reconRate <= range[i])] <- cols[i] } 
        if (i > 1) {reconColors[which((reconRate <= range[i]) & (reconRate > range[i-1]))] <- cols[i]}
        }

	plot(1:100, col = rev(colorRampPalette(c('red','orange','yellow','green','deepskyblue1','blue','black'))(100)), pch = 16, cex = 3)
	plot(1:100, col = palette(rich.colors(100)), pch = 16, cex = 3)
	plot(1:100, col = palette(viridis(100)), pch = 16, cex = 3)
	
	#plot(1:100, col = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100), pch = 16, cex = 3)
   	#plot(1:100, col = palette(rich.colors(100)), pch = 16, cex = 3)
   	#plot(1:100, col = topo.colors(100), pch = 16, cex = 3)

    #range <- quantile(reconRate, seq(0,1, 0.01))[2:101]
    #for (i in 1:100) {
    #    if (i==100) {range[i] <- range[i]+0.1}
    #    if (i==1) {reconColors[which(reconRate <= range[i])] <- palette(rich.colors(100))[i] } 
    #    if (i > 1) {reconColors[which((reconRate <= range[i]) & (reconRate > range[i-1]))] <- palette(rich.colors(100))[i]}
    #    }

    # plotting & tabling
    th <- max(branching.times(plottree))
    gr.col <- gray.colors(2, start=0.90, end=0.95)

#library(phylobase)
#library(phylosignal)
#tree4d<-phylo4d(plottree, cladesDR)
#
#pdf(file="DR_onMamPhy_BDvr_pcsFIXED_NDexp_5911taxa_MCCtarget_wDR_asbarsPhylo4d.pdf", width=25, height=25, onefile=TRUE)
#multiplot.phylo4d(p4d=tree4d, trait="harmMeans", show.tip=FALSE,tree.type = "fan", tree.open.angle = 45)
#dev.off()
missing<-read.table("MamPhy_FIN4_1813sp_missing_LIST.txt", header=FALSE)
sampled<-as.vector(setdiff(plottree$tip.label,missing$V1))
DR_sampled<-na.omit(as.data.frame(DR)[match(names(DR),sampled),])[1:length(sampled)]

tipColors <- rep("black",5911)
tipColors[match(missing$V1,plottree$tip.label)] <- grey(0.6) 


# NO TIPS... DR Dens in UPPER CENTER -- BAMM shifts in LOW LEFT
# or WITH tips...
###
#pdf(file="DR_onMamPhy_BDvr_pcsFIXED_NDexp_5911taxa_MCCtarget_w24shifts_UpperDRdens_higherOrdNumCladeNum_richCol.pdf", width=25, height=25, onefile=TRUE)

pdf(file="DR_onMamPhy_BDvr_pcsFIXED_NDexp_5911taxa_MCCtarget_w24shifts_UpperDRdens_higherOrdNumCladeNum_richCol_wTips.pdf", width=25, height=25, onefile=TRUE)
#pdf(file="DR_onMamPhy_BDvr_pcsFIXED_NDexp_5911taxa_MCCtarget_w24shifts_UpperDRdens_ordClade_viridis.pdf", width=25, height=25, onefile=TRUE)
#pdf(file="DR_onMamPhy_BDvr_pcsFIXED_NDexp_5911taxa_MCCtarget_w24shifts_UpperDRdens_ordClade_richCol.pdf", width=25, height=25, onefile=TRUE)

XX1=-220
XX2=220
##FBD
#XX1=-300
#XX2=300

plot(plottree, show.tip.label=FALSE, type="f", edge.width=0.9, cex=0.04, open.angle=45,no.margin=FALSE, root.edge=TRUE, edge.color="white", plot=FALSE, x.lim=c(XX1,XX2), y.lim=c(XX1,XX2))

draw.circle(0,0, radius=th, col =gr.col[2], border=gr.col[2])
draw.circle(0,0, radius=th-23.03, col=gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-66, col =gr.col[2], border=gr.col[2])
draw.circle(0,0, radius=th-66, col =NA, lty=2, lwd=3, border="black")
#draw.circle(0,0, radius=th-145, col =gr.col[1], border=gr.col[1])
lines<-rev(seq(0,170,10))
for (i in 1:length(lines)){
	draw.circle(0,0, radius=lines[i],col=NA,lty=2, lwd=2, border=grey(0.8,alpha=0.6))
}

legend(x=-60,y=-15, legend = c("up, 4x", "down, 2x", "up/down, ~1x"), title="Rate shifts", title.adj=0, y.intersp=1.5, bty = "n", box.lwd=2, box.col="black",lwd=2, cex=1.6, pt.cex = c(8,4,2), lty=c(NA,NA,NA), pt.bg = c(up, down, upDown), pch = c(21, 21, 21))

par(new=T)
#plot(plottree, show.tip.label=FALSE, type="f", edge.width=0.8, cex=0.05, open.angle=45,no.margin=FALSE, root.edge=TRUE, edge.color=as.matrix(reconColors), plot=TRUE, x.lim=c(XX1,XX2), y.lim=c(XX1,XX2))
obj <- plot2.phylo(plottree, show.tip.label=TRUE, tip.color=tipColors, cex=0.05, label.offset=0.08, type="f", edge.width=0.9, no.margin=FALSE, root.edge=TRUE, edge.color=as.matrix(reconColors),x.lim=c(XX1,XX2), y.lim=c(XX1,XX2))

nodelabels(pch=21, node=toPlotFIN$node, frame = "circle", bg = bgcols, cex= 2*bgsizes)
nodelabels(text=toPlotFIN$ID, node=toPlotFIN$node, adj=c(1.5,1.5), font=2, frame = "none", col = "black", cex= 2)
#nodelabels(text=toPlotFIN$CLADE, node=toPlotFIN$node, adj=c(1.5,2.5), font=1, frame = "none", col = "black", cex= 1.5)

par(new=T)

#group.label.tip.rad3(obj, lab=famNumAll, c( "light grey", grey(0.7) ), "black", offset.bar=1, offset.lab=2, cex=1, lwd=4, arc.bar.width=1.05)
group.label.tip.rad4(obj, lab=higherTax, c( "black",grey(0.4) ), "white", offset.bar=11, offset.lab=15, cex=1.8,font=1, lwd=4, arc.bar.width=1.04)
group.label.tip.rad4(obj, lab=cladeOrdNums, c( "black", grey(0.4) ), "black", offset.bar=27.5, offset.lab=38, cex=1.8, font=1, lwd=4, arc.bar.width=1.034)
group.label.tip.rad3(obj, lab=ordNumsSelect, c( "light grey", grey(0.6) ), "black", offset.bar=19.3, offset.lab=7, cex=1.8, font=1,lwd=4, arc.bar.width=1.037)


#group.label.tip.rad4(obj, lab=ords1, c( "black",grey(0.4) ), "white", offset.bar=2, offset.lab=6, cex=1.2,font=1, lwd=4, arc.bar.width=1.04)
#group.label.tip.rad3(obj, lab=ords2, c( NA,NA ), "black", offset.bar=2, offset.lab=6, cex=1,font=1, lwd=4, arc.bar.width=1.04)
#group.label.tip.rad3(obj, lab=clades, c( "light grey", grey(0.7) ), "black", offset.bar=12, offset.lab=18, cex=1, lwd=4, arc.bar.width=1.02)

#group.label.tip.rad3(obj, lab=fams, c( "light grey", grey(0.7) ), "black", offset.bar=12, offset.lab=18, cex=1, lwd=4, arc.bar.width=1.02)
#group.label.tip.rad3(obj, lab=famNums, c( "white", "black" ), "black", offset.bar=12, offset.lab=18, cex=1, lwd=4, arc.bar.width=1.02)

# for ND normal ##
#cols2<-palette(rich.colors(100))
color.legend(-46.3, 34.5, 32.7, 29.5, legend=NULL, rect.col= cols, gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey") #xl, yb, xr, yt
#color.legend(-75.1, -26, -30.4, -32, legend=NULL, rect.col= cols, gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey") #xl, yb, xr, yt
text(x=-4, y=15.5, "Tip diversification rate (species/Ma)", font=1, cex=1.8) 

par(fig=c(0.40, 0.655, 0.57, 0.67),mar=c(0,0,0,0), mar=c(0,0,0,0), new=T) #c(x1, x2, y1, y2)

plot(density(DR), col="dark grey", main="", bty="n", xlab="", ylab="",axes=F, xlim=range(DR))
polygon(density(DR), col="light grey", border="black", bty="n",main="")
x.tick <- quantile(DR, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick), labels=c(NA,round(x.tick,2)), side=1, line=1.3, cex=2.5, lwd=1, tck=-0.05, cex.axis=1.5, mgp=c(1,1,0))
dens.rate <- density(DR)$y
axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=2.5, las=1, lwd=1, cex.axis=1.5, tck=-0.05, mgp=c(1,1,0))
seg.tick <- quantile(DR, c(0.01,0.5,0.99))
segments(seg.tick[[1]],0,seg.tick[[1]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")
segments(seg.tick[[2]],0,seg.tick[[2]],max(dens.rate)*1, lty=2, lwd=2,col="black")
segments(seg.tick[[3]],0,seg.tick[[3]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")

#par(new=T)
#plot(density(DR_sampled), col="light grey", main="", bty="n", axes=F, xlim=range(DR))
#polygon(density(DR_sampled), col="light grey", border="black", bty="n")

dev.off()

#########
#####
# SAME plotting strategy but with the PRUNED mammal tree and DR results based on that...
###

## THE SAME, but for the pruned data -- 4098 species
####

# read in the tree on which data will be plotted
missing<-read.table("MamPhy_FIN4_1813sp_missing_LIST.txt", header=FALSE)

# plotting tree
bbone <- "NDexp" #"FBD"
mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")
plottree <- ladderize(mamMCC)

plottree <- drop.tip(plottree,as.character(missing$V1))

# read in the data for colouring branches
data <- read.table("DR-SUMMARY_MamPhy_BDvarRates_17Exp_all10ktrees_pruned4098sp.txt", row.names=1, header=TRUE)
data <- data[2:4099,]


>> Ahhhh, but these values were calculated on the BD-varRates trees-- 
i.e., with WRONG topologies-- need to re-calc on the FIXED topologies you did based on the global ML tree.



grp <- read.table("MamPhy_plotting_groups_prune4098.txt", header=TRUE)

fams <- grp[,2][match(plottree$tip.label,grp[,1])]
ords <- grp[,3][match(plottree$tip.label,grp[,1])]








source("DR_circleTree_functions.R")

# NO TIPS... DR Dens in LOW LEFT -- BAMM shifts in LOW RIGHT
###
pdf(file="DR_onMamPhy_BDvr_pcsFIXED_NDexp_5911taxa_MCCtarget_w35BAMMshifts_ordFam_difCol_10lines_LegR_DRdens.pdf", width=25, height=25, onefile=TRUE)

XX1=-220
XX2=220
##FBD
#XX1=-300
#XX2=300

plot(plottree, show.tip.label=FALSE, type="f", edge.width=0.8, cex=0.05, open.angle=45,no.margin=FALSE, root.edge=TRUE, edge.color=as.matrix(reconColors), plot=FALSE, x.lim=c(XX1,XX2), y.lim=c(XX1,XX2))
draw.circle(0,0, radius=th, col =gr.col[2], border=gr.col[2])
draw.circle(0,0, radius=th-23.03, col=gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-66, col =gr.col[2], border=gr.col[2])
draw.circle(0,0, radius=th-66, col =NA, lty=2, lwd=3, border="black")
#draw.circle(0,0, radius=th-145, col =gr.col[1], border=gr.col[1])
lines<-rev(seq(0,170,10))
for (i in 1:length(lines)){
	draw.circle(0,0, radius=lines[i],col=NA,lty=2, lwd=2, border=grey(0.8,alpha=0.6))
}

legend(x=20,y=-20, legend = c("up, 4x", "down, 2x", "up/down, ~1x"), title="BAMM rate shifts", title.adj=0, y.intersp=1.5, bty = "n", box.lwd=2, box.col="black",lwd=2, cex=1.4, pt.cex = c(8,4,2), lty=c(NA,NA,NA), pt.bg = c(up, down, upDown), pch = c(21, 21, 21))

par(new=T)
#plot(plottree, show.tip.label=FALSE, type="f", edge.width=0.8, cex=0.05, open.angle=45,no.margin=FALSE, root.edge=TRUE, edge.color=as.matrix(reconColors), plot=TRUE, x.lim=c(XX1,XX2), y.lim=c(XX1,XX2))
obj <- plot2.phylo(plottree, show.tip.label=FALSE, cex=0.05, label.offset=0.05, type="f", edge.width=0.8, no.margin=FALSE, root.edge=TRUE, edge.color=as.matrix(reconColors),x.lim=c(XX1,XX2), y.lim=c(XX1,XX2))

nodelabels(pch=21, node=toPlotFIN$node, frame = "circle", bg = bgcols, cex= 2*bgsizes)
nodelabels(text=toPlotFIN$ID, node=toPlotFIN$node, adj=c(-0.5,-0.5), frame = "none", col = "black", cex= 1.4)

par(new=T)

group.label.tip.rad4(obj, lab=ords, c( "black",grey(0.4) ), "white", offset.bar=2, offset.lab=6, cex=1.2,font=1, lwd=4, arc.bar.width=1.04)
group.label.tip.rad3(obj, lab=famFirst, c( "light grey", grey(0.7) ), "black", offset.bar=12, offset.lab=18, cex=1, lwd=4, arc.bar.width=1.02)


# for ND normal ##
#cols2<-palette(rich.colors(100))
#color.legend(-34.3, 55.5, 57.7, 50.5, legend=NULL, rect.col= cols2, gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey") #xl, yb, xr, yt
color.legend(-75.1, -26, -30.4, -32, legend=NULL, rect.col= cols, gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey") #xl, yb, xr, yt
text(x=-42, y=-45.5, "Diversification rate (species/Ma)", font=1, cex=1.5) 

par(fig=c(0.35, 0.495, 0.45, 0.50),mar=c(0,0,0,0), mar=c(0,0,0,0), new=T) #c(x1, x2, y1, y2)

plot(density(DR), col="dark grey", main="", bty="n", axes=F, xlim=range(DR))
polygon(density(DR), col="dark grey", border="black", bty="n")
x.tick <- quantile(DR, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick), labels=c(NA,round(x.tick,2)), side=1, line=1.3, cex=2.5, lwd=1, tck=-0.05, cex.axis=1.2, mgp=c(1,1,0))
dens.rate <- density(DR)$y
axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=2.5, las=1, lwd=1, cex.axis=1.2, tck=-0.05, mgp=c(1,1,0))
par(new=T)
plot(density(DR_sampled), col="light grey", main="", bty="n", axes=F, xlim=range(DR))
polygon(density(DR_sampled), col="light grey", border="black", bty="n")

dev.off()




###
# NO TIPS -- LINEAR with DR

pdf(file="DR_onMamPhy_BDvr_pcsFIXED_NDexp_5911taxa_MCCtarget_linear.pdf", width=2, height=6, onefile=TRUE)
obj <- plot2.phylo(plottree, show.tip.label=FALSE, cex=0.05, label.offset=0.05, type="phylogram", edge.width=0.3, no.margin=TRUE, root.edge=TRUE, edge.color=as.matrix(reconColors))
#nodelabels(pch=21, node=toPlotFIN$node, frame = "circle", bg = bgcols, cex= 2*bgsizes)
#nodelabels(text=toPlotFIN$ID, node=toPlotFIN$node, adj=c(-0.5,-0.5), frame = "none", col = "black", cex= 1.4)
#color.legend(-20,1,-30,15, legend=NULL, rect.col= cols, gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey") #xl, yb, xr, yt

#par(new=TRUE)
par(new=TRUE, fig=c(0.1, 0.6,0,0.5),mar=c(0,0,0,0), mar=c(0,0,0,0)) #c(x1, x2, y1, y2)

plot(x=seq(0.0,0.8,length.out=100),y=rep(0,100), col = rich.colors(100)[1], pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")
points(x=seq(0.4,0.8,length.out=100),y=rep(0,100), col = rich.colors(100)[100], pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")
points(x=seq(0.04,0.55,length.out=100),y=rep(0,100), col = rich.colors(100), pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")
x.tick <- quantile(DR, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick), labels=c(NA,round(x.tick,2)), side=1, line=-7, cex=2.5, lwd=1, tck=-0.05, cex.axis=0.4, mgp=c(1,0,0))
text(x=0.4, y=-0.25, "Diversification rate (species/Ma)", font=1, cex=0.4) 

dev.off()


###
# Now with TIPLABELS
###

missing<-read.table("MamPhy_FIN4_1813sp_missing_LIST.txt", header=FALSE)

sampled<-setdiff(plottree$tip.label,missing$V1)

tipColors <- rep("black",5911)
tipColors[match(missing$V1,plottree$tip.label)] <- "red" 

pdf(file="DR_onMamPhy_BDvarRates_Exp_5910taxa_all10k_MCC-tH_harmMean_colTips.pdf", width=25, height=25, onefile=TRUE)
pdf(file="DR_onMamPhy_BDvarRates_Exp_5910taxa_all10k_MCC-tH_CV_colTips.pdf", width=25, height=25, onefile=TRUE)

pdf(file="DR_onMamPhy_BDvr_pcsFIXED_FBD_5911taxa_sample1_harmMean_colTips.pdf", width=25, height=25, onefile=TRUE)

pdf(file="DR_onMamPhy_BDvr_pcsFIXED_NDexp_5911taxa_sample1_harmMean_colTips.pdf", width=25, height=25, onefile=TRUE)

XX1=-220
XX2=220
##FBD
#XX1=-300
#XX2=300

plot(plottree, show.tip.label=FALSE, type="f", edge.width=1, no.margin=TRUE, root.edge=TRUE, edge.color="white", plot=FALSE, x.lim=c(XX1,XX2), y.lim=c(XX1,XX2))
draw.circle(0,0, radius=th, col =gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-23.03, col=gr.col[2], border=gr.col[2])
draw.circle(0,0, radius=th-66, col =gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-145, col =gr.col[2], border=gr.col[2])

par(new=T)
obj <- plot2.phylo(plottree, show.tip.label=TRUE, tip.color=tipColors, cex=0.05, label.offset=0.05, type="f", edge.width=0.65, no.margin=TRUE, root.edge=TRUE, edge.color=as.matrix(reconColors),x.lim=c(XX1,XX2), y.lim=c(XX1,XX2))

par(new=T)

# for ND normal ##
color.legend(-27.3, -56.5, 57.7, -51.5, legend=NULL, rect.col= palette(rich.colors(100)), gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey")
text(x=26, y=-71, "Diversification rate (species per Myr)", cex=2) 
par(new=T)
group.label.tip.rad3(obj, lab=ords, c( "black","light grey"), "black", offset.bar=12, offset.lab=19, cex=1, lwd=4, arc.bar.width=1.02)
par(fig=c(0.42, 0.72, 0.39, 0.49),mar=c(0,0,0,0), mar=c(0,0,0,0), new=T)

##for FBD
#color.legend(-39.3, -76.5, 72.7, -70, legend=NULL, rect.col= palette(rich.colors(100)), gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey")
#text(x=26, y=-98, "Diversification rate (species per Myr)", cex=2) #"ED (Millions of Years)", cex=1.5)
##text(x=26, y=-65, "Coefficient of variation in 10,000 trees", cex=2) #"ED (Millions of Years)", cex=1.5)
#par(new=T)
#group.label.tip.rad3(obj, lab=ords, c( "black","light grey"), "black", offset.bar=12, offset.lab=19, cex=1, lwd=4, arc.bar.width=1.02)
#par(fig=c(0.42, 0.72, 0.39, 0.49),mar=c(0,0,0,0), mar=c(0,0,0,0), new=T)

plot(density(DR_median), col="dark grey", main="", bty="n", axes=F, xlim=range(DR_median))
polygon(density(DR_median), col="dark grey", border="dark grey", bty="n")

x.tick <- quantile(DR_median, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick), labels=c(0,round(x.tick,2)), side=1, line=1.3, cex=2.5, lwd=1, tck=-0.05, cex.axis=1.5, mgp=c(1,1,0))

dens.rate <- density(DR_median)$y
axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=2.5, las=1, lwd=1, cex.axis=1.5, tck=-0.05, mgp=c(1,1,0))

dev.off()






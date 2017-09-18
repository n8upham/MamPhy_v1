# Code to repeat the DR analyses on OTHER MAMMAL TREES...
#####
library(ape); library(picante); library(phytools); library(geiger)
setwd("/mnt/data/personal/nateu/Nate_Backup/DR_calcOnOther_mamPhys/")
source("DR_functions.R")

#=============
# ES and DR on per-tip basis
#######

# Do the HEDGES DR CALCS
##
# Load:: Hedges et al 2015 >> is ONE tree, no posterior !!!
##### #
#Hedges_phyUnsmooth<-read.tree(file="8.TTOL_mammals_unsmoothed.tre") # 1 tree-- 3738 tips and 3737 internal nodes
#Hedges_phySmooth<-read.tree(file="9.TTOL_mammals_smoothed_interpolated.tre") # 1 tree-- 5364 tips and 5363 internal nodes.
Hedges_phyUnsmooth_1<-scan("8.TTOL_mammals_unsmoothed.tre", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 
Hedges_phyUnsmooth<-strsplit(Hedges_phyUnsmooth_1,"[;]")  #trees now indexed #tree1=trees[[1]]  #first tree

Hedges_phySmooth<-scan("9.TTOL_mammals_smoothed_interpolated.tre", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 
Hedges_phySmooth<-strsplit(Hedges_phySmooth_1,"[;]")  #trees now indexed #tree1=trees[[1]]  #first tree

# gives pairwise clade matrix from CAIC function
clade_matrix = readCAIC(Hedges_phySmooth)

ES = ES_v2(clade_matrix)
DR = 1/ES
res = cbind.data.frame(DR,ES)
res1 = res[order(rownames(res)),]

write.table(res1, file="Hedges_phySmooth_TTOL_tipCalcs_1tree_DR-ES.txt")


#########
# Now looping through the calcs of 1k posteriors...
library(foreach);library(doSNOW)
cl2 = makeCluster(15, type = 'SOCK', outfile="")
registerDoSNOW(cl2)
source("DR_functions.R")


# LOAD:: Faurby and Svenning 2015-- 1k trees
##
##### #setwd("/Users/Nate/Desktop/JETZ-pdfs/Faurby_Svenning2015_SUPP/AppendixD_producedPhylogenies_Rcode")
#FS2015_phy1k<-read.tree(file="FS2015_Fully_resolved_phylogeny.nex") # NOT in nexus, is a NEWICK formated file 1000 trees, 5747 taxa each
#FS2015_phy1k_i<-scan("FS2015_Fully_resolved_phylogeny.nex", what="list",sep="\n",quiet=TRUE,n=1,skip=i+1,comment.char="#") 
FS2015_phy1k<-scan("FS2015_Fully_resolved_phylogeny.nex", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") # is newick format
#FS2015_phy1k<-strsplit(FS2015_phy1k_1,"[;]")  #trees now indexed #tree1=trees[[1]]  #first tree
trees<-FS2015_phy1k

# match all names with the followng rownames
cm = readCAIC(trees[[1]])

#run parallell
ntrees = length(trees)

g_es_dr= foreach(i=1:ntrees, .combine = cbind) %dopar% {   
   ES = ES_v2(readCAIC(trees[[i]]))
   DR = 1/ES
   res1 = cbind.data.frame(DR,ES)
   res2 = as.matrix(res1[match(cm$tip.label, rownames(res1)),]) 
   return(res2)
}

write.table(g_es_dr, file="FS2015_resolved_tipCalcs_1kTrees_DR-ES.txt")


# Load:: Kuhn et al. 2011 1k trees POSTERIOR of the Fritz et al. 2009 tree which is an update of Bininda-Emonds et al. 2007 tree
##
##### 
#Kuhn_phy1k_i<-scan("Fritz.Resolved.Normal.cr1r5_v2b.nwk", what="list",sep="\n",quiet=TRUE,n=1,skip=i+1,comment.char="#") 
Kuhn_phy1k<-scan("Fritz.Resolved.Normal.cr1r5_v2b.nwk", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 
#Kuhn_phy1k<-strsplit(Kuhn_phy1k_1,"[;]")  #trees now indexed #tree1=trees[[1]]  #first tree
trees<-Kuhn_phy1k

# match all names with the followng rownames
cm = readCAIC(trees[[1]])

#run parallell
ntrees = length(trees)

g_es_dr= foreach(i=1:ntrees, .combine = cbind) %dopar% {   
   ES = ES_v2(readCAIC(trees[[i]]))
   DR = 1/ES
   res1 = cbind.data.frame(DR,ES)
   res2 = as.matrix(res1[match(cm$tip.label, rownames(res1)),]) 
   return(res2)
}

write.table(g_es_dr, file="Kuhn2011_resolved_tipCalcs_1kTrees_DR-ES.txt")


#================
# SUMMARIZE the DR matrices into harmonic means of the posterior.
library(ape); library(picante); library(phytools); library(geiger)
setwd("/mnt/data/personal/nateu/Nate_Backup/DR_calcOnOther_mamPhys/")

# KUHN first.
Kuhn_DR1k<-read.table(file="Kuhn2011_resolved_tipCalcs_1kTrees_DR-ES.txt", header=TRUE)

DR_table<-vector("list",length=1001)
for (i in 0:1000){
	DR_table[[i+1]]<-Kuhn_DR1k[,i*2+1]
}
DR_tableAll<-do.call(cbind,DR_table)
rownames(DR_tableAll)<-rownames(Kuhn_DR1k)	

ES_table<-vector("list",length=1001)
for (i in 1:1001){
	ES_table[[i]]<-Kuhn_DR1k[,i*2]
}
ES_tableAll<-do.call(cbind,ES_table)
rownames(ES_tableAll)<-rownames(Kuhn_DR1k)	

#table<-DR_tableAll
table<-ES_tableAll

# Will want the mean, median, mode values across the ROWS for each tip in the tree >> 
harmMeans<-data.frame(matrix(NA, nrow = length(table[,1]), ncol = 1), row.names=rownames(table))
medians<-data.frame(matrix(NA, nrow = length(table[,1]), ncol = 1), row.names=rownames(table))
means<-data.frame(matrix(NA, nrow = length(table[,1]), ncol = 1), row.names=rownames(table))

for (i in 1:length(rownames(table))){
	harmMeans[i,] <- 1/(mean(1/(unlist(table[i,]))))
	medians[i,] <- median(unlist(table[i,]))
	means[i,] <- mean(unlist(table[i,]))
}
summary <- cbind(harmMeans,medians,means, deparse.level=1)
colnames(summary)<-c("harmMeans","medians","means")

#write.table(summary, "DR-SUMMARY_Kuhn2011_resolved_1kTrees_5020species.txt")
write.table(summary, "ES-SUMMARY_Kuhn2011_resolved_1kTrees_5020species.txt")


####
# Faurby and Svenning 2015 now.
##
FS2015_DR1k<-read.table(file="FS2015_resolved_tipCalcs_1kTrees_DR-ES.txt", header=TRUE)

DR_table<-vector("list",length=1000)
for (i in 0:999){
	DR_table[[i+1]]<-FS2015_DR1k[,i*2+1]
}
DR_tableAll<-do.call(cbind,DR_table)
rownames(DR_tableAll)<-rownames(FS2015_DR1k)	

ES_table<-vector("list",length=1000)
for (i in 1:1000){
	ES_table[[i]]<-FS2015_DR1k[,i*2]
}
ES_tableAll<-do.call(cbind,ES_table)
rownames(ES_tableAll)<-rownames(FS2015_DR1k)	

table<-DR_tableAll
#table<-ES_tableAll

# Will want the mean, median, mode values across the ROWS for each tip in the tree >> 
harmMeans<-data.frame(matrix(NA, nrow = length(table[,1]), ncol = 1), row.names=rownames(table))
medians<-data.frame(matrix(NA, nrow = length(table[,1]), ncol = 1), row.names=rownames(table))
means<-data.frame(matrix(NA, nrow = length(table[,1]), ncol = 1), row.names=rownames(table))

for (i in 1:length(rownames(table))){
	harmMeans[i,] <- 1/(mean(1/(unlist(table[i,]))))
	medians[i,] <- median(unlist(table[i,]))
	means[i,] <- mean(unlist(table[i,]))
}
summary <- cbind(harmMeans,medians,means, deparse.level=1)
colnames(summary)<-c("harmMeans","medians","means")

write.table(summary, "DR-SUMMARY_FS2015_resolved_1kTrees_5747species.txt")
#write.table(summary, "ES-SUMMARY_FS2015_resolved_1kTrees_5747species.txt")



# ==================
# now PLOT to SUMMARIZE
###
# Plot the PAIRWISE COMPARISONS too...
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/DR_calcsOnOther_mamPhys")
source("DR_circleTree_functions.R")
library(RColorBrewer)

bbone<-"NDexp"
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""))
colnames(cladesDR)<-c("tiplabel","gen","fam","famLabel","famNumLabel","famNumAll","ord","ordNums","ordLabel1","ordLabel2","clade","cladeCommonAll","cladeCommonSubs","cladeCommonSubs2","cladeCommonSubs2Nums","cladeCombo","higher","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")
head(cladesDR)

compareTable<-read.table("MamPhy_5911sp_clades_DR-COMPARE_vs_Kuhn-FS-Hedges.txt")
colnames(compareTable)<-c("tiplabel","ord","clade","higher","mamPhy","Kuhn2011","FS2015","Hedges2015")

# how many sp common across all trees?
noNAs<-na.omit(compareTable)
	# > 4403 species

length(na.omit(compareTable$Kuhn2011)) # 4670
length(na.omit(compareTable$FS2015)) # 5335
length(na.omit(compareTable$Hedges2015)) # 5033

# set up point colors
# by higherTaxon:
##
library(viridis)
cats<-names(sort(table(compareTable$higher),decreasing=FALSE))

#cols<-viridis(length(cats), alpha=0.7)
colsOrder<-rich.colors(length(cats), alpha=1)
cols<-c(colsOrder[6],colsOrder[4],colsOrder[2],colsOrder[5],colsOrder[3],colsOrder[1])

pointCols<-rep("black",length(noNAs[,1]))

for(j in 1:length(cats)){
pointCols[which(as.vector(noNAs$higher)==cats[j])]<-cols[j]
}

# by DR scale:
##
#cols<-rich.colors(100)
cols<-viridis(100)

pointCols<-rep("black",length(noNAs[,1]))

	range <- quantile(noNAs$mamPhy, seq(0,1, 0.01))[2:101]

    for (i in 1:100) {
        if (i==100) {range[i] <- range[i]+0.1}
        if (i==1) {pointCols[which(noNAs$mamPhy <= range[i])] <- cols[i] } 
        if (i > 1) {pointCols[which((noNAs$mamPhy <= range[i]) & (noNAs$mamPhy > range[i-1]))] <- cols[i]}
        }


# PLOT
#pdf(file="DR_pairwise_COMPARE_vs_Kuhn-FS-Hedges_richCol.pdf", height=3, width=8, onefile=TRUE)
#pdf(file="DR_pairwise_COMPARE_vs_Kuhn-FS-Hedges_viridis-by-higherTax.pdf", height=3, width=8, onefile=TRUE)
pdf(file="DR_pairwise_COMPARE_vs_Kuhn-FS-Hedges_richCol-by-higherTax_relative.pdf", height=3, width=10, onefile=TRUE)

#quartz(height=3,width=10)
par(oma = rep(3,4) + 0.1,
    mar = rep(1,4) + 0.1)

layout(matrix(1:4,1,4, byrow=TRUE))#,widths=rep(0.25,4),heights=rep(,4))
#layout(matrix(1:4,1,4, byrow=FALSE),widths=rep(0.25,4),heights=rep(0.75,4))
dat<-noNAs
maxDR<-max(dat[5:8])
	
	plot(mamPhy ~ Kuhn2011, data=dat, xlim=c(0,xMax),ylim=c(0,max(dat[5])), type="n", xlab="",ylab="",col=pointCols, bty="n", xaxt="n", yaxt="n", cex.lab=1.2,cex.axis=1.3)
	legend(x=0.1,y=0.6,legend=c("Monotremata", "Xenarthra", "Afrotheria","Marsupialia","Laurasiatheria", "Euarchontoglires"),fill=cols)

	xMax<-max(dat[6])
	plot(mamPhy ~ Kuhn2011, data=dat, xlim=c(0,xMax),ylim=c(0,max(dat[5])), xlab="",ylab="This study",col=pointCols, cex.lab=1.2,cex.axis=1.3)
	#corr<-cor(y=dat$mamPhy, x=dat$Kuhn2011, method="spearman")
	model<-summary(lm(mamPhy ~ Kuhn2011, data=dat))
	a<-round(model$coef[[1]], 2)
	b<-round(model$coef[[2]], 2)
	r2<-round(model$r.squared,digits=2)
	X=seq(min(dat$Kuhn2011),max(dat$Kuhn2011),length.out=20)
	lines(x=X,y=b*X+a,col="red", lwd=3,lty=1) #x1=min(dat$MRCA),x2=max(dat$MRCA),y1=min(log(dat$richness)),y2=max(log(dat$richness)),untf=FALSE,

	abline(a=0, b=max(dat[5])/xMax, lty=2, lwd=2, col="black")
	#text(labels=bquote("y = "*.(a)*" + "*.(b)*"x; "*r^2*" = "*.(r2)), x=xMax-(xMax/2.5), y=0.02, col="black", cex=1.1)
	mtext(text=bquote("y = "*.(a)*" + "*.(b)*"x; "*r^2*" = "*.(r2)), side=3, col="black", cex=0.8)
	mtext(text="Kuhn et al. 2011",side=1,line=3, font=2)
	mtext(text="This study",side=2,line=3, font=2)

	xMax<-max(dat[7])
	plot(mamPhy ~ FS2015, data=dat, xlim=c(0,max(dat[7])),ylim=c(0,max(dat[5])),  xlab="",ylab="",col=pointCols, yaxt="n", cex.axis=1.3)
	axis(side=2,labels=FALSE)
	mtext(text="Faurby and Svenning 2015",side=1,line=3, font=2)

	model<-summary(lm(mamPhy ~ FS2015, data=dat))
	a<-round(model$coef[[1]], 2)
	b<-round(model$coef[[2]], 2)
	r2<-round(model$r.squared,digits=2)
	X=seq(min(dat$FS2015),max(dat$FS2015),length.out=20)
	lines(x=X,y=b*X+a,col="red", lwd=3,lty=1) #x1=min(dat$MRCA),x2=max(dat$MRCA),y1=min(log(dat$richness)),y2=max(log(dat$richness)),untf=FALSE,

	abline(a=0, b=max(dat[5])/xMax, lty=2, lwd=2, col="black")
	#text(labels=bquote("y = "*.(a)*" + "*.(b)*"x; "*r^2*" = "*.(r2)), x=xMax-(xMax/2.5), y=0.02, col="black", cex=1.1)
	mtext(text=bquote("y = "*.(a)*" + "*.(b)*"x; "*r^2*" = "*.(r2)), side=3, col="black", cex=0.8)

	xMax<-max(dat[8])
	plot(mamPhy ~ Hedges2015, data=dat, xlim=c(0,max(dat[8])),ylim=c(0,max(dat[5])),  xlab="",ylab="",col=pointCols, yaxt="n", cex.axis=1.3)
	axis(side=2,labels=FALSE)
	mtext(text="Hedges et al. 2015",side=1,line=3, font=2)

	model<-summary(lm(mamPhy ~ Hedges2015, data=dat))
	a<-round(model$coef[[1]], 2)
	b<-round(model$coef[[2]], 2)
	r2<-round(model$r.squared,digits=2)
	X=seq(min(dat$Hedges2015),max(dat$Hedges2015),length.out=20)
	lines(x=X,y=b*X+a,col="red", lwd=3,lty=1) #x1=min(dat$MRCA),x2=max(dat$MRCA),y1=min(log(dat$richness)),y2=max(log(dat$richness)),untf=FALSE,

	abline(a=0, b=max(dat[5])/xMax, lty=2, lwd=2, col="black")
	#text(labels=bquote("y = "*.(a)*" + "*.(b)*"x; "*r^2*" = "*.(r2)), x=xMax-(xMax/2.5), y=0.02, col="black", cex=1.1)
	mtext(text=bquote("y = "*.(a)*" + "*.(b)*"x; "*r^2*" = "*.(r2)), side=3, col="black", cex=0.8)


dev.off()






# NO TIPS -- LINEAR with DR
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/DR_calcsOnOther_mamPhys")
source("DR_circleTree_functions.R")
library(ape); library(phytools)

##
# Trees TOGETHER in a loop, same timescale and color scale.

# Load trees
mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")
plottree1 <- ladderize(mamMCC)
r1<-max(nodeHeights(plottree1))

Kuhn_phy<-read.tree(file="Fritz.Resolved.Normal.cr1r5_v2b_firstTree.nwk")
plottree2<-ladderize(Kuhn_phy)
r2<-max(nodeHeights(plottree2))

FS2015_phy<-read.tree(file="FS2015_Fully_resolved_phylogeny_firstTree.nwk")
plottree3<-ladderize(FS2015_phy)
r3<-max(nodeHeights(plottree3))

Hedges_phySmooth<-read.tree(file="9.TTOL_mammals_smoothed_interpolated.tre") # 1 tree-- 5364 tips and 5363 internal nodes.
plottree4<-ladderize(Hedges_phySmooth)
r4<-max(nodeHeights(plottree3))

allPlotPhys<-list(plottree1,plottree2,plottree3,plottree4)
roots<-c(r1,r2,r3,r4)


# Load DRs
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""))
colnames(cladesDR)<-c("tiplabel","gen","fam","famLabel","famNumLabel","famNumAll","ord","ordNums","ordLabel1","ordLabel2","clade","cladeCommonAll","cladeCommonSubs","cladeCommonSubs2","cladeCommonSubs2Nums","cladeCombo","higher","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")
#head(cladesDR)
DR_mamPhy <- cladesDR[,"harmMeans"]
names(DR_mamPhy) <- cladesDR$tiplabel

DR_table3<-read.table("DR-SUMMARY_Kuhn2011_resolved_1kTrees_5020species.txt")
DR2 <- DR_table3$harmMeans
names(DR2) <- rownames(DR_table3)

DR_table4<-read.table("DR-SUMMARY_FS2015_resolved_1kTrees_5747species.txt")
DR3 <- DR_table4$harmMeans
names(DR3) <- rownames(DR_table4)

DR_ES<-read.table("Hedges_phySmooth_TTOL_tipCalcs_1tree_DR-ES.txt")
DR4<- DR_ES[,"DR"]
names(DR4) <- rownames(DR_ES)

allDRmeans<-list(DR_mamPhy,DR2,DR3,DR4)

Names<-c("This study","Kuhn et al. 2011","Faurby and Svenning 2015","Hedges et al. 2015")

# Make the color scale for ALL:
## basing this on the DR tip range, not the reconRate range, which is buffered by the BM ancestral recon...
	range <- quantile(DR_mamPhy, seq(0,1, 0.01))[2:101]
	range[100] <- range[100]+0.1
	#cols<-rev(colorRampPalette(c('red','orange','yellow','green','deepskyblue1','blue','black'))(100))
	#library(viridis)
	#cols<-viridis(100)
	cols<-rich.colors(100)

	x.tick <- quantile(DR_mamPhy, c(0.01,0.5,0.99,1))

# DO THE LOOP
##


#pdf(file="DR_on4phylosCompared_linear_richCol_justScale_ownColors.pdf", width=8, height=6, onefile=TRUE)
#pdf(file="DR_on4phylosCompared_linear_viridis_justScale_ownColors.pdf", width=8, height=6, onefile=TRUE)
pdf(file="DR_on4phylosCompared_linear_viridis_justScale_ownColors_withTips_80in.pdf", width=6, height=80, onefile=TRUE)

#quartz(width=8, height=6)
par(oma = rep(3,4) + 0.1, mar = rep(2,4) + 0.1)

#layout(matrix(1:2,2,1, byrow=FALSE),widths=rep(0.5,2),heights=rep(c(0.9,1),1))
#layout(matrix(1:8,2,4, byrow=FALSE),widths=rep(0.25,8),heights=rep(c(0.75,0.25),4))
#layout(matrix(1:4,1,4, byrow=FALSE),widths=rep(0.25,4),heights=rep(0.75,4))

for (j in 1:length(allPlotPhys)){
	plottree<-allPlotPhys[[j]]
	DR<-allDRmeans[[j]]

	DR_ordered <- DR[match(plottree$tip.label,names(DR))]
	DR_anc <- ace(DR_ordered, plottree, method="pic", scaled=TRUE)$ace

	# Match ancestors totree edges
    match.anc <- match(as.numeric(names(DR_anc)), plottree$edge[,2])[-1]

    # Assign rates to each internal node
    reconRate <- vector(mode="numeric", length=length(plottree$edge[,2]))
    reconRate[match.anc] <- DR_anc[2:(length(plottree$tip.label)-1)] #[2:7237]

    # Assign rates to tips
    tip.idx <- sort(plottree$tip.label, index.return=TRUE)$ix

    reconRate[match(tip.idx, plottree$edge[,2])] <- DR[sort(names(DR))]

    # Create colour palette
    reconColors <- reconRate

	range <- quantile(DR, seq(0,1, 0.01))[2:101]

    for (i in 1:100) {
        if (i==100) {range[i] <- range[i]+0.1}
        if (i==1) {reconColors[which(reconRate <= range[i])] <- cols[i] } 
        if (i > 1) {reconColors[which((reconRate <= range[i]) & (reconRate > range[i-1]))] <- cols[i]}
        }

    # Plot the tree 
	obj <- plot2.phylo(plottree, show.tip.label=TRUE, cex=0.05, tip.color="black",x.lim=c(0,roots[j]),label.offset=0.1, type="phylogram", edge.width=0.3, no.margin=TRUE, root.edge=TRUE, edge.color=as.matrix(reconColors))
	axisPhylo()#at=c(0,50,100,150,200),labels=TRUE)
	mtext(side=3,text=Names[j], cex=0.8, font=2)
	#par(new=TRUE, fig=c(0.0, 0.15,0.05,0.25),mar=rep(1,4)) #c(x1, x2, y1, y2)
	#par(fig=c(0.1, 0.9,0,0.9))#,mar=c(0,0,0,0), mar=c(0,0,0,0)) #c(x1, x2, y1, y2)
#	plot(density(DR), add=TRUE, col="dark grey", main="", bty="n", xlab="", ylab="",axes=F, xlim=range(DR))
#	polygon(density(DR), col="light grey", border="black", bty="n",main="")
#	dens.rate <- density(DR)$y
#	axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=1, las=1, lwd=1, cex.axis=1, tck=-0.05, mgp=c(1,1,0))

	x.tick_other <- quantile(DR, c(0.01,0.5,0.99,1))

## USE THIS #
####
#	Y=0
#	plot(x=seq(0.0,x.tick_other[[4]],length.out=100),y=rep(Y,100), col = cols[1], pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")
#	points(x=seq(x.tick_other[[3]],x.tick_other[[4]],length.out=100),y=rep(Y,100), col = cols[100], pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")
#	points(x=seq(x.tick_other[[1]],x.tick_other[[3]],length.out=100),y=rep(Y,100), col = cols, pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")
#
#	axis(at=c(0,x.tick_other), side=1, line=-4, labels=FALSE, cex=2.5, lwd=1, tck=-0.05, cex.axis=1, mgp=c(1,0,0))
#	text(x= c(0,x.tick_other), y = rep(-0.4,5), cex=0.8,font=2,srt = 45, labels = c(NA,round(x.tick_other,2)), xpd = TRUE, adj=c(0.95,0.05))


#	Y=0
#	points(x=seq(0.0,x.tick[[4]],length.out=100),y=rep(Y,100), col = cols[1], pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")
#	points(x=seq(x.tick[[3]],x.tick[[4]],length.out=100),y=rep(Y,100), col = cols[100], pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")
#	points(x=seq(x.tick[[1]],x.tick[[3]],length.out=100),y=rep(Y,100), col = cols, pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")
	
	#x.tick_other <- quantile(DR, c(0.01,0.5,0.99,1))
	#axis(at=c(0,x.tick_other), labels=FALSE, side=1, line=-7, cex=2.5, lwd=1, tck=-0.05, cex.axis=0.4, mgp=c(1,0,0))
	#text(x= c(0,x.tick_other)+0.03, y = rep(-0.15,5), cex=0.45,font=2,srt = 45, labels = c(NA,round(x.tick_other,2)), xpd = TRUE, adj=c(0.95,0.05))
	#text(x=0.4, y=-0.25, "Diversification rate (species/Ma)", font=1, cex=0.4) 

}

dev.off()


>> basically it appears the polygon part is not so interesting!

# compare Ctenomys
DR1<-allDRmeans[[1]]
DR1_gen<-cladesDR[which(names(DR1)==cladesDR$tiplabel),"gen"]
DR1_clade<-DR1[which(as.vector(DR1_gen)=="Ctenomys")]
DR1_clade_mean<-mean(DR1_clade)
DR1_clade_mean

compareTableTips<-cbind.data.frame(sort(cladesDR$tiplabel),cladesDR$gen[order(cladesDR$tiplabel)],compareTable)

colnames(compareTableTips)<-c("tiplabel","gen","sciName",colnames(compareTable)[2:8])

# Get Ctenomys SUBSET from the compare table
allDR_sub<-compareTableTips[which(compareTableTips$gen=="Ctenomys"),]

allDR_sub_means<-colMeans(allDR_sub[,7:10],na.rm=TRUE)

allDR_sub_means_inQuantile<-rep(NA,4)
up95s<-vector("list",length=4)
for(i in 1:4){
up95s[[i]]<-quantile(compareTableTips[,6+i],seq(0.01,0.99,by=0.01),na.rm=TRUE)
	res<-rep(NA,length(up95s[[i]]))
	names(res)<-names(up95s[[i]])
	for(j in 1:length(up95s[[i]])){
		if(up95s[[i]][[j]] <= allDR_sub_means[i] && up95s[[i]][[j+1]] >= allDR_sub_means[i]) {res[j]<-1} else {res[j]<-0}
	}
	allDR_sub_means_inQuantile[i]<-names(res[which(res==1)])
}
allRes<-rbind(round(allDR_sub_means,3),allDR_sub_means_inQuantile)
rownames(allRes)<-c("mean","inQuant")
#        mamPhy  Kuhn2011 FS2015  Hedges2015
#mean    "0.526" "0.147"  "0.282" "0.203"   
#inQuant "98%"   "65%"    "69%"   "80%"     



	seg.tick <- quantile(DR_mamPhy, c(0.01,0.5,0.99))
	segments(seg.tick[[1]],0,seg.tick[[1]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")
	segments(seg.tick[[2]],0,seg.tick[[2]],max(dens.rate)*1, lty=2, lwd=2,col="black")
	segments(seg.tick[[3]],0,seg.tick[[3]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")


	color.legend(0, 0, 0.5, 0.3, legend=NULL, rect.col= cols, gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey") #xl, yb, xr, yt
	#color.legend(-75.1, -26, -30.4, -32, legend=NULL, rect.col= cols, gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey") #xl, yb, xr, yt
	text(x=-4, y=15.5, "Tip diversification rate (species/Ma)", font=1, cex=1.8) 

	par(fig=c(0.40, 0.655, 0.57, 0.67),mar=c(0,0,0,0), mar=c(0,0,0,0), new=T) #c(x1, x2, y1, y2)



	x.tick <- quantile(DR, c(0.01,0.5,0.99,1))
	axis(at=c(0,x.tick), labels=c(NA,round(x.tick,2)), side=1, line=1.3, cex=2.5, lwd=1, tck=-0.05, cex.axis=1.5, mgp=c(1,1,0))









# ===================
##
# Each tree SEPARATELY

#=====
# on our MamPhy NDexp
#=====
mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")
plottree <- ladderize(mamMCC)

#max(nodeHeights(mamMCC))
#[1] 180.8472

cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""))
colnames(cladesDR)<-c("tiplabel","gen","fam","famLabel","famNumLabel","famNumAll","ord","ordNums","ordLabel1","ordLabel2","clade","cladeCommonAll","cladeCommonSubs","cladeCommonSubs2","cladeCommonSubs2Nums","cladeCombo","higher","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")
head(cladesDR)

# load DR
DR_mamPhy <- cladesDR[,"harmMeans"]
names(DR_mamPhy) <- cladesDR$tiplabel

DR_ordered <- DR_mamPhy[match(plottree$tip.label,names(DR_mamPhy))]
DR_anc <- ace(DR_ordered, plottree, method="pic", scaled=TRUE)$ace

# Match ancestors totree edges
    match.anc <- match(as.numeric(names(DR_anc)), plottree$edge[,2])[-1]

    # Assign rates to each internal node
    reconRate_mamPhy <- vector(mode="numeric", length=length(plottree$edge[,2]))
    reconRate_mamPhy[match.anc] <- DR_anc[2:(length(plottree$tip.label)-1)] #[2:7237]

    # Assign rates to tips
    tip.idx <- sort(plottree$tip.label, index.return=TRUE)$ix

    reconRate_mamPhy[match(tip.idx, plottree$edge[,2])] <- DR_mamPhy[sort(names(DR_mamPhy))]

    # Create colour palette
    reconColors_mamPhy <- reconRate_mamPhy

	#cols<-rev(colorRampPalette(c('red','orange','yellow','green','deepskyblue1','blue','black'))(100))
	
	#library(viridis)
	#cols<-viridis(100)
	cols<-rich.colors(100)
	range <- quantile(reconRate_mamPhy, seq(0,1, 0.01))[2:101]
    for (i in 1:100) {
        if (i==100) {range[i] <- range[i]+0.1}
        if (i==1) {reconColors_mamPhy[which(reconRate_mamPhy <= range[i])] <- cols[i] } 
        if (i > 1) {reconColors_mamPhy[which((reconRate_mamPhy <= range[i]) & (reconRate_mamPhy > range[i-1]))] <- cols[i]}
        }

pdf(file="DR_onMamPhy_BDvr_pcsFIXED_NDexp_5911taxa_MCCtarget_linear_richCol.pdf", width=2, height=6, onefile=TRUE)
obj <- plot2.phylo(plottree, show.tip.label=FALSE, cex=0.05, label.offset=0.05, type="phylogram", edge.width=0.3, no.margin=TRUE, root.edge=TRUE, edge.color=as.matrix(reconColors_mamPhy))
#nodelabels(pch=21, node=toPlotFIN$node, frame = "circle", bg = bgcols, cex= 2*bgsizes)
#nodelabels(text=toPlotFIN$ID, node=toPlotFIN$node, adj=c(-0.5,-0.5), frame = "none", col = "black", cex= 1.4)
#color.legend(-20,1,-30,15, legend=NULL, rect.col= cols, gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey") #xl, yb, xr, yt

#par(new=TRUE)
par(new=TRUE, fig=c(0.1, 0.6,0,0.5),mar=c(0,0,0,0), mar=c(0,0,0,0)) #c(x1, x2, y1, y2)

x.tick <- quantile(DR_mamPhy, c(0.01,0.5,0.99,1))
plot(x=seq(0.0,x.tick[[4]],length.out=100),y=rep(0,100), col = cols[1], pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")
points(x=seq(x.tick[[3]],x.tick[[4]],length.out=100),y=rep(0,100), col = cols[100], pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")
points(x=seq(x.tick[[1]],x.tick[[3]],length.out=100),y=rep(0,100), col = cols, pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")
axis(at=c(0,x.tick), labels=FALSE, side=1, line=-7, cex=2.5, lwd=1, tck=-0.05, cex.axis=0.4, mgp=c(1,0,0))
text(x= c(0,x.tick)+0.03, y = rep(-0.15,5), cex=0.45,font=2,srt = 45, labels = c(NA,round(x.tick,2)), xpd = TRUE, adj=c(0.95,0.05))

text(x=0.4, y=-0.25, "Diversification rate (species/Ma)", font=1, cex=0.4) 

dev.off()


#=====
# on Hedges et al. 2015
#=====
DR_ES<-read.table("Hedges_phySmooth_TTOL_tipCalcs_1tree_DR-ES.txt")
Hedges_phySmooth<-read.tree(file="9.TTOL_mammals_smoothed_interpolated.tre") # 1 tree-- 5364 tips and 5363 internal nodes.
plottree<-ladderize(Hedges_phySmooth)

#max(nodeHeights(Hedges_phySmooth))
#[1] 179.1419

# load DR
DR <- DR_ES[,"DR"]
names(DR) <- rownames(DR_ES)

DR_ordered <- DR[match(plottree$tip.label,names(DR))]
DR_anc <- ace(DR_ordered, plottree, method="pic", scaled=TRUE)$ace

# Match ancestors totree edges
    match.anc <- match(as.numeric(names(DR_anc)), plottree$edge[,2])[-1]

    # Assign rates to each internal node
    reconRate <- vector(mode="numeric", length=length(plottree$edge[,2]))
    reconRate[match.anc] <- DR_anc[2:(length(plottree$tip.label)-1)] #[2:7237]

    # Assign rates to tips
    tip.idx <- sort(plottree$tip.label, index.return=TRUE)$ix

    reconRate[match(tip.idx, plottree$edge[,2])] <- DR[sort(names(DR))]

    # USE THE SAME COLORS AS FOR MAMPHY

    #Create colour palette
    reconColors <- reconRate

    for (i in 1:100) {
        if (i==1) {reconColors[which(reconRate <= range[i])] <- cols[i] } 
        if (i > 1) {reconColors[which((reconRate <= range[i]) & (reconRate > range[i-1]))] <- cols[i]}
        }


#pdf(file="DR_onHedges2015_phySmooth_5364taxa_linear_viridis.pdf", width=2, height=6, onefile=TRUE)
pdf(file="DR_onHedges2015_phySmooth_5364taxa_linear_richCol.pdf", width=2, height=6, onefile=TRUE)
obj <- plot2.phylo(plottree, show.tip.label=FALSE, cex=0.05, label.offset=0.05, type="phylogram", edge.width=0.3, no.margin=TRUE, root.edge=TRUE, edge.color=as.matrix(reconColors))
#nodelabels(pch=21, node=toPlotFIN$node, frame = "circle", bg = bgcols, cex= 2*bgsizes)
#nodelabels(text=toPlotFIN$ID, node=toPlotFIN$node, adj=c(-0.5,-0.5), frame = "none", col = "black", cex= 1.4)
#color.legend(-20,1,-30,15, legend=NULL, rect.col= cols, gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey") #xl, yb, xr, yt

#par(new=TRUE)
par(new=TRUE, fig=c(0.1, 0.6,0,0.5),mar=c(0,0,0,0), mar=c(0,0,0,0)) #c(x1, x2, y1, y2)

plot(x=seq(0.0,x.tick[[4]],length.out=100),y=rep(0,100), col = cols[1], pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")
points(x=seq(x.tick[[3]],x.tick[[4]],length.out=100),y=rep(0,100), col = cols[100], pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")
points(x=seq(x.tick[[1]],x.tick[[3]],length.out=100),y=rep(0,100), col = cols, pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")

x.tick_other <- quantile(DR, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick_other), labels=FALSE, side=1, line=-7, cex=2.5, lwd=1, tck=-0.05, cex.axis=0.4, mgp=c(1,0,0))
text(x= c(0,x.tick_other)+0.03, y = rep(-0.15,5), cex=0.45,font=2,srt = 45, labels = c(NA,round(x.tick_other,2)), xpd = TRUE, adj=c(0.95,0.05))

text(x=0.4, y=-0.25, "Diversification rate (species/Ma)", font=1, cex=0.4) 

dev.off()



#=====
# on KUHN et al. 2011
#=====
DR_table<-read.table("DR-SUMMARY_Kuhn2011_resolved_1kTrees_5020species.txt")

#Kuhn_phy_1<-scan("Fritz.Resolved.Normal.cr1r5_v2b.nwk", what="list",sep="\n",n=1,quiet=TRUE,skip=0,comment.char="#") 
#tre<-read.tree(text=Kuhn_phy_1)
#write.tree(tre,file="Fritz.Resolved.Normal.cr1r5_v2b_firstTree.nwk")
Kuhn_phy<-read.tree(file="Fritz.Resolved.Normal.cr1r5_v2b_firstTree.nwk")

plottree<-ladderize(Kuhn_phy)

#max(nodeHeights(Kuhn_phy))
#[1] 166.18

# load DR
DR <- DR_table$harmMeans
names(DR) <- rownames(DR_table)

DR_ordered <- DR[match(plottree$tip.label,names(DR))]
DR_anc <- ace(DR_ordered, plottree, method="pic", scaled=TRUE)$ace

# Match ancestors totree edges
    match.anc <- match(as.numeric(names(DR_anc)), plottree$edge[,2])[-1]

    # Assign rates to each internal node
    reconRate <- vector(mode="numeric", length=length(plottree$edge[,2]))
    reconRate[match.anc] <- DR_anc[2:(length(plottree$tip.label)-1)] #[2:7237]

    # Assign rates to tips
    tip.idx <- sort(plottree$tip.label, index.return=TRUE)$ix

    reconRate[match(tip.idx, plottree$edge[,2])] <- DR[sort(names(DR))]

    # USE THE SAME COLORS AS FOR MAMPHY

    #Create colour palette
    reconColors <- reconRate

    for (i in 1:100) {
        if (i==1) {reconColors[which(reconRate <= range[i])] <- cols[i] } 
        if (i > 1) {reconColors[which((reconRate <= range[i]) & (reconRate > range[i-1]))] <- cols[i]}
        }


pdf(file="DR_onKuhn2011_5020species_linear_richCol.pdf", width=2, height=6, onefile=TRUE)
obj <- plot2.phylo(plottree, show.tip.label=FALSE, cex=0.05, label.offset=0.05, type="phylogram", edge.width=0.3, no.margin=TRUE, root.edge=TRUE, edge.color=as.matrix(reconColors))
#nodelabels(pch=21, node=toPlotFIN$node, frame = "circle", bg = bgcols, cex= 2*bgsizes)
#nodelabels(text=toPlotFIN$ID, node=toPlotFIN$node, adj=c(-0.5,-0.5), frame = "none", col = "black", cex= 1.4)
#color.legend(-20,1,-30,15, legend=NULL, rect.col= cols, gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey") #xl, yb, xr, yt

#par(new=TRUE)
par(new=TRUE, fig=c(0.1, 0.6,0,0.5),mar=c(0,0,0,0), mar=c(0,0,0,0)) #c(x1, x2, y1, y2)

plot(x=seq(0.0,x.tick[[4]],length.out=100),y=rep(0,100), col = cols[1], pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")
points(x=seq(x.tick[[3]],x.tick[[4]],length.out=100),y=rep(0,100), col = cols[100], pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")
points(x=seq(x.tick[[1]],x.tick[[3]],length.out=100),y=rep(0,100), col = cols, pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")

x.tick_other <- quantile(DR, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick_other), labels=FALSE, side=1, line=-7, cex=2.5, lwd=1, tck=-0.05, cex.axis=0.4, mgp=c(1,0,0))
text(x= c(0,x.tick_other)+0.03, y = rep(-0.15,5), cex=0.45,font=2,srt = 45, labels = c(NA,round(x.tick_other,2)), xpd = TRUE, adj=c(0.95,0.05))

text(x=0.4, y=-0.25, "Diversification rate (species/Ma)", font=1, cex=0.4) 

dev.off()


#=====
# on FS2015
#=====
DR_table<-read.table("DR-SUMMARY_FS2015_resolved_1kTrees_5747species.txt")

#FS2015_phy_1<-scan("FS2015_Fully_resolved_phylogeny.nex", what="list",sep="\n",n=1,quiet=TRUE,skip=0,comment.char="#") 
#tre<-read.tree(text=FS2015_phy_1)
#write.tree(tre,file="FS2015_Fully_resolved_phylogeny_firstTree.nwk")
FS2015_phy<-read.tree(file="FS2015_Fully_resolved_phylogeny_firstTree.nwk")

plottree<-ladderize(FS2015_phy)

#max(nodeHeights(FS2015_phy))
#[1] 217.885

# load DR
DR <- DR_table$harmMeans
names(DR) <- rownames(DR_table)

DR_ordered <- DR[match(plottree$tip.label,names(DR))]
DR_anc <- ace(DR_ordered, plottree, method="pic", scaled=TRUE)$ace

# Match ancestors totree edges
    match.anc <- match(as.numeric(names(DR_anc)), plottree$edge[,2])[-1]

    # Assign rates to each internal node
    reconRate <- vector(mode="numeric", length=length(plottree$edge[,2]))
    reconRate[match.anc] <- DR_anc[2:(length(plottree$tip.label)-1)] #[2:7237]

    # Assign rates to tips
    tip.idx <- sort(plottree$tip.label, index.return=TRUE)$ix

    reconRate[match(tip.idx, plottree$edge[,2])] <- DR[sort(names(DR))]

    # USE THE SAME COLORS AS FOR MAMPHY

    #Create colour palette
    reconColors <- reconRate

    for (i in 1:100) {
        if (i==1) {reconColors[which(reconRate <= range[i])] <- cols[i] } 
        if (i > 1) {reconColors[which((reconRate <= range[i]) & (reconRate > range[i-1]))] <- cols[i]}
        }


pdf(file="DR_onFS2015_5747species_linear_richCol.pdf", width=2, height=6, onefile=TRUE)
obj <- plot2.phylo(plottree, show.tip.label=FALSE, cex=0.05, label.offset=0.05, type="phylogram", edge.width=0.3, no.margin=TRUE, root.edge=TRUE, edge.color=as.matrix(reconColors))
#nodelabels(pch=21, node=toPlotFIN$node, frame = "circle", bg = bgcols, cex= 2*bgsizes)
#nodelabels(text=toPlotFIN$ID, node=toPlotFIN$node, adj=c(-0.5,-0.5), frame = "none", col = "black", cex= 1.4)
#color.legend(-20,1,-30,15, legend=NULL, rect.col= cols, gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey") #xl, yb, xr, yt

#par(new=TRUE)
par(new=TRUE, fig=c(0.1, 0.6,0,0.5),mar=c(0,0,0,0), mar=c(0,0,0,0)) #c(x1, x2, y1, y2)

plot(x=seq(0.0,x.tick[[4]],length.out=100),y=rep(0,100), col = cols[1], pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")
points(x=seq(x.tick[[3]],x.tick[[4]],length.out=100),y=rep(0,100), col = cols[100], pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")
points(x=seq(x.tick[[1]],x.tick[[3]],length.out=100),y=rep(0,100), col = cols, pch = 16, cex = 1, bty="n", xlab="",ylab="",xaxt="n",yaxt="n")

x.tick_other <- quantile(DR, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick_other), labels=FALSE, side=1, line=-7, cex=2.5, lwd=1, tck=-0.05, cex.axis=0.4, mgp=c(1,0,0))
text(x= c(0,x.tick_other)+0.03, y = rep(-0.15,5), cex=0.45,font=2,srt = 45, labels = c(NA,round(x.tick_other,2)), xpd = TRUE, adj=c(0.95,0.05))

text(x=0.4, y=-0.25, "Diversification rate (species/Ma)", font=1, cex=0.4) 

dev.off()


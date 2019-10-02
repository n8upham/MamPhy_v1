#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy MS1 -- Upham et al. 2019 -- PLOS Biology
###
# Figure 4 - comparing divergence times among Mammalia backbone-level studies
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# DIVTIME comparison across BACKBONE-LEVEL TREES
####
setwd("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
library(ape); library(phytools); library(viridis); library(phangorn); library(geiger); library(magick); library(phyloch)
library(plotrix)

# load CSV of times
divTimes1<-read.csv("compare_mamPhy_cladeDivTimes_toPlot.csv")
#divTimes<-divTimes1#[as.numeric(rownames(na.omit(divTimes1[,1:8]))),]
#divTimes<-divTimes1[c(which(divTimes1$LEVEL < 4), c(29:32,38:39,50,51,53,57:59)),]
divTimes<-read.csv("compare_mamPhy_cladeDivTimes_toPlot_reducedNew.csv")

colBE2007<-plasma(10)[1]
colM2011<-plasma(10)[3]
colDR2012<-plasma(10)[5]
colR2016<-plasma(10)[7]

pdf(file="compare_mamPhy_cladeDivTimes_plot95CI_4studies_wFossils_OK.pdf",height=18,width=9,onefile=TRUE)

	plotCI(y=rev(1:length(divTimes$Est)), x= -divTimes$Est, ui= -divTimes$Low, li= -divTimes$Up, err="x",
	       sfrac=0, add=FALSE, scol="white", pch=20, yaxt="n",bty="n", ylab="",xlab="", xpd=NA, xlim=c(-140,0), lwd=2, cex=1.5)

		for(j in (2*(1:(length(divTimes$Est)/2))) ){
			#abline(h=j+0.1, lwd=0.3, col=grey(0.7,alpha=0.8))	
			#lines(x=c(-125,0), y=c(j+0.2,j+0.2), lwd=0.3, col=grey(0.85,alpha=1))	
			rect(angle=90, xleft= -155, xright=0, ybottom=j-0.8, ytop=j+0.2, border=NA, col=grey(0.95,alpha=1), xpd=NA)	
		}
		text(labels=(divTimes$Taxon), x=(-155+(2*divTimes$LEVEL)), y=rev(1:length(divTimes$Est))-0.2, xpd=NA, adj=0, cex=1) #
		lines(x=c(-66,-66), y=c(-2,55), lwd=1, col=grey(0.4,alpha=1), lty=2)	

	plotCI(y=rev(1:length(divTimes$Est)), x= -divTimes$Est, ui= -divTimes$Low, li= -divTimes$Up, err="x",
	       sfrac=0, add=TRUE, scol="black", pch=20, yaxt="n",bty="n", ylab="",xlab="", xpd=NA, xlim=c(-130,0), lwd=2.3, cex=1.4)

	# BE2007
	plotCI(y=rev(1:length(divTimes$Est))-0.15, x= -divTimes$BE2007_Est, ui= -as.numeric(as.vector(divTimes$BE2007_Low)), li= -as.numeric(as.vector(divTimes$BE2007_Up)), err="x",
	       sfrac=0, add=TRUE, scol=colBE2007, col=colBE2007, pch=20, yaxt="n",bty="n", ylab="",xlab="", xpd=NA, lwd=1.2)

	# M2011
	plotCI(y=rev(1:length(divTimes$Est))-0.3, x= -divTimes$M2011_Est, ui= -as.numeric(as.vector(divTimes$M2011_Low)), li= -as.numeric(as.vector(divTimes$M2011_Up)), err="x",
	       sfrac=0, add=TRUE, scol=colM2011, col=colM2011, pch=20, yaxt="n",bty="n", ylab="",xlab="", xpd=NA, lwd=1.2)

	# DR2012
	plotCI(y=rev(1:length(divTimes$Est))-0.45, x= -divTimes$DR2012_Est, ui= -as.numeric(as.vector(divTimes$DR2012_Low)), li= -as.numeric(as.vector(divTimes$DR2012_Up)), err="x",
	       sfrac=0, add=TRUE, scol=colDR2012, col=colDR2012, pch=20, yaxt="n",bty="n", ylab="",xlab="", xpd=NA, lwd=1.2)

	# R2016
	plotCI(y=rev(1:length(divTimes$Est))-0.6, x= -divTimes$R2016_Est, ui= -as.numeric(as.vector(divTimes$R2016_Low)), li= -as.numeric(as.vector(divTimes$R2016_Up)), err="x",
	       sfrac=0, add=TRUE, scol=colR2016, col=colR2016, pch=20, yaxt="n",bty="n", ylab="",xlab="", xpd=NA, lwd=1.2)

	# MIN CROWN AGE to MAX
	#fossilDat<-na.omit(divTimes[,c("Foley2016_crownMIN","PBDB_maxStemAge")])
	middleDates<- (as.numeric(as.vector(divTimes$Foley2016_crownMAX))-as.numeric(as.vector(divTimes$Foley2016_crownMIN)))/2+ -(as.numeric(as.vector(divTimes$Foley2016_crownMAX)) )

	plotCI(y=rev(1:length(divTimes$Est))-0.75, x= middleDates, ui= -as.numeric(as.vector(divTimes$Foley2016_crownMIN)), li= -as.numeric(as.vector(divTimes$Foley2016_crownMAX)), err="x",
		col="goldenrod1", add=TRUE, sfrac=0.002, pch=NA, lwd=2)

	# PBDB max fossil AS STAR
	points(y=rev(1:length(divTimes$Est))-0.73, x= -as.numeric(as.vector(divTimes$PBDB_maxStemAge)),
		bg="goldenrod1", pch=22, cex=1.1, col="goldenrod1")
	#text(y=rev(1:length(divTimes$Est))-0.65, x= -as.numeric(as.vector(divTimes$PBDB_maxStemAge)),
	#		labels="t", cex=0.4, font=2, col="black")
#	points(y=rev(1:length(divTimes$Est))-0.65, x= -as.numeric(as.vector(divTimes$Foley2016_crownMIN)),
#		bg="chartreuse3", pch=22, cex=1.2, col="black")
#

dev.off()


# Get the LEGEND
####
pdf(file="legend_compare_mamPhy_cladeDivTimes_plot95CI.pdf")
	plot(1:10, ylab="",xlab="", xaxt="n",yaxt="n",pch=NA)
	legend(x=1, y=10, legend=c(NA, "This study", "Bininda-Emonds et al. 2007", "Meredith et al. 2011", "dos Reis et al. 2012", "Ronquist et al. 2016", "", "Crown min-to-max (Foley et al. 2016)", "Stem max (Paleobiology Database)"), 
		col=c(NA, "black", colBE2007, colM2011, colDR2012, colR2016, NA, "goldenrod1","goldenrod1"), bty="n",lty=1, lwd=c(NA, 2.3, rep(1.2,4), NA, 2, NA), pch=c(NA, rep(20,5),NA,NA,15), pt.cex=c(NA, 1.4,rep(1,4),NA,NA,1.1) )

	text(x= 1, y= 9.65, label="Fossil-calibrated molecular ages", font=2, cex=1, adj=0)
	text(x= 1, y= 7.3, label="Fossil ages directly", font=2, cex=1, adj=0)
dev.off()


legend(x= -100, y= 6500, legend=c(" ","Bininda-Emonds et al. 2007 (one tree)", "Kuhn et al. 2011 (set of trees)"," ","Faurby & Svenning 2015 (set of trees)", 
	" ", "Hedges et al. 2015 (one tree)", " ",
	"This study (set of trees)", "This study (MCC tree)"), 
	col=c(NA, colBE, colKuhn,NA, colFS,NA, colHedges,NA, "goldenrod1","black"), #col=c(NA, "purple", NA, "blue",NA,"red",NA,"goldenrod1",NA,grey(0.5),"black"), 
	lwd=2, lty=1, cex=0.7)
text(x= -98, y= 4900, label="MRP supertree", font=2, cex=0.7, adj=0)
text(x= -98, y= 1850, label="DNA supertree", font=2, cex=0.7, adj=0)
text(x= -98, y= 1030, label="Consensus timetree", font=2, cex=0.7, adj=0)
text(x= -98, y= 520, label="Backbone-and-patch", font=2, cex=0.7, adj=0)




# GET TREE TO JOIN
####
mamMCC1<-read.nexus(file="DR_calcsOnOther_mamPhys/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_target.tre")
mamMCC <- ladderize(drop.tip(mamMCC1, "_Anolis_carolinensis"))

mamData<-read.csv(file="traits_datReady_mamPhy_5911species.csv",header=TRUE)[,c(1:11)]
	
ordNames<-names(table(mamData$ord))
ordReps<-as.vector(mamData[match(ordNames,mamData$ord),"tiplabel"])
ordTable<-cbind.data.frame(ordNames,ordReps)

toDrop<-setdiff(mamMCC$tip.label,ordReps)
mamMCC_ords<-drop.tip(mamMCC,toDrop)
newTipNames<-do.call(rbind,strsplit(mamMCC_ords$tip.label,"_"))[,4]
mamMCC_ords$tip.label<-newTipNames

pdf(file="ordTree_forDivTimeCompare.pdf",height=18,width=8)
plot(mamMCC_ords, show.tip.label=TRUE, edge.width=4, cex=1)
dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy MS1 -- Upham et al. 2019 -- PLOS Biology
###
# Figure 1 - large circular phylogeny of Mammalia
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plot phylogeny with tip DR colored on branches
###
setwd("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
dir<-"_R-CODE/source_functions/"
source(paste0(dir,"circleTree_plotting_functions.R"))
library(RColorBrewer); library(ape); library(phangorn)

# plotting tree
bbone <- "NDexp" #"FBDasZhouEtAl"
mamMCC<-drop.tip(read.nexus(file=paste0("MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_",bbone,"_MCC_v2_target.tre")),"_Anolis_carolinensis")
plottree <- ladderize(mamMCC, right=TRUE)

rootAge=max(node.depth.edgelength(plottree))

# read in the data for colouring branches
tipDR_dat<-read.table("DR-SUMMARY_MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_all10k_v2_expanded.txt", header=TRUE)
cladesDR<-read.csv(paste0("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.csv"), header=TRUE)
head(cladesDR)

#===
# TAXON labels
# ===
# if all fams as nums
famNumAll <- cladesDR$famNumAll[match(plottree$tip.label,cladesDR$tiplabel)]

# ords, some full, some numbers
ordNums<- cladesDR$ordNums[match(plottree$tip.label,cladesDR$tiplabel)]
ordNumsSelect<- cladesDR$ordNumsSelect[match(plottree$tip.label,cladesDR$tiplabel)]

# major clades with numbers
cladeOrdNums<- cladesDR$cladeOrdNums[match(plottree$tip.label,cladesDR$tiplabel)]

# higher taxa
higherTax<- as.vector(cladesDR$higher[match(plottree$tip.label,cladesDR$tiplabel)])
higherTax[which(higherTax=="Afroth.")]<-"Afro."
higherTax[which(higherTax=="Xenart.")]<-"X."

#===
# Tip DR branch colors
# ===
# load DR
DR <- tipDR_dat[-1,"harmMeans"]
names(DR) <- rownames(tipDR_dat[-1,])

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
	
	cols<-rich.colors(100)
	range <- quantile(reconRate, seq(0,1, 0.01))[2:101]
    for (i in 1:100) {
        if (i==100) {range[i] <- range[i]+0.1}
        if (i==1) {reconColors[which(reconRate <= range[i])] <- cols[i] } 
        if (i > 1) {reconColors[which((reconRate <= range[i]) & (reconRate > range[i-1]))] <- cols[i]}
        }

    # plotting & tabling
    th <- max(branching.times(plottree))
    gr.col <- gray.colors(2, start=0.90, end=0.95)

# color tips for missing species
# ===
missing<-cladesDR[which(cladesDR$genes==0),"tiplabel"]
sampled<-as.vector(setdiff(plottree$tip.label,missing))
tipColors <- rep("black",5911)
tipColors[match(missing,plottree$tip.label)] <- grey(0.6) 


# WITH tips... DR Dens in UPPER CENTER
###
#pdf(file="DR_onMamPhy_BDvr_pcsFIXED_NDexp_5911taxa_MCCtarget_UpperDRdens_richCol_topoConsNew_v2_lagosOK_higher.pdf", width=25, height=25, onefile=TRUE)
#pdf(file="DR_onMamPhy_BDvr_pcsFIXED_NDexp_5911taxa_MCCtarget_UpperDRdens_richCol_topoConsNew_v2_lagosOK_higherAll.pdf", width=25, height=25, onefile=TRUE)
pdf(file="DR_onMamPhy_BDvr_pcsFIXED_NDexp_5911taxa_MCCtarget_UpperDRdens_richCol_topoConsNew_v2_lagosOK_higherAll_noLab.pdf", width=25, height=25, onefile=TRUE)

	XX1=-190
	XX2=190
	par(xpd=NA)

	plot(plottree, show.tip.label=FALSE, type="f", edge.width=0.9, cex=0.04, open.angle=45,no.margin=FALSE, root.edge=TRUE, edge.color="white", plot=FALSE, x.lim=c(XX1,XX2), y.lim=c(XX1,XX2))

	draw.circle(0,0, radius=th, col =gr.col[2], border=gr.col[2])
	draw.circle(0,0, radius=th-23.03, col=gr.col[1], border=gr.col[1])
	draw.circle(0,0, radius=th-66, col =gr.col[2], border=gr.col[2])
	draw.circle(0,0, radius=th-66, col =NA, lty=2, lwd=3, border="black")

	lines<-rev(seq(0,170,10))
	for (i in 1:length(lines)){
		draw.circle(0,0, radius=lines[i],col=NA,lty=2, lwd=2, border=grey(0.8,alpha=0.6))
	}

	text(y=-th+95.5, x=3, labels="Crown\nPlacentalia", font=2, cex=1.5, srt=0)
	text(y=-th+170, x=80, labels="Crown\nMarsupialia", font=2, cex=1.5, srt=0)

	text(y=3, x=seq(20,160,20), labels=rev(seq(20,160,20)), font=2, cex=1.5)
	text(y=10, x=th-66, labels="K-Pg\nboundary", font=2, cex=1.5)

	par(new=T)
	obj <- plot2.phylo(plottree, show.tip.label=TRUE, tip.color=tipColors, cex=0.05, label.offset=0.08, type="f", edge.width=0.9, no.margin=FALSE, root.edge=TRUE, edge.color=as.matrix(reconColors),x.lim=c(XX1,XX2), y.lim=c(XX1,XX2))

	par(new=T)

		#group.label.tip.rad3(obj, lab=famNumAll, c( "light grey", grey(0.7) ), "black", offset.bar=1, offset.lab=2, cex=1, lwd=4, arc.bar.width=1.05)
#		group.label.tip.rad4(obj, lab=higherTax, c( "black",grey(0.4) ), "white", offset.bar=11, offset.lab=15, cex=1.8,font=1, lwd=4, arc.bar.width=1.04)
#		group.label.tip.rad4(obj, lab=cladeOrdNums, c( "black", grey(0.4) ), "black", offset.bar=27.5, offset.lab=38, cex=1.8, font=1, lwd=4, arc.bar.width=1.034)
#		group.label.tip.rad3(obj, lab=ordNumsSelect, c( "light grey", grey(0.6) ), "black", offset.bar=19.3, offset.lab=20, cex=1.8, font=1,lwd=4, arc.bar.width=1.037)

		group.label.tip.rad4(obj, lab=higherTax, c( "black",grey(0.4) ), col.lab=NA, offset.bar=11, offset.lab=15, cex=1.8,font=1, lwd=4, arc.bar.width=1.04)
		group.label.tip.rad4(obj, lab=cladeOrdNums, c( "black", grey(0.4) ), col.lab=NA, offset.bar=27.5, offset.lab=38, cex=1.8, font=1, lwd=4, arc.bar.width=1.034)
		group.label.tip.rad3(obj, lab=ordNumsSelect, c( "light grey", grey(0.6) ), col.lab=NA, offset.bar=19.3, offset.lab=20, cex=1.8, font=1,lwd=4, arc.bar.width=1.037)

	# for ND normal ##
	color.legend(-51.5, 30.5, 0, 25.5, legend=NULL, rect.col= cols[1], gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey") #xl, yb, xr, yt
	color.legend(0, 30.5, 88, 25.5, legend=NULL, rect.col= cols[100], gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey") #xl, yb, xr, yt

	color.legend(-44.5, 30.5, 45.2, 25.5, legend=NULL, rect.col= cols, gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey") #xl, yb, xr, yt
	text(x=4, y=13.5, "Tip speciation rate (species/Ma)", font=1, cex=2) 

	par(fig=c(0.385, 0.72, 0.57, 0.695),mar=c(0,0,0,0), mar=c(0,0,0,0), new=T) #c(x1, x2, y1, y2)

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

dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy MS1 -- Upham et al. 2019 -- PLOS Biology
###
# Figure 5 - comparing LTT plots among Mammalia species-level studies
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# PLOTTING OTHER BIG TREES
####
setwd("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
source("source_functions/DR_circleTree_functions.R")
library(ape); library(phytools); library(viridis); library(phangorn); library(geiger); library(magick); library(phyloch)
library(dplyr)
##
# All BIG TREES TOGETHER in a loop, same timescale and color scale.

# Load trees
dir<-"_DATA/"
mamMCC1<-read.nexus(file=paste0(dir,"MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_target.tre"))
mamMCC <- ladderize(drop.tip(mamMCC1, "_Anolis_carolinensis"))
	#r1<-max(nodeHeights(plottree1))

BE_etAl_phy1<-read.nexus(file=paste0(dir,"bininda-emonds-et-al2007_treeFile_Supp1.nex"))
BE_etAl_phy <- ladderize(BE_etAl_phy1[[1]])
	#rX<-max(nodeHeights(plottreeX))

KuhnEtAl_100<-read.nexus(file=paste0(dir,"FritzTree.rs200k.100trees.nex"))

mamPhy_100<-read.nexus(file=paste0(dir,"MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_sample100_nexus.trees"))

FS2015_100<-read.nexus(file=paste0(dir,"FS2015_Fully_resolved_phylogeny_random100trees.nex")) # 100 trees, 5747 taxa each

Hedges_phySmooth<-read.tree(file=paste0(dir,"9.TTOL_mammals_smoothed_interpolated.tre"))


###
# Make the LTT comparisons plot
########
# get PER-STUDY colors...
colKuhn<-viridis(10, alpha=0.4)[1]
colBE<-viridis(10, alpha=1)[2]
colFS<-viridis(10, alpha=0.4)[6]
colHedges<-viridis(10, alpha=1)[9]

colThisStudy<-rgb((as.vector(col2rgb("goldenrod1"))/255)[1],(as.vector(col2rgb("goldenrod1"))/255)[2],(as.vector(col2rgb("goldenrod1"))/255)[3],alpha=0.5)

RESDIR<-"Fig5_compare_LTTs/"

png(file=paste0(RESDIR,"comparing_MamPhys_LTT_all_BE-Kuhn-FS-Hedges_viridis_legByType_withMamPhy.png"), width=6, height=6, res=400, units="in")
#	ltt.plot(KuhnEtAl_100[[1]], log="y", col=rgb(0,0,1,alpha=0.5), lwd=1, xlim=c(-100,0), bty="n", ylab="", xlab="", xaxt="n",yaxt="n")#, ylim=c(4,5500))
	# KUHN
	ltt.plot(KuhnEtAl_100[[1]], log="y", col=rgb(0,0,1,alpha=0.5), lwd=1, xlim=c(-100,0), bty="n", ylab="log (Lineages)", xlab="Time before present (Ma)")#, ylim=c(4,5500))
	#axis(side=1, at=c(0,-20,-40,-60,-80,-100),labels=c(0,-20,-40,-60,-80,-100))
	for(i in 2:length(KuhnEtAl_100)){
		ltt.lines(KuhnEtAl_100[[i]], col=colKuhn, lwd=1) #rgb(0,0,1,alpha=0.5)
	}
	ltt.lines(BE_etAl_phy, col=colBE, lwd=2) #"purple"
	# MAMPHY
	for(i in 1:length(mamPhy_100)){
		ltt.lines(mamPhy_100[[i]], col=colThisStudy, lwd=1) #grey(0.5,alpha=0.4)
	}
	ltt.lines(mamMCC, col="black", lwd=2)
	# FS2015
	for(i in 1:length(FS2015_100)){
		ltt.lines(FS2015_100[[i]], col=colFS, lwd=1) #rgb(1,0,0,alpha=0.4)
	}
	par(new=TRUE)
	#HEDGES
	ltt.lines(Hedges_phySmooth, col=colHedges, lwd=2) #"goldenrod1"
	rect(xleft= 0, ybottom=1, xright=10, ytop=6000,  angle = 90, col = "white", border = "white", lwd=3)

#ltt.plot(KuhnEtAl_100[[1]], log="y", col=rgb(0,0,1,alpha=0.5), lwd=1, xlim=c(-100,0), bty="n",ylab="log ( Lineages )")#, ylim=c(4,5500))
#legend(x= -100, y= 6500, legend=c(" ","Supertree w/ polytomies", "","Supertree resolved credible set"," ", 
#	"Supertree resolved credible set", " ", "Consensus timetree smoothed", " ",
#	"Backbone-and-patch credible set", "Backbone-and-patch MCC"), 
#	col=c(NA, colBE, NA, colKuhn,NA, colFS,NA, colHedges,NA, "goldenrod1","black"), #col=c(NA, "purple", NA, "blue",NA,"red",NA,"goldenrod1",NA,grey(0.5),"black"), 
#	lwd=2, lty=1, cex=0.7)
#text(x= -98, y= 4900, label="Bininda-Emonds et al. 2007", font=2, cex=0.7, adj=0)
#text(x= -98, y= 2550, label="Kuhn et al. 2011", font=2, cex=0.7, adj=0)
#text(x= -98, y= 1350, label="Faurby and Svenning 2015", font=2, cex=0.7, adj=0)
#text(x= -98, y= 720, label="Hedges et al. 2015", font=2, cex=0.7, adj=0)
#text(x= -98, y= 400, label="This study", font=2, cex=0.7, adj=0)

legend(x= -100, y= 6500, legend=c(" ","Bininda-Emonds et al. 2007 (one tree)", "Kuhn et al. 2011 (set of trees)"," ","Faurby & Svenning 2015 (set of trees)", 
	" ", "Hedges et al. 2015 (one tree)", " ",
	"This study (set of trees)", "This study (MCC tree)"), 
	col=c(NA, colBE, colKuhn,NA, colFS,NA, colHedges,NA, "goldenrod1","black"), #col=c(NA, "purple", NA, "blue",NA,"red",NA,"goldenrod1",NA,grey(0.5),"black"), 
	lwd=2, lty=1, cex=0.7)
text(x= -98, y= 4900, label="MRP supertree", font=2, cex=0.7, adj=0)
text(x= -98, y= 1850, label="DNA supertree", font=2, cex=0.7, adj=0)
text(x= -98, y= 1030, label="Consensus timetree", font=2, cex=0.7, adj=0)
text(x= -98, y= 520, label="Backbone-and-patch", font=2, cex=0.7, adj=0)

dev.off()


##
# Plot just the RODENTIA portions...
###
# plot each big tree with node labels numbered and TIPS


pdf(file=paste0(RESDIR,"comparing_MamPhys_fullplots_nodeNumsTips.pdf"), width=120, height=100, onefile=TRUE)

par(oma = rep(3,4) + 0.1,
    mar = c(4,1.5,4,1) + 0.1) #‘c(bottom, left, top, right)’

layout(matrix(1:6,1,6, byrow=TRUE))#,widths=rep(0.25,4),heights=rep(,4))

	plot(ladderize(KuhnEtAl_100[[1]]), cex=0.2, label.offset=1)
	nodelabels(cex=0.2)
	mtext(side=3, "Kuhn 2011 tree 1 of 100", cex=5)

	plot(ladderize(BE_etAl_phy), cex=0.2, label.offset=1)
	nodelabels(cex=0.2)
	mtext(side=3, "BE2007 1 tree ", cex=5)

	plot(ladderize(FS2015_100[[1]]), cex=0.2, label.offset=1)
	nodelabels(cex=0.2)
	mtext(side=3, "FS2015 tree 1 of 100", cex=5)

	plot(ladderize(Hedges_phySmooth), cex=0.2, label.offset=1)
	nodelabels(cex=0.2)
	mtext(side=3, "Hedges 1 tree", cex=5)

	plot(ladderize(mamPhy_100[[1]]), cex=0.2, label.offset=1)
	nodelabels(cex=0.2)
	mtext(side=3, "MamPhy tree 1 of 100", cex=5)

	plot(ladderize(mamMCC), cex=0.2, label.offset=1)
	nodelabels(cex=0.2)
	mtext(side=3, "MamPhy MCC", cex=5)

dev.off()

# GET TAXON LISTS & PRUNE
###
	# BE_etAl_phy -- 
		# Rodentia -- node 4520
		toKeep<-as.vector(na.omit(BE_etAl_phy$tip.label[getDescendants((BE_etAl_phy), node=4520)]))
		toDrop<-setdiff(BE_etAl_phy$tip.label,toKeep)
		BE_Rodentia<-drop.tip(BE_etAl_phy,toDrop)
		# Chiroptera -- node 5842
		toKeep<-as.vector(na.omit(BE_etAl_phy$tip.label[getDescendants((BE_etAl_phy), node=5842)]))
		toDrop<-setdiff(BE_etAl_phy$tip.label,toKeep)
		BE_Chiroptera<-drop.tip(BE_etAl_phy,toDrop)
		# Primates -- node 5259
		toKeep<-as.vector(na.omit(BE_etAl_phy$tip.label[getDescendants((BE_etAl_phy), node=5259)]))
		toDrop<-setdiff(BE_etAl_phy$tip.label,toKeep)
		BE_Primates<-drop.tip(BE_etAl_phy,toDrop)
	
	# Hedges_phySmooth -- 
		# Rodentia -- node 5371
		toKeep<-as.vector(na.omit(Hedges_phySmooth$tip.label[getDescendants((Hedges_phySmooth), node=5371)]))
		toDrop<-setdiff(Hedges_phySmooth$tip.label,toKeep)
		Hedges_Rodentia<-drop.tip(Hedges_phySmooth,toDrop)
		# Chiroptera -- node 8027
		toKeep<-as.vector(na.omit(Hedges_phySmooth$tip.label[getDescendants((Hedges_phySmooth), node=8027)]))
		toDrop<-setdiff(Hedges_phySmooth$tip.label,toKeep)
		Hedges_Chiroptera<-drop.tip(Hedges_phySmooth,toDrop)
		# Primates -- node 7581
		toKeep<-as.vector(na.omit(Hedges_phySmooth$tip.label[getDescendants((Hedges_phySmooth), node=7581)]))
		toDrop<-setdiff(Hedges_phySmooth$tip.label,toKeep)
		Hedges_Primates<-drop.tip(Hedges_phySmooth,toDrop)

	# MamPhy -- 
		# Rodentia tag
		toKeep<-mamMCC$tip.label[grep(pattern="RODENTIA", x=mamMCC$tip.label, ignore.case=FALSE)]
		toDrop<-setdiff(mamMCC$tip.label,toKeep)
		MamMCC_Rodentia<-drop.tip(mamMCC,toDrop)
		MamPhy100_Rodentia<-list()
		for(i in 1:length(mamPhy_100)){
			MamPhy100_Rodentia[[i]]<-drop.tip(mamPhy_100[[i]],toDrop)
		}
		class(MamPhy100_Rodentia)<-"multiPhylo"
		# Chiroptera tag
		toKeep<-mamMCC$tip.label[grep(pattern="CHIROPTERA", x=mamMCC$tip.label, ignore.case=FALSE)]
		toDrop<-setdiff(mamMCC$tip.label,toKeep)
		MamMCC_Chiroptera<-drop.tip(mamMCC,toDrop)
		MamPhy100_Chiroptera<-list()
		for(i in 1:length(mamPhy_100)){
			MamPhy100_Chiroptera[[i]]<-drop.tip(mamPhy_100[[i]],toDrop)
		}
		class(MamPhy100_Chiroptera)<-"multiPhylo"
		# Primates tag
		toKeep<-mamMCC$tip.label[grep(pattern="PRIMATES", x=mamMCC$tip.label, ignore.case=FALSE)]
		toDrop<-setdiff(mamMCC$tip.label,toKeep)
		MamMCC_Primates<-drop.tip(mamMCC,toDrop)
		MamPhy100_Primates<-list()
		for(i in 1:length(mamPhy_100)){
			MamPhy100_Primates[[i]]<-drop.tip(mamPhy_100[[i]],toDrop)
		}
		class(MamPhy100_Primates)<-"multiPhylo"

	# Kuhn 2011
		# RODENTIA
		toKeep<-as.vector(unlist(read.table(file=paste0(dir,"KuhnEtAl2011_TAXA_Rodentia.txt"))))
		toDrop<-setdiff(KuhnEtAl_100[[1]]$tip.label,toKeep)
		Kuhn100_Rodentia<-list()
		for(i in 1:100){
			Kuhn100_Rodentia[[i]]<-drop.tip(KuhnEtAl_100[[i]],toDrop)
		}
		class(Kuhn100_Rodentia)<-"multiPhylo"

		# CHIROPTERA
		toKeep<-as.vector(unlist(read.table(file=paste0(dir,"KuhnEtAl2011_TAXA_Chiroptera.txt"))))
		toDrop<-setdiff(KuhnEtAl_100[[1]]$tip.label,toKeep)
		Kuhn100_Chiroptera<-list()
		for(i in 1:100){
			Kuhn100_Chiroptera[[i]]<-drop.tip(KuhnEtAl_100[[i]],toDrop)
		}
		class(Kuhn100_Chiroptera)<-"multiPhylo"

		# PRIMATES
		toKeep<-as.vector(unlist(read.table(file=paste0(dir,"KuhnEtAl2011_TAXA_Primates.txt"))))
		toDrop<-setdiff(KuhnEtAl_100[[1]]$tip.label,toKeep)
		Kuhn100_Primates<-list()
		for(i in 1:100){
			Kuhn100_Primates[[i]]<-drop.tip(KuhnEtAl_100[[i]],toDrop)
		}
		class(Kuhn100_Primates)<-"multiPhylo"
	

	# FS2015
		# RODENTIA
		toKeep<-as.vector(unlist(read.table(file=paste0(dir,"FS2015_TAXA_Rodentia.txt"))))
		toDrop<-setdiff(FS2015_100[[1]]$tip.label,toKeep)
		FS100_Rodentia<-list()
		for(i in 1:100){
			FS100_Rodentia[[i]]<-drop.tip(FS2015_100[[i]],toDrop)
		}
		class(FS100_Rodentia)<-"multiPhylo"

		# CHIROPTERA
		toKeep<-as.vector(unlist(read.table(file=paste0(dir,"FS2015_TAXA_Chiroptera.txt"))))
		toDrop<-setdiff(FS2015_100[[1]]$tip.label,toKeep)
		FS100_Chiroptera<-list()
		for(i in 1:100){
			FS100_Chiroptera[[i]]<-drop.tip(FS2015_100[[i]],toDrop)
		}
		class(FS100_Chiroptera)<-"multiPhylo"

		# PRIMATES
			# Fixing the FS2015 PRIMATES...NEEDED because TARSIERS are non-ultrametric in all their trees! (need to email them the correction)
			library(phangorn)
		toKeep<-as.vector(unlist(read.table(file=paste0(dir,"FS2015_TAXA_Primates.txt"))))
		toDrop<-setdiff(FS2015_100[[1]]$tip.label,toKeep)
		FS100_Primates<-list()
		for(i in 1:100){
			toBeScaled<-drop.tip(FS2015_100[[i]],toDrop)
			FS100_Primates[[i]]<-nnls.tree(cophenetic(toBeScaled),toBeScaled,rooted=TRUE)
		}
		class(FS100_Primates)<-"multiPhylo"


# PLOT the PER-STUDY-based LTTs by 3 taxa...
	colRod<-plasma(10, alpha=0.3)[2]
	colChir<-plasma(10, alpha=0.3)[5]
	colPri<-plasma(10, alpha=0.3)[8]
	colRod_2<-plasma(10, alpha=0.9)[2]
	colChir_2<-plasma(10, alpha=0.9)[5]
	colPri_2<-plasma(10, alpha=0.9)[8]

# COMBO
png(file=paste0(RESDIR,"comparing_MamPhys_LTT_all_COMBO_ROD-CHIR-PRI_plasma.png"), width=6, height=6, res=600, units="in")
#pdf(file=paste0(RESDIR,"comparing_MamPhys_LTT_all_COMBO_ROD-CHIR-PRI_plasma.pdf"), width=6, height=6)

	par(oma = c(4,4,0.1,0.1) + 0.1,
	    mar = c(0.5,0.5,1,0.5) + 0.1) #‘c(bottom, left, top, right)’

	layout(matrix(1:4,2,2, byrow=TRUE))#,widths=rep(0.25,4),heights=rep(,4))

	# THIS STUDY (MamPhy)
	###
		#png(file=paste0(RESDIR,"comparing_MamPhys_LTT_all_THISSTUDY_ROD-CHIR-PRI_plasma.png"), width=6, height=6, res=600, units="in")
		#pdf(file=paste0(RESDIR,"comparing_MamPhys_LTT_all_THISSTUDY_ROD-CHIR-PRI_plasma.pdf"), width=6, height=6, onefile=TRUE)#, res=400, units="in")
		#	ltt.plot(KuhnEtAl_100[[1]], log="y", col=rgb(0,0,1,alpha=0.5), lwd=1, xlim=c(-100,0), bty="n", ylab="", xlab="", xaxt="n",yaxt="n")#, ylim=c(4,5500))
			# RODENTIA
			#ltt.plot(MamPhy100_Rodentia[[1]], log="y", col=colRod, lwd=1, xlim=c(-100,0), bty="n", ylab="log (Lineages)", xlab="Time before present (Ma)")#, ylim=c(4,5500))
			ltt.plot(MamPhy100_Rodentia[[1]], log="y", col=colRod, lwd=1, xlim=c(-100,0), bty="n", xaxt="n",ylab="", xlab="")#, ylim=c(4,5500))
			axis(side=1, at=c(0,-20,-40,-60,-80,-100),labels=NA)#c(0,-20,-40,-60,-80,-100))
			for(i in 2:100){
				ltt.lines(MamPhy100_Rodentia[[i]], col=colRod, lwd=1)
			}	
			# CHIROP
			for(i in 1:100){
				ltt.lines(MamPhy100_Chiroptera[[i]], col=colChir, lwd=1)
			}
			# PRIMATES
			for(i in 1:100){
				ltt.lines(MamPhy100_Primates[[i]], col=colPri, lwd=1)
			}
			#	ltt.lines(MamMCC_Rodentia, col=colRod_2, lwd=3)
			#	ltt.lines(MamMCC_Chiroptera, col=colChir_2, lwd=3)
			#	ltt.lines(MamMCC_Primates, col=colPri_2, lwd=3)
			#par(new=TRUE)
			#rect(xleft= -80, ybottom=1, xright=10, ytop=2, angle = 90, col = "white", border = "white", lwd=0.1)

		#dev.off()

	# FS2015
	###
	#png(file=paste0(RESDIR,"comparing_MamPhys_LTT_all_FS2015_ROD-CHIR-PRI_plasma.png"), width=6, height=6, res=600, units="in")
	#pdf(file=paste0(RESDIR,"comparing_MamPhys_LTT_all_FS2015_ROD-CHIR-PRI_plasma.pdf"), width=6, height=6)#, res=400, units="in")
	#	ltt.plot(KuhnEtAl_100[[1]], log="y", col=rgb(0,0,1,alpha=0.5), lwd=1, xlim=c(-100,0), bty="n", ylab="", xlab="", xaxt="n",yaxt="n")#, ylim=c(4,5500))
		# RODENTIA
		#ltt.plot(FS100_Rodentia[[1]], log="y", col=colRod, lwd=1, xlim=c(-100,0), bty="n",ylab="log (Lineages)", xlab="Time before present (Ma)")#, ylim=c(4,5500))
		ltt.plot(FS100_Rodentia[[1]], log="y", col=colRod, lwd=1, xlim=c(-100,0), bty="n", xaxt="n",yaxt="n",ylab="", xlab="")#, ylim=c(4,5500))
		axis(side=1, at=c(0,-20,-40,-60,-80,-100),labels=NA)#c(0,-20,-40,-60,-80,-100))
		axis(side=2, at=c(1,5,10,50,100,500,1000),labels=NA)#c(0,-20,-40,-60,-80,-100))
		
		for(i in 2:length(FS100_Rodentia)){
			ltt.lines(FS100_Rodentia[[i]], col=colRod, lwd=1)
		}	# CHIROP
		for(i in 1:length(FS100_Chiroptera)){
			ltt.lines(FS100_Chiroptera[[i]], col=colChir, lwd=1)
		}

		# PRIMATES
		for(i in 1:length(FS100_Primates)){
			ltt.lines(FS100_Primates[[i]], col=colPri, lwd=1)
		}

		par(new=TRUE)
		rect(xleft= 0, ybottom=1, xright=10, ytop=6000,  angle = 90, col = "white", border = "white", lwd=3)

	#dev.off()

	# HEDGES2015
	###
	#png(file=paste0(RESDIR,"comparing_MamPhys_LTT_all_HEDGES2015_ROD-CHIR-PRI_plasma.png"), width=6, height=6, res=600, units="in")
	#pdf(file=paste0(RESDIR,"comparing_MamPhys_LTT_all_FS2015_ROD-CHIR-PRI_plasma.pdf"), width=6, height=6)#, res=400, units="in")
	#	ltt.plot(KuhnEtAl_100[[1]], log="y", col=rgb(0,0,1,alpha=0.5), lwd=1, xlim=c(-100,0), bty="n", ylab="", xlab="", xaxt="n",yaxt="n")#, ylim=c(4,5500))
		# RODENTIA
		#ltt.plot(Hedges_Rodentia, log="y", col=colRod_2, lwd=3, xlim=c(-100,0), bty="n", ylab="log (Lineages)", xlab="Time before present (Ma)")#, ylim=c(4,5500))
		ltt.plot(Hedges_Rodentia, log="y", col=colRod_2, lwd=2, xlim=c(-100,0), bty="n", ylab="", xlab="")#, ylim=c(4,5500))
			ltt.lines(Hedges_Chiroptera, col=colChir_2, lwd=2)
			ltt.lines(Hedges_Primates, col=colPri_2, lwd=2)
	#dev.off()

	# KUHN2011 and BE2007
	###
	#png(file=paste0(RESDIR,"comparing_MamPhys_LTT_all_KUHN2011_ROD-CHIR-PRI_plasma.png"), width=6, height=6, res=600, units="in")
	#pdf(file=paste0(RESDIR,"comparing_MamPhys_LTT_all_FS2015_ROD-CHIR-PRI_plasma.pdf"), width=6, height=6)#, res=400, units="in")
	#	ltt.plot(KuhnEtAl_100[[1]], log="y", col=rgb(0,0,1,alpha=0.5), lwd=1, xlim=c(-100,0), bty="n", ylab="", xlab="", xaxt="n",yaxt="n")#, ylim=c(4,5500))
		# RODENTIA
		#ltt.plot(Kuhn100_Rodentia[[1]], log="y", col=colRod, lwd=1, xlim=c(-100,0), bty="n", ylab="log (Lineages)", xlab="Time before present (Ma)")#, ylim=c(4,5500))
		ltt.plot(Kuhn100_Rodentia[[1]], log="y", col=colRod, lwd=1, xlim=c(-100,0), bty="n", yaxt="n", ylab="", xlab="")#, ylim=c(4,5500))
		#axis(side=1, at=c(0,-20,-40,-60,-80,-100),labels=NA)#c(0,-20,-40,-60,-80,-100))
		axis(side=2, at=c(1,5,10,50,100,500,1000),labels=NA)#c(0,-20,-40,-60,-80,-100))
			ltt.lines(BE_Rodentia, col="grey", lwd=3)
			ltt.lines(BE_Chiroptera, col="grey", lwd=3)
			ltt.lines(BE_Primates, col="grey", lwd=3)
		for(i in 2:100){
			ltt.lines(Kuhn100_Rodentia[[i]], col=colRod, lwd=1)
		}	
		# CHIROP
		for(i in 1:100){
			ltt.lines(Kuhn100_Chiroptera[[i]], col=colChir, lwd=1)
		}
		# PRIMATES
		for(i in 1:100){
			ltt.lines(Kuhn100_Primates[[i]], col=colPri, lwd=1)
		}
	
	title(main="",sub="", xlab="Time before present (Ma)", ylab="log (Lineages)", line=2, 
		outer=TRUE, cex.axis=1,cex.lab=1.2,font.axis=1, font.lab=1)
	#dev.off()
dev.off()







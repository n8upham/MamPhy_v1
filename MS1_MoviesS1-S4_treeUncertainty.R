#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy MS1 -- Upham et al. 2019 -- PLOS Biology
###
# Movies S1-S4 - animated gifs showing tree uncertainty int eh credible sets of Mammalia
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####
# LOAD the 100 mammal trees of interest.
###
dirname<-"/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/11_makingFullPosteriors/10k_tree_andDR_distributions_FINAL"
setwd(dirname)
library(ape); library(phyloch); library(phytools)

#mamPhy_100<-read.nexus(file="MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_v2_sample100_nexus.trees")
mamPhy_100<-read.nexus(file="MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_v2_sample100_nexus.trees")
#mamPhy_100<-read.nexus(file="MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_v2_sample100_nexus.trees")
#mamPhy_100<-read.nexus(file="MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_v2_sample100_nexus.trees")
#TYPE<-"Completed-NDexp"
#TYPE<-"DNA-only-NDexp"
#TYPE<-"Completed-FBDasZhouEtAl"
TYPE<-"DNA-only-FBDasZhouEtAl"

cladesDR_all<-read.table("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_NDexp_DRstats_DRtreeLABELS.txt",header=TRUE)
#	cladesDR<-cladesDR_all[which(cladesDR_all$genes > 0),]
	cladesDR<-cladesDR_all

# which tips are missing?
missingSp<-as.vector(cladesDR_all[which(cladesDR_all$genes == 0),"tiplabel"])

# get a list of species per higher taxon
higher<-names(table(cladesDR$higher))[c(5,4,1,6,2,3)]
spHigher<-list()
for(j in 1:length(higher)){
	spHigher[[j]]<-as.vector(cladesDR[which(cladesDR$higher==higher[j] & cladesDR$genes > 0),"tiplabel"])
}

# prune trees and get root heights
treesToPlot<-vector("list", length(mamPhy_100))
subRoots<-c()
for(i in 1:length(mamPhy_100)){
	treesToPlot[[i]]<-drop.tip(mamPhy_100[[i]], "_Anolis_carolinensis")
	subRoots[i]<-max(node.depth.edgelength(phy=treesToPlot[[i]]))
}
maxOfAll<-max(subRoots)
maxOfAll_ID<-match(maxOfAll,subRoots)


# set node colors
library(viridis); library(gplots); 
source("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/_R-CODE/source_functions/circleTree_plotting_functions.R")
#nodeCols<-cividis(6)#10)[c(0,2,4,6,8,10)]
#nodeCols<-rich.colors(6)#10)[c(0,2,4,6,8,10)]
nodeCols<-magma(6)#10)[c(0,2,4,6,8,10)]

RESDIR<-"gif_100trees_v2"

# NOW ANIMATE
install.packages("magick")
library(magick)

img <- image_graph(height=800, width=600, res = 96, bg="white")
#pdf(file=paste0(RESDIR,"/credibleSet_mamPhy_",TYPE,"_100trees_all_higherNodeCols_tipLabel_v2.pdf"),height=8,width=6, onefile = TRUE)

for(i in 1:100){

	tree_i<-ladderize(treesToPlot[[i]])
	# get the nodes to mark
	nodeSet<-c()
	for(j in 1:length(spHigher)){
		nodeSet[j]<-getMRCA(phy=treesToPlot[[i]], tip=spHigher[[j]])
	}
	#plot oldest tree and axis
	plot(ladderize(treesToPlot[[maxOfAll_ID]]), edge.width=0.5, show.tip.label=FALSE, cex=0.3, edge.color=NA, type="phylogram")
	axis(side=1,cex.lab=0.4,at=(maxOfAll-c(0,50,100,150,200)),labels=(c(0,50,100,150,200)))

	subRoot_BASE<-max(node.depth.edgelength(treesToPlot[[maxOfAll_ID]]))

#NDexp
#	rootTime<-200
#FBD
	rootTime<-300

	# plot legend
# COMPLETED
#	legend(x=subRoot_BASE-rootTime, y=2000, legend=c("Monotremata","Marsupialia","Afrotheria","Xenarthra","Euarchontoglires","Laurasiatheria","","Imputed species"),
#		pch=c(rep(21,6),NA,45), pt.bg=c(nodeCols,"white",grey(0.6)), col=c(rep("black",6),"white",grey(0.6)), cex=c(rep(1,6),1,1),pt.cex=c(rep(2,6),1,2),text.col=c(rep("black",6),"white",grey(0.4)) )
# DNA-only
	legend(x=subRoot_BASE-rootTime, y=2000, legend=c("Monotremata","Marsupialia","Afrotheria","Xenarthra","Euarchontoglires","Laurasiatheria"),
		pch=c(rep(21,6)), pt.bg=c(nodeCols), col=c(rep("black",6)), cex=c(rep(1,6)),pt.cex=c(rep(2,6)),text.col=c(rep("black",6)) )
	text(x=subRoot_BASE-rootTime, y=2200, labels=paste0("tree ",i), adj=0)

	#plot target tree
	## get last phylo plot parameters
	obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	subRoot_TARGET<-max(node.depth.edgelength(tree_i))
	par(new=TRUE)
	plot2.phylo((tree_i), edge.width=0.5, show.tip.label=FALSE, cex=0.3, edge.color="black", type="phylogram", 
		x.lim=c(subRoot_TARGET-obj$x.lim[2],subRoot_TARGET) )

	## get last phylo plot parameters
	obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	x.tip<-obj$xx[1:obj$Ntip]
	y.tip<-obj$yy[1:obj$Ntip]
	#tipsMissing<-as.vector(na.omit(match(missingSp,tree_i$tip.label)))#c(1:5911)[tree_i$tip.label %in% missingSp]
	#tipsMissing<-y.tip[tree_i$tip.label %in% missingSp]
#	tipsMissing<-c(1:obj$Ntip)[tree_i$tip.label %in% missingSp]

	#label nodes on that tree
	nodelabels(node=nodeSet, pch=21, col="black", bg=nodeCols, cex=3 )

	#label the TIPS for MISSING
	#tiplabels(cex=0.1)
#	tiplabels(tip=tipsMissing, frame="none",text="-----------",col=gray(0.2), cex=0.1, offset=5) # 
	#lines(y=cbind(tipsMissing,tipsMissing), x=cbind(rep(x.tip[1],length(tipsMissing)), rep(x.tip[1]+5,length(tipsMissing)) ), col="red", lty=1) # 
	#lines(y=tipsMissing, x=rep(x.tip[1],length(tipsMissing)), col="red", lty=1) # 

	#rug(side=4, x=tipsMissing, lwd=0.00001, col="black")#grey(0.6))
#	mtext(side=3, text="Node-dating backbone, 5911 species (completed)", font=2)
#	mtext(side=3, text="Node-dating backbone, 4098 species (DNA-only)", font=2)
#	mtext(side=3, text="Fossilized birth-death backbone, 5911 species + 76 fossils (completed)", font=2, cex=0.9)
	mtext(side=3, text="Fossilized birth-death backbone, 4098 species + 76 fossils (DNA-only)", font=2, cex=0.9)
	mtext(side=1, text="Millions of years before present (Ma)", cex=1, line=3)

}
dev.off()

animation <- image_animate(img, fps = 10)
print(animation)
image_write(animation, paste0(RESDIR,"/credibleSet_mamPhy_",TYPE,"_100trees_all_higherNodeCols_tipLabel_v2_tipSeq.gif"))




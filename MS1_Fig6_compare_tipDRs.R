#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy MS1 -- Upham et al. 2019 -- PLOS Biology
###
# Figure 6 - comparing tip DR calculations among studies
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Code to repeat the DR analyses on OTHER MAMMAL TREES...
#####
library(ape); library(picante); library(phytools); library(geiger)
setwd("/mnt/data/personal/nateu/Nate_Backup/DR_calcOnOther_mamPhys/")
source("DR_functions.R")

dir<-"_DATA/"

#=============
# ES and DR on per-tip basis
#######

# Do the HEDGES DR CALCS
##
# Load:: Hedges et al 2015 >> is ONE tree, no posterior !!!
##### #
#Hedges_phyUnsmooth<-read.tree(file=paste0(dir,"8.TTOL_mammals_unsmoothed.tre")) # 1 tree-- 3738 tips and 3737 internal nodes
#Hedges_phySmooth<-read.tree(file=paste0(dir,"9.TTOL_mammals_smoothed_interpolated.tre")) # 1 tree-- 5364 tips and 5363 internal nodes.
Hedges_phyUnsmooth_1<-scan(paste0(dir,"8.TTOL_mammals_unsmoothed.tre"), what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 
Hedges_phyUnsmooth<-strsplit(Hedges_phyUnsmooth_1,"[;]")  #trees now indexed #tree1=trees[[1]]  #first tree

Hedges_phySmooth<-scan(paste0(dir,"9.TTOL_mammals_smoothed_interpolated.tre"), what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 
Hedges_phySmooth<-strsplit(Hedges_phySmooth_1,"[;]")  #trees now indexed #tree1=trees[[1]]  #first tree

# gives pairwise clade matrix from CAIC function
clade_matrix = readCAIC(Hedges_phySmooth)

ES = ES_v2(clade_matrix)
DR = 1/ES
res = cbind.data.frame(DR,ES)
res1 = res[order(rownames(res)),]

write.table(res1, file=paste0(dir,"Hedges_phySmooth_TTOL_tipCalcs_1tree_DR-ES.txt"))


#########
# Now looping through the calcs of 1k posteriors...
library(foreach);library(doSNOW)
cl2 = makeCluster(15, type = 'SOCK', outfile="")
registerDoSNOW(cl2)
source("DR_functions.R")


# LOAD:: Faurby and Svenning 2015-- 1k trees
##
##### #setwd("/Users/Nate/Desktop/JETZ-pdfs/Faurby_Svenning2015_SUPP/AppendixD_producedPhylogenies_Rcode")
FS2015_phy1k<-read.tree(file=paste0(dir,"FS2015_Fully_resolved_phylogeny.nex")) # NOT in nexus, is a NEWICK formated file 1000 trees, 5747 taxa each
#FS2015_phy1k_i<-scan(paste0(dir,"FS2015_Fully_resolved_phylogeny.nex"), what="list",sep="\n",quiet=TRUE,n=1,skip=i+1,comment.char="#") 
FS2015_phy1k<-scan(paste0(dir,"FS2015_Fully_resolved_phylogeny.nex"), what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") # is newick format
#FS2015_phy1k<-strsplit(FS2015_phy1k_1,"[;]")  #trees now indexed #tree1=trees[[1]]  #first tree
trees<-FS2015_phy1k

# PRUNE these FS2015 trees to extant and recently extinct taxa (using the same recently extinct ones as Upham et al. 2019)
	FS_traits<-read.csv(paste0(dir,"AdditionalFile1_edited.csv"))
		FS_extinct<-as.character(FS_traits[which(FS_traits$Extinct.or.extant=="Extinct"), "Binomial.name"])

	mamData<-read.csv(file=paste0(dir,"taxonomy_mamPhy_5911species_toPublish.csv"),header=TRUE)
		mam_extinct<-as.character(mamData[which(mamData$extinct.==1), "Species_Name"])

	toDrop<-setdiff(FS_extinct,mam_extinct)

	FS2015_phy1k_dropped<-list()
	for(i in 1:length(FS2015_phy1k)){
		FS2015_phy1k_dropped[[i]]<-drop.tip(FS2015_phy1k[[i]],toDrop)
	}
	class(FS2015_phy1k_dropped)<-"multiPhylo"
	write.tree(FS2015_phy1k_dropped,file=paste0(dir,"FS2015_Fully_resolved_phylogeny_dropPleistoExtinct_5502species.nex"))

FS2015_phy1k_noPleisto<-scan(paste0(dir,"FS2015_Fully_resolved_phylogeny_dropPleistoExtinct_5502species.nex"), what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") # is newick format
trees<-FS2015_phy1k_noPleisto


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

write.table(g_es_dr, file=paste0(dir,"FS2015_resolved_tipCalcs_1kTrees_DR-ES_noPleistoExtinct.txt"))



# Load:: Kuhn et al. 2011 1k trees POSTERIOR of the Fritz et al. 2009 tree which is an update of Bininda-Emonds et al. 2007 tree
##
##### 
Kuhn_phy1k<-scan(paste0(dir,"Fritz.Resolved.Normal.cr1r5_v2b.nwk"), what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 
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

write.table(g_es_dr, file=paste0(dir,"Kuhn2011_resolved_tipCalcs_1kTrees_DR-ES.txt"))


#================
# SUMMARIZE the DR matrices into harmonic means of the posterior.
library(ape); library(picante); library(phytools); library(geiger)
setwd("/mnt/data/personal/nateu/Nate_Backup/DR_calcOnOther_mamPhys/")

# KUHN first.
Kuhn_DR1k<-read.table(file=paste0(dir,"Kuhn2011_resolved_tipCalcs_1kTrees_DR-ES.txt"), header=TRUE)

DR_table<-vector("list",length=1001)
for (i in 0:1000){
	DR_table[[i+1]]<-Kuhn_DR1k[,i*2+1]
}
DR_tableAll<-do.call(cbind,DR_table)
rownames(DR_tableAll)<-rownames(Kuhn_DR1k)	

table<-DR_tableAll

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

write.table(summary, paste0(dir,"DR-SUMMARY_Kuhn2011_resolved_1kTrees_5020species.txt"))


####
# Faurby and Svenning 2015 now.
##
FS2015_DR1k<-read.table(file=paste0(dir,"FS2015_resolved_tipCalcs_1kTrees_DR-ES_noPleistoExtinct.txt"), header=TRUE)

DR_table<-vector("list",length=1000)
for (i in 0:999){
	DR_table[[i+1]]<-FS2015_DR1k[,i*2+1]
}
DR_tableAll<-do.call(cbind,DR_table)
rownames(DR_tableAll)<-rownames(FS2015_DR1k)	

table<-DR_tableAll

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

write.table(summary, paste0(dir,"DR-SUMMARY_FS2015_resolved_1kTrees_5747species.txt"))



# ==================
# now PLOT to SUMMARIZE
###
# Plot the PAIRWISE COMPARISONS too...
setwd("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/DR_calcsOnOther_mamPhys")
source("DR_circleTree_functions.R")
library(RColorBrewer)

RESDIR<-"Fig6_compare_tipDRs/"

bbone<-"NDexp"
cladesDR<-read.table(paste0(dir,"MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt"))
head(cladesDR)

## Get the ACTUAL 95% CIs for all of these calculations:
####
library(data.table)
	# MAMPHY
	##
	all10kDR<-fread(file="/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/11_makingFullPosteriors/10k_tree_andDR_distributions_AugDec2018/DR-matrix_MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_all10k.txt")
	nRows<-5912
	# FS2015
	##
	all10kDR<-fread(file="/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/DR_calcsOnOther_mamPhys/FS2015_resolved_tipCalcs_1kTrees_DR-ES_noPleistoExtinct.txt")
	#nRows<-5747
		nRows<-5502

	# KUHN2011
	##
	all10kDR<-fread(file="/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/DR_calcsOnOther_mamPhys/Kuhn2011_resolved_tipCalcs_1kTrees_DR-ES.txt")
	nRows<-5020

	# Will want the mean, median, mode values across the ROWS for each tip in the tree >> 
	harmMeans<-data.frame(matrix(NA, nrow = nRows, ncol = 1), row.names=all10kDR$V1)
	medians<-data.frame(matrix(NA, nrow = nRows, ncol = 1), row.names=all10kDR$V1)
	means<-data.frame(matrix(NA, nrow = nRows, ncol = 1), row.names=all10kDR$V1)
	low95<-data.frame(matrix(NA, nrow = nRows, ncol = 1), row.names=all10kDR$V1)
	high95<-data.frame(matrix(NA, nrow = nRows, ncol = 1), row.names=all10kDR$V1)

	for (i in 1:length(rownames(all10kDR))){
		x <- as.numeric(unlist(all10kDR[i,])[(c(1:1000)*2)])
	#	x <- as.numeric(unlist(all10kDR[i,])[2:10001])
		harmMeans[i,] <- 1/(mean(1/x))
		medians[i,] <- median(x)
		means[i,] <- mean(x)
		low95[i,] <- as.numeric(quantile(x, 0.025))
		high95[i,] <- as.numeric(quantile(x, 0.975))
	}
	CIwidth<-high95-low95

	all10kDR_summary <- cbind(harmMeans,medians,means,low95, high95, CIwidth)#,deparse.level=1)
	colnames(all10kDR_summary) <- c("harmMeans","medians","means","low95","high95", "CIwidth")
	rownames(all10kDR_summary) <- all10kDR$V1
	tail(all10kDR_summary[order(all10kDR_summary$CIwidth),],100)

	#write.table(all10kDR_summary, paste0(dir,"DR-SUMMARY_MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_all10k_with95CI.txt"), col.names=TRUE, row.names=TRUE)
	write.table(all10kDR_summary, paste0(dir,"DR-SUMMARY_FS2015_5502sp_noPleistoExtinct_all1k_with95CI.txt"), col.names=TRUE, row.names=TRUE)
	#write.table(all10kDR_summary, paste0(dir,"DR-SUMMARY_KUHN2011_5020sp_all1k_with95CI.txt"), col.names=TRUE, row.names=TRUE)

	# Now JOIN those summary tables...
	###
	library(dplyr)
	setwd("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/DR_calcsOnOther_mamPhys")
	# load each:
	RESDIR<-"/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/DR_calcsOnOther_mamPhys/"
	mamPhy_DR<-read.table(file=paste0(RESDIR,"DR-SUMMARY_MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_all10k_with95CI.txt"), header=TRUE)
	FS2015_DR<-read.table(file=paste0(RESDIR,"DR-SUMMARY_FS2015_5502sp_noPleistoExtinct_all1k_with95CI.txt"), header=TRUE)
	K2011_DR<-read.table(file=paste0(RESDIR,"DR-SUMMARY_KUHN2011_5020sp_all1k_with95CI.txt"), header=TRUE)
	H2015_DR<-read.table(file=paste0(RESDIR,"Hedges_phySmooth_TTOL_tipCalcs_1tree_DR-ES.txt"), header=TRUE)
		# make first col the rownames
		mamPhy_DR<-cbind.data.frame(tiplabel=rownames(mamPhy_DR), mamPhy_DR)
		FS2015_DR<-cbind.data.frame(Species_Name=rownames(FS2015_DR),FS2015_DR)
		K2011_DR<-cbind.data.frame(Species_Name=rownames(K2011_DR),K2011_DR)
		H2015_DR<-cbind.data.frame(Species_Name=rownames(H2015_DR),harmMeans=H2015_DR$DR)

	# join MamPhy to taxonomy
	mamData<-read.csv(file="taxonomy_mamPhy_5911species_toPublish.csv",header=TRUE)
	mamPhy_DR_wTax<-left_join(x=mamData, y=mamPhy_DR, by="tiplabel")

	# join FS2015 to MamPhy
	step1 <- left_join(x=mamPhy_DR_wTax,y=FS2015_DR, by="Species_Name", suffix=c("_mamPhy","_FS"))
	# join KUHN2011 to MamPhy
	step2 <- left_join(x=mamPhy_DR_wTax,y=K2011_DR, by="Species_Name", suffix=c("_mamPhy","_KUHN"))
	# join H2015 to MamPhy
	step3 <- left_join(x=mamPhy_DR_wTax,y=H2015_DR, by="Species_Name", suffix=c("_mamPhy","_HEDGES"))

	#combine together
	compareTable<-cbind.data.frame(step1,step2[,18:23],harmMeans_HEDGES=step3[,18])
	write.table(compareTable,file=paste0("tipDR_compareTable_mamPhy-and-3supertrees.txt"))


# how many sp common across all trees?
noNAs<-na.omit(compareTable)
	# > 4403 species

length(na.omit(compareTable$harmMeans_KUHN)) # 4670
length(na.omit(compareTable$harmMeans_FS)) # 5329
length(na.omit(compareTable$harmMeans_HEDGES)) # 5033

sp_Overlaps<-c( length(na.omit(compareTable$Kuhn2011)), length(na.omit(compareTable$FS2015)), length(na.omit(compareTable$Hedges2015)) )

# how many sp in the top 1% of tip DR are common between each data set?
dataSets<-c("mamPhy","Kuhn2011","FS2015","Hedges2015")
dataSets_topPercentile<-vector("list",length(dataSets))
dataSets_topPercentile_gen<-vector("list",length(dataSets))
for(j in 1:length(dataSets)){
	
	dat<-compareTable_wGen[,c("tiplabel","gen", dataSets[j])]
	top1p<-quantile(dat[,3], 0.95, na.rm=TRUE)

	dataSets_topPercentile[[j]]<-dat[which(dat[,3] >= top1p),1]
	print(length(dat[which(dat[,3] >= top1p),1]))
	dataSets_topPercentile_gen[[j]]<-unique(dat[which(dat[,3] >= top1p),2])
	print( length( unique(dat[which(dat[,3] >= top1p),2])  ))

}

intersect(dataSets_topPercentile[[1]], dataSets_topPercentile[[2]])
intersect(dataSets_topPercentile[[1]], dataSets_topPercentile[[3]])
intersect(dataSets_topPercentile[[1]], dataSets_topPercentile[[4]])


dataSets_topPercentile_gen[[1]]
# [1] Crocidura      Ctenomys       Delphinus      Lepus          Petrogale     
# [6] Pteropus       Rattus         Sorex          Sousa          Spermophilus  
#[11] Stenella       Tamias         Trachypithecus
intersect(dataSets_topPercentile_gen[[1]], dataSets_topPercentile_gen[[2]])
#[1] "Spermophilus"   "Trachypithecus"
intersect(dataSets_topPercentile_gen[[1]], dataSets_topPercentile_gen[[3]])
#[1] "Pteropus" "Rattus"  
intersect(dataSets_topPercentile_gen[[1]], dataSets_topPercentile_gen[[4]])
#[1] "Spermophilus"   "Trachypithecus"

top_Overlaps_sp<-c( length(intersect(dataSets_topPercentile[[1]], dataSets_topPercentile[[2]])), length(intersect(dataSets_topPercentile[[1]], dataSets_topPercentile[[3]])), length(intersect(dataSets_topPercentile[[1]], dataSets_topPercentile[[4]])) )
top_Overlaps_gen<-c( length(intersect(dataSets_topPercentile_gen[[1]], dataSets_topPercentile_gen[[2]])), length(intersect(dataSets_topPercentile_gen[[1]], dataSets_topPercentile_gen[[3]])), length(intersect(dataSets_topPercentile_gen[[1]], dataSets_topPercentile_gen[[4]])) )

totalSp<-length(dataSets_topPercentile[[1]])
totalGen<-length(dataSets_topPercentile_gen[[1]])

top_Overlaps_sp/totalSp*100
# SPECIES-LEVEL
# top 1%: [1] 1.666667 6.666667 1.666667
# top 5%: [1] 11.14865 19.25676 15.87838
	#if na.omit: [1] 12.21719 20.81448 16.74208

top_Overlaps_gen/totalGen*100
# GENUS-LEVEL
# top 1%: [1] 15.38462 15.38462 15.38462
# top 5%: [1] 20.00000 27.27273 38.18182
	#if na.omit: [1] 21.73913 28.26087 41.30435

#[1] 296 sp		[1] 221 sp na.omit
#[1] 55 gen		[1] 46  gen na.omit
#[1] 234		[1] 221
#[1] 29			[1] 28
#[1] 267		[1] 221
#[1] 30			[1] 26
#[1] 253		[1] 221
#[1] 61			[1] 53
base<-dataSets_topPercentile_gen[[2]]
totalGen<-length(base)
top_Overlaps_gen<-c( length(intersect(base, dataSets_topPercentile_gen[[2]])), length(intersect(base, dataSets_topPercentile_gen[[3]])), length(intersect(base, dataSets_topPercentile_gen[[4]])) )

top_Overlaps_gen/totalGen*100
				#	BE2007	  #FS2015  #Hedges2015
# for BE2007: [1] 100.00000  34.48276  68.96552
# for FS2015: [1]  33.33333 100.00000  56.66667
# for He2015: [1]  32.78689  27.86885 100.00000
	#if.na.omit
				#	BE2007	  #FS2015  #Hedges2015
# for BE2007: [1]  100.00000  25.00000  67.85714
# for FS2015: [1]  26.92308 100.00000  57.69231
# for He2015: [1]  35.84906  28.30189 100.00000



base<-dataSets_topPercentile[[2]]
totalSp<-length(base)
top_Overlaps_sp<-c( length(intersect(base, dataSets_topPercentile[[2]])), length(intersect(base, dataSets_topPercentile[[3]])), length(intersect(base, dataSets_topPercentile[[4]])) )

top_Overlaps_sp/totalSp*100
				#	BE2007	  #FS2015  #Hedges2015
# for BE2007: [1] 100.000000   7.692308  35.294118
# for FS2015: [1]   7.692308 100.000000  14.932127
# for He2015: [1]  35.29412  14.93213 100.00000



# set up point colors
# by higherTaxon:
##
library(viridis)
cats<-names(sort(table(compareTable$higher),decreasing=FALSE))[c(4,3,2,5,6)]

#cols<-viridis(length(cats), alpha=0.7)
colsOrder<-rich.colors(length(cats)+1, alpha=1)
#cols<-c(colsOrder[6],colsOrder[4],colsOrder[2],colsOrder[5],colsOrder[3])
cols<-c(colsOrder[6],colsOrder[3],colsOrder[3],colsOrder[5],"blueviolet")
boxCols<-c(colsOrder[6],colsOrder[3],colsOrder[5],"blueviolet")

#	pointCols<-rep("black",length(noNAs[,1]))
#	
#	for(j in 1:length(cats)){
#	pointCols[which(as.vector(noNAs$higher)==cats[j])]<-cols[j]
#	}
#	
#	# by DR scale:
#	##
#	#cols<-rich.colors(100)
#	cols<-viridis(100)
#	
#	pointCols<-rep("black",length(noNAs[,1]))
#	
#		range <- quantile(noNAs$mamPhy, seq(0,1, 0.01))[2:101]
#	
#	    for (i in 1:100) {
#	        if (i==100) {range[i] <- range[i]+0.1}
#	        if (i==1) {pointCols[which(noNAs$mamPhy <= range[i])] <- cols[i] } 
#	        if (i > 1) {pointCols[which((noNAs$mamPhy <= range[i]) & (noNAs$mamPhy > range[i-1]))] <- cols[i]}
#	        }

# get data ready
datReady<-compareTable[,c("tiplabel","higher","harmMeans_mamPhy","low95_mamPhy","high95_mamPhy","harmMeans_FS","low95_FS","high95_FS","harmMeans_KUHN","low95_KUHN","high95_KUHN","harmMeans_HEDGES")]
	#maxDR<-max(dat[5:8])

RESDIR<-"Fig6_compare_tipDRs/"

# PLOT -- NEW all pairwise, but with CIs
####
png(file=paste0(RESDIR,"DR_pairwise_COMPARE_vs_Kuhn-FS-Hedges_viridis-by-higherTax_wSpearman_justTop3_wCI_greyCI_noPleistoExtinct.png"), height=4, width=11, units="in",res=600)
#pdf(file=paste0(RESDIR,"DR_pairwise_COMPARE_vs_Kuhn-FS-Hedges_viridis-by-higherTax_wSpearman_justTop3_wCI_greyCI_noPleistoExtinct.pdf"), height=4, width=11, onefile=TRUE)

#quartz(height=9,width=10)
par(oma = rep(3,4) + 0.1,
    mar = c(4,1.5,4,1) + 0.1) #‘c(bottom, left, top, right)’

layout(matrix(1:4,1,4, byrow=TRUE))#,widths=rep(0.25,4),heights=rep(,4))

# legend
	plot(harmMeans_mamPhy ~ harmMeans_FS, data=datReady, xlim=c(0,10),ylim=c(0,1), type="n", xlab="",ylab="", bty="n", xaxt="n", yaxt="n", cex.lab=1.2,cex.axis=1.3)
	legend(x=-1,y=0.95,legend=c( "Marsupialia","Afrotheria / Xenarthra", "Laurasiatheria", "Euarchontoglires"),fill=boxCols, cex=1.2, xpd="NA")
	par(new=TRUE)
	legend(x=-1,y=0.45,legend=c("Relative 1-to-1 line", "Regression line", "95% CI (this study)","95% CI (previous)"),lty=c(2,1,1,1), lwd=c(2,2,1,1),col=c("black","black","blue","grey"), cex=1.2, xpd="NA")

# mamPhy ~ Kuhn2011
	dat<-na.omit(datReady[,c(1:5,9:11)])
	yMax<-max(dat[,5])
	xMax<-max(dat[,8])
	form<-as.formula(harmMeans_mamPhy ~ harmMeans_KUHN)

	pointCols<-rep("black",length(dat[,1]))
	for(j in 1:length(cats)){
	pointCols[which(as.vector(dat$higher)==cats[j])]<-cols[j]
	}
	otherCI<-grey(0.5,alpha=0.1) #rgb(1,0,0,alpha=0.1)

	plot(form, data=dat, xlim=c(0,xMax),ylim=c(0,yMax), xlab="",ylab="",col=NA, cex.lab=1.2,cex.axis=1.3)
    plotCI(add=TRUE, y=dat[,3], x=dat[,6], err="y", li=dat[,4],ui=dat[,5], sfrac=0, scol=rgb(0,0,1,alpha=0.1), pch=NA)
    plotCI(add=TRUE, y=dat[,3], x=dat[,6], err="x", li=dat[,7],ui=dat[,8], sfrac=0, scol=otherCI, pch=NA)
 	points(form, data=dat, xlim=c(0,xMax),ylim=c(0,yMax), xlab="",ylab="",col=pointCols, cex.lab=1.2,cex.axis=1.3, pch=".")

	mtext(text="Kuhn et al. 2011",side=1,line=3, font=2)
	mtext(text="1,000 trees",side=1,line=4.5, font=1, cex=1)

	mtext(text="This study",side=2,line=4, font=2)
	mtext(text="10,000 trees",side=2,line=2.5, font=1, cex=1)

	model<-summary(lm(form, data=dat))
	a<-round(model$coef[[1]], 2)
	b<-round(model$coef[[2]], 2)
	r2<-round(model$r.squared,digits=2)
	X=seq(min(dat[,6]),max(dat[,6]),length.out=20)
	lines(x=X,y=b*X+a,col="black", lwd=3,lty=1) #x1=min(dat$MRCA),x2=max(dat$MRCA),y1=min(log(dat$richness)),y2=max(log(dat$richness)),untf=FALSE,

	abline(a=0, b=max(dat[,3])/max(dat[,6]), lty=2, lwd=2, col="black")
	#text(labels=bquote("y = "*.(a)*" + "*.(b)*"x; "*r^2*" = "*.(r2)), x=xMax-(xMax/2.5), y=0.02, col="black", cex=1.1)
	r<-round(cor(y=dat$harmMeans_mamPhy, x=dat[,6], method="spearman"),digits=2)
	mtext(text=bquote("spearman "*r*" = "*.(r) ), side=3, col="black", cex=0.8, line =1.5 )
	mtext(text=bquote("y = "*.(a)*" + "*.(b)*"x, "*R^2*" = "*.(r2) ), side=3, col="black", cex=0.8)
#dev.off()

# mamPhy ~ FS2015
	dat<-na.omit(datReady[,c(1:5,6:8)])

	xMax<-max(dat[,8])
	form<-as.formula(harmMeans_mamPhy ~ harmMeans_FS)
	pointCols<-rep("black",length(dat[,1]))
	for(j in 1:length(cats)){
	pointCols[which(as.vector(dat$higher)==cats[j])]<-cols[j]
	}

	plot(form, data=dat, xlim=c(0,xMax),ylim=c(0,yMax), xlab="",ylab="",col=NA, cex.lab=1.2,cex.axis=1.3, yaxt="n")
    plotCI(add=TRUE, y=dat[,3], x=dat[,6], err="y", li=dat[,4],ui=dat[,5], sfrac=0, scol=rgb(0,0,1,alpha=0.1), pch=NA)
    plotCI(add=TRUE, y=dat[,3], x=dat[,6], err="x", li=dat[,7],ui=dat[,8], sfrac=0, scol=otherCI, pch=NA)
 	points(form, data=dat, xlim=c(0,xMax),ylim=c(0,yMax), xlab="",ylab="",col=pointCols, cex.lab=1.2,cex.axis=1.3, pch=".")

	axis(side=2,labels=FALSE)
	mtext(text="Faurby and Svenning 2015",side=1,line=3, font=2)
	mtext(text="1,000 trees",side=1,line=4.5, font=1, cex=1)

	model<-summary(lm(form, data=dat))
	a<-round(model$coef[[1]], 2)
	b<-round(model$coef[[2]], 2)
	r2<-round(model$r.squared,digits=2)
	X=seq(min(dat[,6]),max(dat[,6]),length.out=20)
	lines(x=X,y=b*X+a,col="black", lwd=3,lty=1) #x1=min(dat$MRCA),x2=max(dat$MRCA),y1=min(log(dat$richness)),y2=max(log(dat$richness)),untf=FALSE,

	abline(a=0, b=max(dat[,3])/max(dat[,6]), lty=2, lwd=2, col="black")
	#text(labels=bquote("y = "*.(a)*" + "*.(b)*"x; "*r^2*" = "*.(r2)), x=xMax-(xMax/2.5), y=0.02, col="black", cex=1.1)
	r<-round(cor(y=dat$harmMeans_mamPhy, x=dat[,6], method="spearman"),digits=2)
	mtext(text=bquote("spearman "*r*" = "*.(r) ), side=3, col="black", cex=0.8, line =1.5 )
	mtext(text=bquote("y = "*.(a)*" + "*.(b)*"x, "*R^2*" = "*.(r2) ), side=3, col="black", cex=0.8)

# mamPhy ~ Hedges2015
	dat<-na.omit(datReady[,c(1:5,12)])

	xMax<-max(dat[,6])
	form<-as.formula(harmMeans_mamPhy ~ harmMeans_HEDGES)
	pointCols<-rep("black",length(dat[,1]))
	for(j in 1:length(cats)){
	pointCols[which(as.vector(dat$higher)==cats[j])]<-cols[j]
	}

	plot(form, data=dat, xlim=c(0,xMax),ylim=c(0,yMax), xlab="",ylab="",col=NA, cex.lab=1.2,cex.axis=1.3, yaxt="n")
    plotCI(add=TRUE, y=dat[,3], x=dat[,6], err="y", li=dat[,4],ui=dat[,5], sfrac=0, scol=rgb(0,0,1,alpha=0.1), pch=NA)
    #plotCI(add=TRUE, y=dat[,3], x=dat[,6], err="x", li=dat[,7],ui=dat[,8], sfrac=0, scol=rgb(1,0,0,alpha=0.1), pch=NA)
 	points(form, data=dat, xlim=c(0,xMax),ylim=c(0,yMax), xlab="",ylab="",col=pointCols, cex.lab=1.2,cex.axis=1.3, pch=".")

	axis(side=2,labels=FALSE)
	mtext(text="Hedges et al. 2015",side=1,line=3, font=2)
	mtext(text="1 tree",side=1,line=4.5, font=1, cex=1)

	model<-summary(lm(form, data=dat))
	a<-round(model$coef[[1]], 2)
	b<-round(model$coef[[2]], 2)
	r2<-round(model$r.squared,digits=2)
	X=seq(min(dat[,6]),max(dat[,6]),length.out=20)
	lines(x=X,y=b*X+a,col="black", lwd=3,lty=1) #x1=min(dat$MRCA),x2=max(dat$MRCA),y1=min(log(dat$richness)),y2=max(log(dat$richness)),untf=FALSE,

	abline(a=0, b=max(dat[,3])/max(dat[,6]), lty=2, lwd=2, col="black")
	#text(labels=bquote("y = "*.(a)*" + "*.(b)*"x; "*r^2*" = "*.(r2)), x=xMax-(xMax/2.5), y=0.02, col="black", cex=1.1)
	r<-round(cor(y=dat$harmMeans_mamPhy, x=dat[,6], method="spearman"),digits=2)
	mtext(text=bquote("spearman "*r*" = "*.(r) ), side=3, col="black", cex=0.8, line =1.5 )
	mtext(text=bquote("y = "*.(a)*" + "*.(b)*"x, "*R^2*" = "*.(r2) ), side=3, col="black", cex=0.8)


dev.off()




# PLOT -- [[old]] all pairwise, just points
####
pdf(file=paste0(RESDIR,"DR_pairwise_COMPARE_vs_Kuhn-FS-Hedges_viridis-by-higherTax_wSpearman_justTop3_wCI.pdf"), height=4, width=10, onefile=TRUE)

#quartz(height=9,width=10)
par(oma = rep(3,4) + 0.1,
    mar = c(4,1,4,1) + 0.1) #‘c(bottom, left, top, right)’

#layout(matrix(1:4,1,4, byrow=TRUE))#,widths=rep(0.25,4),heights=rep(,4))
layout(matrix(1:12,3,4, byrow=TRUE))#,widths=rep(0.25,4),heights=rep(,4))
#layout(matrix(1:4,1,4, byrow=FALSE),widths=rep(0.25,4),heights=rep(0.75,4))


# legend
	plot(mamPhy ~ Kuhn2011, data=dat, xlim=c(0,0),ylim=c(0,4000), type="n", xlab="",ylab="", bty="n", xaxt="n", yaxt="n", cex.lab=1.2,cex.axis=1.3)
	legend(x=0.1,y=0.6,legend=c("Monotremata", "Xenarthra", "Afrotheria","Marsupialia","Laurasiatheria", "Euarchontoglires"),fill=cols)

# mamPhy ~ Kuhn2011
	xMax<-max(dat[6])
	plot(mamPhy ~ Kuhn2011, data=dat, xlim=c(0,xMax),ylim=c(0,max(dat[5])), xlab="",ylab="",col=pointCols, cex.lab=1.2,cex.axis=1.3)
	#corr<-cor(y=dat$mamPhy, x=dat$Kuhn2011, method="spearman")
	mtext(text="Kuhn et al. 2011",side=1,line=3, font=2)
	mtext(text="This study",side=2,line=3, font=2)

	model<-summary(lm(mamPhy ~ Kuhn2011, data=dat))
	a<-round(model$coef[[1]], 2)
	b<-round(model$coef[[2]], 2)
	r2<-round(model$r.squared,digits=2)
	X=seq(min(dat$Kuhn2011),max(dat$Kuhn2011),length.out=20)
	lines(x=X,y=b*X+a,col="red", lwd=3,lty=1) #x1=min(dat$MRCA),x2=max(dat$MRCA),y1=min(log(dat$richness)),y2=max(log(dat$richness)),untf=FALSE,

	abline(a=0, b=max(dat[5])/xMax, lty=2, lwd=2, col="black")
	#text(labels=bquote("y = "*.(a)*" + "*.(b)*"x; "*r^2*" = "*.(r2)), x=xMax-(xMax/2.5), y=0.02, col="black", cex=1.1)
	r<-round(cor(y=dat$mamPhy, x=dat$Kuhn2011, method="spearman"),digits=2)
	mtext(text=bquote("spearman "*r*" = "*.(r) ), side=3, col="black", cex=0.8, line =1.5 )
	mtext(text=bquote("y = "*.(a)*" + "*.(b)*"x, "*R^2*" = "*.(r2) ), side=3, col="black", cex=0.8)

# mamPhy ~ FS2015
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
	r<-round(cor(y=dat$mamPhy, x=dat$FS2015, method="spearman"),digits=2)
	mtext(text=bquote("spearman "*r*" = "*.(r) ), side=3, col="black", cex=0.8, line =1.5 )
	mtext(text=bquote("y = "*.(a)*" + "*.(b)*"x, "*R^2*" = "*.(r2) ), side=3, col="black", cex=0.8)

# mamPhy ~ Hedges2015
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
	r<-round(cor(y=dat$mamPhy, x=dat$Hedges2015, method="spearman"),digits=2)
	mtext(text=bquote("spearman "*r*" = "*.(r) ), side=3, col="black", cex=0.8, line =1.5 )
	mtext(text=bquote("y = "*.(a)*" + "*.(b)*"x, "*R^2*" = "*.(r2) ), side=3, col="black", cex=0.8)


# fill 6
	plot(mamPhy ~ Kuhn2011, data=dat, xlim=c(0,xMax),ylim=c(0,max(dat[5])), type="n", xlab="",ylab="",col=pointCols, bty="n", xaxt="n", yaxt="n", cex.lab=1.2,cex.axis=1.3)
# fill 7
	plot(mamPhy ~ Kuhn2011, data=dat, xlim=c(0,xMax),ylim=c(0,max(dat[5])), type="n", xlab="",ylab="",col=pointCols, bty="n", xaxt="n", yaxt="n", cex.lab=1.2,cex.axis=1.3)

# Kuhn2011 ~ FS2015
	xMax<-max(dat[7])
	plot(Kuhn2011 ~ FS2015, data=dat, xlim=c(0,max(dat[7])),ylim=c(0,max(dat[6])),  xlab="",ylab="",col=pointCols, yaxt="n", cex.axis=1.3)
	axis(side=2,labels=FALSE)
	mtext(text="Faurby and Svenning 2015",side=1,line=3, font=2)
	mtext(text="Kuhn et al. 2011",side=2,line=3, font=2)

	model<-summary(lm(Kuhn2011 ~ FS2015, data=dat))
	a<-round(model$coef[[1]], 2)
	b<-round(model$coef[[2]], 2)
	r2<-round(model$r.squared,digits=2)
	X=seq(min(dat$FS2015),max(dat$FS2015),length.out=20)
	lines(x=X,y=b*X+a,col="red", lwd=3,lty=1) #x1=min(dat$MRCA),x2=max(dat$MRCA),y1=min(log(dat$richness)),y2=max(log(dat$richness)),untf=FALSE,

	abline(a=0, b=max(dat[6])/xMax, lty=2, lwd=2, col="black")
	#text(labels=bquote("y = "*.(a)*" + "*.(b)*"x; "*r^2*" = "*.(r2)), x=xMax-(xMax/2.5), y=0.02, col="black", cex=1.1)
	r<-round(cor(y=dat$Kuhn2011, x=dat$FS2015, method="spearman"),digits=2)
	mtext(text=bquote("spearman "*r*" = "*.(r) ), side=3, col="black", cex=0.8, line =1.5 )
	mtext(text=bquote("y = "*.(a)*" + "*.(b)*"x, "*R^2*" = "*.(r2) ), side=3, col="black", cex=0.8)

# Kuhn2011 ~ Hedges2015
	xMax<-max(dat[8])
	plot(Kuhn2011 ~ Hedges2015, data=dat, xlim=c(0,max(dat[8])),ylim=c(0,max(dat[6])),  xlab="",ylab="",col=pointCols, yaxt="n", cex.axis=1.3)
	axis(side=2,labels=FALSE)
	mtext(text="Hedges et al. 2015",side=1,line=3, font=2)
	#mtext(text="Kuhn et al. 2011",side=2,line=3, font=2)

	model<-summary(lm(Kuhn2011 ~ Hedges2015, data=dat))
	a<-round(model$coef[[1]], 2)
	b<-round(model$coef[[2]], 2)
	r2<-round(model$r.squared,digits=2)
	X=seq(min(dat$Hedges2015),max(dat$Hedges2015),length.out=20)
	lines(x=X,y=b*X+a,col="red", lwd=3,lty=1) #x1=min(dat$MRCA),x2=max(dat$MRCA),y1=min(log(dat$richness)),y2=max(log(dat$richness)),untf=FALSE,

	abline(a=0, b=max(dat[6])/xMax, lty=2, lwd=2, col="black")
	#text(labels=bquote("y = "*.(a)*" + "*.(b)*"x; "*r^2*" = "*.(r2)), x=xMax-(xMax/2.5), y=0.02, col="black", cex=1.1)
	r<-round(cor(y=dat$Kuhn2011, x=dat$Hedges2015, method="spearman"),digits=2)
	mtext(text=bquote("spearman "*r*" = "*.(r) ), side=3, col="black", cex=0.8, line =1.5 )
	mtext(text=bquote("y = "*.(a)*" + "*.(b)*"x, "*R^2*" = "*.(r2) ), side=3, col="black", cex=0.8)

# fill 9
	plot(mamPhy ~ Kuhn2011, data=dat, xlim=c(0,xMax),ylim=c(0,max(dat[5])), type="n", xlab="",ylab="",col=pointCols, bty="n", xaxt="n", yaxt="n", cex.lab=1.2,cex.axis=1.3)
# fill 10
	plot(mamPhy ~ Kuhn2011, data=dat, xlim=c(0,xMax),ylim=c(0,max(dat[5])), type="n", xlab="",ylab="",col=pointCols, bty="n", xaxt="n", yaxt="n", cex.lab=1.2,cex.axis=1.3)
# fill 11
	plot(mamPhy ~ Kuhn2011, data=dat, xlim=c(0,xMax),ylim=c(0,max(dat[5])), type="n", xlab="",ylab="",col=pointCols, bty="n", xaxt="n", yaxt="n", cex.lab=1.2,cex.axis=1.3)

# FS2015 ~ Hedges2015
	xMax<-max(dat[8])
	plot(FS2015 ~ Hedges2015, data=dat, xlim=c(0,max(dat[8])),ylim=c(0,max(dat[7])),  xlab="",ylab="",col=pointCols, yaxt="n", cex.axis=1.3)
	axis(side=2,labels=FALSE)
	mtext(text="Hedges et al. 2015",side=1,line=3, font=2)
	mtext(text="Faurby and Svenning 2015",side=2,line=3, font=2)

	model<-summary(lm(FS2015 ~ Hedges2015, data=dat))
	a<-round(model$coef[[1]], 2)
	b<-round(model$coef[[2]], 2)
	r2<-round(model$r.squared,digits=2)
	X=seq(min(dat$Hedges2015),max(dat$Hedges2015),length.out=20)
	lines(x=X,y=b*X+a,col="red", lwd=3,lty=1) #x1=min(dat$MRCA),x2=max(dat$MRCA),y1=min(log(dat$richness)),y2=max(log(dat$richness)),untf=FALSE,

	abline(a=0, b=max(dat[7])/xMax, lty=2, lwd=2, col="black")
	#text(labels=bquote("y = "*.(a)*" + "*.(b)*"x; "*r^2*" = "*.(r2)), x=xMax-(xMax/2.5), y=0.02, col="black", cex=1.1)
	r<-round(cor(y=dat$FS2015, x=dat$Hedges2015, method="spearman"),digits=2)
	mtext(text=bquote("spearman "*r*" = "*.(r) ), side=3, col="black", cex=0.8, line =1.5 )
	mtext(text=bquote("y = "*.(a)*" + "*.(b)*"x, "*R^2*" = "*.(r2) ), side=3, col="black", cex=0.8)


dev.off()



# NO TIPS -- LINEAR with DR
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")
setwd("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/DR_calcsOnOther_mamPhys")
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

pdf(file=paste0(RESDIR,"DR_on4phylosCompared_linear_richCol_justScale_ownColors.pdf"), width=8, height=6, onefile=TRUE)
#pdf(file=paste0(RESDIR,"DR_on4phylosCompared_linear_richCol_justScale_ownColors_withTips_80in.pdf"), width=6, height=80, onefile=TRUE)

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



########
# BAMMtools analyses...
###
library(BAMMtools)
library(coda)

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/_FBD_10trees")

folders<-c(1:10)
# load in tables
MSC_all10<-vector("list",length(folders))
for (i in 1:length(folders)){

MSC_i<-read.table(paste("results_MSC_",i,".txt",sep=""), header=TRUE)

MSC_all10[[i]]<-MSC_i[ order(MSC_i[,"node"]), ]

}

numShifts<-vector()
for (i in 1:10) {
	numShifts[i]<-nrow(MSC_all10[[i]])
}
#numShifts
# 26 28 32 27 24 28 20 29 32 24
mean(numShifts)
# 27
min(numShifts)
# 20
max(numShifts)
# 32

library(plyr)

allByNode_to1<-join_all(dfs=MSC_all10, by='node', type='full')

intersect(MSC_all10[[1]][,1],MSC_all10[[2]][,1])

reducedAll<-Reduce(function(...) merge(..., all=TRUE), MSC_all10)

tableAll<-table(reducedAll$node)

tableAll<-as.data.frame(tableAll[ order(tableAll, decreasing=TRUE) ]) # 205 nodes identified in at least one tree

tableAll_byNode<-as.data.frame(tableAll[order(as.numeric(rownames(tableAll))),])
names(tableAll_byNode)<-"num" 

tableAll_2orMore<-as.data.frame(tableAll[which(tableAll$num > 1),])
names(tableAll_2orMore)<-c("num") # 30 nodes ID'ed in 2 or more trees

      num
5912   10
11461   7
5918    6
8667    6
6043    4
7831    4
8288    4
9608    4
6279    3
6550    3
8523    3
8980    3
9341    3
11716   3
6044    2
6124    2
6173    2
6549    2
6728    2
6910    2
7089    2
7598    2
8185    2
8192    2
8524    2
8547    2
9773    2
9853    2
9970    2
10166   2

reducedAll[ with(reducedAll$node==names(tableAll_2orMore))]
reducedAll[ which(reducedAll$node %in% as.numeric(rownames(tableAll_2orMore))),]

	>>this was giving 




####
# Comparing BAMM and DR tip rates across runs...
###
# for FBD
folders<-c(1:10)
# load in tables
tipRates_all10<-vector("list",length(folders))
for (i in 1:length(folders)){

#setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/grace_2.5_1week_",folders[i],sep=""))
setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/grace_2.5_1week_NDexp_",folders[i],sep=""))

tipRates_i<-read.table(paste("meanTipRates_sample10_tree",i,".txt",sep=""), header=TRUE)

tipRates_all10[[i]]<-tipRates_i[ order(row.names(tipRates_i)), ]

}

#tipRates_all10[[1]][order(tipRates_all10[[1]][,"harmMeans_NetDiv"]),]
#tipMeans_FBD[ order(tipMeans_FBD[,"harmMeans_NetDiv"]), ]

# summarize across tables
globalMeans<-data.frame(matrix(NA, nrow = length(tipRates_all10[[1]][,1]), ncol = 8), row.names=rownames(tipRates_all10[[1]]))
colnames(globalMeans)<-c("harmMeans_Lam","arithMeans_Lam","harmMeans_Mu","arithMeans_Mu","harmMeans_NetDiv","arithMeans_NetDiv","harmMeans_NetDiv_sub","arithMeans_NetDiv_sub")

for (i in 1:length(tipRates_all10[[1]][,1])){
	harmLams<-vector()
	arithLams<-vector()
	harmMus<-vector()
	arithMus<-vector()
	harmNetDivs<-vector()
	arithNetDivs<-vector()
	harmNetDivs_sub<-vector()
	arithNetDivs_sub<-vector()

	for (j in 1:length(folders)){
		harmLams[j]<-tipRates_all10[[j]][i,"harmMeans_Lam"]
		arithLams[j]<-tipRates_all10[[j]][i,"arithMeans_Lam"]
		harmMus[j]<-tipRates_all10[[j]][i,"harmMeans_Mu"]
		arithMus[j]<-tipRates_all10[[j]][i,"arithMeans_Mu"]
		harmNetDivs[j]<-tipRates_all10[[j]][i,"harmMeans_NetDiv"]
		arithNetDivs[j]<-tipRates_all10[[j]][i,"arithMeans_NetDiv"]
		harmNetDivs_sub[j]<-tipRates_all10[[j]][i,"harmMeans_Lam"]-tipRates_all10[[j]][i,"harmMeans_Mu"]
		arithNetDivs_sub[j]<-tipRates_all10[[j]][i,"arithMeans_Lam"]-tipRates_all10[[j]][i,"arithMeans_Mu"]

		globalMeans[i,1]<-mean(harmLams)
		globalMeans[i,2]<-mean(arithLams)
		globalMeans[i,3]<-mean(harmMus)
		globalMeans[i,4]<-mean(arithMus)
		globalMeans[i,5]<-mean(harmNetDivs)
		globalMeans[i,6]<-mean(arithNetDivs)
		globalMeans[i,7]<-mean(harmNetDivs_sub)
		globalMeans[i,8]<-mean(arithNetDivs_sub)
}
}
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/")
#write.table(globalMeans, "globalMeanTipRates_sample10_tree1-10_FBD_wNetDiv.txt")
write.table(globalMeans, "globalMeanTipRates_sample10_tree1-10_NDexp_wNetDiv.txt")

# reload
tipMeans_FBD<-read.table("globalMeanTipRates_sample10_tree1-10_FBD_wNetDiv.txt")
tipMeans_ND<-read.table("globalMeanTipRates_sample10_tree1-10_NDexp_wNetDiv.txt")

# compare-- 
pdf(file="BAMM_tipMeanRates_FBD-ND_compare.pdf", width=5, height=5)
# within BAMM FBD
x="harmMeans_Lam"
y="arithMeans_Lam"
plot(tipMeans_FBD[,x], tipMeans_FBD[,y], main="within FBD", xlab=x, ylab=y)

x="harmMeans_Mu"
y="arithMeans_Mu"
plot(tipMeans_FBD[,x], tipMeans_FBD[,y], main="within FBD", xlab=x, ylab=y)

x="harmMeans_NetDiv_sub"
y="arithMeans_NetDiv_sub"
plot(tipMeans_FBD[,x], tipMeans_FBD[,y], main="within FBD", xlab=x, ylab=y)


# within NDexp
x="harmMeans_Lam"
y="arithMeans_Lam"
plot(tipMeans_ND[,x], tipMeans_ND[,y], main="within NDexp", xlab=x, ylab=y)

x="harmMeans_Mu"
y="arithMeans_Mu"
plot(tipMeans_ND[,x], tipMeans_ND[,y], main="within NDexp", xlab=x, ylab=y)

x="harmMeans_NetDiv_sub"
y="arithMeans_NetDiv_sub"
plot(tipMeans_ND[,x], tipMeans_ND[,y], main="within NDexp", xlab=x, ylab=y)


# between FBD and ND
x="FBD"
y="ND"
var="harmMeans_Lam"
plot(tipMeans_FBD[,var], tipMeans_ND[,var], main="harmMeans_Lam (FBD vs NDexp)", xlab=x, ylab=y)

x="FBD"
y="ND"
var="arithMeans_Lam"
plot(tipMeans_FBD[,var], tipMeans_ND[,var], main="arithMeans_Lam (FBD vs NDexp)", xlab=x, ylab=y)

x="FBD"
y="ND"
var="harmMeans_Mu"
plot(tipMeans_FBD[,var], tipMeans_ND[,var], main="harmMeans_Mu (FBD vs NDexp)", xlab=x, ylab=y)

x="FBD"
y="ND"
var="arithMeans_Mu"
plot(tipMeans_FBD[,var], tipMeans_ND[,var], main="arithMeans_Mu (FBD vs NDexp)", xlab=x, ylab=y)

x="FBD"
y="ND"
var="harmMeans_NetDiv_sub"
plot(tipMeans_FBD[,var], tipMeans_ND[,var], main="harmMeans_NetDiv_sub (FBD vs NDexp)", xlab=x, ylab=y)

x="FBD"
y="ND"
var="arithMeans_NetDiv_sub"
plot(tipMeans_FBD[,var], tipMeans_ND[,var], main="arithMeans_NetDiv_sub (FBD vs NDexp)", xlab=x, ylab=y)

dev.off()

# now ONLY the HARMONIC mean comparisons
# within
pdf(file="BAMM_tipMeanRates_withinFBD-ND_harmMean.pdf", width=5, height=5)
# within BAMM FBD, speciation vs extinction
x="harmMeans_Lam"
y="harmMeans_Mu"
plot(tipMeans_FBD[,x], tipMeans_FBD[,y], main="within FBD", xlab=x, ylab=y)

x="harmMeans_Lam"
y="harmMeans_NetDiv_sub"
plot(tipMeans_FBD[,x], tipMeans_FBD[,y], main="within FBD", xlab=x, ylab=y)

x="harmMeans_Mu"
y="harmMeans_NetDiv_sub"
plot(tipMeans_FBD[,x], tipMeans_FBD[,y], main="within FBD", xlab=x, ylab=y)


# within NDexp
x="harmMeans_Lam"
y="harmMeans_Mu"
plot(tipMeans_ND[,x], tipMeans_ND[,y], main="within NDexp", xlab=x, ylab=y)

x="harmMeans_Lam"
y="harmMeans_NetDiv_sub"
plot(tipMeans_ND[,x], tipMeans_ND[,y], main="within NDexp", xlab=x, ylab=y)

x="harmMeans_Mu"
y="harmMeans_NetDiv_sub"
plot(tipMeans_ND[,x], tipMeans_ND[,y], main="within NDexp", xlab=x, ylab=y)

dev.off()

# now ONLY the HARMONIC mean comparisons
# between FBD and ND
pdf(file="BAMM_tipMeanRates_btwnFBD-ND_harmMean.pdf", width=5, height=5)

x="FBD"
y="ND"
var="harmMeans_Lam"
plot(tipMeans_FBD[,var], tipMeans_ND[,var], main="harmMeans_Lam (FBD vs NDexp)", xlab=x, ylab=y)

x="FBD"
y="ND"
var="harmMeans_Mu"
plot(tipMeans_FBD[,var], tipMeans_ND[,var], main="harmMeans_Mu (FBD vs NDexp)", xlab=x, ylab=y)

x="FBD"
y="ND"
var="harmMeans_NetDiv_sub"
plot(tipMeans_FBD[,var], tipMeans_ND[,var], main="harmMeans_NetDiv_sub (FBD vs NDexp)", xlab=x, ylab=y)

dev.off()


# load in the TIP-LEVEL DR data...
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")

cladesDR_FBD<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_DR-SUMMARY-expanded_FBD.txt")
head(cladesDR_FBD)
colnames(cladesDR_FBD)<-c("tiplabel","gen","fam","ord","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")

cladesDR_ND<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_DR-SUMMARY-expanded_NDexp.txt")
head(cladesDR_ND)
colnames(cladesDR_ND)<-c("tiplabel","gen","fam","ord","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")

DR_FBD<-cladesDR_FBD$harmMeans
names(DR_FBD)<-cladesDR_FBD$tiplabel
FBD_DR_BAMM<-cbind(DR_FBD, tipMeans_FBD[which(rownames(tipMeans_FBD)==names(DR_FBD)),"harmMeans_Lam"], tipMeans_FBD[which(rownames(tipMeans_FBD)==names(DR_FBD)),"harmMeans_NetDiv_sub"], tipMeans_FBD[which(rownames(tipMeans_FBD)==names(DR_FBD)),"harmMeans_Mu"])
colnames(FBD_DR_BAMM)<-c("DR", "Lam", "netDiv", "Mu")

DR_ND<-cladesDR_ND$harmMeans
names(DR_ND)<-cladesDR_ND$tiplabel
ND_DR_BAMM<-cbind(DR_ND, tipMeans_ND[which(rownames(tipMeans_ND)==names(DR_ND)),"harmMeans_Lam"], tipMeans_ND[which(rownames(tipMeans_ND)==names(DR_ND)),"harmMeans_NetDiv_sub"], tipMeans_ND[which(rownames(tipMeans_ND)==names(DR_ND)),"harmMeans_Mu"])
colnames(ND_DR_BAMM)<-c("DR", "Lam", "netDiv", "Mu")

# Look at DR > or < than BAMM LAMBDA
##
subtract=FBD_DR_BAMM[,"Lam"] - FBD_DR_BAMM[,"DR"]

LamLessDR<-names(subtract[which(subtract < 0)]) # 2291 species
LamGreaterDR<-names(subtract[which(subtract > 0)]) # 3620 species

plot(x=FBD_DR_BAMM[ which(rownames(FBD_DR_BAMM) %in% LamGreaterDR),"DR"],y=FBD_DR_BAMM[which(rownames(FBD_DR_BAMM) %in% LamGreaterDR),"Lam"], xlim=c(0,1),ylim=c(0,1))
plot(x=FBD_DR_BAMM[ which(rownames(FBD_DR_BAMM) %in% LamLessDR),"DR"],y=FBD_DR_BAMM[which(rownames(FBD_DR_BAMM) %in% LamLessDR),"Lam"], xlim=c(0,1),ylim=c(0,1))

tree<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target.tre")
tree2 <- drop.tip(tree,"_Anolis_carolinensis")
tree3 <- ladderize(tree2, right=TRUE)

tipColors <- rep("black",5911)
tipColors[match(LamLessDR,tree3$tip.label)] <- "red" 
missing<-read.table("MamPhy_FIN4_1813sp_missing_LIST.txt", header=FALSE)

missingAndLamGrThanDR<-setdiff(missing$V1,as.vector(LamGreaterDR)) #  728 sp of 1813
missingAndLamLsThanDR<-setdiff(missing$V1,as.vector(LamLessDR)) #  1085 sp of 1813



setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/")

pdf(file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_5911taxa_radial_BAMM_bdLam_2291sp_lessThanDR-OK.pdf", width=30, height=30, onefile=TRUE)

plot(tree3, cex=0.15, tip.color=tipColors, type="fan", label.offset=0.4)

dev.off()



pdf(file="BAMM_tipMeanRates_FBD-ND_compareTo_DR_withDRGrLessLam_DRresponse.pdf", width=10, height=5, onefile=TRUE)

pdf(file="BAMM_tipMeanRates_FBD-ND_compareTo_DRresponse.pdf", width=10, height=5, onefile=TRUE)

pdf(file="BAMM_tipMeanRates_FBD-ND_compareTo_DRresponse_smoothScatterBl.pdf", width=10, height=5, onefile=TRUE)
margins= c(5, 4, 2, 2)

layout(matrix(c(1:2), 1, 2, byrow = TRUE))

# DR vs BAMM 
##
# LAMBDA
x="Tip-level BAMM"
y="Tip-level DR"
form=FBD_DR_BAMM[,"DR"] ~ FBD_DR_BAMM[,"Lam"]
#hexbinplot(form, cex.lab=1.3, cex.main=1.5, xlim=c(0,1.2), ylim=c(0,1.2), main="DR ~ BAMM speciation rate", xlab=x, ylab=y)#, mar=margins)
smoothScatter(form, nrpoints=0, colramp = colorRampPalette(c("white", blues9)), cex.lab=1.3, cex.main=1.5, xlim=c(0,1.1), ylim=c(0,1.1), main="DR ~ BAMM speciation rate", xlab=x, ylab=y, mar=margins)
#plot(form, col="#00000033", cex.lab=1.3, cex.main=1.5, xlim=c(0,1), ylim=c(0,1), main="DR ~ BAMM speciation rate", xlab=x, ylab=y, mar=margins)
dd<-lm(form)
abline(coef(dd)[1],coef(dd)[2])
ss<-summary(dd)
mtext(bquote( plain("FBD backbone; ")* plain("R") ^ plain("2") == .(round(ss$r.squared, digits=3)) * plain("; Y = ") * .(round(coef(dd)[1], digits=3)) * plain(" + ") * .(round(coef(dd)[2], digits=3)) ~ plain("X")))

x="Tip-level BAMM"
y=""
form=ND_DR_BAMM[,"DR"] ~ ND_DR_BAMM[,"Lam"]
#plot(form, cex.lab=1.4, xlab=x, ylab=y, xlim=c(0,1), ylim=c(0,1), mar=margins)
smoothScatter(form, nrpoints=0, colramp = colorRampPalette(c("white", blues9)), cex.lab=1.3, cex.main=1.5, xlim=c(0,1.1), ylim=c(0,1.1), main="", xlab=x, ylab=y, mar=margins)
dd<-lm(form)
abline(coef(dd)[1],coef(dd)[2])
ss<-summary(dd)
mtext(bquote( plain("NDexp backbone; ")* plain("R") ^ plain("2") == .(round(ss$r.squared, digits=3)) * plain("; Y = ") * .(round(coef(dd)[1], digits=3)) * plain(" + ") * .(round(coef(dd)[2], digits=3)) ~ plain("X")))

subtract=FBD_DR_BAMM[,"Lam"] - FBD_DR_BAMM[,"DR"]
LamGreaterDR<-names(subtract[which(subtract < 0)]) # 2291 species
LamLessDR<-names(subtract[which(subtract > 0)]) # 3620 species

subtract2=ND_DR_BAMM[,"Lam"] - ND_DR_BAMM[,"DR"]
LamGreaterDR2<-names(subtract[which(subtract2 < 0)]) # 2253 species
LamLessDR2<-names(subtract[which(subtract2 > 0)]) # 3658 species

## DIFF > 0
#layout(matrix(c(1:2), 1, 2, byrow = TRUE))
#
#x="Tip-level BAMM"
#y="Tip-level DR"
#form=FBD_DR_BAMM[which(rownames(FBD_DR_BAMM) %in% LamLessDR),"DR"] ~ FBD_DR_BAMM[ which(rownames(FBD_DR_BAMM) %in% LamLessDR),"Lam"]
#
#plot(form, cex.lab=1.3, cex.main=1.5, xlim=c(0,1), ylim=c(0,1), main="BD Lam > yule DR", xlab=x, ylab=y, mar=margins)
#dd<-lm(form)
#abline(coef(dd)[1],coef(dd)[2])
#ss<-summary(dd)
#mtext(bquote( plain("FBD backbone; ")* plain("R") ^ plain("2") == .(round(ss$r.squared, digits=3)) * plain("; Y = ") * .(round(coef(dd)[1], digits=3)) * plain(" + ") * .(round(coef(dd)[2], digits=3)) ~ plain("X")))
#
#x="Tip-level DR"
#y=""
#form=ND_DR_BAMM[which(rownames(ND_DR_BAMM) %in% LamLessDR2),"DR"] ~ ND_DR_BAMM[ which(rownames(ND_DR_BAMM) %in% LamLessDR2),"Lam"]
#
#plot(form, cex.lab=1.4, xlab=x, ylab=y, xlim=c(0,1), ylim=c(0,1), mar=margins)
#dd<-lm(form)
#abline(coef(dd)[1],coef(dd)[2])
#ss<-summary(dd)
#mtext(bquote( plain("NDexp backbone; ")* plain("R") ^ plain("2") == .(round(ss$r.squared, digits=3)) * plain("; Y = ") * .(round(coef(dd)[1], digits=3)) * plain(" + ") * .(round(coef(dd)[2], digits=3)) ~ plain("X")))
#
## DIFF < 0
#layout(matrix(c(1:2), 1, 2, byrow = TRUE))
#
#x="Tip-level BAMM"
#y="Tip-level DR"
#form=FBD_DR_BAMM[which(rownames(FBD_DR_BAMM) %in% LamGreaterDR),"DR"] ~ FBD_DR_BAMM[ which(rownames(FBD_DR_BAMM) %in% LamGreaterDR),"Lam"]
#
#plot(form, cex.lab=1.3, cex.main=1.5, xlim=c(0,1), ylim=c(0,1), main="BD Lam < yule DR", xlab=x, ylab=y, mar=margins)
#dd<-lm(form)
#abline(coef(dd)[1],coef(dd)[2])
#ss<-summary(dd)
#mtext(bquote( plain("FBD backbone; ")* plain("R") ^ plain("2") == .(round(ss$r.squared, digits=3)) * plain("; Y = ") * .(round(coef(dd)[1], digits=3)) * plain(" + ") * .(round(coef(dd)[2], digits=3)) ~ plain("X")))
#
#x="Tip-level BAMM"
#y=""
#form=ND_DR_BAMM[which(rownames(ND_DR_BAMM) %in% LamGreaterDR2),"Lam"] ~ ND_DR_BAMM[ which(rownames(ND_DR_BAMM) %in% LamGreaterDR2),"DR"]
#
#plot(form, cex.lab=1.4, xlab=x, ylab=y, xlim=c(0,1), ylim=c(0,1), mar=margins)
#dd<-lm(form)
#abline(coef(dd)[1],coef(dd)[2])
#ss<-summary(dd)
#mtext(bquote( plain("NDexp backbone; ")* plain("R") ^ plain("2") == .(round(ss$r.squared, digits=3)) * plain("; Y = ") * .(round(coef(dd)[1], digits=3)) * plain(" + ") * .(round(coef(dd)[2], digits=3)) ~ plain("X")))

# NETDIV
##
layout(matrix(c(1:2), 1, 2, byrow = TRUE))

x="Tip-level BAMM"
y="Tip-level DR"
form=FBD_DR_BAMM[,"DR"] ~ FBD_DR_BAMM[,"netDiv"]
#plot(form, cex.lab=1.4, cex.main=1.5, xlim=c(0,1), ylim=c(0,1), main="DR ~ BAMM net div rate", xlab=x, ylab=y, mar=margins)
smoothScatter(form, nrpoints=0, colramp = colorRampPalette(c("white", blues9)), cex.lab=1.3, cex.main=1.5, xlim=c(0,1.1), ylim=c(0,1.1), main="DR ~ BAMM net div rate", xlab=x, ylab=y, mar=margins)
dd<-lm(form)
abline(coef(dd)[1],coef(dd)[2])
ss<-summary(dd)
mtext(bquote( plain("FBD backbone; ")* plain("R") ^ plain("2") == .(round(ss$r.squared, digits=3)) * plain("; Y = ") * .(round(coef(dd)[1], digits=3)) * plain(" + ") * .(round(coef(dd)[2], digits=3)) ~ plain("X")))

x="Tip-level BAMM"
y=""
form=ND_DR_BAMM[,"DR"] ~ ND_DR_BAMM[,"netDiv"]
#plot(form, cex.lab=1.4, xlab=x, ylab=y, xlim=c(0,1), ylim=c(0,1), mar=margins)
smoothScatter(form, nrpoints=0, colramp = colorRampPalette(c("white", blues9)), cex.lab=1.3, cex.main=1.5, xlim=c(0,1.1), ylim=c(0,1.1), main="", xlab=x, ylab=y, mar=margins)
dd<-lm(form)
abline(coef(dd)[1],coef(dd)[2])
ss<-summary(dd)
mtext(bquote( plain("NDexp backbone; ")* plain("R") ^ plain("2") == .(round(ss$r.squared, digits=3)) * plain("; Y = ") * .(round(coef(dd)[1], digits=3)) * plain(" + ") * .(round(coef(dd)[2], digits=3)) ~ plain("X")))

## DIFF > 0
#layout(matrix(c(1:2), 1, 2, byrow = TRUE))
#
#x="Tip-level BAMM"
#y="Tip-level DR"
#form=FBD_DR_BAMM[which(rownames(FBD_DR_BAMM) %in% LamLessDR),"DR"] ~ FBD_DR_BAMM[ which(rownames(FBD_DR_BAMM) %in% LamLessDR),"netDiv"]
#
#plot(form, cex.lab=1.3, cex.main=1.5, xlim=c(0,1), ylim=c(0,1), main="BD Lam > yule DR", xlab=x, ylab=y, mar=margins)
#dd<-lm(form)
#abline(coef(dd)[1],coef(dd)[2])
#ss<-summary(dd)
#mtext(bquote( plain("FBD backbone; ")* plain("R") ^ plain("2") == .(round(ss$r.squared, digits=3)) * plain("; Y = ") * .(round(coef(dd)[1], digits=3)) * plain(" + ") * .(round(coef(dd)[2], digits=3)) ~ plain("X")))
#
#x="Tip-level BAMM"
#y=""
#form=ND_DR_BAMM[which(rownames(ND_DR_BAMM) %in% LamLessDR2),"DR"] ~ ND_DR_BAMM[ which(rownames(ND_DR_BAMM) %in% LamLessDR2),"netDiv"]
#
#plot(form, cex.lab=1.4, xlab=x, ylab=y, xlim=c(0,1), ylim=c(0,1), mar=margins)
#dd<-lm(form)
#abline(coef(dd)[1],coef(dd)[2])
#ss<-summary(dd)
#mtext(bquote( plain("NDexp backbone; ")* plain("R") ^ plain("2") == .(round(ss$r.squared, digits=3)) * plain("; Y = ") * .(round(coef(dd)[1], digits=3)) * plain(" + ") * .(round(coef(dd)[2], digits=3)) ~ plain("X")))
#
## DIFF < 0
#layout(matrix(c(1:2), 1, 2, byrow = TRUE))
#
#x="Tip-level BAMM"
#y="Tip-level DR"
#form=FBD_DR_BAMM[which(rownames(FBD_DR_BAMM) %in% LamGreaterDR),"DR"] ~ FBD_DR_BAMM[ which(rownames(FBD_DR_BAMM) %in% LamGreaterDR),"netDiv"]
#
#plot(form, cex.lab=1.3, cex.main=1.5, xlim=c(0,1), ylim=c(0,1), main="BD Lam < yule DR", xlab=x, ylab=y, mar=margins)
#dd<-lm(form)
#abline(coef(dd)[1],coef(dd)[2])
#ss<-summary(dd)
#mtext(bquote( plain("FBD backbone; ")* plain("R") ^ plain("2") == .(round(ss$r.squared, digits=3)) * plain("; Y = ") * .(round(coef(dd)[1], digits=3)) * plain(" + ") * .(round(coef(dd)[2], digits=3)) ~ plain("X")))
#
#x="Tip-level DR"
#y=""
#form=ND_DR_BAMM[which(rownames(ND_DR_BAMM) %in% LamGreaterDR2),"DR"] ~ ND_DR_BAMM[ which(rownames(ND_DR_BAMM) %in% LamGreaterDR2),"netDiv"]
#
#plot(form, cex.lab=1.4, xlab=x, ylab=y, xlim=c(0,1), ylim=c(0,1), mar=margins)
#dd<-lm(form)
#abline(coef(dd)[1],coef(dd)[2])
#ss<-summary(dd)
#mtext(bquote( plain("NDexp backbone; ")* plain("R") ^ plain("2") == .(round(ss$r.squared, digits=3)) * plain("; Y = ") * .(round(coef(dd)[1], digits=3)) * plain(" + ") * .(round(coef(dd)[2], digits=3)) ~ plain("X")))

# MU
##
layout(matrix(c(1:2), 1, 2, byrow = TRUE))

x="Tip-level BAMM"
y="Tip-level DR"
form=FBD_DR_BAMM[,"DR"] ~ FBD_DR_BAMM[,"Mu"]
#plot(form, cex.lab=1.4, cex.main=1.5, xlim=c(0,1), ylim=c(0,1), main="DR ~ BAMM extinction rate", xlab=x, ylab=y, mar=margins)
smoothScatter(form, nrpoints=0, colramp = colorRampPalette(c("white", blues9)), cex.lab=1.3, cex.main=1.5, xlim=c(0,1.1), ylim=c(0,1.1), main="DR ~ BAMM extinction rate", xlab=x, ylab=y, mar=margins)
dd<-lm(form)
abline(coef(dd)[1],coef(dd)[2])
ss<-summary(dd)
mtext(bquote( plain("FBD backbone; ")* plain("R") ^ plain("2") == .(round(ss$r.squared, digits=3)) * plain("; Y = ") * .(round(coef(dd)[1], digits=3)) * plain(" + ") * .(round(coef(dd)[2], digits=3)) ~ plain("X")))

x="Tip-level BAMM"
y=""
form=ND_DR_BAMM[,"DR"] ~ ND_DR_BAMM[,"Mu"]
#plot(form, cex.lab=1.4, xlab=x, ylab=y, xlim=c(0,1), ylim=c(0,1), mar=margins)
smoothScatter(form, nrpoints=0, colramp = colorRampPalette(c("white", blues9)), cex.lab=1.3, cex.main=1.5, xlim=c(0,1.1), ylim=c(0,1.1), main="", xlab=x, ylab=y, mar=margins)
dd<-lm(form)
abline(coef(dd)[1],coef(dd)[2])
ss<-summary(dd)
mtext(bquote( plain("NDexp backbone; ")* plain("R") ^ plain("2") == .(round(ss$r.squared, digits=3)) * plain("; Y = ") * .(round(coef(dd)[1], digits=3)) * plain(" + ") * .(round(coef(dd)[2], digits=3)) ~ plain("X")))

## DIFF > 0
#layout(matrix(c(1:2), 1, 2, byrow = TRUE))
#
#x="Tip-level BAMM"
#y="Tip-level DR"
#form=FBD_DR_BAMM[which(rownames(FBD_DR_BAMM) %in% LamLessDR),"DR"] ~ FBD_DR_BAMM[ which(rownames(FBD_DR_BAMM) %in% LamLessDR),"Mu"]
#
#plot(form, cex.lab=1.3, cex.main=1.5, xlim=c(0,1), ylim=c(0,1), main="BD Lam > yule DR", xlab=x, ylab=y, mar=margins)
#dd<-lm(form)
#abline(coef(dd)[1],coef(dd)[2])
#ss<-summary(dd)
#mtext(bquote( plain("FBD backbone; ")* plain("R") ^ plain("2") == .(round(ss$r.squared, digits=3)) * plain("; Y = ") * .(round(coef(dd)[1], digits=3)) * plain(" + ") * .(round(coef(dd)[2], digits=3)) ~ plain("X")))
#
#x="Tip-level BAMM"
#y=""
#form=ND_DR_BAMM[which(rownames(ND_DR_BAMM) %in% LamLessDR2),"DR"] ~ ND_DR_BAMM[ which(rownames(ND_DR_BAMM) %in% LamLessDR2),"Mu"]
#
#plot(form, cex.lab=1.4, xlab=x, ylab=y, xlim=c(0,1), ylim=c(0,1), mar=margins)
#dd<-lm(form)
#abline(coef(dd)[1],coef(dd)[2])
#ss<-summary(dd)
#mtext(bquote( plain("NDexp backbone; ")* plain("R") ^ plain("2") == .(round(ss$r.squared, digits=3)) * plain("; Y = ") * .(round(coef(dd)[1], digits=3)) * plain(" + ") * .(round(coef(dd)[2], digits=3)) ~ plain("X")))
#
## DIFF < 0
#layout(matrix(c(1:2), 1, 2, byrow = TRUE))
#
#x="Tip-level BAMM"
#y="Tip-level DR"
#form=FBD_DR_BAMM[which(rownames(FBD_DR_BAMM) %in% LamGreaterDR),"DR"] ~ FBD_DR_BAMM[ which(rownames(FBD_DR_BAMM) %in% LamGreaterDR),"Mu"]
#
#plot(form, cex.lab=1.3, cex.main=1.5, xlim=c(0,1), ylim=c(0,1), main="BD Lam < yule DR", xlab=x, ylab=y, mar=margins)
#dd<-lm(form)
#abline(coef(dd)[1],coef(dd)[2])
#ss<-summary(dd)
#mtext(bquote( plain("FBD backbone; ")* plain("R") ^ plain("2") == .(round(ss$r.squared, digits=3)) * plain("; Y = ") * .(round(coef(dd)[1], digits=3)) * plain(" + ") * .(round(coef(dd)[2], digits=3)) ~ plain("X")))
#
#x="Tip-level BAMM"
#y=""
#form=ND_DR_BAMM[which(rownames(ND_DR_BAMM) %in% LamGreaterDR2),"DR"] ~ ND_DR_BAMM[ which(rownames(ND_DR_BAMM) %in% LamGreaterDR2),"Mu"]
#
#plot(form, cex.lab=1.4, xlab=x, ylab=y, xlim=c(0,1), ylim=c(0,1), mar=margins)
#dd<-lm(form)
#abline(coef(dd)[1],coef(dd)[2])
#ss<-summary(dd)
#mtext(bquote( plain("NDexp backbone; ")* plain("R") ^ plain("2") == .(round(ss$r.squared, digits=3)) * plain("; Y = ") * .(round(coef(dd)[1], digits=3)) * plain(" + ") * .(round(coef(dd)[2], digits=3)) ~ plain("X")))
#
dev.off()




#setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/grace_2.5_1week_",folders[i],sep=""))
setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/grace_2.5_1week_NDexp_",folders[i],sep=""))

#tree <- read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample10_",folders[i],".tre",sep=""))
#edata <- getEventData(tree, eventdata = "mamPhy_FBD_event_data.txt", burnin=0.33, nsamples=1000)

assign(paste("tree_ND_",folders[i],sep=""), read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample10_",folders[i],".tre",sep="")))

assign(paste("edata_ND_",folders[i],sep=""), getEventData(tree, eventdata = "mamPhy_NDexp_event_data.txt", burnin=0.33, nsamples=1000))


shift_probs <- summary(edata)

write.table(shift_probs, file="results.txt", append=TRUE)

bfmat <- computeBayesFactors(postburn, expectedNumberOfShifts=1, burnin=0.0) #burnin already done

write.table(bfmat, file="results.txt", append=TRUE)

pdf(file="prior-posterior_compare.pdf")
plotPrior(postburn, expectedNumberOfShifts = 1, burnin = 0.0)
dev.off()

# phylorate

pdf(file="mean-phylorate.pdf", width=8.5, height=200)
plot.bammdata(edata, lwd=2, legend=TRUE, labels=TRUE,cex=0.2)
dev.off()

pdf(file="mean-phylorate_polar.pdf", width=35, height=35)
plot.bammdata(edata, lwd=2, legend=TRUE, method="polar", labels=TRUE,cex=0.11)
dev.off()


pdf(file="mean-phylorate_25th_wShifts.pdf")
index <- 25
e2 <- subsetEventData(edata, index = index)
plot.bammdata(e2, lwd=2, legend=TRUE)
addBAMMshifts(e2, cex=2)
dev.off()

# Credible shifts
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)

write.table(summary(css), file="results.txt", append=TRUE)

pdf(file="credibleShiftSet.pdf")
plot.credibleshiftset(css)
dev.off()


best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1)
pdf(file="bestShiftSet.pdf")
plot.bammdata(best, lwd = 2)
addBAMMshifts(best, cex=2.5)
dev.off()

# maximum shift sets...

msc.set <- maximumShiftCredibility(edata, maximize='product')

msc.config <- subsetEventData(edata, index = msc.set$sampleindex)

pdf(file="MAX_ShiftSet.pdf")
plot.bammdata(msc.config, lwd=2)
addBAMMshifts(msc.config, cex = 2)
dev.off()

pdf(file="MAX_ShiftSet_largePolar.pdf", width=35, height=35)
plot.bammdata(msc.config, lwd=2, legend=TRUE, method="polar", labels=TRUE,cex=0.11)
addBAMMshifts(msc.config, cex = 2, method="polar")
dev.off()

write.table(msc.config$eventData, file="results.txt", append=TRUE)

dsc <- distinctShiftConfigurations(edata, expectedNumberOfShifts=1, threshold=5)

# Here is one random sample with the BEST shift configuration
pdf(file="random_bestShift1.pdf")
plot.bammshifts(dsc, edata, rank=1, legend=F)
dev.off()

# Here is another (read the label text):
pdf(file="random_bestShift2.pdf")
plot.bammshifts(dsc, edata, rank=1, legend=F)
dev.off()

# Here is one random sample with the 2nd BEST shift configuration
pdf(file="random_2ndbestShift1.pdf")
plot.bammshifts(dsc, edata, rank=2, legend=F)
dev.off()

# Here is another (read the label text):
pdf(file="random_2ndbestShift2.pdf")
plot.bammshifts(dsc, edata, rank=2, legend=F)
dev.off()

##
# Now using APE functions...

mysample <- 25  # this is the sample we'll plot

nrow(edata$eventData[[ mysample ]])

shiftnodes <- getShiftNodesFromIndex(edata, index = mysample)

pdf(file="apePlot_sample25.pdf", width=8.5,height=150)
plot.phylo(tree, cex=0.15)
nodelabels(node = shiftnodes, pch=21, col="red", cex=1.5)
dev.off()

# marginal shifts...

marg_probs <- marginalShiftProbsTree(edata)

pdf(file="marginalProbs.pdf", width=8.5,height=150)
plot.phylo(marg_probs, cex=0.15)
dev.off()

# clade-specific rates (sp and ex...)
allrates <- getCladeRates(edata)

write.table("ALL - LAMBDA", file="results.txt", append=TRUE)

write.table(mean(allrates$lambda), file="results.txt", append=TRUE)

write.table(quantile(allrates$lambda, c(0.05, 0.95)), file="results.txt", append=TRUE)

write.table("ALL - MU", file="results.txt", append=TRUE)

write.table(mean(allrates$mu), file="results.txt", append=TRUE)

write.table(quantile(allrates$mu, c(0.05, 0.95)), file="results.txt", append=TRUE)

# DOLPHIN compare
dolphinrates <- getCladeRates(edata, node= 8185) #141)

write.table("DOLPHINS", file="results.txt", append=TRUE)

write.table(mean(dolphinrates$lambda), file="results.txt", append=TRUE)

write.table(quantile(dolphinrates$lambda, c(0.05, 0.95)), file="results.txt", append=TRUE)

nondolphinrate <- getCladeRates(edata, node = 8185, nodetype = "exclude") #141, nodetype = "exclude")

write.table("NON-DOLPHINS", file="results.txt", append=TRUE)

write.table(mean(nondolphinrate$lambda), file="results.txt", append=TRUE)

write.table(quantile(nondolphinrate$lambda, c(0.05, 0.95)), file="results.txt", append=TRUE)

# per-branch rates
# components edata$meanTipLambda and edata$meanTipMu are the relevant model-averaged mean rates of speciation and extinction at the tips of the tree.
## How does this compare to DR??

BAMM_harmMeans_Lam<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=edata$tip.label)
BAMM_arithMeans_Lam<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=edata$tip.label)
BAMM_harmMeans_Mu<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=edata$tip.label)
BAMM_arithMeans_Mu<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=edata$tip.label)

for (k in 1:length(edata$tipLambda[[1]])){

tipLam<-vector()
tipMu<-vector()
for (j in 1:length(edata$tipLambda)){
	tipLam[j]<-edata$tipLambda[[j]][k]
	tipMu[j]<-edata$tipMu[[j]][k]
}
	BAMM_harmMeans_Lam[k,]<-1/(mean(1/tipLam))
	BAMM_arithMeans_Lam[k,]<-mean(tipLam)
	BAMM_harmMeans_Mu[k,]<-1/(mean(1/tipMu))
	BAMM_arithMeans_Mu[k,]<-mean(tipMu)
}

res<-cbind(BAMM_harmMeans_Lam,BAMM_arithMeans_Lam,BAMM_harmMeans_Mu,BAMM_arithMeans_Mu)

colnames(res)<-c("harmMeans_Lam","arithMeans_Lam","harmMeans_Mu","arithMeans_Mu")

write.table(res, paste("meanTipRates_sample10_tree",folders[i],".txt", sep=""))

#plot(density(res$harmMeans_Lam), ylim=c(0,8))
#plot(density(res$arithMeans_Lam), ylim=c(0,8))
#plot(density(edata$meanTipLambda), ylim=c(0,8)) # It is identical to the arithmetic means, good.

#####

# RATES through time...
pdf(file="rateThroughTime_spec.pdf")
plotRateThroughTime(edata, ratetype="speciation")
dev.off()

pdf(file="rateThroughTime_ext.pdf")
plotRateThroughTime(edata, ratetype="extinction")
dev.off()

pdf(file="rateThroughTime_div.pdf")
plotRateThroughTime(edata, ratetype="netdiv")
dev.off()

# cohort matrix...
cmat <- getCohortMatrix(edata)
pdf(file="cohortMatrixPlot.pdf")
cohorts(cmat, edata)
dev.off()

# cumulative shift marg_probs
cst <- cumulativeShiftProbsTree(edata)
pdf(file="cumShiftProbs_tree.pdf", width=8.5,height=150)
plot.phylo(cst, cex=0.15)
dev.off()


cst <- cumulativeShiftProbsTree(edata)
edgecols <- rep('black', length(tree$edge.length))
is_highprobshift <- cst$edge.length >= 0.95
edgecols[ is_highprobshift ] <- "red"

pdf(file="cumShiftProbs_tree_cols.pdf", width=35,height=55)
plot.phylo(tree, edge.color = edgecols, cex=0.05, type="phylogram")
dev.off()

pdf(file="cumShiftProbs_tree_cols_fan.pdf", width=35,height=35)
plot.phylo(tree, edge.color = edgecols, cex=0.15, type="fan")
dev.off()



}




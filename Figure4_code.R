# Figure 4 - R Code - MamPhy -- Upham et al. 2017
###
library(ape); library(phytools)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"

tipDataAll<-read.table("MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)
head(tipDataAll)

mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")

agesAllClades_origALL<-read.table(file=paste("divTime_24BAMMnodes_",bbone,"_MCC_target_DETAILS.txt", sep=""), header=TRUE)
agesAllClades_ALL<-agesAllClades_origALL[which(agesAllClades_origALL$Indep==1),]

nodeNums<-agesAllClades_ALL$nodeMCC
#nodeNumsAll<-agesAllClades_ALL$nodeMCC
#firstShift<-c(1,3,7,8,10,11,13,15,17:29,31,32) #5 == rm Placentalia
#nodeNums<-nodeNumsAll[firstShift]
#shiftDivDat1<-agesAllClades_ALL[firstShift,c("CLADE_label","ID_label","nodeMCC","numTaxa","avgFactor","avgIncDec","avgCladeDiv","cladeDiv_low95","cladeDiv_high95","mean","lower","upper")]
#colnames(shiftDivDat1)<-c("clade","ID","nodeMCC","numTaxa","factor","shift","cladeDiv","cladeDiv_low","cladeDiv_high","mean","lower","upper")
#shiftDivDat<-shiftDivDat1[order(shiftDivDat1$cladeDiv),]
#nodeNums<-shiftDivDat$nodeMCC # re-ordered

# gather the BAMM rate-shift data... 
shiftDivDat<-agesAllClades_ALL[order(agesAllClades_ALL$avgCladeDiv),]
nodeNums<-shiftDivDat$nodeMCC # re-ordered

shiftedNames<-shiftDivDat$clade
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






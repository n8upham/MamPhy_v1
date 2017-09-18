library(ape)
library(phytools)
library(phyloch)
library(moments)
library(nlme)
library(geiger)
library(phylolm)
library(RPANDA)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")

#===============
# Doing this across 100 tree samples
# FREQUENCY--- getting results for EACH of the 100 trees as a result
# ================
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/_NDexp_divModels_100_final")
#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_nexus-and-newickTrees")
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")

library(ape); library(picante); library(phytools); library(geiger); library(DDD); library(pspline);library(TreePar)

library(foreach);library(doSNOW)
cl = makeCluster(35, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees=100

foreach(i=1:ntrees, .packages=c('ape', 'picante', 'phytools','geiger','DDD','pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar% {

# which backbone?
bbone<- "NDexp" #"FBD" # 

# load tree
mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_nexus.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 
#rootTime<-max(node.age(mamPhy)$ages)

# tip categories
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""))
colnames(cladesDR)<-c("tiplabel","gen","fam","famLabel","famNumLabel","famNumAll","ord","ordLabel1","ordLabel2","clade","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")
head(cladesDR)

ordNames<-names(which(table(cladesDR$ord) > 4))
famNames<-names(which(table(cladesDR$fam) > 4))
patchNames<-names(which(table(cladesDR$PC) > 4))
cladeNames<-names(which(table(cladesDR$clade) > 4))

redoNames<-c("Ruminantia", "Yangochiroptera", "OldWorldMonkey2", "NewWorldMonkey2")

ruminants<-cladesDR[which(cladesDR$fam=="BOVIDAE" | cladesDR$fam=="MOSCHIDAE" | cladesDR$fam=="ANTILOCAPRIDAE" | cladesDR$fam=="GIRAFFIDAE" | cladesDR$fam=="TRAGULIDAE" | cladesDR$fam=="CERVIDAE"),"tiplabel"]
yango<-cladesDR[which(cladesDR$PC=="PC17_VespertilionidRelated" | cladesDR$PC=="PC18_EmballonuridRelated" | cladesDR$PC=="PC16_PhyllostomidRelated"),"tiplabel"]
OWM<-cladesDR[which(cladesDR$fam=="CERCOPITHECIDAE" | cladesDR$fam=="HYLOBATIDAE" |  cladesDR$fam=="HOMINIDAE"),"tiplabel"] 
NWM<-cladesDR[which(cladesDR$fam=="PITHECIIDAE" | cladesDR$fam=="ATELIDAE" |  cladesDR$fam=="AOTIDAE" |  cladesDR$fam=="CEBIDAE" |  cladesDR$fam=="CALLITRICHIDAE"),"tiplabel"] 

redoTipNames<-list(ruminants,yango,OWM,NWM)

# REDOS breakout
mamPhy_REDOS<-vector("list",length(redoNames))
for (j in 1:length(redoTipNames)){
    toDrop<-setdiff(mamPhy$tip.label,redoTipNames[[j]])
    mamPhy_REDOS[[j]]<-drop.tip(mamPhy,toDrop)
    }

source("run_diversification_analyses.R")

# Morlon 
# ON CLADES-- REDOS 
for (j in 1:length(redoNames)){
    run_Morlon_models(tree_file=mamPhy_REDOS[[j]], clade=paste(redoNames[[j]],"_",bbone,"_tree",i,sep=""),sampling_fraction=1, number_of_trees=1)    
}

}


# just redo NWM
ntrees=100

foreach(i=1:ntrees, .packages=c('ape', 'picante', 'phytools','geiger','DDD','pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar% {

# which backbone?
bbone<- "NDexp" #"FBD" # 

# load tree
mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_nexus.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 

cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""))
colnames(cladesDR)<-c("tiplabel","gen","fam","famLabel","famNumLabel","famNumAll","ord","ordLabel1","ordLabel2","clade","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")
head(cladesDR)

NWM<-cladesDR[which(cladesDR$fam=="PITHECIIDAE" | cladesDR$fam=="ATELIDAE" |  cladesDR$fam=="AOTIDAE" |  cladesDR$fam=="CEBIDAE" |  cladesDR$fam=="CALLITRICHIDAE"),"tiplabel"] 
toDrop<-setdiff(mamPhy$tip.label,NWM)
mamPhy_NWM<-drop.tip(mamPhy,toDrop)

run_Morlon_models(tree_file=mamPhy_NWM, clade=paste("NewWorldMonkey2","_",bbone,"_tree",i,sep=""),sampling_fraction=1, number_of_trees=1)    

}


# CLADES breakout
cladeTipNames<-vector("list",length(cladeNames))
for (j in 1:length(cladeNames)){
    cladeTipNames[[j]]<-cladesDR[which(cladesDR$clade==cladeNames[j]),"tiplabel"]
}

mamPhy_clades<-vector("list",length(cladeNames))
for (j in 1:length(cladeTipNames)){
    toDrop<-setdiff(mamPhy$tip.label,cladeTipNames[[j]])
    mamPhy_clades[[j]]<-drop.tip(mamPhy,toDrop)
    }

source("run_diversification_analyses.R")

# Morlon 
# ON CLADES 
for (j in 1:length(cladeNames)){
    run_Morlon_models(tree_file=mamPhy_clades[[j]], clade=paste(cladeNames[[j]],"_",bbone,"_tree",i,sep=""),sampling_fraction=1, number_of_trees=1)    
}

}



# ORDERS breakout
ordTipNames<-vector("list",length(ordNames))
for (j in 1:length(ordNames)){
    ordTipNames[[j]]<-cladesDR[which(cladesDR$ord==ordNames[j]),"tiplabel"]
}

mamPhy_ords<-vector("list",length(ordNames))
for (j in 1:length(ordTipNames)){
    toDrop<-setdiff(mamPhy$tip.label,ordTipNames[[j]])
    mamPhy_ords[[j]]<-drop.tip(mamPhy,toDrop)
    }

# PATCHES breakout
patchTipNames<-vector("list",length(patchNames))
for (j in 1:length(patchNames)){
    patchTipNames[[j]]<-cladesDR[which(cladesDR$PC==patchNames[j]),"tiplabel"]
}

mamPhy_pcs<-vector("list",length(patchNames))
for (j in 1:length(patchTipNames)){
    toDrop<-setdiff(mamPhy$tip.label,patchTipNames[[j]])
    mamPhy_pcs[[j]]<-drop.tip(mamPhy,toDrop)
    }

# MAMMAL SUBCLADE BREAKOUTS
Monotremata<-cladesDR[which(cladesDR$PC=="PC23_Monotremata"),"tiplabel"]
Marsupalia<-cladesDR[which(cladesDR$PC=="PC1_Marsupials"),"tiplabel"]
Placentalia<-setdiff(mamPhy$tip.label,(cladesDR[which(cladesDR$PC=="PC1_Marsupials" | cladesDR$PC=="PC23_Monotremata"),"tiplabel"]))
Afrotheria<-cladesDR[which(cladesDR$PC=="PC2_Afrotheria"),"tiplabel"]
Xenarthra<-cladesDR[which(cladesDR$PC=="PC3_Xenarthra"),"tiplabel"]
Euarchontoglires<-cladesDR[which(cladesDR$ord=="RODENTIA" | cladesDR$ord=="LAGOMORPHA" | cladesDR$ord=="PRIMATES" | cladesDR$ord=="SCANDENTIA" | cladesDR$ord=="DERMOPTERA"),"tiplabel"]
Laurasiatheria<-cladesDR[which(cladesDR$ord=="EULIPOTYPHLA" | cladesDR$ord=="CHIROPTERA" | cladesDR$ord=="PHOLIDOTA" | cladesDR$ord=="CARNIVORA" | cladesDR$ord=="PERISSODACTYLA" | cladesDR$ord=="CETARTIODACTYLA"),"tiplabel"]

mamSubsTips<-list(Monotremata,Marsupalia,Placentalia,Afrotheria,Xenarthra,Euarchontoglires,Laurasiatheria)
mamSubsNames<-c("Monotremata","Marsupalia","Placentalia","Afrotheria","Xenarthra","Euarchontoglires","Laurasiatheria")

mamSubs<-vector("list",length(mamSubsTips))
for (j in 1:length(mamSubsTips)){
    toDrop<-setdiff(mamPhy$tip.label,mamSubsTips[[j]])
    mamSubs[[j]]<-drop.tip(mamPhy,toDrop)
    }

#######
# CONDAMINE code for running 
# MORLON: all Yule, BD con + 4 BD timeDep Expo + 4 BD timeDepLin
# COND: BD con + 4 BD tempDep Expo + 4 BD tempDep Lin
# ETTIENE: 6 DD models (linear and expo combos)
# STADLER: TreePar model of RATE SHIFTS (up to 4) through time

source("run_diversification_analyses.R")

# Morlon 
# ON ORDERS 
for (j in 1:length(ordNames)){
    run_Morlon_models(tree_file=mamPhy_ords[[j]], clade=paste(ordNames[[j]],"_",bbone,"_tree",i,sep=""),sampling_fraction=1, number_of_trees=1)    
}
# ON PATCHES
for (j in 1:length(patchNames)){
    run_Morlon_models(tree_file=mamPhy_pcs[[j]], clade=paste(patchNames[[j]],"_",bbone,"_tree",i,sep=""),sampling_fraction=1, number_of_trees=1)    
}

# ON MAMMAL SUBCLADES
for (j in 1:length(mamSubsNames)){
    run_Morlon_models(tree_file=mamSubs[[j]], clade=paste(mamSubsNames[[j]],"_",bbone,"_tree",i,sep=""),sampling_fraction=1, number_of_trees=1)    
}

}

########
# DDD
# ON ORDERS
for (j in 1:length(ordNames)){
    run_DDD(tree_file=mamPhy_ords[[j]], clade=paste(ordNames[[j]],"_",bbone,"_tree",i,sep=""),number_of_trees=1)    
}

# ON PATCHES
for (j in 1:length(patchNames)){
    run_DDD(tree_file=mamPhy_pcs[[j]], clade=paste(patchNames[[j]],"_",bbone,"_tree",i,sep=""),number_of_trees=1)    
}

# ON MAMMAL SUBCLADES
for (j in 1:length(mamSubsNames)){
    run_DDD(tree_file=mamSubs[[j]], clade=paste(mamSubsNames[[j]],"_",bbone,"_tree",i,sep=""),number_of_trees=1)    
}

######
# on FULL trees
run_Morlon_models(tree_file=mamPhy, clade=paste("allMammalia_",bbone,"_tree",i,sep=""), sampling_fraction=1, number_of_trees=1)

run_DDD(tree_file=mamPhy, clade=paste("allMammalia_",bbone,"_tree",i,sep=""),number_of_trees=1)    

############
# ==============
# SUMMARIZE the ML modeling results across the 100 trees
# ==============
#######
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/_NDexp_divModels_100_final")
#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_nexus-and-newickTrees")
library(ape); library(picante); library(phytools); library(geiger); library(DDD); library(pspline);library(TreePar)

ntrees=100

setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_nexus-and-newickTrees")

# which backbone?
bbone<- "NDexp" #"FBD" # 
i=1
mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_nexus.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 

setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_divModels_100_final")
# tip categories
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
#colnames(cladesDR)<-c("tiplabel","gen","fam","famLabel","famNumLabel","famNumAll","ord","ordLabel1","ordLabel2","clade","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")
head(cladesDR)

ordNames<-names(which(table(cladesDR$ord) > 4))
famNames<-names(which(table(cladesDR$fam) > 4))
patchNames<-names(which(table(cladesDR$PC) > 4))
cladeNames<-names(which(table(cladesDR$clade) > 4))

redoNames<-c("Ruminantia", "Yangochiroptera", "OldWorldMonkey2", "NewWorldMonkey2")

# MAMMAL SUBCLADE BREAKOUTS
Monotremata<-cladesDR[which(cladesDR$PC=="PC23_Monotremata"),"tiplabel"]
Marsupalia<-cladesDR[which(cladesDR$PC=="PC1_Marsupials"),"tiplabel"]
Placentalia<-setdiff(mamPhy$tip.label,(cladesDR[which(cladesDR$PC=="PC1_Marsupials" | cladesDR$PC=="PC23_Monotremata"),"tiplabel"]))
Afrotheria<-cladesDR[which(cladesDR$PC=="PC2_Afrotheria"),"tiplabel"]
Xenarthra<-cladesDR[which(cladesDR$PC=="PC3_Xenarthra"),"tiplabel"]
Euarchontoglires<-cladesDR[which(cladesDR$ord=="RODENTIA" | cladesDR$ord=="LAGOMORPHA" | cladesDR$ord=="PRIMATES" | cladesDR$ord=="SCANDENTIA" | cladesDR$ord=="DERMOPTERA"),"tiplabel"]
Laurasiatheria<-cladesDR[which(cladesDR$ord=="EULIPOTYPHLA" | cladesDR$ord=="CHIROPTERA" | cladesDR$ord=="PHOLIDOTA" | cladesDR$ord=="CARNIVORA" | cladesDR$ord=="PERISSODACTYLA" | cladesDR$ord=="CETARTIODACTYLA"),"tiplabel"]

mamSubsTips<-list(Monotremata,Marsupalia,Placentalia,Afrotheria,Xenarthra,Euarchontoglires,Laurasiatheria)
mamSubsNames<-c("Monotremata","Marsupalia","Placentalia","Afrotheria","Xenarthra","Euarchontoglires","Laurasiatheria")

# LOAD data back in
# to summarize per clade across the 100 trees, want all 100 trees for a given unit
###
# REDOS
redoNames<-c("Ruminantia", "Yangochiroptera", "OldWorldMonkey2", "NewWorldMonkey2")
cladeNames<-redoNames

# CLADES
Unit<-vector("list",length=ntrees)
bestModels<-vector("list",length=ntrees)
bestModelsGrThan1<-vector("list",length=ntrees)
bestModelsGrThan2<-vector("list",length=ntrees)
deltaAIC_RCs<-vector("list",length=ntrees)
cladeRes_all100<-vector("list",length(cladeNames))
cladeRes_bestModels<-vector("list",length(cladeNames))
cladeRes_bestModelsGrThan1<-vector("list",length(cladeNames))
cladeRes_bestModelsGrThan2<-vector("list",length(cladeNames))
cladeRes_deltaAIC_RCs<-vector("list",length(cladeNames))
for (j in 1:length(cladeNames)){
    for (i in 1:ntrees){
    res<-read.table(paste(cladeNames[j],"_",bbone,"_tree",i,"_results_Morlon.txt",sep=""), header=TRUE)  
    res2<-cbind(res[,1:8],res$means.AICc-min(res$means.AICc),rep(i,10))
    colnames(res2)<-c(colnames(res)[1:8],"deltaAICc","tree")
    Unit[[i]]<-res2
    bestModels[[i]]<-res2[which(res2$means.AICc==min(res2$means.AICc)),"means.Models"]
    sortResCON<-sort(res2$deltaAICc[1:2], decreasing=FALSE)
    sortResVAR<-sort(res2$deltaAICc[3:10], decreasing=FALSE)
    deltaAIC_RCs[[i]]<-sortResCON[1]-sortResVAR[1]
    if (abs(deltaAIC_RCs[[i]]) > 1) { 
        bestModelsGrThan1[[i]]<-res2[which(res2$means.AICc==min(res2$means.AICc)),"means.Models"]
        } else {
        bestModelsGrThan1[[i]]<-as.factor("NA")    
        }
    if (abs(deltaAIC_RCs[[i]]) > 2) { 
        bestModelsGrThan2[[i]]<-res2[which(res2$means.AICc==min(res2$means.AICc)),"means.Models"]
        } else {
        bestModelsGrThan2[[i]]<-as.factor("NA")    
        }
    }
    cladeRes_all100[[j]]<-do.call(rbind,Unit)
    cladeRes_deltaAIC_RCs[[j]]<-sort(unlist(deltaAIC_RCs))
    cladeRes_bestModels[[j]]<-table(unlist(bestModels))
    cladeRes_bestModelsGrThan1[[j]]<-as.data.frame(table(unlist(bestModelsGrThan1)))
    cladeRes_bestModelsGrThan2[[j]]<-as.data.frame(table(unlist(bestModelsGrThan2)))
}
modelNames<-res$means.Models

cladeRes_sum<-vector("list",length(cladeNames))
for (j in 1:length(cladeNames)){
    res<-cladeRes_all100[[j]]
    summaries<-vector("list",length(modelNames))    
    for (k in 1:length(modelNames)){
        model<-as.vector(modelNames[k])
        res2<-res[which(res$means.Models==model),]
        res3<-res2[,3:9]
        sumRes<-data.frame(matrix(NA, nrow = length(res3)*3+1, ncol = 1),row.names=c("means.logL","q1.logL","q2.logL", "means.AICc","q1.AICc","q2.AICc", "means.Lambda","q1.Lambda","q2.Lambda", "means.AlphaTime","q1.AlphaTime","q2.AlphaTime", "means.Mu","q1.Mu","q2.Mu", "means.BetaTime","q1.BetaTime","q2.BetaTime", "means.deltaAICc","q1.deltaAICc","q2.deltaAICc","numTimesBEST"))
        colnames(sumRes)<-model
        for (l in 0:(length(3:9)-1)){
            sumRes[3*l+1,]<-mean(res3[,1+l])
            quants<-quantile(res3[,1+l],c(0.025,0.975),na.rm=TRUE)
            sumRes[3*l+2,]<-quants[[1]]
            sumRes[3*l+3,]<-quants[[2]]
            sumRes[22,]<-cladeRes_bestModels[[j]][[k]]
        }
        summaries[[k]]<-sumRes
    }
    cladeRes_sum[[j]]<-do.call(cbind,summaries)
    write.table(cladeRes_sum[[j]],paste(cladeNames[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees.txt",sep=""))
    write.table(cladeRes_bestModelsGrThan1[[j]],paste(cladeNames[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_grThan1.txt",sep=""))
    write.table(cladeRes_deltaAIC_RCs[[j]],paste(cladeNames[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_deltaAIC-RCs.txt",sep=""))
}


# ORDERS
Unit<-vector("list",length=ntrees)
bestModels<-vector("list",length=ntrees)
bestModelsGrThan1<-vector("list",length=ntrees)
bestModelsGrThan2<-vector("list",length=ntrees)
deltaAIC_RCs<-vector("list",length=ntrees)
ordRes_all100<-vector("list",length(ordNames))
ordRes_bestModels<-vector("list",length(ordNames))
ordRes_bestModelsGrThan1<-vector("list",length(ordNames))
ordRes_bestModelsGrThan2<-vector("list",length(ordNames))
ordRes_deltaAIC_RCs<-vector("list",length(ordNames))
for (j in 1:length(ordNames)){
    for (i in 1:ntrees){
    res<-read.table(paste(ordNames[j],"_",bbone,"_tree",i,"_results_Morlon.txt",sep=""), header=TRUE)  
    res2<-cbind(res[,1:8],res$means.AICc-min(res$means.AICc),rep(i,10))
    colnames(res2)<-c(colnames(res)[1:8],"deltaAICc","tree")
    Unit[[i]]<-res2
    bestModels[[i]]<-res2[which(res2$means.AICc==min(res2$means.AICc)),"means.Models"]
    sortResCON<-sort(res2$deltaAICc[1:2], decreasing=FALSE)
    sortResVAR<-sort(res2$deltaAICc[3:10], decreasing=FALSE)
    deltaAIC_RCs[[i]]<-sortResCON[1]-sortResVAR[1]
    if (abs(deltaAIC_RCs[[i]]) > 1) { 
        bestModelsGrThan1[[i]]<-res2[which(res2$means.AICc==min(res2$means.AICc)),"means.Models"]
        } else {
        bestModelsGrThan1[[i]]<-as.factor("NA")    
        }
    if (abs(deltaAIC_RCs[[i]]) > 2) { 
        bestModelsGrThan2[[i]]<-res2[which(res2$means.AICc==min(res2$means.AICc)),"means.Models"]
        } else {
        bestModelsGrThan2[[i]]<-as.factor("NA")    
        }
    }
    ordRes_all100[[j]]<-do.call(rbind,Unit)
    ordRes_deltaAIC_RCs[[j]]<-sort(unlist(deltaAIC_RCs))
    ordRes_bestModels[[j]]<-table(unlist(bestModels))
    ordRes_bestModelsGrThan1[[j]]<-as.data.frame(table(unlist(bestModelsGrThan1)))
    ordRes_bestModelsGrThan2[[j]]<-as.data.frame(table(unlist(bestModelsGrThan2)))
}
modelNames<-res$means.Models

ordRes_sum<-vector("list",length(ordNames))
for (j in 1:length(ordNames)){
    res<-ordRes_all100[[j]]
    summaries<-vector("list",length(modelNames))    
    for (k in 1:length(modelNames)){
        model<-as.vector(modelNames[k])
        res2<-res[which(res$means.Models==model),]
        res3<-res2[,3:9]
        sumRes<-data.frame(matrix(NA, nrow = length(res3)*3+1, ncol = 1),row.names=c("means.logL","q1.logL","q2.logL", "means.AICc","q1.AICc","q2.AICc", "means.Lambda","q1.Lambda","q2.Lambda", "means.AlphaTime","q1.AlphaTime","q2.AlphaTime", "means.Mu","q1.Mu","q2.Mu", "means.BetaTime","q1.BetaTime","q2.BetaTime", "means.deltaAICc","q1.deltaAICc","q2.deltaAICc","numTimesBEST"))
        colnames(sumRes)<-model
        for (l in 0:(length(3:9)-1)){
            sumRes[3*l+1,]<-mean(res3[,1+l])
            quants<-quantile(res3[,1+l],c(0.025,0.975),na.rm=TRUE)
            sumRes[3*l+2,]<-quants[[1]]
            sumRes[3*l+3,]<-quants[[2]]
            sumRes[22,]<-ordRes_bestModels[[j]][[k]]
        }
        summaries[[k]]<-sumRes
    }
    ordRes_sum[[j]]<-do.call(cbind,summaries)
    write.table(ordRes_sum[[j]],paste(ordNames[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees.txt",sep=""))
    write.table(ordRes_bestModelsGrThan1[[j]],paste(ordNames[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_grThan1.txt",sep=""))
    write.table(ordRes_deltaAIC_RCs[[j]],paste(ordNames[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_deltaAIC-RCs.txt",sep=""))
}


# PATCHES
Unit<-vector("list",length=ntrees)
bestModels<-vector("list",length=ntrees)
bestModelsGrThan1<-vector("list",length=ntrees)
bestModelsGrThan2<-vector("list",length=ntrees)
deltaAIC_RCs<-vector("list",length=ntrees)
patchRes_all100<-vector("list",length(patchNames))
patchRes_bestModels<-vector("list",length(patchNames))
patchRes_bestModelsGrThan1<-vector("list",length(patchNames))
patchRes_bestModelsGrThan2<-vector("list",length(patchNames))
patchRes_deltaAIC_RCs<-vector("list",length(patchNames))
for (j in 1:length(patchNames)){
    for (i in 1:ntrees){
    res<-read.table(paste(patchNames[j],"_",bbone,"_tree",i,"_results_Morlon.txt",sep=""), header=TRUE)  
    res2<-cbind(res[,1:8],res$means.AICc-min(res$means.AICc),rep(i,10))
    colnames(res2)<-c(colnames(res)[1:8],"deltaAICc","tree")
    Unit[[i]]<-res2
    bestModels[[i]]<-res2[which(res2$means.AICc==min(res2$means.AICc)),"means.Models"]
    sortResCON<-sort(res2$deltaAICc[1:2], decreasing=FALSE)
    sortResVAR<-sort(res2$deltaAICc[3:10], decreasing=FALSE)
    deltaAIC_RCs[[i]]<-sortResCON[1]-sortResVAR[1]
    if (abs(deltaAIC_RCs[[i]]) > 1) { 
        bestModelsGrThan1[[i]]<-res2[which(res2$means.AICc==min(res2$means.AICc)),"means.Models"]
        } else {
        bestModelsGrThan1[[i]]<-as.factor("NA")    
        }
    if (abs(deltaAIC_RCs[[i]]) > 2) { 
        bestModelsGrThan2[[i]]<-res2[which(res2$means.AICc==min(res2$means.AICc)),"means.Models"]
        } else {
        bestModelsGrThan2[[i]]<-as.factor("NA")    
        }
    }
    patchRes_all100[[j]]<-do.call(rbind,Unit)
    patchRes_deltaAIC_RCs[[j]]<-sort(unlist(deltaAIC_RCs))
    patchRes_bestModels[[j]]<-table(unlist(bestModels))
    patchRes_bestModelsGrThan1[[j]]<-as.data.frame(table(unlist(bestModelsGrThan1)))
    patchRes_bestModelsGrThan2[[j]]<-as.data.frame(table(unlist(bestModelsGrThan2)))
}
modelNames<-res$means.Models

patchRes_sum<-vector("list",length(patchNames))
for (j in 1:length(patchNames)){
    res<-patchRes_all100[[j]]
    summaries<-vector("list",length(modelNames))    
    for (k in 1:length(modelNames)){
        model<-as.vector(modelNames[k])
        res2<-res[which(res$means.Models==model),]
        res3<-res2[,3:9]
        sumRes<-data.frame(matrix(NA, nrow = length(res3)*3+1, ncol = 1),row.names=c("means.logL","q1.logL","q2.logL", "means.AICc","q1.AICc","q2.AICc", "means.Lambda","q1.Lambda","q2.Lambda", "means.AlphaTime","q1.AlphaTime","q2.AlphaTime", "means.Mu","q1.Mu","q2.Mu", "means.BetaTime","q1.BetaTime","q2.BetaTime", "means.deltaAICc","q1.deltaAICc","q2.deltaAICc","numTimesBEST"))
        colnames(sumRes)<-model
        for (l in 0:(length(3:9)-1)){
            sumRes[3*l+1,]<-mean(res3[,1+l])
            quants<-quantile(res3[,1+l],c(0.025,0.975),na.rm=TRUE)
            sumRes[3*l+2,]<-quants[[1]]
            sumRes[3*l+3,]<-quants[[2]]
            sumRes[22,]<-patchRes_bestModels[[j]][[k]]
        }
        summaries[[k]]<-sumRes
    }
    patchRes_sum[[j]]<-do.call(cbind,summaries)
    write.table(patchRes_sum[[j]],paste(patchNames[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees.txt",sep=""))
    write.table(patchRes_bestModelsGrThan1[[j]],paste(patchNames[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_grThan1.txt",sep=""))
    write.table(patchRes_deltaAIC_RCs[[j]],paste(patchNames[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_deltaAIC-RCs.txt",sep=""))
}



# SUBCLADES
Unit<-vector("list",length=ntrees)
bestModels<-vector("list",length=ntrees)
bestModelsGrThan1<-vector("list",length=ntrees)
bestModelsGrThan2<-vector("list",length=ntrees)
deltaAIC_RCs<-vector("list",length=ntrees)
mamSubsRes_all100<-vector("list",length(mamSubsNames))
mamSubsRes_bestModels<-vector("list",length(mamSubsNames))
mamSubsRes_bestModelsGrThan1<-vector("list",length(mamSubsNames))
mamSubsRes_bestModelsGrThan2<-vector("list",length(mamSubsNames))
mamSubsRes_deltaAIC_RCs<-vector("list",length(mamSubsNames))
for (j in 1:length(mamSubsNames)){
    for (i in 1:ntrees){
    res<-read.table(paste(mamSubsNames[j],"_",bbone,"_tree",i,"_results_Morlon.txt",sep=""), header=TRUE)  
    res2<-cbind(res[,1:8],res$means.AICc-min(res$means.AICc),rep(i,10))
    colnames(res2)<-c(colnames(res)[1:8],"deltaAICc","tree")
    Unit[[i]]<-res2
    bestModels[[i]]<-res2[which(res2$means.AICc==min(res2$means.AICc)),"means.Models"]
    sortResCON<-sort(res2$means.AICc[1:2], decreasing=FALSE)
    sortResVAR<-sort(res2$means.AICc[3:10], decreasing=FALSE)
    deltaAIC_RCs[[i]]<-sortResCON[1]-sortResVAR[1]
    if (abs(deltaAIC_RCs[[i]]) > 1) { 
        bestModelsGrThan1[[i]]<-res2[which(res2$means.AICc==min(res2$means.AICc)),"means.Models"]
        } else {
        bestModelsGrThan1[[i]]<-as.factor("NA")    
        }
    if (abs(deltaAIC_RCs[[i]]) > 2) { 
        bestModelsGrThan2[[i]]<-res2[which(res2$means.AICc==min(res2$means.AICc)),"means.Models"]
        } else {
        bestModelsGrThan2[[i]]<-as.factor("NA")    
        }
    }
    mamSubsRes_all100[[j]]<-do.call(rbind,Unit)
    mamSubsRes_deltaAIC_RCs[[j]]<-sort(unlist(deltaAIC_RCs))
    mamSubsRes_bestModels[[j]]<-table(unlist(bestModels))
    mamSubsRes_bestModelsGrThan1[[j]]<-as.data.frame(table(unlist(bestModelsGrThan1)))
    mamSubsRes_bestModelsGrThan2[[j]]<-as.data.frame(table(unlist(bestModelsGrThan2)))
}
modelNames<-res$means.Models

mamSubsRes_sum<-vector("list",length(mamSubsNames))
for (j in 1:length(mamSubsNames)){
    res<-mamSubsRes_all100[[j]]
    summaries<-vector("list",length(modelNames))    
    for (k in 1:length(modelNames)){
        model<-as.vector(modelNames[k])
        res2<-res[which(res$means.Models==model),]
        res3<-res2[,3:9]
        sumRes<-data.frame(matrix(NA, nrow = length(res3)*3+1, ncol = 1),row.names=c("means.logL","q1.logL","q2.logL", "means.AICc","q1.AICc","q2.AICc", "means.Lambda","q1.Lambda","q2.Lambda", "means.AlphaTime","q1.AlphaTime","q2.AlphaTime", "means.Mu","q1.Mu","q2.Mu", "means.BetaTime","q1.BetaTime","q2.BetaTime", "means.deltaAICc","q1.deltaAICc","q2.deltaAICc","numTimesBEST"))
        colnames(sumRes)<-model
        for (l in 0:(length(3:9)-1)){
            sumRes[3*l+1,]<-mean(res3[,1+l])
            quants<-quantile(res3[,1+l],c(0.025,0.975),na.rm=TRUE)
            sumRes[3*l+2,]<-quants[[1]]
            sumRes[3*l+3,]<-quants[[2]]
            sumRes[22,]<-mamSubsRes_bestModels[[j]][[k]]
        }
        summaries[[k]]<-sumRes
    }
    mamSubsRes_sum[[j]]<-do.call(cbind,summaries)
    write.table(mamSubsRes_sum[[j]],paste(mamSubsNames[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees.txt",sep=""))
    write.table(mamSubsRes_bestModelsGrThan1[[j]],paste(mamSubsNames[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_grThan1.txt",sep=""))
    write.table(mamSubsRes_deltaAIC_RCs[[j]],paste(mamSubsNames[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_deltaAIC-RCs.txt",sep=""))
}

# toPLOT
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_divModels_100_final")
orderToPlot<-c("MouseRelated", "SquirrelRelated", "GuineaPigRelated", "Yinpterochiroptera", "Yangochiroptera", "Soricidae", "Talpidae","Erinaceidae","OldWorldMonkey2", "NewWorldMonkey2", "Strepsirrhini", "Ruminantia", "Cetaceans", "AFROSORICIDA")#, "SCANDENTIA", "MACROSCELIDEA", "PILOSA", "CARNIVORA","CatRelated","DogRelated", "LAGOMORPHA", "PERISSODACTYLA", "CINGULATA", "PHOLIDOTA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA")

Unit<-vector("list",length=ntrees)
bestModels<-vector("list",length=ntrees)
bestModelsGrThan1<-vector("list",length=ntrees)
bestModelsGrThan2<-vector("list",length=ntrees)
deltaAIC_RCs<-vector("list",length=ntrees)
toPlotFig2_all100<-vector("list",length(orderToPlot))
toPlotFig2_bestModels<-vector("list",length(orderToPlot))
toPlotFig2_bestModelsGrThan1<-vector("list",length(orderToPlot))
toPlotFig2_bestModelsGrThan2<-vector("list",length(orderToPlot))
toPlotFig2_deltaAIC_RCs<-vector("list",length(orderToPlot))
for (j in 1:length(orderToPlot)){
    for (i in 1:ntrees){
    res<-read.table(paste(orderToPlot[j],"_",bbone,"_tree",i,"_results_Morlon.txt",sep=""), header=TRUE)  
    res2<-cbind(res[,1:8],res$means.AICc-min(res$means.AICc),rep(i,10))
    colnames(res2)<-c(colnames(res)[1:8],"deltaAICc","tree")
    Unit[[i]]<-res2
    bestModels[[i]]<-res2[which(res2$means.AICc==min(res2$means.AICc)),"means.Models"]
    sortResCON<-sort(res2$means.AICc[1:2], decreasing=FALSE)
    sortResVAR<-sort(res2$means.AICc[3:10], decreasing=FALSE)
    deltaAIC_RCs[[i]]<-sortResCON[1]-sortResVAR[1]
    if (abs(deltaAIC_RCs[[i]]) > 1) { 
        bestModelsGrThan1[[i]]<-res2[which(res2$means.AICc==min(res2$means.AICc)),"means.Models"]
        } else {
        bestModelsGrThan1[[i]]<-as.factor("NA")    
        }
    if (abs(deltaAIC_RCs[[i]]) > 2) { 
        bestModelsGrThan2[[i]]<-res2[which(res2$means.AICc==min(res2$means.AICc)),"means.Models"]
        } else {
        bestModelsGrThan2[[i]]<-as.factor("NA")    
        }
    }
    toPlotFig2_all100[[j]]<-do.call(rbind,Unit)
    toPlotFig2_deltaAIC_RCs[[j]]<-sort(unlist(deltaAIC_RCs))
    toPlotFig2_bestModels[[j]]<-table(unlist(bestModels))
    toPlotFig2_bestModelsGrThan1[[j]]<-as.data.frame(table(unlist(bestModelsGrThan1)))
    toPlotFig2_bestModelsGrThan2[[j]]<-as.data.frame(table(unlist(bestModelsGrThan2)))
}
modelNames<-res$means.Models

toPlotFig2_sum<-vector("list",length(orderToPlot))
for (j in 1:length(orderToPlot)){
    res<-toPlotFig2_all100[[j]]
    summaries<-vector("list",length(modelNames))    
    for (k in 1:length(modelNames)){
        model<-as.vector(modelNames[k])
        res2<-res[which(res$means.Models==model),]
        res3<-res2[,3:9]
        sumRes<-data.frame(matrix(NA, nrow = length(res3)*3+1, ncol = 1),row.names=c("means.logL","q1.logL","q2.logL", "means.AICc","q1.AICc","q2.AICc", "means.Lambda","q1.Lambda","q2.Lambda", "means.AlphaTime","q1.AlphaTime","q2.AlphaTime", "means.Mu","q1.Mu","q2.Mu", "means.BetaTime","q1.BetaTime","q2.BetaTime", "means.deltaAICc","q1.deltaAICc","q2.deltaAICc","numTimesBEST"))
        colnames(sumRes)<-model
        for (l in 0:(length(3:9)-1)){
            sumRes[3*l+1,]<-mean(res3[,1+l])
            quants<-quantile(res3[,1+l],c(0.025,0.975),na.rm=TRUE)
            sumRes[3*l+2,]<-quants[[1]]
            sumRes[3*l+3,]<-quants[[2]]
            sumRes[22,]<-toPlotFig2_bestModels[[j]][[k]]
        }
        summaries[[k]]<-sumRes
    }
    toPlotFig2_sum[[j]]<-do.call(cbind,summaries)
    write.table(toPlotFig2_sum[[j]],paste(orderToPlot[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees.txt",sep=""))
    write.table(toPlotFig2_bestModelsGrThan1[[j]],paste(orderToPlot[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_grThan1.txt",sep=""))
    write.table(toPlotFig2_deltaAIC_RCs[[j]],paste(orderToPlot[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_deltaAIC-RCs.txt",sep=""))
}




###
# Re-load the SUM tables:
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_divModels_100_final/res2_parameterSummaries_10models-100trees")
i=100

ordRes_sum<-vector("list",length(ordNames))
for (j in 1:length(ordRes_sum)){
    ordRes_sum[[j]]<-read.table(paste(ordNames[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees.txt",sep=""))
}

patchRes_sum<-vector("list",length(patchNames))
for (j in 1:length(patchRes_sum)){
    patchRes_sum[[j]]<-read.table(paste(patchNames[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees.txt",sep=""))
}

cladeRes_sum<-vector("list",length(cladeNames))
for (j in 1:length(cladeRes_sum)){
    cladeRes_sum[[j]]<-read.table(paste(cladeNames[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees.txt",sep=""))
}

subRes_sum<-vector("list",length(mamSubsNames))
for (j in 1:length(subRes_sum)){
    subRes_sum[[j]]<-read.table(paste(mamSubsNames[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees.txt",sep=""))
}

# load in....
i=100
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_divModels_100_final/res2_parameterSummaries_10models-100trees")
orderToPlot<-c("MouseRelated", "SquirrelRelated", "GuineaPigRelated", "Yinpterochiroptera", "Yangochiroptera", "Soricidae", "Talpidae","Erinaceidae","OldWorldMonkey2", "NewWorldMonkey2", "Strepsirrhini", "Ruminantia", "Cetaceans", "AFROSORICIDA")#, "SCANDENTIA", "MACROSCELIDEA", "PILOSA", "CARNIVORA","CatRelated","DogRelated", "LAGOMORPHA", "PERISSODACTYLA", "CINGULATA", "PHOLIDOTA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA")
toPlotFig2_sum<-vector("list",length(orderToPlot))
for (j in 1:length(toPlotFig2_sum)){
    toPlotFig2_sum[[j]]<-read.table(paste(orderToPlot[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees.txt",sep=""))
}

setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_divModels_100_final/res2_grThan1")
orderToPlot<-c("MouseRelated", "SquirrelRelated", "GuineaPigRelated", "Yinpterochiroptera", "Yangochiroptera", "Soricidae", "Talpidae","Erinaceidae","OldWorldMonkey2", "NewWorldMonkey2", "Strepsirrhini", "Ruminantia", "Cetaceans", "AFROSORICIDA")#, "SCANDENTIA", "MACROSCELIDEA", "PILOSA", "CARNIVORA","CatRelated","DogRelated", "LAGOMORPHA", "PERISSODACTYLA", "CINGULATA", "PHOLIDOTA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA")
toPlotFig2_bestModelsGrThan1<-vector("list",length(orderToPlot))
for (j in 1:length(toPlotFig2_sum)){
    toPlotFig2_bestModelsGrThan1[[j]]<-read.table(paste(orderToPlot[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_grThan1.txt",sep=""))
}

setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_divModels_100_final/res2_deltaAICrcs")
orderToPlot<-c("MouseRelated", "SquirrelRelated", "GuineaPigRelated", "Yinpterochiroptera", "Yangochiroptera", "Soricidae", "Talpidae","Erinaceidae","OldWorldMonkey2", "NewWorldMonkey2", "Strepsirrhini", "Ruminantia", "Cetaceans", "AFROSORICIDA")#, "SCANDENTIA", "MACROSCELIDEA", "PILOSA", "CARNIVORA","CatRelated","DogRelated", "LAGOMORPHA", "PERISSODACTYLA", "CINGULATA", "PHOLIDOTA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA")
toPlotFig2_deltaAIC_RCs<-vector("list",length(orderToPlot))
for (j in 1:length(toPlotFig2_sum)){
    toPlotFig2_deltaAIC_RCs[[j]]<-read.table(paste(orderToPlot[j],"_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_deltaAIC-RCs.txt",sep=""))
}

# SUM as time-CONSTANT vs. time-VARIABLE -- but KEEPING the best model info
# toPlotFig2 
dat<-data.frame(matrix(NA, nrow = length(orderToPlot), ncol = 17),row.names=orderToPlot)
colnames(dat)<-c("timeCon","timeVar","%timeCon","richness","GrThan1_TC","GrThan1_TV","GrThan1_NA", "GrThan2_TC","GrThan2_TV","GrThan2_NA", "mean_dAICrc","low95_dAICrc","up95_dAICrc","median_dAICrc","low50_dAICrc","up50_dAICrc","diffZero")
for (j in 1:length(toPlotFig2_sum)){
    res<-toPlotFig2_sum[[j]]["numTimesBEST",]
    GrThan1<-toPlotFig2_bestModelsGrThan1[[j]]
    #GrThan2<-cladeRes_bestModelsGrThan2[[j]]
    dAICrc<-toPlotFig2_deltaAIC_RCs[[j]][[1]]
    dat[j,1]<-sum(res[,1:2])
    dat[j,2]<-sum(res[,3:10])
    dat[j,3]<-round(sum(res[,1:2])/100,digits=2)
    dat[j,4]<-length(cladesDR[which(cladesDR$clade==orderToPlot[j]),"tiplabel"])
    dat[j,5]<-sum(GrThan1[which(GrThan1$Var1=="BCST" | GrThan1$Var1=="BCSTDCST"),"Freq"])
    NA1<-GrThan1[which(GrThan1$Var1=="NA"),"Freq"]
    dat[j,7]<-if(length(NA1)>0) {NA1} else {0}
    dat[j,6]<-100-(dat[j,5]+dat[j,7])
    #dat[j,8]<-sum(GrThan2[which(GrThan2$Var1=="BCST" | GrThan2$Var1=="BCSTDCST"),"Freq"])
    #NA2<-GrThan2[which(GrThan2$Var1=="NA"),"Freq"]
    #dat[j,10]<-if(length(NA2)>0) {NA2} else {0}
    #dat[j,9]<-100-(dat[j,8]+dat[j,10])
    dat[j,11]<-mean(dAICrc)
    dat[j,12]<-quantile(dAICrc,c(0.025,0.975))[[1]]
    dat[j,13]<-quantile(dAICrc,c(0.025,0.975))[[2]]
    dat[j,14]<-median(dAICrc)
    dat[j,15]<-quantile(dAICrc,c(0.25,0.75))[[1]]
    dat[j,16]<-quantile(dAICrc,c(0.25,0.75))[[2]]
    dat[j,17]<-if (dat[j,12] < 0 && dat[j,13] < 0 ) {"smaller"} else if (dat[j,12] > 0 && dat[j,13] > 0) {"bigger"} else {"NS"}

}
write.table(dat,paste("toPlotFig2_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_COUNTS-n-more.txt",sep=""))

# CLADES 
dat<-data.frame(matrix(NA, nrow = length(cladeNames), ncol = 13),row.names=cladeNames)
colnames(dat)<-c("timeCon","timeVar","%timeCon","richness","GrThan1_TC","GrThan1_TV","GrThan1_NA", "GrThan2_TC","GrThan2_TV","GrThan2_NA", "mean_dAICrc","low95_dAICrc","up95_dAICrc")
for (j in 1:length(cladeRes_sum)){
    res<-cladeRes_sum[[j]]["numTimesBEST",]
    GrThan1<-cladeRes_bestModelsGrThan1[[j]]
    GrThan2<-cladeRes_bestModelsGrThan2[[j]]
    dAICrc<-cladeRes_deltaAIC_RCs[[j]]
    dat[j,1]<-sum(res[,1:2])
    dat[j,2]<-sum(res[,3:10])
    dat[j,3]<-round(sum(res[,1:2])/100,digits=2)
    dat[j,4]<-length(cladesDR[which(cladesDR$clade==cladeNames[j]),"tiplabel"])
    dat[j,5]<-sum(GrThan1[which(GrThan1$Var1=="BCST" | GrThan1$Var1=="BCSTDCST"),"Freq"])
    NA1<-GrThan1[which(GrThan1$Var1=="NA"),"Freq"]
    dat[j,7]<-if(length(NA1)>0) {NA1} else {0}
    dat[j,6]<-100-(dat[j,5]+dat[j,7])
    dat[j,8]<-sum(GrThan2[which(GrThan2$Var1=="BCST" | GrThan2$Var1=="BCSTDCST"),"Freq"])
    NA2<-GrThan2[which(GrThan2$Var1=="NA"),"Freq"]
    dat[j,10]<-if(length(NA2)>0) {NA2} else {0}
    dat[j,9]<-100-(dat[j,8]+dat[j,10])
    dat[j,11]<-mean(dAICrc)
    dat[j,12]<-quantile(dAICrc,c(0.025,0.975))[[1]]
    dat[j,13]<-quantile(dAICrc,c(0.025,0.975))[[2]]
}
write.table(dat,paste("allClades_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_COUNTS-n-more.txt",sep=""))
#write.table(dat,paste("allRedos_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_COUNTS-n-more.txt",sep=""))


# ORDERS 
dat<-data.frame(matrix(NA, nrow = length(ordNames), ncol = 13),row.names=ordNames)
colnames(dat)<-c("timeCon","timeVar","%timeCon","richness","GrThan1_TC","GrThan1_TV","GrThan1_NA", "GrThan2_TC","GrThan2_TV","GrThan2_NA", "mean_dAICrc","low95_dAICrc","up95_dAICrc")
for (j in 1:length(ordRes_sum)){
    res<-ordRes_sum[[j]]["numTimesBEST",]
    GrThan1<-ordRes_bestModelsGrThan1[[j]]
    GrThan2<-ordRes_bestModelsGrThan2[[j]]
    dAICrc<-ordRes_deltaAIC_RCs[[j]]
    dat[j,1]<-sum(res[,1:2])
    dat[j,2]<-sum(res[,3:10])
    dat[j,3]<-round(sum(res[,1:2])/100,digits=2)
    dat[j,4]<-length(cladesDR[which(cladesDR$ord==ordNames[j]),"tiplabel"])
    dat[j,5]<-sum(GrThan1[which(GrThan1$Var1=="BCST" | GrThan1$Var1=="BCSTDCST"),"Freq"])
    NA1<-GrThan1[which(GrThan1$Var1=="NA"),"Freq"]
    dat[j,7]<-if(length(NA1)>0) {NA1} else {0}
    dat[j,6]<-100-(dat[j,5]+dat[j,7])
    dat[j,8]<-sum(GrThan2[which(GrThan2$Var1=="BCST" | GrThan2$Var1=="BCSTDCST"),"Freq"])
    NA2<-GrThan2[which(GrThan2$Var1=="NA"),"Freq"]
    dat[j,10]<-if(length(NA2)>0) {NA2} else {0}
    dat[j,9]<-100-(dat[j,8]+dat[j,10])
    dat[j,11]<-mean(dAICrc)
    dat[j,12]<-quantile(dAICrc,c(0.025,0.975))[[1]]
    dat[j,13]<-quantile(dAICrc,c(0.025,0.975))[[2]]
}
write.table(dat,paste("allOrds_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_COUNTS-n-more.txt",sep=""))


# PCS 
dat<-data.frame(matrix(NA, nrow = length(patchNames), ncol = 13),row.names=patchNames)
colnames(dat)<-c("timeCon","timeVar","%timeCon","richness","GrThan1_TC","GrThan1_TV","GrThan1_NA", "GrThan2_TC","GrThan2_TV","GrThan2_NA", "mean_dAICrc","low95_dAICrc","up95_dAICrc")
for (j in 1:length(patchRes_sum)){
    res<-patchRes_sum[[j]]["numTimesBEST",]
    GrThan1<-patchRes_bestModelsGrThan1[[j]]
    GrThan2<-patchRes_bestModelsGrThan2[[j]]
    dAICrc<-patchRes_deltaAIC_RCs[[j]]
    dat[j,1]<-sum(res[,1:2])
    dat[j,2]<-sum(res[,3:10])
    dat[j,3]<-round(sum(res[,1:2])/100,digits=2)
    dat[j,4]<-length(cladesDR[which(cladesDR$patch==patchNames[j]),"tiplabel"])
    dat[j,5]<-sum(GrThan1[which(GrThan1$Var1=="BCST" | GrThan1$Var1=="BCSTDCST"),"Freq"])
    NA1<-GrThan1[which(GrThan1$Var1=="NA"),"Freq"]
    dat[j,7]<-if(length(NA1)>0) {NA1} else {0}
    dat[j,6]<-100-(dat[j,5]+dat[j,7])
    dat[j,8]<-sum(GrThan2[which(GrThan2$Var1=="BCST" | GrThan2$Var1=="BCSTDCST"),"Freq"])
    NA2<-GrThan2[which(GrThan2$Var1=="NA"),"Freq"]
    dat[j,10]<-if(length(NA2)>0) {NA2} else {0}
    dat[j,9]<-100-(dat[j,8]+dat[j,10])
    dat[j,11]<-mean(dAICrc)
    dat[j,12]<-quantile(dAICrc,c(0.025,0.975))[[1]]
    dat[j,13]<-quantile(dAICrc,c(0.025,0.975))[[2]]
}
write.table(dat,paste("allPatches_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_COUNTS-n-more.txt",sep=""))


#SUBCLADES
dat<-data.frame(matrix(NA, nrow = length(mamSubsNames), ncol = 13),row.names=mamSubsNames)
colnames(dat)<-c("timeCon","timeVar","%timeCon","richness","GrThan1_TC","GrThan1_TV","GrThan1_NA", "GrThan2_TC","GrThan2_TV","GrThan2_NA", "mean_dAICrc","low95_dAICrc","up95_dAICrc")
for (j in 1:length(mamSubsRes_sum)){
    res<-mamSubsRes_sum[[j]]["numTimesBEST",]
    GrThan1<-mamSubsRes_bestModelsGrThan1[[j]]
    GrThan2<-mamSubsRes_bestModelsGrThan2[[j]]
    dAICrc<-mamSubsRes_deltaAIC_RCs[[j]]
    dat[j,1]<-sum(res[,1:2])
    dat[j,2]<-sum(res[,3:10])
    dat[j,3]<-round(sum(res[,1:2])/100,digits=2)
    dat[j,4]<-length(cladesDR[which(cladesDR$mamSubs==mamSubsNames[j]),"tiplabel"])
    dat[j,5]<-sum(GrThan1[which(GrThan1$Var1=="BCST" | GrThan1$Var1=="BCSTDCST"),"Freq"])
    NA1<-GrThan1[which(GrThan1$Var1=="NA"),"Freq"]
    dat[j,7]<-if(length(NA1)>0) {NA1} else {0}
    dat[j,6]<-100-(dat[j,5]+dat[j,7])
    dat[j,8]<-sum(GrThan2[which(GrThan2$Var1=="BCST" | GrThan2$Var1=="BCSTDCST"),"Freq"])
    NA2<-GrThan2[which(GrThan2$Var1=="NA"),"Freq"]
    dat[j,10]<-if(length(NA2)>0) {NA2} else {0}
    dat[j,9]<-100-(dat[j,8]+dat[j,10])
    dat[j,11]<-mean(dAICrc)
    dat[j,12]<-quantile(dAICrc,c(0.025,0.975))[[1]]
    dat[j,13]<-quantile(dAICrc,c(0.025,0.975))[[2]]
}
write.table(dat,paste("allMamSubs_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_COUNTS-n-more.txt",sep=""))


# Load back in just the summary tables...
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_divModels_100_final")
i=100
finalSum_clade<-read.table(paste("allClades_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_COUNTS-n-more.txt",sep=""))
finalSum_ord<-read.table(paste("allOrds_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_COUNTS-n-more.txt",sep=""))
finalSum_PC<-read.table(paste("allPatches_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_COUNTS-n-more.txt",sep=""))
finalSum_SUB<-read.table(paste("allMamSubs_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_COUNTS-n-more.txt",sep=""))

finALL<-rbind(finalSum_SUB,finalSum_clade,finalSum_ord,finalSum_PC)

write.table(finALL,paste("allCOMBINED_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_COUNTS-n-more.txt",sep=""))

## RE-LOAD for a RE-PLOT of the model results for Fig2 LTTs
library(plotrix);library(viridis)
setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_100trees/trees_1to100_brokenOut/_NDexp_divModels_100_final")
i=100
bbone<-"NDexp"
#finALL<-read.table(paste("allCOMBINED_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_COUNTS-n-more.txt",sep=""))
finALL<-read.table(paste("toPlotFig2_",bbone,"_tree",i,"_results_Morlon_SUMMARIZED_100trees_COUNTS-n-more.txt",sep=""))

finALL[which(finALL$richness >10),]

AFROSORICIDA, CHIROPTERA, DIPROTODONTIA, MACROSCELIDEA << Equivocal


# OK, awesome, got these.
# make ALL the pie charts at once...
allClades<-rownames(finALL)

#pdf(file="LTT_pieCharts_ofMLdivModels_all86clades_NDexp_blue_2cat.pdf", onefile=TRUE, width=20,height=20)
pdf(file="LTT_pieCharts_ofMLdivModels_REDOS_NDexp_blue_2cat.pdf", onefile=TRUE, width=20,height=20)
op <- par(mfrow = c(9,10),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1, lwd=2)

for (j in 1:length(allClades)){
    percent<-finALL[allClades[j],"%timeCon"]#"X.timeCon"]
    pie(x=c(percent,1-percent),labels=NA,col=c("black","deepskyblue2"))
    mtext(allClades[j],cex=0.6)
}
dev.off()


#pdf(file="LTT_pieCharts_ofMLdivModels_all86clades_NDexp_blue_3cat_grThan1.pdf", onefile=TRUE, width=20,height=20)
pdf(file="LTT_pieCharts_ofMLdivModels_REDOS_NDexp_blue_3cat_grThan1.pdf", onefile=TRUE, width=20,height=20)
op <- par(mfrow = c(9,10),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1, lwd=2)

for (j in 1:length(allClades)){
    perCON<-finALL[allClades[j],"GrThan1_TC"]/100
    perVAR<-finALL[allClades[j],"GrThan1_TV"]/100
    perNA<-finALL[allClades[j],"GrThan1_NA"]/100
    pie(x=c(perCON,perVAR,perNA),labels=NA,col=c("black","deepskyblue2","white"))
    mtext(allClades[j],cex=0.6)
}
dev.off()

#pdf(file="LTT_pieCharts_ofMLdivModels_all86clades_NDexp_blue_3cat_grThan2.pdf", onefile=TRUE, width=20,height=20)
pdf(file="LTT_pieCharts_ofMLdivModels_REDOS_NDexp_blue_3cat_grThan2.pdf", onefile=TRUE, width=20,height=20)
op <- par(mfrow = c(9,10),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1, lwd=2)

for (j in 1:length(allClades)){
    perCON<-finALL[allClades[j],"GrThan2_TC"]/100
    perVAR<-finALL[allClades[j],"GrThan2_TV"]/100
    perNA<-finALL[allClades[j],"GrThan2_NA"]/100
    pie(x=c(perCON,perVAR,perNA),labels=NA,col=c("black","deepskyblue2","white"))
    mtext(allClades[j],cex=0.6)
}
dev.off()

########
# NOW get the 95% CI of each of the delta AIC rcs plotted together...
# plot a SERIES of CIs together... all of them

#orderToPlot<-c("DIPROTODONTIA", "DIDELPHIMORPHIA", "DASYUROMORPHIA", "PERAMELEMORPHIA", "PAUCITUBERCULATA", "RODENTIA", "MouseRelated", "SquirrelRelated", "GuineaPigRelated", "PRIMATES", "OldWorldMonkey", "NewWorldMonkey", "Strepsirrhini", "LAGOMORPHA", "SCANDENTIA", "CHIROPTERA", "EmballonuridRelated", "VespertilionidRelated", "Yinpterochiroptera", "PhyllostomidRelated", "EULIPOTYPHLA", "CETARTIODACTYLA", "TerrestrialCetartios", "Cetaceans", "CARNIVORA","CatRelated","DogRelated", "PERISSODACTYLA", "PHOLIDOTA")
#orderToPlot<-redoNames
#orderToPlot<-c("RODENTIA", "MouseRelated", "SquirrelRelated", "GuineaPigRelated", "CHIROPTERA", "Yinpterochiroptera", "Yangochiroptera", "EULIPOTYPHLA", "PRIMATES", "OldWorldMonkey2", "NewWorldMonkey2", "Strepsirrhini", "CETARTIODACTYLA", "Ruminantia", "Cetaceans", "AFROSORICIDA", "SCANDENTIA", "MACROSCELIDEA", "PILOSA", "CARNIVORA","CatRelated","DogRelated", "LAGOMORPHA", "PERISSODACTYLA", "CINGULATA", "PHOLIDOTA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA")
orderToPlot<-c("MouseRelated", "SquirrelRelated", "GuineaPigRelated", "Yinpterochiroptera", "Yangochiroptera", "Soricidae", "Talpidae","Erinaceidae","OldWorldMonkey2", "NewWorldMonkey2", "Strepsirrhini", "Ruminantia", "Cetaceans", "AFROSORICIDA")#, "SCANDENTIA", "MACROSCELIDEA", "PILOSA", "CARNIVORA","CatRelated","DogRelated", "LAGOMORPHA", "PERISSODACTYLA", "CINGULATA", "PHOLIDOTA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA")

allToPlot<-finALL[orderToPlot,]
allToPlot[which(allToPlot$X.timeCon >= 0.95),] # RC favored signif -- 5 clades, richnesses are 21,12,8,5,5
    #           GrThan2_TC GrThan2_TV GrThan2_NA mean_dAICrc low95_dAICrc
    #PILOSA             65          0         35    -2.17057    -2.924525
    #CINGULATA          45          0         55    -1.62382    -2.453000
    #PHOLIDOTA          60          0         40    -2.30030    -3.721675
    #HYRACOIDEA         98          0          2    -5.58128    -6.655525
    #SIRENIA           100          0          0    -5.72266    -6.625850

allToPlot[which(allToPlot$X.timeCon <= 0.05),] # RV favored signif
meanRV<-allToPlot[which(allToPlot$X.timeCon <= 0.05),"mean_dAICrc"] # RV favored signif
    # 14 clades

allToPlot[which(allToPlot$X.timeCon > 0.05 & allToPlot$X.timeCon < 0.95),] # Equivocal
    # 10 clades, including some speciose ones -- all bats (and sep for Yin and Yango), GP-rel clade, Strepsirrhini, NW monkeys, 

marsups<-c("DIPROTODONTIA", "DIDELPHIMORPHIA", "DASYUROMORPHIA", "PERAMELEMORPHIA", "PAUCITUBERCULATA")
marsupDat<-finALL[marsups,]
    # DIDELPHIMORPHIA: -1.13170

    #RC signif and > 10 sp: -2.17, -1.62, -1.13
    meanRC<-c(-2.17, -1.62, -1.13)
    mean(meanRC) # [1] -1.64


###
#orderToPlot<-c("RODENTIA", "MouseRelated", "SquirrelRelated", "GuineaPigRelated", "CHIROPTERA", "Yinpterochiroptera", "Yangochiroptera", "EULIPOTYPHLA", "PRIMATES", "OldWorldMonkey2", "NewWorldMonkey2", "Strepsirrhini", "CETARTIODACTYLA", "Ruminantia", "Cetaceans", "AFROSORICIDA", "SCANDENTIA", "MACROSCELIDEA", "PILOSA", "CARNIVORA","CatRelated","DogRelated", "LAGOMORPHA", "PERISSODACTYLA", "CINGULATA", "PHOLIDOTA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA")
orderToPlot<-c("MouseRelated", "SquirrelRelated", "GuineaPigRelated", "Yinpterochiroptera", "Yangochiroptera", "Soricidae", "Talpidae","Erinaceidae","OldWorldMonkey2", "NewWorldMonkey2", "Strepsirrhini", "Ruminantia", "Cetaceans", "AFROSORICIDA")#, "SCANDENTIA", "MACROSCELIDEA", "PILOSA", "CARNIVORA","CatRelated","DogRelated", "LAGOMORPHA", "PERISSODACTYLA", "CINGULATA", "PHOLIDOTA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA")

allToPlot<-finALL[orderToPlot,]

means<-finALL[orderToPlot,"mean_dAICrc"]
low95s<-finALL[orderToPlot,"low95_dAICrc"]
up95s<-finALL[orderToPlot,"up95_dAICrc"]

#pdf(file="LTT_deltaAICrc_95dists_forUpdatedLTT.pdf", onefile=TRUE, width=3,height=4.2)
pdf(file="LTT_deltaAICrc_95dists_forUpdatedLTT_finalFor6ords.pdf", onefile=TRUE, width=3,height=4.2)
op <- par(oma = c(2,2,2,2) + 0.1,
          mar = c(0,0,0,0) + 0.1, lwd=2)

plotCI(x=rev(means),y=1:length(orderToPlot), ui=rev(up95s),li=rev(low95s), cex=0.8, err="x",sfrac=0, xlim=c(-55,55),bty="n",xaxt="n",yaxt="n",xlab="delta AICrc", col="white",ylab="",lwd=2, pch=NA)
axis(side=1,at=c(-40,-20,0,20,40),labels=TRUE,cex.axis=0.7)
text(x=rep(-65,length(orderToPlot)),y=1:length(orderToPlot),labels=rev(orderToPlot), pos=4,cex=0.5)
for(i in 1:length(orderToPlot)){
    abline(h=i,lty=1, lwd=0.5, col=grey(0.6, alpha=0.5))
}
abline(v=0,lty=2, lwd=2, col=rgb(1,0,0, alpha=0.5)) #rgb(1,0,0,alpha=0.5))#
plotCI(add=TRUE,x=rev(means),y=1:length(orderToPlot), ui=rev(up95s),li=rev(low95s), cex=0.8, err="x",sfrac=0, xlim=c(-55,55),bty="n",xaxt="n",yaxt="n",xlab="delta AICrc", ylab="",lwd=2, pch=20)
dev.off()

# ADD THE TIP DR STUFF HERE TOO.

bbone<- "NDexp" # "FBD"
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
head(cladesDR)

#orderToPlot<-c("RODENTIA", "MouseRelated", "SquirrelRelated", "GuineaPigRelated", "CHIROPTERA", "Yinpterochiroptera", "Yangochiroptera", "EULIPOTYPHLA", "PRIMATES", "OldWorldMonkey2", "NewWorldMonkey2", "Strepsirrhini", "CETARTIODACTYLA", "Ruminantia", "Cetaceans", "AFROSORICIDA", "SCANDENTIA", "MACROSCELIDEA", "PILOSA", "CARNIVORA","CatRelated","DogRelated", "LAGOMORPHA", "PERISSODACTYLA", "PHOLIDOTA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA")
orderToPlot_ords<-c("RODENTIA", "CHIROPTERA", "EULIPOTYPHLA", "PRIMATES", "CETARTIODACTYLA", "AFROSORICIDA", "SCANDENTIA", "MACROSCELIDEA", "PILOSA", "CARNIVORA","LAGOMORPHA", "PERISSODACTYLA", "CINGULATA","PHOLIDOTA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA")
orderToPlot_clades<-c("Mouse-related", "Squirrel-related", "Guinea_pig-related", "Yinpterochiroptera", "Yangochiroptera", "Catarrhini", "Playtrrhini", "Strepsirrhini", "Ruminantia", "Whippomorpha", "Feliformes","Caniformes")

#orderToPlot<-c("RODENTIA", "Mouse-related", "Squirrel-related", "Guinea_pig-related", "CHIROPTERA", "Yinpterochiroptera", "Yangochiroptera", "EULIPOTYPHLA", "PRIMATES", "Catarrhini", "Playtrrhini", "Strepsirrhini", "CETARTIODACTYLA", "Ruminantia", "Whippomorpha", "AFROSORICIDA", "SCANDENTIA", "MACROSCELIDEA", "PILOSA", "CARNIVORA","Feliformes","Caniformes", "LAGOMORPHA", "PERISSODACTYLA", "CINGULATA", "PHOLIDOTA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA")
orderToPlot<-c("Mouse-related", "Squirrel-related", "Guinea_pig-related", "Yinpterochiroptera", "Yangochiroptera", "Soricidae", "Talpidae","Erinaceidae","Catarrhini", "Platyrrhini", "Strepsirrhini", "Ruminantia", "Whippomorpha", "AFROSORICIDA")#, "SCANDENTIA", "MACROSCELIDEA", "PILOSA", "CARNIVORA","CatRelated","DogRelated", "LAGOMORPHA", "PERISSODACTYLA", "CINGULATA", "PHOLIDOTA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA")

cols<-rep("black",length(orderToPlot))
#cols[match(orderToPlot_clades,orderToPlot)]<-grey(0.5)

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

cladeDR_allMams<-cladesDR[,"harmMeans"]

meanDat<-rep(NA,length(cladeDR_all)); lowDat<-rep(NA,length(cladeDR_all)); upDat<-rep(NA,length(cladeDR_all)); pVals<-rep(NA,length(cladeDR_all)); overlapZero<-rep(NA,length(cladeDR_all))
for (i in 1:length(cladeDR_all)){
    dat<-cladeDR_all[[i]]
        test<-wilcox.test(x=cladeDR_allMams, y=dat, alternative="two.sided", conf.int=TRUE)
        meanDat[i]<-test$estimate[[1]]
        lowDat[i]<-test$conf.int[1]
        upDat[i]<-test$conf.int[2]
        pVals[i]<-test$p.value
        overlapZero[i]<-if (test$conf.int[1] < 0 && test$conf.int[2] < 0 ) {"bigger"} else if (test$conf.int[1] > 0 && test$conf.int[2] > 0) {"smaller"} else {"NS"}
}
cladeDR_diffMedians<-cbind.data.frame(meanDat,lowDat,upDat,pVals,overlapZero)
rownames(cladeDR_diffMedians)<-orderToPlot

median<-rep(NA,length(cladeDR_all)); low95<-rep(NA,length(cladeDR_all)); up95<-rep(NA,length(cladeDR_all)); low50<-rep(NA,length(cladeDR_all)); up50<-rep(NA,length(cladeDR_all))
for (i in 1:length(cladeDR_all)){
    dat<-cladeDR_all[[i]]
        int95<-quantile(dat,c(0.025,0.975))
        int50<-quantile(dat,c(0.25,0.75))
        median[i]<-median(dat)
        low95[i]<-int95[[1]]
        up95[i]<-int95[[2]]
        low50[i]<-int50[[1]]
        up50[i]<-int50[[2]]
}
cladeDR_int95<-cbind.data.frame(median,low95,up95,low50,up50)
rownames(cladeDR_int95)<-orderToPlot

# Get COLORS
cols50<-c(viridis(5,alpha=1)[1:4],"gold","darkgoldenrod3")

ID<-c(1,1,1,2,2,3,3,3,4,4,4,5,5,6)

cols50All<-cols50[ID]

last2ords<-c("gold","darkgoldenrod3")
last2ordsAlpha<-c()
for(q in 1:length(last2ords)){
    asRGB<-col2rgb(last2ords[q])/255
    asAlpha<-rgb(red=asRGB[[1]],green=asRGB[[2]],blue=asRGB[[3]],alpha=0.5)
    last2ordsAlpha[q]<-asAlpha
}

cols95<-c(viridis(5,alpha=0.5)[1:4],last2ordsAlpha)

cols95All<-cols95[ID]

ptCols_DR<-cols50All
ptCols_DR[which(cladeDR_diffMedians$overlapZero=="NS")]<-"white"

ptCols_ML<-cols50All
ptCols_ML[which(finALL$diffZero=="NS")]<-"white"

orderToPlot<-c("Mouse-related", "Squirrel-related", "Guinea_pig-related", "Yinpterochiroptera", "Yangochiroptera", "Soricidae", "Talpidae","Erinaceidae","Catarrhini", "Platyrrhini", "Strepsirrhini", "Ruminantia", "Whippomorpha", "AFROSORICIDA")#, "SCANDENTIA", "MACROSCELIDEA", "PILOSA", "CARNIVORA","CatRelated","DogRelated", "LAGOMORPHA", "PERISSODACTYLA", "CINGULATA", "PHOLIDOTA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA")

#pdf(file="LTT_deltaAICrc_95dists_forUpdatedLTT_dualWith_tipDR.pdf", onefile=TRUE, width=4,height=5.3)
pdf(file="LTT_deltaAICrc_95dists_forUpdatedLTT_dualWith_tipDR_finalFor6ords.pdf", onefile=TRUE, width=4,height=5.3)
op <- par(mfrow = c(1,2),
            oma = c(2,2,2,2) + 0.1,
            mar = c(0,0,0,0) + 0.1, lwd=2)

# TipDR statistical differences
#plotCI(x=rev(cladeDR_diffMedians$meanDat),y=1:length(orderToPlot), ui=rev(cladeDR_diffMedians$upDat),li=rev(cladeDR_diffMedians$lowDat), cex=0.8, err="x",sfrac=0, xlim=c(-0.4,0.4),bty="n",xaxt="n",yaxt="n",xlab="delta AICrc", col="white",ylab="",lwd=2, pch=NA)
#axis(side=1,at=c(-0.4,-0.2,0,0.2,0.4),labels=c(-0.4,NA,0,NA,0.4),cex.axis=0.7, padj=-1)
#text(x=rep(-65,length(orderToPlot)),y=1:length(orderToPlot),labels=rev(orderToPlot), pos=4,cex=0.5)
#for(i in 1:length(orderToPlot)){
#    abline(h=i,lty=1, lwd=0.5, col=grey(0.6, alpha=0.5))
#}
#abline(v=0,lty=2, lwd=2, col=rgb(1,0,0, alpha=0.5)) #rgb(1,0,0,alpha=0.5))#
#plotCI(add=TRUE,x=rev(cladeDR_diffMedians$meanDat),y=1:length(orderToPlot), ui=rev(cladeDR_diffMedians$upDat),li=rev(cladeDR_diffMedians$lowDat), cex=0.8, err="x",sfrac=0, xlim=c(-55,55),bty="n",xaxt="n",yaxt="n",xlab="delta AICrc", col="black",ylab="",lwd=2, pch=20)

# TipDR 95% CI per clade
plotCI(x=rev(cladeDR_int95$median),y=1:length(orderToPlot), ui=rev(cladeDR_int95$up95),li=rev(cladeDR_int95$low95), cex=0.8, err="x",sfrac=0, xlim=c(0,0.8),bty="n",xaxt="n",yaxt="n",xlab="delta AICrc", col="white",ylab="",lwd=2, pch=NA)
axis(side=1,at=c(0,0.2,0.4,0.6,0.8),labels=c(0,NA,0.4,NA,0.8),cex.axis=0.7, padj=-1)
int95_mam<-quantile(cladeDR_allMams,c(0.025,0.975))
rect(xleft=int95_mam[[1]],xright=int95_mam[[2]],ybottom=0,ytop=30, col=grey(0.95), border=NA)
for(i in 1:length(orderToPlot)){
    abline(h=i,lty=1, lwd=0.5, col=grey(0.6, alpha=0.5))
}
abline(v=median(cladeDR_allMams),lty=1, lwd=2, col="dark grey") #rgb(1,0,0,alpha=0.5))#

plotCI(add=TRUE, x=rev(cladeDR_int95$median),y=1:length(orderToPlot), ui=rev(cladeDR_int95$up95),li=rev(cladeDR_int95$low95),cex=0.8, err="x", sfrac=0, xlim=c(-55,55),bty="n",xaxt="n",yaxt="n",xlab="delta AICrc", col=rev(cols95All),scol=rev(cols95All),ylab="",lwd=3, pch=NA)
plotCI(add=TRUE, x=rev(cladeDR_int95$median),y=1:length(orderToPlot), ui=rev(cladeDR_int95$up50),li=rev(cladeDR_int95$low50),cex=0.8, err="x", sfrac=0, xlim=c(-55,55),bty="n",xaxt="n",yaxt="n",xlab="delta AICrc", col=rev(ptCols_DR),scol=rev(cols50All),ylab="",lwd=6, pch=21)


# Models for RC/RV

plotCI(x=rev(means),y=1:length(orderToPlot), ui=rev(up95s),li=rev(low95s), cex=0.8, err="x",sfrac=0, xlim=c(-20,55),bty="n",xaxt="n",yaxt="n",xlab="delta AICrc", col="white",ylab="",lwd=2, pch=NA)
axis(side=1,at=c(-20,0,20,40),labels=c(-20,0,NA,40),cex.axis=0.7, padj=-1)
for(i in 1:length(orderToPlot)){
    abline(h=i,lty=1, lwd=0.5, col=grey(0.6, alpha=0.5))
}
abline(v=0,lty=1, lwd=2, col="salmon") #rgb(1,0,0,alpha=0.5))#
text(x=rep(-30,length(orderToPlot)),y=1:length(orderToPlot),labels=rev(orderToPlot), pos=4,cex=0.5)

#plotCI(add=TRUE,x=rev(means),y=1:length(orderToPlot), ui=rev(up95s),li=rev(low95s), cex=0.8, err="x",sfrac=0, xlim=c(-55,55),bty="n",xaxt="n",yaxt="n",xlab="delta AICrc", ylab="",col=rev(ptCols_ML),scol=rev(cols50All),lwd=4, pch=21)

plotCI(add=TRUE, x=rev(finALL$median_dAICrc),y=1:length(orderToPlot), ui=rev(finALL$up95_dAICrc),li=rev(finALL$low95_dAICrc), cex=0.8, err="x",sfrac=0, xlim=c(-55,55),bty="n",xaxt="n",yaxt="n",xlab="delta AICrc", ylab="", col=rev(cols95All),scol=rev(cols95All),ylab="",lwd=3, pch=NA)
plotCI(add=TRUE, x=rev(finALL$median_dAICrc),y=1:length(orderToPlot), ui=rev(finALL$up50_dAICrc),li=rev(finALL$low50_dAICrc), cex=0.8, err="x",sfrac=0, xlim=c(-55,55),bty="n",xaxt="n",yaxt="n",xlab="delta AICrc", ylab="",col=rev(ptCols_ML),scol=rev(cols50All),ylab="",lwd=6, pch=21)


dev.off()





########
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



# =================================
# GLOBAL TREE - 10 samples, FBD
# TreePar RATE SHIFTS
###
# Load results back in (FBD)
###
#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_10trees")
#load("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample10_1_nexus_complete_results_TreePar.Rdata")

# Set up to run TREEPAR-- on 10 trees, and then subsetting those trees...
# ready the trees and subsets first.
library(ape); library(picante); library(phytools); library(geiger); library(DDD); library(pspline);library(TreePar)
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/TreePar_analyses_MamPhy")
#source("run_diversification_analyses.R")
#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_10trees")

library(foreach);library(doSNOW)
cl = makeCluster(10, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees=10

foreach(i=1:ntrees, .packages=c('ape', 'picante', 'phytools','geiger','DDD','pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar% {

# which backbone?
bbone<- "NDexp" #"FBD" # 

# load tree
mamPhy<-ladderize(drop.tip(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample10_",i,".tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 

# tip categories
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
head(cladesDR)

#ordNames<-names(which(table(cladesDR$ord) > 25))[c(1:4,8:11)] # [8:11] # excluding "DASYUROMORPHIA"  "DIDELPHIMORPHIA" "DIPROTODONTIA" bc crowns are ~ 20, 21, 30.1 Ma < 30 Ma needed for 7 splits of 5 Ma...
#setNames<-ordNames
#setCat<-"ord"
higherNames<-names(which(table(cladesDR$higher) > 25))
setNames<-higherNames
setCat<-"higher"

# BREAKOUTS
mamPhy_cladeSets<-vector("list",length(setNames))
for (j in 1:length(setNames)){
    setTipNames<-cladesDR[which(as.vector(cladesDR[,setCat])==setNames[j]),"tiplabel"]
    toDrop<-setdiff(mamPhy$tip.label,as.vector(setTipNames))
    mamPhy_cladeSets[[j]]<-drop.tip(mamPhy,toDrop)
    }

#lapply(lapply(mamPhy_cladeSets,getx),max)

setNames<-"globalTree"

for(q in 1:length(setNames)){
#foreach(q=1:length(setNames), .packages=c('ape', 'picante', 'phytools','geiger','DDD','pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar% { 

#tree_file<-mamPhy_cladeSets[[q]]
tree_file<-mamPhy

#num_shifts=7; grid=1; sampling_fraction=1; number_of_trees=1 # looking at 7 shifts over -- now 1 Ma-- 5 Ma increments >> means clades must > 30 Ma !
grid=1; sampling_fraction=1; number_of_trees=1 # looking for 1 mass extinction time only...  1 Ma grid

        brTimes<-getx(tree_file)

#        BD_shifts<-bd.shifts.optim(x=brTimes,c(sampling_fraction,rep(1,num_shifts)),grid=grid,start=0,end=ceiling(max(brTimes)),yule=FALSE,ME=FALSE,all=FALSE,posdiv=FALSE)
        BD_shifts<-bd.shifts.optim(x=brTimes,sampling=c(sampling_fraction,0.05),grid=grid,start=0,end=ceiling(max(brTimes)),yule=FALSE,ME=TRUE,all=FALSE,posdiv=FALSE)
    
        res<-BD_shifts[[2]]

############# RESULTS ###########################################
results<-data.frame(matrix(NA,(1+num_shifts),5+(3*(1+num_shifts)-1) ))
colnames(results)<-c("Model","NP","logL","AICc","deltaAICc",paste("DivRate",1:(num_shifts+1),sep=""),paste("Turnover",1:(num_shifts+1),sep=""),paste("ShiftTime",1:num_shifts,sep=""))
results[,1]<-c("NoShiftTime",paste(1:num_shifts,"ShiftTime",sep=""))

    #NP
        for(j in 1:(num_shifts+1)){
            results[j,2]<-length(res[[j]])-1
        }

    #logL
        for(j in 1:(num_shifts+1)){
            results[j,3]<-round(-res[[j]][1],4)
        }
        
    #AICc
        for(j in 1:(num_shifts+1)){
            NP<-length(res[[j]])-1
            results[j,4]<-round((2*(-round(-res[[j]][1],4))+2*NP+(2*NP*(NP+1))/((length(brTimes)+1)-NP-1)),3)
        }

    #deltaAICc
        for(j in 1:(num_shifts+1)){
            MIN<-min(results[,4])
            results[j,5]<-results[j,4]-MIN
        }

    #DivRates
        for(j in 1:(num_shifts+1)){ # down rows first
            for(k in 0:(num_shifts)){ # going across columns
                if( (3+(3*k)) > length(res[[j]])) { break } else { 
                results[j,(k+6)]<-round(res[[j]][j+2+k],4)
                }
            }
        }

    #Turnovers
        for(j in 1:(num_shifts+1)){ # down rows first
            for(k in 0:(num_shifts)){ # going across columns
                if( (3+(3*k)) > length(res[[j]])) { break } else { 
                results[j,(k+(num_shifts+7))]<-round(res[[j]][k+2],4)
                }
            }
        }

    #ShiftTimes
        for(j in 1:(num_shifts+1)){ # down rows first           
            for(k in 0:(num_shifts-1)){ # going across columns
                if( (4+(3*k)) > length(res[[j]])) { break } else { 
                pos<-length(res[[j]])-(j-2)
                results[j,(k+((2*num_shifts)+8))]<-round(res[[j]][pos+k],4)
                }
            }
        }

write.table(results, file=paste("TreePar_shifts",num_shifts,"_grid",grid,"Ma_clade_",setNames[q],"_tree",i,".txt",sep=""))

} # end per clade loop

} # end 10 tree loop


######
# Load back in data to summarize. Do model averaging across 10 trees.
####
#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_10trees")
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/TreePar_results")

# which backbone?
bbone<- "NDexp" #"FBD" # 

# tip categories
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
head(cladesDR)

num_shifts=7; grid=1

ordNames<-names(which(table(cladesDR$ord) > 25))[c(1:4,8:11)] # [8:11] # excluding "DASYUROMORPHIA"  "DIDELPHIMORPHIA" "DIPROTODONTIA" bc crowns are ~ 20, 21, 30.1 Ma < 30 Ma needed for 7 splits of 5 Ma...
setNames<-ordNames
setCat<-"ord"
#higherNames<-c(names(which(table(cladesDR$higher) > 25)),"globalTree")
#setNames<-higherNames
#setCat<-"higher"

numTrees<-10

# SUMMARIZE
for(q in 1:length(setNames)){

results<-vector("list", length=numTrees)
for(i in 1:numTrees){
    results[[i]]<-read.table(file=paste("TreePar_shifts",num_shifts,"_grid",grid,"Ma_clade_",setNames[q],"_tree",i,".txt",sep=""))
}
resAll<-do.call(rbind,results)

bestModels<-resAll[which(resAll$deltaAICc==0),]
bestOf10trees<-names(sort(table(bestModels$Model),decreasing=TRUE)[1])
bestData<-resAll[which(as.vector(resAll$Model)==bestOf10trees),]

sumBestData<-rbind(sapply(bestData[,-1],mean,names=FALSE),sapply(bestData[,-1],quantile,c(0.5,0.025,0.975),na.rm=TRUE,names=FALSE))
rownames(sumBestData)<-paste(setNames[q],".",c("mean","median","low95","up95"),sep="")

write.table(sumBestData,file=paste("TreePar_shifts",num_shifts,"_grid",grid,"Ma_clade_",setNames[q],"_10trees_sumBestData.txt",sep=""))

#}

## PLOT the data
### Want a line for the mean, and then a polygon for the 95% CI...
#meanShifts<- -sumBestData[1,][21:27]
#meanRates<-sumBestData[1,][5:11]
#plot(type="s",x= meanShifts, y= meanRates)
    # >> issue here with taking the mean of a rate category that's not the same across models...

# What if just plot the 10 models...
# Best models (any shift number)

#toPlot<-bestData
toPlot<-bestModels

pdf(file=paste("TreePar_shifts",num_shifts,"_grid",grid,"Ma_clade_",setNames[q],"_10trees_Plot_10bestModels.pdf",sep=""),onefile=TRUE,height=5,width=8)
#pdf(file=paste("TreePar_shifts",num_shifts,"_grid",grid,"Ma_clade_",setNames[q],"_10trees_Plot_bestModel10trees.pdf",sep=""),onefile=TRUE,height=5,width=8)
Lwd<-2
Jit<-NULL

if(toPlot[1,"Model"]=="NoShiftTime"){
    plot(type="l",x= c(-100,0), 
        y= c(toPlot[1,][6:12][c(which(toPlot[1,][22:28]!="NA"),(length(which(toPlot[1,][22:28]!="NA"))+1))],toPlot[1,][6:12][c(which(toPlot[1,][22:28]!="NA"),(length(which(toPlot[1,][22:28]!="NA"))+1))]), 
        lwd=Lwd, ylab="divRate",xlab="Time before present (Ma)")
} else {
    Xvals<- toPlot[1,][22:28][which(toPlot[1,][22:28]!="NA")]
    plot(type="S",x= c(-Xvals,-100), 
        y= toPlot[1,][6:12][c(which(toPlot[1,][22:28]!="NA"),(length(which(toPlot[1,][22:28]!="NA"))+1))], 
        lwd=Lwd, ylab="divRate",xlab="Time before present (Ma)")
    rug(x= jitter(as.numeric(-toPlot[1,][22:28][which(toPlot[1,][22:28]!="NA")]),amount=Jit),lwd=1, side=3)
    
    }


for(j in 2:length(toPlot[,1])){
    par(new=TRUE)
    
    if(toPlot[j,"Model"]=="NoShiftTime"){
    points(type="l",x= c(-100,0), 
        y= c(toPlot[j,][6:12][c(which(toPlot[j,][22:28]!="NA"),(length(which(toPlot[j,][22:28]!="NA"))+1))],toPlot[j,][6:12][c(which(toPlot[j,][22:28]!="NA"),(length(which(toPlot[j,][22:28]!="NA"))+1))]), 
        lwd=Lwd)
    } else {

    Yvals<-which(toPlot[j,][22:28]!="NA")
    
    if(length(Yvals) < 7) {
    points(type="S",x= c(-toPlot[j,][22:28][which(toPlot[j,][22:28]!="NA")],-100), 
        y= toPlot[j,][6:12][c(Yvals,(length(Yvals)+1))], 
        lwd=Lwd)
    rug(x= jitter(as.numeric(-toPlot[j,][22:28][Yvals]),amount=Jit),lwd=1, side=3)
    } else {
    points(type="S",x= -toPlot[j,][22:28][which(toPlot[j,][22:28]!="NA")], 
        y= toPlot[j,][6:12][Yvals], 
        lwd=Lwd)
    rug(x= jitter(as.numeric(-toPlot[j,][22:28]),amount=Jit),lwd=1, side=3)
    }
    }
}
mtext(text=paste(setNames[q]," -- up to ",num_shifts," shifts, ",grid,"Ma grid"," -- 10 best models",sep=""))
#mtext(text=paste(setNames[q]," -- up to ",num_shifts," shifts, ",grid,"Ma grid"," -- best model 10 trees",sep=""))

dev.off()


}



##
# REDO- the same thing but for CoMET model-- across all clades
#####
# =================================
# GLOBAL TREE - 10 trees, ND
# CoMET model RATE SHIFTS
###
# ready the trees and subsets first.
library(ape); library(picante); library(phytools); library(geiger); library(DDD); library(pspline);library(TreePar); library(TESS); library(phangorn)
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")

#setwd("/Volumes/MercurySSD-240GB/Users/NateSSD/Vertlife_CURRENT/Diversification_analyses-Condamine/_mamPhy_global_10trees")

library(foreach);library(doSNOW)
cl = makeCluster(10, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees=10

foreach(i=1:ntrees, .packages=c('TESS','phangorn','geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

# set tree dir
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")

# which backbone?
bbone<- "NDexp" #"FBD" # 

# load tree
#mamPhy<-ladderize(drop.tip(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample10_",i,".tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 
mamPhy_orig<-ladderize(drop.tip(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample10_",i,".tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 

is.ultrametric(mamPhy_orig) ## fails

mamPhy<-nnls.tree(cophenetic(mamPhy_orig),mamPhy_orig,rooted=TRUE)
is.ultrametric(mamPhy) ## should pass

# tip categories
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
head(cladesDR)

ordNames<-names(which(table(cladesDR$ord) > 25))[c(1:4,8:11)] # [8:11] # excluding "DASYUROMORPHIA"  "DIDELPHIMORPHIA" "DIPROTODONTIA" bc crowns are ~ 20, 21, 30.1 Ma < 30 Ma needed for 7 splits of 5 Ma...
higherNames<-names(which(table(cladesDR$higher) > 25))
setNames<-c(higherNames, ordNames)

# BREAKOUTS
mamPhy_cladeSets<-vector("list",length(setNames))
for (j in 1:length(setNames)){
    if(j < 6){ setCat<-"higher"} else {setCat<-"ord"}
    setTipNames<-cladesDR[which(as.vector(cladesDR[,setCat])==setNames[j]),"tiplabel"]
    toDrop<-setdiff(mamPhy$tip.label,as.vector(setTipNames))
    mamPhy_cladeSets[[j]]<-drop.tip(mamPhy,toDrop)
    }

## Run CoMet
# CPP on Mass-Extinction Times (CoMET) model (May et al., 2015)
 
# set results dir
#setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/CoMET_analyses_MamPhy")
#####
# set parameters
samplingFraction <- 1.0 # all tips are eincluded in the full MamPhy
numExpectedMassExtinctions <- 0 # not including reference here-- make analogous to TREEPAR analyses...
numExpectedRateChanges <- 7 # seven possible shift points.
param<-"7shift-noME"

for(q in 1:length(setNames)){
#foreach(q=1:length(setNames), .packages=c('ape', 'picante', 'phytools','geiger','DDD','pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar% { 

tree_file<-mamPhy_cladeSets[[q]]
posterior_directories <- paste("MamPhy_CoMET_estimHyper_",param,"_",setNames[q],"_10trees_tree",i,"_minESS500_run",1:4,sep="")


for(j in 1:length(posterior_directories)){

tess.analysis(tree_file,
  empiricalHyperPriors = TRUE,
  samplingProbability = samplingFraction,
  numExpectedRateChanges = numExpectedRateChanges,
  numExpectedMassExtinctions = numExpectedMassExtinctions,
  estimateNumberMassExtinctions = FALSE,
  estimateNumberRateChanges = TRUE,
  #pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
  #pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
  MAX_ITERATIONS = 100000000,
  MAX_TIME = 24*60*60,
  MIN_ESS = 500,
  dir = posterior_directories[j])

}

} # close loop on 13 clades

} # close multiple trees loop (parallelized)

##
# SUMMARIZE the CoMET results...
######

library(ape); library(phytools); library(picante); library(plotrix); library(phangorn); library(phyloch)
library(TESS); library(matrixStats)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/")
source("tess.plot.output2.R")

#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_MamPhy")
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_MamPhy_7shift-noME")
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_MamPhy_subclades")
previousTree<-"MamPhy"

# set parameters
samplingFraction <- 1.0 # all tips are eincluded in the full MamPhy
numExpectedMassExtinctions <- 0 # not including reference here-- make analogous to TREEPAR analyses...
numExpectedRateChanges <- 7 # seven possible shift points.
param<-"7shift-noME"

for (q in 1:length(setNames)){

for (i in c(1:10)) {  #start multiple trees loop
posterior_directories <- paste("MamPhy_CoMET_estimHyper_",param,"_",setNames[q],"_10trees_tree",i,"_minESS500_run",1:4,sep="")

setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_MamPhy_subclades/",posterior_directories[1],sep=""))
tre<-read.nexus("input_tree.tre")
root<-max(branching.times(tre))
setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_MamPhy_subclades",sep=""))

for(j in 1:length(posterior_directories)){

assign(paste("output_",j,sep=""), tess.process.output(posterior_directories[j],
  numExpectedRateChanges = numExpectedRateChanges,
  numExpectedMassExtinctions = numExpectedMassExtinctions, 
  burnin=0.25,
  numIntervals=round((root/5),0),
  criticalBayesFactors=c(2,6,10)))
 
}

output_list <- list(output_1,output_2,output_3,output_4)

for(j in 1:length(output_list)){

pdf(file=paste("MamPhy_CoMET_estimHyper_",param,"_",setNames[q],"_10trees_tree",i,"_minESS500_run",j,"_results_5Ma.pdf",sep=""),width=8,height=8)

layout.mat <- matrix(1:6,nrow=3,ncol=2,byrow=TRUE)
layout(layout.mat)
  tess.plot.output(output_list[[j]],
  fig.types = c("speciation rates","speciation shift times","extinction rates","extinction shift times","net-diversification rates","relative-extinction rates"),#"mass extinction Bayes factors","mass extinction times"),
  las=2)
 
dev.off()

}

pdf(file=paste("MamPhy_CoMET_estimHyper_",param,"_",setNames[q],"_10trees_tree",i,"_minESS500_all4runsCombined_chainAssessment.pdf",sep=""), width=8,height=16)

layout.mat <- matrix(1:5,nrow=5,ncol=1,byrow=TRUE)
layout(layout.mat)
tess.plot.multichain.diagnostics(output_list, parameters = c("speciation rates"
    ,"speciation shift times","extinction rates","net-diversification rates","relative-extinction rates"),#"mass extinction times"), 
    las=2)
# "extinction shift times"
dev.off()


} # close multiple trees loop

} # close 13 clade loop


# Combine runs TOGETHER
##
# 4 runs of each trees...
###

for (q in 1:length(setNames)){

for (i in c(1:10)) {  #start multiple trees loop

posterior_directories <- paste("MamPhy_CoMET_estimHyper_",param,"_",setNames[q],"_10trees_tree",i,"_minESS500_run",1:4,sep="")

setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_MamPhy_subclades/",posterior_directories[1],sep=""))
tre<-read.nexus("input_tree.tre")
root<-max(branching.times(tre))
setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_MamPhy_subclades",sep=""))
binSize<-1 #5

for(j in 1:length(posterior_directories)){

assign(paste("output_",j,sep=""), tess.process.output(posterior_directories[j],
  numExpectedRateChanges = numExpectedRateChanges,
  numExpectedMassExtinctions = numExpectedMassExtinctions, 
  burnin=0.25,
  numIntervals=round((root/binSize),0),
  criticalBayesFactors=c(2,6,10)))
 
}

output_list <- list(output_1,output_2,output_3,output_4)

keysAll <- unique(unlist(lapply(output_list, names)))

join4runs<-vector("list",(length(keysAll)))
names(join4runs)<-keysAll
for(j in c(1:18)){
  key<-keysAll[j]
  if(j %in% c(1:4)){
    join4runs[[j]]<-c(output_list[[1]][[key]],output_list[[2]][[key]],output_list[[3]][[key]],output_list[[4]][[key]])
  } else {
  join4runs[[j]]<-rbind(output_list[[1]][[key]],output_list[[2]][[key]],output_list[[3]][[key]],output_list[[4]][[key]])
  }
}
join4runs[[19]]<-output_list[[1]]$tree
join4runs[[20]]<-output_list[[1]][[20]]

for(j in c(7,8,11,12,16,17,18)){
join4runs[[j]]<-colMeans(join4runs[[j]])
}

assign(paste("join4runs_tree",i,sep=""),join4runs)

pdf(file=paste("MamPhy_CoMET_estimHyper_",param,"_",setNames[q],"_10trees_tree",i,"_minESS500_results_join4runs_",binSize,"Ma.pdf",sep=""),width=8,height=8)

layout.mat <- matrix(1:6,nrow=3,ncol=2,byrow=TRUE)
layout(layout.mat)
  tess.plot.output(join4runs, #join4runs_ALL, #join4runs
  fig.types = c("speciation rates","speciation shift times","extinction rates","extinction shift times","net-diversification rates","relative-extinction rates"),#"mass extinction Bayes factors","mass extinction times"),
  las=2)
 
dev.off()

} # close multiple trees loop


###
# Combine 10 trees together

allFileNames<-paste("join4runs_tree",1:10,sep="")

output_list_ALL <- list(join4runs_tree1, join4runs_tree2, join4runs_tree3, join4runs_tree4, join4runs_tree5, join4runs_tree6, join4runs_tree7, join4runs_tree8, join4runs_tree9, join4runs_tree10)

# get root to use
treeAge_ALL<-c()
for(j in 1:length(output_list_ALL)){
  treeAge_ALL[j] <- max(branching.times(output_list_ALL[[j]]$tree))
}
treeAgeMed<-median(treeAge_ALL)
treeAgeMax<-max(treeAge_ALL)

# get bins to use
allCols<-c(); 
for(j in c(1:7,9:10)){
    allCols[j]<-dim(output_list_ALL[[j]][[5]])[2]
}
allCols # this is the number of 5 Ma bins per tree
maxBin<-max(na.omit(allCols))

intToUse<-which(allCols==maxBin)[1]
treeAge<-treeAge_ALL[intToUse] # using tree 4, 41 intervals
output_list_ALL[[intToUse]][["tree"]]$root.age<-treeAge


for(j in c(1:7,9:10)){
  colToAdd<-maxBin-allCols[j] #dim(output_list_ALL[[j]][[k]])[2]

  for(k in c(5:7,9:11,13:16)){ # the variables based on bins where NEED to even them out!
    if(k %in% c(7,11,16)){
      output_list_ALL[[j]][[k]]<-c(rep(output_list_ALL[[j]][[k]][1],colToAdd),output_list_ALL[[j]][[k]])
    } else {

      if(colToAdd==0){ addCol<-NULL } else if(colToAdd==1) { 
      rowToAdd<-dim(output_list_ALL[[j]][[k]])[1]
      addCol<-rep(output_list_ALL[[j]][[k]][1,1],rowToAdd)
  #    names(addCol)<-paste("col",dim(output_list_ALL[[j]][[k]])[2]+1,sep="")
    } else {
      rowToAdd<-dim(output_list_ALL[[j]][[k]])[1]
      addCol<-rep(output_list_ALL[[j]][[k]][1,1],rowToAdd)
      for(z in 1:(colToAdd-1)){
        addCol<-cbind(addCol,rep(output_list_ALL[[j]][[k]][1,1],rowToAdd))
      }
      colnames(addCol)<-paste("col",(allCols[j]+(1:colToAdd)),sep="")
    }
      output_list_ALL[[j]][[k]]<-cbind(addCol,output_list_ALL[[j]][[k]])
    }
  }
}

allCols2<-c(); 
for(j in c(1:7,9:10)){
#    allCols2[j]<-dim(output_list_ALL[[j]][[5]])[2]
    allCols2[j]<-length(output_list_ALL[[j]][[7]])
}
allCols2 # this is the number of 5 Ma bins per tree
maxBin2<-max(na.omit(allCols2))


Q<-18
join4runs_ALL<-vector("list",(length(keysAll)))
names(join4runs_ALL)<-keysAll
for(j in c(1:Q)){
  key<-keysAll[j]
  if(j %in% c(1:4)){
    join4runs_ALL[[j]]<-c(output_list_ALL[[1]][[key]], output_list_ALL[[2]][[key]], output_list_ALL[[3]][[key]], output_list_ALL[[4]][[key]], 
        output_list_ALL[[5]][[key]], output_list_ALL[[6]][[key]], output_list_ALL[[7]][[key]], 
        #output_list_ALL[[8]][[key]], 
        output_list_ALL[[9]][[key]], 
        output_list_ALL[[10]][[key]])
    #(join4runs_tree1[[key]], join4runs_tree2[[key]], join4runs_tree3[[key]], join4runs_tree4[[key]], join4runs_tree5[[key]], join4runs_tree6[[key]], join4runs_tree7[[key]], join4runs_tree8[[key]], join4runs_tree9[[key]], join4runs_tree10[[key]])
  } else {
  join4runs_ALL[[j]]<-rbind(output_list_ALL[[1]][[key]], output_list_ALL[[2]][[key]], output_list_ALL[[3]][[key]], output_list_ALL[[4]][[key]], 
        output_list_ALL[[5]][[key]], output_list_ALL[[6]][[key]], output_list_ALL[[7]][[key]], 
        #output_list_ALL[[8]][[key]], 
        output_list_ALL[[9]][[key]], 
        output_list_ALL[[10]][[key]])
}
}

join4runs_ALL[[19]]<-output_list_ALL[[intToUse]][["tree"]]
join4runs_ALL[[20]]<-output_list_ALL[[intToUse]][[20]] #using the deepest rooted tree to set the intervals

for(j in c(7,8,11,12,16,17,18)){
  m<-join4runs_ALL[[j]]
  #m[!is.finite(m)] <- 0
  join4runs_ALL[[j]]<-colMeans(m)
}


pdf(file=paste("MamPhy_CoMET_estimHyper_",param,"_",setNames[q],"_10trees_minESS500_results_join4runs_ALL10trees_",binSize,"Ma_treeAgeUsing.pdf",sep=""),width=8,height=8)

op <- par(oma = c(1,1,5,1) + 0.1,   #c(bottom, left, top, right)
          mar = rep(4,4) + 0.1)

layout.mat <- matrix(1:6,nrow=3,ncol=2,byrow=TRUE)
layout(layout.mat)
  tess.plot.output2(join4runs_ALL, #join4runs
  fig.types = c("speciation rates","speciation shift times","extinction rates","extinction shift times","net-diversification rates","relative-extinction rates"),#"mass extinction Bayes factors","mass extinction times"),
  las=2)

title(main=paste(previousTree," ", setNames[q],"(10 trees, 4 runs each)",sep=""), xlab = "",
      ylab = "",
      outer = TRUE, line = 1.5,cex.main=1.5,font.main=2)

dev.off()

save(join4runs_ALL,file=paste("MAMPHY_10tree_CoMET_summary_5Ma_7shiftnoME_",setNames[q],".Rda",sep=""))

} # close 13 clade loop













###########
###############
# ready the data
trees10<-read.tree("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample10.tre")
for(i in 1:length(trees10)){
    write.nexus(trees10[[i]],file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample10_",i,"_nexus.tre",sep=""))
}

# ready nodes
library(foreach);library(doSNOW)
cl = makeCluster(10, type = 'SOCK', outfile="")
registerDoSNOW(cl)

# run parallel
ntrees = 10

treeParMamPhy10 <- foreach(i=1:ntrees, .packages=c('ape', 'DDD', 'picante', 'pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar%
    run_TreePar(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample10_",i,"_nexus.tre",sep=""), sampling_fraction=1, grid=0.1, number_of_trees=1)


# GLOBAL TREE - 10 samples, FBD
# Condamine envData correlations
###
# ready nodes
library(foreach);library(doSNOW)
cl3 = makeCluster(10, type = 'SOCK', outfile="")
registerDoSNOW(cl3)

# run parallel
ntrees = 10

paleoEnvMamPhy10 <- foreach(i=1:ntrees, .packages=c('ape', 'DDD', 'picante', 'pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar% { run_PaleoEnv(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample10_",i,"_nexus.tre",sep=""), env_data_file="./PaleoEnv/PastTemperature.txt", sampling_fraction=1, number_of_trees=1)
    return(final_table_tree_file)
}


# GLOBAL TREE - all 100 samples, FBD
###
# ready the MAMPHY data, 1 tree per file...
mamFBD_100<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_nexus.trees")
for(i in 1:length(mamFBD_100)){
    write.nexus(mamFBD_100[[i]],file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_",i,"_nexus.tre",sep=""))
}

mamNDexp_100<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_nexus.trees")
for(i in 1:length(mamNDexp_100)){
    write.nexus(mamNDexp_100[[i]],file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_",i,"_nexus.tre",sep=""))
}

# BD, Morlon models
# ready nodes
library(foreach);library(doSNOW)
cl2 = makeCluster(15, type = 'SOCK', outfile="")
registerDoSNOW(cl2)

# run parallel
ntrees = 100
morlonParallel <- foreach(i=1:ntrees, .packages=c('ape', 'DDD', 'picante', 'pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar% run_Morlon_models(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_",i,"_nexus.tre",sep=""), sampling_fraction=1, number_of_trees=1)


# Condamine envData correlations
# ready nodes
library(foreach);library(doSNOW)
cl4 = makeCluster(15, type = 'SOCK', outfile="")
registerDoSNOW(cl4)

# run parallel
ntrees = 100

paleoEnvMamPhy100 <- foreach(i=1:ntrees, .packages=c('ape', 'DDD', 'picante', 'pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar%run_PaleoEnv(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_",i,"_nexus.tre",sep=""), env_data_file="./PaleoEnv/PastTemperature.txt", sampling_fraction=1, number_of_trees=1)

# PLACENTALS - all 100 samples, FBD
###
# ready the MAMPHY data, 1 tree per file...
mamFBD_100<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_nexus.trees")

# drop all MARSUP and MONOTREMEs... those patches PC1 and PC23
patchNames<-names(which(table(cladesDR$PC) > 4))

patchTipNames<-vector("list",length(patchNames))
for (i in 1:length(patchNames)){
    patchTipNames[[i]]<-cladesDR[which(cladesDR$PC==patchNames[i]),"tiplabel"]
}

toDrop<-as.vector(unlist(c(patchTipNames[11],patchTipNames[15])))
mamFBD_100_placentals<-vector("list",length(mamFBD_100))
for (i in 1:length(mamFBD_100)){
     mamFBD_100_placentals[[i]] <- drop.tip(mamFBD_100[[i]], toDrop)
}
write.nexus(mamFBD_100_placentals, file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_PLACENTALS_nexus.trees")

# Now BREAK APART to 1 file per tree...
#mamFBD_100_placentals<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_PLACENTALS_nexus.trees")
for(i in 1:length(mamFBD_100_placentals)){
    write.nexus(mamFBD_100_placentals[[i]],file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_PLACENTALS_",i,"_nexus.tre",sep=""))
}

# BD, Morlon models
# ready nodes
library(foreach);library(doSNOW)
cl2 = makeCluster(15, type = 'SOCK', outfile="")
registerDoSNOW(cl2)

# run parallel
ntrees = 100
morlonParallel <- foreach(i=1:ntrees, .packages=c('ape', 'DDD', 'picante', 'pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar% run_Morlon_models(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_PLACENTALS_",i,"_nexus.tre",sep=""), sampling_fraction=1, number_of_trees=1)










###################

#===============
# Doing this across 100 tree samples
# MEANS--- getting the MEANS as a result
# ================

# ML models of lineage diversification
#######
# ORDER AND FAMS AND GLOBAL-- fit MANY models across these
####
# Go...

# Load in data
mamFBD_100<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_nexus.trees")
mamNDexp_100<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_nexus.trees")

mamFBD_MCC<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target.tre")
mamNDexp_MCC<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target.tre")

cladesDR<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_DR-SUMMARY-expanded_FBD.txt")
head(cladesDR)
colnames(cladesDR)<-c("tiplabel","gen","fam","ord","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")

ordNames<-names(which(table(cladesDR$ord) > 4))
famNames<-names(which(table(cladesDR$fam) > 4))
patchNames<-names(which(table(cladesDR$PC) > 4))

#ordPhy<-read.tree("MamPhy_BDvr_pcsFIXED_FBD_MCC_target_27ORDERS.tre")
#famPhy<-read.tree("MamPhy_BDvr_pcsFIXED_FBD_MCC_target_162FAMILIES.tre")
#patchPhy<-read.tree("MamPhy_BDvr_pcsFIXED_FBD_MCC_target_28PATCHES.tre")

# ORDERS breakout
ordTipNames<-vector("list",length(ordNames))
for (i in 1:length(ordNames)){
	ordTipNames[[i]]<-cladesDR[which(cladesDR$ord==ordNames[i]),"tiplabel"]
}

mamFBD_100_ord_FBD<-vector("list",length(mamFBD_100))
mamFBD_100_ord_NDexp<-vector("list",length(mamNDexp_100))

for (k in 1:length(ordTipNames)){
	toDrop<-setdiff(mamFBD_100[[1]]$tip.label,ordTipNames[[k]])
	
	for (i in 1:length(mamFBD_100_ord_FBD)){
		mamFBD_100_ord_FBD[[i]]<-drop.tip(mamFBD_100[[i]],toDrop)
	}
	write.nexus(mamFBD_100_ord_FBD, file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_ORDS_",ordNames[k],".trees",sep=""))

	for (i in 1:length(mamFBD_100_ord_NDexp)){
		mamFBD_100_ord_NDexp[[i]]<-drop.tip(mamNDexp_100[[i]],toDrop)
	}
	write.nexus(mamFBD_100_ord_NDexp, file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_ORDS_",ordNames[k],".trees",sep=""))
}

mamFBD_MCC_ord<-vector("list",length(ordTipNames))
mamNDexp_MCC_ord<-vector("list",length(ordTipNames))

for (k in 1:length(ordTipNames)){
	toDrop<-setdiff(mamFBD_MCC$tip.label,ordTipNames[[k]])
	mamFBD_MCC_ord[[i]]<-drop.tip(mamFBD_MCC,toDrop)
	write.nexus(mamFBD_MCC_ord[[i]], file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_ORDS_",ordNames[k],".tre",sep=""))

	mamNDexp_MCC_ord[[i]]<-drop.tip(mamNDexp_MCC,toDrop)
	write.nexus(mamNDexp_MCC_ord[[i]], file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_ORDS_",ordNames[k],".tre",sep=""))
}

# ORDERS re-load
# for FBD
for (i in 1:length(ordNames))
    {
    assign(paste(ordNames[i],".FBD.trees",sep=""), read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_ORDS_",ordNames[i],".trees",sep=""))) 
    }

allORDS_FBD <- vector("list",length(ordNames))
for (i in 1:length(ordNames)){
	allORDS_FBD[[i]]<-get(paste(ordNames[i],".FBD.trees",sep=""))
}

# for NDexp
for (i in 1:length(ordNames))
    {
    assign(paste(ordNames[i],".ND.trees",sep=""), read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_ORDS_",ordNames[i],".trees",sep=""))) 
    }

allORDS_ND <- vector("list",length(ordNames))
for (i in 1:length(ordNames)){
	allORDS_ND[[i]]<-get(paste(ordNames[i],".ND.trees",sep=""))
}

# PATCHES breakout
patchTipNames<-vector("list",length(patchNames))
for (i in 1:length(patchNames)){
	patchTipNames[[i]]<-cladesDR[which(cladesDR$PC==patchNames[i]),"tiplabel"]
}

mamFBD_100_pc_FBD<-vector("list",length(mamFBD_100))
mamFBD_100_pc_NDexp<-vector("list",length(mamNDexp_100))

for (k in 1:length(patchTipNames)){
	toDrop<-setdiff(mamFBD_100[[1]]$tip.label,patchTipNames[[k]])
	
	for (i in 1:length(mamFBD_100_pc_FBD)){
		mamFBD_100_pc_FBD[[i]]<-drop.tip(mamFBD_100[[i]],toDrop)
	}
	write.nexus(mamFBD_100_pc_FBD, file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_PCS_",patchNames[k],".trees",sep=""))

	for (i in 1:length(mamFBD_100_pc_NDexp)){
		mamFBD_100_pc_NDexp[[i]]<-drop.tip(mamNDexp_100[[i]],toDrop)
	}
	write.nexus(mamFBD_100_pc_NDexp, file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_PCS_",patchNames[k],".trees",sep=""))
}


# PATCHES re-load
# for FBD
for (i in 1:length(patchNames))
    {
    assign(paste(patchNames[i],".FBD.trees",sep=""), read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_PCS_",patchNames[i],".trees",sep=""))) 
    }

allPCS_FBD <- vector("list",length(patchNames))
for (i in 1:length(patchNames)){
	allPCS_FBD[[i]]<-get(paste(patchNames[i],".FBD.trees",sep=""))
}

# for NDexp
for (i in 1:length(patchNames))
    {
    assign(paste(patchNames[i],".ND.trees",sep=""), read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_PCS_",patchNames[i],".trees",sep=""))) 
    }

allPCS_ND <- vector("list",length(patchNames))
for (i in 1:length(patchNames)){
	allPCS_ND[[i]]<-get(paste(patchNames[i],".ND.trees",sep=""))
}


#######
# CONDAMINE code for running 
# MORLON: all Yule, BD con + 4 BD timeDep Expo + 4 BD timeDepLin
# COND: BD con + 4 BD tempDep Expo + 4 BD tempDep Lin
# ETTIENE: 6 DD models (linear and expo combos)
# STADLER: TreePar model of RATE SHIFTS (up to 4) through time

#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/Diversification_analyses-Condamine/")

##
# run this on LITORIA
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/")

cladesDR<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_DR-SUMMARY-expanded_FBD.txt")
head(cladesDR)
colnames(cladesDR)<-c("tiplabel","gen","fam","ord","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")

ordNames<-names(which(table(cladesDR$ord) > 4))
famNames<-names(which(table(cladesDR$fam) > 4))
patchNames<-names(which(table(cladesDR$PC) > 4))

setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/Diversification_analyses-Condamine/")
source("run_diversification_analyses.R")

# Parallelized
##
# stopCluster(cl) 
# PATCHES -- MORLON MODELS -- FBD
###
library(foreach);library(doSNOW)
cl = makeCluster(16, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees = length(patchNames)

morlonParallel <- foreach(i=1:ntrees, .packages=c('ape', 'DDD', 'picante', 'pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar% run_Morlon_models(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_PCS_", patchNames[i] ,".trees",sep=""), sampling_fraction=1, number_of_trees=100)

# ORDERS -- MORLON MODELS -- FBD, in parallel starting from the 9th order... = Eulipto.
###
library(foreach);library(doSNOW)
cl = makeCluster(15, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees = length(ordNames)

morlonParallel <- foreach(i=9:ntrees, .packages=c('ape', 'DDD', 'picante', 'pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar% run_Morlon_models(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_ORDS_", ordNames[i] ,".trees",sep=""), sampling_fraction=1, number_of_trees=100)

# ORDERS -- MORLON MODELS -- ND, in parallel starting from the 9th order... = Eulipto.
###
library(foreach);library(doSNOW)
cl = makeCluster(15, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees = length(ordNames)

morlonParallel <- foreach(i=9:ntrees, .packages=c('ape', 'DDD', 'picante', 'pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar% run_Morlon_models(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_ORDS_", ordNames[i] ,".trees",sep=""), sampling_fraction=1, number_of_trees=100)



# GLOBAL TREE - 10 samples, FBD
# TreePar RATE SHIFTS
###
# ready the data
trees10<-read.tree("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample10.tre")
for(i in 1:length(trees10)){
	write.nexus(trees10[[i]],file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample10_",i,"_nexus.tre",sep=""))
}

# ready nodes
library(foreach);library(doSNOW)
cl = makeCluster(10, type = 'SOCK', outfile="")
registerDoSNOW(cl)

# run parallel
ntrees = 10

treeParMamPhy10 <- foreach(i=1:ntrees, .packages=c('ape', 'DDD', 'picante', 'pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar%
	run_TreePar(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample10_",i,"_nexus.tre",sep=""), sampling_fraction=1, grid=0.1, number_of_trees=1)


# GLOBAL TREE - 10 samples, FBD
# Condamine envData correlations
###
# ready nodes
library(foreach);library(doSNOW)
cl3 = makeCluster(10, type = 'SOCK', outfile="")
registerDoSNOW(cl3)

# run parallel
ntrees = 10

paleoEnvMamPhy10 <- foreach(i=1:ntrees, .packages=c('ape', 'DDD', 'picante', 'pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar% { run_PaleoEnv(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample10_",i,"_nexus.tre",sep=""), env_data_file="./PaleoEnv/PastTemperature.txt", sampling_fraction=1, number_of_trees=1)
	return(final_table_tree_file)
}


# GLOBAL TREE - all 100 samples, FBD
###
# ready the MAMPHY data, 1 tree per file...
mamFBD_100<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_nexus.trees")
for(i in 1:length(mamFBD_100)){
    write.nexus(mamFBD_100[[i]],file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_",i,"_nexus.tre",sep=""))
}

mamNDexp_100<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_nexus.trees")
for(i in 1:length(mamNDexp_100)){
    write.nexus(mamNDexp_100[[i]],file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_",i,"_nexus.tre",sep=""))
}

# BD, Morlon models
# ready nodes
library(foreach);library(doSNOW)
cl2 = makeCluster(15, type = 'SOCK', outfile="")
registerDoSNOW(cl2)

# run parallel
ntrees = 100
morlonParallel <- foreach(i=1:ntrees, .packages=c('ape', 'DDD', 'picante', 'pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar% run_Morlon_models(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_",i,"_nexus.tre",sep=""), sampling_fraction=1, number_of_trees=1)


# Condamine envData correlations
# ready nodes
library(foreach);library(doSNOW)
cl4 = makeCluster(15, type = 'SOCK', outfile="")
registerDoSNOW(cl4)

# run parallel
ntrees = 100

paleoEnvMamPhy100 <- foreach(i=1:ntrees, .packages=c('ape', 'DDD', 'picante', 'pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar%run_PaleoEnv(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_",i,"_nexus.tre",sep=""), env_data_file="./PaleoEnv/PastTemperature.txt", sampling_fraction=1, number_of_trees=1)

# PLACENTALS - all 100 samples, FBD
###
# ready the MAMPHY data, 1 tree per file...
mamFBD_100<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_nexus.trees")

# drop all MARSUP and MONOTREMEs... those patches PC1 and PC23
patchNames<-names(which(table(cladesDR$PC) > 4))

patchTipNames<-vector("list",length(patchNames))
for (i in 1:length(patchNames)){
    patchTipNames[[i]]<-cladesDR[which(cladesDR$PC==patchNames[i]),"tiplabel"]
}

toDrop<-as.vector(unlist(c(patchTipNames[11],patchTipNames[15])))
mamFBD_100_placentals<-vector("list",length(mamFBD_100))
for (i in 1:length(mamFBD_100)){
     mamFBD_100_placentals[[i]] <- drop.tip(mamFBD_100[[i]], toDrop)
}
write.nexus(mamFBD_100_placentals, file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_PLACENTALS_nexus.trees")

# Now BREAK APART to 1 file per tree...
#mamFBD_100_placentals<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_PLACENTALS_nexus.trees")
for(i in 1:length(mamFBD_100_placentals)){
    write.nexus(mamFBD_100_placentals[[i]],file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_PLACENTALS_",i,"_nexus.tre",sep=""))
}

# BD, Morlon models
# ready nodes
library(foreach);library(doSNOW)
cl2 = makeCluster(15, type = 'SOCK', outfile="")
registerDoSNOW(cl2)

# run parallel
ntrees = 100
morlonParallel <- foreach(i=1:ntrees, .packages=c('ape', 'DDD', 'picante', 'pspline', 'TreePar'), .combine=cbind, .verbose=TRUE) %dopar% run_Morlon_models(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_PLACENTALS_",i,"_nexus.tre",sep=""), sampling_fraction=1, number_of_trees=1)








###
# Visualizing the model results, once have a best ML model:
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/Diversification_analyses-Condamine")


# Cetartiodactyla

BTimeVar_LIN

means.Lambda	means.AlphaTime
0.318064	-0.0109548

f.lamb<- y[1]+(y[2]*x)}


plot_dtt(fit.bd, tot_time, N0)

# BTimeVar LIN
f.lamb<-function(x,y){y[1]+(y[2]*x)}
f.mu<-function(x,y){0}
lamb_par<-c(0.01,0.001)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

phyloi <- read.nexus(file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_ORDS_CETARTIODACTYLA.tre")
tot_time<-max(node.age(phyloi)$ages)
f<-1
cond="crown"

BTimeVar_LIN<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond,dt=1e-3)

print(BTimeVar_LIN)

class(res) <- "fit.bd"

lamb_par :
[1]  0.295730781 -0.008963531

pdf(file="Cetartio_test_DTT.pdf", width=5, height=5, onefile=TRUE)

plot_dtt(BTimeVar_LIN, tot_time, N0=length(phyloi$tip.label))

dev.off()


plot_dtt<-function (fit.bd, tot_time, N0) 
{
    if (!inherits(fit.bd, "fit.bd")) 
        stop("object \"fit.bd\" is not of class \"fit.bd\"")
    t <- seq(0, tot_time, length.out = 100)
    if ("f.mu" %in% attributes(fit.bd)$names) {
        r <- function(t) {
            -fit.bd$f.lamb(t) + fit.bd$f.mu(t)
        }
        R <- function(s) {
            .Integrate(Vectorize(r), 0, s)
        }
        N <- N0 * exp(Vectorize(R)(t))
        dev.new()
        plot(-t, N, type = "l", xlab = "time", ylab = "Number of species", 
            main = "Diversity Through Time")
    }
    else {
        r <- function(t) {
            -fit.bd$f.lamb(t)
        }
        R <- function(s) {
            .Integrate(Vectorize(r), 0, s)
        }
        N <- N0 * exp(Vectorize(R)(t))
        dev.new()
        plot(-t, N, type = "l", xlab = "time", ylab = "Number of species", 
            main = "Diversity Through Time")
    }
}


plotdtt <- function (fit.bd, tot_time, N0, col=1, add=FALSE, div.time, xlim, ylim)
{
    if (!inherits(fit.bd, "fit.bd"))
        stop("object \"fit.bd\" is not of class \"fit.bd\"")
    t <- seq(tot_time-div.time, tot_time, 0.01)
    if ("f.mu" %in% attributes(fit.bd)$names) {
        r <- function(t) {
            -fit.bd$f.lamb(t) + fit.bd$f.mu(t)
        }
        R <- function(s) {
            RPANDA:::.Integrate(Vectorize(r), 0, s)
        }
        N <- N0 * exp(Vectorize(R)(t))
                                        #dev.new()
        if(add==FALSE)
            {
        plot(-t, N, type = "l", xlab = "time", ylab = "Number of species",
             main = "Diversity Through Time", col=col, xlim=xlim, ylim=ylim)
    }
        else
            {
                lines(-t, N, type = "l", xlab = "time", ylab = "Number of species",
                     main = "Diversity Through Time", col=col, xlim=xlim, ylim=ylim)
            }
    }
    else {
        r <- function(t) {
            -fit.bd$f.lamb(t)
        }
        R <- function(s) {
            RPANDA:::.Integrate(Vectorize(r), 0, s)
        }
        N <- N0 * exp(Vectorize(R)(t))
                                        #dev.new()
        if(add==FALSE)
            {
        plot(-t, N, type = "l", xlab = "time", ylab = "Number of species",
             main = "Diversity Through Time",col=col, xlim=xlim, ylim=ylim)
    }
        else
            {
                lines(-t, N, type = "l", xlab = "time", ylab = "Number of species",
                     main = "Diversity Through Time",col=col, xlim=xlim, ylim=ylim)
            }
    }
}
plotdtt(results$balaenopteridae$bcstdcst,div.times[1],N0=Ntip(balaenopteridae.tree),xlim=c(-max(div.times),0),ylim=c(0,150),div.time=div.times[1])
plotdtt(results$delphinidae$bvardcst,div.times[2],N0=Ntip(delphinidae.tree),col=6,add=TRUE,xlim=c(-max(div.times),0),ylim=c(0,150),div.time=div.times[2])
plotdtt(results$phocoenidae$bcstdcst,div.times[3],N0=Ntip(phocoenidae.tree),col="goldenrod",add=TRUE,xlim=c(-max(div.times),0),ylim=c(0,150),div.time=div.times[3])
plotdtt(results$ziphidae$bcstdcst,div.times[4],N0=Ntip(ziphidae.tree),col=4,add=TRUE,xlim=c(-max(div.times),0),ylim=c(0,150),div.time=div.times[4])
plotdtt(results$othercetaceans$bcstdvar,div.times[5],N0=Ntip(othercetaceans.tree),col="darkred",add=TRUE,xlim=c(-max(div.times),0),ylim=c(0,150),div.time=div.times[5])
legend("topleft",legend=c("Balaenopteridae","Delphinidae","Phocoenidae","Ziphidae","Other Cetaceans"),text.col=c(1,6,"goldenrod",4,"darkred"))



#########
# Non-parallel
##
# Morlon-- FBD-- orders, then Global 
for (i in 1:length(ordNames)){
	run_Morlon_models(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_ORDS_", ordNames[i] ,".trees",sep=""), sampling_fraction=1, number_of_trees=100)
}
run_Morlon_models(tree_file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_nexus.trees", sampling_fraction=1, number_of_trees=100)

# Morlon-- NDexp-- orders, then Global 
for (i in 1:length(ordNames)){
	run_Morlon_models(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_ORDS_", ordNames[i] ,".trees",sep=""), sampling_fraction=1, number_of_trees=100)
}
run_Morlon_models(tree_file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_nexus.trees", sampling_fraction=1, number_of_trees=100)

# PaleoTEMP-- FBD-- orders, then Global 
for (i in 1:length(ordNames)){
	run_PaleoEnv(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_ORDS_", ordNames[i] ,".trees",sep=""), env_data_file="PastTemperature.txt", sampling_fraction=1, number_of_trees=100)
}
run_PaleoEnv(tree_file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_nexus.trees", env_data_file="PastTemperature.txt", sampling_fraction=1, number_of_trees=100)

# PaleoTEMP-- NDexp-- orders, then Global 
for (i in 1:length(ordNames)){
	run_PaleoEnv(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_ORDS_", ordNames[i] ,".trees",sep=""), env_data_file="PastTemperature.txt", sampling_fraction=1, number_of_trees=100)
}
run_PaleoEnv(tree_file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_nexus.trees", env_data_file="PastTemperature.txt", sampling_fraction=1, number_of_trees=100)


# TREEPAR-- FBD-- orders, then Global
for (i in 1:length(ordNames)){
	run_TreePar(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_ORDS_", ordNames[i] ,".trees",sep=""), sampling_fraction=1, grid=0.1, number_of_trees=100)
}
run_TreePar(tree_file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_.trees", sampling_fraction=1, grid=0.1, number_of_trees=100)

# TREEPAR-- NDexp-- orders, then Global
for (i in 1:length(ordNames)){
	run_TreePar(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_ORDS_", ordNames[i] ,".trees",sep=""), sampling_fraction=1, grid=0.1, number_of_trees=100)
}
run_TreePar(tree_file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_.trees", sampling_fraction=1, grid=0.1, number_of_trees=100)


## DDD models-- FBD-- orders, then Global -- ONLY 1 tree (does not scale across posterior)
for (i in 1:length(ordNames)){
	tt<-read.nexus(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_ORDS_", ordNames[i] ,".tre",sep=""))
	Ntips<-length(tt$tip.label)	
	run_DDD(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_ORDS_", ordNames[i] ,".tre",sep=""), total_richness=Ntips, number_of_trees=1)
}
	tt<-read.nexus(file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_5911sp.tre")
	Ntips<-length(tt$tip.label)	
	run_DDD(tree_file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_5911sp.tre", total_richness=Ntips, number_of_trees=1)


## DDD models-- NDexp-- orders, then Global -- ONLY 1 tree (does not scale across posterior)
for (i in 1:length(ordNames)){
	tt<-read.nexus(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_ORDS_", ordNames[i] ,".tre",sep=""))
	Ntips<-length(tt$tip.label)	
	run_DDD(tree_file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_ORDS_", ordNames[i] ,".tre",sep=""), total_richness=Ntips, number_of_trees=1)
}
	tt<-read.nexus(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp.tre")
	Ntips<-length(tt$tip.label)	
	run_DDD(tree_file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp.tre", total_richness=Ntips, number_of_trees=1)




# TESTS
run_TreePar(tree_file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_ORDS_AFROSORICIDA.trees", sampling_fraction=1, grid=0.1, number_of_trees=100)

	tt<-read.nexus(file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target.tre")
	write.nexus(drop.tip(tt,"_Anolis_carolinensis"), file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_5911sp.tre")

	tt<-read.nexus(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target.tre")
	write.nexus(drop.tip(tt,"_Anolis_carolinensis"), file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp.tre")












####
res60Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees60Ma.txt", header=TRUE)
res50Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees50Ma.txt", header=TRUE)
res40Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees40Ma.txt", header=TRUE)
res30Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees30Ma.txt", header=TRUE)
res20Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees20Ma.txt", header=TRUE)
res10Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees10Ma.txt", header=TRUE)
res5Ma <-read.table("MamPhy_5911sp_FBD_cladeLevel_DRstats_trees5Ma.txt", header=TRUE)

results<-list(res5Ma, res10Ma, res20Ma, res30Ma, res40Ma, res50Ma, res60Ma)

allCladeSetNames<-c("5Ma","10Ma","20Ma","30Ma","40Ma","50Ma","60Ma")

for (i in 1:length(allCladeSetNames))
    {
    assign(paste("trees",allCladeSetNames[[i]],sep=""), read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_slice",allCladeSetNames[[i]],"_newick.trees",sep=""))) 
    }

allCladeSets<-list(trees5Ma,trees10Ma,trees20Ma,trees30Ma,trees40Ma,trees50Ma,trees60Ma)

slicePhys<-read.tree("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target_slicePhy-5-10-20-30-40-50-60.trees")

###
# MY code, tests...
####
# Yule (empirical), BD (sp + ext, net div), DD 
###
ymle = function(tree) (.subset2(tree,2)-1L)/sum(.subset2(tree,4)) # this take the # of number of nodes in a tree / sum of branch lengths.

# Loop through trees and record the results
ymle(allCladeSets[[1]][[1]])

##
# BD, using RPANDA
##
#fit_bd(phylo, tot_time, f.lamb, f.mu, lamb_par, mu_par, f = 1,
#            meth = "Nelder-Mead", cst.lamb = FALSE, cst.mu = FALSE,
#            expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE,
#            dt=0, cond = "crown")
##
# altered function for multi-plotting...
plot_fit_bd<-function (fit.bd, tot_time) 
{
    if (!inherits(fit.bd, "fit.bd")) 
        stop("object \"fit.bd\" is not of class \"fit.bd\"")
    t <- seq(0, tot_time, length.out = 100)
    plot(-t, fit.bd$f.lamb(t), type = "l", xlab = "time", ylab = "speciation rate", 
        main = "Fitted speciation rate")
    if ("f.mu" %in% attributes(fit.bd)$names) {
        plot(-t, fit.bd$f.mu(t), type = "l", xlab = "time", ylab = "extinction rate", 
            main = "Fitted extinction rate")
        r <- function(t) {
            fit.bd$f.lamb(t) - fit.bd$f.mu(t)
        }
        plot(-t, r(t), type = "l", xlab = "time", ylab = "net diversification rate", 
            main = "Fitted net diversification rate")
    }
    else {
        plot(-t, fit.bd$f.lamb(t), type = "l", xlab = "time", 
            ylab = "net diversification rate", main = "Fitted net diversification rate")
    }
}

plot_fit_env<-function (fit.env, env_data, tot_time) 
{
    if (!inherits(fit.env, "fit.env")) 
        stop("object is not of class \"fit.env\"")
    t <- seq(0, tot_time, length.out = 100)
    dev.new()
    plot(-t, fit.env$f.lamb(t), type = "l", xlab = "time", ylab = "speciation rate", 
        main = "Fitted speciation rate")
    df <- smooth.spline(env_data[, 1], env_data[, 2])$df
    spline_result <- sm.spline(env_data[, 1], env_data[, 2], 
        df = df)
    env_func <- function(t) {
        predict(spline_result, t)
    }
    dev.new()
    plot(env_func(t), fit.env$f.lamb(t), type = "l", xlab = "Environmental data", 
        ylab = "speciation rate", main = "Fitted speciation rate")
    if ("f.mu" %in% attributes(fit.env)$names) {
        dev.new()
        plot(-t, fit.env$f.mu(t), type = "l", xlab = "time", 
            ylab = "extinction rate", main = "Fitted extinction rate")
        plot(env_func(t), fit.env$f.mu(t), type = "l", xlab = "Environmental data", 
            ylab = "extinction rate", main = "Fitted extinction rate")
        r <- function(t) {
            fit.env$f.lamb(t) - fit.env$f.mu(t)
        }
        dev.new()
        plot(-t, r(t), type = "l", xlab = "time", ylab = "net diversification rate", 
            main = "Fitted net diversification rate")
        dev.new()
        plot(env_func(t), r(t), type = "l", xlab = "Environmental data", 
            ylab = "net diversification rate", main = "Fitted net diversification rate")
    }
    else {
        dev.new()
        plot(-t, fit.env$f.lamb(t), type = "l", xlab = "time", 
            ylab = "net diversification rate", main = "Fitted net diversification rate")
        dev.new()
        plot(env_func(t), fit.env$f.lamb(t), type = "l", xlab = "Environmental data", 
            ylab = "net diversification rate", main = "Fitted net diversification rate")
    }
}

#####
# BD constant rates
##
tree<-allCladeSets[[7]][[6]]
tot_time<-max(node.age(tree)$ages)
rate_conT<-function(t,y){y[1]} # BOTH constant rate through time
f.lamb<-rate_conT
f.mu<-rate_conT
lamb_par<-c(0.01)
mu_par<-c(0.005)

bd<-fit_bd(phylo=tree, tot_time, f.lamb, f.mu, lamb_par, mu_par, f = 1, meth = "Nelder-Mead", cst.lamb = TRUE, cst.mu = TRUE, expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE, dt=0, cond = "crown")

bd$lamb_par
bd$mu_par
bd$aicc
quartz(width=10, height=3.5)
par(mfrow=c(1,3))
plot_fit_bd(bd,tot_time)


#####
# BD rates vary LINEARLY with time
rate_linT<-function(t,y){y[1]+y[2]*t} # speciation: linear variation in rate through time
f.lamb<-rate_linT
f.mu<-rate_linT
lamb_par<-c(0.05,0.01)
mu_par<-c(0.005,0.001)

bd<-fit_bd(phylo=tree, tot_time, f.lamb, f.mu, lamb_par, mu_par, f = 1, meth = "Nelder-Mead", cst.lamb = FALSE, cst.mu = FALSE, expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE, dt=2, cond = "crown")

bd$lamb_par
bd$mu_par
bd$aicc
quartz(width=10, height=3.5)
par(mfrow=c(1,3))
plot_fit_bd(bd,tot_time)

#####
# BD rates vary EXPONENTIALLY with time
rate_expT<-function(t,y){y[1]*exp(y[2]*t)} # speciation: exponential variation in rate through time
f.lamb<-rate_expT
f.mu<-rate_expT
lamb_par<-c(0.05,0.01)
mu_par<-c(0.005,0.001)

bd<-fit_bd(phylo=tree, tot_time, f.lamb, f.mu, lamb_par, mu_par, f = 1, meth = "Nelder-Mead", cst.lamb = FALSE, cst.mu = FALSE, expo.lamb = TRUE, expo.mu = TRUE, fix.mu = FALSE, dt=2, cond = "crown")

bd$lamb_par
bd$mu_par
bd$aicc
quartz(width=10, height=3.5)
par(mfrow=c(1,3))
plot_fit_bd(bd,tot_time)



# BD coalescent equilibrium model -- models 1 and 2 from Morlon et al. 2010
##
fit_coal_cst(phylo, tau0 = 1e-2, gamma = 1, cst.rate = FALSE,
                 meth = "Nelder-Mead", N0 = 0)


# BD coalescent rate var models -- models 3-6 from Morlon et al. 2010
##
fit_coal_var(phylo, lamb0 = 0.1, alpha = 1, mu0 = 0.01, beta = 0,
                 meth = "Nelder-Mead", N0 = 0, cst.lamb = FALSE, cst.mu = FALSE,
                 fix.eps = FALSE, mu.0 = FALSE, pos = TRUE)


# Environmental BD model
##
#fit_env(phylo, env_data, tot_time, f.lamb, f.mu, lamb_par, mu_par, df= NULL, f = 1, meth = "Nelder-Mead", cst.lamb = FALSE, cst.mu = FALSE, expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE, dt=0, cond = "crown")
tree<-allCladeSets[[7]][[6]]
tot_time<-max(node.age(tree)$ages)
data(InfTemp)
# plot(InfTemp) # AWESOME.
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)} # exp dependence of speciation rate on the env variable, no time dep, no ext
f.mu<-function(t,x,y){0}
lamb_par<-c(0.10,0.01)
mu_par<-c()


envBD<-fit_env(phylo=tree, env_data=InfTemp, tot_time, f.lamb, f.mu, lamb_par, mu_par, df= NULL, f = 1, meth = "Nelder-Mead", cst.lamb = FALSE, cst.mu = TRUE, expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE, dt=1e-3, cond = "crown")

quartz(width=10, height=3.5)
par(mfrow=c(1,3))

plot_fit_env(envBD,InfTemp,tot_time)





########
# HARMON code for FITTING these model simultaneously
# http://lukejharmon.github.io/ilhabela/instruction/2015/07/02/diversification-analysis-bamm-rpanda/#rpanda
# returning the AICcs also...
##
lambda.cst <- function(x,y){y}
lambda.var <- function(x,y){y[1]*exp(y[2]*x)}
mu.cst <- function(x,y){y}
mu.var <- function(x,y){y[1]*exp(y[2]*x)}

fit.multi.rpanda <- function(tree,par)
    {
        bcstdcst <- fit_bd(tree, max(branching.times(tree)), f.lamb=lambda.cst, f.mu=mu.cst, lamb_par=par[[1]][1],mu_par=par[[1]][2],cst.lamb=TRUE,cst.mu=TRUE,cond="crown",f=1,dt=1e-3)
        bvardcst <- fit_bd(tree, max(branching.times(tree)), f.lamb=lambda.var, f.mu=mu.cst, lamb_par=par[[2]][c(1,2)],mu_par=par[[2]][3],expo.lamb=TRUE,cst.mu=TRUE,cond="crown",f=1,dt=1e-3)
        bcstdvar <- fit_bd(tree, max(branching.times(tree)), f.lamb=lambda.cst, f.mu=mu.var, lamb_par=par[[3]][1],mu_par=par[[3]][c(2,3)],cst.lamb=TRUE,expo.mu=TRUE,cond="crown",f=1,dt=1e-3)
        bvardvar <- fit_bd(tree, max(branching.times(tree)), f.lamb=lambda.var, f.mu=mu.var, lamb_par=par[[4]][c(1,2)],mu_par=par[[4]][c(3,4)],expo.lamb=TRUE,expo.mu=TRUE,cond="crown",f=1,dt=1e-3)
        return(list("bcstdcst"=bcstdcst,"bvardcst"=bvardcst,"bcstdvar"=bcstdvar,"bvardvar"=bvardvar))
    }

inits <- list(c(0.4,0),c(0.4,-0.05,0),c(0.4,0.1,0.05),c(0.4,-0.05,0.1,0.05))

###
tree6<-allCladeSets[[7]][[6]]
tree7<-allCladeSets[[7]][[7]]
tree8<-allCladeSets[[7]][[8]]
tree9<-allCladeSets[[7]][[9]]
tree10<-allCladeSets[[7]][[10]]

results<-list(fit.multi.rpanda(tree6,par=inits), fit.multi.rpanda(tree7,par=inits), fit.multi.rpanda(tree8,par=inits), fit.multi.rpanda(tree9,par=inits), fit.multi.rpanda(tree10,par=inits))

aic.table <- matrix(nrow=4,ncol=length(results),NA)
for(i in 1:length(results))
    {
        for(j in 1:4)
            {
                aic.table[j,i] <- results[[i]][[j]]$aicc
            }
    }
colnames(aic.table) <- c("6","7","8","9","10")
rownames(aic.table) <- c("bcstdcst","bvardcst","bcstdvar","bvardvar")
aic.table

aic.table <- matrix(nrow=4,ncol=length(results),NA)
for(i in 1:length(results))
    {
        for(j in 1:4)
            {
                aic.table[j,i] <- results[[i]][[j]]$aicc
            }
    }

par.table <- data.frame(

	"Balaenopteridae"=c(results[[1]]$bcstdcst$lamb_par[1:2],results[[1]]$bcstdcst$mu_par[1:2]),

	"Delphinidae"=c(results[[2]]$bvardcst$lamb_par[1:2],results[[2]]$bvardcst$mu_par[1:2]),"Phocoenidae"=c(results[[3]]$bcstdcst$lamb_par[1:2],results[[3]]$bcstdcst$mu_par[1:2]),"Ziphidae"=c(results[[4]]$bcstdcst$lamb_par[1:2],results[[4]]$bcstdcst$mu_par[1:2]),"Other Cetaceans"=c(results[[5]]$bcstdvar$lamb_par[1:2],results[[5]]$bcstdvar$mu_par[1:2]))
par.table


##          Balaenopteridae Delphinidae Phocoenidae Ziphidae Other Cetaceans
## bcstdcst        56.73336    171.4150    31.79207 132.6440        119.4690
## bvardcst        58.44456    170.0338    34.13348 130.1843        117.5562
## bcstdvar        58.32515    170.6453    41.79207 127.6094        115.3723
## bvardvar        65.67244    171.3258    64.13348 133.2251        118.9789

##   Balaenopteridae   Delphinidae   Phocoenidae     Ziphidae Other.Cetaceans
## 1    7.343080e-02  1.405102e-01  1.410234e-01 9.485149e-02       0.1857764
## 2              NA  1.257642e-01            NA           NA              NA
## 3   -6.391602e-10 -1.521115e-07 -6.012454e-08 7.545982e-08       0.8313274
## 4              NA            NA            NA           NA      -0.1742352

# Function to calculate species richness in a given point in time
div.times <- c(max(branching.times(balaenopteridae.tree)),max(branching.times(delphinidae.tree)),max(branching.times(phocoenidae.tree)),max(branching.times(ziphidae.tree)),max(branching.times(othercetaceans.tree)))

# Function modified from plot_dtt from RPANDA package
plotdtt <- function (fit.bd, tot_time, N0, col=1, add=FALSE, div.time, xlim, ylim)
{
    if (!inherits(fit.bd, "fit.bd"))
        stop("object \"fit.bd\" is not of class \"fit.bd\"")
    t <- seq(tot_time-div.time, tot_time, 0.01)
    if ("f.mu" %in% attributes(fit.bd)$names) {
        r <- function(t) {
            -fit.bd$f.lamb(t) + fit.bd$f.mu(t)
        }
        R <- function(s) {
            RPANDA:::.Integrate(Vectorize(r), 0, s)
        }
        N <- N0 * exp(Vectorize(R)(t))
                                        #dev.new()
        if(add==FALSE)
            {
        plot(-t, N, type = "l", xlab = "time", ylab = "Number of species",
             main = "Diversity Through Time", col=col, xlim=xlim, ylim=ylim)
    }
        else
            {
                lines(-t, N, type = "l", xlab = "time", ylab = "Number of species",
                     main = "Diversity Through Time", col=col, xlim=xlim, ylim=ylim)
            }
    }
    else {
        r <- function(t) {
            -fit.bd$f.lamb(t)
        }
        R <- function(s) {
            RPANDA:::.Integrate(Vectorize(r), 0, s)
        }
        N <- N0 * exp(Vectorize(R)(t))
                                        #dev.new()
        if(add==FALSE)
            {
        plot(-t, N, type = "l", xlab = "time", ylab = "Number of species",
             main = "Diversity Through Time",col=col, xlim=xlim, ylim=ylim)
    }
        else
            {
                lines(-t, N, type = "l", xlab = "time", ylab = "Number of species",
                     main = "Diversity Through Time",col=col, xlim=xlim, ylim=ylim)
            }
    }
}

plotdtt(results$balaenopteridae$bcstdcst,div.times[1],N0=Ntip(balaenopteridae.tree),xlim=c(-max(div.times),0),ylim=c(0,150),div.time=div.times[1])
plotdtt(results$delphinidae$bvardcst,div.times[2],N0=Ntip(delphinidae.tree),col=6,add=TRUE,xlim=c(-max(div.times),0),ylim=c(0,150),div.time=div.times[2])
plotdtt(results$phocoenidae$bcstdcst,div.times[3],N0=Ntip(phocoenidae.tree),col="goldenrod",add=TRUE,xlim=c(-max(div.times),0),ylim=c(0,150),div.time=div.times[3])
plotdtt(results$ziphidae$bcstdcst,div.times[4],N0=Ntip(ziphidae.tree),col=4,add=TRUE,xlim=c(-max(div.times),0),ylim=c(0,150),div.time=div.times[4])
plotdtt(results$othercetaceans$bcstdvar,div.times[5],N0=Ntip(othercetaceans.tree),col="darkred",add=TRUE,xlim=c(-max(div.times),0),ylim=c(0,150),div.time=div.times[5])
legend("topleft",legend=c("Balaenopteridae","Delphinidae","Phocoenidae","Ziphidae","Other Cetaceans"),text.col=c(1,6,"goldenrod",4,"darkred"))



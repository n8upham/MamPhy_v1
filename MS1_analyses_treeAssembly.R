# ASSEMBLE the full Mammalia posterior, 10,000 trees
###

# packages
library(ape); library(phytools); library(phyloch); library(phangorn); library(paleotree)

setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors")

# LOAD FUNCTIONS
# =========
rescaleTree<-function(tree,scale)
{
	tree$edge.length<-tree$edge.length/max(nodeHeights(tree))*scale
	#tree$edge.length<-tree$edge.length/max(nodeHeights(tree)[,2])*scale
	return(tree)
}

# Code for one PALEOTREE FUNCTION:: https://github.com/dwbapst/paleotree/blob/master/R/dateNodes.R
# For use with the FBD backbone assembly only (see below)
####
dateNodes <- function (tree, rootAge = tree$root.time, labelDates = FALSE, 
    tolerance = 0.001) 
{
    if (!inherits(tree, "phylo")) {
        stop("tree must be of class 'phylo'")
    }
    if (is.null(tree$edge.length)) {
        stop("tree does not appear to have edge lengths?")
    }
    nodeRelTimes <- node.depth.edgelength(tree)
    if (is.null(rootAge)) {
        rootAge <- max(nodeRelTimes)
        message("Root age not given; treating tree as if latest tip was at modern day (time=0)")
    }
    res <- rootAge - nodeRelTimes
    if (any(res < (-tolerance))) {
        message("Warning: Some dates are negative? rootAge may be incorrectly defined or you are using a time-scaling method that warps the tree, like aba or zbla.")
    }
    if (labelDates) {
        names(res) <- sapply(Descendants(tree), function(x) paste0(sort(tree$tip.label[x]), 
            collapse = " "))
    }
    else {
        names(res) <- 1:(Ntip(tree) + Nnode(tree))
    }
    return(res)
}

# Static code

#clade<-read.table("clade_inRep_outgroup_BDvr_topoCons.txt", header=TRUE)
clade<-read.table("clade_inRep_outgroup_BDvr_topoFreeOut.txt", header=TRUE)
clades <- as.vector(clade[,"PC"])
ingroup <- as.vector(clade[,"Ingroup"])
outgroup <- as.vector(clade[,"Outgroup_MarsMono"])

#Changing code

#Backbone run, NOW with 10,000 trees
backbones <- read.nexus("postburn10k_Bbone_ND_60taxa_17calExp_BDallDNA.inMYR.trees")
#backbones <- read.nexus("postburn10k_Bbone_FBD_136taxa_topoAsZhou_FIN.inMYR.trees")

	# for ND BACKBONE
	# Drop taxa to just the INGROUP-OUTGROUP reps ...
	toKeep <- c("_Anolis_carolinensis", "Anomalurus_beecrofti_ANOMALURIDAE_RODENTIA", "Bison_bison_BOVIDAE_CETARTIODACTYLA", "Callithrix_jacchus_CALLITRICHIDAE_PRIMATES", "Calomyscus_baluchi_CALOMYSCIDAE_RODENTIA", "Castor_canadensis_CASTORIDAE_RODENTIA", "Cricetomys_gambianus_NESOMYIDAE_RODENTIA", "Cricetulus_barabensis_CRICETIDAE_RODENTIA", "Dasypus_novemcinctus_DASYPODIDAE_CINGULATA", "Eptesicus_fuscus_VESPERTILIONIDAE_CHIROPTERA", "Equus_caballus_EQUIDAE_PERISSODACTYLA", "Erethizon_dorsatum_ERETHIZONTIDAE_RODENTIA", "Erinaceus_europaeus_ERINACEIDAE_EULIPOTYPHLA", "Felis_catus_FELIDAE_CARNIVORA", "Galeopterus_variegatus_CYNOCEPHALIDAE_DERMOPTERA", "Ictidomys_tridecemlineatus_SCIURIDAE_RODENTIA", "Jaculus_jaculus_DIPODIDAE_RODENTIA", "Manis_pentadactyla_MANIDAE_PHOLIDOTA", "Ornithorhynchus_anatinus_ORNITHORHYNCHIDAE_MONOTREMATA", "Oryctolagus_cuniculus_LEPORIDAE_LAGOMORPHA", "Pteronotus_parnellii_MORMOOPIDAE_CHIROPTERA", "Pteropus_alecto_PTEROPODIDAE_CHIROPTERA", "Rattus_norvegicus_MURIDAE_RODENTIA", "Saccopteryx_bilineata_EMBALLONURIDAE_CHIROPTERA", "Spalax_ehrenbergi_SPALACIDAE_RODENTIA", "Trichechus_manatus_TRICHECHIDAE_SIRENIA", "Tupaia_belangeri_TUPAIIDAE_SCANDENTIA", "Typhlomys_cinereus_PLATACANTHOMYIDAE_RODENTIA", 
				"Vombatus_ursinus_VOMBATIDAE_DIPROTODONTIA", "Caenolestes_fuliginosus_CAENOLESTIDAE_PAUCITUBERCULATA","Tachyglossus_aculeatus_TACHYGLOSSIDAE_MONOTREMATA","Ochotona_princeps_OCHOTONIDAE_LAGOMORPHA")
	toDrop <- setdiff(backbones[[1]]$tip.label,toKeep)

	# for FBD+fossils BACKBONE
	# Drop taxa to just the INGROUP-OUTGROUP reps ...
	toKeep <- c("_Anolis_carolinensis", "Anomalurus_beecrofti_ANOMALURIDAE_RODENTIA", "Bison_bison_BOVIDAE_CETARTIODACTYLA", "Callithrix_jacchus_CALLITRICHIDAE_PRIMATES", "Calomyscus_baluchi_CALOMYSCIDAE_RODENTIA", "Castor_canadensis_CASTORIDAE_RODENTIA", "Cricetomys_gambianus_NESOMYIDAE_RODENTIA", "Cricetulus_barabensis_CRICETIDAE_RODENTIA", "Dasypus_novemcinctus_DASYPODIDAE_CINGULATA", "Eptesicus_fuscus_VESPERTILIONIDAE_CHIROPTERA", "Equus_caballus_EQUIDAE_PERISSODACTYLA", "Erethizon_dorsatum_ERETHIZONTIDAE_RODENTIA", "Erinaceus_europaeus_ERINACEIDAE_EULIPOTYPHLA", "Felis_catus_FELIDAE_CARNIVORA", "Galeopterus_variegatus_CYNOCEPHALIDAE_DERMOPTERA", "Ictidomys_tridecemlineatus_SCIURIDAE_RODENTIA", "Jaculus_jaculus_DIPODIDAE_RODENTIA", "Manis_pentadactyla_MANIDAE_PHOLIDOTA", "Ornithorhynchus_anatinus_ORNITHORHYNCHIDAE_MONOTREMATA", "Oryctolagus_cuniculus_LEPORIDAE_LAGOMORPHA", "Pteronotus_parnellii_MORMOOPIDAE_CHIROPTERA", "Pteropus_alecto_PTEROPODIDAE_CHIROPTERA", "Rattus_norvegicus_MURIDAE_RODENTIA", "Saccopteryx_bilineata_EMBALLONURIDAE_CHIROPTERA", "Spalax_ehrenbergi_SPALACIDAE_RODENTIA", "Trichechus_manatus_TRICHECHIDAE_SIRENIA", "Tupaia_belangeri_TUPAIIDAE_SCANDENTIA", "Typhlomys_cinereus_PLATACANTHOMYIDAE_RODENTIA", 
				"Vombatus_ursinus_VOMBATIDAE_DIPROTODONTIA", "Caenolestes_fuliginosus_CAENOLESTIDAE_PAUCITUBERCULATA","Tachyglossus_aculeatus_TACHYGLOSSIDAE_MONOTREMATA","Ochotona_princeps_OCHOTONIDAE_LAGOMORPHA",
				"X_Adelobasileus", "X_Aegialodon", "X_Akidolestes", "X_Albertatherium", "X_Ambondro", "X_Amphilestes", "X_Amphitherium", "X_Anchistodelphys", "X_Andinodelphys", "X_Asfaltomylos", "X_Asiatherium", "X_Asioryctes", "X_Atokatheridium", "X_Ausktribosphenos", "X_Bishops", "X_Castorocauda", "X_Cimolodontidae", "X_Daulestes", "X_Deltatheridium", "X_Didelphodon", "X_Dryolestes", "X_Eomaia", "X_Fruitafossor", "X_Gobiconodon", "X_Hadrocodium", "X_Haldanodon", "X_Haramiyavia", "X_Henkelotherium", "X_Holoclemensia", "X_Jeholodens", "X_Juramaia", "X_Kennalestes", "X_Kielantherium", "X_Kokopellia", "X_Kuehneodon", "X_Maotherium", "X_Massetognathus", "X_Mayulestes", "X_Megaconus", "X_Megazostrodon", "X_Montanalestes", "X_Morganucodon", "X_Murtoilestes", "X_Nanolestes", "X_Obdurodon", "X_Pachygenelus", "X_Pediomys", "X_Peramus", "X_Plagiaulacidae", "X_Priacodon", "X_Probainognathus", "X_Prokennalestes", "X_Pseudotribos", "X_Pucadelphys", "X_Repenomamus", "X_Rugosodon", "X_Shuotherium", "X_Sineleutherus", "X_Sinobaatar", "X_Sinoconodon", "X_Sinodelphys", "X_Spalacotherium", "X_Steropodon", "X_Sulestes", "X_Teinolophos", "X_Thomasia", "X_Thrinaxodon", "X_Tinodon", "X_Trioracodon", "X_Tritylodontidae", "X_Turgidodon", "X_Ukhaatherium", "X_Vincelestes", "X_Yanoconodon", "X_Zalambdalestes", "X_Zhangheotherium")
	toDrop <- setdiff(backbones[[1]]$tip.label,toKeep)

# drop non-IN-OUT reps, NOW for 10,000 trees
bboneDrop<-vector("list",10000)
for (i in 1:length(backbones)){
	bboneDrop[[i]]<-drop.tip(backbones[[i]],toDrop)
	#bboneDrop[[i]]$tip.label<-gsub("__", "_", bboneDrop[[i]]$tip.label) #correcting tip names
	}

#Clade trees
for (i in 1:length(clades))
    {
    #assign(paste(clades[i], ".trees", sep=""), read.nexus(paste("postburn10k_BDvr_FIXED_",clades[i], ".trees", sep=""))) 
    assign(paste(clades[i], ".trees", sep=""), read.nexus(paste("postburn10k_BDvr_topoFreeOut_",clades[i], ".trees", sep=""))) 
    }

# confirm patches load as 10k trees...
d<-vector()
for(i in 1:28){
	d[i]<-length(get(paste(clades[i], ".trees", sep="")))
	}	


#Prune out extra taxa from PC25_Dermoptera and PC28_Platacanthomyidae (there as placeholders)
PC25_DermopteraOK<-vector("list",length(PC25_Dermoptera.trees))
PC28_PlatacanthomyidaeOK<-vector("list",length(PC28_Platacanthomyidae.trees))

for (i in 1:length(PC25_Dermoptera.trees)){
	PC25_DermopteraOK[[i]]<-drop.tip(PC25_Dermoptera.trees[[i]],c("Callithrix_jacchus_CALLITRICHIDAE_PRIMATES","Gorilla_gorilla_HOMINIDAE_PRIMATES"))
	PC28_PlatacanthomyidaeOK[[i]]<-drop.tip(PC28_Platacanthomyidae.trees[[i]],c("Rattus_norvegicus_MURIDAE_RODENTIA","Spalax_ehrenbergi_SPALACIDAE_RODENTIA"))
}
PC25_Dermoptera.trees<-vector("list",length(PC25_DermopteraOK))
PC28_Platacanthomyidae.trees<-vector("list",length(PC28_PlatacanthomyidaeOK))
for (i in 1:length(PC25_Dermoptera.trees)){
	PC25_Dermoptera.trees[[i]]<-PC25_DermopteraOK[[i]]
	PC28_Platacanthomyidae.trees[[i]]<-PC28_PlatacanthomyidaeOK[[i]]
}


######
#Loop to generate 10,000 trees 

	# list to keep track of how many trees need branch shortening
	trees_cladesNeedingAdHoc<-list()

for (z in 0:9999)
{   

#Backbone
backbone <- bboneDrop[[10000 - z]]
# ND
	ages <- branching.times(backbone) 
# FBD
	#ages <- dateNodes(backbone) 
# both
nodes <- mrca(backbone)

tips <- vector()
root <- vector()
for(i in 1:length(clades)) # 1:28
{
tips[i] <- match(ingroup[i],backbone$tip.label)
root[i] <- nodes[ingroup[i],outgroup[i]]
}

# for the bats
	# find only sister-species pairs
	tree<-backbone
	sisters <- matrix(NA, ncol=2, nrow=length(tree$tip.label)) # empty matrix
	sisters[,1] <- tree$tip.label # populate first column with tip.labels
	for(i in 1:length(tree$tip.label)){
		tmp <- phytools::getSisters(tree, tree$tip.label[i], mode="label")
		if(!is.null(tmp$tips)){sisters[i,2] <- tmp$tips}
	}
	sisters <- sisters[-which(is.na(sisters[,2])),] # prune away tips that do not have a labelled tip as sister

	# find the bat sister-species pairs
	batSisters<-sisters[match(c("Pteronotus_parnellii_MORMOOPIDAE_CHIROPTERA","Eptesicus_fuscus_VESPERTILIONIDAE_CHIROPTERA", "Saccopteryx_bilineata_EMBALLONURIDAE_CHIROPTERA"),sisters[,1]),]
	batPCs<-c(16,17,18)
	batPCs_sister<-batPCs[!is.na(batSisters[,1])]
	batPCs_nonSister<-batPCs[is.na(batSisters[,1])]

	# change the root nodes for the sister pairs
	sisterRoot<-nodes[ingroup[batPCs_sister[1]],ingroup[batPCs_sister[2]]]
	root[ batPCs_sister[1] ]<-sisterRoot
	root[ batPCs_sister[2] ]<-sisterRoot
	# and the non-sister...
	root[ batPCs_nonSister ]<-root[ batPCs_nonSister ]+1
		# test:
		#cbind.data.frame(clades,root)

#Clades
clade.trees <- list()
for (i in 1:length(clades)) # 1:28
{
#Clade trees
clade.trees[[i]]<-get(paste(clades[i],".trees",sep=""))[[length(get(paste(clades[i],".trees",sep=""))) - z]]
}

#Branching times of clades, unscaled
clade.times <- lapply(clade.trees, branching.times)

#Scale of outgroup to ingroup
#clade.scales <- lapply(clade.times, function(x) {x[2]/x[1]})
clade.scales <- list()
getStemBranch<-function(x) {x[2]/x[1]}
for (i in 1:length(clades)) # 1:28
{
#Scale of outgroup to ingroup
	 if(i==28){ 
	 clade.scales[[i]]<-clade.times[[i]][1] # because only 2 species in the Platacanthomyidae PC
	 } else 
	if(i==1 | i==6 | i==23) {
	clade.scales[[i]]<- 1 # because want to attach to root (pruned of outgroup) for Marsupialia, Monotremata, & Lagomorpha
	} else {
	clade.scales[[i]]<-getStemBranch(clade.times[[i]])
	}
}

#Clade trees unscaled without outgroup
clade.trees.add <- lapply(clade.trees, drop.tip, 1)

#Root ages of trees
clade.root.ages <- ages[as.character(root)]

#Tip.ages of trees
clade.tip.ages <- setNames(backbone$edge.length[unlist(lapply(tips,function(x){which.edge(backbone,x)}))],tips)
	#rootsToConfirm<-cbind.data.frame(PC=clades,rootAge=backbone$edge.length[unlist(lapply(tips,function(x){which.edge(backbone,x)}))])

#Scaled root age of tree
clade.root.scale <- unlist(clade.scales) * clade.root.ages

#Clade trees rescaled, ready to paste
clade.trees.scaled <- list()
for (i in 1:length(clade.trees))
{
clade.trees.scaled[[i]] <- rescaleTree(clade.trees.add[[i]], clade.root.scale[i])
}


#Loop to bind trees
pasted_fulltree <- backbone 

	# vector to keep track of how many trees need branch shortening
	cladesNeedingAdHoc<-rep(NA, length(clades))

for (i in 1:length(clades))
{

#Tip edge
clade.tip <- match(ingroup[i],pasted_fulltree$tip.label)
tip.edge <- which.edge(pasted_fulltree, clade.tip)

if (clade.root.scale[i] < clade.tip.ages[i])
{
#Shorten branch
pasted_fulltree$edge.length[tip.edge] <- as.vector(clade.tip.ages[i] - clade.root.scale[i])

#Bind tree
pasted_fulltree <- bind.tree(pasted_fulltree, clade.trees.scaled[[i]], clade.tip, 0)

} else {
	# This is *just in case* the patch clade scale was greater than the bacbone scale (not applicable to our mammal trees)
	
	# Mark which clades need this shortening
	cladesNeedingAdHoc[i]<-"1"
	#Make tree fit tip age
	re.tree <- rescaleTree(clade.trees.scaled[[i]], as.vector(pasted_fulltree$edge.length[tip.edge] * 0.999999))

	#Shorten branch
	pasted_fulltree$edge.length[tip.edge] <- as.vector(pasted_fulltree$edge.length[tip.edge] * 0.000001)

	#Bind tree
	pasted_fulltree <- bind.tree(pasted_fulltree, re.tree, clade.tip, 0)
}

if(i==1){
pasted_fulltree<-drop.tip(pasted_fulltree,match("Caenolestes_fuliginosus_CAENOLESTIDAE_PAUCITUBERCULATA",pasted_fulltree$tip.label))
} else if(i==6){
pasted_fulltree<-drop.tip(pasted_fulltree,match("Ochotona_princeps_OCHOTONIDAE_LAGOMORPHA",pasted_fulltree$tip.label))
} else if(i==23){
pasted_fulltree<-drop.tip(pasted_fulltree,match("Tachyglossus_aculeatus_TACHYGLOSSIDAE_MONOTREMATA",pasted_fulltree$tip.label))
} else { NULL }

} # end clade loop

#write.tree(pasted_fulltree, "MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_topoAsZhou_all10k.trees", append = TRUE)
#write.tree(pasted_fulltree, "MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_all10k.trees", append = TRUE)

write.tree(pasted_fulltree, "MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_all10k_v2.trees", append = TRUE)
#write.tree(pasted_fulltree, "MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_all10k_v2.trees", append = TRUE)

#write.tree(pasted_fulltree, "MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_all10k_v2.trees", append = TRUE)
#write.tree(pasted_fulltree, "MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_all10k_v2.trees", append = TRUE)

trees_cladesNeedingAdHoc[[z+1]] <- cladesNeedingAdHoc
}


# Summarize the number of clades needing the adHoc 0.000001 scale adjustment::
###

#res<-as.data.frame(do.call(cbind,trees_cladesNeedingAdHoc))
#
#SUMS<-c()
#for(j in 1:length(res[,1])){
#	x<-as.numeric(na.omit( as.numeric(res[j,]) ))
#	SUMS[j]<-sum(x)
#}
#cbind.data.frame(clades,SUMS)
#                       clades SUMS
#1              PC1_Marsupials 6296
#2              PC2_Afrotheria    0
#3               PC3_Xenarthra   24
#4              PC4_Scandentia    9
#5                PC5_Primates    0
#6              PC6_Lagomorpha    0
#7           PC7_Castorimorpha    0
#8               PC8_Dipodidae    0
#9              PC9_Spalacidae    0
#10            PC10_Nesomyidae  172
#11               PC11_Muridae  672
#12            PC12_Cricetidae    0
#13       PC13_SquirrelRelated    0
#14      PC14_GuineaPigRelated    0
#15          PC15_Eulipotyphla    0
#16   PC16_PhyllostomidRelated    0
#17 PC17_VespertilionidRelated    0
#18   PC18_EmballonuridRelated    0
#19    PC19_Yinpterochiroptera    0
#20       PC20_Cetartiodactyla 2273
#21        PC21_Perissodactyla    0
#22             PC22_Carnivora    0
#23           PC23_Monotremata 5003
#24             PC24_Pholidota    0
#25            PC25_Dermoptera    0
#26       PC26_Anomaluromorpha    0
#27          PC27_Calomyscidae   31
#28     PC28_Platacanthomyidae   12







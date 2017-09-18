library(ape)
library(phytools)

#Static code
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors")

rescaleTree<-function(tree,scale)
{
tree$edge.length<-tree$edge.length/max(nodeHeights(tree)[,2])*scale
return(tree)
}

#clade <- read.csv("clades.csv")
#clades <- as.vector(clade[,1])
#ingroup <- as.vector(clade[,2])
#outgroup <- as.vector(clade[,3])

clade<-read.table("clade_inRep_outgroup_BDvr_Fixed-n-Free.txt")
clades <- as.vector(clade[,1])
ingroup <- as.vector(clade[,2])
outgroup <- as.vector(clade[,3])
PC_out <- as.vector(clade[,4])


#Changing code

#Backbone run, NOW with 10,000 trees
backbones <- read.nexus("postburn10k_Bbone_FBD_60taxa_BDallDNA.inMYR.trees")


# Drop taxa to just the INGROUP-OUTGROUP reps ...
toKeep <- c("_Anolis_carolinensis", "Anomalurus_beecrofti__ANOMALURIDAE__RODENTIA", "Bison_bison__BOVIDAE__CETARTIODACTYLA", "Callithrix_jacchus__CALLITRICHIDAE__PRIMATES", "Calomyscus_baluchi__CALOMYSCIDAE__RODENTIA", "Castor_canadensis__CASTORIDAE__RODENTIA", "Cricetomys_gambianus__NESOMYIDAE__RODENTIA", "Cricetulus_barabensis__CRICETIDAE__RODENTIA", "Dasypus_novemcinctus__DASYPODIDAE__CINGULATA", "Eptesicus_fuscus__VESPERTILIONIDAE__CHIROPTERA", "Equus_caballus__EQUIDAE__PERISSODACTYLA", "Erethizon_dorsatum__ERETHIZONTIDAE__RODENTIA", "Erinaceus_europaeus__ERINACEIDAE__EULIPOTYPHLA", "Felis_catus__FELIDAE__CARNIVORA", "Galeopterus_variegatus__CYNOCEPHALIDAE__DERMOPTERA", "Ictidomys_tridecemlineatus__SCIURIDAE__RODENTIA", "Jaculus_jaculus__DIPODIDAE__RODENTIA", "Manis_pentadactyla__MANIDAE__PHOLIDOTA", "Ornithorhynchus_anatinus__ORNITHORHYNCHIDAE__MONOTREMATA", "Oryctolagus_cuniculus__LEPORIDAE__LAGOMORPHA", "Pteronotus_parnellii__MORMOOPIDAE__CHIROPTERA", "Pteropus_alecto__PTEROPODIDAE__CHIROPTERA", "Rattus_norvegicus__MURIDAE__RODENTIA", "Saccopteryx_bilineata__EMBALLONURIDAE__CHIROPTERA", "Spalax_ehrenbergi__SPALACIDAE__RODENTIA", "Trichechus_manatus__TRICHECHIDAE__SIRENIA", "Tupaia_belangeri__TUPAIIDAE__SCANDENTIA", "Typhlomys_cinereus__PLATACANTHOMYIDAE__RODENTIA", "Vombatus_ursinus__VOMBATIDAE__DIPROTODONTIA")
toDrop <- setdiff(backbones[[1]]$tip.label,toKeep)

# drop non-IN-OUT reps, NOW for 10,000 trees
bboneDrop<-vector("list",10000)
for (i in 1:length(backbones)){
	bboneDrop[[i]]<-drop.tip(backbones[[i]],toDrop)
	bboneDrop[[i]]$tip.label<-gsub("__", "_", bboneDrop[[i]]$tip.label) #correcting tip names
	}

#Clade trees
for (i in 1:length(clades))
    {
    assign(paste(clades[i], ".trees", sep=""), read.nexus(paste("postburn10k_BDvr_FREE_",clades[i], ".trees", sep=""))) 
    }

# confirm patches load as 10k trees...
d<-vector()
for(i in 1:28){
	d[i]<-length(get(paste(clades[i], ".trees", sep="")))
	}	
d

# RE-root all clade trees
for (i in 1:length(clades)){
	temp<-vector("list",10000)
    tt<-get(paste(clades[i],".trees",sep=""))
    for (j in 1:length(get(paste(clades[i],".trees",sep="")))){
    	uu<-root(tt[[j]],PC_out[i],resolve.root=TRUE)
    	makeChronosCalib(uu, node = "root", age.min = 1, age.max = age.min, interactive = FALSE, soft.bounds = FALSE)


    	temp[[j]]<-uu
    }

    assign(paste(clades[i], ".trees", sep=""), temp) 
}







#Prune out extra taxa from PC25_Dermoptera and PC28_Platacanthomyidae (there as placeholders)
PC25_DermopteraOK<-vector("list",length(PC25_Dermoptera.trees))
PC28_PlatacanthomyidaeOK<-vector("list",length(PC28_Platacanthomyidae.trees))
for (i in 1:length(PC25_Dermoptera.trees)){
	PC25_DermopteraOK[[i]]<-drop.tip(PC25_Dermoptera.trees[[i]],c("Callithrix_jacchus_CALLITRICHIDAE_PRIMATES","Gorilla_gorilla_HOMINIDAE_PRIMATES"))
	PC28_PlatacanthomyidaeOK[[i]]<-drop.tip(PC28_Platacanthomyidae.trees[[i]],c("Spalax_ehrenbergi_SPALACIDAE_RODENTIA"))
}

PC25_Dermoptera.trees<-vector("list",length(PC25_DermopteraOK))
PC28_Platacanthomyidae.trees<-vector("list",length(PC28_PlatacanthomyidaeOK))
for (i in 1:length(PC25_Dermoptera.trees)){
	PC25_Dermoptera.trees[[i]]<-PC25_DermopteraOK[[i]]
	PC28_Platacanthomyidae.trees[[i]]<-PC28_PlatacanthomyidaeOK[[i]]
}


######
#Loop to generate 10,000 trees 

for (z in 0:9999)
{   

#Backbone
backbone <- bboneDrop[[10000 - z]]
ages <- branching.times(backbone) 
nodes <- mrca(backbone)

tips <- vector()
root <- vector()
for(i in 1:length(clades)) 
{
tips[i] <- match(ingroup[i],backbone$tip.label)
root[i] <- nodes[ingroup[i],outgroup[i]]
}

#Clades
clade.trees <- list()
for (i in 1:length(clades))
{
#Clade trees
clade.trees[[i]] <- get(paste(clades[i],".trees",sep=""))[[length(get(paste(clades[i],".trees",sep=""))) - z]]
}

#Branching times of clades, unscaled
clade.times <- lapply(clade.trees, branching.times)

#Scale of outgroup to ingroup
clade.scales <- lapply(clade.times, function(x) {x[2]/x[1]})

#Clade trees unscaled without outgroup
clade.trees.add <- lapply(clade.trees, drop.tip, PC_out)

#Root ages of trees
clade.root.ages <- ages[as.character(root)]

#Tip.ages of trees
clade.tip.ages <- setNames(backbone$edge.length[unlist(lapply(tips,function(x){which.edge(backbone,x)}))],tips)

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

#Make tree fit tip age
re.tree <- rescaleTree(clade.trees.scaled[[i]], as.vector(pasted_fulltree$edge.length[tip.edge] * 0.999999))

#Shorten branch
pasted_fulltree$edge.length[tip.edge] <- as.vector(pasted_fulltree$edge.length[tip.edge] * 0.000001)

#Bind tree
pasted_fulltree <- bind.tree(pasted_fulltree, re.tree, clade.tip, 0)
}

}

write.tree(pasted_fulltree, "MamPhy_fullPosterior_BDvr_pcsFREE_FBD_all10k.trees", append = TRUE)

}


q()

n


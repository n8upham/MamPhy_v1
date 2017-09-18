library(ape)
library(phytools)

#Static code
setwd("~/MamPhy_mtDNA_fullPosterior")

rescaleTree<-function(tree,scale)
{
tree$edge.length<-tree$edge.length/max(nodeHeights(tree)[,2])*scale
return(tree)
}

clade<-read.table("clade_inRep_outgroup.txt")
clades <- as.vector(clade[,1])
ingroup <- as.vector(clade[,2])
outgroup <- as.vector(clade[,3])


#Changing code

#Backbone run, 10,000 trees
backbones <- read.nexus("postburn10k_Bbone_ND_76taxa_17Exp_mtDNA_inMYR.trees")

# Drop taxa to just the IN-OUT reps and the PASTIS taxa...
toDrop <- c("Aepyprymnus_rufescens__POTOROIDAE__DIPROTODONTIA", "Caenolestes_fuliginosus__CAENOLESTIDAE__PAUCITUBERCULATA", "Dromiciops_gliroides__MICROBIOTHERIIDAE__MICROBIOTHERIA", "Isoodon_macrourus__PERAMELIDAE__PERAMELEMORPHIA", "Macropus_eugenii__MACROPODIDAE__DIPROTODONTIA", "Myrmecobius_fasciatus__MYRMECOBIIDAE__DASYUROMORPHIA", "Notoryctes_typhlops__NOTORYCTIDAE__NOTORYCTEMORPHIA", "Vombatus_ursinus__VOMBATIDAE__DIPROTODONTIA", "Chrysochloris_asiatica__CHRYSOCHLORIDAE__AFROSORICIDA", "Myrmecophaga_tridactyla__MYRMECOPHAGIDAE__PILOSA", "Gorilla_gorilla__HOMINIDAE__PRIMATES", "Lemur_catta__LEMURIDAE__PRIMATES", "Nycticebus_coucang__LORISIDAE__PRIMATES", "Tarsius_syrichta__TARSIIDAE__PRIMATES", "Ochotona_princeps__OCHOTONIDAE__LAGOMORPHA", "Fukomys_damarensis__BATHYERGIDAE__RODENTIA", "Solenodon_paradoxus__SOLENODONTIDAE__EULIPOTYPHLA", "Lama_glama__CAMELIDAE__CETARTIODACTYLA", "Mephitis_mephitis__MEPHITIDAE__CARNIVORA")

# for all 10,000
bboneDrop<-vector("list",10000)
for (i in 1:length(backbones)){
	bboneDrop[[i]]<-drop.tip(backbones[[i]],toDrop)
	}

#Clade trees
for (i in 1:length(clades))
    {
    assign(paste(clades[i], ".trees", sep=""), read.nexus(paste("postburn10k_YuleCR_",clades[i], ".trees", sep=""))) 
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
clade.trees.add <- lapply(clade.trees, drop.tip, 1)

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

write.tree(pasted_fulltree, "MamPhy_fullPosterior_YuleCR_all10k.trees", append = TRUE)

}

q()

n



##################
###########
#######    functions

#the CAIC function from Iain Martyn August 2011
#this produces a clade.matrix from a newick file very quickly, bypassing ape.  This speeds up ED production
#Iain says he now has ED even 10x faster than the previous 10x faster...
#Updated by Gavin Thomas April 11 2013 (readCAICnew)
#Now assigns names to rows and columns of clade.matrix and gives names to each element in the list

treeES<-readCAIC(readLines("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/9_PASTIS-ing/BayesPASTIS_EchimCaproMyo_168_115tax/combined_MCC_180Mgens_18ktrees_FigTree_newick.tre"))

#########

readCAIC = function(tree) { 

	tpc = unlist(strsplit(tree, "[\\(\\),;]"))
	tpc=tpc[grep(":",tpc)]

	# find the tip and edge labels
	tiplabels=tpc[grep(".:",tpc)]
	edgelabels=tpc[-grep(".:",tpc)]
	edgelabels=c(edgelabels,":0")

	# locate the clusters and edges
	tree <- unlist(strsplit(tree, NULL))
	x=which(tree=="(")
	y=which(tree==")")
	v=which(tree==":")
	#these are the locations of the tips
	w=setdiff(v,y+1) 

	##Pass through string from left to right locating the paired parenthesis, thus
	## allowing for easy assigment of edge to subtending tips
	#initialize objects (M is the actual clade matrix, while E is vector with associated edge weights)
	j=2 
	k=length(w)+1
	M=matrix(0,length(w)*2-1,length(w))
	E=as.vector(matrix(0,1,length(w)*2-1))
	x=c(x,y[length(y)]+1)
	
	# Main Pass
	while (length(x)>1) {	
		if (x[j+1]<y[1]) j=j+1 
		else { 
			M[k,which(x[j]<w & w<y[1])]=1
			E[k]=strsplit(edgelabels[k-length(w)],"[:]")[[1]][2]
			k=k+1
			y=y[-1]
			x=x[-j]
			j=j-1
		}
	}

	# Assign branch lengths and finished tip names to the tips
	for (i in 1:length(w)) {
		M[i,i]=1
		tmp=strsplit(tiplabels[i],"[:]")[[1]]
		E[i]=tmp[2]
		tiplabels[i]=tmp[1]
	}
	
	rownames(M) <- 1:dim(M)[1]
	colnames(M) <- tiplabels
	M=list(clade.matrix=M,edge.length=as.numeric(E),tip.label=tiplabels)
	attr(M, "class") <- "clade.matrix"
	return(M)
}


### Fast ES (Equal Splits) scores script
## from functions for faster ES.r, dated 31 Oct 2011
## older approach !
## clade.matrix is the output from readCAIC
## debugged by Ignacio-Sept/14

ES_v2 = function(CM) {
	# second way
	caicmatrix = CM$clade.matrix
	caicmatrix[caicmatrix>1] = 1
	caicmatrix[(1:dim(caicmatrix)[2]),(1:dim(caicmatrix)[2])] = diag(dim(caicmatrix)[2])
	lambda = CM$edge.length*caicmatrix
	#tmp=lambda/(2^caicmatrix)
	tmp = lambda/(caicmatrix*rowSums(caicmatrix))
	rm(caicmatrix)
	rm(lambda)
	tmp[!is.finite(tmp)] = 0
	ESS=as.matrix(colSums(tmp))
	rm(tmp)
	rownames(ESS)=CM$tip.label
	return(ESS)
}

# Say Ignacio:
#   Attached is some code I debugged some time ago from what Walter had. The first function reads the tree in newick (bypassing read.tree from ape: much faster), and then use the second function to get the tip rates.
# Just do:
#    ES_v2(readCAIC(“~path_to_newick.tre”))
# should work.

library(ape)
library(phyloch)
library(TreeSim)
library(TreePar)
library(geiger)

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/9_PASTIS-ing/BayesPASTIS_EchimCaproMyo_168_115tax")

Echi10k<-read.nexus("infile.nex.trprobs_10ktrees.nex")

Echi10k_man<-read.nexus("combined_manual-180M.t") ## Yes, these now have EDGE LENGTHS

Echi18k_man<-Echi10k_man

Echi18k_Times<-vector("list",18000)
for (i in 1:18000){
  Echi18k_Times[[i]]<-branching.times(Echi18k_man[[i]])
}

LTT.plot(c(Echi18k_man,Echi18k_Times))

mltt.plot(Echi18k_man, backward=TRUE, legend=FALSE)

Echi18k_drop<-vector("list",18000)
for (i in 1:18000){
  Echi18k_drop[[i]] <- drop.tip(Echi18k_man[[i]], drop);
} ## this works.

Echi18k_drop[[1]] 

mltt.plot(Echi18k_drop, backward=TRUE, legend=FALSE)#, ylim=c(0,114))

ltt.plot(Echi18k_drop[[1]],ylim=c(0,114))#, backward=TRUE, legend=FALSE)#, 
for (i in 2:18000){
  ltt.lines(Echi18k_drop[[i]]);
} ## thi

ltt.plot(Echi18k_man[[1]],ylim=c(0,114))#, backward=TRUE, legend=FALSE)#, 
for (i in 2:18000){
  ltt.lines(Echi18k_man[[i]]);
} ## thi



## And now the ES calculation...

ES_v2(readCAIC(“~path_to_newick.tre”))

ES_v2(readCAIC(readLines("combined_MCC_180Mgens_18ktrees_FigTree_newick.tre")))

## Just need to get into NEWICK now...

for (i in 1:1000){
  write.tree(Echi18k_man[[i]], file=paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/9_PASTIS-ing/BayesPASTIS_EchimCaproMyo_168_115tax/Echi18k_man_newicks/Echi18k_man_newick_",i,".trees",sep=""), append=FALSE)
}

ES_1000_Echi_man<-vector("list",1000)
for (i in 1:1000){
  ES_1000_Echi_man[[i]]<-ES_v2(readCAIC(readLines(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/9_PASTIS-ing/BayesPASTIS_EchimCaproMyo_168_115tax/Echi18k_man_newicks/Echi18k_man_newick_",i,".trees",sep=""))))
}

ES_1000_Echi_man_combo = apply( cbind( ES_1000_Echi_man ) , 1 , unlist )

hist(ES_1000_Echi_man_combo, breaks=1000)

library(moments)
man<-c(ES_1000_Echi_man_combo)
drop<-c(ES_1000_Echi_drop_combo)
ks.test(man, drop)
hist(dd, breaks=1000)


for (i in 1:1000){
  write.tree(Echi18k_drop[[i]], file=paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/9_PASTIS-ing/BayesPASTIS_EchimCaproMyo_168_115tax/Echi18k_man_newicks/Echi18k_drop_newick_",i,".trees",sep=""), append=FALSE)
}

ES_1000_Echi_drop<-vector("list",1000)
for (i in 1:1000){
  ES_1000_Echi_drop[[i]]<-ES_v2(readCAIC(readLines(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/9_PASTIS-ing/BayesPASTIS_EchimCaproMyo_168_115tax/Echi18k_man_newicks/Echi18k_drop_newick_",i,".trees",sep=""))))
}
ES_1000_Echi_drop_combo = apply( cbind( ES_1000_Echi_drop ) , 1 , unlist )

hist(ES_1000_Echi_drop_combo, breaks=1000, ylim=c(0,400))
abline(v=mean(ES_1000_Echi_drop_combo), lty=2)
abline(v=median(ES_1000_Echi_drop_combo), lty=1)

hist(ES_1000_Echi_man_combo, breaks=1000)
abline(v=mean(ES_1000_Echi_man_combo), lty=2)
abline(v=median(ES_1000_Echi_man_combo), lty=1)


library(plotrix)
plotCI(ES_1000_Echi_drop_combo)

## Need to get a DROPSET for this

cat EchimCaproMyo_sampling.txt | awk '{if ($7=="0") print $5 "__" $2 "__" $1 }' > ECM_unSampled.txt

cat EchimCaproMyo_sampling.txt | awk '{if ($7=="1") print $5 "__" $2 "__" $1 }' > ECM_Sampled.txt

toDrop<-read.table("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/9_PASTIS-ing/ECM_unSampled.txt", quote="")

toKeep<-read.table("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/9_PASTIS-ing/ECM_Sampled.txt", quote="")

drop<-c("Boromys_offella__ECHIMYIDAE__RODENTIA", "Boromys_torrei__ECHIMYIDAE__RODENTIA", "Brotomys_voratus__ECHIMYIDAE__RODENTIA", "Clyomys_bishopi__ECHIMYIDAE__RODENTIA", "Diplomys_caniceps__ECHIMYIDAE__RODENTIA", "Echimys_saturnus__ECHIMYIDAE__RODENTIA", "Echimys_vieirai__ECHIMYIDAE__RODENTIA", "Heteropsomys_insulans__ECHIMYIDAE__RODENTIA", "Makalata_obscura__ECHIMYIDAE__RODENTIA", "Mesomys_leniceps__ECHIMYIDAE__RODENTIA", "Olallamys_edax__ECHIMYIDAE__RODENTIA", "Phyllomys_kerri__ECHIMYIDAE__RODENTIA", "Phyllomys_medius__ECHIMYIDAE__RODENTIA", "Phyllomys_thomasi__ECHIMYIDAE__RODENTIA", "Phyllomys_unicolor__ECHIMYIDAE__RODENTIA", "Proechimys_canicollis__ECHIMYIDAE__RODENTIA", "Proechimys_chrysaeolus__ECHIMYIDAE__RODENTIA", "Proechimys_decumanus__ECHIMYIDAE__RODENTIA", "Proechimys_echinothrix__ECHIMYIDAE__RODENTIA", "Proechimys_gardneri__ECHIMYIDAE__RODENTIA", "Proechimys_goeldii__ECHIMYIDAE__RODENTIA", "Proechimys_guairae__ECHIMYIDAE__RODENTIA", "Proechimys_kulinae__ECHIMYIDAE__RODENTIA", "Proechimys_magdalenae__ECHIMYIDAE__RODENTIA", "Proechimys_mincae__ECHIMYIDAE__RODENTIA", "Proechimys_oconnelli__ECHIMYIDAE__RODENTIA", "Proechimys_pattoni__ECHIMYIDAE__RODENTIA", "Proechimys_poliopus__ECHIMYIDAE__RODENTIA", "Proechimys_semispinosus__ECHIMYIDAE__RODENTIA", "Proechimys_trinitatus__ECHIMYIDAE__RODENTIA", "Proechimys_urichi__ECHIMYIDAE__RODENTIA", "Trinomys_mirapitanga__ECHIMYIDAE__RODENTIA", "Trinomys_myosuros__ECHIMYIDAE__RODENTIA", "Geocapromys_columbianus__CAPROMYIDAE__RODENTIA", "Geocapromys_thoracatus__CAPROMYIDAE__RODENTIA", "Hexolobodon_phenax__CAPROMYIDAE__RODENTIA", "Isolobodon_montanus__CAPROMYIDAE__RODENTIA", "Isolobodon_portoricensis__CAPROMYIDAE__RODENTIA", "Mesocapromys_nanus__CAPROMYIDAE__RODENTIA", "Mesocapromys_sanfelipensis__CAPROMYIDAE__RODENTIA", "Mysateles_garridoi__CAPROMYIDAE__RODENTIA", "Mysateles_gundlachi__CAPROMYIDAE__RODENTIA", "Mysateles_meridionalis__CAPROMYIDAE__RODENTIA", "Plagiodontia_ipnaeum__CAPROMYIDAE__RODENTIA")


####
# Issue w getting branch lengths in the Bayes sample of trees at the end...

cd /file/path/to

#gets the taxon block from one of the files (in this case the first one) and makes a new tree file called "combined.t"

grep -v -e "tree gen" -e "end;" infile.nex.run1.t > combined_manual-180M.t

#searches for the trees using the term "tree rep" in both tree files and dumps them into alternate lines in the combined file

paste -d"\n" <(grep "tree gen" infile.nex.run1_500burnin.t) <(grep "tree gen" infile.nex.run2_500burnin.t) <(grep "tree gen" infile.nex.run3_500burnin.t) <(grep "tree gen" infile.nex.run4_500burnin.t) >> combined_manual-180M.t

#adds the closing nexus code
echo "end;" >> combined_manual-180M.t






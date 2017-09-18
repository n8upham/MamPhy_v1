library(ape)

####

####
# PRUNE the MAMPHY from 5910 taxa to 4098 mammals (5911 to 4099)
###

setwd("~/MamPhy_analyses")

# scan in trees
trees=scan("MamPhy_fullPosterior_YuleCR_all10k.trees", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

# List of taxa to DROP
ss<-as.vector(read.table("MamPhy_FIN4_1812sp_missing_LIST.txt", header=FALSE))

# global -- with imputed-- 5910 species
#library(foreach);library(doSNOW)
#cl = makeCluster(20, type = 'SOCK')
#registerDoSNOW(cl)

# For each tree, read in as PHYLO and then prune and write to file (append)
ntrees = 10000

#treesDrop = foreach(i=1:ntrees, .combine = cbind) %dopar% {
#	phy <- read.tree(text=trees[[i]])
#	phy_drop <- drop.tip(phy, phy$tip.label[match(ss[,1], phy$tip.label)])
#	return(phy_drop)
#	write.tree(phy_drop, "MamPhy_fullPosterior_YuleCR_all10k_pruned_4098spp.trees", append = TRUE)
#}


for (i in 1:ntrees) {
	phy <- read.tree(text=trees[[i]])
	phy_drop <- drop.tip(phy, phy$tip.label[match(ss[,1], phy$tip.label)])
	write.tree(phy_drop, "MamPhy_fullPosterior_YuleCR_all10k_pruned_4098spp.trees", append = TRUE)
}



##################
###########
#######    functions

#the CAIC function from Iain Martyn August 2011
#this produces a clade.matrix from a newick file very quickly, bypassing ape.  This speeds up ED production
#Iain says he now has ED even 10x faster than the previous 10x faster...
#Updated by Gavin Thomas April 11 2013 (readCAICnew)
#Now assigns names to rows and columns of clade.matrix and gives names to each element in the list

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

#####

trees=scan("MamPhy_fullPosterior_YuleCR_all10k_pruned_4098spp.trees", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

trees=strsplit(trees,"[;]")  #trees now indexed #tree1=trees[[1]]  #first tree

# match all names with the followng rownames
cm = readCAIC(trees[[1]])

#run parallell
ntrees = length(trees)

g_es_dr_DROP= foreach(i=1:ntrees, .combine = cbind) %dopar% {   
   res1 = 1/ES_v2(readCAIC(trees[[i]]))
   res2 = as.matrix(res1[match(cm$tip.label, rownames(res1)),]) 
   rm(res1)
   return(res2)
}

write.table(g_es_dr_DROP, "DR-matrix_MamPhy_YuleCR_all10ktrees_DROP_to4098spp.txt")

pdf(file = "DRplot_4098spp.pdf", height=8, width=8, title="DR - MamPhy 4098 spp, Yule CR, 10k trees")
hist(unlist(g_es_dr_DROP), breaks=1000)
dev.off()

all10kDR<-read.table("DR-matrix_MamPhy_YuleCR_all10ktrees.txt", header=TRUE)

pdf(file = "DRplot_5910spp.pdf", height=8, width=8, title="DR - MamPhy 5910 spp, Yule CR, 10k trees")
hist(unlist(all10kDR), breaks=1000)
dev.off()



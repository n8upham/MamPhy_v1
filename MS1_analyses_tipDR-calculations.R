

# Calculate tip DR across the 10k trees
### 

library(ape)

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

setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors")


# IF FBD samples (need to prune extinct taxa first)::dev.size
###
# load in as NEXUS
#MamPhy_orig<-read.nexus(file="MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_all10k_v2_nexus.trees")
MamPhy_orig<-read.nexus(file="MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_all10k_v2_nexus.trees")

# prune the unneeded tips...
MamPhy<-vector("list",length(MamPhy_orig))
for(i in 1: length(MamPhy_orig)){
MamPhy[[i]]<-drop.tip(MamPhy_orig[[i]], c(as.vector(MamPhy_orig[[1]]$tip.label[grep("X_",MamPhy_orig[[1]]$tip.label)]),"_Anolis_carolinensis") )

# save it as a NEWICK
#write.tree(MamPhy[[i]], file="MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_all10k_v2_prune5911.trees", append=TRUE)
write.tree(MamPhy[[i]], file="MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_all10k_v2_prune4098.trees", append=TRUE)
}


#==================
# DR
#==================

#######
library(foreach);library(doSNOW)
cl = makeCluster(5, type = 'SOCK', outfile="")
registerDoSNOW(cl)


# scan in NEWICK FORMATTED trees
#####

#trees1=scan("MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_all10k_v2.trees", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 
#trees1=scan("MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_all10k_v2.trees", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 
#trees1=scan("MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_all10k_v2_prune5911.trees", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 
trees1=scan("MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_all10k_v2_prune4098.trees", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

trees=strsplit(trees1,"[;]")  #trees now indexed #tree1=trees[[1]]  #first tree


# match all names with the followng rownames
cm = readCAIC(trees[[1]])

#run parallell
ntrees = length(trees)

g_es_dr= foreach(i=1:ntrees, .combine = cbind) %dopar% {   
   res1 = 1/ES_v2(readCAIC(trees[[i]]))
   res2 = as.matrix(res1[match(cm$tip.label, rownames(res1)),]) 
   rm(res1)
   return(res2)
}


#write.table(g_es_dr, "DR-matrix_MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_all10k_v2.txt")
#write.table(g_es_dr, "DR-matrix_MamPhy_BDvr_DNAonly_4098sp_topoFree_NDexp_all10k_v2.txt")
#write.table(g_es_dr, "DR-matrix_MamPhy_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_all10k_v2_prune5911.txt")
write.table(g_es_dr, "DR-matrix_MamPhy_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_all10k_v2_prune4098.txt")




####
# Now loading the matrix to calculate HARM MEANS

library(data.table)

#all10kDR<-fread("DR-matrix_MamPhy_BDvr_DNAonly_4098sp_topoFree_NDexp_all10k_v2.txt", header="auto")
#all10kDR<-fread("DR-matrix_MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_all10k_v2.txt", header="auto")
#all10kDR<-fread("DR-matrix_MamPhy_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_all10k_v2_prune5911.txt", header="auto")
all10kDR<-fread("DR-matrix_MamPhy_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_all10k_v2_prune4098.txt", header="auto")

nRows<-4099
#nRows<-4100
#nRows<-5912
#nRows<-5911

# Will want the mean, median, mode values across the ROWS for each tip in the tree >> 
harmMeans<-data.frame(matrix(NA, nrow = nRows, ncol = 1), row.names=all10kDR$V1)
medians<-data.frame(matrix(NA, nrow = nRows, ncol = 1), row.names=all10kDR$V1)
means<-data.frame(matrix(NA, nrow = nRows, ncol = 1), row.names=all10kDR$V1)
Count.of.NA<-data.frame(matrix(NA, nrow = nRows, ncol = 1), row.names=all10kDR$V1)
range<-data.frame(matrix(NA, nrow = nRows, ncol = 1), row.names=all10kDR$V1)
variance<-data.frame(matrix(NA, nrow = nRows, ncol = 1), row.names=all10kDR$V1)
stdev<-data.frame(matrix(NA, nrow = nRows, ncol = 1), row.names=all10kDR$V1)
cv<-data.frame(matrix(NA, nrow = nRows, ncol = 1), row.names=all10kDR$V1)
sterror<-data.frame(matrix(NA, nrow = nRows, ncol = 1), row.names=all10kDR$V1)
	low95<-data.frame(matrix(NA, nrow = nRows, ncol = 1), row.names=all10kDR$V1)
	high95<-data.frame(matrix(NA, nrow = nRows, ncol = 1), row.names=all10kDR$V1)

for (i in 1:length(rownames(all10kDR))){
	x <- as.numeric(unlist(all10kDR[i,])[2:10001])
	harmMeans[i,] <- 1/(mean(1/x))
	medians[i,] <- median(x)
	means[i,] <- mean(x)
	Count.of.NA[i,] <- length(x[is.na(x)])
	range[i,] <- max(x) - min(x)
	variance[i,] <- var(x)
	stdev[i,] <- sd(x)
	cv[i,] <- sd(x)/mean(x)
	sterror[i,] <- sd(x)/sqrt(length(x[!is.na(x)]))
		low95[i,] <- as.numeric(quantile(x, 0.025))
		high95[i,] <- as.numeric(quantile(x, 0.975))
}


all10kDR_summary <- cbind(harmMeans,medians,means,Count.of.NA,range,variance,stdev,cv,sterror,low95, high95, deparse.level=1)

colnames(all10kDR_summary) <- c("harmMeans","medians","means","Count.of.NA","range","variance","stdev","cv","sterror","low95", "high95")
rownames(all10kDR_summary) <- all10kDR$V1

#write.table(all10kDR_summary, "DR-SUMMARY_MamPhy_BDvr_DNAonly_4098sp_topoFree_NDexp_all10k_v2_expanded.txt", col.names=TRUE, row.names=TRUE)
#write.table(all10kDR_summary, "DR-SUMMARY_MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_all10k_v2_expanded.txt", col.names=TRUE, row.names=TRUE)
#write.table(all10kDR_summary, "DR-SUMMARY_MamPhy_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_all10k_v2_expanded.txt", col.names=TRUE, row.names=TRUE)
write.table(all10kDR_summary, "DR-SUMMARY_MamPhy_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_all10k_v2_expanded.txt", col.names=TRUE, row.names=TRUE)







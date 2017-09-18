
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



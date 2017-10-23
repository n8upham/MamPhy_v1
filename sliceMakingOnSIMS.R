##
# Make the time slices at 5 Ma for each of the three SIMS sets.
###
library(moments); library(nlme); library(ape); library(picante); library(phytools); library(geiger)
library(foreach);library(doSNOW)

cl = makeCluster(100, type = 'MPI', outfile="")
registerDoSNOW(cl)

ntrees<-100
bbone<-"NDexp"
#for(i in 1:ntrees){
foreach(i=1:ntrees, .packages=c('geiger','moments', 'nlme', 'ape', 'picante', 'phytools'), .combine=cbind, .verbose=TRUE) %dopar% {

sims<-c("mamPhyE","lowE_0p2","highE_0p8")

lengthsALL_SIMS<-vector("list",length(sims))
sumDat_SIMS<-vector("list",length(sims))
for(z in 1:length(sims)){

## read in 1 of 100 sumulated full trees
simPhy<-read.tree(paste("MamPhy_SIMS_",sims[z],"_",bbone,"_tree",i,".tre",sep=""))

#=======================================
# Make time slices every 5 Ma to create clades
#=======================================
sliceEvery = 5 # million years
upTo = 70 # million years
numSlices = upTo/sliceEvery

root=max(node.age(simPhy)$ages)

allCladeSets<-vector("list",length=numSlices)
slicePhys<-vector("list",length=numSlices)
allCladeSetNames<-vector("list",length=numSlices)
for (j in 1:numSlices){
	allCladeSetNames[[j]]<-paste(sliceEvery*j,"Ma",sep="")
	# make tipward slice clades
	allCladeSets[[j]]<-treeSlice(simPhy, slice=root-(sliceEvery*j), trivial=TRUE) # now keeping trivial (single tips)
		# write clades to per-slice files
		write.tree(allCladeSets[[j]],file=paste("MamPhy_SIMS_",sims[z],"_",bbone,"_sample100_",i,"_timeslice_cladeTrees_",(sliceEvery*j),"Ma.trees",sep=""))
	# make rootward slice *backbones*
	slicePhys[[j]]<-treeSlice(simPhy, slice=root-(sliceEvery*j), trivial=TRUE, orientation="rootwards") # toward root, slicePhys
}

# re-label the slicePhys with standard names
# Note: need to match tips with nodes to re-label these properly! 
for (j in 1:numSlices){
	cladeSet<-allCladeSets[[j]]
	newTipNames<-paste(i,"_",j,"_",c(1:length(cladeSet)),sep="")
	
	for (k in 1:length(cladeSet)){
		cladeSp<-cladeSet[[k]]$tip.label
		if(length(cladeSp)==1){ 
			slicePhys[[j]]$tip.label[which(slicePhys[[j]]$tip.label==cladeSp)]<-newTipNames[k]
		} else {
			node <- getMRCA(simPhy, cladeSp) # find the MRCA node of those species
			slicePhys[[j]]$tip.label[which(slicePhys[[j]]$tip.label==node)]<-newTipNames[k]
		}
	}
}

# write slice phys
for(j in 1:length(slicePhys)){
	write.tree(slicePhys[[j]], file=paste("MamPhy_SIMS_",sims[z],"_",bbone,"_tree_",i,"_CORRECT_slicePhy-5to70Ma_by5.trees",sep=""), append=TRUE)
}

# record clade sizes per slice
allLengths_i<-vector("list",length(allCladeSets))
for (j in 1:length(allCladeSets)){
	lengths<-vector()
	for(k in 1:length(allCladeSets[[j]])){
		lengths[k]<-length(allCladeSets[[j]][[k]]$tip.label)
	}
allLengths_i[[j]]<-lengths
}
names(allLengths_i)<-allCladeSetNames 
save(allLengths_i, file=paste("MamPhy_SIMS_",sims[z],"_",bbone,"_sample100_",i,"_timeslice_cladeRichnesses_wSingletons.Rda",sep=""))


} # end 3 sim loop

} # end 100 tree loop

q()

n


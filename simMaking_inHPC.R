
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy v1 -- Upham et al. 2017
###
# Figure 3 - time slices to explain clade richness
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# For SIMULATIONS, do time slice analyses
######
# - Generate simulations of the full MamPhy for testing null expectations.


# ~~~~~~~~
# Simulate as LamMu CONSTANT / N as 5911 / root age SET as T == REJECTION or DIRECT sampling
# This was the FIRST ROUND, Dec 2016.
# ====================================
# initialize
library(moments); library(nlme); library(ape); library(picante); library(phytools); library(geiger)
library(foreach);library(doSNOW)

cl = makeCluster(100, type = 'MPI', outfile="")
registerDoSNOW(cl)

ntrees=100
foreach(i=1:ntrees, .packages=c('moments', 'nlme', 'ape', 'picante', 'phytools','geiger'), .combine=cbind, .verbose=TRUE) %dopar% {

# which backbone?
bbone<- "NDexp" #"FBD" # 

mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_nexus.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 
root<-max(node.age(mamPhy)$ages)

bd<-birthdeath(mamPhy)

extFrac<-bd$para[[1]] # d/b
netDiv<-bd$para[[2]] # b-d
lam<-netDiv/(1-extFrac)
mu<-lam-netDiv
rateTable<-cbind.data.frame(extFrac,netDiv,lam,mu)
write.table(rateTable,paste("MamPhy_SIMS_mamPhyE_",bbone,"_tree",i,"newRoundOfSims_RATES.txt",sep=""))

# sim for mamPhy rates
start.time <- Sys.time()
	simMam<-pbtree(b=lam,d=mu,n=5911,t=root, scale=NULL,nsim=1,type="continuous", extant.only=FALSE, max.count=1e5, method="direct", ape=FALSE)
#	simMam<-pbtree(b=lam,d=mu,n=5911,t=root, scale=NULL,nsim=1,type="continuous", extant.only=TRUE, max.count=1e5, method="rejection", ape=FALSE)
end.time <- Sys.time()
time.taken <- end.time - start.time

write.table(time.taken,paste("MamPhy_SIMS_mamPhyE_",bbone,"_tree",i,"_timeTakenFor1Sim_Direct.txt",sep=""))
write.tree(simMam, paste("MamPhy_SIMS_mamPhyE_",bbone,"_tree",i,"_RC_N5911_Troot_Direct.tre",sep=""))

#write.table(time.taken,paste("MamPhy_SIMS_mamPhyE_",bbone,"_tree",i,"_timeTakenFor1Sim_Rejection.txt",sep=""))
#write.tree(simMam, paste("MamPhy_SIMS_mamPhyE_",bbone,"_tree",i,"_RC_N5911_Troot_Rejection.tre",sep=""))

} # end 100 tree loop


#start.time <- Sys.time()
#	simMam<-pbtree(b=lam,d=mu,n=5911,t=root, scale=NULL,nsim=1,type="continuous", extant.only=TRUE, max.count=1e5, method="rejection", ape=FALSE)
#end.time <- Sys.time()
#time.taken <- end.time - start.time
#time.taken

stopCluster(cl)

q()

n



## ~~~~~~~~
## Simulate as LamMu CONSTANT / N as 5911 / root age SCALED after simulation
## This was the FIRST ROUND, Dec 2016.
# ====================================
#
## which backbone?
#bbone<- "NDexp" #"FBD" # 
#
#mamPhy<-ladderize(drop.tip(read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",i,"_nexus.tre",sep="")), "_Anolis_carolinensis")) ## Use the MCC target node heights one... 
#root<-max(node.age(mamPhy)$ages)
#
#bd<-birthdeath(mamPhy)
#
#extFrac<-bd$para[[1]] # d/b
#netDiv<-bd$para[[2]] # b-d
#lam<-netDiv/(1-extFrac)
#mu<-lam-netDiv
#rateTable<-cbind.data.frame(extFrac,netDiv,lam,mu)
#write.table(rateTable,paste("MamPhy_SIMS_mamPhyE_",bbone,"_tree",i,"_RATES.txt",sep=""))
## sim for mamPhy rates
#	for (j in 1:10){
#	simMam<-pbtree(b=lam,d=mu,n=5911,t=NULL, scale=root,nsim=1,type="continuous", extant.only=TRUE)
#	if (class(simMam) == "NULL"){
#		cat("another sim...")
#		simMam<-pbtree(b=lam,d=mu,n=5911,t=NULL, scale=root,nsim=1,type="continuous", extant.only=TRUE)
#		} else break
#		cat("got one, simulating...")
#	}
#write.tree(simMam, paste("MamPhy_SIMS_mamPhyE_",bbone,"_tree",redos[i],".tre",sep=""))
#
## sim for e=0.8, same lam
#muHiE<-0.8*lam
#	for (j in 1:20){
#	simHiE<-pbtree(b=lam,d=muHiE,n=5911,t=NULL, scale=root,nsim=1,type="continuous", extant.only=TRUE)
#	if (class(simHiE) == "NULL"){
#		cat("another sim...")
#		simHiE<-pbtree(b=lam,d=muHiE,n=5911,t=NULL, scale=root,nsim=1,type="continuous", extant.only=TRUE)
#		} else break
#		cat("got one, simulating...")
#	}
#write.tree(simHiE, paste("MamPhy_SIMS_highE_0p8_",bbone,"_tree",redos[i],".tre",sep=""))
#
## sim for e=0.2, same net div 
#muLoE<-0.2*lam
#	simLoE<-pbtree(b=lam,d=muLoE,n=5911,t=NULL, scale=root,nsim=1,type="continuous", extant.only=TRUE)
#
#write.tree(simLoE, paste("MamPhy_SIMS_lowE_0p2_",bbone,"_tree",redos[i],".tre",sep=""))
#
#
#}
#
#}
#




library(ape)
library(TreeSim)

setwd("/mnt/data2/scratch/Nate_Backup")

mamExp <- read.nexus("MamPhy_fullPosterior_BDvarRates_17Exp_all10k_nexus.trees")

mamUni <- read.nexus("MamPhy_fullPosterior_BDvarRates_17Uni_all10k_nexus.trees")

jpeg(file="LTT_MamPhy_BDvarRates_Exp-v-Uni_all10k_noOut.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamExp[[1]],"_Anolis_carolinensis"), log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:length(mamExp)) {
	ltt.lines(drop.tip(mamExp[[i]],"_Anolis_carolinensis"), col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}
for (i in 1:length(mamUni)) {
	ltt.lines(drop.tip(mamUni[[i]],"_Anolis_carolinensis"), col=hsv(0.9708995,0.2470588,1,alpha=0.5), lty = 2, lwd=1)
}

title(main="MamPhy, BD varRates, 5910 spp, all10k - Exp (blue) vs Uni (pink)")

dev.off()

#######
# With SIMULATIONS

mamSims<-sim.bd.taxa.age(n=5910, numbsim=10, lambda=0.09234, mu=0.000000000001378, frac = 1.0, age=187.3, mrca = TRUE)
# lambda=0.09234, mu=0.000000000001378 >> values from EXP backbone post-burnin... but that doesn't make sense either bc not species level.

lamMu<-data.frame(matrix(NA, nrow = length(mamExp), ncol = 2))
colnames(lamMu)<-c("lam","mu")
for(i in 1:length(mamExp)){
	RES<-bd.shifts.optim(getx(mamExp[[i]]),sampling=1, grid=180, start=180, end=0, maxitk = 5, yule = FALSE, ME = FALSE, all = FALSE, posdiv = FALSE, miniall = c(0), survival = 1,groups=0)
	lamMu[i,1]<-RES[[1]][[1]][[1]]$par[1]
	lamMu[i,2]<-RES[[1]][[1]][[1]]$par[2]
}

# >>> This will provide speciation and extinction rates across 100 full trees
meanLam<-mean(lamMu[,1]) # [1] 0.6788691
meanMu<-mean(lamMu[,2]) # [1] 0.08429332

mamSims<-sim.bd.taxa.age(n=5910, numbsim=100, lambda=meanLam, mu=meanMu, frac = 1.0, age=187.3, mrca = TRUE)

# NOW PLOT:

jpeg(file="LTT_MamPhy_BDvarRates_Exp-v-Uni_all10k_noOut_vSIMs.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamExp[[1]],"_Anolis_carolinensis"), log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:length(mamExp)) {
	ltt.lines(drop.tip(mamExp[[i]],"_Anolis_carolinensis"), col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}
for (i in 1:length(mamUni)) {
	ltt.lines(drop.tip(mamUni[[i]],"_Anolis_carolinensis"), col=hsv(0.9708995,0.2470588,1,alpha=0.5), lty = 2, lwd=1)
}
for (i in 1:length(mamSims)) {
	ltt.lines(mamSims[[i]], col=hsv(0.4,0.8,1,alpha=0.5), lty = 2, lwd=1)
}

title(main="MamPhy, BD varRates, 5910 spp v SIMs - Exp (blue) vs Uni (pink)")

dev.off()
# >>> Those SIMS end up being mad weird ! With constant rates at that level it needs to be very rapid speciation...



##
jpeg(file="LTT_MamPhy_BDvarRates_Exp-v-Uni_all10k_noOut_wLINE.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamExp[[1]],"_Anolis_carolinensis"), log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:length(mamExp)) {
	ltt.lines(drop.tip(mamExp[[i]],"_Anolis_carolinensis"), col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}
for (i in 1:length(mamUni)) {
	ltt.lines(drop.tip(mamUni[[i]],"_Anolis_carolinensis"), col=hsv(0.9708995,0.2470588,1,alpha=0.5), lty = 2, lwd=1)
}
segments(-max(nodeHeights(drop.tip(mamUni[[1]],"_Anolis_carolinensis"))),log(2),0,5910,lty=2, col="red")

title(main="MamPhy, BD varRates, 5910 spp wLINE - Exp (blue) vs Uni (pink)")

dev.off()








mamExp2 <- read.nexus("MamPhy_fullPosterior_BDvarRates_17Exp_all10k_pruned_4098spp_nexus.trees")

mamUni2 <- read.nexus("MamPhy_fullPosterior_BDvarRates_17Uni_all10k_pruned_4098spp_nexus.trees")

jpeg(file="LTT_MamPhy_BDvarRates_Exp-v-Uni_all10k_pruned4098_noOut.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamExp2[[1]],"_Anolis_carolinensis"), log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:length(mamExp2)) {
	ltt.lines(drop.tip(mamExp2[[i]],"_Anolis_carolinensis"), col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}
for (i in 1:length(mamUni2)) {
	ltt.lines(drop.tip(mamUni2[[i]],"_Anolis_carolinensis"), col=hsv(0.9708995,0.2470588,1,alpha=0.5), lty = 2, lwd=1)
}

title(main="MamPhy, BD varRates, pruned4098 spp, all10k - Exp (blue) vs Uni (pink)")

dev.off()


jpeg(file="LTT_MamPhy_BDvarRates_Exp-v-Uni_all10k_5910-vs-pruned4098_noOut.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamExp[[1]],"_Anolis_carolinensis"), log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:length(mamExp)) {
	ltt.lines(drop.tip(mamExp[[i]],"_Anolis_carolinensis"), col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}
for (i in 1:length(mamUni)) {
	ltt.lines(drop.tip(mamUni[[i]],"_Anolis_carolinensis"), col=hsv(0.9708995,0.2470588,1,alpha=0.5), lty = 2, lwd=1)
}
for (i in 2:length(mamExp2)) {
	ltt.lines(drop.tip(mamExp2[[i]],"_Anolis_carolinensis"), col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}
for (i in 1:length(mamUni2)) {
	ltt.lines(drop.tip(mamUni2[[i]],"_Anolis_carolinensis"), col=hsv(0.9708995,0.2470588,1,alpha=0.5), lty = 2, lwd=1)
}

title(main="MamPhy, BD varRates, 5910 vs 4098 spp, all10k - Exp (blue) vs Uni (pink)")

dev.off()




q()

n
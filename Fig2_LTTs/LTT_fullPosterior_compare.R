library(ape)

setwd("/mnt/data2/scratch/Nate_Backup")

mamExp <- read.nexus("MamPhy_fullPosterior_BDvarRates_17Exp_all10k_nexus.trees")

mamUni <- read.nexus("MamPhy_fullPosterior_BDvarRates_17Uni_all10k_nexus.trees")

jpeg(file="LTT_MamPhy_BDvarRates_Exp-v-Uni_all10k.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(mamExp[[1]], log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:length(mamExp)) {
	ltt.lines(mamExp[[i]], col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}
for (i in 1:length(mamUni)) {
	ltt.lines(mamUni[[i]], col=hsv(0.9708995,0.2470588,1,alpha=0.5), lty = 2, lwd=1)
}

title(main="MamPhy, BD varRates, 5911 spp, all10k - Exp (blue) vs Uni (pink)")

dev.off()


mamExp2 <- read.nexus("MamPhy_fullPosterior_BDvarRates_17Exp_all10k_pruned_4098spp_nexus.trees")

mamUni2 <- read.nexus("MamPhy_fullPosterior_BDvarRates_17Uni_all10k_pruned_4098spp_nexus.trees")

jpeg(file="LTT_MamPhy_BDvarRates_Exp-v-Uni_all10k_pruned4098.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(mamExp2[[1]], log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:length(mamExp2)) {
	ltt.lines(mamExp2[[i]], col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}
for (i in 1:length(mamUni2)) {
	ltt.lines(mamUni2[[i]], col=hsv(0.9708995,0.2470588,1,alpha=0.5), lty = 2, lwd=1)
}

title(main="MamPhy, BD varRates, pruned4099 spp, all10k - Exp (blue) vs Uni (pink)")

dev.off()


jpeg(file="LTT_MamPhy_BDvarRates_Exp-v-Uni_all10k_5911-vs-pruned4098.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(mamExp2[[1]], log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:length(mamExp2)) {
	ltt.lines(mamExp2[[i]], col=hsv(0.65,1,1,alpha=0.25), lty = 2, lwd=1)
}
for (i in 1:length(mamUni2)) {
	ltt.lines(mamUni2[[i]], col=hsv(0.9708995,0.2470588,1,alpha=0.25), lty = 2, lwd=1)
}
for (i in 2:length(mamExp)) {
	ltt.lines(mamExp[[i]], col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}
for (i in 1:length(mamUni)) {
	ltt.lines(mamUni[[i]], col=hsv(0.9708995,0.2470588,1,alpha=0.5), lty = 2, lwd=1)
}

title(main="MamPhy, BD varRates, 5911 vs 4099 spp, all10k - Exp (blue) vs Uni (pink)")

dev.off()




q()

n
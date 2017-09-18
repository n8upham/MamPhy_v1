########
# BAMMtools analyses...
###
library(BAMMtools)
library(coda)

folders<-c(1:10)

#folders<-1
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/grace_2.5_1week_1")

library(foreach);library(doSNOW)
cl = makeCluster(10, type = 'SOCK', outfile="")
registerDoSNOW(cl)

ntrees = length(folders)

runs <- foreach(i=1:ntrees, .packages=c('BAMMtools', 'coda'), .combine=cbind, .verbose=TRUE) %dopar% {

#for (i in 1:length(folders)){

#setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/grace_2.5_1week_",folders[i],sep=""))
setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/grace_2.5_1week_NDexp_",folders[i],sep=""))

#tree <- read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample10_",folders[i],".tre",sep=""))
#edata <- getEventData(tree, eventdata = "mamPhy_FBD_event_data.txt", burnin=0.33, nsamples=1000)
tree <- ladderize(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample10_",folders[i],".tre",sep="")))
edata <- getEventData(tree, eventdata = "mamPhy_NDexp_event_data.txt", burnin=0.33, nsamples=1000)

# convergence

#mcmcout <- read.csv("mamPhy_FBD_mcmc_out.txt", header=T)
#pdf(file="mamPhy_FBD_mcmc_out.pdf")
mcmcout <- read.csv("mamPhy_NDexp_mcmc_out.txt", header=T)
pdf(file="mamPhy_NDexp_mcmc_out.pdf")
plot(mcmcout$logLik ~ mcmcout$generation)
dev.off()

burnstart <- floor(0.33 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

#write.table(paste("FBD_sample10_",folders[i],sep=""), file="results.txt", append=TRUE)
write.table(paste("NDexp_sample10_",folders[i],sep=""), file="results.txt", append=TRUE)

write.table(effectiveSize(postburn$N_shifts), file="results.txt", append=TRUE)

write.table(effectiveSize(postburn$logLik), file="results.txt", append=TRUE)

post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)

shift_probs <- summary(edata)

write.table(shift_probs, file="results.txt", append=TRUE)

bfmat <- computeBayesFactors(postburn, expectedNumberOfShifts=1, burnin=0.0) #burnin already done

write.table(bfmat, file="results.txt", append=TRUE)

pdf(file="prior-posterior_compare.pdf")
plotPrior(postburn, expectedNumberOfShifts = 1, burnin = 0.0)
dev.off()

# phylorate

pdf(file="mean-phylorate.pdf", width=8.5, height=200)
plot.bammdata(edata, lwd=2, legend=TRUE, labels=TRUE,cex=0.2)
dev.off()

pdf(file="mean-phylorate_polar.pdf", width=35, height=35)
plot.bammdata(edata, lwd=2, legend=TRUE, method="polar", labels=TRUE,cex=0.11)
dev.off()


pdf(file="mean-phylorate_25th_wShifts.pdf")
index <- 25
e2 <- subsetEventData(edata, index = index)
plot.bammdata(e2, lwd=2, legend=TRUE)
addBAMMshifts(e2, cex=2)
dev.off()

# Credible shifts
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)

write.table(summary(css), file="results.txt", append=TRUE)

pdf(file="credibleShiftSet.pdf")
plot.credibleshiftset(css)
dev.off()


best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1)
pdf(file="bestShiftSet.pdf")
plot.bammdata(best, lwd = 2)
addBAMMshifts(best, cex=2.5)
dev.off()

# maximum shift sets...

msc.set <- maximumShiftCredibility(edata, maximize='product')

msc.config <- subsetEventData(edata, index = msc.set$sampleindex)

pdf(file="MAX_ShiftSet.pdf")
plot.bammdata(msc.config, lwd=2)
addBAMMshifts(msc.config, cex = 2)
dev.off()

pdf(file="MAX_ShiftSet_largePolar.pdf", width=35, height=35)
plot.bammdata(msc.config, lwd=2, legend=TRUE, method="polar", labels=TRUE,cex=0.11)
addBAMMshifts(msc.config, cex = 2, method="polar")
dev.off()

write.table(msc.config$eventData, file="results.txt", append=TRUE)

dsc <- distinctShiftConfigurations(edata, expectedNumberOfShifts=1, threshold=5)

# Here is one random sample with the BEST shift configuration
pdf(file="random_bestShift1.pdf")
plot.bammshifts(dsc, edata, rank=1, legend=F)
dev.off()

# Here is another (read the label text):
pdf(file="random_bestShift2.pdf")
plot.bammshifts(dsc, edata, rank=1, legend=F)
dev.off()

# Here is one random sample with the 2nd BEST shift configuration
pdf(file="random_2ndbestShift1.pdf")
plot.bammshifts(dsc, edata, rank=2, legend=F)
dev.off()

# Here is another (read the label text):
pdf(file="random_2ndbestShift2.pdf")
plot.bammshifts(dsc, edata, rank=2, legend=F)
dev.off()

##
# Now using APE functions...

mysample <- 25  # this is the sample we'll plot

nrow(edata$eventData[[ mysample ]])

shiftnodes <- getShiftNodesFromIndex(edata, index = mysample)

pdf(file="apePlot_sample25.pdf", width=8.5,height=150)
plot.phylo(tree, cex=0.15)
nodelabels(node = shiftnodes, pch=21, col="red", cex=1.5)
dev.off()

# marginal shifts...

marg_probs <- marginalShiftProbsTree(edata)

pdf(file="marginalProbs.pdf", width=8.5,height=150)
plot.phylo(marg_probs, cex=0.15)
dev.off()

# clade-specific rates (sp and ex...)
allrates <- getCladeRates(edata)

write.table("ALL - LAMBDA", file="results.txt", append=TRUE)

write.table(mean(allrates$lambda), file="results.txt", append=TRUE)

write.table(quantile(allrates$lambda, c(0.05, 0.95)), file="results.txt", append=TRUE)

write.table("ALL - MU", file="results.txt", append=TRUE)

write.table(mean(allrates$mu), file="results.txt", append=TRUE)

write.table(quantile(allrates$mu, c(0.05, 0.95)), file="results.txt", append=TRUE)

# DOLPHIN compare
dolphinrates <- getCladeRates(edata, node= 8185) #141)

write.table("DOLPHINS", file="results.txt", append=TRUE)

write.table(mean(dolphinrates$lambda), file="results.txt", append=TRUE)

write.table(quantile(dolphinrates$lambda, c(0.05, 0.95)), file="results.txt", append=TRUE)

nondolphinrate <- getCladeRates(edata, node = 8185, nodetype = "exclude") #141, nodetype = "exclude")

write.table("NON-DOLPHINS", file="results.txt", append=TRUE)

write.table(mean(nondolphinrate$lambda), file="results.txt", append=TRUE)

write.table(quantile(nondolphinrate$lambda, c(0.05, 0.95)), file="results.txt", append=TRUE)

# per-branch rates
# components edata$meanTipLambda and edata$meanTipMu are the relevant model-averaged mean rates of speciation and extinction at the tips of the tree.
## How does this compare to DR??

BAMM_harmMeans_Lam<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=edata$tip.label)
BAMM_arithMeans_Lam<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=edata$tip.label)
BAMM_harmMeans_Mu<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=edata$tip.label)
BAMM_arithMeans_Mu<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=edata$tip.label)
BAMM_harmMeans_NetDiv<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=edata$tip.label)
BAMM_arithMeans_NetDiv<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=edata$tip.label)

for (k in 1:length(edata$tipLambda[[1]])){

tipLam<-vector()
tipMu<-vector()
tipDiv<-vector()
for (j in 1:length(edata$tipLambda)){
	tipLam[j]<-edata$tipLambda[[j]][k]
	tipMu[j]<-edata$tipMu[[j]][k]
	tipDiv[j]<-edata$tipLambda[[j]][k]-edata$tipMu[[j]][k]

}
	BAMM_harmMeans_Lam[k,]<-1/(mean(1/tipLam))
	BAMM_arithMeans_Lam[k,]<-mean(tipLam)
	BAMM_harmMeans_Mu[k,]<-1/(mean(1/tipMu))
	BAMM_arithMeans_Mu[k,]<-mean(tipMu)
	BAMM_harmMeans_NetDiv[k,]<-1/(mean(1/tipDiv))
	BAMM_arithMeans_NetDiv[k,]<-mean(tipDiv)

}

res<-cbind(BAMM_harmMeans_Lam,BAMM_arithMeans_Lam,BAMM_harmMeans_Mu,BAMM_arithMeans_Mu,BAMM_harmMeans_NetDiv,BAMM_arithMeans_NetDiv)

colnames(res)<-c("harmMeans_Lam","arithMeans_Lam","harmMeans_Mu","arithMeans_Mu","harmMeans_NetDiv","arithMeans_NetDiv")

write.table(res, paste("meanTipRates_sample10_tree",folders[i],".txt", sep=""))

#plot(density(res$harmMeans_Lam), ylim=c(0,8))
#plot(density(res$arithMeans_Lam), ylim=c(0,8))
#plot(density(edata$meanTipLambda), ylim=c(0,8)) # It is identical to the arithmetic means, good.

#####

# RATES through time...
pdf(file="rateThroughTime_spec.pdf")
plotRateThroughTime(edata, ratetype="speciation")
dev.off()

pdf(file="rateThroughTime_ext.pdf")
plotRateThroughTime(edata, ratetype="extinction")
dev.off()

pdf(file="rateThroughTime_div.pdf")
plotRateThroughTime(edata, ratetype="netdiv")
dev.off()

# cohort matrix...
cmat <- getCohortMatrix(edata)
pdf(file="cohortMatrixPlot.pdf")
cohorts(cmat, edata)
dev.off()

# cumulative shift marg_probs
cst <- cumulativeShiftProbsTree(edata)
pdf(file="cumShiftProbs_tree.pdf", width=8.5,height=150)
plot.phylo(cst, cex=0.15)
dev.off()


cst <- cumulativeShiftProbsTree(edata)
edgecols <- rep('black', length(tree$edge.length))
is_highprobshift <- cst$edge.length >= 0.95
edgecols[ is_highprobshift ] <- "red"

pdf(file="cumShiftProbs_tree_cols.pdf", width=35,height=55)
plot.phylo(tree, edge.color = edgecols, cex=0.05, type="phylogram")
dev.off()

pdf(file="cumShiftProbs_tree_cols_fan.pdf", width=35,height=35)
plot.phylo(tree, edge.color = edgecols, cex=0.15, type="fan")
dev.off()



}

#######
# for FBD again (with ladderized)

folders<-c(1:10)

#folders<-1
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/grace_2.5_1week_1")

library(foreach);library(doSNOW)
cl2 = makeCluster(10, type = 'SOCK', outfile="")
registerDoSNOW(cl2)

ntrees = length(folders)

runs <- foreach(i=1:ntrees, .packages=c('BAMMtools', 'coda'), .combine=cbind, .verbose=TRUE) %dopar% {

#for (i in 1:length(folders)){

setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/grace_2.5_1week_",folders[i],sep=""))
#setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/grace_2.5_1week_NDexp_",folders[i],sep=""))

tree <- ladderize(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample10_",folders[i],".tre",sep="")))
edata <- getEventData(tree, eventdata = "mamPhy_FBD_event_data.txt", burnin=0.33, nsamples=1000)
#tree <- ladderize(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample10_",folders[i],".tre",sep="")))
#edata <- getEventData(tree, eventdata = "mamPhy_NDexp_event_data.txt", burnin=0.33, nsamples=1000)

# convergence

mcmcout <- read.csv("mamPhy_FBD_mcmc_out.txt", header=T)
pdf(file="mamPhy_FBD_mcmc_out.pdf")
#mcmcout <- read.csv("mamPhy_NDexp_mcmc_out.txt", header=T)
#pdf(file="mamPhy_NDexp_mcmc_out.pdf")
plot(mcmcout$logLik ~ mcmcout$generation)
dev.off()

burnstart <- floor(0.33 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

write.table(paste("FBD_sample10_",folders[i],sep=""), file="results.txt", append=TRUE)
#write.table(paste("NDexp_sample10_",folders[i],sep=""), file="results.txt", append=TRUE)

write.table(effectiveSize(postburn$N_shifts), file="results.txt", append=TRUE)

write.table(effectiveSize(postburn$logLik), file="results.txt", append=TRUE)

post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)

shift_probs <- summary(edata)

write.table(shift_probs, file="results.txt", append=TRUE)

bfmat <- computeBayesFactors(postburn, expectedNumberOfShifts=1, burnin=0.0) #burnin already done

write.table(bfmat, file="results.txt", append=TRUE)

pdf(file="prior-posterior_compare.pdf")
plotPrior(postburn, expectedNumberOfShifts = 1, burnin = 0.0)
dev.off()

# phylorate

pdf(file="mean-phylorate.pdf", width=8.5, height=200)
plot.bammdata(edata, lwd=2, legend=TRUE, labels=TRUE,cex=0.2)
dev.off()

pdf(file="mean-phylorate_polar.pdf", width=35, height=35)
plot.bammdata(edata, lwd=2, legend=TRUE, method="polar", labels=TRUE,cex=0.11)
dev.off()


pdf(file="mean-phylorate_25th_wShifts.pdf")
index <- 25
e2 <- subsetEventData(edata, index = index)
plot.bammdata(e2, lwd=2, legend=TRUE)
addBAMMshifts(e2, cex=2)
dev.off()

# Credible shifts
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)

write.table(summary(css), file="results.txt", append=TRUE)

pdf(file="credibleShiftSet.pdf")
plot.credibleshiftset(css)
dev.off()


best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1)
pdf(file="bestShiftSet.pdf")
plot.bammdata(best, lwd = 2)
addBAMMshifts(best, cex=2.5)
dev.off()

# maximum shift sets...

msc.set <- maximumShiftCredibility(edata, maximize='product')

msc.config <- subsetEventData(edata, index = msc.set$sampleindex)

pdf(file="MAX_ShiftSet.pdf")
plot.bammdata(msc.config, lwd=2)
addBAMMshifts(msc.config, cex = 2)
dev.off()

pdf(file="MAX_ShiftSet_largePolar.pdf", width=35, height=35)
plot.bammdata(msc.config, lwd=2, legend=TRUE, method="polar", labels=TRUE,cex=0.11)
addBAMMshifts(msc.config, cex = 2, method="polar")
dev.off()

write.table(msc.config$eventData, file="results.txt", append=TRUE)

dsc <- distinctShiftConfigurations(edata, expectedNumberOfShifts=1, threshold=5)

# Here is one random sample with the BEST shift configuration
pdf(file="random_bestShift1.pdf")
plot.bammshifts(dsc, edata, rank=1, legend=F)
dev.off()

# Here is another (read the label text):
pdf(file="random_bestShift2.pdf")
plot.bammshifts(dsc, edata, rank=1, legend=F)
dev.off()

# Here is one random sample with the 2nd BEST shift configuration
pdf(file="random_2ndbestShift1.pdf")
plot.bammshifts(dsc, edata, rank=2, legend=F)
dev.off()

# Here is another (read the label text):
pdf(file="random_2ndbestShift2.pdf")
plot.bammshifts(dsc, edata, rank=2, legend=F)
dev.off()

##
# Now using APE functions...

mysample <- 25  # this is the sample we'll plot

nrow(edata$eventData[[ mysample ]])

shiftnodes <- getShiftNodesFromIndex(edata, index = mysample)

pdf(file="apePlot_sample25.pdf", width=8.5,height=150)
plot.phylo(tree, cex=0.15)
nodelabels(node = shiftnodes, pch=21, col="red", cex=1.5)
dev.off()

# marginal shifts...

marg_probs <- marginalShiftProbsTree(edata)

pdf(file="marginalProbs.pdf", width=8.5,height=150)
plot.phylo(marg_probs, cex=0.15)
dev.off()

# clade-specific rates (sp and ex...)
allrates <- getCladeRates(edata)

write.table("ALL - LAMBDA", file="results.txt", append=TRUE)

write.table(mean(allrates$lambda), file="results.txt", append=TRUE)

write.table(quantile(allrates$lambda, c(0.05, 0.95)), file="results.txt", append=TRUE)

write.table("ALL - MU", file="results.txt", append=TRUE)

write.table(mean(allrates$mu), file="results.txt", append=TRUE)

write.table(quantile(allrates$mu, c(0.05, 0.95)), file="results.txt", append=TRUE)

# DOLPHIN compare
dolphinrates <- getCladeRates(edata, node= 8185) #141)

write.table("DOLPHINS", file="results.txt", append=TRUE)

write.table(mean(dolphinrates$lambda), file="results.txt", append=TRUE)

write.table(quantile(dolphinrates$lambda, c(0.05, 0.95)), file="results.txt", append=TRUE)

nondolphinrate <- getCladeRates(edata, node = 8185, nodetype = "exclude") #141, nodetype = "exclude")

write.table("NON-DOLPHINS", file="results.txt", append=TRUE)

write.table(mean(nondolphinrate$lambda), file="results.txt", append=TRUE)

write.table(quantile(nondolphinrate$lambda, c(0.05, 0.95)), file="results.txt", append=TRUE)

# per-branch rates
# components edata$meanTipLambda and edata$meanTipMu are the relevant model-averaged mean rates of speciation and extinction at the tips of the tree.
## How does this compare to DR??

BAMM_harmMeans_Lam<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=edata$tip.label)
BAMM_arithMeans_Lam<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=edata$tip.label)
BAMM_harmMeans_Mu<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=edata$tip.label)
BAMM_arithMeans_Mu<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=edata$tip.label)
BAMM_harmMeans_NetDiv<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=edata$tip.label)
BAMM_arithMeans_NetDiv<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=edata$tip.label)

for (k in 1:length(edata$tipLambda[[1]])){

tipLam<-vector()
tipMu<-vector()
tipDiv<-vector()
for (j in 1:length(edata$tipLambda)){
	tipLam[j]<-edata$tipLambda[[j]][k]
	tipMu[j]<-edata$tipMu[[j]][k]
	tipDiv[j]<-edata$tipLambda[[j]][k]-edata$tipMu[[j]][k]

}
	BAMM_harmMeans_Lam[k,]<-1/(mean(1/tipLam))
	BAMM_arithMeans_Lam[k,]<-mean(tipLam)
	BAMM_harmMeans_Mu[k,]<-1/(mean(1/tipMu))
	BAMM_arithMeans_Mu[k,]<-mean(tipMu)
	BAMM_harmMeans_NetDiv[k,]<-1/(mean(1/tipDiv))
	BAMM_arithMeans_NetDiv[k,]<-mean(tipDiv)

}

res<-cbind(BAMM_harmMeans_Lam,BAMM_arithMeans_Lam,BAMM_harmMeans_Mu,BAMM_arithMeans_Mu,BAMM_harmMeans_NetDiv,BAMM_arithMeans_NetDiv)

colnames(res)<-c("harmMeans_Lam","arithMeans_Lam","harmMeans_Mu","arithMeans_Mu","harmMeans_NetDiv","arithMeans_NetDiv")

write.table(res, paste("meanTipRates_sample10_tree",folders[i],".txt", sep=""))


#plot(density(res$harmMeans_Lam), ylim=c(0,8))
#plot(density(res$arithMeans_Lam), ylim=c(0,8))
#plot(density(edata$meanTipLambda), ylim=c(0,8)) # It is identical to the arithmetic means, good.

#####

# RATES through time...
pdf(file="rateThroughTime_spec.pdf")
plotRateThroughTime(edata, ratetype="speciation")
dev.off()

pdf(file="rateThroughTime_ext.pdf")
plotRateThroughTime(edata, ratetype="extinction")
dev.off()

pdf(file="rateThroughTime_div.pdf")
plotRateThroughTime(edata, ratetype="netdiv")
dev.off()

# cohort matrix...
cmat <- getCohortMatrix(edata)
pdf(file="cohortMatrixPlot.pdf")
cohorts(cmat, edata)
dev.off()

# cumulative shift marg_probs
cst <- cumulativeShiftProbsTree(edata)
pdf(file="cumShiftProbs_tree.pdf", width=8.5,height=150)
plot.phylo(cst, cex=0.15)
dev.off()


cst <- cumulativeShiftProbsTree(edata)
edgecols <- rep('black', length(tree$edge.length))
is_highprobshift <- cst$edge.length >= 0.95
edgecols[ is_highprobshift ] <- "red"

pdf(file="cumShiftProbs_tree_cols.pdf", width=35,height=55)
plot.phylo(tree, edge.color = edgecols, cex=0.05, type="phylogram")
dev.off()

pdf(file="cumShiftProbs_tree_cols_fan.pdf", width=35,height=35)
plot.phylo(tree, edge.color = edgecols, cex=0.15, type="fan")
dev.off()

}




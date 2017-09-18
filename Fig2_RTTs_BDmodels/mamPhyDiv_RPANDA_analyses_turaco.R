library(ape)
library(RPANDA)


######
# with MAMPHY:

setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors")

mamPhy10k<-read.tree("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_all10k.trees")

trees1=scan("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_all10k.trees", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

#get 1 random tree
phy1<-read.tree(text=sample(trees1,size=1))

phy1a<-drop.tip(phy1,"_Anolis_carolinensis")

write.tree(phy1a, "MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample1.tre")

# get 10 random trees
phy10<-list()
phy_text10<-sample(trees1,size=10)

for (i in 1:length(phy_text10)){
	phy<-drop.tip(read.tree(text=phy_text10[[i]]),"_Anolis_carolinensis")
	write.tree(phy, "MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample10.tre", append=TRUE)
}




# load 1 random tree
phy<-read.tree("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample1.tre")

# load 10 random trees
phy10<-read.tree("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample10.tre")

for (i in 1:length(phy10)){
	res1<-spectR(phy10[[i]])
	pdf(file=paste("spectR_mamPhy_FBD_10trees_",i,".pdf",sep=""))
	plot_spectR(res1)
	dev.off()
	write.table(res1[2:6],file="spectR_mamPhy_FBD_10trees_TABLE.txt", append=TRUE)
}
	


res2<-BICompare(phy10[[10]],res1$eigengap)

pdf(file="spectR_BIC_mamPhy_FBD_tree10_fan.pdf")
plot_BICompare(phy10[[10]],res2) # shows branches colored by given modalities, default to FAN
dev.off()

plot_BICompare <- function (phylo, BICompare) 
{
    if (!inherits(BICompare, "BICompare")) 
        stop("object \"BICompare\" is not of class \"BICompare\"")
    t <- max(BICompare[[2]])
    col_edge <- rainbow(t)[BICompare[[2]][phylo$edge[, 2]]]
    col_tip <- rainbow(t)[BICompare[[2]][1:length(phylo$tip.label)]]
    plot(phylo, edge.color = col_edge, tip.color = col_tip, type = "phylogram", 
        cex = 0.4)
} # changes to normal phylogram

pdf(file="spectR_BIC_mamPhy_FBD_tree10_phylogram.pdf")
plot_BICompare(phy,res2) # shows branches colored by given modalities, default to FAN
dev.off()

q()

n
###############
########
# BAMMtools analyses...
###
library(BAMMtools)
library(coda)

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/grace_2.5_1week_4")

tree <- read.tree("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample10_4.tre")
edata <- getEventData(tree, eventdata = "mamPhy_FBD_event_data.txt", burnin=0.33, nsamples=1000)

# convergence

mcmcout <- read.csv("mamPhy_FBD_mcmc_out.txt", header=T)
pdf(file="mamPhy_FBD_mcmc_out.pdf")
plot(mcmcout$logLik ~ mcmcout$generation)
dev.off()

burnstart <- floor(0.33 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
effectiveSize(postburn$N_shifts)
    var1 
630.5419 

effectiveSize(postburn$logLik)
    var1 
578.2652 

# All good now!

# shifts
post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)

shift_probs <- summary(edata)

   shifts  prob
1      20 0.001
2      21 0.001
3      23 0.001
4      24 0.004
5      25 0.007
6      26 0.014
7      27 0.021
8      28 0.023
9      29 0.026
10     30 0.045
11     31 0.044
12     32 0.059
13     33 0.044
14     34 0.052
15     35 0.057
16     36 0.063
17     37 0.067 #<< 37 shifts it seems.
18     38 0.055
19     39 0.058
20     40 0.028
21     41 0.040
22     42 0.031
23     43 0.026
24     44 0.036
25     45 0.028
26     46 0.026
27     47 0.024
28     48 0.020
29     49 0.009
30     50 0.017
31     51 0.014
32     52 0.006
33     53 0.012
34     54 0.005
35     55 0.004
36     56 0.008
37     57 0.004
38     58 0.004
39     59 0.005
40     60 0.004
41     63 0.002
42     67 0.001
43     69 0.002
44     71 0.001
45     77 0.001

bfmat <- computeBayesFactors(postburn, expectedNumberOfShifts=1, burnin=0.0) #burnin already done

## >> returns progressively larger numbers for Bayes Factors it seems...

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

   rank     probability cumulative  Core_shifts
         1      0.001      0.001         39
         2      0.001      0.002         40
         3      0.001      0.003         39
         4      0.001      0.004         34
         5      0.001      0.005         51
         6      0.001      0.006         36
         7      0.001      0.007         45
         8      0.001      0.008         44
         9      0.001      0.009         45

...omitted 941 additional distinct shift configurations
from the credible set. You can access the full set from your 
credibleshiftset object

css$number.distinct
[1] 950

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

msc.config$eventData
[[1]]
    node    time      lam1         lam2        mu1 mu2 index
1   5912   0.000 0.0256736  0.006451350 0.04109040   0     1
2   5919 168.182 0.1980130 -0.003422330 0.13049300   0     2
3   6043 174.353 0.2277210  0.006885620 0.21152700   0     3
4   6044 180.763 0.0977820  0.008906230 0.01759670   0     4
5  11461 194.504 0.1535530  0.003765290 0.05689840   0     5
6   9863 202.268 0.1229780 -0.003404850 0.01184120   0     6
7   6934 211.189 0.1936960  0.003878350 0.05519230   0     7
8   7831 215.403 0.4093570 -0.005069780 0.16206000   0     8
9   8667 218.655 0.3836340  0.008103060 0.25796900   0     9
10  6357 222.040 0.3273920  0.000541586 0.10093900   0    10
11 10685 222.894 0.4090490  0.003663970 0.22511800   0    11
12 10830 224.219 0.2614610 -0.000803713 0.01912310   0    12
13  9958 228.080 0.3493960 -0.007359080 0.00252133   0    13
14 10113 231.086 0.1979800 -0.002901760 0.02662170   0    14
15  6550 231.379 0.4600620  0.000651859 0.10264900   0    15
16  8185 233.009 0.4573310  0.000239422 0.11592200   0    16
17  6124 234.134 0.3550330 -0.002738230 0.05140980   0    17
18  8288 236.931 0.4941480 -0.003266330 0.10815100   0    18
19  4873 237.597 0.1845350  0.011280400 0.36167600   0    19
20  7663 237.780 0.4722560 -0.005489950 0.11747300   0    20
21 10292 240.549 0.6567390 -0.003602960 0.26664200   0    21
22 11043 240.772 0.5561480  0.000801693 0.03836340   0    22
23  6787 241.311 0.6115640 -0.004338310 0.06718830   0    23
24 10861 241.439 0.4918560  0.008705580 0.02551090   0    24
25 11716 241.682 0.3529800  0.004534200 0.00486131   0    25
26  9608 242.842 0.9350640 -0.005551850 0.02301210   0    26
27 11255 245.929 0.2640060 -0.002922520 0.16640400   0    27

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

mean(allrates$lambda)
[1] 0.251775
quantile(allrates$lambda, c(0.05, 0.95))
       5%       95% 
0.2350588 0.2669584 

mean(allrates$mu)
[1] 0.08200702
quantile(allrates$mu, c(0.05, 0.95))
        5%        95% 
0.05593678 0.10467251 

# could do... and compare to the DOLPHINS in your tree !!
pdf(file="mean-phylorate_wNodelabels.pdf", width=8.5, height=200)
plot.bammdata(edata, lwd=2, legend=TRUE, labels=TRUE,cex=0.2)
nodelabels(cex=0.2)
dev.off()


dolphinrates <- getCladeRates(edata, node= 8185) #141)
mean(dolphinrates$lambda)
[1] 0.3448141
quantile(dolphinrates$lambda, c(0.05, 0.95))
       5%       95% 
0.2603402 0.4355796 

nondolphinrate <- getCladeRates(edata, node = 8185, nodetype = "exclude") #141, nodetype = "exclude")
[1] 0.2508217
quantile(nondolphinrate$lambda, c(0.05, 0.95))
       5%       95% 
0.2341132 0.2661783 

# per-branch rates
# To pull out the mean rates for individual branches, you can use the function getMeanBranchLengthTree (see the ?getMeanBranchLengthTree for help on this function). The function generates a copy of your original phylogenetic tree, but where each branch length is replaced by the mean of the marginal distribution of evolutionary rates on each branch. The function can be used to extract branch-specific mean rates of speciation, extinction, net diversification, and trait evolution.
# components edata$meanTipLambda and edata$meanTipMu are the relevant model-averaged mean rates of speciation and extinction at the tips of the tree.

tree_rate<-getMeanBranchLengthTree(edata, rate="speciation") #takes a lil time

edata$meanTipLambda ## How does this compare to DR??

write.table(edata$meanTipLambda, "meanTipLambda_sample10_tree4.txt")


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



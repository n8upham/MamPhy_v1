library(ape)
library(data.table)

#setwd("~/MamPhy_analyses")
#setwd("~/MamPhy_BDvarRatesALL_fullPosterior")
setwd("/mnt/data2/scratch/Nate_Backup")

all10kDR<-fread("DR-matrix_MamPhy_BDvarRates_17Exp_all10ktrees.txt", header="auto") # 5911 obs. of  10000 variables

# Will want the mean, median, mode values across the ROWS for each tip in the tree >> 
harmMeans<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=rownames(all10kDR))
medians<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=rownames(all10kDR))
means<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=rownames(all10kDR))
Count.of.NA<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=rownames(all10kDR))
range<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=rownames(all10kDR))
variance<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=rownames(all10kDR))
stdev<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=rownames(all10kDR))
cv<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=rownames(all10kDR))
sterror<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=rownames(all10kDR))

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
}


all10kDR_summary <- cbind(harmMeans,medians,means,Count.of.NA,range,variance,stdev,cv,sterror,deparse.level=1)

colnames(all10kDR_summary) <- c("harmMeans","medians","means","Count.of.NA","range","variance","stdev","cv","sterror")
rownames(all10kDR_summary) <- all10kDR$V1

write.table(all10kDR_summary, "DR-SUMMARY_MamPhy_BDvarRates_17Exp_all10ktrees_5910sp_expanded.txt", col.names=TRUE, row.names=TRUE)


###
all10kDR<-fread("DR-matrix_MamPhy_BDvarRates_17Exp_all10ktrees_pruned4098.txt", header="auto") # 5911 obs. of  10000 variables

# Will want the mean, median, mode values across the ROWS for each tip in the tree >> 
harmMeans<-data.frame(matrix(NA, nrow = 4099, ncol = 1), row.names=rownames(all10kDR))
medians<-data.frame(matrix(NA, nrow = 4099, ncol = 1), row.names=rownames(all10kDR))
means<-data.frame(matrix(NA, nrow = 4099, ncol = 1), row.names=rownames(all10kDR))
Count.of.NA<-data.frame(matrix(NA, nrow = 4099, ncol = 1), row.names=rownames(all10kDR))
range<-data.frame(matrix(NA, nrow = 4099, ncol = 1), row.names=rownames(all10kDR))
variance<-data.frame(matrix(NA, nrow = 4099, ncol = 1), row.names=rownames(all10kDR))
stdev<-data.frame(matrix(NA, nrow = 4099, ncol = 1), row.names=rownames(all10kDR))
cv<-data.frame(matrix(NA, nrow = 4099, ncol = 1), row.names=rownames(all10kDR))
sterror<-data.frame(matrix(NA, nrow = 4099, ncol = 1), row.names=rownames(all10kDR))

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
}


all10kDR_summary <- cbind(harmMeans,medians,means,Count.of.NA,range,variance,stdev,cv,sterror,deparse.level=1)

colnames(all10kDR_summary) <- c("harmMeans","medians","means","Count.of.NA","range","variance","stdev","cv","sterror")
rownames(all10kDR_summary) <- all10kDR$V1

write.table(all10kDR_summary, "DR-SUMMARY_MamPhy_BDvarRates_17Exp_all10ktrees_pruned4098sp_expanded.txt", col.names=TRUE, row.names=TRUE)

######

all10kDR<-fread("DR-matrix_MamPhy_BDvarRates_17Uni_all10ktrees.txt", header="auto") # 5911 obs. of  10000 variables

# Will want the mean, median, mode values across the ROWS for each tip in the tree >> 
harmMeans<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=rownames(all10kDR))
medians<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=rownames(all10kDR))
means<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=rownames(all10kDR))
Count.of.NA<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=rownames(all10kDR))
range<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=rownames(all10kDR))
variance<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=rownames(all10kDR))
stdev<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=rownames(all10kDR))
cv<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=rownames(all10kDR))
sterror<-data.frame(matrix(NA, nrow = 5911, ncol = 1), row.names=rownames(all10kDR))

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
}


all10kDR_summary <- cbind(harmMeans,medians,means,Count.of.NA,range,variance,stdev,cv,sterror,deparse.level=1)

colnames(all10kDR_summary) <- c("harmMeans","medians","means","Count.of.NA","range","variance","stdev","cv","sterror")
rownames(all10kDR_summary) <- all10kDR$V1

write.table(all10kDR_summary, "DR-SUMMARY_MamPhy_BDvarRates_17Uni_all10ktrees_5910sp_expanded.txt", col.names=TRUE, row.names=TRUE)


###
all10kDR<-fread("DR-matrix_MamPhy_BDvarRates_17Uni_all10ktrees_pruned4098.txt", header="auto") # 5911 obs. of  10000 variables

# Will want the mean, median, mode values across the ROWS for each tip in the tree >> 
harmMeans<-data.frame(matrix(NA, nrow = 4099, ncol = 1), row.names=rownames(all10kDR))
medians<-data.frame(matrix(NA, nrow = 4099, ncol = 1), row.names=rownames(all10kDR))
means<-data.frame(matrix(NA, nrow = 4099, ncol = 1), row.names=rownames(all10kDR))
Count.of.NA<-data.frame(matrix(NA, nrow = 4099, ncol = 1), row.names=rownames(all10kDR))
range<-data.frame(matrix(NA, nrow = 4099, ncol = 1), row.names=rownames(all10kDR))
variance<-data.frame(matrix(NA, nrow = 4099, ncol = 1), row.names=rownames(all10kDR))
stdev<-data.frame(matrix(NA, nrow = 4099, ncol = 1), row.names=rownames(all10kDR))
cv<-data.frame(matrix(NA, nrow = 4099, ncol = 1), row.names=rownames(all10kDR))
sterror<-data.frame(matrix(NA, nrow = 4099, ncol = 1), row.names=rownames(all10kDR))

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
}


all10kDR_summary <- cbind(harmMeans,medians,means,Count.of.NA,range,variance,stdev,cv,sterror,deparse.level=1)

colnames(all10kDR_summary) <- c("harmMeans","medians","means","Count.of.NA","range","variance","stdev","cv","sterror")
rownames(all10kDR_summary) <- all10kDR$V1

write.table(all10kDR_summary, "DR-SUMMARY_MamPhy_BDvarRates_17Uni_all10ktrees_pruned4098sp_expanded.txt", col.names=TRUE, row.names=TRUE)

######


q()

n


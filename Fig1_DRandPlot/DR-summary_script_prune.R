library(ape)
library(data.table)

#setwd("~/MamPhy_analyses")
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors")


all10kDR<-read.table("DR-matrix_MamPhy_BDvrFIXED_FBD_all10k.txt", header=TRUE)

# Will want the mean, median, mode values across the ROWS for each tip in the tree >> 
harmMeans<-data.frame(matrix(NA, nrow = 5912, ncol = 1), row.names=rownames(all10kDR))
medians<-data.frame(matrix(NA, nrow = 5912, ncol = 1), row.names=rownames(all10kDR))
means<-data.frame(matrix(NA, nrow = 5912, ncol = 1), row.names=rownames(all10kDR))

for (i in 1:length(rownames(all10kDR))){
	harmMeans[i,] <- 1/(mean(1/(unlist(all10kDR[i,]))))
	medians[i,] <- median(unlist(all10kDR[i,]))
	means[i,] <- mean(unlist(all10kDR[i,]))
}

all10kDR_summary <- cbind(harmMeans,medians,means, deparse.level=1)
colnames(all10kDR_summary)<-c("harmMeans","medians","means")

write.table(all10kDR_summary, "DR-SUMMARY_MamPhy_BDvrFIXED_FBD_all10k_5912species.txt")


#######
######
# And the same for NDexp
##

all10kDR<-read.table("DR-matrix_MamPhy_BDvrFIXED_NDexp_all10k.txt", header=TRUE)

# Will want the mean, median, mode values across the ROWS for each tip in the tree >> 
harmMeans<-data.frame(matrix(NA, nrow = 5912, ncol = 1), row.names=rownames(all10kDR))
medians<-data.frame(matrix(NA, nrow = 5912, ncol = 1), row.names=rownames(all10kDR))
means<-data.frame(matrix(NA, nrow = 5912, ncol = 1), row.names=rownames(all10kDR))

for (i in 1:length(rownames(all10kDR))){
	harmMeans[i,] <- 1/(mean(1/(unlist(all10kDR[i,]))))
	medians[i,] <- median(unlist(all10kDR[i,]))
	means[i,] <- mean(unlist(all10kDR[i,]))
}

all10kDR_summary <- cbind(harmMeans,medians,means, deparse.level=1)
colnames(all10kDR_summary)<-c("harmMeans","medians","means")

write.table(all10kDR_summary, "DR-SUMMARY_MamPhy_BDvrFIXED_NDexp_all10k_5912species.txt")

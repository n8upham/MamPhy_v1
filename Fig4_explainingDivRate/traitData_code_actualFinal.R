





######
# Load data
library(ape); library(phytools)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"

mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")

tipDataAll_wComments_DETAILS<-read.table("MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass.txt", header=TRUE)
tipDataAll_1<-tipDataAll_wComments_DETAILS[,c("GENUS","FAMILY","ORDER","perVert","perInvert","perPlant","trophicCat","CarnOrNot","HerbOrNot","OmniOrNot","lifemodeCat6","AquaOrNot","ArboOrNot","FlysOrNot","SubtOrNot","TerrOrNot","activityCat3","NoctOrNot","CathOrNot","DiurOrNot","tipDR_mean10k","geoRange_numGrid","geoRange_km2","BM_final_kg","BM_final_g")]
tipDataAll_2<-tipDataAll_wComments_DETAILS[,c(52:81,83:101)]
tipDataAll<-cbind.data.frame(tipDataAll_1,tipDataAll_2)
rownames(tipDataAll)<-tipDataAll_wComments_DETAILS$tiplabel

# Re-calculating the home range to mass relationship in Jetz et al. 2005

datHR_all<-as.data.frame(na.omit(tipDataAll[,c("trophicCat","X22.2_HomeRange_Indiv_km2","BM_final_kg")])) # exclude "X22.1_HomeRange_km2" not corrected for group size
colnames(datHR_all)<-c("trophicCat","HomeRange_Indiv_km2","BodyMass_kg")
	# gives 606 entries

pdf(file="reCalc_of_log-log_HR-vs-BM_mammalsTrophic_4way.pdf", onefile=TRUE, height=8, width=8)

cats<-c("All mammals:", "Herbivores:", "Omnivores:", "Carnivores:")
dats<-list(datHR_all,datHR_all[which(datHR_all$trophicCat=="Herb"),],datHR_all[which(datHR_all$trophicCat=="Omni"),],datHR_all[which(datHR_all$trophicCat=="Carn"),])

HR_equations<-data.frame(matrix(NA, nrow = length(cats), ncol = 5))
colnames(HR_equations)<-c("log10_Int_km2","log10_Slope_kg","nonLog_Int_ha","nonLog_Slope_Kg","R2")
rownames(HR_equations)<-cats

layout(matrix(c(1:4), 2, 2, byrow = TRUE))#, widths=4, heights=4)
par(oma = c(5,4,5,3) + 0.1, mar = rep(1.7,4) + 0.1)

for(i in 1:length(cats)){
	cat<-cats[i]
	dat<-dats[[i]]

	form<-log10(HomeRange_Indiv_km2) ~ log10(BodyMass_kg)
	plot(form, data=dat)
	fit<-lm(form,data=dat)
	sum<-summary(fit)

	n=length(dat[,1])
	a=round(sum$coef[1],2)
	b=round(sum$coef[2],2)
	R2=round(sum$r.squared,2)
	X=seq(min(log10(dat$BodyMass_kg)),max(log10(dat$BodyMass_kg)),length.out=20)
	lines(x=X,y=b*X+a,col="black", lwd=4,lty=1) 

	mtext(side=3, text=bquote(bold(.(cat) ~ n ~ "=" ~ .(n) ~~~~ R^2 ~ "=" ~ .(R2))),adj=0, cex=0.7)

	Y=log10(dat$HomeRange_Indiv_km2)
	text(min(X),max(Y),label=bquote(log10 ~ km^2 ~ bold(y ~ "=" ~ .(a)+.(b)*x)~log10 ~ kg),adj=0, cex=0.8)

	text(min(X),max(Y)-(max(Y)-min(Y))/12,label=bquote(ha ~ bold(y ~ "=" ~ .(round(100*10^a,2))+.(round(10^b,2))*x)~kg),adj=0, cex=0.8)

	HR_equations[i,1]<-a
	HR_equations[i,2]<-b
	HR_equations[i,3]<-round(100*10^a,2)
	HR_equations[i,4]<-round(10^b,2)
	HR_equations[i,5]<-R2

}
title(main="", xlab = "log10(Body Mass in kg)",
      ylab = "log10(Home Range per indiv in km2)",
      outer = TRUE, line = 3,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.5,font.axis=1, font.lab=2)


dev.off()

write.table(HR_equations,file="reCalc_of_log-log_HR-vs-BM_mammalsTrophic_4way_equations.txt")

## Now do those HR calcs
datHR_all<-as.data.frame(na.omit(tipDataAll[,c("trophicCat","BM_final_kg")])) # exclude "X22.1_HomeRange_km2" not corrected for group size
colnames(datHR_all)<-c("trophicCat","BodyMass_kg")
	# gives 4988 entries

cats<-c("All mammals:", "Herbivores:", "Omnivores:", "Carnivores:")
dats<-list(datHR_all,datHR_all[which(datHR_all$trophicCat=="Herb"),],datHR_all[which(datHR_all$trophicCat=="Omni"),],datHR_all[which(datHR_all$trophicCat=="Carn"),])

HR_ext_vals<-data.frame(matrix(NA, nrow = length(tipDataAll[,1]), ncol = 1))
rownames(HR_ext_vals)<-tipDataAll_wComments_DETAILS$tiplabel

for(i in 2:length(cats)){
	cat<-cats[i]
	dat<-dats[[i]]

	a=HR_equations[i,1]
	b=HR_equations[i,2]

	for(j in 1:length(dat[,1])){
		HR_j <- a+(log(dat[j,"BodyMass_kg"])*b)
		names(HR_j)<-rownames(dat[j,])
		HR_ext_vals[which(rownames(HR_ext_vals)==names(HR_j)),]<-HR_j
	}
}
colnames(HR_ext_vals)<-"log10homeRange_km2_ext"

homeRange_km2_ext<-10^HR_ext_vals
colnames(homeRange_km2_ext)<-"homeRange_km2_ext"

HR_only<-cbind(HR_ext_vals,homeRange_km2_ext)
write.table(HR_only,file="MamPhy_5911sp_homeRangeExt_only.txt")

tipDataAllHR<-cbind(tipDataAll,HR_ext_vals,homeRange_km2_ext)
write.table(tipDataAllHR,file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_homeRangeExt.txt")

####
# Now do the same extrapolation for the DISPERSAL DISTANCE
library(ape); library(phytools)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"

tipDataAllHR<-read.table(file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_homeRangeExt.txt")

# BIN the home range data also...
lnHomeRange_km2_ext<-as.vector(log(tipDataAllHR$homeRange_km2_ext))
names(lnHomeRange_km2_ext)<-rownames(tipDataAllHR)
#binsEven
	bins<-quantile(lnHomeRange_km2_ext,probs=seq(0,1,by=0.1),na.rm=TRUE)
	lnHomeRange_rank1to10_ext<-data.frame(matrix(NA, nrow = length(tipDataAllHR[,1]), ncol = 1))
	rownames(lnHomeRange_rank1to10_ext)<-rownames(tipDataAllHR)
	colnames(lnHomeRange_rank1to10_ext)<-"lnHomeRange_rank1to10"

	toBin<-as.matrix(na.omit(lnHomeRange_km2_ext))

	for(i in 1:length(toBin)){
		for(j in 1:(length(bins)-1)){
			if(toBin[[i]] >= bins[[j]] && toBin[[i]] <= bins[[j+1]]){
			lnHomeRange_rank1to10_ext[which(rownames(lnHomeRange_rank1to10_ext)==names(toBin[i,])),]<-j} else {next}
		}
	}

#binsEqual
	range<-range(na.omit(lnHomeRange_km2_ext))
	binsEqual<-seq(from=range[1],to=range[2],length.out=11)
	lnHomeRange_rank1to10_binsEqual<-data.frame(matrix(NA, nrow = length(tipDataAllHR[,1]), ncol = 1))
	rownames(lnHomeRange_rank1to10_binsEqual)<-rownames(tipDataAllHR)
	colnames(lnHomeRange_rank1to10_binsEqual)<-"lnHomeRange_rank1to10_binsEqual"

	toBin<-as.matrix(na.omit(lnHomeRange_km2_ext))

	for(i in 1:length(toBin)){
		for(j in 1:(length(binsEqual)-1)){
			if(toBin[[i]] >= binsEqual[j] && toBin[[i]] <= binsEqual[j+1]){
			lnHomeRange_rank1to10_binsEqual[which(rownames(lnHomeRange_rank1to10_binsEqual)==names(toBin[i,])),]<-j} else {next}
		}
	}

#hist(lnDispDistAll_rank1to10[,1],breaks=30)
pdf(file="lnHomeRange_km2_ext_histoAsBinned.pdf")
hist(toBin,breaks=30, col="grey")
for(j in 1:(length(bins)-1)){
	abline(v=bins[[j+1]], col="blue",lwd=2)
}
dev.off()

pdf(file="lnHomeRange_rank1to10_binsEqual_histoAsBinned.pdf")
hist(toBin,breaks=30, col="grey")
for(j in 1:(length(binsEqual)-1)){
	abline(v=binsEqual[j+1], col="blue",lwd=2)
}
dev.off()



##########
# calc dispersal distance...
datForDispersal<-as.data.frame(tipDataAllHR[,c("BM_final_g","homeRange_km2_ext","geoRange_km2")])
	# gives 4961 values

# RE-calc for dispersal
library(ape); library(phytools)
setwd("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"

tipDataAllHR<-read.table(file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)
datForDispersal<-as.data.frame(tipDataAllHR[,c("BM_final_g","homeRange_km2_ext","geoArea_km2")])

# Dispersal distance = -1.153 + 0.315*bodyMass(g) + 0.220*homeRange(km2) + 0.252*geoRange(km2)

lnBM<-log(datForDispersal[,"BM_final_g"])
lnHR<-log(datForDispersal[,"homeRange_km2_ext"])
lnGR<-log(datForDispersal[,"geoArea_km2"])
lnDispDistAll_km_ext = -1.153 + 0.315*lnBM + 0.220*lnHR + 0.252*lnGR
names(lnDispDistAll_km_ext)<-tipDataAllHR$tiplabel

DispDistAll_km_ext <- exp(lnDispDistAll_km_ext)

#binsEven
	bins<-quantile(lnDispDistAll_km_ext,probs=seq(0,1,by=0.1),na.rm=TRUE)
	lnDispDistAll_rank1to10_ext<-data.frame(matrix(NA, nrow = length(lnDispDistAll_km_ext), ncol = 1))
	rownames(lnDispDistAll_rank1to10_ext)<-names(lnDispDistAll_km_ext)
	colnames(lnDispDistAll_rank1to10_ext)<-"lnDispDistAll_rank1to10"

	toBin<-na.omit(lnDispDistAll_km_ext)

	for(i in 1:length(toBin)){
		for(j in 1:(length(bins)-1)){
			if(toBin[[i]] >= bins[[j]] && toBin[[i]] <= bins[[j+1]]){
			lnDispDistAll_rank1to10_ext[which(rownames(lnDispDistAll_rank1to10_ext)==names(toBin[i])),]<-j} else {next}
		}
	}

#binsEqual
	range<-range(na.omit(lnDispDistAll_km_ext))
	binsEqual<-seq(from=range[1],to=range[2],length.out=11)
	lnDispDistAll_rank1to10_binsEqual<-data.frame(matrix(NA, nrow = length(tipDataAllHR[,1]), ncol = 1))
	rownames(lnDispDistAll_rank1to10_binsEqual)<-tipDataAllHR$tiplabel
	colnames(lnDispDistAll_rank1to10_binsEqual)<-"lnDispDistAll_rank1to10_binsEqual"

	toBin<-as.matrix(na.omit(lnDispDistAll_km_ext))

	for(i in 1:length(toBin)){
		for(j in 1:(length(binsEqual)-1)){
			if(toBin[[i]] >= binsEqual[j] && toBin[[i]] <= binsEqual[j+1]){
			lnDispDistAll_rank1to10_binsEqual[which(rownames(lnDispDistAll_rank1to10_binsEqual)==names(toBin[i,])),]<-j} else {next}
		}
	}

#hist(lnDispDistAll_rank1to10[,1],breaks=30)
pdf(file="lnDispDistAll_km_ext_histoAsBinned.pdf")
hist(toBin,breaks=30, col="grey")
for(j in 1:(length(bins)-1)){
	abline(v=bins[[j+1]], col="blue",lwd=2)
}
dev.off()

pdf(file="lnDispDistAll_rank1to10_binsEqual_histoAsBinned.pdf")
hist(toBin,breaks=30, col="grey")
for(j in 1:(length(binsEqual)-1)){
	abline(v=binsEqual[j+1], col="blue",lwd=2)
}
dev.off()


DispDist_only<-cbind(lnBM,lnHR,lnGR,lnDispDistAll_km_ext,DispDistAll_km_ext,lnDispDistAll_rank1to10_ext,lnDispDistAll_rank1to10_binsEqual)
write.table(DispDist_only,file="MamPhy_5911sp_dispersalDistance_Ext_only_corrected.txt")

DispDist_only<-cbind(lnBM,lnHR,lnGR,lnDispDistAll_km_ext,DispDistAll_km_ext,lnDispDistAll_rank1to10_ext,lnDispDistAll_rank1to10_binsEqual,lnHomeRange_rank1to10_ext,lnHomeRange_rank1to10_binsEqual)
write.table(DispDist_only,file="MamPhy_5911sp_dispersalDistance_Ext_only.txt")

tipDataAll_ext<-cbind(tipDataAll,HR_ext_vals,homeRange_km2_ext,lnDispDistAll_km_ext,DispDistAll_km_ext,lnDispDistAll_rank1to10_ext)
write.table(tipDataAll_ext,file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp.txt")

#####
# validate the home range calcs...
library(ape); library(phytools)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"

tipDataAll_ext<-read.table(file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp.txt")

dat<-na.omit(tipDataAll_ext[,c("homeRange_km2_ext","X22.2_HomeRange_Indiv_km2")])
	# gives 606 species, as expected.
form<-log(homeRange_km2_ext) ~ log(X22.2_HomeRange_Indiv_km2)
plot(form, data=dat)
fit<-lm(form,data=dat)
sum<-summary(fit)
	# Adjusted R-squared:  0.734 >> that makes sense, follows from the 3 prediction equations used.  Cool.


###
# Look at PAIRWISE with tip DR relationships with these new variables
##
library(nlme); library(geiger)
tipDataAll<-read.table(file="MamPhy_5911sp_tipDR-range-Pantheria-EltonTraits-mass_extended_HR_Disp_ALL-comments.txt", header=TRUE)
mamMCC<-read.tree(file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target_5911sp_newick.tre")


# as geoArea
dat1<-tipDataAll[which(tipDataAll$TerrOrMar=="Terrestrial"),c("tipDR_mean10k","geoArea_km2")]
rownames(dat1)<-tipDataAll[which(tipDataAll$TerrOrMar=="Terrestrial"),"tiplabel"]
cladeData<-treedata(mamMCC,na.omit(dat1)) 
dat<-as.data.frame(cladeData$data)

		#dat<-na.omit(tipDataAll[,c("tipDR_mean10k","GenerationLength_d")])
form<-(tipDR_mean10k ~ log(geoArea_km2))
plot(form, data=tipDataAll)
fit<-gls(form, data=dat)
sum<-summary(fit)

# as Lat of centroid
dat1<-tipDataAll[,c("tipDR_mean10k","Lat_centroid")]
rownames(dat1)<-tipDataAll$tiplabel
cladeData<-treedata(mamMCC,na.omit(dat1)) 
dat<-as.data.frame(cladeData$data)

pdf(file="tipDR_mean10k_vs_Lat_centroid_plot.pdf")
form<-(tipDR_mean10k ~ (Lat_centroid))
plot(form, data=tipDataAll)
fit<-gls(form, data=dat)
sum<-summary(fit)
dev.off()


# as dispersalDistance
dat1<-tipDataAll[,c("tipDR_mean10k","DispDistAll_km_ext")]
rownames(dat1)<-tipDataAll$tiplabel
cladeData<-treedata(mamMCC,na.omit(dat1)) 
dat<-as.data.frame(cladeData$data)

pdf(file="logTipDR_mean10k_vs_LogDispDistAll_plot.pdf")
form<-(log(tipDR_mean10k) ~ log(DispDistAll_km_ext))
plot(form, data=dat)
fit<-gls(form, data=dat)
sum<-summary(fit)
dev.off()


fitPagelModel<-gls(form, data=dat, correlation=corPagel(value=1,phy=cladeData$phy))
sumPagel<-summary(fitPagelModel)



# as continuous...
dat1<-tipDataAll[,c("tipDR_mean10k","GenerationLength_d")]
rownames(dat1)<-tipDataAll$tiplabel
cladeData<-treedata(mamMCC,na.omit(dat1)) 
dat<-as.data.frame(cladeData$data)

		#dat<-na.omit(tipDataAll[,c("tipDR_mean10k","GenerationLength_d")])
form<-(tipDR_mean10k ~ log(GenerationLength_d))
plot(form, data=tipDataAll)

fit<-gls(form, data=dat, correlation=corPagel(value=1,phy=cladeData$phy))
sum<-summary(fit)


dat<-na.omit(tipDataAll[,c("tipDR_mean10k","lnHomeRange_rank1to10")])
form<-(tipDR_mean10k ~ (lnHomeRange_rank1to10))
plot(form, data=tipDataAll)

fit<-lm(form, data=dat)
sum<-summary(fit)


dat<-na.omit(tipDataAll[,c("tipDR_mean10k","trophicCat")])
form<-(tipDR_mean10k ~ trophicCat)
plot(form, data=tipDataAll)

fit<-gls(form, data=dat)
sum<-summary(fit)



	a=round(sum$coef[1],2)
	b=round(sum$coef[2],2)
	R2=round(sum$r.squared,2)
	X=seq(min(log10(dat$BodyMass_kg)),max(log10(dat$BodyMass_kg)),length.out=20)
	lines(x=X,y=b*X+a,col="black", lwd=4,lty=1) 

	mtext(side=3, text=bquote(bold(.(cat) ~ n ~ "=" ~ .(n) ~~~~ R^2 ~ "=" ~ .(R2))),adj=0, cex=0.7)

	Y=log10(dat$HomeRange_Indiv_km2)
	text(min(X),max(Y),label=bquote(log10 ~ km^2 ~ bold(y ~ "=" ~ .(a)+.(b)*x)~log10 ~ kg),adj=0, cex=0.8)



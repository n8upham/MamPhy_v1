#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy MS1 -- Upham et al. 2019 -- PLOS Biology
###
# Figure 3 - comparing divergence times of node-and tip-dated Mammalia backbones
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Now, COMPARE the FBD and ND divergences in a pairwise fashion
########
setwd("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")
library(plotrix)

div59<-read.table("compare_backbones_ND-FBD_59taxa.txt",header=TRUE)
div28<-read.table("compare_backbones_ND-FBD_28taxa.txt",header=TRUE)

whichRows<-c(3:(nrow(div28)))
yMegaMin<-min(div28[whichRows,"FBD_min"], na.rm=TRUE)
#whichRows<-c(1:nrow(div28))
#yMegaMin<-0

yMegaMax<-max(div28[whichRows,"FBD_max"], na.rm=TRUE)
xEst<-div28[whichRows,"NDexp_mean"]
yEst<-div28[whichRows,"FBD_mean"]
xMin<-div28[whichRows,"NDexp_min"]
xMax<-div28[whichRows,"NDexp_max"]
yMin<-div28[whichRows,"FBD_min"]
yMax<-div28[whichRows,"FBD_max"]

#pdf(file="compare_backbones_ND-FBD_pairwiseDIVs_28taxa_theria.pdf")
#pdf(file="compare_backbones_ND-FBD_pairwiseDIVs_28taxa.pdf")
pdf(file="compare_backbones_ND-FBD_pairwiseDIVs_28taxa_theria_wOutline_2.pdf")

plot(data=div28, FBD_mean ~ NDexp_mean, xlim=c(yMegaMin,yMegaMax), ylim=c(yMegaMin,yMegaMax), pch=NA, bty="l", xlab="Node-dated", ylab="Tip-dated" )
abline(a=0,b=1,lty=2)

plotCI(x=xEst,y=yEst, ui=xMax,li=xMin,err="x",
	sfrac=0.0,gap=0,slty=par("lty"),
	add=TRUE,scol="black")

plotCI(x=xEst,y=yEst, ui=yMax,li=yMin,err="y",
	sfrac=0.0,gap=0,slty=par("lty"),
	add=TRUE,scol="black")

dev.off()


plotCI(x=xEst,y=1:length(xEst), ui=xMax,li=xMin,err="x", xlim=c(0,150),
	sfrac=0.0,gap=0)

plotCI(x=yEst,y=1:length(yEst)+0.4, ui=yMax,li=yMin,err="x",
	sfrac=0.0,gap=0, col="red",
	add=TRUE)


####
# Then compare the TIP DR values of the FBD and ND...
########
###
setwd("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/11_makingFullPosteriors/10k_tree_andDR_distributions_AugDec2018")
library(dplyr)
	
drFBD<-read.table("DR-SUMMARY_MamPhy_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_all10k_expanded2.txt",header=TRUE)
drNDe<-read.table("DR-SUMMARY_MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_all10k_expanded2.txt",header=TRUE)[c(2:5912),]

setwd("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")
mamData<-read.csv(file="traits_datReady_mamPhy_5911species.csv",header=TRUE)[,c(1:11)]

tally<-rep(NA,5911)
for(i in 1:5911){
if(rownames(drFBD)[i] == rownames(drNDe)[i]){
	tally[i]<-0
} else
	tally[i]<-1
} # many different, so need to JOIN

drJOINED<-left_join(x=cbind.data.frame(tiplabel=rownames(drNDe), drNDe),y=cbind.data.frame(tiplabel=rownames(drFBD),drFBD),by="tiplabel",
	suffix= c("_NDexp","_FBD"))
allJOINED<-left_join(x=mamData,y=drJOINED,by="tiplabel")
write.table(allJOINED, file="compare_backbones_ND-FBD_tipDR_5911sp_TABLE.txt")




# get COLORS
library(gplots)
cats<-names(sort(table(allJOINED$higher),decreasing=FALSE))[c(4,3,2,5,6)]

colsOrder<-rich.colors(length(cats)+1, alpha=1)
#cols<-c(colsOrder[6],colsOrder[4],colsOrder[2],colsOrder[5],colsOrder[3])
cols<-c(colsOrder[6],colsOrder[3],colsOrder[3],colsOrder[5],"blueviolet")
boxCols<-c(colsOrder[6],colsOrder[3],colsOrder[5],"blueviolet")

# get data ready
datReady<-allJOINED


png(file="compare_backbones_ND-FBD_tipDR_5911sp_1to1Line_w95CI_noLeg.png", height=6, width=5.5, units="in",res=600)
#pdf(file="compare_backbones_ND-FBD_tipDR_5911sp_1to1Line_w95CI_forLeg.pdf", height=6, width=5.5, onefile=TRUE)

#par(oma = rep(3,4) + 0.1,
#    mar = c(4,1.5,4,1) + 0.1) #‘c(bottom, left, top, right)’

#layout(matrix(1:2,1,2, byrow=TRUE),widths=c(0.25,0.75),heights=c(1,1))

yMegaMin<-0
yMegaMax<-max(datReady[,"high95_FBD"])-0.3

xEst<-datReady[,"harmMeans_NDexp"]
yEst<-datReady[,"harmMeans_FBD"]

xMin<-datReady[,"low95_NDexp"]
xMax<-datReady[,"high95_NDexp"]
yMin<-datReady[,"low95_FBD"]
yMax<-datReady[,"high95_FBD"]


# FBD ~ NDexp
	dat<-datReady
	form<-as.formula(harmMeans_FBD ~ harmMeans_NDexp)

	pointCols<-rep("black",length(dat[,1]))
	for(j in 1:length(cats)){
	pointCols[which(as.vector(dat$higher)==cats[j])]<-cols[j]
	}
	otherCI<-grey(0.5,alpha=0.1) #rgb(1,0,0,alpha=0.1)

	plot(form, data=dat, xlim=c(0,yMegaMax),ylim=c(0,yMegaMax), xlab="",ylab="",col=NA, cex.lab=1,cex.axis=1, bty="l",xpd=NA)
    plotCI(add=TRUE, y=yEst, x=xEst, err="y", li=yMin, ui=yMax, sfrac=0, gap=0,col=otherCI, pch=NA)
    plotCI(add=TRUE,y=yEst, x=xEst, err="x", li=xMin, ui=xMax, sfrac=0, gap=0,col=rgb(0,0,1,alpha=0.05), pch=NA)
 	abline(a=0,b=1,lty=2, lwd=2)
 	points(form, data=dat,col=pointCols, cex.lab=1,cex.axis=1, pch=".")

	mtext(text="Tip-dated",side=2,line=3, font=1)
	#mtext(text="1,000 trees",side=1,line=4.5, font=1, cex=1)

	mtext(text="Node-dated",side=1,line=3, font=1)
	#mtext(text="10,000 trees",side=2,line=2.5, font=1, cex=1)

#	model<-summary(lm(form, data=dat))
#	a<-round(model$coef[[1]], 2)
#	b<-round(model$coef[[2]], 2)
#	r2<-round(model$r.squared,digits=2)
#	X=seq(min(dat[,"harmMeans_NDexp"]),max(dat[,"harmMeans_NDexp"]),length.out=20)
#	lines(x=X,y=b*X+a,col="black", lwd=3,lty=1) #x1=min(dat$MRCA),x2=max(dat$MRCA),y1=min(log(dat$richness)),y2=max(log(dat$richness)),untf=FALSE,
#
#	abline(a=0, b=max(dat[,"harmMeans_FBD"])/max(dat[,"harmMeans_NDexp"]), lty=2, lwd=2, col="black")
#	#text(labels=bquote("y = "*.(a)*" + "*.(b)*"x; "*r^2*" = "*.(r2)), x=xMax-(xMax/2.5), y=0.02, col="black", cex=1.1)
#	r<-round(cor(y=dat$harmMeans_FBD, x=dat[,"harmMeans_NDexp"], method="spearman"),digits=2)
#	mtext(text=bquote("spearman "*r*" = "*.(r) ), side=3, col="black", cex=0.8, line =1.5 )
#	mtext(text=bquote("y = "*.(a)*" + "*.(b)*"x, "*R^2*" = "*.(r2) ), side=3, col="black", cex=0.8)

# legend
	#plot(harmMeans_FBD ~ harmMeans_NDexp, data=datReady, xlim=c(0,10),ylim=c(0,1), type="n", xlab="",ylab="", bty="n", xaxt="n", yaxt="n", cex.lab=1.2,cex.axis=1.3)
#	legend(box.col="white",x=0.01,y=1.5,legend=c( "Marsupialia","Afrotheria / Xenarthra", "Laurasiatheria", "Euarchontoglires"),fill=boxCols, cex=0.8, xpd="NA")
#	par(new=TRUE)
#	legend(box.col="white",x=0.57,y=1.5,legend=c("1-to-1 line", #"Regression line", 
#				"95% CI (node-dated)","95% CI (tip-dated)"),
#				lty=c(2,1,#1,
#					1), lwd=c(2,1,#1,
#					1),col=c("black","blue",#"blue",
#					"grey"), cex=0.8, xpd="NA")

dev.off()









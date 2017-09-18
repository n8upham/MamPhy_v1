library(ape)
library(TreeSim)

setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors")

mamFBD10k<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_all10k_nexus.trees")

mamNDexp10k<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_all10k_nexus.trees")

mamNDexp_1000<-sample(x=mamNDexp10k,size=1000)

write.tree(mamNDexp_1000,file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample1000.trees")
write.nexus(mamNDexp_1000,file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample1000_nexus.trees")


##
# plot the FULL 10k versions
###
jpeg(file="LTT_MamPhy_BDvr_FBD_5911spp_all10k_noOut.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamFBD10k[[1]],"_Anolis_carolinensis"), log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:length(mamFBD10k)) {
	ltt.lines(drop.tip(mamFBD10k[[i]],"_Anolis_carolinensis"), col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}

title(main="MamPhy, BDvr, 5911 spp, all10k - FBD (blue)")

dev.off()


jpeg(file="LTT_MamPhy_BDvr_NDexp_5911spp_all10k_noOut.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamNDexp10k[[1]],"_Anolis_carolinensis"), log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:length(mamNDexp10k)) {
	ltt.lines(drop.tip(mamNDexp10k[[i]],"_Anolis_carolinensis"), col=hsv(0.9708995,0.2470588,1,alpha=0.5), lty = 2, lwd=1)
}

title(main="MamPhy, BDvr, 5911 spp, all10k - ND (pink)")

dev.off()


jpeg(file="LTT_MamPhy_BDvr_FBD-vs-NDexp_5911spp_all10k_noOut.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamFBD10k[[1]],"_Anolis_carolinensis"), log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:length(mamFBD10k)) {
	ltt.lines(drop.tip(mamFBD10k[[i]],"_Anolis_carolinensis"), col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}
for (i in 1:length(mamNDexp10k)) {
	ltt.lines(drop.tip(mamNDexp10k[[i]],"_Anolis_carolinensis"), col=hsv(0.9708995,0.2470588,1,alpha=0.5), lty = 2, lwd=1)
}

title(main="MamPhy, BDvr, 5911 spp, all10k - FBD (blue) vs ND (pink)")

dev.off()


jpeg(file="LTT_MamPhy_BDvr_NDexp-vs-FBD_5911spp_all10k_noOut.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamNDexp10k[[1]],"_Anolis_carolinensis"), log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col=hsv(0.9708995,0.2470588,1,alpha=0.5))
for (i in 2:length(mamNDexp10k)) {
	ltt.lines(drop.tip(mamNDexp10k[[i]],"_Anolis_carolinensis"), col=hsv(0.9708995,0.2470588,1,alpha=0.5), lty = 2, lwd=1)
}
for (i in 1:length(mamFBD10k)) {
	ltt.lines(drop.tip(mamFBD10k[[i]],"_Anolis_carolinensis"), col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}

title(main="MamPhy, BDvr, 5911 spp, all10k - FBD (blue) vs ND (pink)")

dev.off()

###
# Now SAMPLE out 100 random trees
##
library(ape)
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors")

# FBD
trees1=scan("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_all10k.trees", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

phy_text100<-sample(trees1,size=100)
ntrees<-length(phy_text100)

mamFBD_100<-vector("list",ntrees)
for (i in 1:length(phy_text100)){
	phy<-drop.tip(read.tree(text=phy_text100[[i]]),"_Anolis_carolinensis")
	write.tree(phy, "MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100.trees", append=TRUE)
	mamFBD_100[[i]]<-phy
}
write.nexus(mamFBD_100, file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_nexus.trees")

# NDexp
trees1=scan("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_all10k.trees", what="",sep="\n",quiet=TRUE,skip=0,comment.char="#") 

phy_text100<-sample(trees1,size=100)
ntrees<-length(phy_text100)

mamNDexp_100<-vector("list",ntrees)
for (i in 1:length(phy_text100)){
	phy<-drop.tip(read.tree(text=phy_text100[[i]]),"_Anolis_carolinensis")
	write.tree(phy, "MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100.trees", append=TRUE)
	mamNDexp_100[[i]]<-phy
}
write.nexus(mamNDexp_100, file="MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_nexus.trees")

# plot the LTT compare of these...
mamFBD_100<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_nexus.trees")
mamNDexp_100<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_nexus.trees")

jpeg(file="LTT_MamPhy_BDvr_FBD-vs-NDexp_5911spp_sample100_noOut.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(mamFBD_100[[1]], log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:length(mamFBD_100)) {
	ltt.lines(mamFBD_100[[i]], col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}
for (i in 1:length(mamNDexp_100)) {
	ltt.lines(mamNDexp_100[[i]], col=hsv(0.9708995,0.2470588,1,alpha=0.5), lty = 2, lwd=1)
}

title(main="MamPhy, BDvr, 5911 spp, sample100 - FBD (blue) vs ND (pink)")

dev.off()


#######
# Now with those 100 random trees...
library(ape)
library(phytools)
setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors")
#setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")

mamFBD_100<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_nexus.trees")
mamNDexp_100<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_nexus.trees")

cladesDR<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_DR-SUMMARY-expanded_FBD.txt")
head(cladesDR)
colnames(cladesDR)<-c("tiplabel","gen","fam","ord","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")

###
# Doing ltt95 curves with filled polygons somehow...

marsupials<-cladesDR[which(cladesDR$PC=="PC1_Marsupials"),"tiplabel"]
monotremes<-cladesDR[which(cladesDR$PC=="PC23_Monotremata"),"tiplabel"]
marMono<-c(as.vector(marsupials),as.vector(monotremes))
placentals<-setdiff(mamFBD_100[[1]]$tip.label,marMono)
placMono<-c(as.vector(placentals),as.vector(monotremes))

mamFBD_100_marsupials<-vector("list",length(mamFBD_100))
xy_all_M<-vector("list", length(mamFBD_100_marsupials))str
for (i in 1:length(mamFBD_100_marsupials)){
	mamFBD_100_marsupials[[i]]<-drop.tip(mamFBD_100[[i]],placMono)
	xy_all_M[[i]] <- ltt.plot.coords(mamFBD_100_marsupials[[i]])
}
xy_all_M_appended<-do.call(rbind, xy_all_M)
write.nexus(mamFBD_100_marsupials, file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_nexus_Marsupials.trees")

mamFBD_100_marsupials<-read.nexus(file="MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_nexus_Marsupials.trees")

###
plot.ltt95<-function(x,...){
	if(hasArg(log)) log<-list(...)$log
	else log<-attr(x,"log")
	if(hasArg(xaxis)) xaxis<-list(...)$xaxis
	else xaxis<-"standard"
	if(attr(x,"method")=="times"){
		n<-max(x[,1])
		if(xaxis=="negative"){ 
			x[,2:4]<-x[,2:4]-max(x[,2:4])
		}
		if(xaxis=="flipped"){
			x[,2:4]<-max(x[,2:4])-x[,2:4]
			x.lim<-c(max(x[,2:4]),min(x[,2:4]))
		} else x.lim<-range(x[,2:4])
		plot(x[,3],x[,1],lwd=2,xlim=x.lim,
			type=if(attr(x,"mode")=="median") "s" else "l",main=NULL,
			xlab="time",ylab="lineages",log=if(log) "y" else "")
		lines(x[,2],x[,1],lty="dashed",type=if(attr(x,"mode")=="median") "s" else "l")
		lines(x[,4],x[,1],lty="dashed",type=if(attr(x,"mode")=="median") "s" else "l")	
	} else if(attr(x,"method")=="lineages"){
		if(xaxis=="negative") x[,1]<-x[,1]-max(x[,1])
		if(xaxis=="flipped"){
			x[,1]<-max(x[,1])-x[,1]
			x.lim<-c(max(x[,1]),min(x[,1]))
		} else x.lim<-range(x[,1])
		plot(x[,1],x[,3],xlim=x.lim,ylim=c(min(x[,2]),max(x[,4])),lwd=2,
			type=if(attr(x,"mode")=="median") "s" else "l",main=NULL,
			xlab="time",ylab="lineages",log=if(log) "y" else "")
		lines(x[,1],x[,2],lty="dashed",type="s")
		lines(x[,1],x[,4],lty="dashed",type="s")
	}
}
##
ltt95_mod<-function (trees, alpha = 0.05, log = TRUE, method = "times", mode = "mean", ...) 
{
    if (!inherits(trees, "multiPhylo")) 
        stop("trees should be an object of class \"multiPhylo\".")
    method <- method[1]
    mode <- mode[1]
	res <- 100
    X <- ltt(trees, plot = FALSE, gamma = FALSE)
        N <- length(X)
        tt <- sapply(X, function(x) max(x$times))
        zz <- max(tt) - tt
        for (i in 1:N) X[[i]]$times <- X[[i]]$times + zz[i]
        n <- sapply(X, function(x) max(x$ltt))
        if (all(n == max(n))) 
            n <- max(n)
        else stop("for method=\"times\" all trees must contain the same numer of lineages")
        LL <- sapply(X, function(x) x$times[1:length(x$times)])
        ii <- floor(alpha/2 * N)
        jj <- ceiling((1 - alpha/2) * N)
        low <- apply(LL, 1, function(x) sort(x)[ii])
        high <- apply(LL, 1, function(x) sort(x)[jj])
        ll <- rowMeans(LL)
        obj <- cbind(c(1:n, n), low, ll, high)
        colnames(obj) <- c("lineages", "low(time)", "time", "high(time)")
        rownames(obj) <- NULL
    attr(obj, "class") <- "ltt95"
    attr(obj, "alpha") <- alpha
    attr(obj, "method") <- method
    attr(obj, "mode") <- mode
    attr(obj, "log") <- log
    invisible(obj)
}

##
dd<-ltt95_mod(mamFBD_100_marsupials, xaxis="negative", alpha=0.05, log=TRUE, method="times", mode="mean")
dd[1:5,]

#      lineages low(time)     time high(time)
# [1,]        1  2.892210 11.44651   19.08623
# [2,]        2  2.892222 11.44652   19.08623
# [3,]        3  6.959139 15.07496   22.20218
# 
# >>The area bounded by low and high is what you want.

x<-dd
n<-max(x[,1])
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
x.lim<-c(min(x[,2:4]),max(x[,2:4]))

plot(x[,3],x[,1],lwd=2,xlim=x.lim,type="l",main=NULL,log="y",font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")

polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col = "grey50", border = NA) #awesome.

lines(x[,3],x[,1],lty=1, lwd=2,type="l")

# Nice!
###
# Now, scale this up to the order-level comparisons, first with all placental orders
###
# LTTs with shaded 95% CIs... - Placentals, ORDER-level
##
ordNames<-names(which(sort(table(cladesDR$ord),decreasing=TRUE) > 2))

ordNames_P<-c(ordNames[1:6],ordNames[9],ordNames[11:12],ordNames[14:18],ordNames[20:21],ordNames[23])

ordTipNames_P<-vector("list",length(ordNames_P))
for (i in 1:length(ordNames_P)){
	ordTipNames_P[[i]]<-cladesDR[which(cladesDR$ord==ordNames_P[i]),"tiplabel"]
}

# get the 100 tree samples -- Placental ORDERS
##
ltt95_allOrds<-vector("list", length(ordNames_P))
mamFBD_100_allOrds<-vector("list", length(ordNames_P))

for (j in 1:length(ordNames_P)){
	mamFBD_100_i<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_ORDS_",ordNames_P[j],".trees",sep=""))
	mamFBD_100_allOrds[[j]]<-mamFBD_100_i
	ltt95_allOrds[[j]] <- ltt95_mod(mamFBD_100_i, alpha=0.05, log=TRUE, method="times", mode="mean")
}

# get the MCCs -- Placental ORDERS
mamFBD_MCC <- drop.tip(read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target.tre"),"_Anolis_carolinensis")

mamFBD_MCC_ord_P<-vector("list",length(ordNames_P))
for (i in 1:length(mamFBD_MCC_ord_P)){
	toDrop<-setdiff(mamFBD_MCC$tip.label,ordTipNames_P[[i]])
	mamFBD_MCC_ord_P[[i]]<-drop.tip(mamFBD_MCC,toDrop)
}

cols<-palette(rainbow(length(ordTipNames_P),alpha=0.3))

# PLOT the polygons for each order - TOGETHER

pdf(file="_test_LTT_95polygon_Placentals_eachOrd_wNames.pdf", width=5, height=5, onefile=TRUE)

xLim<-c(-115,0)

x<-ltt95_allOrds[[1]]
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[1], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
text(x=0,y=nrow(ltt95_allOrds[[1]][,])-1,labels=ordNames_P[1],pos=2,cex=0.3,offset=0)

for (i in 2:length(ltt95_allOrds)){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=ordNames_P[i],pos=2, cex=0.3,offset=0)
}
dev.off()

######
# PLOT the polygons for each order - TOGETHER -- With PLACENTALS LTT...
# prep the placentals LTT...
mamFBD_100_placentals<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_PLACENTALS_nexus.trees")
ltt95_placentals <- ltt95_mod(mamFBD_100_placentals, alpha=0.05, log=TRUE, method="times", mode="mean")

# plot
jpeg(file="_test_LTT_95polygon_PlacentalsLTTpolygon_eachOrd_wNames_quadMarsupRodBat.jpg", width=10, height=10,units="in", res=450, quality=100)
xLim<-c(-110,13)

cols<-palette(rainbow(length(ordTipNames_P),alpha=0.3))
cols<-palette(rainbow(length(ordTipNames_P),alpha=0.3))

#layout(matrix(c(1:2), 1, 2, byrow = TRUE))
op <- par(mfrow = c(2,2),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)

#smoothScatter(xy_all_P_appended, bandwidth=c(3,1)/3, log="y", xlim=xLim, colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=1000, font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
x<-ltt95_placentals
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="", ylab="",xaxt="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
text(x=0,y=nrow(ltt95_placentals[,])-1,labels="PLACENTALIA",cex=0.3,offset=0, font=2,adj=c(0,3.5))#pos=4)
text(x=-115,y=5000,labels="(a) Placental mammals",cex=1,font=2,pos=4)

#ltt.lines(drop.tip(mamFBD_MCC,marMono), col="black", lty = 1, lwd=1)

x<-ltt95_allOrds[[1]]
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[1], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
text(x=0,y=nrow(ltt95_allOrds[[1]][,])-1,labels=ordNames_P[1],cex=0.3,offset=0, adj=c(0,0))#pos=4)

for (i in 2:3){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=ordNames_P[i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}
for (i in 4:4){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=ordNames_P[i],cex=0.3,offset=0,adj=c(0,2))#pos=4)
}
for (i in 5:10){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=ordNames_P[i],cex=0.3,offset=0,adj=c(0,1))#pos=4)
}
for (i in 11:11){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=ordNames_P[i],cex=0.3,offset=0,adj=c(0,2))#pos=4)
}
for (i in 12:12){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=ordNames_P[i],cex=0.3,offset=0,adj=c(0,3))#pos=4)
}
for (i in 13:15){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=ordNames_P[i],cex=0.3,offset=0,adj=c(0,1))#pos=4)
}
for (i in 16:16){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=ordNames_P[i],cex=0.3,offset=0,adj=c(0,-1))#pos=4)
}
for (i in 17:17){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=ordNames_P[i],cex=0.3,offset=0,adj=c(0,1))#pos=4)
}

#dev.off()
cols<-palette(rainbow(length(ordTipNames_M),alpha=0.3))
cols<-palette(rainbow(length(ordTipNames_M),alpha=0.3))


x<-ltt95_placentals
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="", ylab="",box=TRUE,yaxt="n",,xaxt="n")
#axis(side=1)

#smoothScatter(xy_all_P_appended, bandwidth=c(3,1)/3, log="y", xlim=xLim, colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=1000, font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
x<-ltt95_marsupials
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
text(x=0,y=nrow(ltt95_marsupials[,])-1,labels="MARSUPALIA",cex=0.3,offset=0, font=2,adj=c(0,3.5))#pos=4)
text(x=-115,y=5000,labels="(b) Marsupial mammals",cex=1,font=2,pos=4)

#ltt.lines(drop.tip(mamFBD_MCC,marMono), col="black", lty = 1, lwd=1)

x<-ltt95_allOrds_M[[1]]
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[1], border = NA) #awesome.
#lines(x[,3],x[,1],lty=2,lwd=1,type="l")
text(x=0,y=nrow(ltt95_allOrds_M[[1]][,])-1,labels=ordNames_M[1],cex=0.3,offset=0, adj=c(0,0))#pos=4)

for (i in 2:5){
	x<-ltt95_allOrds_M[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds_M[[i]][,])-1,labels=ordNames_M[i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}

cols<-palette(rainbow(length(patchNames_Rod),alpha=0.3))
cols<-palette(rainbow(length(patchNames_Rod),alpha=0.3))


x<-ltt95_placentals
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="", ylab="")

#smoothScatter(xy_all_P_appended, bandwidth=c(3,1)/3, log="y", xlim=xLim, colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=1000, font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
x<-ltt95_rodents
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
text(x=0,y=nrow(ltt95_rodents[,])-1,labels="RODENTIA",cex=0.3,offset=0, font=2,adj=c(0,3.5))#pos=4)
text(x=-115,y=5000,labels="(c) Rodents",cex=1,font=2,pos=4)

#ltt.lines(drop.tip(mamFBD_MCC,marMono), col="black", lty = 1, lwd=1)

for (i in 10:4){
	x<-ltt95_allPCs_R[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allPCs_R[[i]][,])-1,labels=patchNames_Rod2[i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}
for (i in 3:3){
	x<-ltt95_allPCs_R[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allPCs_R[[i]][,])-1,labels=patchNames_Rod2[i],cex=0.3,offset=0,adj=c(0,-1))#pos=4)
}
for (i in 2:2){
	x<-ltt95_allPCs_R[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allPCs_R[[i]][,])-1,labels=patchNames_Rod2[i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}
for (i in 1:1){
	x<-ltt95_allPCs_R[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allPCs_R[[i]][,])-1,labels=patchNames_Rod2[i],cex=0.3,offset=0,adj=c(0,-1))#pos=4)
}
#dev.off()

cols<-palette(rainbow(length(famTipNames_Bat),alpha=0.3))
cols<-palette(rainbow(length(famTipNames_Bat),alpha=0.3))

x<-ltt95_placentals
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="", ylab="",yaxt="n")

#smoothScatter(xy_all_P_appended, bandwidth=c(3,1)/3, log="y", xlim=xLim, colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=1000, font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
x<-ltt95_bats
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
text(x=0,y=nrow(ltt95_bats[,])-1,labels="CHIROPTERA",cex=0.3,offset=0, font=2,adj=c(0,3.5))#pos=4)
text(x=-115,y=5000,labels="(d) Bats",cex=1,font=2,pos=4)

#ltt.lines(drop.tip(mamFBD_MCC,marMono), col="black", lty = 1, lwd=1)

for (i in 1:2){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}
for (i in 3:3){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,2))#pos=4)
}
for (i in 4:5){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}
for (i in 6:6){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,2))#pos=4)
}
for (i in 7:11){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}
for (i in 12:12){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,1.5))#pos=4)
}
for (i in 13:13){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,3))#pos=4)
}

title(xlab = "Time before present (Ma)",
      ylab = "(log) Number of lineages",
      outer = TRUE, line = 3,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2)
par(op)

dev.off()


####
# Same now, with MARSUPIALS
mamFBD_100_marsupials<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_nexus_Marsupials.trees")
ltt95_marsupials <- ltt95_mod(mamFBD_100_marsupials, alpha=0.05, log=TRUE, method="times", mode="mean")

ordNames<-names(which(sort(table(cladesDR$ord),decreasing=TRUE) > 2))

ordNames_M<-c(ordNames[7:8],ordNames[10],ordNames[13],ordNames[19])

ordTipNames_M<-vector("list",length(ordNames_M))
for (i in 1:length(ordNames_M)){
	ordTipNames_M[[i]]<-cladesDR[which(cladesDR$ord==ordNames_M[i]),"tiplabel"]
}

ltt95_allOrds_M<-vector("list", length(ordNames_M))
mamFBD_100_allOrds_M<-vector("list", length(ordNames_M))

for (j in 1:length(ordNames_M)){
	mamFBD_100_i<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_ORDS_",ordNames_M[j],".trees",sep=""))
	mamFBD_100_allOrds_M[[j]]<-mamFBD_100_i
	ltt95_allOrds_M[[j]] <- ltt95_mod(mamFBD_100_i, alpha=0.05, log=TRUE, method="times", mode="mean")
}

cols<-palette(rainbow(length(ordTipNames_M),alpha=0.3))

# PLOT

jpeg(file="_test_LTT_95polygon_MarsupialsLTTpolygon_eachOrd_wNames_sameY.jpg", width=5, height=5,units="in", res=450, quality=100)
xLim<-c(-110,13)

x<-ltt95_placentals
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")

#smoothScatter(xy_all_P_appended, bandwidth=c(3,1)/3, log="y", xlim=xLim, colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=1000, font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
x<-ltt95_marsupials
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
text(x=0,y=nrow(ltt95_marsupials[,])-1,labels="MARSUPALIA",cex=0.3,offset=0, font=2,adj=c(0,3.5))#pos=4)

#ltt.lines(drop.tip(mamFBD_MCC,marMono), col="black", lty = 1, lwd=1)

x<-ltt95_allOrds_M[[1]]
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[1], border = NA) #awesome.
#lines(x[,3],x[,1],lty=2,lwd=1,type="l")
text(x=0,y=nrow(ltt95_allOrds_M[[1]][,])-1,labels=ordNames_M[1],cex=0.3,offset=0, adj=c(0,0))#pos=4)

for (i in 2:5){
	x<-ltt95_allOrds_M[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds_M[[i]][,])-1,labels=ordNames_M[i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}

dev.off()

##
# Now RODENTS
mamFBD_100_rodentia<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_ORDS_RODENTIA.trees")
ltt95_rodents <- ltt95_mod(mamFBD_100_rodentia, alpha=0.05, log=TRUE, method="times", mode="mean")

patchNames<-names(which(sort(table(cladesDR$PC),decreasing=TRUE) > 2))
patchNames_Rod<-c(patchNames[1:2],patchNames[9:10],patchNames[13],patchNames[17:18],patchNames[21],patchNames[23],patchNames[25])
patchNames_Rod2<-c("MURIDAE", "CRICETIDAE", "SQUIRREL_RELATED", "GUINEAPIG_RELATED", "CASTORIMORPHA", "NESOMYIDAE", "DIPODIDAE", "SPALACIDAE", "ANOMALUROMORPHA", "CALOMYSCIDAE")

ltt95_allPCs_R<-vector("list", length(patchNames_Rod))
mamFBD_100_allPCs_R<-vector("list", length(patchNames_Rod))

for (j in 1:length(patchNames_Rod)){
	mamFBD_100_i<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_PCS_",patchNames_Rod[j],".trees",sep=""))
	mamFBD_100_allPCs_R[[j]]<-mamFBD_100_i
	ltt95_allPCs_R[[j]] <- ltt95_mod(mamFBD_100_i, alpha=0.05, log=TRUE, method="times", mode="mean")
}

cols<-palette(rainbow(length(patchNames_Rod),start=0,end=0.5,alpha=0.3))

# PLOT
jpeg(file="_test_LTT_95polygon_RodentiaLTTpolygon_eachPC_wNames_sameY_dualChiroptera.jpg", width=10, height=5,units="in", res=450, quality=100)
xLim<-c(-110,13)

cols<-palette(rainbow(length(patchNames_Rod),alpha=0.3))
cols<-palette(rainbow(length(patchNames_Rod),alpha=0.3))

#layout(matrix(c(1:2), 1, 2, byrow = TRUE))
op <- par(mfrow = c(1,2),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)

x<-ltt95_placentals
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")

#smoothScatter(xy_all_P_appended, bandwidth=c(3,1)/3, log="y", xlim=xLim, colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=1000, font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
x<-ltt95_rodents
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
text(x=0,y=nrow(ltt95_rodents[,])-1,labels="RODENTIA",cex=0.3,offset=0, font=2,adj=c(0,3.5))#pos=4)
text(x=-115,y=5000,labels="(c) Rodents",cex=1,font=2,pos=4)

#ltt.lines(drop.tip(mamFBD_MCC,marMono), col="black", lty = 1, lwd=1)

for (i in 10:4){
	x<-ltt95_allPCs_R[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allPCs_R[[i]][,])-1,labels=patchNames_Rod2[i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}
for (i in 3:3){
	x<-ltt95_allPCs_R[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allPCs_R[[i]][,])-1,labels=patchNames_Rod2[i],cex=0.3,offset=0,adj=c(0,-1))#pos=4)
}
for (i in 2:2){
	x<-ltt95_allPCs_R[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allPCs_R[[i]][,])-1,labels=patchNames_Rod2[i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}
for (i in 1:1){
	x<-ltt95_allPCs_R[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allPCs_R[[i]][,])-1,labels=patchNames_Rod2[i],cex=0.3,offset=0,adj=c(0,-1))#pos=4)
}
#dev.off()

cols<-palette(rainbow(length(famTipNames_Bat),alpha=0.3))
cols<-palette(rainbow(length(famTipNames_Bat),alpha=0.3))

x<-ltt95_placentals
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="",yaxt="n")

#smoothScatter(xy_all_P_appended, bandwidth=c(3,1)/3, log="y", xlim=xLim, colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=1000, font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
x<-ltt95_bats
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
text(x=0,y=nrow(ltt95_bats[,])-1,labels="CHIROPTERA",cex=0.3,offset=0, font=2,adj=c(0,3.5))#pos=4)
text(x=-115,y=5000,labels="(d) Bats",cex=1,font=2,pos=4)

#ltt.lines(drop.tip(mamFBD_MCC,marMono), col="black", lty = 1, lwd=1)

for (i in 1:2){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}
for (i in 3:3){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,2))#pos=4)
}
for (i in 4:5){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}
for (i in 6:6){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,2))#pos=4)
}
for (i in 7:11){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}
for (i in 12:12){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,1.5))#pos=4)
}
for (i in 13:13){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,3))#pos=4)
}

title(xlab = "Time before present (Ma)",
      ylab = "(log) Number of lineages",
      outer = TRUE, line = 3,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2)
par(op)

dev.off()




###
# RODENTS with FAMILIES
###
rodents<-cladesDR[ which(cladesDR$ord=="RODENTIA"),]
famNames_Rod<-names(which(sort(table(rodents$fam),decreasing=TRUE) >2))

famTipNames_Rod<-vector("list",length(famNames_Rod))
for (i in 1:length(famTipNames_Rod)){
	famTipNames_Rod[[i]]<-cladesDR[which(cladesDR$fam==famNames_Rod[i]),"tiplabel"]
}

for (j in 1:length(famNames_Rod)){

mamFBD_100_rodFams_i<-vector("list",length(mamFBD_100))
for (i in 1:length(mamFBD_100)){
	toDrop<-setdiff(mamFBD_100[[1]]$tip.label,famTipNames_Rod[[j]])
	mamFBD_100_rodFams_i[[i]]<-drop.tip(mamFBD_100[[i]],toDrop)
}
write.nexus(mamFBD_100_rodFams_i,file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_rodFAMS_",famNames_Rod[j],".trees",sep=""))
}

ltt95_allFAMs_R<-vector("list", length(famNames_Rod))
mamFBD_100_allFAMs_R<-vector("list", length(famNames_Rod))

for (j in 1:length(famNames_Rod)){
	mamFBD_100_i<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_rodFAMS_",famNames_Rod[j],".trees",sep=""))
	mamFBD_100_allFAMs_R[[j]]<-mamFBD_100_i
	ltt95_allFAMs_R[[j]] <- ltt95_mod(mamFBD_100_i, alpha=0.05, log=TRUE, method="times", mode="mean")
}

cols<-palette(rainbow(length(famNames_Rod),alpha=0.3))

cols<-palette(rainbow(length(famNames_Rod),start=0,end=0.5,alpha=0.3))


# PLOT
jpeg(file="_test_LTT_95polygon_RodentiaLTTpolygon_eachFAM_wNames_sameY_rainbow0p5.jpg", width=5, height=5,units="in", res=450, quality=100)
xLim<-c(-110,13)

x<-ltt95_placentals
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")

#smoothScatter(xy_all_P_appended, bandwidth=c(3,1)/3, log="y", xlim=xLim, colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=1000, font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
x<-ltt95_rodents
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
text(x=0,y=nrow(ltt95_rodents[,])-1,labels="RODENTIA",cex=0.3,offset=0, font=2,adj=c(0,3.5))#pos=4)

#ltt.lines(drop.tip(mamFBD_MCC,marMono), col="black", lty = 1, lwd=1)

for (i in 1:23){
	x<-ltt95_allFAMs_R[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_R[[i]][,])-1,labels=famNames_Rod[i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}

dev.off()

###
# BATS - with PCs
mamFBD_100_bats<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_ORDS_CHIROPTERA.trees")
ltt95_bats <- ltt95_mod(mamFBD_100_bats, alpha=0.05, log=TRUE, method="times", mode="mean")

patchNames<-names(which(sort(table(cladesDR$PC),decreasing=TRUE) > 2))
patchNames_Bat<-c(patchNames[3],patchNames[6],patchNames[12],patchNames[16])
patchNames_Bat2<-c("VESPERTILIONID_RELATED", "YINPTEROCHIROPTERA", "PHYLLOSTOMID_RELATED", "EMBALLONURID_RELATED")

ltt95_allPCs_B<-vector("list", length(patchNames_Bat))
mamFBD_100_allPCs_B<-vector("list", length(patchNames_Bat))

for (j in 1:length(patchNames_Bat2)){
	mamFBD_100_i<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_PCS_",patchNames_Bat[j],".trees",sep=""))
	mamFBD_100_allPCs_B[[j]]<-mamFBD_100_i
	ltt95_allPCs_B[[j]] <- ltt95_mod(mamFBD_100_i, alpha=0.05, log=TRUE, method="times", mode="mean")
}

cols<-palette(rainbow(length(patchNames_Bat),start=0,end=0.5,alpha=0.3))

jpeg(file="_test_LTT_95polygon_ChiropteraLTTpolygon_eachPC_wNames_sameY.jpg", width=5, height=5,units="in", res=450, quality=100)
xLim<-c(-110,13)

cols<-palette(rainbow(length(patchNames_Bat),alpha=0.3))
cols<-palette(rainbow(length(patchNames_Bat),alpha=0.3))

x<-ltt95_placentals
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")

#smoothScatter(xy_all_P_appended, bandwidth=c(3,1)/3, log="y", xlim=xLim, colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=1000, font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
x<-ltt95_bats
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
text(x=0,y=nrow(ltt95_bats[,])-1,labels="CHIROPTERA",cex=0.3,offset=0, font=2,adj=c(0,3.5))#pos=4)

#ltt.lines(drop.tip(mamFBD_MCC,marMono), col="black", lty = 1, lwd=1)

for (i in 1:4){
	x<-ltt95_allPCs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allPCs_B[[i]][,])-1,labels=patchNames_Bat2[i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}

dev.off()


###
# BATS, fam-level
bats<-cladesDR[ which(cladesDR$ord=="CHIROPTERA"),]
famNames_Bat<-names(which(sort(table(bats$fam),decreasing=TRUE) >2))

famTipNames_Bat<-vector("list",length(famNames_Bat))
for (i in 1:length(famTipNames_Bat)){
	famTipNames_Bat[[i]]<-cladesDR[which(cladesDR$fam==famNames_Bat[i]),"tiplabel"]
}

for (j in 1:length(famNames_Bat)){

mamFBD_100_batFams_i<-vector("list",length(mamFBD_100))
for (i in 1:length(mamFBD_100)){
	toDrop<-setdiff(mamFBD_100[[1]]$tip.label,famTipNames_Bat[[j]])
	mamFBD_100_batFams_i[[i]]<-drop.tip(mamFBD_100[[i]],toDrop)
}
write.nexus(mamFBD_100_batFams_i,file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_batFAMS_",famNames_Bat[j],".trees",sep=""))
}

ltt95_allFAMs_B<-vector("list", length(famNames_Bat))
mamFBD_100_allFAMs_B<-vector("list", length(famNames_Bat))

for (j in 1:length(famNames_Bat)){
	mamFBD_100_i<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_batFAMS_",famNames_Bat[j],".trees",sep=""))
	mamFBD_100_allFAMs_B[[j]]<-mamFBD_100_i
	ltt95_allFAMs_B[[j]] <- ltt95_mod(mamFBD_100_i, alpha=0.05, log=TRUE, method="times", mode="mean")
}


jpeg(file="_test_LTT_95polygon_ChiropteraLTTpolygon_eachFAM_wNames_sameY.jpg", width=5, height=5,units="in", res=450, quality=100)
xLim<-c(-110,13)

cols<-palette(rainbow(length(famTipNames_Bat),alpha=0.3))
cols<-palette(rainbow(length(famTipNames_Bat),alpha=0.3))

x<-ltt95_placentals
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")

#smoothScatter(xy_all_P_appended, bandwidth=c(3,1)/3, log="y", xlim=xLim, colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=1000, font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
x<-ltt95_bats
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
text(x=0,y=nrow(ltt95_bats[,])-1,labels="CHIROPTERA",cex=0.3,offset=0, font=2,adj=c(0,3.5))#pos=4)

#ltt.lines(drop.tip(mamFBD_MCC,marMono), col="black", lty = 1, lwd=1)

for (i in 1:2){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}
for (i in 3:3){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,2))#pos=4)
}
for (i in 4:5){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}
for (i in 6:6){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,2))#pos=4)
}
for (i in 7:11){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}
for (i in 12:12){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,1.5))#pos=4)
}
for (i in 13:13){
	x<-ltt95_allFAMs_B[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allFAMs_B[[i]][,])-1,labels=famNames_Bat[i],cex=0.3,offset=0,adj=c(0,3))#pos=4)
}
dev.off()



######
# SMOOTH SCATTER
# combine all 100 trees branching times into a single data frame
xy_all<-vector("list", length(mamFBD_100))
for (i in 1:length(mamFBD_100)){
	xy_all[[i]] <- ltt.plot.coords(mamFBD_100[[i]])
}
xy_all_appended<-do.call(rbind, xy_all)

mamFBD_MCC <- drop.tip(read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target.tre"),"_Anolis_carolinensis")

# jpeg(file="_test_LTT_smoothScatter_MCC_3kBin.jpg", width=8, height=8, units="in", res=450, quality=100)
# smoothScatter(xy_all_appended, bandwidth=c(3,1)/3, log="y", colramp = colorRampPalette(c("white", "blue")), nrpoints=0, nbin=3000, font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
# ltt.lines(mamFBD_MCC, col="red", lty = 1, lwd=2)
# title(main="MamPhy, BDvr, 5911 spp, FBD sample100")
# dev.off()
# 
# jpeg(file="_test_LTT_smoothScatter_MCC_4kBin_9blue.jpg", width=8, height=8, units="in", res=450, quality=100)
# smoothScatter(xy_all_appended, bandwidth=c(3,1)/3, log="y", colramp = colorRampPalette(c("white", blues9)), nrpoints=0, nbin=4000, font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
# ltt.lines(mamFBD_MCC, col="red", lty = 1, lwd=2.5)
# title(main="MamPhy, BDvr, 5911 spp, FBD sample100")
# dev.off()
# 

# Drop the MCC tree to just each ORDER, for plotting on top
##
cladesDR<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_DR-SUMMARY-expanded_FBD.txt")
head(cladesDR)
colnames(cladesDR)<-c("tiplabel","gen","fam","ord","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")

ordNames<-names(which(sort(table(cladesDR$ord),decreasing=TRUE) > 2))

ordTipNames<-vector("list",length(ordNames))
for (i in 1:length(ordNames)){
	ordTipNames[[i]]<-cladesDR[which(cladesDR$ord==ordNames[i]),"tiplabel"]
}

mamFBD_MCC_ord<-vector("list",length(ordNames))
for (i in 1:length(mamFBD_MCC_ord)){
	toDrop<-setdiff(mamFBD_MCC$tip.label,ordTipNames[[i]])
	mamFBD_MCC_ord[[i]]<-drop.tip(mamFBD_MCC,toDrop)
}
cols<-palette(rainbow(length(ordTipNames)))

## The MCC of each ORDER in the same LTT plot as the MamPhy...

jpeg(file="_test_LTT_smoothScatter_MCC_5kBin_4blue_withORDs.jpg", width=8, height=8, units="in", res=450, quality=100)

smoothScatter(xy_all_appended, bandwidth=c(3,1)/3, log="y", colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=5000, font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
ltt.lines(mamFBD_MCC, col="black", lty = 2, lwd=2.5)

for(i in 1:length(mamFBD_MCC_ord)){
	ltt.lines(mamFBD_MCC_ord[[i]], col=cols[i], lty = 1, lwd=2)
}
title(main="MamPhy, BDvr, 5911 spp, FBD sample100, MCC of 23 ORDERS")

dev.off()

# focus on PLACENTALS
##
marsupials<-cladesDR[which(cladesDR$PC=="PC1_Marsupials"),"tiplabel"]
monotremes<-cladesDR[which(cladesDR$PC=="PC23_Monotremata"),"tiplabel"]
marMono<-c(as.vector(marsupials),as.vector(monotremes))

mamFBD_100_placentals<-vector("list",length(mamFBD_100))
xy_all_P<-vector("list", length(mamFBD_100_placentals))
for (i in 1:length(mamFBD_100_placentals)){
	mamFBD_100_placentals[[i]]<-drop.tip(mamFBD_100[[i]],marMono)
	xy_all_P[[i]] <- ltt.plot.coords(mamFBD_100_placentals[[i]])
}
xy_all_P_appended<-do.call(rbind, xy_all_P)

#######
# MCCs - Placentals, ORDER-level
##
# removing the 6 marMono orders...
# DIPROTODONTIA
# DIDELPHIMORPHIA
# DASYUROMORPHIA
# PERAMELEMORPHIA
# PAUCITUBERCULATA
# MONOTREMATA

ordNames<-names(which(sort(table(cladesDR$ord),decreasing=TRUE) > 2))

ordNames_P<-c(ordNames[1:6],ordNames[9],ordNames[11:12],ordNames[14:18],ordNames[20:21],ordNames[23])

ordTipNames_P<-vector("list",length(ordNames_P))
for (i in 1:length(ordNames_P)){
	ordTipNames_P[[i]]<-cladesDR[which(cladesDR$ord==ordNames_P[i]),"tiplabel"]
}

mamFBD_MCC <- drop.tip(read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target.tre"),"_Anolis_carolinensis")

mamFBD_MCC_ord_P<-vector("list",length(ordNames_P))
for (i in 1:length(mamFBD_MCC_ord_P)){
	toDrop<-setdiff(mamFBD_MCC$tip.label,ordTipNames_P[[i]])
	mamFBD_MCC_ord_P[[i]]<-drop.tip(mamFBD_MCC,toDrop)
}
cols<-palette(rainbow(length(ordTipNames_P)))

##
# Plotting, placentals with ORDS
jpeg(file="_test_LTT_smoothScatter_MCC-Placentals_1kBin_4blue_withORDs.jpg", width=8, height=8, units="in", res=450, quality=100)

smoothScatter(xy_all_P_appended, bandwidth=c(3,1)/3, log="y", colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=1000, font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
ltt.lines(drop.tip(mamFBD_MCC,marMono), col="black", lty = 2, lwd=2.5)

for(i in 1:length(mamFBD_MCC_ord_P)){
	ltt.lines(mamFBD_MCC_ord_P[[i]], col=cols[i], lty = 1, lwd=2)
}
title(main="MamPhy, BDvr, 5911 spp, FBD sample100, Placentals-- MCC of 17 ORDERS")
dev.off()


jpeg(file="_test_LTT_smoothScatter_MCC-Placentals_5kBin_4blue_withORDs.jpg", width=8, height=8, units="in", res=450, quality=100)

smoothScatter(xy_all_P_appended, bandwidth=c(3,1)/3, log="y", colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=5000, font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
ltt.lines(drop.tip(mamFBD_MCC,marMono), col="black", lty = 2, lwd=2.5)

for(i in 1:length(mamFBD_MCC_ord_P)){
	ltt.lines(mamFBD_MCC_ord_P[[i]], col=cols[i], lty = 1, lwd=2)

}
title(main="MamPhy, BDvr, 5911 spp, FBD sample100, Placentals-- MCC of 17 ORDERS")

dev.off()

###
# focus on MARSUPIALS
##
marsupials<-cladesDR[which(cladesDR$PC=="PC1_Marsupials"),"tiplabel"]
monotremes<-cladesDR[which(cladesDR$PC=="PC23_Monotremata"),"tiplabel"]
marMono<-c(as.vector(marsupials),as.vector(monotremes))
placentals<-setdiff(mamFBD_100[[1]]$tip.label,marMono)
placMono<-c(as.vector(placentals),as.vector(monotremes))

mamFBD_100_marsupials<-vector("list",length(mamFBD_100))
xy_all_M<-vector("list", length(mamFBD_100_marsupials))
for (i in 1:length(mamFBD_100_marsupials)){
	mamFBD_100_marsupials[[i]]<-drop.tip(mamFBD_100[[i]],placMono)
	xy_all_M[[i]] <- ltt.plot.coords(mamFBD_100_marsupials[[i]])
}
xy_all_M_appended<-do.call(rbind, xy_all_M)

#######
# MCCs - Marsupials, ORDER-level
##
# removing the 17 + 1 placental Mono orders...
# DIPROTODONTIA
# DIDELPHIMORPHIA
# DASYUROMORPHIA
# PERAMELEMORPHIA
# PAUCITUBERCULATA
# MONOTREMATA

ordNames<-names(which(sort(table(cladesDR$ord),decreasing=TRUE) > 2))

ordNames_M<-c(ordNames[7:8],ordNames[10],ordNames[13],ordNames[19])

ordTipNames_M<-vector("list",length(ordNames_M))
for (i in 1:length(ordNames_M)){
	ordTipNames_M[[i]]<-cladesDR[which(cladesDR$ord==ordNames_M[i]),"tiplabel"]
}

mamFBD_MCC <- drop.tip(read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target.tre"),"_Anolis_carolinensis")

mamFBD_MCC_ord_M<-vector("list",length(ordNames_M))
for (i in 1:length(mamFBD_MCC_ord_M)){
	toDrop<-setdiff(mamFBD_MCC$tip.label,ordTipNames_M[[i]])
	mamFBD_MCC_ord_M[[i]]<-drop.tip(mamFBD_MCC,toDrop)
}
cols<-palette(rainbow(length(ordTipNames_M)))

##
# Plotting, MARSUPIALS with ORDS
jpeg(file="_test_LTT_smoothScatter_MCC-Marsupials_1kBin_4blue_withORDs.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamFBD_MCC,marMono), col="white", lty = 2, lwd=2.5, log="y", font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
smoothScatter(xy_all_M_appended, bandwidth=c(3,1)/3, add=TRUE, log="y", colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=1000) #font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
ltt.lines(drop.tip(mamFBD_MCC,placMono), col="black", lty = 2, lwd=2.5)

for(i in 1:length(mamFBD_MCC_ord_M)){
	ltt.lines(mamFBD_MCC_ord_M[[i]], col=cols[i], lty = 1, lwd=2)
}
title(main="MamPhy, BDvr, 5911 spp, FBD sample100, Marsupials-- MCC of 5 ORDERS")
dev.off()


jpeg(file="_test_LTT_smoothScatter_MCC-Marsupials_5kBin_4blue_withORDs.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamFBD_MCC,marMono), col="white", lty = 2, lwd=2.5, log="y", font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
smoothScatter(xy_all_M_appended, bandwidth=c(3,1)/3, add=TRUE, log="y", colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=5000) #font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
ltt.lines(drop.tip(mamFBD_MCC,placMono), col="black", lty = 2, lwd=2.5)

for(i in 1:length(mamFBD_MCC_ord_M)){
	ltt.lines(mamFBD_MCC_ord_M[[i]], col=cols[i], lty = 1, lwd=2)
}
title(main="MamPhy, BDvr, 5911 spp, FBD sample100, Marsupials-- MCC of 5 ORDERS")
dev.off()

######
###
# focus on RODENTIA
##
# Want to plot the Rodent PATCHES separately...
rodentNames<-cladesDR[which(cladesDR$ord=="RODENTIA"),"tiplabel"]
nonRodents<-setdiff(mamFBD_100[[1]]$tip.label,rodentNames)

mamFBD_100_rodents<-vector("list",length(mamFBD_100))
xy_all_R<-vector("list", length(mamFBD_100_rodents))
for (i in 1:length(mamFBD_100_rodents)){
	mamFBD_100_rodents[[i]]<-drop.tip(mamFBD_100[[i]],nonRodents)
	xy_all_R[[i]] <- ltt.plot.coords(mamFBD_100_rodents[[i]])
}
xy_all_R_appended<-do.call(rbind, xy_all_R)

#######
# MCCs - Rodentia, PC-level
##

patchNames<-names(which(sort(table(cladesDR$PC),decreasing=TRUE) > 2))

patchNames_Rod<-c(patchNames[1:2],patchNames[9:10],patchNames[13],patchNames[17:18],patchNames[21],patchNames[23],patchNames[25])

patchTipNames_Rod<-vector("list",length(patchNames_Rod))
for (i in 1:length(patchNames_Rod)){
	patchTipNames_Rod[[i]]<-cladesDR[which(cladesDR$PC==patchNames_Rod[i]),"tiplabel"]
}

mamFBD_MCC <- drop.tip(read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target.tre"),"_Anolis_carolinensis")

mamFBD_MCC_PC_Rod<-vector("list",length(patchTipNames_Rod))
for (i in 1:length(mamFBD_MCC_PC_Rod)){
	toDrop<-setdiff(mamFBD_MCC$tip.label,patchTipNames_Rod[[i]])
	mamFBD_MCC_PC_Rod[[i]]<-drop.tip(mamFBD_MCC,toDrop)
}
cols<-palette(rainbow(length(patchTipNames_Rod)))

###
# MCCs - Rodentia, fam-level
rodents<-cladesDR[ which(cladesDR$ord=="RODENTIA"),]

famNames_Rod<-names(which(sort(table(rodents$fam),decreasing=TRUE) >2))

famTipNames_Rod<-vector("list",length(famNames_Rod))
for (i in 1:length(famTipNames_Rod)){
	famTipNames_Rod[[i]]<-cladesDR[which(cladesDR$fam==famNames_Rod[i]),"tiplabel"]
}

mamFBD_MCC <- drop.tip(read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target.tre"),"_Anolis_carolinensis")

mamFBD_MCC_fam_Rod<-vector("list",length(famTipNames_Rod))
for (i in 1:length(mamFBD_MCC_fam_Rod)){
	toDrop<-setdiff(mamFBD_MCC$tip.label,famTipNames_Rod[[i]])
	mamFBD_MCC_fam_Rod[[i]]<-drop.tip(mamFBD_MCC,toDrop)
}
cols<-palette(rainbow(length(famTipNames_Rod)))

##
# Plotting, RODENTIA with rodent PCs
jpeg(file="_test_LTT_smoothScatter_MCC-Rodentia_1kBin_4blue_withPCs.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamFBD_MCC,marMono), col="white", lty = 2, lwd=2.5, log="y", font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
smoothScatter(xy_all_R_appended, bandwidth=c(3,1)/3, add=TRUE, log="y", colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=1000) #font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
ltt.lines(drop.tip(mamFBD_MCC,nonRodents), col="black", lty = 2, lwd=2.5)

for(i in 1:length(mamFBD_MCC_PC_Rod)){
	ltt.lines(mamFBD_MCC_PC_Rod[[i]], col=cols[i], lty = 1, lwd=2)
}
title(main="MamPhy, BDvr, 5911 spp, FBD sample100, Rodentia-- MCC of 10 PCs")
dev.off()


jpeg(file="_test_LTT_smoothScatter_MCC-Rodentia_5kBin_4blue_withORDs.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamFBD_MCC,marMono), col="white", lty = 2, lwd=2.5, log="y", font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
smoothScatter(xy_all_R_appended, bandwidth=c(3,1)/3, add=TRUE, log="y", colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=5000) #font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
ltt.lines(drop.tip(mamFBD_MCC,nonRodents), col="black", lty = 2, lwd=2.5)

for(i in 1:length(mamFBD_MCC_PC_Rod)){
	ltt.lines(mamFBD_MCC_PC_Rod[[i]], col=cols[i], lty = 1, lwd=2)
}
title(main="MamPhy, BDvr, 5911 spp, FBD sample100, Rodentia-- MCC of 10 PCs")
dev.off()

# Plotting, RODENTIA with rodent FAMILIES
jpeg(file="_test_LTT_smoothScatter_MCC-Rodentia_1kBin_4blue_withFAMs.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamFBD_MCC,marMono), col="white", lty = 2, lwd=2.5, log="y", font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
smoothScatter(xy_all_R_appended, bandwidth=c(3,1)/3, add=TRUE, log="y", colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=1000) #font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
ltt.lines(drop.tip(mamFBD_MCC,nonRodents), col="black", lty = 2, lwd=2.5)

for(i in 1:length(mamFBD_MCC_fam_Rod)){
	ltt.lines(mamFBD_MCC_fam_Rod[[i]], col=cols[i], lty = 1, lwd=2)
}
title(main="MamPhy, BDvr, 5911 spp, FBD sample100, Rodentia-- MCC of 23 FAMILIES")
dev.off()


jpeg(file="_test_LTT_smoothScatter_MCC-Rodentia_5kBin_4blue_withFAMs.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamFBD_MCC,marMono), col="white", lty = 2, lwd=2.5, log="y", font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
smoothScatter(xy_all_R_appended, bandwidth=c(3,1)/3, add=TRUE, log="y", colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=5000) #font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
ltt.lines(drop.tip(mamFBD_MCC,nonRodents), col="black", lty = 2, lwd=2.5)

for(i in 1:length(mamFBD_MCC_fam_Rod)){
	ltt.lines(mamFBD_MCC_fam_Rod[[i]], col=cols[i], lty = 1, lwd=2)
}
title(main="MamPhy, BDvr, 5911 spp, FBD sample100, Rodentia-- MCC of 23 FAMILIES")
dev.off()


######
###
# focus on BATS
##
# Want to plot the bat PATCHES separately...
batNames<-cladesDR[which(cladesDR$ord=="CHIROPTERA"),"tiplabel"]
nonBats<-setdiff(mamFBD_100[[1]]$tip.label,batNames)

mamFBD_100_bats<-vector("list",length(mamFBD_100))
xy_all_B<-vector("list", length(mamFBD_100_bats))
for (i in 1:length(mamFBD_100_bats)){
	mamFBD_100_bats[[i]]<-drop.tip(mamFBD_100[[i]],nonBats)
	xy_all_B[[i]] <- ltt.plot.coords(mamFBD_100_bats[[i]])
}
xy_all_B_appended<-do.call(rbind, xy_all_B)

#######
# MCCs - BATS, PC-level
##

patchNames<-names(which(sort(table(cladesDR$PC),decreasing=TRUE) > 2))

patchNames_Bat<-c(patchNames[3],patchNames[6],patchNames[12],patchNames[16])

patchTipNames_Bat<-vector("list",length(patchNames_Bat))
for (i in 1:length(patchTipNames_Bat)){
	patchTipNames_Bat[[i]]<-cladesDR[which(cladesDR$PC==patchNames_Bat[i]),"tiplabel"]
}

mamFBD_MCC <- drop.tip(read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target.tre"),"_Anolis_carolinensis")

mamFBD_MCC_PC_Bat<-vector("list",length(patchTipNames_Bat))
for (i in 1:length(mamFBD_MCC_PC_Bat)){
	toDrop<-setdiff(mamFBD_MCC$tip.label,patchTipNames_Bat[[i]])
	mamFBD_MCC_PC_Bat[[i]]<-drop.tip(mamFBD_MCC,toDrop)
}
cols<-palette(rainbow(length(patchTipNames_Bat)))

###
# MCCs - BATS, fam-level
bats<-cladesDR[ which(cladesDR$ord=="CHIROPTERA"),]

famNames_Bat<-names(which(sort(table(bats$fam),decreasing=TRUE) >2))

famTipNames_Bat<-vector("list",length(famNames_Bat))
for (i in 1:length(famTipNames_Bat)){
	famTipNames_Bat[[i]]<-cladesDR[which(cladesDR$fam==famNames_Bat[i]),"tiplabel"]
}

mamFBD_MCC <- drop.tip(read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_MCC_target.tre"),"_Anolis_carolinensis")

mamFBD_MCC_fam_Bat<-vector("list",length(famTipNames_Bat))
for (i in 1:length(mamFBD_MCC_fam_Bat)){
	toDrop<-setdiff(mamFBD_MCC$tip.label,famTipNames_Bat[[i]])
	mamFBD_MCC_fam_Bat[[i]]<-drop.tip(mamFBD_MCC,toDrop)
}
cols<-palette(rainbow(length(famTipNames_Bat)))



##
# Plotting, BATS with bat PCs
jpeg(file="_test_LTT_smoothScatter_MCC-Chiroptera_1kBin_4blue_withPCs.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamFBD_MCC,marMono), col="white", lty = 2, lwd=2.5, log="y", font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
smoothScatter(xy_all_B_appended, bandwidth=c(3,1)/3, add=TRUE, log="y", colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=1000) #font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
ltt.lines(drop.tip(mamFBD_MCC,nonBats), col="black", lty = 2, lwd=2.5)

for(i in 1:length(mamFBD_MCC_PC_Bat)){
	ltt.lines(mamFBD_MCC_PC_Bat[[i]], col=cols[i], lty = 1, lwd=2)
}
title(main="MamPhy, BDvr, 5911 spp, FBD sample100, Chiroptera-- MCC of 10 PCs")
dev.off()


jpeg(file="_test_LTT_smoothScatter_MCC-Chiroptera_5kBin_4blue_withPCs.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamFBD_MCC,marMono), col="white", lty = 2, lwd=2.5, log="y", font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
smoothScatter(xy_all_B_appended, bandwidth=c(3,1)/3, add=TRUE, log="y", colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=5000) #font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
ltt.lines(drop.tip(mamFBD_MCC,nonBats), col="black", lty = 2, lwd=2.5)

for(i in 1:length(mamFBD_MCC_PC_Bat)){
	ltt.lines(mamFBD_MCC_PC_Bat[[i]], col=cols[i], lty = 1, lwd=2)
}
title(main="MamPhy, BDvr, 5911 spp, FBD sample100, Chiroptera-- MCC of 10 PCs")
dev.off()

##
# Plotting, BATS with bat FAMILIES
jpeg(file="_test_LTT_smoothScatter_MCC-Chiroptera_1kBin_4blue_withFAMs.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamFBD_MCC,marMono), col="white", lty = 2, lwd=2.5, log="y", font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
smoothScatter(xy_all_B_appended, bandwidth=c(3,1)/3, add=TRUE, log="y", colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=1000) #font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
ltt.lines(drop.tip(mamFBD_MCC,nonBats), col="black", lty = 2, lwd=2.5)

for(i in 1:length(mamFBD_MCC_fam_Bat)){
	ltt.lines(mamFBD_MCC_fam_Bat[[i]], col=cols[i], lty = 1, lwd=2)
}
title(main="MamPhy, BDvr, 5911 spp, FBD sample100, Chiroptera-- MCC of 13 FAMs")
dev.off()


jpeg(file="_test_LTT_smoothScatter_MCC-Chiroptera_5kBin_4blue_withFAMs.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(drop.tip(mamFBD_MCC,marMono), col="white", lty = 2, lwd=2.5, log="y", font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
smoothScatter(xy_all_B_appended, bandwidth=c(3,1)/3, add=TRUE, log="y", colramp = colorRampPalette(c("white", "lightblue1", "deepskyblue","royalblue1","royalblue4")), nrpoints=0, nbin=5000) #font.axis=2, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
ltt.lines(drop.tip(mamFBD_MCC,nonBats), col="black", lty = 2, lwd=2.5)

for(i in 1:length(mamFBD_MCC_fam_Bat)){
	ltt.lines(mamFBD_MCC_fam_Bat[[i]], col=cols[i], lty = 1, lwd=2)
}
title(main="MamPhy, BDvr, 5911 spp, FBD sample100, Chiroptera-- MCC of 13 FAMs")
dev.off()







# plot the 100 tree LTT with DENSITY of the lines
# with ALPHA 0.1
jpeg(file="_test_LTT_density.jpg", width=8, height=8, units="in", res=450, quality=100)

ltt.plot(mamFBD_100[[1]], log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:length(mamFBD_100)) {
	ltt.lines(mamFBD_100[[i]], col=hsv(0.65,1,1,alpha=0.1), lty = 1, lwd=1)
}
title(main="MamPhy, BDvr, 5911 spp, sample100 - FBD ")

dev.off()

# 95 CI...
dd<-ltt95(trees=mamFBD_100, alpha=0.05, log=TRUE, method="lineages")

# ways to make as density...
xy <- ltt.plot.coords(mamFBD_100[[1]])
plot(xy, log="y")
smoothScatter(xy, log="y")

ltt.lines.DENS<-function (phy, backward = TRUE, tol = 1e-06, ...) 
{
    xy <- ltt.plot.coords(phy, backward, tol)
    smoothScatter(xy, log="y", add=TRUE, nrpoints=0, colramp = colorRampPalette(c("white", blues9)), ...)
}
for (i in 2:length(mamFBD_100)) {
	ltt.lines.DENS(mamFBD_100[[i]], nbin=100, add=TRUE)#, col=hsv(0.65,1,1,alpha=0.1), lty = 1, lwd=1)
}
dev.off()

# try with the old LASER function -- plotLtt
library(laser)

Btimes
x<-branching.times(mamFBD_100[[1]])


plotLtt<-function (x) 
{
    if (!is.numeric(x)) 
        stop("object x not of class 'numeric'")
    x <- rev(sort(x))
    x <- c(-x, 0)
    n <- length(x)
    z <- max(x) - (x)
    vn <- 1:n
    plot(log(vn + 1) ~ z, xlab = "Time From Basal Divergence", 
        ylab = "Log Lineages", main = "Log-Lineages Through Time")
    print(log(vn + 1))
}





#######
# Now with those 100 random trees...
library(ape)
library(phytools)
#setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors")
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")

mamFBD_100<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample100_nexus.trees")
mamNDexp_100<-read.nexus("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample100_nexus.trees")

cladesDR<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_DR-SUMMARY-expanded_FBD.txt")
head(cladesDR)
colnames(cladesDR)<-c("tiplabel","gen","fam","ord","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")

ordNames<-names(which(table(cladesDR$ord) > 2))
famNames<-names(which(table(cladesDR$fam) > 2))
patchNames<-names(which(table(cladesDR$PC) > 2))

ordTipNames<-vector("list",length(ordNames))
for (i in 1:length(ordNames)){
	ordTipNames[[i]]<-cladesDR[which(cladesDR$ord==ordNames[i]),"tiplabel"]
}

cols<-palette(rainbow(length(ordTipNames)))

## Each ORDER in own LTT plot... (OWN axis)
pdf(file="LTT_MamPhy_BDvr_FBD_27ords_EACH_sample100.pdf", onefile=TRUE, width=8, height=8) #units="in", res=450, quality=100)

for (k in 1:length(ordTipNames)){
	toDrop<-setdiff(mamFBD_100[[1]]$tip.label,ordTipNames[[k]])
	
	mamFBD_100_ord<-vector("list",length(mamFBD_100))
	for (i in 1:length(mamFBD_100_ord)){
		mamFBD_100_ord[[i]]<-drop.tip(mamFBD_100[[i]],toDrop)
	}
	ltt.plot(mamFBD_100_ord[[1]], log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col=cols[k])
	for (j in 2:length(mamFBD_100_ord)) {
		ltt.lines(mamFBD_100_ord[[j]], col=cols[k], lty = 2, lwd=1)
	}
}
dev.off()

## Each ORDER in own LTT plot... (SAME axis--BOTH)
pdf(file="LTT_MamPhy_BDvr_FBD_27ords_EACH-sameAXIS_sample100.pdf", onefile=TRUE, width=8, height=8) #units="in", res=450, quality=100)

for (k in 1:length(ordTipNames)){
	toDrop<-setdiff(mamFBD_100[[1]]$tip.label,ordTipNames[[k]])
	
	mamFBD_100_ord<-vector("list",length(mamFBD_100))
	for (i in 1:length(mamFBD_100_ord)){
		mamFBD_100_ord[[i]]<-drop.tip(mamFBD_100[[i]],toDrop)
	}
	ltt.plot(mamFBD_100_ord[[1]], log="y", xlim=c(-150,0), ylim=c(1,2000), main=ordNames[k], xlab="Time before present (Ma)", ylab="(log) Number of lineages", col=cols[k])
	for (j in 2:length(mamFBD_100_ord)) {
		ltt.lines(mamFBD_100_ord[[j]], col=cols[k], lty = 2, lwd=1)
	}
}
dev.off()


## All ORDERS in SAME LTT plot...
pdf(file="LTT_MamPhy_BDvr_FBD_27ords_ALL-in-one_sample100.pdf", onefile=TRUE, width=8, height=8) #units="in", res=450, quality=100)

ltt.plot(mamFBD_100[[1]], log="y", xlim=c(-150,0), ylim=c(1,2000), xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="white")

for (k in 1:length(ordTipNames)){
	toDrop<-setdiff(mamFBD_100[[1]]$tip.label,ordTipNames[[k]])
	
	mamFBD_100_ord<-vector("list",length(mamFBD_100))
	for (i in 1:length(mamFBD_100_ord)){
		mamFBD_100_ord[[i]]<-drop.tip(mamFBD_100[[i]],toDrop)
	}
	for (j in 1:length(mamFBD_100_ord)) {
		ltt.lines(mamFBD_100_ord[[j]], col=cols[k], lty = 2, lwd=1)
	}
}
dev.off()

###
# FAMS, same...

famTipNames<-vector("list",length(famNames))
for (i in 1:length(famNames)){
	famTipNames[[i]]<-cladesDR[which(cladesDR$fam==famNames[i]),"tiplabel"]
}

cols<-palette(rainbow(length(famTipNames)))

## Each FAM in own LTT plot... (SAME axis)
pdf(file="LTT_MamPhy_BDvr_FBD_162fams_EACH-sameAXIS_sample100.pdf", onefile=TRUE, width=8, height=8) #units="in", res=450, quality=100)

for (k in 1:length(famTipNames)){
	toDrop<-setdiff(mamFBD_100[[1]]$tip.label,famTipNames[[k]])
	
	mamFBD_100_fam<-vector("list",length(mamFBD_100))
	for (i in 1:length(mamFBD_100_fam)){
		mamFBD_100_fam[[i]]<-drop.tip(mamFBD_100[[i]],toDrop)
	}
	ltt.plot(mamFBD_100_fam[[1]], log="y", xlim=c(-100,0), ylim=c(1,800), main=famNames[k], xlab="Time before present (Ma)", ylab="(log) Number of lineages", col=cols[k])
	for (j in 2:length(mamFBD_100_fam)) {
		ltt.lines(mamFBD_100_fam[[j]], col=cols[k], lty = 2, lwd=1)
	}
}
dev.off()


## All FAMS in SAME LTT plot...
pdf(file="LTT_MamPhy_BDvr_FBD_162fams_ALL-in-one_sample100.pdf", onefile=TRUE, width=8, height=8) #units="in", res=450, quality=100)

ltt.plot(mamFBD_100[[1]], log="y", xlim=c(-100,0), ylim=c(1,800), xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="white")

for (k in 1:length(famTipNames)){
	toDrop<-setdiff(mamFBD_100[[1]]$tip.label,famTipNames[[k]])
	
	mamFBD_100_fam<-vector("list",length(mamFBD_100))
	for (i in 1:length(mamFBD_100_fam)){
		mamFBD_100_fam[[i]]<-drop.tip(mamFBD_100[[i]],toDrop)
	}
	for (j in 1:length(mamFBD_100_fam)) {
		ltt.lines(mamFBD_100_fam[[j]], col=cols[k], lty = 2, lwd=1)
	}
}
dev.off()







## JUST PRUNING the MAMPHY TO THE higher taxa....
#### NOT really interesting...
# LTTs for each...
ordSp<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_1spMostGenes_perORD.txt")
ordSpOrd<-cbind(ordSp[1],ordSp[4])
colnames(ordSpOrd)<-c("tip","ord")

famSp<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_1spMostGenes_perFAM.txt")
famSpFam<-cbind(famSp[1],famSp[3])
colnames(famSpFam)<-c("tip","fam")

patchSp<-read.table("MamPhy_5911sp_tipGenFamOrdGenesSampPC_1spMostGenes_perPC.txt")
patchSpPatch<-cbind(patchSp[1],patchSp[7])
colnames(patchSpPatch)<-c("tip","PC")

# ORDERS
toDrop<-setdiff(mamFBD_100[[1]]$tip.label,ordSpOrd$tip)

mamFBD_100_ord1<-vector("list",length(mamFBD_100))
for (i in 1:length(mamFBD_100)){
	mamFBD_100_ord1[[i]]<-drop.tip(mamFBD_100[[i]],toDrop)
}
mamFBD_100_ord2<-vector("list",length(mamNDexp_100))
for (i in 1:length(mamNDexp_100)){
	mamFBD_100_ord2[[i]]<-drop.tip(mamNDexp_100[[i]],toDrop)
}

jpeg(file="LTT_MamPhy_BDvr_FBD-vs-NDexp_27ORDS_sample100.jpg", width=8, height=8, units="in", res=450, quality=100)
ltt.plot(mamFBD_100_ord1[[1]], log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:length(mamFBD_100_ord1)) {
	ltt.lines(mamFBD_100_ord1[[i]], col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}
for (i in 1:length(mamFBD_100_ord2)) {
	ltt.lines(mamFBD_100_ord2[[i]], col=hsv(0.9708995,0.2470588,1,alpha=0.5), lty = 2, lwd=1)
}
title(main="MamPhy, BDvr, 27 ORDERS, sample100 - FBD (blue) vs ND (pink)")
dev.off()

# FAMS
toDrop<-setdiff(mamFBD_100[[1]]$tip.label,famSpFam$tip)

mamFBD_100_fam1<-vector("list",length(mamFBD_100))
for (i in 1:length(mamFBD_100)){
	mamFBD_100_fam1[[i]]<-drop.tip(mamFBD_100[[i]],toDrop)
}
mamFBD_100_fam2<-vector("list",length(mamNDexp_100))
for (i in 1:length(mamNDexp_100)){
	mamFBD_100_fam2[[i]]<-drop.tip(mamNDexp_100[[i]],toDrop)
}

jpeg(file="LTT_MamPhy_BDvr_FBD-vs-NDexp_162FAMS_sample100.jpg", width=8, height=8, units="in", res=450, quality=100)
ltt.plot(mamFBD_100_fam1[[1]], log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:length(mamFBD_100_fam1)) {
	ltt.lines(mamFBD_100_fam1[[i]], col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}
for (i in 1:length(mamFBD_100_fam2)) {
	ltt.lines(mamFBD_100_fam2[[i]], col=hsv(0.9708995,0.2470588,1,alpha=0.5), lty = 2, lwd=1)
}
title(main="MamPhy, BDvr, 162 FAMILIES, sample100 - FBD (blue) vs ND (pink)")
dev.off()

# PCs
toDrop<-setdiff(mamFBD_100[[1]]$tip.label,patchSpPatch$tip)

mamFBD_100_pc1<-vector("list",length(mamFBD_100))
for (i in 1:length(mamFBD_100)){
	mamFBD_100_pc1[[i]]<-drop.tip(mamFBD_100[[i]],toDrop)
}
mamFBD_100_pc2<-vector("list",length(mamNDexp_100))
for (i in 1:length(mamNDexp_100)){
	mamFBD_100_pc2[[i]]<-drop.tip(mamNDexp_100[[i]],toDrop)
}

jpeg(file="LTT_MamPhy_BDvr_FBD-vs-NDexp_28PCs_sample100.jpg", width=8, height=8, units="in", res=450, quality=100)
ltt.plot(mamFBD_100_pc1[[1]], log="y", xlab="Time before present (Ma)", ylab="(log) Number of lineages", col="blue")
for (i in 2:length(mamFBD_100_pc1)) {
	ltt.lines(mamFBD_100_pc1[[i]], col=hsv(0.65,1,1,alpha=0.5), lty = 2, lwd=1)
}
for (i in 1:length(mamFBD_100_pc2)) {
	ltt.lines(mamFBD_100_pc2[[i]], col=hsv(0.9708995,0.2470588,1,alpha=0.5), lty = 2, lwd=1)
}
title(main="MamPhy, BDvr, 28 PATCHCLADES, sample100 - FBD (blue) vs ND (pink)")
dev.off()


## THINK about this... your goal was not just to prune the tree to those higher clade reps... 
# and then plot LTT... all you get is major SATURATION
# REALLY-- wanted to plot all 28 patches, for example on the same LTT plot...





cladeReps<-vector()
for (i in 1:length(allCladeSets[[k]])){
	cladeSp<-allCladeSets[[k]][[i]]$tip.label
	cladeReps[i]<-cladeSp[1]
	}
toDrop<-setdiff(mamPhy$tip.label,cladeReps)
assign(paste("slice",allCladeSetNames[[k]],"Phy",sep=""),drop.tip(mamPhy,toDrop))
}

slicePhys<-list(slice5MaPhy,slice10MaPhy,slice20MaPhy,slice30MaPhy,slice40MaPhy,slice50MaPhy,slice60MaPhy)







for (k in 1:length(allCladeSets)){
cladeReps<-vector()
for (i in 1:length(allCladeSets[[k]])){
	cladeSp<-allCladeSets[[k]][[i]]$tip.label
	cladeReps[i]<-cladeSp[1]
	}
toDrop<-setdiff(mamPhy$tip.label,cladeReps)
assign(paste("slice",allCladeSetNames[[k]],"Phy",sep=""),drop.tip(mamPhy,toDrop))
}

slicePhys<-list(slice5MaPhy,slice10MaPhy,slice20MaPhy,slice30MaPhy,slice40MaPhy,slice50MaPhy,slice60MaPhy)








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
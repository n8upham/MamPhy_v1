# author: N Upham
# 7 Nov 2016
###
# packages
library(ape)
library(phytools)

# modified ltt95() function [phytools] >> faster, doesn't plot, returns upper and lower 95% CI + mean.
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

# directory
setwd("set/to/directory")

# sample of 100 trees
fullTrees_100<-read.nexus("filename.trees")

# species taxonomic assignments, arranged in a table
speciesClades<-read.table("speciesClades.txt")
head(speciesClades)
colnames(speciesClades)<-c("tiplabel","gen","fam","ord","higher")

# 95% confidence intervals around LTTs, as colored transparent polygons
###

##
# Specific to MAMMALIA, plotting out placentals, marsupials, rodents, bats
###
# specify backbone
bbone<-"FBD" # "NDexp"

# DATA PREP - get 95% CIs of the 100 trees
# generate placentals and marsupials
marsupials<-speciesClades[which(speciesClades$higher=="MARSUPIALS"),"tiplabel"]
monotremes<-speciesClades[which(speciesClades$higher=="MONOTREMES"),"tiplabel"]
placentals<-speciesClades[which(speciesClades$higher=="PLACENTALS"),"tiplabel"]
marMono<-c(as.vector(marsupials),as.vector(monotremes))
placMono<-c(as.vector(placentals),as.vector(monotremes))

fullTrees_100_placentals<-vector("list",length(fullTrees_100))
for (i in 1:length(fullTrees_100_placentals)){
	fullTrees_100_placentals[[i]]<-drop.tip(fullTrees_100[[i]],marMono)
}
write.nexus(fullTrees_100_placentals, file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_PLACENTALS.trees",sep=""))

fullTrees_100_marsupials<-vector("list",length(fullTrees_100))
for (i in 1:length(fullTrees_100_marsupials)){
	fullTrees_100_marsupials[[i]]<-drop.tip(fullTrees_100[[i]],placMono)
}
write.nexus(fullTrees_100_marsupials, file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_MARSUPIALS.trees",sep=""))


# load data
# Placentalia
fullTrees_100_placentals<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_PLACENTALS.trees",sep=""))
ltt95_placentals <- ltt95_mod(fullTrees_100_placentals, alpha=0.05, log=TRUE, method="times", mode="mean")

# placental orders (17)
placentalSp<-speciesClades[ which(speciesClades$higher=="PLACENTALS"),]
ordNames_P<-names(which(sort(table(placentalSp$ord),decreasing=TRUE) > 2))

ordTipNames_P<-vector("list",length(ordNames_P))
for (i in 1:length(ordNames_P)){
	ordTipNames_P[[i]]<-speciesClades[which(speciesClades$ord==ordNames_P[i]),"tiplabel"]
}

ltt95_allOrds<-vector("list", length(ordNames_P))
fullTrees_100_allOrds<-vector("list", length(ordNames_P))

for (j in 1:length(ordNames_P)){
	fullTrees_100_i<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_ORDS_",ordNames_P[j],".trees",sep=""))
	fullTrees_100_allOrds[[j]]<-fullTrees_100_i
	ltt95_allOrds[[j]] <- ltt95_mod(fullTrees_100_i, alpha=0.05, log=TRUE, method="times", mode="mean")
}

# Marsupalia
fullTrees_100_marsupials<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_MARSUPIALS.trees",sep=""))
ltt95_marsupials <- ltt95_mod(fullTrees_100_marsupials, alpha=0.05, log=TRUE, method="times", mode="mean")

# marsupial orders (5)
ordNames<-names(which(sort(table(speciesClades$ord),decreasing=TRUE) > 2))
ordNames_M<-c(ordNames[7:8],ordNames[10],ordNames[13],ordNames[19])

ordTipNames_M<-vector("list",length(ordNames_M))
for (i in 1:length(ordNames_M)){
	ordTipNames_M[[i]]<-speciesClades[which(speciesClades$ord==ordNames_M[i]),"tiplabel"]
}

ltt95_allOrds_M<-vector("list", length(ordNames_M))
fullTrees_100_allOrds_M<-vector("list", length(ordNames_M))

for (j in 1:length(ordNames_M)){
	fullTrees_100_i<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_ORDS_",ordNames_M[j],".trees",sep=""))
	fullTrees_100_allOrds_M[[j]]<-fullTrees_100_i
	ltt95_allOrds_M[[j]] <- ltt95_mod(fullTrees_100_i, alpha=0.05, log=TRUE, method="times", mode="mean")
}

# Rodentia
fullTrees_100_rodentia<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_ORDS_RODENTIA.trees",sep=""))
ltt95_rodents <- ltt95_mod(fullTrees_100_rodentia, alpha=0.05, log=TRUE, method="times", mode="mean")

# rodent patch-clade level (=families + multi-family clades)
patchNames<-names(which(sort(table(speciesClades$PC),decreasing=TRUE) > 2))
patchNames_Rod<-c(patchNames[1:2],patchNames[9:10],patchNames[13],patchNames[17:18],patchNames[21],patchNames[23],patchNames[25])
patchNames_Rod2<-c("MURIDAE", "CRICETIDAE", "SQUIRREL_RELATED", "GUINEAPIG_RELATED", "CASTORIMORPHA", "NESOMYIDAE", "DIPODIDAE", "SPALACIDAE", "ANOMALUROMORPHA", "CALOMYSCIDAE")

ltt95_allPCs_R<-vector("list", length(patchNames_Rod))
fullTrees_100_allPCs_R<-vector("list", length(patchNames_Rod))

for (j in 1:length(patchNames_Rod)){
	fullTrees_100_i<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_PCS_",patchNames_Rod[j],".trees",sep=""))
	fullTrees_100_allPCs_R[[j]]<-fullTrees_100_i
	ltt95_allPCs_R[[j]] <- ltt95_mod(fullTrees_100_i, alpha=0.05, log=TRUE, method="times", mode="mean")
}

# Chiroptera
fullTrees_100_bats<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_ORDS_CHIROPTERA.trees",sep=""))
ltt95_bats <- ltt95_mod(fullTrees_100_bats, alpha=0.05, log=TRUE, method="times", mode="mean")

# bats, fam-level
bats<-speciesClades[ which(speciesClades$ord=="CHIROPTERA"),]
famNames_Bat<-names(which(sort(table(bats$fam),decreasing=TRUE) >2))

famTipNames_Bat<-vector("list",length(famNames_Bat))
for (i in 1:length(famTipNames_Bat)){
	famTipNames_Bat[[i]]<-speciesClades[which(speciesClades$fam==famNames_Bat[i]),"tiplabel"]
}

# generate per-family files
for (j in 1:length(famNames_Bat)){
fullTrees_100_batFams_i<-vector("list",length(fullTrees_100))

for (i in 1:length(fullTrees_100)){
	toDrop<-setdiff(fullTrees_100[[1]]$tip.label,famTipNames_Bat[[j]])
	fullTrees_100_batFams_i[[i]]<-drop.tip(fullTrees_100[[i]],toDrop)
}
write.nexus(fullTrees_100_batFams_i,file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_batFAMS_",famNames_Bat[j],".trees",sep=""))
}

ltt95_allFAMs_B<-vector("list", length(famNames_Bat))
fullTrees_100_allFAMs_B<-vector("list", length(famNames_Bat))

for (j in 1:length(famNames_Bat)){
	fullTrees_100_i<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_batFAMS_",famNames_Bat[j],".trees",sep=""))
	fullTrees_100_allFAMs_B[[j]]<-fullTrees_100_i
	ltt95_allFAMs_B[[j]] <- ltt95_mod(fullTrees_100_i, alpha=0.05, log=TRUE, method="times", mode="mean")
}



# PLOT the four-part figure.

# plot
jpeg(file="LTT_95polygon_PlacentalsLTTpolygon_eachOrd_wNames_quadMarsupRodBat_NDexp.jpg", width=10, height=10,units="in", res=600, quality=100)

#png(file="LTT_95polygon_PlacentalsLTTpolygon_eachOrd_wNames_quadMarsupRodBat_NDexp.png", width=10, height=10,units="in", res=600)
#pdf(file="LTT_95polygon_PlacentalsLTTpolygon_eachOrd_wNames_quadMarsupRodBat_NDexp.pdf", width=10, height=10,onefile=TRUE)
#tiff(file="LTT_95polygon_PlacentalsLTTpolygon_eachOrd_wNames_quadMarsupRodBat_NDexp.tiff", width=10, height=10,units="in", res=450)
xLim<-c(-110,13)

#layout(matrix(c(1:2), 1, 2, byrow = TRUE))
op <- par(mfrow = c(2,2),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)

# PLACENTALS
cols<-palette(rainbow(length(ordTipNames_P),alpha=0.3))
cols<-palette(rainbow(length(ordTipNames_P),alpha=0.3))

x<-ltt95_placentals
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="", ylab="",xaxt="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
text(x=0,y=nrow(ltt95_placentals[,])-1,labels="PLACENTALIA",cex=0.3,offset=0, font=2,adj=c(0,3.5))#pos=4)
text(x=-115,y=5000,labels="(a) Placental mammals",cex=1,font=2,pos=4)

for (i in 1:3){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
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
# MARSUPIALS
cols<-palette(rainbow(length(ordTipNames_M),alpha=0.3))
cols<-palette(rainbow(length(ordTipNames_M),alpha=0.3))


x<-ltt95_placentals # plot white to placeholder the Y axis
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="", ylab="",yaxt="n",xaxt="n")

x<-ltt95_marsupials
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
text(x=0,y=nrow(ltt95_marsupials[,])-1,labels="MARSUPALIA",cex=0.3,offset=0, font=2,adj=c(0,3.5))#pos=4)
text(x=-115,y=5000,labels="(b) Marsupial mammals",cex=1,font=2,pos=4)

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

#dev.off()
# RODENTS
cols<-palette(rainbow(length(patchNames_Rod),alpha=0.3))
cols<-palette(rainbow(length(patchNames_Rod),alpha=0.3))

x<-ltt95_placentals # plot white to placeholder the Y axis
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="", ylab="")

x<-ltt95_rodents
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
text(x=0,y=nrow(ltt95_rodents[,])-1,labels="RODENTIA",cex=0.3,offset=0, font=2,adj=c(0,3.5))#pos=4)
text(x=-115,y=5000,labels="(c) Rodents",cex=1,font=2,pos=4)

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
# BATS
cols<-palette(rainbow(length(famTipNames_Bat),alpha=0.3))
cols<-palette(rainbow(length(famTipNames_Bat),alpha=0.3))

x<-ltt95_placentals # plot white to placeholder the Y axis
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="", ylab="",yaxt="n")

x<-ltt95_bats
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
text(x=0,y=nrow(ltt95_bats[,])-1,labels="CHIROPTERA",cex=0.3,offset=0, font=2,adj=c(0,3.5))#pos=4)
text(x=-115,y=5000,labels="(d) Bats",cex=1,font=2,pos=4)

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





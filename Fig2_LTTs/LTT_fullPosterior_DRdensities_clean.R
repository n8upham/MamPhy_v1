# packages
library(ape)
library(phytools)

#load(file="LTT_workspace.Rdata")
#save.image(file="LTT_workspace.Rdata")

# modified ltt95() function [phytools] >> faster, doesn't plot.
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
#setwd("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors")
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses")

# specify the backbone to use - FBD or NDexp
bbone <- "NDexp" # "FBD"

# sample of 100 trees (different backbones)
mamPhy_100<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_nexus.trees",sep=""))

# species 'trait' data
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdGenesSampPC_DR-SUMMARY-expanded_",bbone,".txt",sep=""))
colnames(cladesDR)<-c("tiplabel","gen","fam","ord","clade","genes","samp","PC", "harmMeans", "medians", "means", "Count.of.NA", "range", "variance", "stdev", "cv", "sterror")
head(cladesDR)

# 95% confidence intervals around LTTs, as colored transparent polygons
###

# DATA PREP - get 95% CIs of the 100 trees
# MAMMAL SUBCLADE BREAKOUTS
Monotremata<-cladesDR[which(cladesDR$PC=="PC23_Monotremata"),"tiplabel"]
Marsupalia<-cladesDR[which(cladesDR$PC=="PC1_Marsupials"),"tiplabel"]
Placentalia<-setdiff(mamPhy$tip.label,(cladesDR[which(cladesDR$PC=="PC1_Marsupials" | cladesDR$PC=="PC23_Monotremata"),"tiplabel"]))
Afrotheria<-cladesDR[which(cladesDR$PC=="PC2_Afrotheria"),"tiplabel"]
Xenarthra<-cladesDR[which(cladesDR$PC=="PC3_Xenarthra"),"tiplabel"]
Euarchontoglires<-cladesDR[which(cladesDR$ord=="RODENTIA" | cladesDR$ord=="LAGOMORPHA" | cladesDR$ord=="PRIMATES" | cladesDR$ord=="SCANDENTIA" | cladesDR$ord=="DERMOPTERA"),"tiplabel"]
Laurasiatheria<-cladesDR[which(cladesDR$ord=="EULIPOTYPHLA" | cladesDR$ord=="CHIROPTERA" | cladesDR$ord=="PHOLIDOTA" | cladesDR$ord=="CARNIVORA" | cladesDR$ord=="PERISSODACTYLA" | cladesDR$ord=="CETARTIODACTYLA"),"tiplabel"]

mamSubsTips<-list(Monotremata,Marsupalia,Placentalia,Afrotheria,Xenarthra,Euarchontoglires,Laurasiatheria)
mamSubsNames<-c("Monotremata","Marsupalia","Placentalia","Afrotheria","Xenarthra","Euarchontoglires","Laurasiatheria")

mamSubs<-vector("list",length(mamSubsTips))

for (j in 1:length(mamSubsTips)){
    toDrop<-setdiff(mamPhy$tip.label,mamSubsTips[[j]])
    trees100<-vector("list",length=100)
    for (k in 1:length(mamPhy_100)){
    	trees100[[k]]<-drop.tip(mamPhy_100[[k]],toDrop)
    	}
    mamSubs[[j]]<-trees100
    }

for (j in 1:length(mamSubs)){
	write.nexus(mamSubs[[j]], file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",mamSubsNames[j],".trees",sep=""))
	}

# Marsupials v Euarchontoglires v Laurasiatheria
####################
# calculate LTT 95 intervals...
j=2 # Marsupalia
mamPhy_100<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",mamSubsNames[j],".trees",sep=""))
ltt95_Marsupalia <- ltt95_mod(mamPhy_100, alpha=0.05, log=TRUE, method="times", mode="mean")

j=6 # Euarchontoglires
mamPhy_100<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",mamSubsNames[j],".trees",sep=""))
ltt95_Euarchontoglires <- ltt95_mod(mamPhy_100, alpha=0.05, log=TRUE, method="times", mode="mean")

j=7 # Laurasiatheria
mamPhy_100<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",mamSubsNames[j],".trees",sep=""))
ltt95_Laurasiatheria <- ltt95_mod(mamPhy_100, alpha=0.05, log=TRUE, method="times", mode="mean")

j=3 # Placentalia
mamPhy_100<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",mamSubsNames[j],".trees",sep=""))
ltt95_Placentalia <- ltt95_mod(mamPhy_100, alpha=0.05, log=TRUE, method="times", mode="mean")

j=4 # Afrotheria
mamPhy_100<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",mamSubsNames[j],".trees",sep=""))
ltt95_Afrotheria <- ltt95_mod(mamPhy_100, alpha=0.05, log=TRUE, method="times", mode="mean")

j=5 # Xenarthra
mamPhy_100<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_",mamSubsNames[j],".trees",sep=""))
ltt95_Xenarthra <- ltt95_mod(mamPhy_100, alpha=0.05, log=TRUE, method="times", mode="mean")


# Marsupial ORDERS (5)
ordNames<-names(which(sort(table(cladesDR$ord),decreasing=TRUE) > 2))
ordNames_M<-c(ordNames[7:8],ordNames[10],ordNames[13],ordNames[19])

mamPhy_100_ords_M<-vector("list", length(ordNames_M))
ltt95_ords_M<-vector("list", length(ordNames_M))

for (j in 1:length(ordNames_M)){
	mamPhy_100_ords_M[[j]]<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_ORDS_",ordNames_M[j],".trees",sep=""))
	ltt95_ords_M[[j]] <- ltt95_mod(mamPhy_100_ords_M[[j]], alpha=0.05, log=TRUE, method="times", mode="mean")
}

# Euarchontoglires ORDERS
ordNames_E<-c("RODENTIA","PRIMATES","LAGOMORPHA","SCANDENTIA","DERMOPTERA")

mamPhy_100_ords_E<-vector("list", length(ordNames_E))
ltt95_ords_E<-vector("list", length(ordNames_E))

for (j in 1:length(ordNames_E)){
	mamPhy_100_ords_E[[j]]<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_ORDS_",ordNames_E[j],".trees",sep=""))
	ltt95_ords_E[[j]] <- ltt95_mod(mamPhy_100_ords_E[[j]], alpha=0.05, log=TRUE, method="times", mode="mean")
}

# Xenarthra + Afrotheria...
ordNames_XA<-c("AFROSORICIDA","CINGULATA","MACROSCELIDEA", "PILOSA","PROBOSCIDEA", "HYRACOIDEA", "SIRENIA")

mamPhy_100_ords_XA<-vector("list", length(ordNames_XA))
ltt95_ords_XA<-vector("list", length(ordNames_XA))

for (j in 1:length(ordNames_XA)){
	mamPhy_100_ords_XA[[j]]<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_ORDS_",ordNames_XA[j],".trees",sep=""))
	ltt95_ords_XA[[j]] <- ltt95_mod(mamPhy_100_ords_XA[[j]], alpha=0.05, log=TRUE, method="times", mode="mean")
}

# Laurasiatheria ORDERS
ordNames_L<-c("CHIROPTERA","EULIPOTYPHLA","CETARTIODACTYLA","CARNIVORA","PERISSODACTYLA","PHOLIDOTA")

mamPhy_100_ords_L<-vector("list", length(ordNames_L))
ltt95_ords_L<-vector("list", length(ordNames_L))

for (j in 1:length(ordNames_L)){
	mamPhy_100_ords_L[[j]]<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_ORDS_",ordNames_L[j],".trees",sep=""))
	ltt95_ords_L[[j]] <- ltt95_mod(mamPhy_100_ords_L[[j]], alpha=0.05, log=TRUE, method="times", mode="mean")
}

# pre/at-KPg ORDERS
ordNames_preKPg<-c("RODENTIA", "CHIROPTERA", "EULIPOTYPHLA", "PRIMATES", "CETARTIODACTYLA", "AFROSORICIDA", "SCANDENTIA", "MACROSCELIDEA", "PILOSA")

mamPhy_100_preKPg<-vector("list", length(ordNames_preKPg))
ltt95_preKPg<-vector("list",length(ordNames_preKPg))

for (j in 1:length(ordNames_preKPg)){
	mamPhy_100_preKPg[[j]]<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_ORDS_",ordNames_preKPg[j],".trees",sep=""))
	ltt95_preKPg[[j]] <- ltt95_mod(mamPhy_100_preKPg[[j]], alpha=0.05, log=TRUE, method="times", mode="mean")
}

# post-KPg ORDERS
ordNames_postKPg<-c("CARNIVORA","LAGOMORPHA", "PERISSODACTYLA","CINGULATA", "PHOLIDOTA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA") #"DERMOPTERA","TUBULIDENTATA")

mamPhy_100_postKPg<-vector("list", length(ordNames_postKPg))
ltt95_postKPg<-vector("list",length(ordNames_postKPg))

for (j in 1:length(ordNames_postKPg)){
	mamPhy_100_postKPg[[j]]<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_ORDS_",ordNames_postKPg[j],".trees",sep=""))
	ltt95_postKPg[[j]] <- ltt95_mod(mamPhy_100_postKPg[[j]], alpha=0.05, log=TRUE, method="times", mode="mean")
}



#save.image(file="LTT_workspace_Marsup-Euarch-Lauras.Rdata")
load(file="LTT_workspace_Marsup-Euarch-Lauras.Rdata")

ordNames_M2<-c("DIPROTODONTIA", "DIDELPHIMORPHIA", "DASYUROMORPHIA", "PERAMELEMORPHIA", "PAUCITUBERCULATA")
ordNames_M<-c("Diprodontia", "Didelphimorphia", "Dasyuromorphia", "Peramelemorphia", "Paucituberculata")
commonOrdM<-c("kangaroos, wombats", "opossums", "quolls, dunnarts", "bandicoots", "shrew opossums")
namesMarsupials<-cbind(ordNames_M,commonOrdM)

ordNames_E2<-c("RODENTIA","PRIMATES","LAGOMORPHA","SCANDENTIA","DERMOPTERA")
ordNames_E<-c("Rodentia","Primates","Lagomorpha","Scandentia","DERMOPTERA")
commonOrdE<-c("rats, mice", "monkeys, apes", "rabbits, pikas","treeshrews", "flying lemurs")
namesEurarch<-cbind(ordNames_E,commonOrdE)

ordNames_L2<-c("CHIROPTERA","EULIPOTYPHLA","CETARTIODACTYLA","CARNIVORA","PERISSODACTYLA","PHOLIDOTA")
ordNames_L<-c("Chiroptera","Eulipotyphla","Artiodactyla","Carnivora","Perissodactyla","Pholidota")
commonOrdL<-c("bats","shrews, moles","deer, cows, whales","cats, dogs, bears","horses, rhinos","pangolins")
namesLauras<-cbind(ordNames_L,commonOrdL)

###
# PLOT as PLACENTALS -- 9 orders pre/at-KPg vs 10 orders post-KPg == 8 with > 2 species
Eulipotyphla, Afrosoricida, Rodentia, Primates, Scandentia, Pilosa, Macroscelidea, Cetartiodactyla, Chiroptera

ordNames_preKPg<-c("EULIPOTYPHLA", "AFROSORICIDA", "RODENTIA", "PRIMATES", "SCANDENTIA", "PILOSA", "MACROSCELIDEA", "CETARTIODACTYLA", "CHIROPTERA")
ordNames_postKPg<-c("CARNIVORA","PERISSODACTYLA","PHOLIDOTA","LAGOMORPHA","CINGULATA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA") #"DERMOPTERA","TUBULIDENTATA")

ltt95_Euarchontoglires
ltt95_Laurasiatheria


pdf(file=paste("LTT_95polygon_Placentals_pre-vs-postKPg.pdf",sep=""), onefile=TRUE, width=5, height=6)

xLim<-c(-110,25)
quartz(width=5,height=8)
#layout(matrix(c(1:2), 1, 2, byrow = TRUE))
op <- par(mfrow = c(3,1),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)

cexLab1<-0.6
cexLab2<-0.9

cols<-viridis(7,alpha=0.3)
cols<-viridis(7,alpha=0.3)

# Euarchontoglires
cols<-viridis(length(ltt95_ords_E),alpha=0.3)
cols<-viridis(length(ltt95_ords_E),alpha=0.3)

x<-ltt95_Placentalia
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=1,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",xaxt="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
axis(side=1,labels=FALSE, at=c(0,-20,-40,-60,-80,-100))

abline(v=-66,lwd=3,lty=2,col="light grey")
for (i in 1:length(ltt95_ords_E)){
	x<-ltt95_ords_E[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_ords_E[[i]][,])-1,labels=ordNames_E[i],font=2, cex=cexLab1,offset=0,adj=c(0,0))#pos=4)
}

# Laurasiatheria
cols<-viridis(length(ltt95_ords_L),alpha=0.3)
cols<-viridis(length(ltt95_ords_L),alpha=0.3)

x<-ltt95_Placentalia
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=1,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",xaxt="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
axis(side=1,labels=FALSE, at=c(0,-20,-40,-60,-80,-100))

text(x=-115,y=3000,labels="(b)",cex=1,font=2,pos=4)

abline(v=-66,lwd=3,lty=2,col="light grey")
for (i in 1:length(ltt95_ords_L)){
	x<-ltt95_ords_L[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_ords_L[[i]][,])-1,labels=ordNames_L[i],font=2, cex=cexLab1,offset=0,adj=c(0,0))#pos=4)
}





# Xenarthra + Afrotheria
x<-ltt95_ords_XA[[1]]
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=1,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",xaxt="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
axis(side=1,labels=FALSE, at=c(0,-20,-40,-60,-80,-100))
#x<-ltt95_Placentalia
#x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
#lines(x[,3],x[,1],lty=1,lwd=1,type="l")

text(x=-115,y=3000,labels="(a)",cex=1,font=2,pos=4)

abline(v=-66,lwd=3,lty=2,col="light grey")
for (i in 1:length(ltt95_ords_XA)){
	x<-ltt95_ords_XA[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_ords_XA[[i]][,])-1,labels=ordNames_XA[i],font=2, cex=cexLab1,offset=0,adj=c(0,0))#pos=4)
}

library(plotrix)
library(viridis)
library(gplots)

divTimesKPg<-read.table(file="LTT_divTimesToPlot.txt",header=TRUE)
ordinal<-divTimesKPg[which(divTimesKPg$LEVEL=="ordinal"),]
superord<-divTimesKPg[which(divTimesKPg$LEVEL=="superordinal"),]

superord_Bor<-superord[which(superord$majorClade1=="Boreoeutheria"),]
superord_Atlan<-superord[which(superord$majorClade1=="Atlanogenata"),]
superord_Euarch<-superord[which(superord$majorClade1=="Euarchontoglires"),]
superord_Lauras<-superord[which(superord$majorClade1=="Laurasiatheria"),]


########
pdf(file=paste("LTT_95polygon_Placentals_pre-vs-postKPg.pdf",sep=""), onefile=TRUE, width=5, height=6)

xLim<-c(-110,25)
quartz(width=5,height=6)
#layout(matrix(c(1:2), 1, 2, byrow = TRUE))
op <- par(mfrow = c(2,1),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)

cexLab1<-0.6
cexLab2<-0.9

# pre/at-KPg Placentalia 
cols<-viridis(length(ordNames_preKPg),alpha=0.3)
cols2<-viridis(length(ordNames_preKPg),alpha=1)

x<-ltt95_Placentalia
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=1,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",xaxt="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
axis(side=1,labels=FALSE, at=c(0,-20,-40,-60,-80,-100))

text(x=-115,y=3000,labels="(a)",cex=1,font=2,pos=4)

abline(v=-66,lwd=3,lty=2,col="light grey")
for (i in 1:length(ltt95_preKPg)){
	x<-ltt95_preKPg[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_preKPg[[i]][,])-1,labels=ordNames_preKPg[i],font=2, cex=cexLab1,offset=0,adj=c(0,0))#pos=4)
}

plotCI(add=TRUE, gap=0, x=rev(-ordinal$means),y=(exp(seq(3,8,length.out=9))), ui=rev(-ordinal$hpd95_mins),li=rev(-ordinal$hpd95_maxs), cex=1, err="x",sfrac=0, col=rev(cols2),lwd=2, pch=16)


# post-KPg Placentalia 
cols<-viridis(length(ordNames_postKPg),alpha=0.3)
cols2<-viridis(length(ordNames_postKPg),alpha=1)
colsSuper<-c("steelblue1","steelblue2","steelblue3","steelblue4") #rich.colors(4)

x<-ltt95_Placentalia
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=1,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",xaxt="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
axis(side=1,labels=FALSE, at=c(0,-20,-40,-60,-80,-100))
x<-ltt95_Placentalia
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")

text(x=-115,y=3000,labels="(b)",cex=1,font=2,pos=4)

abline(v=-66,lwd=2,lty=2,col="dark grey")
for (i in 1:length(ltt95_postKPg)){
	x<-ltt95_postKPg[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_postKPg[[i]][,])-1,labels=ordNames_postKPg[i],font=2, cex=cexLab1,offset=0,adj=c(0,0))#pos=4)
}

lwdBar<-2
plotCI(add=TRUE, gap=0, x=rev(-superord_Atlan$means),y=(exp(seq(log(1000),log(4500),length.out=6))), ui=rev(-superord_Atlan$hpd95_mins),li=rev(-superord_Atlan$hpd95_maxs), cex=1, err="x",sfrac=0, col=colsSuper[1],lwd=lwdBar, pch=16)
plotCI(add=TRUE, gap=0, x=rev(-superord_Bor$means),y=(500), ui=rev(-superord_Bor$hpd95_mins),li=rev(-superord_Bor$hpd95_maxs), cex=1, err="x",sfrac=0, col=colsSuper[2],lwd=lwdBar, pch=16)
plotCI(add=TRUE, gap=0, x=rev(-superord_Euarch$means),y=(exp(seq(log(300),log(150),length.out=3))), ui=rev(-superord_Euarch$hpd95_mins),li=rev(-superord_Euarch$hpd95_maxs), cex=1, err="x",sfrac=0, col=colsSuper[3],lwd=lwdBar, pch=16)
plotCI(add=TRUE, gap=0, x=rev(-superord_Lauras$means),y=(exp(seq(log(80),log(20),length.out=5))), ui=rev(-superord_Lauras$hpd95_mins),li=rev(-superord_Lauras$hpd95_maxs), cex=1, err="x",sfrac=0, col=colsSuper[4],lwd=lwdBar, pch=16)

dev.off()



superord_Atlan
superord_Bor

superord_Euarch
superord_Lauras



# Euarchontoglires 
cols<-palette(rainbow(length(ordNames_E[1:4]),alpha=0.3))
cols<-palette(rainbow(length(ordNames_E[1:4]),alpha=0.3))

x<-ltt95_Euarchontoglires
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=1,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",xaxt="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
axis(side=1,labels=FALSE, at=c(0,-20,-40,-60,-80,-100))
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
#text(x=0,y=nrow(ltt95_Euarchontoglires[,])-1,labels="Euarchontoglires",cex=cexLab,offset=0, font=2,adj=c(0,1))#pos=4)

text(x=-115,y=3000,labels="(b) Euarchontoglires",cex=1,font=2,pos=4)
for(i in 0:5){
	abline(v=-100+(20*i),lty=2,col="light grey")
}

for (i in 1:4){
	x<-ltt95_ords_E[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_ords_E[[i]][,])-1,labels=namesEurarch[,1][i],font=2,cex=cexLab1,offset=0,adj=c(0,0))#pos=4)
}
for (i in 1:4){
	text(x=0,y=nrow(ltt95_ords_E[[i]][,])-1,labels=namesEurarch[,2][i],cex=cexLab2,offset=0,adj=c(-0.05,1))#pos=4)
}


# Laurasiatheria
cols<-palette(rainbow(length(ordNames_L),alpha=0.3))
cols<-palette(rainbow(length(ordNames_L),alpha=0.3))

x<-ltt95_Euarchontoglires
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=1,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",xaxt="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
axis(side=1,labels=TRUE, at=c(0,-20,-40,-60,-80,-100))
x<-ltt95_Laurasiatheria
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
#text(x=0,y=nrow(ltt95_Laurasiatheria[,])-1,labels="Laurasiatheria",cex=cexLab,offset=0, font=2,adj=c(0,1))#pos=4)

text(x=-115,y=3000,labels="(c) Laurasiatheria",cex=1,font=2,pos=4)
for(i in 0:5){
	abline(v=-100+(20*i),lty=2,col="light grey")
}

for (i in 1:6){
	x<-ltt95_ords_L[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_ords_L[[i]][,])-1,labels=namesLauras[,1][i],font=2,cex=cexLab1,offset=0,adj=c(0,0))#pos=4)
}
for (i in 1:6){
	text(x=0,y=nrow(ltt95_ords_L[[i]][,])-1,labels=namesLauras[,2][i],cex=cexLab2,offset=0,adj=c(-0.05,1))#pos=4)
}


title(xlab = "Time before present (Ma)",
      ylab = "(log) Number of lineages",
      outer = TRUE, line = 3,cex.axis=1,cex.lab=1,font.axis=1, font.lab=2)
par(op)

dev.off()








#######
# PLOT the THREE-part figure. -- Marsupial v Euarchontoglires v Laurasiatheria

# plot
jpeg(file=paste("LTT_95polygon_Marsup-Euarch-Lauras_eachOrd_triVert_",bbone,".jpg",sep=""), width=5, height=9,units="in", res=600, quality=100)

pdf(file=paste("LTT_95polygon_Marsup-Euarch-Lauras_eachOrd_triVert_",bbone,".pdf",sep=""), onefile=TRUE, width=5, height=9)

pdf(file=paste("LTT_95polygon_Marsup-Euarch-Lauras_eachOrd_triVert_",bbone,"withPlacentals.pdf",sep=""), onefile=TRUE, width=5, height=9)

pdf(file=paste("LTT_95polygon_Marsup-Euarch-Lauras_eachOrd_triVert_",bbone,"withPlacentals_withTC-TVsymbols.pdf",sep=""), onefile=TRUE, width=5, height=9)

xLim<-c(-110,25)

#layout(matrix(c(1:2), 1, 2, byrow = TRUE))
op <- par(mfrow = c(3,1),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)

cexLab1<-1
cexLab2<-0.9

# MARSUPIALS 
cols<-palette(rainbow(length(ordNames_M),alpha=0.3))
cols<-palette(rainbow(length(ordNames_M),alpha=0.3))

x<-ltt95_Euarchontoglires
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=1,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",xaxt="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
axis(side=1,labels=FALSE, at=c(0,-20,-40,-60,-80,-100))
x<-ltt95_Placentalia
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
x<-ltt95_Marsupalia
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
#text(x=0,y=nrow(ltt95_Marsupalia[,])-1,labels="All marsupials",cex=cexLab,offset=0, font=2,adj=c(0,3.5))#pos=4)


text(x=-115,y=3000,labels="(a) Marsupials",cex=1,font=2,pos=4)

for(i in 0:5){
	abline(v=-100+(20*i),lty=2,col="light grey")
}
for (i in 1:5){
	x<-ltt95_ords_M[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_ords_M[[i]][,])-1,labels=namesMarsupials[,1][i],font=2, cex=cexLab1,offset=0,adj=c(0,0))#pos=4)
}
for (i in 1:5){
	text(x=0,y=nrow(ltt95_ords_M[[i]][,])-1,labels=namesMarsupials[,2][i],cex=cexLab2,offset=0,adj=c(-0.05,1))#pos=4)
}


# Euarchontoglires 
cols<-palette(rainbow(length(ordNames_E[1:4]),alpha=0.3))
cols<-palette(rainbow(length(ordNames_E[1:4]),alpha=0.3))

x<-ltt95_Euarchontoglires
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=1,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",xaxt="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
axis(side=1,labels=FALSE, at=c(0,-20,-40,-60,-80,-100))
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
#text(x=0,y=nrow(ltt95_Euarchontoglires[,])-1,labels="Euarchontoglires",cex=cexLab,offset=0, font=2,adj=c(0,1))#pos=4)

text(x=-115,y=3000,labels="(b) Euarchontoglires",cex=1,font=2,pos=4)
for(i in 0:5){
	abline(v=-100+(20*i),lty=2,col="light grey")
}

for (i in 1:4){
	x<-ltt95_ords_E[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_ords_E[[i]][,])-1,labels=namesEurarch[,1][i],font=2,cex=cexLab1,offset=0,adj=c(0,0))#pos=4)
}
for (i in 1:4){
	text(x=0,y=nrow(ltt95_ords_E[[i]][,])-1,labels=namesEurarch[,2][i],cex=cexLab2,offset=0,adj=c(-0.05,1))#pos=4)
}


# Laurasiatheria
cols<-palette(rainbow(length(ordNames_L),alpha=0.3))
cols<-palette(rainbow(length(ordNames_L),alpha=0.3))

x<-ltt95_Euarchontoglires
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=1,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",xaxt="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
axis(side=1,labels=TRUE, at=c(0,-20,-40,-60,-80,-100))
x<-ltt95_Laurasiatheria
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
#text(x=0,y=nrow(ltt95_Laurasiatheria[,])-1,labels="Laurasiatheria",cex=cexLab,offset=0, font=2,adj=c(0,1))#pos=4)

text(x=-115,y=3000,labels="(c) Laurasiatheria",cex=1,font=2,pos=4)
for(i in 0:5){
	abline(v=-100+(20*i),lty=2,col="light grey")
}

for (i in 1:6){
	x<-ltt95_ords_L[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_ords_L[[i]][,])-1,labels=namesLauras[,1][i],font=2,cex=cexLab1,offset=0,adj=c(0,0))#pos=4)
}
for (i in 1:6){
	text(x=0,y=nrow(ltt95_ords_L[[i]][,])-1,labels=namesLauras[,2][i],cex=cexLab2,offset=0,adj=c(-0.05,1))#pos=4)
}


title(xlab = "Time before present (Ma)",
      ylab = "(log) Number of lineages",
      outer = TRUE, line = 3,cex.axis=1,cex.lab=1,font.axis=1, font.lab=2)
par(op)

dev.off()





#===========================================
# Old code with Marsupials v Placentals
###################
# load data
# Placentalia
mamPhy_100_placentals<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_PLACENTALS.trees",sep=""))
ltt95_placentals <- ltt95_mod(mamPhy_100_placentals, alpha=0.05, log=TRUE, method="times", mode="mean")

# placental orders (17)
ordNames<-names(which(sort(table(cladesDR$ord),decreasing=TRUE) > 2))
ordNames_P<-c(ordNames[1:6],ordNames[9],ordNames[11:12],ordNames[14:18],ordNames[20:21],ordNames[23])

ordTipNames_P<-vector("list",length(ordNames_P))
for (i in 1:length(ordNames_P)){
	ordTipNames_P[[i]]<-cladesDR[which(cladesDR$ord==ordNames_P[i]),"tiplabel"]
}

ltt95_allOrds<-vector("list", length(ordNames_P))
mamPhy_100_allOrds<-vector("list", length(ordNames_P))

for (j in 1:length(ordNames_P)){
	mamPhy_100_i<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_ORDS_",ordNames_P[j],".trees",sep=""))
	mamPhy_100_allOrds[[j]]<-mamPhy_100_i
	ltt95_allOrds[[j]] <- ltt95_mod(mamPhy_100_i, alpha=0.05, log=TRUE, method="times", mode="mean")
}

# Marsupalia
mamPhy_100_marsupials<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_MARSUPIALS.trees",sep=""))
ltt95_marsupials <- ltt95_mod(mamPhy_100_marsupials, alpha=0.05, log=TRUE, method="times", mode="mean")

# marsupial orders (5)
ordNames<-names(which(sort(table(cladesDR$ord),decreasing=TRUE) > 2))
ordNames_M<-c(ordNames[7:8],ordNames[10],ordNames[13],ordNames[19])

ordTipNames_M<-vector("list",length(ordNames_M))
for (i in 1:length(ordNames_M)){
	ordTipNames_M[[i]]<-cladesDR[which(cladesDR$ord==ordNames_M[i]),"tiplabel"]
}

ltt95_allOrds_M<-vector("list", length(ordNames_M))
mamPhy_100_allOrds_M<-vector("list", length(ordNames_M))

for (j in 1:length(ordNames_M)){
	mamPhy_100_i<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_ORDS_",ordNames_M[j],".trees",sep=""))
	mamPhy_100_allOrds_M[[j]]<-mamPhy_100_i
	ltt95_allOrds_M[[j]] <- ltt95_mod(mamPhy_100_i, alpha=0.05, log=TRUE, method="times", mode="mean")
}

# Rodentia
mamPhy_100_rodentia<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_ORDS_RODENTIA.trees",sep=""))
ltt95_rodents <- ltt95_mod(mamPhy_100_rodentia, alpha=0.05, log=TRUE, method="times", mode="mean")

# rodent patch-clade level (=families + multi-family clades)
patchNames<-names(which(sort(table(cladesDR$PC),decreasing=TRUE) > 2))
patchNames_Rod<-c(patchNames[1:2],patchNames[9:10],patchNames[13],patchNames[17:18],patchNames[21],patchNames[23],patchNames[25])
patchNames_Rod2<-c("MURIDAE", "CRICETIDAE", "SQUIRREL_RELATED", "GUINEAPIG_RELATED", "CASTORIMORPHA", "NESOMYIDAE", "DIPODIDAE", "SPALACIDAE", "ANOMALUROMORPHA", "CALOMYSCIDAE")

ltt95_allPCs_R<-vector("list", length(patchNames_Rod))
mamPhy_100_allPCs_R<-vector("list", length(patchNames_Rod))

for (j in 1:length(patchNames_Rod)){
	mamPhy_100_i<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_PCS_",patchNames_Rod[j],".trees",sep=""))
	mamPhy_100_allPCs_R[[j]]<-mamPhy_100_i
	ltt95_allPCs_R[[j]] <- ltt95_mod(mamPhy_100_i, alpha=0.05, log=TRUE, method="times", mode="mean")
}

# Chiroptera
mamPhy_100_bats<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_ORDS_CHIROPTERA.trees",sep=""))
ltt95_bats <- ltt95_mod(mamPhy_100_bats, alpha=0.05, log=TRUE, method="times", mode="mean")

# bats, fam-level
bats<-cladesDR[ which(cladesDR$ord=="CHIROPTERA"),]
famNames_Bat<-names(which(sort(table(bats$fam),decreasing=TRUE) >2))

famTipNames_Bat<-vector("list",length(famNames_Bat))
for (i in 1:length(famTipNames_Bat)){
	famTipNames_Bat[[i]]<-cladesDR[which(cladesDR$fam==famNames_Bat[i]),"tiplabel"]
}

# generate per-family files
#for (j in 1:length(famNames_Bat)){
#mamPhy_100_batFams_i<-vector("list",length(mamPhy_100))
#
#for (i in 1:length(mamPhy_100)){
#	toDrop<-setdiff(mamPhy_100[[1]]$tip.label,famTipNames_Bat[[j]])
#	mamPhy_100_batFams_i[[i]]<-drop.tip(mamPhy_100[[i]],toDrop)
#}
#write.nexus(mamPhy_100_batFams_i,file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_batFAMS_",famNames_Bat[j],".trees",sep=""))
#}
#
#for (j in 1:length(famNames_Bat)){
#mamPhy_100_batFams_i<-read.nexus(file=paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_batFAMS_",famNames_Bat[j],".trees",sep=""))
#}

ltt95_allFAMs_B<-vector("list", length(famNames_Bat))
mamPhy_100_allFAMs_B<-vector("list", length(famNames_Bat))

for (j in 1:length(famNames_Bat)){
	mamPhy_100_i<-read.nexus(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample100_batFAMS_",famNames_Bat[j],".trees",sep=""))
	mamPhy_100_allFAMs_B[[j]]<-mamPhy_100_i
	ltt95_allFAMs_B[[j]] <- ltt95_mod(mamPhy_100_i, alpha=0.05, log=TRUE, method="times", mode="mean")
}


###
# NEW LABELS for curves
##
ordNames_P
commonOrdP<-c("rodents", "bats", "shrews", "primates", "deer, cows, whales", "bears, seals, cats", "rabbits, pikas", "tenrecs, golden moles", "horses, rhinos", "armadillos", "treeshrews", "elephant shrews", "sloths, anteaters", "pangolins", "elephants", "hyraxes", "manatees")
namesPlacentals<-cbind(ordNames_P,commonOrdP)

ordNames_M
commonOrdM<-c("kangaroos, wombats", "opossums", "quolls, dunnarts", "bandicoots", "shrew opossums")
namesMarsupials<-cbind(ordNames_M,commonOrdM)

patchNames_Rod 
commonPC_R<-c(
namesRodents<-cbind(patchNames_Rod,commonPC_R)

famNames_Bat
commonFam_B<-c(
namesBats<-cbind(famNames_Bat,commonFam_B)



# PLOT the TWO-part figure. -- Placental-Marsupial

# plot
jpeg(file=paste("LTT_95polygon_PlacentalsLTTpolygon_eachOrd_wNames_dualPlacentMarsup_",bbone,"_commonN_vert.jpg",sep=""), width=5, height=9,units="in", res=600, quality=100)

#png(file="LTT_95polygon_PlacentalsLTTpolygon_eachOrd_wNames_quadMarsupRodBat_NDexp.png", width=10, height=10,units="in", res=600)
#pdf(file="LTT_95polygon_PlacentalsLTTpolygon_eachOrd_wNames_quadMarsupRodBat_NDexp.pdf", width=10, height=10,onefile=TRUE)
#tiff(file="LTT_95polygon_PlacentalsLTTpolygon_eachOrd_wNames_quadMarsupRodBat_NDexp.tiff", width=10, height=10,units="in", res=450)
xLim<-c(-110,13)

#layout(matrix(c(1:2), 1, 2, byrow = TRUE))
op <- par(mfrow = c(2,1),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)

cexLab<-0.4
# PLACENTALS
cols<-palette(rainbow(length(ordTipNames_P),alpha=0.3))
cols<-palette(rainbow(length(ordTipNames_P),alpha=0.3))

x<-ltt95_placentals
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="", ylab="",xaxt="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
axis(side=1,labels=FALSE)
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,tRodype="l")
text(x=0,y=nrow(ltt95_placentals[,])-1,labels="All placentals",cex=cexLab,offset=0, font=2,adj=c(0,3.5))#pos=4)
text(x=-115,y=5000,labels="(a) Placental mammals",cex=1,font=2,pos=4)
for(i in 0:5){
	abline(v=-100+(20*i),lty=2,col="light grey")
}

for (i in 1:3){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=namesPlacentals[,2][i],cex=cexLab,offset=0,adj=c(0,0))#pos=4)
}
for (i in 4:4){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=namesPlacentals[,2][i],cex=cexLab,offset=0,adj=c(0,2))#pos=4)
}
for (i in 5:10){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=namesPlacentals[,2][i],cex=cexLab,offset=0,adj=c(0,1))#pos=4)
}
for (i in 11:11){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=namesPlacentals[,2][i],cex=cexLab,offset=0,adj=c(0,2))#pos=4)
}
for (i in 12:12){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=namesPlacentals[,2][i],cex=cexLab,offset=0,adj=c(0,3))#pos=4)
}
for (i in 13:15){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=namesPlacentals[,2][i],cex=cexLab,offset=0,adj=c(0,1))#pos=4)
}
for (i in 16:16){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=namesPlacentals[,2][i],cex=cexLab,offset=0,adj=c(0,-1))#pos=4)
}
for (i in 17:17){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=namesPlacentals[,2][i],cex=cexLab,offset=0,adj=c(0,1))#pos=4)
}

#dev.off()
# MARSUPIALS
cols<-palette(rainbow(length(ordTipNames_M),alpha=0.3))
cols<-palette(rainbow(length(ordTipNames_M),alpha=0.3))


x<-ltt95_placentals # plot white to placeholder the Y axis
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="", ylab="")#,yaxt="n")#,xaxt="n")
#axis(side=1,labels=FALSE)

x<-ltt95_marsupials
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
text(x=0,y=nrow(ltt95_marsupials[,])-1,labels="All marsupials",cex=cexLab,offset=0, font=2,adj=c(0,3.5))#pos=4)
text(x=-115,y=5000,labels="(b) Marsupial mammals",cex=1,font=2,pos=4)
for(i in 0:5){
	abline(v=-100+(20*i),lty=2,col="light grey")
}

x<-ltt95_allOrds_M[[1]]
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[1], border = NA) #awesome.
#lines(x[,3],x[,1],lty=2,lwd=1,type="l")
text(x=0,y=nrow(ltt95_allOrds_M[[1]][,])-1,labels=namesMarsupials[,2][1],cex=cexLab,offset=0, adj=c(0,0))#pos=4)

for (i in 2:5){
	x<-ltt95_allOrds_M[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds_M[[i]][,])-1,labels=namesMarsupials[,2][i],cex=cexLab,offset=0,adj=c(0,0))#pos=4)
}

title(xlab = "Time before present (Ma)",
      ylab = "(log) Number of lineages",
      outer = TRUE, line = 3,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2)
par(op)

dev.off()




##################
# PLOT the four-part figure.

# plot
jpeg(file=paste("LTT_95polygon_PlacentalsLTTpolygon_eachOrd_wNames_quadMarsupRodBat_",bbone,"_commonN.jpg",sep=""), width=10, height=10,units="in", res=600, quality=100)

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
axis(side=1,labels=FALSE)
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
text(x=0,y=nrow(ltt95_placentals[,])-1,labels="All placentals",cex=0.3,offset=0, font=2,adj=c(0,3.5))#pos=4)
text(x=-115,y=5000,labels="(a) Placental mammals",cex=1,font=2,pos=4)
for(i in 0:5){
	abline(v=-100+(20*i),lty=2,col="light grey")
}

for (i in 1:3){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=namesPlacentals[,2][i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
}
for (i in 4:4){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=namesPlacentals[,2][i],cex=0.3,offset=0,adj=c(0,2))#pos=4)
}
for (i in 5:10){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=namesPlacentals[,2][i],cex=0.3,offset=0,adj=c(0,1))#pos=4)
}
for (i in 11:11){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=namesPlacentals[,2][i],cex=0.3,offset=0,adj=c(0,2))#pos=4)
}
for (i in 12:12){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=namesPlacentals[,2][i],cex=0.3,offset=0,adj=c(0,3))#pos=4)
}
for (i in 13:15){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=namesPlacentals[,2][i],cex=0.3,offset=0,adj=c(0,1))#pos=4)
}
for (i in 16:16){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=namesPlacentals[,2][i],cex=0.3,offset=0,adj=c(0,-1))#pos=4)
}
for (i in 17:17){
	x<-ltt95_allOrds[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds[[i]][,])-1,labels=namesPlacentals[,2][i],cex=0.3,offset=0,adj=c(0,1))#pos=4)
}

#dev.off()
# MARSUPIALS
cols<-palette(rainbow(length(ordTipNames_M),alpha=0.3))
cols<-palette(rainbow(length(ordTipNames_M),alpha=0.3))


x<-ltt95_placentals # plot white to placeholder the Y axis
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="", ylab="",yaxt="n",xaxt="n")
axis(side=1,labels=FALSE)

x<-ltt95_marsupials
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=0.8,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) Number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.4,alpha=0.2), border = NA) #awesome.
lines(x[,3],x[,1],lty=1,lwd=1,type="l")
text(x=0,y=nrow(ltt95_marsupials[,])-1,labels="All marsupials",cex=0.3,offset=0, font=2,adj=c(0,3.5))#pos=4)
text(x=-115,y=5000,labels="(b) Marsupial mammals",cex=1,font=2,pos=4)
for(i in 0:5){
	abline(v=-100+(20*i),lty=2,col="light grey")
}

x<-ltt95_allOrds_M[[1]]
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
#plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,font.axis=1, font.lab=2, xlab="Time before present (Ma)", ylab="(log) number of lineages")
polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[1], border = NA) #awesome.
#lines(x[,3],x[,1],lty=2,lwd=1,type="l")
text(x=0,y=nrow(ltt95_allOrds_M[[1]][,])-1,labels=namesMarsupials[,2][1],cex=0.3,offset=0, adj=c(0,0))#pos=4)

for (i in 2:5){
	x<-ltt95_allOrds_M[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_allOrds_M[[i]][,])-1,labels=namesMarsupials[,2][i],cex=0.3,offset=0,adj=c(0,0))#pos=4)
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
for(i in 0:5){
	abline(v=-100+(20*i),lty=2,col="light grey")
}

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
for(i in 0:5){
	abline(v=-100+(20*i),lty=2,col="light grey")
}

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





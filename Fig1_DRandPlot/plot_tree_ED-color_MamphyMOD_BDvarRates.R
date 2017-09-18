setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
	
########################################################
#### Libraries and functions 
########################################################

# Required libraries
library(ape)
library(caper)
# gplots is needed for the following functions: rich.colors
library(gplots) ## need this one

# plotrix is needed for the following functions: draw.circle
library(plotrix)

# diversitree is needed for the internal plotting functions extracted below
library(diversitree) # last run using version 0.9-3
plot2.phylo <- diversitree:::plot2.phylo
radial.text <- diversitree:::radial.text

radial.text <- function (r, theta, labels, cex = 1, col = "black", font = 1, 
    ...) 
{
    n <- length(labels)
    col <- rep(col, length.out = n)
    x <- r * cos(theta)
    y <- r * sin(theta)
    srt <- theta/(2 * pi) * 360
    adj <- rep(0, n)
    i <- srt > 90 & srt < 270
    adj[i] <- 1
    srt[i] <- (srt[i] + 180)%%360
    if (length(font) > 1) 
        font <- font[1]
    if (length(cex) > 1) 
        cex <- cex[1]
    for (i in seq_len(n)) text(x[i], y[i], labels[i], cex = cex, 
        col = col[i], font = font, srt = srt[i], adj = adj[i], 
        ...)
}

group.label.tip.rad3 <- function (obj, lab, col.bar, col.lab, lwd = 1, offset.bar = 0, 
    offset.lab = 0, cex = 1, font = 1, check = FALSE, quiet = FALSE, lend="butt", arc.bar.width=5,...) 
{
    n.taxa <- obj$n.taxa
    np <- obj$Ntip
    n <- obj$n.spp
    if (is.null(n.taxa)) 
        dt <- 1/6/n * 2 * pi
    else dt <- (n.taxa/2 - 0.5 + 1/6)/n * 2 * pi
    theta <- obj$xy$theta[seq_len(obj$Ntip)]
    t0 <- tapply(theta - dt, lab, min)
    t1 <- tapply(theta + dt, lab, max)
    str <- names(t0)
    if (check) {
        i <- order(t0)
        t0 <- t0[i]
        t1 <- t1[i]
        str <- str[i]
        g <- integer(length(t0))
        g[1] <- j <- 1
        end <- t1[1]
        for (i in seq_along(g)[-1]) {
            if (t0[i] > end) {
                j <- j + 1
                end <- t1[i]
            }
            else {
                end <- max(end, t1[i])
            }
            g[i] <- j
        }
        tg <- table(g)
        if (any(tg > 1)) {
            if (!quiet) {
                err <- sapply(which(tg != 1), function(x) paste(str[g == 
                  x], collapse = ", "))
                warn <- c("Collapsing non-monophyletic groups:", 
                  sprintf("\t%s", err))
                warning(paste(warn, collapse = "\n"))
            }
            t0 <- tapply(t0, g, min)
            t1 <- tapply(t1, g, max)
            str <- as.character(tapply(str, g, collapse))
        }
    }
    tm <- (t0 + t1)/2
    r.bar <- rep(max(obj$xx) + offset.bar, length(t0))
    r.lab <- rep(max(obj$xx) + offset.lab, length(t0))
    arcs(t0, t1, r.bar, col = col.bar, lwd = lwd, lend=lend, np=np, arc.bar.width=arc.bar.width)
    if (any(!is.na(col.lab))) 
        radial.text(r.lab, tm, str, col = col.lab, font = font, 
            cex = cex, ...)
}


arcs <- function (theta0, theta1, r, col, lty = 1, lwd = par("lwd"), 
    np=np, lend="butt", arc.bar.width) 
{
    if (length(r) != length(theta0) || length(r) != length(theta1)) 
        stop("theta0, theta1 and r must be of same length")
    if (any(lty[!is.na(lty)] != 1)) 
        warning("lwd != 1 will probably be ugly for arcs")
    theta0 <- sort(theta0)
    theta1 <- sort(theta1) 
    dx <- max(r) * 2 * pi/np
    nn <- pmax(2, ceiling((theta1 - theta0) * r/2/dx))
    tmp <- lapply(seq_along(nn), function(i) cbind(seq(theta0[i], 
        theta1[i], length = nn[i]), rep(r[i], nn[i])))
    tmp0 <- do.call(rbind, lapply(seq_along(nn), function(i) rbind(tmp[[i]][-nn[i], 
        ])))
    tmp1 <- do.call(rbind, lapply(seq_along(nn), function(i) rbind(tmp[[i]][-1, 
        ])))
    if (length(lwd) > 1) 
        lwd <- rep(rep(lwd, length = length(theta0)), nn - 1)
    if (length(col) > 1)
        col <- rep(rep(col, length = length(theta0)), nn - 1)
     {
	segments(x0=tmp0[, 2] * cos(tmp0[, 1]), y0=tmp0[, 2] * sin(tmp0[, 
        1]), x1=tmp0[, 2] * cos(tmp0[, 1])*arc.bar.width, y1=tmp0[, 2] * sin(tmp0[, 1])*arc.bar.width, lwd = lwd, col = col, lend=lend)
        }
}

color.legend <- function (xl, yb, xr, yt, legend, rect.col, cex = 1, align = "lt", 
gradient = "x", border="grey", lwd=0.1,...) {

	oldcex <- par("cex")
	par(xpd = TRUE, cex = cex)
	gradient.rect(xl, yb, xr, yt, col = rect.col, nslices = length(rect.col), gradient = gradient, border=border, lwd=lwd)
	if (gradient == "x") {
		xsqueeze <- (xr - xl)/(2 * length(rect.col))
		textx <- seq(xl + xsqueeze, xr - xsqueeze, length.out = length(legend))
			if (match(align, "rb", 0)) {
				texty <- yb - 0.2 * strheight("O")
				textadj <- c(0.5, 1)
				}
			else {
				texty <- yt + 0.2 * strheight("O")
				textadj <- c(0.5, 0)
			}
		}
	else {
	ysqueeze <- (yt - yb)/(2 * length(rect.col))
	texty <- seq(yb + ysqueeze, yt - ysqueeze, length.out = length(legend))
	if (match(align, "rb", 0)) {
	textx <- xr + 0.2 * strwidth("O")
	textadj <- c(0, 0.5)
	}
	else {
	textx <- xl - 0.2 * strwidth("O")
	textadj <- c(1, 0.5)
	}
	}
	text(textx, texty, labels = legend, adj = textadj, ...)
	par(xpd = FALSE, cex = oldcex)
}

gradient.rect <- function (xleft, ybottom, xright, ytop, reds, greens, blues, col = NULL, nslices = 50, gradient = "x", border = par("fg"), lwd) {
	if (is.null(col)) 
	col <- color.gradient(reds, greens, blues, nslices)
	else nslices <- length(col)
	nrect <- max(unlist(lapply(list(xleft, ybottom, xright, ytop), 
	length)))
	if (nrect > 1) {
		if (length(xleft) < nrect) 
		xleft <- rep(xleft, length.out = nrect)
		if (length(ybottom) < nrect) 
		ybottom <- rep(ybottom, length.out = nrect)
		if (length(xright) < nrect) 
		xright <- rep(xright, length.out = nrect)
		if (length(ytop) < nrect) 
		ytop <- rep(ytop, length.out = nrect)
		for (i in 1:nrect) gradient.rect(xleft[i], ybottom[i], 
		xright[i], ytop[i], reds, greens, blues, col, nslices, 
		gradient, border = border)
		}
		else {
			if (gradient == "x") {
				xinc <- (xright - xleft)/nslices
				xlefts <- seq(xleft, xright - xinc, length = nslices)
				xrights <- xlefts + xinc
				rect(xlefts, ybottom, xrights, ytop, col = col, lty = 0, lwd=lwd)
				rect(xlefts[1], ybottom, xrights[nslices], ytop, 
				border = border, lwd=lwd)
				}
				else {
					yinc <- (ytop - ybottom)/nslices
					ybottoms <- seq(ybottom, ytop - yinc, length = nslices)
					ytops <- ybottoms + yinc
					rect(xleft, ybottoms, xright, ytops, col = col, lty = 0, lwd=lwd)
					rect(xleft, ybottoms[1], xright, ytops[nslices], border = border, lwd=lwd)
					}
					}
					invisible(col)
					}




	
########################################################
#### Preparing the data for plotting
########################################################

# read in the tree on which data will be plotted
plottree <- drop.tip(read.nexus("MamPhy_fullPosterior_BDvarRates_17Exp_all10k_nexus_MCC-targetHeights.tre"),"_Anolis_carolinensis")
plottree <- ladderize(plottree, right=FALSE)


plottree <- drop.tip(read.tree("MamPhy_fullPosterior_BDvr_pcsFIXED_FBD_sample1.tre"),"_Anolis_carolinensis")
plottree <- ladderize(plottree, right=FALSE)


plottree <- drop.tip(read.tree("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_sample1.tre"),"_Anolis_carolinensis")
plottree <- ladderize(plottree, right=FALSE)


# read in the data for colouring branches
data <- read.table("DR-SUMMARY_MamPhy_BDvrFIXED_NDexp_all10k_5912species.txt", row.names=1, header=TRUE)
data <- data[2:5912,]

grp <- read.table("5911sp_mammals_spFamOrd.txt", header=TRUE)

fams <- grp[,2][match(plottree$tip.label,grp[,1])]
ords <- grp[,3][match(plottree$tip.label,grp[,1])]


## THE SAME, but for the pruned data -- 4098 species
####

# read in the tree on which data will be plotted
missing<-read.table("MamPhy_FIN4_1813sp_missing_LIST.txt", header=FALSE)

plottree <- drop.tip(plottree,as.character(missing$V1))

plottree <- drop.tip(plottree,"_Anolis_carolinensis")
plottree <- ladderize(plottree, right=FALSE)

# read in the data for colouring branches
data <- read.table("DR-SUMMARY_MamPhy_BDvarRates_17Exp_all10ktrees_pruned4098sp.txt", row.names=1, header=TRUE)
data <- data[2:4099,]

data <- read.table("DR-SUMMARY_MamPhy_BDvarRates_17Uni_all10ktrees_pruned4098sp.txt", row.names=1, header=TRUE)
data <- data[2:4099,]

grp <- read.table("MamPhy_plotting_groups_prune4098.txt", header=TRUE)

fams <- grp[,2][match(plottree$tip.label,grp[,1])]
ords <- grp[,3][match(plottree$tip.label,grp[,1])]


########################################################
#### Ancestral states for each trait
########################################################

DR_median <- data[,"medians"]

DR_median <- data[,"means"]


DR_median <- data[,"cv"]*100

DR_median <- data[,"range"]

DR_median <- data[,"variance"]


DR_median <- data[,"harmMeans"]

names(DR_median) <- rownames(data)


DR_median_ordered <- DR_median[match(plottree$tip.label,names(DR_median))]

DR_median_anc <- ace(DR_median_ordered, plottree, method="pic", scaled=TRUE)$ace



####
# DO a BOXPLOT of the missing vs sampled species in the tip values
##
missing<-read.table("MamPhy_FIN4_1812sp_missing_LIST.txt", header=FALSE)

DR_median_missing <- DR_median[match(missing$V1,names(DR_median))]

DR_median_sampled <- DR_median[match(setdiff(plottree$tip.label,missing$V1),names(DR_median))]

miss<-cbind(DR_median_missing,rep(1,length(DR_median_missing)))
colnames(miss)<-c("DR","factor")
rownames(miss)<-NULL

samp<-cbind(DR_median_sampled,rep(2,length(DR_median_sampled)))
colnames(samp)<-c("DR","factor")
rownames(samp)<-NULL

all<-rbind(miss,samp)

boxplot(DR ~ factor, data=all, plot=TRUE, names=c("Imputed species","Sampled species"), ylab="Coefficient of variation (CV)", main="Per-species variation in DR across 10k trees (5910 sp., Exp)")

boxplot(DR ~ factor, data=all, plot=TRUE, names=c("Imputed species","Sampled species"), ylab="Range in DR", main="Per-species range in DR across 10k trees (5910 sp., Exp)")

boxplot(DR ~ factor, data=all, plot=TRUE, names=c("Imputed species","Sampled species"), ylab="Variance in DR", main="Per-species variance in DR across 10k trees (5910 sp., Exp)")

boxplot(DR ~ factor, data=all, plot=TRUE, names=c("Imputed species","Sampled species"), ylab="Harmonic mean", main="Per-species harmonicMean in DR across 10k trees (5910 sp., Exp)")

boxplot(DR ~ factor, data=all, plot=TRUE, names=c("Imputed species","Sampled species"), ylab="Median", main="Per-species median in DR across 10k trees (5910 sp., Exp)")

boxplot(DR ~ factor, data=all, plot=TRUE, names=c("Imputed species","Sampled species"), ylab="Mean (arithmetic)", main="Per-species mean in DR across 10k trees (5910 sp., Exp)")

########################################################
# Plot 1 - HFPMedian
########################################################

# FOR FULL TREE OF 5911 SPECIES
#####
# Match ancestors totree edges
    match.anc <- match(as.numeric(names(DR_median_anc)), plottree$edge[,2])[-1]

    # Assign rates to each internal node
    reconRate <- vector(mode="numeric", length=length(plottree$edge[,2]))
    reconRate[match.anc] <- DR_median_anc[2:5910] #[2:7237]

    # Assign rates to tips
    tip.idx <- sort(plottree$tip.label, index.return=TRUE)$ix

    reconRate[match(tip.idx, plottree$edge[,2])] <- DR_median[sort(names(DR_median))]

    # Create colour palette
    reconColors <- reconRate

    range <- quantile(reconRate, seq(0,1, 0.01))[2:101]
    for (i in 1:100) {
        if (i==100) {range[i] <- range[i]+0.1}
        if (i==1) {reconColors[which(reconRate <= range[i])] <- palette(rich.colors(100))[i] } 
        if (i > 1) {reconColors[which((reconRate <= range[i]) & (reconRate > range[i-1]))] <- palette(rich.colors(100))[i]}
        }

    # plotting & tabling
    th <- max(branching.times(plottree))
    gr.col <- gray.colors(2, start=0.90, end=0.95)

# FOR PRUNED TREE OF 4098 SPECIES
#####
# Match ancestors totree edges
     match.anc <- match(as.numeric(names(DR_median_anc)), plottree$edge[,2])[-1]

    # Assign rates to each internal node
    reconRate2 <- vector(mode="numeric", length=length(plottree$edge[,2]))
    reconRate2[match.anc] <- DR_median_anc[2:4097] #[2:7237]

    # Assign rates to tips
    tip.idx <- sort(plottree$tip.label, index.return=TRUE)$ix

    reconRate2[match(tip.idx, plottree$edge[,2])] <- DR_median[sort(names(DR_median))]

    # Create colour palette
    reconColors2 <- reconRate2

   range2 <- quantile(reconRate2, seq(0,1, 0.01))[2:101]
    for (i in 1:100) {
    	if (i==100) {range2[i] <- range2[i]+0.1}
    	if (i==1) {reconColors2[which(reconRate2 <= range2[i])] <- palette(rich.colors(100))[i] } 
    	if (i > 1) {reconColors2[which((reconRate2 <= range2[i]) & (reconRate2 > range2[i-1]))] <- palette(rich.colors(100))[i]}
    	}

    # plotting & tabling
    th <- max(branching.times(plottree))
    gr.col <- gray.colors(2, start=0.90, end=0.95)

####
# NO TIPS...
###

pdf(file="DR_onMamPhy_BDvarRates_5910_all10k_LAST_median.pdf", width=25, height=25, onefile=TRUE)
pdf(file="DR_onMamPhy_BDvarRates_Exp_5910taxa_all10k_MCC-tH_harmMean_noTips.pdf", width=25, height=25, onefile=TRUE)

pdf(file="DR_onMamPhy_BDvr_pcsFIXED_FBD_5911taxa_sample1_harmMean_noTips.pdf", width=25, height=25, onefile=TRUE)

pdf(file="DR_onMamPhy_BDvr_pcsFIXED_NDexp_5911taxa_sample1_harmMean_noTips_sizeFBD.pdf", width=25, height=25, onefile=TRUE)

pdf(file="DR_onMamPhy_BDvr_pcsFIXED_NDexp_5911taxa_sample1_harmMean_noTips.pdf", width=25, height=25, onefile=TRUE)

XX1=-220
XX2=220
##FBD
#XX1=-300
#XX2=300

plot(plottree, show.tip.label=FALSE, type="f", edge.width=1, no.margin=FALSE, root.edge=TRUE, edge.color="white", plot=FALSE, x.lim=c(XX1,XX2), y.lim=c(XX1,XX2))
draw.circle(0,0, radius=th, col =gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-23.03, col=gr.col[2], border=gr.col[2])
draw.circle(0,0, radius=th-66, col =gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-145, col =gr.col[2], border=gr.col[2])

par(new=T)
obj <- plot2.phylo(plottree, show.tip.label=FALSE, cex=0.05, label.offset=0.05, type="f", edge.width=0.65, no.margin=FALSE, root.edge=TRUE, edge.color=as.matrix(reconColors),x.lim=c(XX1,XX2), y.lim=c(XX1,XX2))

par(new=T)

## for ND small
#color.legend(-61.0, -50.5, 65.7, -43.5, legend=NULL, rect.col= palette(rich.colors(100)), gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey")
#text(x=26, y=-66, "Diversification rate (species per Myr)", cex=2) #"ED (Millions of Years)", cex=1.5)
#par(new=T)
#group.label.tip.rad3(obj, lab=ords, c( "black","light grey"), "black", offset.bar=2, offset.lab=8, cex=1, lwd=4, arc.bar.width=1.02)
#par(fig=c(0.40, 0.70, 0.44, 0.54),mar=c(0,0,0,0), mar=c(0,0,0,0), new=T)

# for ND normal ##
color.legend(-34.3, -62.5, 57.7, -57.5, legend=NULL, rect.col= palette(rich.colors(100)), gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey")
text(x=26, y=-76, "Diversification rate (species per Myr)", cex=2) 
par(new=T)
group.label.tip.rad3(obj, lab=ords, c( "black","light grey"), "black", offset.bar=2, offset.lab=8, cex=1, lwd=4, arc.bar.width=1.02)
par(fig=c(0.42, 0.72, 0.39, 0.49),mar=c(0,0,0,0), mar=c(0,0,0,0), new=T)

## for FBD ##
#color.legend(-48.3, -85.5, 71.7, -78.5, legend=NULL, rect.col= palette(rich.colors(100)), gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey")
#text(x=26, y=-106, "Diversification rate (species per Myr)", cex=2) 
#par(new=T)
#group.label.tip.rad3(obj, lab=ords, c( "black","light grey"), "black", offset.bar=2, offset.lab=8, cex=1, lwd=4, arc.bar.width=1.02)
#par(fig=c(0.42, 0.72, 0.39, 0.49),mar=c(0,0,0,0), mar=c(0,0,0,0), new=T)

plot(density(DR_median), col="dark grey", main="", bty="n", axes=F, xlim=range(DR_median))
polygon(density(DR_median), col="dark grey", border="dark grey", bty="n")

x.tick <- quantile(DR_median, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick), labels=c(0,round(x.tick,2)), side=1, line=1.3, cex=2.5, lwd=1, tck=-0.05, cex.axis=1.5, mgp=c(1,1,0))

dens.rate <- density(DR_median)$y
axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=2.5, las=1, lwd=1, cex.axis=1.5, tck=-0.05, mgp=c(1,1,0))

dev.off()


###
# Now with TIPLABELS
###

missing<-read.table("MamPhy_FIN4_1813sp_missing_LIST.txt", header=FALSE)

sampled<-setdiff(plottree$tip.label,missing$V1)

tipColors <- rep("black",5911)
tipColors[match(missing$V1,plottree$tip.label)] <- "red" 

pdf(file="DR_onMamPhy_BDvarRates_Exp_5910taxa_all10k_MCC-tH_harmMean_colTips.pdf", width=25, height=25, onefile=TRUE)
pdf(file="DR_onMamPhy_BDvarRates_Exp_5910taxa_all10k_MCC-tH_CV_colTips.pdf", width=25, height=25, onefile=TRUE)

pdf(file="DR_onMamPhy_BDvr_pcsFIXED_FBD_5911taxa_sample1_harmMean_colTips.pdf", width=25, height=25, onefile=TRUE)

pdf(file="DR_onMamPhy_BDvr_pcsFIXED_NDexp_5911taxa_sample1_harmMean_colTips.pdf", width=25, height=25, onefile=TRUE)

XX1=-220
XX2=220
##FBD
#XX1=-300
#XX2=300

plot(plottree, show.tip.label=FALSE, type="f", edge.width=1, no.margin=TRUE, root.edge=TRUE, edge.color="white", plot=FALSE, x.lim=c(XX1,XX2), y.lim=c(XX1,XX2))
draw.circle(0,0, radius=th, col =gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-23.03, col=gr.col[2], border=gr.col[2])
draw.circle(0,0, radius=th-66, col =gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-145, col =gr.col[2], border=gr.col[2])

par(new=T)
obj <- plot2.phylo(plottree, show.tip.label=TRUE, tip.color=tipColors, cex=0.05, label.offset=0.05, type="f", edge.width=0.65, no.margin=TRUE, root.edge=TRUE, edge.color=as.matrix(reconColors),x.lim=c(XX1,XX2), y.lim=c(XX1,XX2))

par(new=T)

# for ND normal ##
color.legend(-27.3, -56.5, 57.7, -51.5, legend=NULL, rect.col= palette(rich.colors(100)), gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey")
text(x=26, y=-71, "Diversification rate (species per Myr)", cex=2) 
par(new=T)
group.label.tip.rad3(obj, lab=ords, c( "black","light grey"), "black", offset.bar=12, offset.lab=19, cex=1, lwd=4, arc.bar.width=1.02)
par(fig=c(0.42, 0.72, 0.39, 0.49),mar=c(0,0,0,0), mar=c(0,0,0,0), new=T)

##for FBD
#color.legend(-39.3, -76.5, 72.7, -70, legend=NULL, rect.col= palette(rich.colors(100)), gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey")
#text(x=26, y=-98, "Diversification rate (species per Myr)", cex=2) #"ED (Millions of Years)", cex=1.5)
##text(x=26, y=-65, "Coefficient of variation in 10,000 trees", cex=2) #"ED (Millions of Years)", cex=1.5)
#par(new=T)
#group.label.tip.rad3(obj, lab=ords, c( "black","light grey"), "black", offset.bar=12, offset.lab=19, cex=1, lwd=4, arc.bar.width=1.02)
#par(fig=c(0.42, 0.72, 0.39, 0.49),mar=c(0,0,0,0), mar=c(0,0,0,0), new=T)

plot(density(DR_median), col="dark grey", main="", bty="n", axes=F, xlim=range(DR_median))
polygon(density(DR_median), col="dark grey", border="dark grey", bty="n")

x.tick <- quantile(DR_median, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick), labels=c(0,round(x.tick,2)), side=1, line=1.3, cex=2.5, lwd=1, tck=-0.05, cex.axis=1.5, mgp=c(1,1,0))

dens.rate <- density(DR_median)$y
axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=2.5, las=1, lwd=1, cex.axis=1.5, tck=-0.05, mgp=c(1,1,0))

dev.off()

#####
# PLOTTING on the 4098 tree but with BOTH polygons in the center... incl the 5910 one in the BACKGROUND.
###

# now the reconRate2 = 4098 tree; reconRate = 5910 tree

pdf(file="DR_onMamPhy_BDvarRates_Exp_4098taxa_all10k_MCC-tH_harmMean_vs5910-2polygons.pdf", width=25, height=25, onefile=TRUE)

pdf(file="DR_onMamPhy_BDvarRates_Uni_4098taxa_all10k_MCC-tH_harmMean_vs5910-2polygons.pdf", width=25, height=25, onefile=TRUE)

XX1=-200
XX2=200

plot(plottree, show.tip.label=FALSE, type="f", edge.width=1, no.margin=TRUE, root.edge=TRUE, edge.color="white", plot=FALSE, x.lim=c(XX1,XX2))#, y.lim=c(YY1,YY2))
draw.circle(0,0, radius=th, col =gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-23.03, col=gr.col[2], border=gr.col[2])
#draw.circle(0,0, radius=th-66, col=gr.col[2], border=gr.col[2])
draw.circle(0,0, radius=th-66, col =gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-145, col =gr.col[2], border=gr.col[2])

par(new=T)
obj <- plot2.phylo(plottree, show.tip.label=TRUE, cex=0.05, label.offset=0.05, type="f", edge.width=0.45, no.margin=TRUE, root.edge=TRUE, edge.color=as.matrix(reconColors2),x.lim=c(XX1,XX2))#), y.lim=c(YY1,YY2))

par(new=T)

color.legend(-24.5, -51.5, 57, -46.5, legend=NULL, rect.col= palette(rich.colors(100)), gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey")

# -34, -49.5, 90, -46.5

text(x=26, y=-65, "Diversification rate (species per Myr)", cex=2) #"ED (Millions of Years)", cex=1.5)

par(new=T)

group.label.tip.rad3(obj, lab=ords, c( "black","light grey"), "black", offset.bar=8, offset.lab=12, cex=1, lwd=4, arc.bar.width=1.02)

par(fig=c(0.42, 0.72, 0.39, 0.49),mar=c(0,0,0,0), mar=c(0,0,0,0), new=T)

plot(density(reconRate), col="dark grey", main="", bty="n", axes=F, xlim=range(reconRate), ylim=range(density(reconRate2)$y))
polygon(density(reconRate), col="dark grey", border="dark grey", bty="n")

col=rgb(0,0,0,0.5)
polygon(density(reconRate2), col=col, border="black", bty="n")

x.tick <- quantile(reconRate2, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick), labels=c(0,round(x.tick,2)), side=1, line=1.3, cex=2.5, lwd=1, tck=-0.05, cex.axis=1.5, mgp=c(1,1,0))

dens.rate <- density(reconRate2)$y
axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=2.5, las=1, lwd=1, cex.axis=1.5, tck=-0.05, mgp=c(1,1,0))

dev.off()





##############
######
# Plotting the 4098 species pruned with DR -- with NO TIPS
#####

pdf(file="DR_onMamPhy_BDvarRates_pruned4098_all10k_LAST_median_noTips.pdf", width=25, height=25, onefile=TRUE)

XX1=-200
XX2=200

plot(plottree, show.tip.label=FALSE, type="f", edge.width=1, no.margin=TRUE, root.edge=TRUE, edge.color="white", plot=FALSE, x.lim=c(XX1,XX2))#, y.lim=c(YY1,YY2))
draw.circle(0,0, radius=th, col =gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-23.03, col=gr.col[2], border=gr.col[2])
#draw.circle(0,0, radius=th-66, col=gr.col[2], border=gr.col[2])
draw.circle(0,0, radius=th-66, col =gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-145, col =gr.col[2], border=gr.col[2])

par(new=T)
obj <- plot2.phylo(plottree, show.tip.label=FALSE, cex=0.05, label.offset=0.05, type="f", edge.width=0.65, no.margin=TRUE, root.edge=TRUE, edge.color=as.matrix(reconColors),x.lim=c(XX1,XX2))#), y.lim=c(YY1,YY2))

par(new=T)

color.legend(-22.3, -51.5, 63.7, -46.5, legend=NULL, rect.col= palette(rich.colors(100)), gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey")

# -34, -49.5, 90, -46.5

text(x=26, y=-65, "Diversification rate (species per Myr)", cex=2) #"ED (Millions of Years)", cex=1.5)

par(new=T)

group.label.tip.rad3(obj, lab=ords, c( "black","light grey"), "black", offset.bar=2, offset.lab=6, cex=1, lwd=4, arc.bar.width=1.02)

par(fig=c(0.42, 0.72, 0.39, 0.49),mar=c(0,0,0,0), mar=c(0,0,0,0), new=T)

plot(density(reconRate), col="dark grey", main="", bty="n", axes=F, xlim=range(reconRate))
polygon(density(reconRate), col="dark grey", border="dark grey", bty="n")

x.tick <- quantile(reconRate, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick), labels=c(0,round(x.tick,2)), side=1, line=1.3, cex=2.5, lwd=1, tck=-0.05, cex.axis=1.5, mgp=c(1,1,0))

dens.rate <- density(reconRate)$y
axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=2.5, las=1, lwd=1, cex.axis=1.5, tck=-0.05, mgp=c(1,1,0))

dev.off()

######
# Plotting the 4098 species pruned with DR -- with YES TIPS
#####

pdf(file="DR_onMamPhy_BDvarRates_pruned4098_all10k_LAST_median_yesTips.pdf", width=25, height=25, onefile=TRUE)

XX1=-200
XX2=200

plot(plottree, show.tip.label=FALSE, type="f", edge.width=1, no.margin=TRUE, root.edge=TRUE, edge.color="white", plot=FALSE, x.lim=c(XX1,XX2))#, y.lim=c(YY1,YY2))
draw.circle(0,0, radius=th, col =gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-23.03, col=gr.col[2], border=gr.col[2])
#draw.circle(0,0, radius=th-66, col=gr.col[2], border=gr.col[2])
draw.circle(0,0, radius=th-66, col =gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-145, col =gr.col[2], border=gr.col[2])

par(new=T)
obj <- plot2.phylo(plottree, show.tip.label=TRUE, cex=0.05, label.offset=0.05, type="f", edge.width=0.65, no.margin=TRUE, root.edge=TRUE, edge.color=as.matrix(reconColors),x.lim=c(XX1,XX2))#), y.lim=c(YY1,YY2))

par(new=T)

color.legend(-22.3, -51.5, 63.7, -46.5, legend=NULL, rect.col= palette(rich.colors(100)), gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey")

# -34, -49.5, 90, -46.5

text(x=26, y=-65, "Diversification rate (species per Myr)", cex=2) #"ED (Millions of Years)", cex=1.5)

par(new=T)

group.label.tip.rad3(obj, lab=ords, c( "black","light grey"), "black", offset.bar=10, offset.lab=14, cex=1, lwd=4, arc.bar.width=1.02)

par(fig=c(0.42, 0.72, 0.39, 0.49),mar=c(0,0,0,0), mar=c(0,0,0,0), new=T)

plot(density(reconRate), col="dark grey", main="", bty="n", axes=F, xlim=range(reconRate))
polygon(density(reconRate), col="dark grey", border="dark grey", bty="n")

x.tick <- quantile(reconRate, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick), labels=c(0,round(x.tick,2)), side=1, line=1.3, cex=2.5, lwd=1, tck=-0.05, cex.axis=1.5, mgp=c(1,1,0))

dens.rate <- density(reconRate)$y
axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=2.5, las=1, lwd=1, cex.axis=1.5, tck=-0.05, mgp=c(1,1,0))

dev.off()











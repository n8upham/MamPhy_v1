setwd("C://Users//Alex//Desktop//amph//amph_shl_new_PASTIS_R")
	
########################################################
#### Libraries and functions 
########################################################

# Required libraries
library(ape)
library(caper)
# gplots is needed for the following functions: rich.colors
library(gplots) 

# plotrix is needed for the following functions: draw.circle
library(plotrix)

# diversitree is needed for the internal plotting functions extracted below
library(diversitree) # last run using version 0.9-3
plot2.phylo <- diversitree:::plot2.phylo
radial.text <- diversitree:::radial.text


group.label.tip.rad3 <- function (obj, lab, col.bar, col.lab, lwd = 1, offset.bar = 0, 
    offset.lab = 0, cex = 1, font = 1, check = FALSE, quiet = FALSE, lend="butt", arc.bar.width=5,
    ...) 
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
plottree <- drop.tip(read.tree("amph_shl_new_Posterior_7238.1.tre"),"Homo_sapiens")
plottree <- ladderize(plottree, right=FALSE)

setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/11_makingFullPosteriors")
last100<-read.tree("MamPhy_fullPosterior_YuleCR_mtDNA_Last100.trees")

plottree <- drop.tip(last100[[1]],"_Anolis_carolinensis")
plottree <- ladderize(plottree, right=FALSE)


#######
# Need to GATHER this data for coloring first...


# read in the data for colouring branches
data <- read.csv("FP.test.1000.medians.csv", row.names=1)
grp <- read.csv("amph_grps.csv", as.is=TRUE, header=TRUE)

fams <- grp[,2][match(plottree$tip.label,grp[,1])]
subs <- grp[,3][match(plottree$tip.label,grp[,1])]
tax <- grp[,4][match(plottree$tip.label,grp[,1])]
lab <- grp[,5][match(plottree$tip.label,grp[,1])]

########################################################
#### Ancestral states for each trait
########################################################

HFPMedian <- log(data[,"x"])[-1]
names(HFPMedian) <- rownames(data)[-1]

HFPMedian_anc <- ace(HFPMedian, plottree, method="pic", scaled=TRUE)$ace

########################################################
# Plot 1 - HFPMedian
########################################################

# Match ancestors totree edges
match.anc <- match(as.numeric(names(HFPMedian_anc)), plottree$edge[,2])[-1]

# Assign rates to each internal node
reconRate <- vector(mode="numeric", length=length(plottree$edge[,2]))
reconRate[match.anc] <- HFPMedian_anc[2:7237]

# Assign rates to tips
tip.idx <- sort(plottree$tip.label, index.return=TRUE)$ix

reconRate[match(tip.idx, plottree$edge[,2])] <- HFPMedian[sort(names(HFPMedian))]

# Create colour palette
reconColors <- reconRate

range <- quantile(reconRate, seq(0,1, 0.01))[2:101]
for (i in 1:100) {
	if (i==100) {range[i] <- range[i]+0.1}
	if (i==1) {reconColors[which(reconRate <= range[i])] <- palette(rich.colors(100))[i] } 
	if (i > 1) {reconColors[which((reconRate <= range[i]) & (reconRate > range[i-1]))] <- palette(rich.colors(100))[i]}
	}

tipColors <- rep(NA, 7238)

# plotting & tabling
th <- max(branching.times(plottree))
gr.col <- gray.colors(2, start=0.90, end=0.95)

pdf(file="FP_amph_1000_1.pdf", width=10, height=10, onefile=TRUE)

plot(plottree, show.tip.label=FALSE, type="f", edge.width=1, no.margin=T, root.edge=TRUE, edge.color="white", plot=FALSE)
draw.circle(0,0, radius=th, col =gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-40, col=gr.col[2], border=gr.col[2])
draw.circle(0,0, radius=th-80, col =gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-120, col=gr.col[2], border=gr.col[2])
draw.circle(0,0, radius=th-160, col =gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-200, col=gr.col[2], border=gr.col[2])

par(new=T)
obj <- plot2.phylo(plottree, show.tip.label=FALSE, type="f", edge.width=0.25, no.margin=TRUE, root.edge=TRUE, edge.color=as.matrix(reconColors))

par(new=T)

color.legend(-51, -113.5, 89.5, -127, legend=NULL, rect.col= palette(rich.colors(100)), gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey")

text(x=22, y=-155, "ED (Millions of Years)", cex=1.5)

par(new=T)

group.label.tip.rad3(obj, lab, c( "black","light grey"), "black", offset.bar=1, offset.lab=10, cex=1.2, lwd=1, arc.bar.width=1.02)

par(fig=c(0.375, 0.74, 0.33, 0.49),mar=c(0,0,0,0), mai=c(0,0,0,0), new=T)

plot(density(reconRate), col="dark grey", main="", bty="n", axes=F, xlim=range(reconRate))
polygon(density(reconRate), col="dark grey", border="dark grey", bty="n")

x.tick <- quantile(reconRate, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick), labels=c(0,round(exp(x.tick), digits=0)), side=1, line=0.8, cex=0.6, lwd=0.5, tck=-0.05, cex.axis=0.8, mgp=c(1,0.25,0))

dens.rate <- density(reconRate)$y
axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=1, las=1, lwd=0.5, cex.axis=0.8, tck=-0.05, mgp=c(1,0.45,0))

dev.off()


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

radial.text2 <- function (r, theta, labels, cex = 1, col = "black", font = 1, 
    ...) 
{
    n <- length(labels)
    col <- rep(col, length.out = n)
    x <- r * cos(theta)
    y <- r * sin(theta)
    srt <- theta/(2 * pi) * 360 + 95
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
group.label.tip.rad4 <- function (obj, lab, col.bar, col.lab, lwd = 1, offset.bar = 0, 
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
        radial.text2(r.lab, tm, str, col = col.lab, font = font, 
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



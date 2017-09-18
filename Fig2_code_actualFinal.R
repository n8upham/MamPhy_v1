#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - MamPhy -- Upham et al. 2017
###
# Figure 2 - lineages and rates through time vs. K-Pg event
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Summarize the LTT data for various clades
###





# Run CoMET on full Mammalia, 10 trees
###



# Run TreePar on full Mammalia, 10 trees
###



# 
###





########
# Plot -- Parts A and C -- LTT of mammal higher taxa and RTT of Mammalia, with BAMM, CoMET, and TreePar rate shifts
########

pdf(file=paste("LTT_95polygon_higherTaxLTT_globalRTT_3part_no2Ma_wBAMM-CoMET-TreePar.pdf",sep=""), onefile=TRUE, width=4, height=8)

xLim<-c(-110,25)
#quartz(width=5,height=8)
layout(matrix(c(1:4), 4, 1, byrow = TRUE), widths=rep(1,4),heights=c(0.25,0.25,0.25,0.25))
op <- par(#mfrow = c(3,1),
          oma = c(5,4,2,0) + 0.1, #c(bottom, left, top, right)
          mar = c(0.5,0,1,1) + 0.1, xpd=TRUE)

LegRCex<-0.9 

# part 1
	# blank set up
x<-ltt95_Placentalia
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",yaxt="n",xaxt="n", bty="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
axis(side=1,labels=FALSE, at=c(0,-20,-40,-60,-80,-100), cex=0.8)
axis(side=2,labels=TRUE, at=c(1,5,50,500,5000), cex=0.8, line=0.25)
text(x=-115,y=8000,labels="(a) Higher taxa",cex=1.2,font=2,pos=4)
abline(v=-66,lwd=2,lty=2,col="dark grey")
mtext(side = 2, text = "(log) lineages", line = 3, cex = 1)

	# Placental and marsupial LTTs
	x<-ltt95_Placentalia
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.7), border = NA) #awesome.
	lines(x[,3],x[,1],lty=1,lwd=1,type="l")
	text(x=-22,y=nrow(ltt95_Placentalia[,])-500,labels="crown\nPlacentalia",font=2, cex=1,offset=0,pos=1)#pos=4)
#bquote(atop("crown","Placentalia"))
	x<-ltt95_Marsupalia
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=grey(0.7), border = NA) #awesome.
	lines(x[,3],x[,1],lty=1,lwd=1,type="l")
	#text(x=0,y=nrow(ltt95_Marsupalia[,])-200,labels="crown\nMarsupialia",font=2, cex=1,offset=0,,pos=1)#pos=4)
	text(x=-13,y=7,labels="crown\nMarsupialia",font=2, cex=1,offset=0,,pos=1)#pos=4)


par(new=TRUE)
colsSuper<-c("steelblue1","royalblue","royalblue4","midnightblue") #"darkgoldenrod1","darkgoldenrod3","darkgoldenrod4","saddlebrown") #c("steelblue1","steelblue2","steelblue3","steelblue4") #rich.colors(4)
lwdBar<-2
cexSuper<-0.9
plotCI(add=TRUE, gap=0, x=rev(-superord_Atlan$means),y=(exp(seq(log(1000),log(4500),length.out=6))), ui=rev(-superord_Atlan$hpd95_mins),li=rev(-superord_Atlan$hpd95_maxs), cex=1, err="x",sfrac=0, col=colsSuper[1],lwd=lwdBar, pch=16)
plotCI(add=TRUE, gap=0, x=rev(-superord_Bor$means),y=(400), ui=rev(-superord_Bor$hpd95_mins),li=rev(-superord_Bor$hpd95_maxs), cex=1, err="x",sfrac=0, col=colsSuper[2],lwd=lwdBar, pch=16)
plotCI(add=TRUE, gap=0, x=rev(-superord_Euarch$means),y=(exp(seq(log(200),log(100),length.out=3))), ui=rev(-superord_Euarch$hpd95_mins),li=rev(-superord_Euarch$hpd95_maxs), cex=1, err="x",sfrac=0, col=colsSuper[3],lwd=lwdBar, pch=16)
plotCI(add=TRUE, gap=0, x=rev(-superord_Lauras$means),y=(exp(seq(log(50),log(15),length.out=5))), ui=rev(-superord_Lauras$hpd95_mins),li=rev(-superord_Lauras$hpd95_maxs), cex=1, err="x",sfrac=0, col=colsSuper[4],lwd=lwdBar, pch=16)
	text(x=-60,y=800,labels=paste("Xenarthra\nand Afrotheria"),font=1, cex=cexSuper,offset=0,adj=c(0,0))#pos=4)
	text(x=-65,y=350,labels="Boreoeutheria",font=1, cex=cexSuper,offset=0,adj=c(0,0))#pos=4)
	text(x=-115,y=125,labels="Euarchontoglires",font=1, cex=cexSuper,offset=0,adj=c(0,0))#pos=4)
	text(x=-106,y=22,labels="Laurasiatheria",font=1, cex=cexSuper,offset=0,adj=c(0,0))#pos=4)


# part 2
	# BAMM RTT of global mammals
#library(BAMMtools)
ratetypes <- c("netdiv","speciation","extinction")
intervalCol <- grey(0.5) #viridis(3) #c("blue","green","darkyellow")
#avgCol <- grey(0,alpha=1) #viridis(3) #c("blue","green","darkyellow")#rep("red",3)
#lineTypes <- c(1,5,4)
lineTypes <- c(1,1,1)
avgCols <- c(grey(0.5),"black","white")

cexRLab<-0.65
cex.lab = 1

for(z in c(2,3,1)){
rmat<-join10runs_ALL
#rmat<-rateMatrices[[1]]

rates <- list((rmat$lambda - rmat$mu),rmat$lambda, rmat$mu)

rate<-rates[[z]][,1:202] # exclude last 2 Ma

ratelabel<-"rate"
ratetype<-ratetypes[z]
lineType<-lineTypes[z]
avgCol<-avgCols[z]

intervals = seq(from = 0, to = 1, by = 0.01)
smooth = FALSE 
opacity = 0.01
plot = TRUE
cex.axis = 1
lwd = 3
xline = 3.5
yline = 3.5
mar = c(6, 6, 1, 1)
xticks = NULL 
yticks = NULL
xlim = "auto"
ylim = "auto"
add = FALSE
axis.labels = TRUE

    maxTime <- max(rmat$times)
    nanCol <- apply(rate, 2, function(x) any(is.nan(x)))
    rate <- rate[, which(nanCol == FALSE)]
    #rmat$times <- rmat$times[which(nanCol == FALSE)]
    rmat$times <- (max(rmat$times) - rmat$times)[1:202] # exclude last 2 Ma

        mm <- apply(rate, 2, quantile, intervals)
        poly <- list()
        q1 <- 1
        q2 <- nrow(mm)
        repeat {
            if (q1 >= q2) {
                break
            }
            a <- as.data.frame(cbind(rmat$times, mm[q1, ]))
            b <- as.data.frame(cbind(rmat$times, mm[q2, ]))
            b <- b[rev(rownames(b)), ]
            colnames(a) <- colnames(b) <- c("x", "y")
            poly[[q1]] <- rbind(a, b)
            q1 <- q1 + 1
            q2 <- q2 - 1
        }

    harmMean<-function(x){
    	res<-(1/mean(1/x))
    	return(res)
    }
	avg <- apply(rate,2,median)#colMeans(rate)
#	offsetBamm<- 0.4
  offsetBamm<- 0

if(z == 2){
       	#plot.new()
        #par(mar = mar)          
 		xMin <- -110 #maxTime
        xMax <- 25
        yMin <- 0 #min(poly[[1]][, 2])
        yMax <- 0.25+offsetBamm #max(poly[[1]][, 2])
        if (yMin >= 0) { yMin <- 0 }
        #plot.window(xlim = xLim, ylim = c(yMin, yMax))

		plot(x = -rmat$times, y = (avg+offsetBamm),lty=2,lwd=1,xlim=c(xMin, xMax),ylim=c(yMin, yMax),type="l", col="white",main=NULL,cex.axis=0.8,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",yaxt="n",xaxt="n",bty="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
        
        axis(at = seq(0,0.2,by=0.05)+offsetBamm, labels=c(0,NA,0.1,NA,0.2),cex.axis = cex.axis, las = 1, side = 2, line=0.25)
        mtext(side = 2, text = "rate", line = 3, cex = cex.lab, adj=0.57+offsetBamm)

		abline(v=-66,lwd=2,lty=2,col="dark grey")
       }

        for (i in 1:length(poly)) {
            polygon(x = -poly[[i]][, 1], y = (poly[[i]][, 2]+offsetBamm), col = transparentColor(intervalCol, opacity), border = NA)
        }
    lines(x = -rmat$times, y = (avg+offsetBamm), lwd = lwd, col = avgCol, lty=lineType)
}
  # plot BAMM shifts... 
    pchs<-rep(21,length(shiftDivDat[,1]))
    #pchs[which(shiftDivDat$Indep==1)]<-16
    alphaCol<-0.7
    ptCols<-rep(rgb(1,0,0,alpha=alphaCol),length(shiftDivDat[,1]))
    ptCols[which(shiftDivDat$avgIncDec=="downShift")]<-rgb(0,0,1,alpha=alphaCol)
    ptCols[which(shiftDivDat$avgIncDec=="up-or-down")]<-grey(0.5,alpha=alphaCol)

#   yVar<-shiftDivDat$avgCladeDiv
    yVar<-shiftDivDat$avgFactor
    #yVar[which(shiftDivDat$avgIncDec=="downShift")]<- -1.7
    xVar<-(-shiftDivDat$mean)
    xLow<-(-shiftDivDat$lower)
    xHigh<-(-shiftDivDat$upper)

    bammShift_yMin<-1 # -2.5
    par(new=TRUE)
    plotCI(x=xVar, y=yVar+(3*offsetBamm),ui=xLow,li=xHigh, xlim=c(xMin, xMax), ylim=c(bammShift_yMin,6), xlab="", ylab="", sfrac=0,err="x", lwd=2,col="black",pt.bg=ptCols, scol="black",pch=pchs,font.lab=2,cex.axis=0.95,cex.lab=1, xaxt="n", yaxt="n",bty="n", cex=1.5)
    axis(side=4,line= -5,at=c(1,2,3,4)+(3*offsetBamm),labels=c(1,2,3,4))
    mtext(side = 4, text = "Shift factor", line = -3, cex = cexRLab, adj=0+offsetBamm)

    text(x=-115, y= 6,labels = "(b) Mammalia-wide rates", cex = 1.2,font=2,pos=4)

  #  rug(xVar[which(shiftDivDat$avgIncDec=="upShift")], col=rgb(1,0,0,alpha=alphaCol), side=1, lwd=2, line=-5-(22*offset))
  #  rug(xVar[which(shiftDivDat$avgIncDec=="downShift")], col=rgb(0,0,1,alpha=alphaCol), side=1, lwd=2, line=-5-(22*offset))
  #  rug(xVar[which(shiftDivDat$avgIncDec=="up-or-down")], col=grey(0.5,alpha=alphaCol), side=1, lwd=2, line=-5-(22*offset))

    legend(x=-115,y=5.5,legend=c("speciation","extinction","net diversification"), col=c(avgCols[2:3],avgCols[1]), lty=c(1,1,1), lwd=lwd, cex=0.9, seg.len=4, bty="o",bg=grey(0.85)) 

    legend(x=3,y=6,legend=c("up","down","up/down"), pt.bg=c(rgb(1,0,0,alpha=alphaCol),rgb(0,0,1,alpha=alphaCol),grey(0.5,alpha=alphaCol)), pch=21, cex=LegRCex, pt.cex=1.5, seg.len=4, bty="n", title="BAMM\nper-lineage\nshifts") 

#   legend(x=-115,y=offset,legend=c("speciation shift","extinction shift"), fill=CoMET_cols, border="black", cex=1, bty="n") 
#   text(x=-115,y=0.25+offset,labels="(c) Mammalia-wide rates",cex=1.2,font=2,pos=4)

#text(x=-115,y=0.22+offset,labels="Rates through time",cex=1,font=2,pos=4)
#text(x=-115,y=-0.03+offset,labels="Tree-wide rate shifts",cex=1,font=2,pos=4)

  #CoMET RTT plots... 
output<-MAMPHY_10tree_CoMET_1Ma_7shiftnoME

ratetypes <- c("netdiv","speciation","extinction")
intervalCol <- grey(0.2) #viridis(3) #c("blue","green","darkyellow")
#avgCol <- grey(0,alpha=1) #viridis(3) #c("blue","green","darkyellow")#rep("red",3)
#lineTypes <- c(1,5,4)
lineTypes <- c(1,1,1)
avgCols <- c(grey(0.5),"black","white")

for(z in c(2,3,1)){
rates <-list((output[["net-diversification rates"]]),(output[["speciation rates"]]),(output[["extinction rates"]]))

rate<-rates[[z]][,1:202] # exclude last 2 Ma
colnames(rate)<-203:2

ratelabel<-"rate"
ratetype<-ratetypes[z]
lineType<-lineTypes[z]
avgCol<-avgCols[z]

    maxTime <- max(output[["intervals"]])
    nanCol <- apply(rate, 2, function(x) any(is.nan(x)))
    rate <- rate[, which(nanCol == FALSE)]
    #plotTimes <- output[["intervals"]][which(nanCol == FALSE)]
    plotTimes <- output[["intervals"]][1:202] # exclude last 2 Ma

        #mm <- apply(rate, 2, quantile, prob = c(0.025, 0.975)) #intervals)

        mm <- apply(rate, 2, quantile, intervals)
        poly <- list()
        q1 <- 1
        q2 <- nrow(mm)
        repeat {
            if (q1 >= q2) {
                break
            }
            a <- as.data.frame(cbind(rmat$times, mm[q1, ]))
            b <- as.data.frame(cbind(rmat$times, mm[q2, ]))
            b <- b[rev(rownames(b)), ]
            colnames(a) <- colnames(b) <- c("x", "y")
            poly[[q1]] <- rbind(a, b)
            q1 <- q1 + 1
            q2 <- q2 - 1
        }

    harmMean<-function(x){
      res<-(1/mean(1/x))
      return(res)
    }
  avg <- apply(rate,2,median)#colMeans(rate)
#  offsetCoMET<- 0.35
  offsetCoMET<- 0.1

if(z == 2){
        #plot.new()
        #par(mar = mar)          
    xMin <- -110 #maxTime
        xMax <- 25
        yMin <- 0 #min(poly[[1]][, 2])

CoMETrate_yMax<-0.38 #0.29
        yMax <- CoMETrate_yMax+offsetCoMET #max(poly[[1]][, 2])
        if (yMin >= 0) { yMin <- 0 }
        #plot.window(xlim = xLim, ylim = c(yMin, yMax))

#    par(new=TRUE)
    plot(x = -plotTimes, y = (avg+offsetCoMET),lty=2,lwd=1,xlim=c(xMin, xMax),ylim=c(yMin, yMax),type="l", col="white",main=NULL,cex.axis=0.8,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",yaxt="n",xaxt="n",bty="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
        
        axis(side=2, at = seq(0,0.2,by=0.05)+offsetCoMET, labels=c(0,NA,0.1,NA,0.2),cex.axis = cex.axis, las = 1, line=0.25)
        mtext(side = 2, text = "rate", line = 3, cex = cex.lab, adj=0+offsetCoMET)

    abline(v=-66,lwd=2,lty=2,col="dark grey")
       }

        for (i in 1:length(poly)) {
            polygon(x = -poly[[i]][, 1], y = (poly[[i]][, 2]+offsetCoMET), col = transparentColor(intervalCol, opacity), border = NA)
        }

  #    polygon(x = -c(0:ncol(mm), ncol(mm):0), 
  #              y = c(c(mm[1, 1], mm[1, 
  #                ]), rev(c(mm[2, 1], mm[2, 
  #                ]))), border = NA, col = transparentColor(intervalCol, 0.5))

    lines(x = -plotTimes, y = (avg+offsetCoMET), lwd = lwd, col = avgCol, lty=lineType)
}


	#CoMET rate shifts...
CoMET_cols<-c(rgb(1,0,0,alpha=0.5),rgb(0,0,1,alpha=0.5))  # c("black",grey(0.7,alpha=0.6))

output<-MAMPHY_10tree_CoMET_1Ma_7shiftnoME
type<-"speciation shift times"
criticalPP <- output[[grep(strsplit(type, " ")[[1]][1], grep("CriticalPosteriorProbabilities", names(output), value = TRUE), value = TRUE)]]

type<-"extinction shift times"
colBars<-CoMET_cols[2]

thisOutput <- output[[type]]
meanThisOutput <- colMeans(thisOutput)
meanThisOutput[202:203]<-NA # excludes last 2 Ma from the plot
par(new=TRUE, lwd=0.1)
barplot(rev(meanThisOutput), space = 0, col = colBars, border = "white", main = "", ylab = "", xlab = "", xaxt = "n", yaxt = "n", xlim=c(-xMin,-xMax),ylim = c(0,0.5))

type<-"speciation shift times"
colBars<-CoMET_cols[1]

thisOutput <- output[[type]]
meanThisOutput <- colMeans(thisOutput)
meanThisOutput[202:203]<-NA # excludes last 2 Ma from the plot
  par(new=TRUE, lwd=0.1)
  barplot(rev(meanThisOutput), space = 0, col = colBars, border = "white", main = "", ylab = "", xlab = "", xaxt = "n", yaxt = "n", xlim=c(-xMin,-xMax),ylim = c(0,0.5))

#axis(2, at = c(0,criticalPP[1],criticalPP[2]), labels = c(0,2,6), line=0.25)#, line= -4)
axis(4, at = c(0,criticalPP[1]), labels = c(0,2), line=-5) #0.25)#, line= -4)
mtext(side = 4, text = "2 ln BF", line = -3, cex = cexRLab, adj=0) 

legend(x=-3,y=0.4,legend=c("speciation","extinction"), pt.bg=c(rgb(1,0,0,alpha=alphaCol),rgb(0,0,1,alpha=alphaCol)), pch=22, cex=LegRCex, pt.cex=1.5, seg.len=4, bty="n", title="CoMET\ntree-wide\nshifts") 

#points(rev(meanThisOutput),type="s",main = "", ylab = "", xlab = "", xaxt = "n", yaxt = "n", xlim=c(-xMin,-xMax),ylim = c(0,1.5), bty="n")

#lines(x=c(115,0), y=c(criticalPP[1],criticalPP[1]), lty = 1, col=grey(0.5))
#lines(x=c(115,0), y=c(criticalPP[2],criticalPP[2]), lty = 1, col=grey(0.5))



 # TreePar RTT plots... 
TP_RTT_all<-list(TP_netdiv_10trees[,length(grid):1], TP_speciation_10trees[,length(grid):1], TP_extinction_10trees[,length(grid):1])

ratetypes <- c("netdiv","speciation","extinction")
intervalCol <- grey(0.2) #viridis(3) #c("blue","green","darkyellow")
#avgCol <- grey(0,alpha=1) #viridis(3) #c("blue","green","darkyellow")#rep("red",3)
#lineTypes <- c(1,5,4)
lineTypes <- c(1,1,1)
avgCols <- c(grey(0.5),"black","white")

for(z in c(2,3,1)){
rates <-TP_RTT_all

rate<-rates[[z]][,1:(length(rates[[z]][1,])-2)] # exclude last 2 Ma
#colnames(rate)<-203:2

ratelabel<-"rate"
ratetype<-ratetypes[z]
lineType<-lineTypes[z]
avgCol<-avgCols[z]

    maxTime <- as.numeric(colnames(rate)[1]) #max(output[["intervals"]])
    nanCol <- apply(rate, 2, function(x) any(is.nan(x)))
    rate <- rate[, which(nanCol == FALSE)]
    #plotTimes <- output[["intervals"]][which(nanCol == FALSE)]
    plotTimes <- as.numeric(colnames(rate)) # output[["intervals"]][1:202] # >> already excludes last 2 Ma

        #mm <- apply(rate, 2, quantile, prob = c(0.025, 0.975)) #intervals)

        mm <- apply(rate, 2, quantile, intervals)
        poly <- list()
        q1 <- 1
        q2 <- nrow(mm)
        repeat {
            if (q1 >= q2) {
                break
            }
            a <- as.data.frame(cbind(plotTimes, mm[q1, ]))
            b <- as.data.frame(cbind(plotTimes, mm[q2, ]))
            b <- b[rev(rownames(b)), ]
            colnames(a) <- colnames(b) <- c("x", "y")
            poly[[q1]] <- rbind(a, b)
            q1 <- q1 + 1
            q2 <- q2 - 1
        }

    harmMean<-function(x){
      res<-(1/mean(1/x))
      return(res)
    }
  avg <- apply(rate,2,median)#colMeans(rate)
#  offsetTreePar<- 0.15
  offsetTreePar<- 0.2

if(z == 2){
        #plot.new()
        #par(mar = mar)          
    xMin <- -110 #maxTime
        xMax <- 25
        yMin <- 0 #min(poly[[1]][, 2])
TreePar_yMax<-0.7
        yMax <- TreePar_yMax+offsetTreePar #max(poly[[1]][, 2])
        if (yMin >= 0) { yMin <- 0 }
        #plot.window(xlim = xLim, ylim = c(yMin, yMax))

#    par(new=TRUE)
    plot(x = -plotTimes, y = (avg+offsetTreePar),lty=2,lwd=1,xlim=c(xMin, xMax),ylim=c(yMin, yMax),type="l", col="white",main=NULL,cex.axis=0.8,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",yaxt="n",xaxt="n",bty="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
        
        axis(side=2, at = seq(0,0.2,by=0.05)+offsetTreePar, labels=c(0,NA,0.1,NA,0.2),cex.axis = cex.axis, las = 1, line=0.25)
        mtext(side = 2, text = "rate", line = 3, cex = cex.lab, adj=0+offsetTreePar)

    abline(v=-66,lwd=2,lty=2,col="dark grey")
    text(x=66,y=TreePar_yMax,labels="K-Pg\nboundary",cex=1,font=2,pos=4, col=grey(0.5))

       }

        for (i in 1:length(poly)) {
            polygon(x = -poly[[i]][, 1], y = (poly[[i]][, 2]+offsetTreePar), col = transparentColor(intervalCol, opacity), border = NA)
        }

  #    polygon(x = -c(0:ncol(mm), ncol(mm):0), 
  #              y = c(c(mm[1, 1], mm[1, 
  #                ]), rev(c(mm[2, 1], mm[2, 
  #                ]))), border = NA, col = transparentColor(intervalCol, 0.5))

    lines(x = -plotTimes, y = (avg+offsetTreePar), lwd = lwd, col = avgCol, lty=lineType)
}

    axis(side=1, at= rev(c(0,-20,-40,-60,-80,-100)), labels=rev(-c(0,-20,-40,-60,-80,-100)), cex=cex.axis, line=0.5)
    mtext(side = 1, text = "time before present (Ma)", line = 3.5, cex = cex.lab, adj=0.4)

  ## TreePar shift HISTOS
shiftUp<-TP_resFIN$shiftMa[which(TP_resFIN$dir=="upShift" & TP_resFIN$shiftMa < -xMin)]
shiftDown<-TP_resFIN$shiftMa[which(TP_resFIN$dir=="downShift" & TP_resFIN$shiftMa < -xMin)]

times<-seq(from=xMin,to=xMax,by=1)
histUp<-hist(-shiftUp, breaks=times, plot=FALSE)$counts
names(histUp)<-times[1:(length(times)-1)]
histDown<-hist(-shiftDown, breaks=times, plot=FALSE)$counts
names(histDown)<-times[1:(length(times)-1)]

par(new=TRUE)
barplot(histUp, space = 0, main="", col=rgb(1,0,0,alpha=0.5), border="white",lwd=0.5, ylim=c(0,10),yaxt="n",xaxt="n")
barplot(add=TRUE,histDown, space = 0,  col=rgb(0,0,1,alpha=0.5), border="white", lwd=0.5, xaxt="n", yaxt="n")

axis(side=4,at=seq(from=0,to=max(histUp),length.out=4),labels=c(0,NA,NA,max(histUp)/10), line=-5) # seq(from=0,to=max(c(histDown,histUp)),length.out=3)
mtext(side=4,text="Shifts per tree", line=-3, cex=cexRLab, adj=-0.2)

legend(x=113,y=7,legend=c("up","down"), pt.bg=c(rgb(1,0,0,alpha=alphaCol),rgb(0,0,1,alpha=alphaCol)), pch=22, cex=LegRCex, pt.cex=1.5, seg.len=4, bty="n", title="TreePar\ntree-wide\nshifts") 

# plot(density(-TP_resFIN$shiftMa), type="n",xlim=c(xMin,xMax), yaxt="n",xaxt="n", bty="n",main="",ylab="",xlab="")
# rug(-shiftDown-1, col=grey(0.3), lwd=2,line=0.5)
# rug(-shiftUp-1, col=rgb(1,0,0,alpha=0.5), lwd=2, line=1.5)
# lines(x=c(-66,-66),y=c(-2,1),col=grey(0.5),lty=2,lwd=3)
# axis(side=1,at=seq(from=-100,to=0,by=20),labels= seq(from=-100,to=0,by=20), line=1.5)


dev.off()






########
# Plot -- Parts B and D -- LTT of K-Pg orders and RTT of 5 most speciose, with BAMM shifts
########





########
# Plot -- Part E -- Tip DR distributions and ML model-fitting results
########



# 
###




#######
# Supplemental aspects of Fig 2 analysis
#######






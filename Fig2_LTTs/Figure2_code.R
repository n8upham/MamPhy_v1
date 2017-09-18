###
# Get the BAMM results back in order to plot net div through time...
# Overall and on a per-order basis...
#####

library(BAMMtools); library(coda); library(phytools); library(ape); library(phangorn)

folders<-c(1:10)
bbone = "NDexp" # "FBD" 
dir = "NDexp_" #"" #

rateMatrices_1Ma<-vector("list",length(folders))
for (i in 1:length(folders)){

#setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/grace_2.5_1week_",dir,folders[i],sep=""))
setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/grace_2.5_1week_NDexp_",folders[i],sep=""))

# load edata
tree <- ladderize(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample10_",folders[i],".tre",sep="")))
edata <- getEventData(tree, eventdata = paste("mamPhy_",bbone,"_event_data.txt",sep=""), burnin=0.33, nsamples=1000)

root<-max(branching.times(tree))

rateMatrices_1Ma[[i]]<-getRateThroughTimeMatrix(edata,nslices=(root/1))

}
setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/",sep=""))
save(rateMatrices_1Ma,file="rateThroughTimeMatrix_10trees_1Ma.Rda")

######
# Now do the same for ** mammal subclades.... **
# Get the BAMM results back in order to plot net div through time...
#########
library(BAMMtools); library(coda); library(phytools); library(ape); library(phangorn)
#setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/",sep=""))
setwd(paste("/mnt/data/personal/nateu/Nate_Backup/MamPhy_BDvr_Fixed-n-Free_fullPosteriors/"))

bbone = "NDexp" # "FBD" 

# tip categories
cladesDR<-read.table(paste("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_",bbone,"_DRstats_DRtreeLABELS.txt",sep=""), header=TRUE)
head(cladesDR)

ordNames<-names(which(table(cladesDR$ord) > 25))[c(1:4,8:11)] # [8:11] # excluding "DASYUROMORPHIA"  "DIDELPHIMORPHIA" "DIPROTODONTIA" bc crowns are ~ 20, 21, 30.1 Ma < 30 Ma needed for 7 splits of 5 Ma...
ordTips<-list()
for(j in 1:length(ordNames)){
	ordTips[[j]]<-cladesDR[which(cladesDR$ord==ordNames[j]),"tiplabel"]
}
higherNames<-names(which(table(cladesDR$higher) > 25))
higherTips<-list()
for(j in 1:length(higherNames)){
	higherTips[[j]]<-cladesDR[which(cladesDR$higher==higherNames[j]),"tiplabel"]
}
setNames<-c(ordNames,higherNames)
setTips<-c(ordTips,higherTips)

folders<-c(1:10)
dir = "NDexp_" #"" #

# go through trees...
# loading edata, subsetting to clade, getting that rate through time matrix
for (i in 1:length(folders)){
setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/grace_2.5_1week_",dir,folders[i],sep=""))
#setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/grace_2.5_1week_NDexp_",folders[i],sep=""))

	# load edata
	tree <- ladderize(read.tree(paste("MamPhy_fullPosterior_BDvr_pcsFIXED_",bbone,"_sample10_",folders[i],".tre",sep="")))
	edata <- getEventData(tree, eventdata = paste("mamPhy_",bbone,"_event_data.txt",sep=""), burnin=0.33, nsamples=1000)

	for(j in 1:length(setNames)){
		subtree_edata<-subtreeBAMM(edata,tips=as.vector(setTips[[j]]))

		toDrop<-setdiff(tree$tip.label,as.vector(setTips[[j]]))
		subtree<-drop.tip(tree,toDrop)
		root<-max(branching.times(subtree))

		rateMat<-getRateThroughTimeMatrix(subtree_edata,nslices=(root/1))

		assign(paste("rateMatrices_1Ma_",setNames[j],"_",i,sep=""),inherits=TRUE,value=rateMat)
	}
}

for(j in 1:length(setNames)){ 
	RM_list<-vector("list",length(folders))
	for(i in 1:length(folders)){
	RM_list[[i]]<-get(paste("rateMatrices_1Ma_",setNames[j],"_",i,sep=""))
	}
	assign(paste("rateMatrices_1Ma_",setNames[j],"_ALL10",sep=""),RM_list)
}

setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/",sep=""))
save.image(file="RTT_workspace_ords-and-higher_bins1Ma_10trees.Rda")







######################

#########

# Load in CoMET results...
library(TESS); library(plotrix)
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMet_inTESS/")
load(file="MAMPHY_10tree_CoMET_summary_5Ma.Rda")
source("tess.plot.output2.R")
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")

##
xMin<- -155
xMax<- 30
labelCex<-0.7

pdf(file="divTime_24BAMMnodes_NDexp_plot95pCIs_byDivTime_Phylo_wExtinctions.pdf", onefile=TRUE, width=6,height=8)
#pdf(file="divTime_24BAMMnodes_NDexp_plot95pCIs_byDivTime_Phylo.pdf", onefile=TRUE, width=6,height=5)
op <- par(oma = c(1,1,1,1) + 0.1,		#c(bottom, left, top, right)
          mar = c(4,4,1,1) + 0.1, lwd=2)
layout(matrix(1:2,2,1), heights=c(0.7,0.3), widths=c(1,1))

# BAMM shifts
plotCI(x=rev(-means),y=1:length(shiftDivDat[,1]), ui=rev(-low95s),li=rev(-up95s), cex=rev(bgsizes), err="x",sfrac=0, col=NA,xlim=c(xMin,xMax),bty="n",xaxt="n",yaxt="n",xlab="Time before present (Ma)", ylab="",lwd=2, pch=NA)
for(i in nameAtY){
    segments(x0=xMin,y0=i,x1=0,y1=i,lty=1, lwd=0.5, col=grey(0.9, alpha=1))
}
	segments(x0=-66.043,y0=0,x1=-66.043,y1=32,lty=2, lwd=2, col=grey(0.5, alpha=1))

plotCI(add=TRUE,x=rev(-means),y=1:length(shiftDivDat[,1]), ui=rev(-low95s),li=rev(-up95s), cex=rev(bgsizes), err="x",sfrac=0, col=rev(bgcols),xlim=c(-110,80),bty="n",xaxt="n",yaxt="n",xlab="Time before present (Ma)", ylab="",lwd=2, pch=20)
axis(side=1,at=seq(-150,0,by=50),labels=c(150,100,50,0), cex.axis=0.85)
text(x=rev(-means)+1,y=(1:length(shiftDivDat[,1]))+0.2,labels=rev(shiftDivDat$ID), cex=0.7, font=2, adj=c(0,0))
text(x=0,y=(1:length(shiftDivDat[,1]))-0.4,labels=rev(numTaxa), cex=labelCex, font=2, adj=c(0,0))
#text(x=rep(0,length(agesAllClades[,1])),y=(1:length(agesAllClades[,1]))+0.2,labels=rev(agesAllClades$ID_label), cex=0.7, font=2, adj=c(0,0))
#par(xpd=TRUE)
text(x=rep(0,length(cladeNames)),y=nameAtY+0.30,labels=cladeNames,cex=labelCex, font=1, adj=0)#, adj=1)

# CoMET shifts
output<-MAMPHY_10tree_CoMET_5Ma
	
	treeAge <- output[["tree"]]$root.age
	treeAge
    numIntervals <- length(output$intervals) - 1
    numIntervals
    plotAt <- 0:numIntervals
    plotAt
    intervalTimes<-seq(0,205,by=5)
    intervalTimes
    intervalSize <- max(intervalTimes)/numIntervals
    intervalSize
 
    type<-"mass extinction times"
    thisOutput <- output[[type]]
    colnames(thisOutput)<-rev(intervalTimes[1:41])

    meanThisOutput <- colMeans(thisOutput)
    toPlot<-round(meanThisOutput,3)[11:41] # plots from 150-0 Mya
   	
   	labels <- pretty(c(0, names(toPlot[1])))
    labels
    labelsAt <- c(41,31,21,11) -10
    labelsAt

    criticalPP<- output[[grep(strsplit(type, " ")[[1]][1], grep("CriticalPosteriorProbabilities", names(output), value = TRUE), value = TRUE)]]

barplot(toPlot, space = 0, xaxt = "n", col = "#4DAF4A", 
        border = "#4DAF4A", main = "", ylab = "posterior probability", 
        xlab = "million years ago", yaxt = "s", pch = 19, ylim=c(0,0.4), xlim=c(0, (185/5)))
axis(1, at = labelsAt, labels = labels)
axis(4, at = c(0,criticalPP[1:2]), labels = c(0,(2 * log(output$criticalBayesFactors[1:2]))), las = 1, tick = TRUE, line = -4)
lines(x=c(max(labelsAt),min(labelsAt)),y=c(criticalPP[1],criticalPP[1]),lty = 2,col="grey")
lines(x=c(max(labelsAt),min(labelsAt)),y=c(criticalPP[2],criticalPP[2]),lty = 2,col="grey")
mtext("Bayes factors",side=4, line=-2)
abline(v = (155/5) - (66.02/5), lty = 2, lwd=3, col="dark grey")

#tess.plot.output2(MAMPHY_10tree_CoMET_5Ma, fig.types = c("mass extinction times"), las=2, plot.tree=FALSE,xlim=c(-xMax,-xMin))

dev.off()


## Calc the MEAN and 95% CI for the inferred extinction divTimes...
#####

    type<-"mass extinction times"
    thisOutput <- output[[type]]
    colnames(thisOutput)<-rev(intervalTimes[1:41])

    meanThisOutput <- colMeans(thisOutput)

    sumThisOutput <- colSums(thisOutput)[11:41] # plots from 150-0 Mya
   
    allVals<-vector("list",length(sumThisOutput))
    for(j in 1:length(sumThisOutput)){
    	val<-as.numeric(names(sumThisOutput[j]))
    	allVals[[j]]<-rep(val,sumThisOutput[[j]])
    }
    freqDist<-do.call(c,allVals)

    mean(freqDist)
    #[1] 99.83792
	quantile(freqDist, c(0.025,0.5,0.975))
	# 2.5%  50%  97.5% 
   	# 80    95   135  

    toPlot<-round(meanThisOutput,3)[11:41] # plots from 150-0 Mya

######



type<-"speciation shift times"
    thisOutput <- output[[type]]
    colnames(thisOutput)<-rev(intervalTimes[1:41])
	meanThisOutput <- colSums(thisOutput)
    toPlot<-round(meanThisOutput/length(thisOutput[,1]),3)[11:41] # plots from 150-0 Mya
   	
barplot(add=TRUE, toPlot, space = 0, xaxt = "n", col = "#4DAF4A", 
        border = "#4DAF4A", main = type, ylab = "posterior probability", 
        xlab = "million years ago", yaxt = "s", pch = 19, ylim=c(0,0.5), xlim=c(0, (180/5)))

            



#####################






#####
# Load BAMM rate through time data back in...

setwd("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/_NDexp_10trees")
#save(rateMatrices_5Ma,file="rateThroughTimeMatrix_10trees_5Ma.Rda")
#load(file="rateThroughTimeMatrix_10trees_5Ma.Rda")
load(file="rateThroughTimeMatrix_10trees_1Ma.Rda")


# I want to join the 10 tree files by column... starting at column 1 (or is it 36, etc?)
allCols<-c(); 
for(j in c(1:7,9:10)){
    allCols[j]<-dim(rateMatrices_1Ma[[j]][["lambda"]])[2]
}
allCols # this is the number of -- 1 Ma -- 5 Ma bins per tree
maxBin<-max(na.omit(allCols)) # using tree 4, 204 intervals-- 41 intervals

rateMatrices<-rateMatrices_1Ma

for(j in c(1:7,9:10)){
  colToAdd<-maxBin-allCols[j] #dim(output_list_ALL[[j]][[k]])[2]

  for(k in c("lambda","mu")){ # the variables based on bins where NEED to even them out!

      if(colToAdd==0){ addCol<-NULL } else if(colToAdd==1) { 
      rowToAdd<-dim(rateMatrices[[j]][[k]])[1]
      addCol<-rep(rateMatrices[[j]][[k]][1,1],rowToAdd)
    } else {
      rowToAdd<-dim(rateMatrices[[j]][[k]])[1]
      addCol<-rep(rateMatrices[[j]][[k]][1,1],rowToAdd)
      for(q in 1:(colToAdd-1)){
        addCol<-cbind(addCol,rep(rateMatrices[[j]][[k]][1,1],rowToAdd))
      }
    }
      rateMatrices[[j]][[k]]<-cbind(addCol,rateMatrices[[j]][[k]])
      colnames(rateMatrices[[j]][[k]])<-203:0 #5*(40:0)
    }
}

allCols2<-c(); 
for(j in c(1:7,9:10)){
#    allCols2[j]<-dim(output_list_ALL[[j]][[5]])[2]
    allCols2[j]<-length(rateMatrices[[j]][["mu"]][1,])
}
allCols2 # this is the number of 5 Ma bins per tree
maxBin2<-max(na.omit(allCols2))


intToUse<-which(allCols==maxBin)

join10runs_ALL<-rateMatrices[[intToUse]]
for(key in c("lambda","mu")){ # the variables based on bins where NEED to even them out!
	join10runs_ALL[[key]]<-rbind(rateMatrices[[1]][[key]], rateMatrices[[2]][[key]], rateMatrices[[3]][[key]], rateMatrices[[4]][[key]], 
        rateMatrices[[5]][[key]], rateMatrices[[6]][[key]], rateMatrices[[7]][[key]], 
        #rateMatrices[[8]][[key]], 
        rateMatrices[[9]][[key]], 
        rateMatrices[[10]][[key]])
}

# Now PLOT it...
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/_NDexp_10trees")
library(viridis)

ratetypes <- c("netdiv","speciation","extinction")
intervalCols <- viridis(3) #c("blue","green","darkyellow")
avgCols <- viridis(3) #c("blue","green","darkyellow")#rep("red",3)


pdf(file=paste("rateThroughTime_allMamm_10trees_combo_bins1Ma_3rates_KPg.pdf",sep=""), onefile=TRUE, width=6, height=4)

for(z in length(ratetypes):1){
rmat<-join10runs_ALL
#rmat<-rateMatrices[[1]]

rates <- list((rmat$lambda - rmat$mu),rmat$lambda, rmat$mu)

rate<-rates[[z]]
ratelabel<-"rate"
ratetype<-ratetypes[z]
intervalCol<-intervalCols[z]
avgCol<-avgCols[z]

intervals = seq(from = 0, to = 1, by = 0.01)
smooth = FALSE 
opacity = 0.01
plot = TRUE
cex.axis = 1
cex.lab = 1.3
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
    rmat$times <- rmat$times[which(nanCol == FALSE)]
    rmat$times <- max(rmat$times) - rmat$times

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

	if(z == 3){
            plot.new()
            par(mar = mar)          
     		xMin <- 150 #maxTime
            xMax <- 0
            yMin <- min(poly[[1]][, 2])
            yMax <- 0.3 #max(poly[[1]][, 2])
            if (yMin >= 0) { yMin <- 0 }
            plot.window(xlim = c(xMin, xMax), ylim = c(yMin, yMax))

			axis(at = c(round(1.2 * xMin), axTicks(1)), cex.axis = cex.axis, side = 1)
            axis(at = seq(0,0.3,by=0.05), labels=c(0,NA,0.1,NA,0.2,NA,0.3),cex.axis = cex.axis, las = 1, side = 2)
            
            mtext(side = 1, text = "time before present", line = xline, cex = cex.lab)
            mtext(side = 2, text = ratelabel, line = yline, cex = cex.lab)

            lines(x = c(66,66), y = c(-0.1,0.4), lwd = 3, lty=2, col = "grey")
           }

            for (i in 1:length(poly)) {
                polygon(x = poly[[i]][, 1], y = poly[[i]][, 2], col = transparentColor(intervalCol, opacity), border = NA)
            }
        lines(x = rmat$time, y = avg, lwd = lwd, col = avgCol)
    }

legend(x=150,y=0.29,legend=c("net diversification","speciation","extinction"), col=viridis(3), lwd=lwd)

dev.off()   

###
# Load back in the TreePar only...

setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/TreePar_results",sep=""))
num_shifts=7; grid=1
numTrees<-10

setNames<-c("globalTree","RODENTIA","CHIROPTERA","EULIPOTYPHLA","PRIMATES","CETARTIODACTYLA","AFROSORICIDA",
        "Euarchontoglires","Laurasiatheria", "Marsupialia","Afroth.","Xenart.","CARNIVORA","LAGOMORPHA")

load(file=paste("TreePar_shifts",num_shifts,"_grid",grid,"Ma_clade--ALL--RES_bestModelShifts.Rda",sep=""))
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")

# BINNING the treePar results...
q<-1
  TP_resFIN<-get(paste("TP_resFIN_",setNames[q],sep=""))
  TP_bestModels<-get(paste("TP_bestModels_",setNames[q],sep=""))

  shifttimes<-TP_bestModels[,paste("ShiftTime",1:7,sep="")]
  netdiv<-TP_bestModels[,paste("DivRate",1:8,sep="")]
  turnover<-TP_bestModels[,paste("Turnover",1:8,sep="")]
  speciation<-(netdiv/(1-turnover))
    colnames(speciation)<-paste("Speciation",1:8,sep="")
  extinction<-((netdiv*turnover)/(1-turnover))
    colnames(extinction)<-paste("Extinction",1:8,sep="")

yVarNames<-c("netdiv","speciation","extinction")
for(z in 1:length(yVarNames)){
yVar<-get(yVarNames[z])
  varByGrid_ALL<-vector("list",length(shifttimes[,1]))

  for(i in 1:length(shifttimes[,1])){
    xVar<-unlist(c(0,shifttimes[i,][which(shifttimes[i,]<116)]))
    dat<-na.omit(cbind.data.frame(x=xVar,y=unlist(yVar[i,][1:length(xVar)])))
      rownames(dat)<-paste("t",1:length(dat[,1]),sep="")

    grid<-(seq(from=0,to=115,by=1))
    varByGrid<-rep(NA,length(grid))
    names(varByGrid)<-grid
    varByGrid[116]<-dat$y[length(dat$y)]

      for(k in 1:length(dat$x)){
        for(j in 1:(length(grid)-1)){
          if(dat$x[k] < grid[j+1]){ varByGrid[j]<-dat$y[k] } else if(dat$x[k] > grid[j+1] & dat$x[k] < grid[j+2]){ varByGrid[j]<-dat$y[k] } else { next }
        }
      }
      varByGrid_ALL[[i]]<-varByGrid
  }
  assign(paste("TP_",yVarNames[z],"_10trees",sep=""),do.call(rbind,varByGrid_ALL))
}

TP_RTT_all<-list(TP_netdiv_10trees[,length(grid):1], TP_speciation_10trees[,length(grid):1], TP_extinction_10trees[,length(grid):1])

# HISTOS for the treePar results...
xMin<- -110
xMax<-25
shiftUp<-TP_resFIN$shiftMa[which(TP_resFIN$dir=="upShift" & TP_resFIN$shiftMa < -xMin)]
shiftDown<-TP_resFIN$shiftMa[which(TP_resFIN$dir=="downShift" & TP_resFIN$shiftMa < -xMin)]

times<-seq(from=xMin,to=xMax,by=1)
histUp<-hist(-shiftUp, breaks=times, plot=FALSE)$counts
names(histUp)<-times[1:(length(times)-1)]
histDown<-hist(-shiftDown, breaks=times, plot=FALSE)$counts
names(histDown)<-times[1:(length(times)-1)]


#     # plot raw 95 CI as points
#          cols<-c(grey(0),"red","blue")
#        pdf(file=paste("TreePar_shifts",num_shifts,"_grid",grid,"Ma_globalMamm_only_bestModelShifts_PLOT_asBinned1Ma_med95CI.pdf",sep=""),onefile=TRUE)
#        z=1
#        plot(type="l",x= -as.numeric(colnames(get(paste("TP_",yVarNames[z],"_10trees",sep="")) )), y=apply(get(paste("TP_",yVarNames[z],"_10trees",sep="")) ,2,median),col=cols[z], ylab="rate",xlab="time",main=yVarNames[z],ylim=c(0,0.5))
#        for(z in 1:length(yVarNames)){
#          points(type="l",x= -as.numeric(colnames(get(paste("TP_",yVarNames[z],"_10trees",sep="")) )), y=apply(get(paste("TP_",yVarNames[z],"_10trees",sep="")) ,2,median),col=cols[z], ylab="rate",xlab="time",main=yVarNames[z],ylim=c(0,0.5))
#          points(x= -as.numeric(colnames(get(paste("TP_",yVarNames[z],"_10trees",sep="")) )), y=apply(get(paste("TP_",yVarNames[z],"_10trees",sep="")) ,2,quantile,c(0.025)),col=cols[z])
#          points(x= -as.numeric(colnames(get(paste("TP_",yVarNames[z],"_10trees",sep="")) )), y=apply(get(paste("TP_",yVarNames[z],"_10trees",sep="")) ,2,quantile,c(0.975)),col=cols[z])
#        }
#        dev.off()

    # TEST
#     # Plot just the TreePar res...
        pdf(file=paste("TreePar_shifts",num_shifts,"_grid",grid,"Ma_globalMamm_only_bestModelShifts_PLOT_upDownRug_new_1Ma.pdf",sep=""),onefile=TRUE)

        q<-1
          TP_resFIN<-get(paste("TP_resFIN_",setNames[q],sep=""))
          TP_bestModels<-get(paste("TP_bestModels_",setNames[q],sep=""))

          xMin<- -110
          xMax<-0
          shiftUp<-TP_resFIN$shiftMa[which(TP_resFIN$dir=="upShift" & TP_resFIN$shiftMa < -xMin)]
          shiftDown<-TP_resFIN$shiftMa[which(TP_resFIN$dir=="downShift" & TP_resFIN$shiftMa < -xMin)]
          
          times<-seq(from=xMin,to=xMax,by=1)
          histUp<-hist(-shiftUp, breaks=times, plot=FALSE)$counts
          names(histUp)<-times[1:(length(times)-1)]
          histDown<-hist(-shiftDown, breaks=times, plot=FALSE)$counts
          names(histDown)<-times[1:(length(times)-1)]

          plot(density(-TP_resFIN$shiftMa), type="n",xlim=c(xMin,xMax), yaxt="n",xaxt="n", bty="n",main="",ylab="",xlab="")
          
          rug(-shiftDown-1, col=grey(0.3), lwd=2,line=0.5)
          rug(-shiftUp-1, col=rgb(1,0,0,alpha=0.5), lwd=2, line=1.5)
          lines(x=c(-66,-66),y=c(-2,1),col=grey(0.5),lty=2,lwd=3)
          axis(side=1,at=seq(from=-100,to=0,by=20),labels= seq(from=-100,to=0,by=20), line=1.5)

          par(new=TRUE)
          barplot(histUp, space = 0, main="", col=rgb(1,0,0,alpha=0.5), border="white",lwd=0.5, xlim=c(0,length(times)), ylim=c(0,max(histUp)),yaxt="n",xaxt="n")
          barplot(add=TRUE,histDown, space = 0,  col=grey(0.3), border="white", lwd=0.5, xaxt="n", yaxt="n")

          axis(side=4,at=seq(from=0,to=max(histUp),length.out=4),labels=seq(from=0,to=max(histUp),length.out=4)/10, line=0) # seq(from=0,to=max(c(histDown,histUp)),length.out=3)
          mtext(side=4,text="Freq shifts per tree")

          par(new=TRUE)
          plot(x=c(0,-TP_bestModels[1,paste("ShiftTime",1:7,sep="")]),y=TP_bestModels[1,paste("DivRate",1:8,sep="")], type="s", xlim=c(xMin,xMax), ylim=c(0,0.5),xaxt="n", bty="n",main="",ylab="",xlab="", col="black", lwd=2, lty=1)
          for(j in 2:length(TP_bestModels[,1])){
            points(x=c(0,-TP_bestModels[j,paste("ShiftTime",1:7,sep="")]),y=TP_bestModels[j,paste("DivRate",1:8,sep="")], type="s", xlim=c(xMin,xMax), xaxt="n", bty="n",main="",ylab="",xlab="", col="black", lwd=2, lty=1)
          }

          legend(x=-110,y=0.4, legend=c("down shift","up shift", "net div"),fill=c(grey(0.3),rgb(1,0,0,alpha=0.5),NA),border=rep(NA,3),lty=c(NA,NA,1),lwd=c(NA,NA,2),col=c(NA,NA,"black"))

          mtext(text=setNames[q],side=3)

        dev.off()





##################
# COMBINED -- GLOBAL TREE data... 
#####
# Now COMBINE the RTT and LTT plots in the same time axis...
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
load(file="LTT_workspace_Marsup-Euarch-Lauras.Rdata")
library(plotrix); library(viridis); library(BAMMtools); library(gplots)

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
# Eulipotyphla, Afrosoricida, Rodentia, Primates, Scandentia, Pilosa, Macroscelidea, Cetartiodactyla, Chiroptera

ordNames_preKPg<-c("RODENTIA", "CHIROPTERA", "EULIPOTYPHLA", "PRIMATES", "ARTIODACTYLA", "AFROSORICIDA", "SCANDENTIA", "MACROSCELIDEA", "PILOSA")
ordNames_postKPg<-c("CARNIVORA","LAGOMORPHA", "PERISSODACTYLA","CINGULATA", "PHOLIDOTA", "PROBOSCIDEA", "HYRACOIDEA", "SIRENIA") #"DERMOPTERA","TUBULIDENTATA")

divTimesKPg<-read.table(file="LTT_divTimesToPlot.txt",header=TRUE)
ordinal<-divTimesKPg[which(divTimesKPg$LEVEL=="ordinal"),]
superord<-divTimesKPg[which(divTimesKPg$LEVEL=="superordinal"),]

superord_Bor<-superord[which(superord$majorClade1=="Boreoeutheria"),]
superord_Atlan<-superord[which(superord$majorClade1=="Atlanogenata"),]
superord_Euarch<-superord[which(superord$majorClade1=="Euarchontoglires"),]
superord_Lauras<-superord[which(superord$majorClade1=="Laurasiatheria"),]

	# BAMM RTT -- get data
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/_NDexp_10trees")
load(file="rateThroughTimeMatrix_10trees_1Ma.Rda")
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
rateMatrices<-get(paste("rateMatrices_1Ma",sep=""))

# I want to join the 10 tree files by column... starting at column 1 (or is it 36, etc?)
allCols<-c(); 
for(j in c(1:7,9:10)){
    allCols[j]<-dim(rateMatrices[[j]][["lambda"]])[2]
}
allCols # this is the number of -- 1 Ma -- 5 Ma bins per tree
maxBin<-max(na.omit(allCols)) # using tree 4, 204 intervals-- 41 intervals

for(j in c(1:7,9:10)){
  colToAdd<-maxBin-allCols[j] #dim(output_list_ALL[[j]][[k]])[2]

  for(k in c("lambda","mu")){ # the variables based on bins where NEED to even them out!

      if(colToAdd==0){ addCol<-NULL } else if(colToAdd==1) { 
      rowToAdd<-dim(rateMatrices[[j]][[k]])[1]
      addCol<-rep(rateMatrices[[j]][[k]][1,1],rowToAdd)
    } else {
      rowToAdd<-dim(rateMatrices[[j]][[k]])[1]
      addCol<-rep(rateMatrices[[j]][[k]][1,1],rowToAdd)
      for(q in 1:(colToAdd-1)){
        addCol<-cbind(addCol,rep(rateMatrices[[j]][[k]][1,1],rowToAdd))
      }
    }
      rateMatrices[[j]][[k]]<-cbind(addCol,rateMatrices[[j]][[k]])
      colnames(rateMatrices[[j]][[k]])<-(maxBin-1):0 #5*(40:0)
    }
}

allCols2<-c(); 
for(j in c(1:7,9:10)){
#    allCols2[j]<-dim(output_list_ALL[[j]][[5]])[2]
    allCols2[j]<-length(rateMatrices[[j]][["mu"]][1,])
}
allCols2 # this is the number of 5 Ma bins per tree
maxBin2<-max(na.omit(allCols2))

intToUse<-which(allCols==maxBin)[1]

join10runs_ALL<-rateMatrices[[intToUse]]
for(key in c("lambda","mu")){ # the variables based on bins where NEED to even them out!
	join10runs_ALL[[key]]<-rbind(rateMatrices[[1]][[key]], rateMatrices[[2]][[key]], rateMatrices[[3]][[key]], rateMatrices[[4]][[key]], 
        rateMatrices[[5]][[key]], rateMatrices[[6]][[key]], rateMatrices[[7]][[key]], 
        #rateMatrices[[8]][[key]], 
        rateMatrices[[9]][[key]], 
        rateMatrices[[10]][[key]])
}

  # BAMM lineage-level shifts -- get data
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"
unloadNamespace("BAMMtools"); unloadNamespace("gplots")

agesAllClades_origALL<-read.table(file=paste("divTime_24BAMMnodes_",bbone,"_MCC_target_DETAILS.txt", sep=""), header=TRUE)
#agesAllClades_ALL<-agesAllClades_origALL[which(agesAllClades_origALL$Indep==1),]
#allDat<-read.table(file=paste("divTime_24BAMMnodes_",bbone,"_MCC_target_cladeDR-mass-range-carn.txt", sep=""), header=TRUE, stringsAsFactors=FALSE)

firstShift<-c(1,3,5,7,8,10,11,13,15,17:29,31,32) #5 == keeping Placentalia
shiftDivDat1<-agesAllClades_origALL[firstShift,c("CLADE_label","ID_label","avgFactor","avgIncDec","mean","lower","upper","avgCladeDiv","cladeDiv_low95","cladeDiv_high95","nodeMCC","numTaxa")]
shiftDivDat<-shiftDivDat1[order(shiftDivDat1$avgFactor),]

  # TEST
    #pdf(file=paste("plotCI_cladeBAMMdiv-vs-cladeAge_withBars_new.pdf",sep=""),onefile=TRUE, width=6,height=4)
#    pdf(file=paste("plotCI_cladeBAMMfactor-vs-cladeAge_timeBars.pdf",sep=""),onefile=TRUE, width=6,height=4)
#    pchs<-rep(21,length(shiftDivDat[,1]))
#    #pchs[which(shiftDivDat$Indep==1)]<-16
#    alphaCol<-1
#    ptCols<-rep(rgb(1,0,0,alpha=alphaCol),length(shiftDivDat[,1]))
#    ptCols[which(shiftDivDat$avgIncDec=="downShift")]<-rgb(0,0,1,alpha=alphaCol)
#    ptCols[which(shiftDivDat$avgIncDec=="up-or-down")]<-grey(0.5,alpha=alphaCol)
#
#    yLab<-"BAMM shift factor"
#    #xLab<-"Shift timing (Ma)"
#
#    yVar<-shiftDivDat$avgFactor
#    #yVar[which(shiftDivDat$avgIncDec=="downShift")]<- -1.7
#    #yLow<-shiftDivDat$cladeDiv_low
#    #yHigh<-shiftDivDat$cladeDiv_high
#    xVar<-(-shiftDivDat$mean)
#    xLow<-(-shiftDivDat$lower)
#    xHigh<-(-shiftDivDat$upper)
#
#    plotCI(x=xVar, y=yVar,ui=xLow,li=xHigh, xlim=c(min(na.omit(xHigh)),max(na.omit(xLow))), ylim=c(min(yVar),max(yVar)), xlab="", ylab=yLab, sfrac=0,err="x", lwd=2,col="black",pt.bg=ptCols, scol="black",pch=pchs,font.lab=2,cex.axis=0.95,cex.lab=1, xaxt="n")
#    axis(side=1,line=2)
#
#    rug(xVar[which(shiftDivDat$avgIncDec=="upShift")], col=rgb(1,0,0,alpha=alphaCol), side=1, lwd=2, line=1)
#    rug(xVar[which(shiftDivDat$avgIncDec=="downShift")], col=rgb(0,0,1,alpha=alphaCol), side=1, lwd=2, line=2)
#    rug(xVar[which(shiftDivDat$avgIncDec=="up-or-down")], col=grey(0.5,alpha=alphaCol), side=1, lwd=2, line=2)
#
#    #plotCI(x=xVar, y=yVar,ui=yHigh,li=yLow, xlim=c(min(na.omit(xHigh)),max(na.omit(xLow))), ylim=c(min(yLow),max(yHigh)), xlab=xLab, ylab=yLab, sfrac=0,err="y", lwd=2,col="black",scol="black",pch=NA,font.lab=2,cex.axis=0.95,cex.lab=1)
#    #plotCI(add=TRUE,x=xVar, y=yVar,ui=xLow,li=xHigh, sfrac=0,err="x", lwd=2,col="black", pt.bg=ptCols, scol="black",pch=pchs,font.lab=2,cex.axis=0.95,cex.lab=1)
#
#    dev.off()


# Now get the global CoMET data...
##
CoMET_cols<-c("black",grey(0.5,alpha=0.5))
setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/",sep=""))
load(file="MAMPHY_10tree_CoMET_summary_1Ma_7shiftnoME.Rda")
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")


# Now get the global TreePAR data...
##
setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/TreePar_results",sep=""))
num_shifts=7; grid=1
numTrees<-10

setNames<-c("globalTree")

#setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/TreePar_results/",sep=""))
#num_shifts=7; grid=1
#numTrees<-10
#
##setNames<-c("globalTree","RODENTIA","CHIROPTERA","EULIPOTYPHLA","PRIMATES","CETARTIODACTYLA","AFROSORICIDA",
##        "Euarchontoglires","Laurasiatheria", "Marsupialia","Afroth.","Xenart.","CARNIVORA","LAGOMORPHA")
#setNames<-c("RODENTIA","CHIROPTERA","EULIPOTYPHLA","PRIMATES","CETARTIODACTYLA","AFROSORICIDA")
#
#
#shiftDeets_ALL<-vector("list",length(setNames))
#for (q in 1:length(setNames)){
#
#  results<-vector("list", length=numTrees)
#  for(i in 1:numTrees){
#      results[[i]]<-read.table(file=paste("TreePar_shifts",num_shifts,"_grid",grid,"Ma_clade_",setNames[q],"_tree",i,".txt",sep=""))
#  }
#  resAll<-do.call(rbind,results)
#    assign(paste("TP_resAll_",setNames[q],sep=""),resAll)
#
#  bestModels<-resAll[which(resAll$deltaAICc==0),]  
#    assign(paste("TP_bestModels_",setNames[q],sep=""),bestModels)
#
#  shiftDirs_all<-vector("list",length(bestModels[,1]))
#  shiftMags_all<-vector("list",length(bestModels[,1]))
#  shiftDivs_all<-vector("list",length(bestModels[,1]))
#  shiftTurns_all<-vector("list",length(bestModels[,1]))
#  for(j in 1:length(bestModels[,1])){
#    shiftDirs<-rep(NA,7)
#    shiftMags<-rep(NA,7)
#    shiftDivs<-list()
#    shiftTurns<-list()
#    for(k in 1:7){
#    if(is.na(bestModels[j,k+6])){ next } else if(bestModels[j,k+6] < bestModels[j,k+5]){shiftDirs[k]<-"upShift"} else {shiftDirs[k]<-"downShift"}  
#    if(is.na(shiftDirs[k])){ next } else if(shiftDirs[k]=="upShift"){ 
#          shiftMags[k]<-bestModels[j,k+5]/bestModels[j,k+6] 
#        } else {
#          shiftMags[k]<- -bestModels[j,k+5]/bestModels[j,k+6]
#        }
#      shiftDivs[[k]]<-c(bestModels[j,k+5],bestModels[j,k+6])
#      shiftTurns[[k]]<-c(bestModels[j,k+13],bestModels[j,k+14])
#    }
#    shiftDirs_all[[j]]<-shiftDirs
#    shiftMags_all[[j]]<-shiftMags
#    shiftDivs_all[[j]]<-do.call(rbind,shiftDivs)
#    shiftTurns_all[[j]]<-do.call(rbind,shiftTurns)
#  }
#
#  shiftDirs_ALL<-do.call(rbind, shiftDirs_all)
#  shiftMags_ALL<-do.call(rbind, shiftMags_all)
#  shiftTimes_ALL<-bestModels[,paste("ShiftTime",1:7,sep="")]
#  shiftDivs_ALL<-do.call(rbind, shiftDivs_all)
#  shiftTurns_ALL<-do.call(rbind, shiftTurns_all)
#
#  shiftDeets_all<-vector("list",length(bestModels[,1]))
#  for(j in 1:length(bestModels[,1])){
#    shiftDeets_all[[j]]<-cbind.data.frame(unlist(shiftTimes_ALL[j,]),shiftDirs_ALL[j,],round(shiftMags_ALL[j,],2))# ,(shiftDivs_ALL[j,]),(shiftTurns_ALL[j,]))
#  }
#  shiftDeets_ALL<-do.call(rbind, shiftDeets_all)
#  colnames(shiftDeets_ALL)<-c("shiftMa","dir","mag")#,"divRates","turnovers")
#  rownames(shiftDeets_ALL)<-c(1:length(shiftDeets_ALL[,1]))
#  shiftDeets_ALL<-na.omit(shiftDeets_ALL)
#
#  shiftDeets_FIN2<-cbind.data.frame(shiftDeets_ALL,shiftDivs_ALL,shiftTurns_ALL,round((shiftDivs_ALL/(1-shiftTurns_ALL)),4),round(((shiftDivs_ALL*shiftTurns_ALL)/(1-shiftTurns_ALL)),4))
#  colnames(shiftDeets_FIN2)<-c("shiftMa","dir","mag","div2","div1","turn2","turn1","sp2","sp1","ex2","ex1")
#  shiftDeets_FIN<-shiftDeets_FIN2[order(shiftDeets_FIN2$shiftMa,decreasing=TRUE),]
#
#  write.table(shiftDeets_FIN,file=paste("TreePar_shifts",num_shifts,"_grid",grid,"Ma_clade_",setNames[q],"_RES_bestModelShifts.txt",sep=""))  
#
#  assign(paste("TP_resFIN_",setNames[q],sep=""),shiftDeets_FIN)
#}
## save.image(file=paste("TreePar_shifts",num_shifts,"_grid",grid,"Ma_clade--ALL--RES_bestModelShifts.Rda",sep=""))


 # # Plot just the TreePar res...
 # pdf(file=paste("TreePar_shifts",num_shifts,"_grid",grid,"Ma_clade--ALL--RES_bestModelShifts_PLOT_upDownRug.pdf",sep=""),onefile=TRUE)
 # for (q in 1:length(setNames)){
 #   if(q==11){next}
 # 
 #   TP_resFIN<-get(paste("TP_resFIN_",setNames[q],sep=""))
 #   TP_bestModels<-get(paste("TP_bestModels_",setNames[q],sep=""))
 # 
 #   xMin<- -110
 #   xMax<-0
 #   times<-seq(from=xMin,to=xMax,by=5)
 #   histDown<-hist(-TP_resFIN$shiftMa[which(TP_resFIN$dir=="downShift" & TP_resFIN$shiftMa < -xMin)], breaks=times, plot=FALSE)$counts
 #   names(histDown)<-times[1:22]
 #   histUp<-hist(-TP_resFIN$shiftMa[which(TP_resFIN$dir=="upShift" & TP_resFIN$shiftMa < -xMin)], breaks=times, plot=FALSE)$counts
 #   names(histUp)<-times[1:22]
 # 
 #   plot(density(-TP_resFIN$shiftMa), type="n",xlim=c(xMin,xMax), yaxt="n",xaxt="n", bty="n",main="",ylab="",xlab="")
 #   
 #   rug(-TP_resFIN$shiftMa[which(TP_resFIN$dir=="downShift" & TP_resFIN$shiftMa < -xMin)], col=grey(0.3), lwd=2,line=0.5)
 #   rug(-TP_resFIN$shiftMa[which(TP_resFIN$dir=="upShift" & TP_resFIN$shiftMa < -xMin)], col=rgb(1,0,0,alpha=0.5), lwd=2, line=1.5)
 #   lines(x=c(-66,-66),y=c(-2,1),col=grey(0.5),lty=2,lwd=3)
 # 
 #   par(new=TRUE)
 #   barplot(histUp, space = 0, main="", col=rgb(1,0,0,alpha=0.5), border="white",lwd=0.5, xlim=c(0,22), ylim=c(0,10),xaxt="n", yaxt="n")
 #   barplot(add=TRUE,histDown, space = 0,  col=grey(0.3), border="white", lwd=0.5, xaxt="n", yaxt="n")
 #   axis(side=1,at=c(2,6,10,14,18,22),labels=c(-100,-80,-60,-40,-20,0), line=1.5)
 #   axis(side=4,at=c(0,2,4,6,8,10),labels=TRUE, line=0) # seq(from=0,to=max(c(histDown,histUp)),length.out=3)
 # 
 #   par(new=TRUE)
 #   plot(x=c(0,-TP_bestModels[1,paste("ShiftTime",1:7,sep="")]),y=TP_bestModels[1,paste("DivRate",1:8,sep="")], type="s", xlim=c(xMin,xMax), ylim=c(0,0.5),xaxt="n", bty="n",main="",ylab="",xlab="", col="black", lwd=2, lty=1)
 #   for(j in 2:length(TP_bestModels[,1])){
 #     points(x=c(0,-TP_bestModels[j,paste("ShiftTime",1:7,sep="")]),y=TP_bestModels[j,paste("DivRate",1:8,sep="")], type="s", xlim=c(xMin,xMax), xaxt="n", bty="n",main="",ylab="",xlab="", col="black", lwd=2, lty=1)
 #   }
 # 
 #   legend(x=-110,y=0.2, legend=c("down shift","up shift", "net div"),fill=c(grey(0.3),rgb(1,0,0,alpha=0.5),NA),border=rep(NA,3),lty=c(NA,NA,1),lwd=c(NA,NA,2),col=c(NA,NA,"black"))
 # 
 #   mtext(text=setNames[q],side=3)
 # }
 # 
 # dev.off()


load(file=paste("TreePar_shifts",num_shifts,"_grid",grid,"Ma_clade--ALL--RES_bestModelShifts.Rda",sep=""))
setwd("/Users/nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")

# BINNING the treePar results...
q<-1
  TP_resFIN<-get(paste("TP_resFIN_",setNames[q],sep=""))
  TP_bestModels<-get(paste("TP_bestModels_",setNames[q],sep=""))

  shifttimes<-TP_bestModels[,paste("ShiftTime",1:7,sep="")]
  netdiv<-TP_bestModels[,paste("DivRate",1:8,sep="")]
  turnover<-TP_bestModels[,paste("Turnover",1:8,sep="")]
  speciation<-(netdiv/(1-turnover))
    colnames(speciation)<-paste("Speciation",1:8,sep="")
  extinction<-((netdiv*turnover)/(1-turnover))
    colnames(extinction)<-paste("Extinction",1:8,sep="")

yVarNames<-c("netdiv","speciation","extinction")
for(z in 1:length(yVarNames)){
yVar<-get(yVarNames[z])
  varByGrid_ALL<-vector("list",length(shifttimes[,1]))

  for(i in 1:length(shifttimes[,1])){
    xVar<-unlist(c(0,shifttimes[i,][which(shifttimes[i,]<116)]))
    dat<-na.omit(cbind.data.frame(x=xVar,y=unlist(yVar[i,][1:length(xVar)])))
      rownames(dat)<-paste("t",1:length(dat[,1]),sep="")

    grid<-(seq(from=0,to=115,by=1))
    varByGrid<-rep(NA,length(grid))
    names(varByGrid)<-grid
    varByGrid[116]<-dat$y[length(dat$y)]

      for(k in 1:length(dat$x)){
        for(j in 1:(length(grid)-1)){
          if(dat$x[k] < grid[j+1]){ varByGrid[j]<-dat$y[k] } else if(dat$x[k] > grid[j+1] & dat$x[k] < grid[j+2]){ varByGrid[j]<-dat$y[k] } else { next }
        }
      }
      varByGrid_ALL[[i]]<-varByGrid
  }
  assign(paste("TP_",yVarNames[z],"_10trees",sep=""),do.call(rbind,varByGrid_ALL))
}

TP_RTT_all<-list(TP_netdiv_10trees[,length(grid):1], TP_speciation_10trees[,length(grid):1], TP_extinction_10trees[,length(grid):1])

# HISTOS for the treePar results...
xMin<- -110
xMax<-25
shiftUp<-TP_resFIN$shiftMa[which(TP_resFIN$dir=="upShift" & TP_resFIN$shiftMa < -xMin)]
shiftDown<-TP_resFIN$shiftMa[which(TP_resFIN$dir=="downShift" & TP_resFIN$shiftMa < -xMin)]

times<-seq(from=xMin,to=xMax,by=1)
histUp<-hist(-shiftUp, breaks=times, plot=FALSE)$counts
names(histUp)<-times[1:(length(times)-1)]
histDown<-hist(-shiftDown, breaks=times, plot=FALSE)$counts
names(histDown)<-times[1:(length(times)-1)]


save.image(file="Fig2_LTT_RTT_workspace.Rda")

load(file="Fig2_LTT_RTT_workspace.Rda")

# PLOT
##
# 4-part-- LTT higher taxa -- Placentals vs marsupials -- and then RTT of global Mammalia (w rate SHIFTS -- 3 different methods)

#pdf(file=paste("LTT_95polygon_pre-KPg_higherTaxLTT_globalRTT.pdf",sep=""), onefile=TRUE, width=5, height=8)
#pdf(file=paste("LTT_95polygon_pre-KPg_higherTaxLTT_globalRTT_3part_reOrder_fin_exclLast2Ma.pdf",sep=""), onefile=TRUE, width=4, height=8)
#pdf(file=paste("LTT_95polygon_pre-KPg_higherTaxLTT_globalRTT_3part_no2Ma_wBAMMshifts.pdf",sep=""), onefile=TRUE, width=4, height=9)
#pdf(file=paste("LTT_95polygon_pre-KPg_higherTaxLTT_globalRTT_3part_no2Ma_wBAMM-CoMET-TreePar.pdf",sep=""), onefile=TRUE, width=4, height=9)
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



##
# LOAD PER CLADE-- Prep all the per-Order data for plotting...
#####
# Get the LTTs...
#####
setwd("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
load(file="LTT_workspace_Marsup-Euarch-Lauras.Rdata")
library(plotrix); library(viridis); library(BAMMtools); library(gplots)

ordNames_preKPg<-c("RODENTIA", "CHIROPTERA", "EULIPOTYPHLA", "PRIMATES", "ARTIODACTYLA", "AFROSORICIDA", "SCANDENTIA", "MACROSCELIDEA", "PILOSA")

divTimesKPg<-read.table(file="LTT_divTimesToPlot.txt",header=TRUE)
ordinal<-divTimesKPg[which(divTimesKPg$LEVEL=="ordinal"),]


# BAMM lineage-level shifts -- get data
#####
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")
bbone<- "NDexp" # "FBD"
#unloadNamespace("BAMMtools"); unloadNamespace("gplots")

agesAllClades_origALL<-read.table(file=paste("divTime_24BAMMnodes_",bbone,"_MCC_target_DETAILS.txt", sep=""), header=TRUE)
#agesAllClades_ALL<-agesAllClades_origALL[which(agesAllClades_origALL$Indep==1),]
#allDat<-read.table(file=paste("divTime_24BAMMnodes_",bbone,"_MCC_target_cladeDR-mass-range-carn.txt", sep=""), header=TRUE, stringsAsFactors=FALSE)

firstShift<-c(1,3,5,7,8,10,11,13,15,17:29,31,32) #5 == keeping Placentalia
shiftDivDat2<-agesAllClades_origALL[firstShift,c("CLADE_label","ID_label","avgFactor","avgIncDec","mean","lower","upper","avgCladeDiv","cladeDiv_low95","cladeDiv_high95","nodeMCC","numTaxa")]
shiftDivDat1<-shiftDivDat2[order(shiftDivDat2$avgFactor),]

orderOfShift<-rep(NA,length(shiftDivDat1[,1]))
for(i in 1: length(orderOfShift)){
  if(shiftDivDat1[i,"ID_label"] %in% c("Q","W","V","X","U","T","S","R")) {orderOfShift[i]<-"RODENTIA"} 
  else if (shiftDivDat1[i,"ID_label"] %in% c("J","L","G","K","I","H") ) {orderOfShift[i]<-"CHIROPTERA"} 
  else if (shiftDivDat1[i,"ID_label"] %in% c("M") ) {orderOfShift[i]<-"EULIPOTYPHLA"} 
  else if (shiftDivDat1[i,"ID_label"] %in% c("O","N") ) {orderOfShift[i]<-"PRIMATES"} 
  else if (shiftDivDat1[i,"ID_label"] %in% c("F","E") ) {orderOfShift[i]<-"CETARTIODACTYLA"} 
  else {orderOfShift[i]<-"notInKPgOrds"} 
}
shiftDivDat<-cbind.data.frame(orderOfShift,shiftDivDat1)


# BAMM rate through time data back in...
####
library(BAMMtools); library(coda); library(phytools); library(ape); library(phangorn); library(viridis)
library(TESS); library(matrixStats)
setwd("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/")
source("tess.plot.output2.R")

#setwd(paste("/mnt/data/personal/nateu/Nate_Backup/BAMM_analyses/",sep=""))
setwd("/Users/Nate/Desktop/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/BAMM_analyses/_COMPLETE_100Mgens/_NDexp_10trees")
load(file="RTT_workspace_ords-and-higher_bins1Ma_10trees.Rda")

setNames<-c("RODENTIA","CHIROPTERA","EULIPOTYPHLA","PRIMATES","CETARTIODACTYLA")#,"AFROSORICIDA")#,
        #"Euarchontoglires","Laurasiatheria", "Marsupialia")#,"Afroth.","Xenart.","CARNIVORA","LAGOMORPHA")

setNames2<-c("RODENTIA","CHIROPTERA","EULIPOTYPHLA","PRIMATES","ARTIODACTYLA")#,"AFROSORICIDA")#,

for(v in 1:length(setNames)){ 

#v=1
# get the BAMM data to plot
rateMatrices<-get(paste("rateMatrices_1Ma_",setNames[v],"_ALL10",sep=""))

# I want to join the 10 tree files by column... starting at column 1 (or is it 36, etc?)
allCols<-c(); 
for(j in c(1:7,9:10)){
    allCols[j]<-dim(rateMatrices[[j]][["lambda"]])[2]
}
allCols # this is the number of -- 1 Ma -- 5 Ma bins per tree
maxBin<-max(na.omit(allCols)) # using tree 4, 204 intervals-- 41 intervals

for(j in c(1:7,9:10)){
  colToAdd<-maxBin-allCols[j] #dim(output_list_ALL[[j]][[k]])[2]

  for(k in c("lambda","mu")){ # the variables based on bins where NEED to even them out!

      if(colToAdd==0){ addCol<-NULL } else if(colToAdd==1) { 
      rowToAdd<-dim(rateMatrices[[j]][[k]])[1]
      addCol<-rep(rateMatrices[[j]][[k]][1,1],rowToAdd)
    } else {
      rowToAdd<-dim(rateMatrices[[j]][[k]])[1]
      addCol<-rep(rateMatrices[[j]][[k]][1,1],rowToAdd)
      for(q in 1:(colToAdd-1)){
        addCol<-cbind(addCol,rep(rateMatrices[[j]][[k]][1,1],rowToAdd))
      }
    }
      rateMatrices[[j]][[k]]<-cbind(addCol,rateMatrices[[j]][[k]])
      colnames(rateMatrices[[j]][[k]])<-(maxBin-1):0 #5*(40:0)
    }
}

allCols2<-c(); 
for(j in c(1:7,9:10)){
#    allCols2[j]<-dim(output_list_ALL[[j]][[5]])[2]
    allCols2[j]<-length(rateMatrices[[j]][["mu"]][1,])
}
allCols2 # this is the number of 5 Ma bins per tree
maxBin2<-max(na.omit(allCols2))

intToUse<-which(allCols==maxBin)[1]

join10runs_ALL<-rateMatrices[[intToUse]]
for(key in c("lambda","mu")){ # the variables based on bins where NEED to even them out!
  join10runs_ALL[[key]]<-rbind(rateMatrices[[1]][[key]], rateMatrices[[2]][[key]], rateMatrices[[3]][[key]], rateMatrices[[4]][[key]], 
        rateMatrices[[5]][[key]], rateMatrices[[6]][[key]], rateMatrices[[7]][[key]], 
        #rateMatrices[[8]][[key]], 
        rateMatrices[[9]][[key]], 
        rateMatrices[[10]][[key]])
}

assign(paste("rateMatrices_",setNames[v],sep=""),join10runs_ALL)
}


# Load in PER CLADE-- CoMET results...
####
setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/CoMET_analyses_MamPhy_subclades/_res_all10trees_all4runs_bin1Ma_Rdata",sep=""))

for (q in 1:length(setNames)){
  load(file=paste("MAMPHY_10tree_CoMET_summary_1Ma_7shiftnoME_",setNames[q],".Rda",sep=""))
  assign(paste("join4runs_ALL_",setNames[q],sep=""),join4runs_ALL)
}

  # Get just the CoMET shifts >2 BF...
for(v in 1:length(setNames)){ 

  v=1

  #CoMET rate shifts...
CoMET_cols<-c(rgb(1,0,0,alpha=0.5),rgb(0,0,1,alpha=0.5))  # c("black",grey(0.7,alpha=0.6))

output<-get(paste("join4runs_ALL_",setNames[v],sep=""))
type<-"speciation shift times"
criticalPP <- output[[grep(strsplit(type, " ")[[1]][1], grep("CriticalPosteriorProbabilities", names(output), value = TRUE), value = TRUE)]]

xMin<- -110
xMax<-25

type<-"extinction shift times"
colBars<-CoMET_cols[2]

thisOutput <- output[[type]]
meanThisOutput <- colMeans(thisOutput)
names(meanThisOutput)<-c((length(meanThisOutput)-1):0)
meanThisOutput[(length(meanThisOutput)-1):length(meanThisOutput)]<-NA # excludes last 2 Ma from the plot
shiftSubstantial_ex<-as.numeric(names(meanThisOutput[which(meanThisOutput >= criticalPP[1])])) # gets shifts >= 2 BF

#par(new=TRUE, lwd=0.1)
par(lwd=0.1)
barplot(rev(meanThisOutput), space = 0, col = colBars, border = "white", main = "", ylab = "", xlab = "", xaxt = "n", yaxt = "n", xlim=c(-xMin,-xMax),ylim = c(0,0.5))
rug(shiftSubstantial_ex, col=colBars, lwd=2,line=0.5)

type<-"speciation shift times"
colBars<-CoMET_cols[1]

thisOutput <- output[[type]]
meanThisOutput <- colMeans(thisOutput)
names(meanThisOutput)<-c((length(meanThisOutput)-1):0)
meanThisOutput[(length(meanThisOutput)-1):length(meanThisOutput)]<-NA # excludes last 2 Ma from the plot
shiftSubstantial_sp<-as.numeric(names(meanThisOutput[which(meanThisOutput >= criticalPP[1])])) # gets shifts >= 2 BF

par(new=TRUE, lwd=0.1)
barplot(rev(meanThisOutput), space = 0, col = colBars, border = "white", main = "", ylab = "", xlab = "", xaxt = "n", yaxt = "n", xlim=c(-xMin,-xMax),ylim = c(0,0.5))
rug(shiftSubstantial_sp, col=colBars, lwd=2,line=1.5)

axis(side=1)




# Load in PER CLADE-- TreePar results...
####
setwd(paste("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/CoMET_and_TreePar/TreePar_results/",sep=""))
num_shifts=7; grid=1
numTrees<-10

load(file=paste("TreePar_shifts",num_shifts,"_grid",grid,"Ma_clade--ALL--RES_bestModelShifts.Rda",sep=""))
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")

#setNames<-c("globalTree","RODENTIA","CHIROPTERA","EULIPOTYPHLA","PRIMATES","CETARTIODACTYLA","AFROSORICIDA",
#        "Euarchontoglires","Laurasiatheria", "Marsupialia","Afroth.","Xenart.","CARNIVORA","LAGOMORPHA")
setNames<-c("RODENTIA","CHIROPTERA","EULIPOTYPHLA","PRIMATES","CETARTIODACTYLA","AFROSORICIDA")

# BINNING the treePar results...
for(q in 1:length(setNames)){

  TP_resFIN<-get(paste("TP_resFIN_",setNames[q],sep=""))
  TP_bestModels<-get(paste("TP_bestModels_",setNames[q],sep=""))

  shifttimes<-TP_bestModels[,paste("ShiftTime",1:7,sep="")]
  netdiv<-TP_bestModels[,paste("DivRate",1:8,sep="")]
  turnover<-TP_bestModels[,paste("Turnover",1:8,sep="")]
  speciation<-(netdiv/(1-turnover))
    colnames(speciation)<-paste("Speciation",1:8,sep="")
  extinction<-((netdiv*turnover)/(1-turnover))
    colnames(extinction)<-paste("Extinction",1:8,sep="")

yVarNames<-c("netdiv","speciation","extinction")
for(z in 1:length(yVarNames)){
yVar<-get(yVarNames[z])
  varByGrid_ALL<-vector("list",length(shifttimes[,1]))

  for(i in 1:length(shifttimes[,1])){
    xVar<-unlist(c(0,shifttimes[i,][which(shifttimes[i,]<116)]))
    dat<-na.omit(cbind.data.frame(x=xVar,y=unlist(yVar[i,][1:length(xVar)])))
      rownames(dat)<-paste("t",1:length(dat[,1]),sep="")

    grid<-(seq(from=0,to=115,by=1))
    varByGrid<-rep(NA,length(grid))
    names(varByGrid)<-grid
    varByGrid[116]<-dat$y[length(dat$y)]

      for(k in 1:length(dat$x)){
        for(j in 1:(length(grid)-1)){
          if(dat$x[k] < grid[j+1]){ varByGrid[j]<-dat$y[k] } else if(dat$x[k] > grid[j+1] & dat$x[k] < grid[j+2]){ varByGrid[j]<-dat$y[k] } else { next }
        }
      }
      varByGrid_ALL[[i]]<-varByGrid
  }
  assign(paste("TP_",yVarNames[z],"_10trees_",setNames[q],sep=""),do.call(rbind,varByGrid_ALL))
}

}

# for global
#TP_RTT_all<-list(TP_netdiv_10trees[,length(grid):1], TP_speciation_10trees[,length(grid):1], TP_extinction_10trees[,length(grid):1])





# PLOT
###
# For each Order, do FULL RTT plots of BAMM-CoMET-TreePar and their shifts, then decide from there...
setNames<-c("RODENTIA","CHIROPTERA","EULIPOTYPHLA","PRIMATES","CETARTIODACTYLA")#,"AFROSORICIDA")
xMin<- -110
xMax<-25
grid<-(seq(from=0,to=115,by=1))
CoMET_cols<-c(rgb(1,0,0,alpha=0.5),rgb(0,0,1,alpha=0.5))  # c("black",grey(0.7,alpha=0.6))


pdf(file=paste("RTT_3part_no2Ma_wBAMM-only_preKPg_5ords_cols.pdf",sep=""), onefile=TRUE, width=4, height=6)
#pdf(file=paste("RTT_3part_no2Ma_wBAMM-CoMET-TreePar_preKPg6Ords.pdf",sep=""), onefile=TRUE, width=4, height=6)

xLim<-c(-110,25)
#quartz(width=5,height=8)
layout(matrix(c(1:5), 5, 1, byrow = TRUE), widths=rep(1,3))#,heights=c(0.333,0.333,0.333))
#layout(matrix(c(1:3), 3, 1, byrow = TRUE), widths=rep(1,3),heights=c(0.333,0.333,0.333))

op <- par(#mfrow = c(3,1),
          oma = c(5,4,2,0) + 0.1, #c(bottom, left, top, right)
          mar = c(0.5,0,1,1) + 0.1, xpd=TRUE)

LegRCex<-0.9 


for(q in 1:length(setNames)){

# LOAD the BAMM RTT results...
join10runs_ALL<-get(paste("rateMatrices_",setNames[q],sep=""))

# LOAD the BAMM shift results...
orderBAMMshifts<-shiftDivDat[which(shiftDivDat$orderOfShift==setNames[q]),]

# # LOAD the CoMET RTT results... is same as the shift results...
# output<-get(paste("join4runs_ALL_",setNames[q],sep=""))
# type<-"speciation shift times"
# criticalPP <- output[[grep(strsplit(type, " ")[[1]][1], grep("CriticalPosteriorProbabilities", names(output), value = TRUE), value = TRUE)]]
# 
# # LOAD the TreePar RTT results...
# TP_RTT_all<-list(get(paste("TP_netdiv_10trees_",setNames[q],sep=""))[,length(grid):1], get(paste("TP_speciation_10trees_",setNames[q],sep=""))[,length(grid):1],get(paste("TP_extinction_10trees_",setNames[q],sep=""))[,length(grid):1])
# 
# # LOAD the TreePar shift histos...
# TP_resFIN<-get(paste("TP_resFIN_",setNames[q],sep=""))
# shiftUp<-TP_resFIN$shiftMa[which(TP_resFIN$dir=="upShift" & TP_resFIN$shiftMa < -xMin)]
# shiftDown<-TP_resFIN$shiftMa[which(TP_resFIN$dir=="downShift" & TP_resFIN$shiftMa < -xMin)]
# 
# times<-seq(from=xMin,to=xMax,by=1)
# histUp<-hist(-shiftUp, breaks=times, plot=FALSE)$counts
# names(histUp)<-times[1:(length(times)-1)]
# histDown<-hist(-shiftDown, breaks=times, plot=FALSE)$counts
# names(histDown)<-times[1:(length(times)-1)]


# now make the plot....
# ~~~~~~
  # BAMM RTT of global mammals
#library(BAMMtools)
ratetypes <- c("netdiv","speciation","extinction")
intervalCol <- grey(0.5) viridis(3) #c("blue","green","darkyellow")
#avgCol <- grey(0,alpha=1) #viridis(3) #c("blue","green","darkyellow")#rep("red",3)
#lineTypes <- c(1,5,4)
lineTypes <- c(1,1,1)
avgCols <- c(intervalCol,"black","white")

cexRLab<-0.65
cex.lab = 1

for(z in c(2,3,1)){
rmat<-join10runs_ALL
rates <- list((rmat$lambda - rmat$mu),rmat$lambda, rmat$mu)

rate<-rates[[z]][,1:(length(rates[[z]][1,])-2)] # exclude last 2 Ma

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
    rmat$times <- (max(rmat$times) - rmat$times)[1:length(rate[1,])] # exclude last 2 Ma

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
# offsetBamm<- 0.4
  offsetBamm<- 0

if(z == 2){
        #plot.new()
        #par(mar = mar)          
    xMin <- -110 #maxTime
        xMax <- 25
        yMin <- 0 #min(poly[[1]][, 2])
        yMax <- 0.4+offsetBamm #max(poly[[1]][, 2])
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
    pchs<-rep(21,length(orderBAMMshifts[,1]))
    #pchs[which(shiftDivDat$Indep==1)]<-16
    alphaCol<-0.7
    ptCols<-rep(rgb(1,0,0,alpha=alphaCol),length(orderBAMMshifts[,1]))
    ptCols[which(orderBAMMshifts$avgIncDec=="downShift")]<-rgb(0,0,1,alpha=alphaCol)
    ptCols[which(orderBAMMshifts$avgIncDec=="up-or-down")]<-grey(0.5,alpha=alphaCol)

    yVar<-orderBAMMshifts$avgFactor
    xVar<-(-orderBAMMshifts$mean)
    xLow<-(-orderBAMMshifts$lower)
    xHigh<-(-orderBAMMshifts$upper)

    bammShift_yMin<-1 # -2.5
    
    if(length(yVar)!=0) { 
    par(new=TRUE)
    plotCI(x=xVar, y=yVar+(3*offsetBamm),ui=xLow,li=xHigh, xlim=c(xMin, xMax), ylim=c(bammShift_yMin,6), xlab="", ylab="", sfrac=0,err="x", lwd=2,col="black",pt.bg=ptCols, scol="black",pch=pchs,font.lab=2,cex.axis=0.95,cex.lab=1, xaxt="n", yaxt="n",bty="n", cex=1.5)
    text(x=xVar, y=yVar+(3*offsetBamm), labels=as.vector(orderBAMMshifts$ID_label), pos=1, offset=0.7)
    
    axis(side=4,line= -5,at=c(1,2,3,4)+(3*offsetBamm),labels=c(1,2,3,4))
    mtext(side = 4, text = "Shift factor", line = -3, cex = cexRLab, adj=0+offsetBamm)
    text(x=-115, y= 6,labels = paste(setNames[q],"-wide rates",sep=""), cex = 1.2,font=2,pos=4)
  legend(x=3,y=6,legend=c("up","down","up/down"), pt.bg=c(rgb(1,0,0,alpha=alphaCol),rgb(0,0,1,alpha=alphaCol),grey(0.5,alpha=alphaCol)), pch=21, cex=LegRCex, pt.cex=1.5, seg.len=4, bty="n", title="BAMM\nper-lineage\nshifts") 
    } else {
    text(x=-115, y= yMax,labels = paste(setNames[q],"-wide rates",sep=""), cex = 1.2,font=2,pos=4)
    }
#    text(x=-115, y= 6,labels = "(b) Mammalia-wide rates", cex = 1.2,font=2,pos=4)

  #  rug(xVar[which(shiftDivDat$avgIncDec=="upShift")], col=rgb(1,0,0,alpha=alphaCol), side=1, lwd=2, line=-5-(22*offset))
  #  rug(xVar[which(shiftDivDat$avgIncDec=="downShift")], col=rgb(0,0,1,alpha=alphaCol), side=1, lwd=2, line=-5-(22*offset))
  #  rug(xVar[which(shiftDivDat$avgIncDec=="up-or-down")], col=grey(0.5,alpha=alphaCol), side=1, lwd=2, line=-5-(22*offset))

#    legend(x=-115,y=5.5,legend=c("speciation","extinction","net diversification"), col=c(avgCols[2:3],avgCols[1]), lty=c(1,1,1), lwd=lwd, cex=0.9, seg.len=4, bty="o",bg=grey(0.85)) 


#   legend(x=-115,y=offset,legend=c("speciation shift","extinction shift"), fill=CoMET_cols, border="black", cex=1, bty="n") 
#   text(x=-115,y=0.25+offset,labels="(c) Mammalia-wide rates",cex=1.2,font=2,pos=4)

#text(x=-115,y=0.22+offset,labels="Rates through time",cex=1,font=2,pos=4)
#text(x=-115,y=-0.03+offset,labels="Tree-wide rate shifts",cex=1,font=2,pos=4)

#   #CoMET RTT plots... 
# ratetypes <- c("netdiv","speciation","extinction")
# intervalCol <- grey(0.2) #viridis(3) #c("blue","green","darkyellow")
# #avgCol <- grey(0,alpha=1) #viridis(3) #c("blue","green","darkyellow")#rep("red",3)
# #lineTypes <- c(1,5,4)
# lineTypes <- c(1,1,1)
# avgCols <- c(grey(0.5),"black","white")
# 
# for(z in c(2,3,1)){
# rates <-list((output[["net-diversification rates"]]),(output[["speciation rates"]]),(output[["extinction rates"]]))
# 
# rate<-rates[[z]][,1:(length(rates[[z]][1,])-2)] # exclude last 2 Ma
# colnames(rate)<-(length(rates[[z]][1,])-1):2
# 
# ratelabel<-"rate"
# ratetype<-ratetypes[z]
# lineType<-lineTypes[z]
# avgCol<-avgCols[z]
# 
#     maxTime <- max(output[["intervals"]])
#     nanCol <- apply(rate, 2, function(x) any(is.nan(x)))
#     rate <- rate[, which(nanCol == FALSE)]
#     #plotTimes <- output[["intervals"]][which(nanCol == FALSE)]
#     plotTimes <- output[["intervals"]][1:(length(rates[[z]][1,])-2)] # exclude last 2 Ma
# 
#         #mm <- apply(rate, 2, quantile, prob = c(0.025, 0.975)) #intervals)
# 
#         mm <- apply(rate, 2, quantile, intervals)
#         poly <- list()
#         q1 <- 1
#         q2 <- nrow(mm)
#         repeat {
#             if (q1 >= q2) {
#                 break
#             }
#             a <- as.data.frame(cbind(rmat$times, mm[q1, ]))
#             b <- as.data.frame(cbind(rmat$times, mm[q2, ]))
#             b <- b[rev(rownames(b)), ]
#             colnames(a) <- colnames(b) <- c("x", "y")
#             poly[[q1]] <- rbind(a, b)
#             q1 <- q1 + 1
#             q2 <- q2 - 1
#         }
# 
#     harmMean<-function(x){
#       res<-(1/mean(1/x))
#       return(res)
#     }
#   avg <- apply(rate,2,median)#colMeans(rate)
# #  offsetCoMET<- 0.35
#   offsetCoMET<- 0.1
# 
# if(z == 2){
#         #plot.new()
#         #par(mar = mar)          
#     xMin <- -110 #maxTime
#         xMax <- 25
#         yMin <- 0 #min(poly[[1]][, 2])
# 
# CoMETrate_yMax<-0.5 #0.29
#         yMax <- CoMETrate_yMax+offsetCoMET #max(poly[[1]][, 2])
#         if (yMin >= 0) { yMin <- 0 }
#         #plot.window(xlim = xLim, ylim = c(yMin, yMax))
# 
# #    par(new=TRUE)
#     plot(x = -plotTimes, y = (avg+offsetCoMET),lty=2,lwd=1,xlim=c(xMin, xMax),ylim=c(yMin, yMax),type="l", col="white",main=NULL,cex.axis=0.8,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",yaxt="n",xaxt="n",bty="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
#         
#         axis(side=2, at = seq(0,0.2,by=0.05)+offsetCoMET, labels=c(0,NA,0.1,NA,0.2),cex.axis = cex.axis, las = 1, line=0.25)
#         mtext(side = 2, text = "rate", line = 3, cex = cex.lab, adj=0+offsetCoMET)
# 
#     abline(v=-66,lwd=2,lty=2,col="dark grey")
#        }
# 
#         for (i in 1:length(poly)) {
#             polygon(x = -poly[[i]][, 1], y = (poly[[i]][, 2]+offsetCoMET), col = transparentColor(intervalCol, opacity), border = NA)
#         }
# 
#   #    polygon(x = -c(0:ncol(mm), ncol(mm):0), 
#   #              y = c(c(mm[1, 1], mm[1, 
#   #                ]), rev(c(mm[2, 1], mm[2, 
#   #                ]))), border = NA, col = transparentColor(intervalCol, 0.5))
# 
#     lines(x = -plotTimes, y = (avg+offsetCoMET), lwd = lwd, col = avgCol, lty=lineType)
# }
# 
# legend(x=-130,y=0.6,legend=c("speciation","extinction","net diversification"), col=c(avgCols[2:3],avgCols[1]), lty=c(1,1,1), lwd=lwd, cex=0.75, seg.len=4, bty="o",bg=grey(0.85)) 
# 
# 
#   #CoMET rate shifts...
# type<-"extinction shift times"
# colBars<-CoMET_cols[2]
# 
# thisOutput <- output[[type]]
# meanThisOutput <- colMeans(thisOutput)
# names(meanThisOutput)<-(length(thisOutput[1,])-1):0
# meanThisOutput[(length(meanThisOutput)-1):(length(meanThisOutput))]<-NA # excludes last 2 Ma from the plot
# par(new=TRUE, lwd=0.1)
# barplot(rev(meanThisOutput), space = 0, col = colBars, border = "white", main = "", ylab = "", xlab = "", xaxt = "n", yaxt = "n", xlim=c(-xMin,-xMax),ylim = c(0,0.5))
# 
# type<-"speciation shift times"
# colBars<-CoMET_cols[1]
# 
# thisOutput <- output[[type]]
# meanThisOutput <- colMeans(thisOutput)
# names(meanThisOutput)<-(length(thisOutput[1,])-1):0
# meanThisOutput[(length(meanThisOutput)-1):(length(meanThisOutput))]<-NA # excludes last 2 Ma from the plot
#   par(new=TRUE, lwd=0.1)
#   barplot(rev(meanThisOutput), space = 0, col = colBars, border = "white", main = "", ylab = "", xlab = "", xaxt = "n", yaxt = "n", xlim=c(-xMin,-xMax),ylim = c(0,0.5))
# 
# #axis(2, at = c(0,criticalPP[1],criticalPP[2]), labels = c(0,2,6), line=0.25)#, line= -4)
# axis(4, at = c(0,criticalPP[1]), labels = c(0,2), line=-5) #0.25)#, line= -4)
# mtext(side = 4, text = "2 ln BF", line = -3, cex = cexRLab, adj=0) 
# 
# legend(x=-3,y=0.4,legend=c("speciation","extinction"), pt.bg=c(rgb(1,0,0,alpha=alphaCol),rgb(0,0,1,alpha=alphaCol)), pch=22, cex=LegRCex, pt.cex=1.5, seg.len=4, bty="n", title="CoMET\ntree-wide\nshifts") 
# 
# #points(rev(meanThisOutput),type="s",main = "", ylab = "", xlab = "", xaxt = "n", yaxt = "n", xlim=c(-xMin,-xMax),ylim = c(0,1.5), bty="n")
# 
# #lines(x=c(115,0), y=c(criticalPP[1],criticalPP[1]), lty = 1, col=grey(0.5))
# #lines(x=c(115,0), y=c(criticalPP[2],criticalPP[2]), lty = 1, col=grey(0.5))
# 
# 
# 
#  # TreePar RTT plots... 
# ratetypes <- c("netdiv","speciation","extinction")
# intervalCol <- grey(0.2) #viridis(3) #c("blue","green","darkyellow")
# #avgCol <- grey(0,alpha=1) #viridis(3) #c("blue","green","darkyellow")#rep("red",3)
# #lineTypes <- c(1,5,4)
# lineTypes <- c(1,1,1)
# avgCols <- c(grey(0.5),"black","white")
# 
# for(z in c(2,3,1)){
# rates <-TP_RTT_all
# 
# crownAge<-116-round(max(output[["intervals"]]),0)
# rate<-rates[[z]][,crownAge:(length(rates[[z]][1,])-2)] # exclude last 2 Ma -- AND that prior to the crown age...
# #colnames(rate)<-203:2
# 
# ratelabel<-"rate"
# ratetype<-ratetypes[z]
# lineType<-lineTypes[z]
# avgCol<-avgCols[z]
# 
#     maxTime <- as.numeric(colnames(rate)[1]) #max(output[["intervals"]])
#     nanCol <- apply(rate, 2, function(x) any(is.nan(x)))
#     rate <- rate[, which(nanCol == FALSE)]
#     #plotTimes <- output[["intervals"]][which(nanCol == FALSE)]
#     plotTimes <- as.numeric(colnames(rate)) # output[["intervals"]][1:202] # >> already excludes last 2 Ma
# 
#         #mm <- apply(rate, 2, quantile, prob = c(0.025, 0.975)) #intervals)
# 
#         mm <- apply(rate, 2, quantile, intervals)
#         poly <- list()
#         q1 <- 1
#         q2 <- nrow(mm)
#         repeat {
#             if (q1 >= q2) {
#                 break
#             }
#             a <- as.data.frame(cbind(plotTimes, mm[q1, ]))
#             b <- as.data.frame(cbind(plotTimes, mm[q2, ]))
#             b <- b[rev(rownames(b)), ]
#             colnames(a) <- colnames(b) <- c("x", "y")
#             poly[[q1]] <- rbind(a, b)
#             q1 <- q1 + 1
#             q2 <- q2 - 1
#         }
# 
#     harmMean<-function(x){
#       res<-(1/mean(1/x))
#       return(res)
#     }
#   avg <- apply(rate,2,median)#colMeans(rate)
# #  offsetTreePar<- 0.15
#   offsetTreePar<- 0.2
# 
# if(z == 2){
#         #plot.new()
#         #par(mar = mar)          
#     xMin <- -110 #maxTime
#         xMax <- 25
#         yMin <- 0 #min(poly[[1]][, 2])
# TreePar_yMax<-0.6
#         yMax <- TreePar_yMax+offsetTreePar #max(poly[[1]][, 2])
#         if (yMin >= 0) { yMin <- 0 }
#         #plot.window(xlim = xLim, ylim = c(yMin, yMax))
# 
# #    par(new=TRUE)
#     plot(x = -plotTimes, y = (avg+offsetTreePar),lty=2,lwd=1,xlim=c(xMin, xMax),ylim=c(yMin, yMax),type="l", col="white",main=NULL,cex.axis=0.8,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",yaxt="n",xaxt="n",bty="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
#         
#         axis(side=2, at = seq(0,0.2,by=0.05)+offsetTreePar, labels=c(0,NA,0.1,NA,0.2),cex.axis = cex.axis, las = 1, line=0.25)
#         mtext(side = 2, text = "rate", line = 3, cex = cex.lab, adj=0+offsetTreePar)
# 
#     abline(v=-66,lwd=2,lty=2,col="dark grey")
#     text(x=66,y=TreePar_yMax,labels="K-Pg\nboundary",cex=1,font=2,pos=4, col=grey(0.5))
# 
#        }
# 
#         for (i in 1:length(poly)) {
#             polygon(x = -poly[[i]][, 1], y = (poly[[i]][, 2]+offsetTreePar), col = transparentColor(intervalCol, opacity), border = NA)
#         }
# 
#   #    polygon(x = -c(0:ncol(mm), ncol(mm):0), 
#   #              y = c(c(mm[1, 1], mm[1, 
#   #                ]), rev(c(mm[2, 1], mm[2, 
#   #                ]))), border = NA, col = transparentColor(intervalCol, 0.5))
# 
#     lines(x = -plotTimes, y = (avg+offsetTreePar), lwd = lwd, col = avgCol, lty=lineType)
# }
# 
#     axis(side=1, at= rev(c(0,-20,-40,-60,-80,-100)), labels=rev(-c(0,-20,-40,-60,-80,-100)), cex=cex.axis, line=0.5)
#     mtext(side = 1, text = "time before present (Ma)", line = 3.5, cex = cex.lab, adj=0.4)
# 
#   ## TreePar shift HISTOS
# par(new=TRUE)
# barplot(histUp, space = 0, main="", col=rgb(1,0,0,alpha=0.5), border="white",lwd=0.5, ylim=c(0,10),yaxt="n",xaxt="n")
# barplot(add=TRUE,histDown, space = 0,  col=rgb(0,0,1,alpha=0.5), border="white", lwd=0.5, xaxt="n", yaxt="n")
# 
# axis(side=4,at=seq(from=0,to=4,length.out=5),labels=c(0,NA,0.2,NA,0.4), line=-5) # seq(from=0,to=max(c(histDown,histUp)),length.out=3)
# mtext(side=4,text="Shifts per tree", line=-3, cex=cexRLab, adj=-0.2)
# 
# legend(x=120,y=7,legend=c("up","down"), pt.bg=c(rgb(1,0,0,alpha=alphaCol),rgb(0,0,1,alpha=alphaCol)), pch=22, cex=LegRCex, pt.cex=1.5, seg.len=4, bty="n", title="TreePar\ntree-wide\nshifts") 
# 
# # plot(density(-TP_resFIN$shiftMa), type="n",xlim=c(xMin,xMax), yaxt="n",xaxt="n", bty="n",main="",ylab="",xlab="")
# # rug(-shiftDown-1, col=grey(0.3), lwd=2,line=0.5)
# # rug(-shiftUp-1, col=rgb(1,0,0,alpha=0.5), lwd=2, line=1.5)
# # lines(x=c(-66,-66),y=c(-2,1),col=grey(0.5),lty=2,lwd=3)
# # axis(side=1,at=seq(from=-100,to=0,by=20),labels= seq(from=-100,to=0,by=20), line=1.5)

}

dev.off()



#########




# PLOT all 6 clades...
setwd("/Users/Nate/Desktop/VertLife-Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/")

#pdf(file=paste("rateThroughTime_10trees_combo_bins1Ma_3rates_KPg_ALL.pdf",sep=""), onefile=TRUE, width=6, height=4)
#pdf(file=paste("rateThroughTime_10trees_combo_bins1Ma_3rates_KPg_ALL_6ordsInOne.pdf",sep=""), onefile=TRUE, width=5, height=10)
pdf(file=paste("rateThroughTime_10trees_combo_bins1Ma_3rates_KPg_ALL_6ordsInOne_wCoMETshifts_sameCol_fin_exclLast2Ma.pdf",sep=""), onefile=TRUE, width=4.5, height=8)
#pdf(file=paste("rateThroughTime_10trees_combo_bins1Ma_3rates_KPg_ALL_6ordsInOne_3higher.pdf",sep=""), onefile=TRUE, width=5, height=10)

#quartz(width=5, height=9)
op <- par(mfrow = c(6,1),
          oma = c(5,0,0,2) + 0.1, #c(bottom, left, top, right)
          mar = c(0,0,0,1) + 0.1)

# Now PLOT these per clade BAMM data...
####
ratetypes <- c("netdiv","speciation","extinction")
#intervalCols <- viridis(3) #c("blue","green","darkyellow")
#avgCols <- viridis(3) #c("blue","green","darkyellow")#rep("red",3)

intervalCols <- c(viridis(5)[1:4],"darkgoldenrod1","darkorange4") #c("blue","green","darkyellow")
#avgCols <- c(viridis(5)[1:4],"gold","darkgoldenrod3") #c("blue","green","darkyellow")#rep("red",3)

#lineTypes <- c(1, 5, 4)
lineTypes <- c(1, 1, 1)

for(v in 1:length(setNames)){ 
join10runs_ALL_v<-get(paste("rateMatrices_",setNames[v],sep=""))

intervalCol<-intervalCols[v]
#avgCol<-avgCols[v]
avgCols <- c(intervalCol,"black","white") #c("blue","green","darkyellow")#rep("red",3)

for(z in c(2,3,1)){
rmat<-join10runs_ALL_v
#rmat<-rateMatrices[[1]]

rates <- list((rmat$lambda - rmat$mu),rmat$lambda, rmat$mu)

rate<-rates[[z]]
rate<-rate[,1:(length(rate[1,])-2)]

ratelabel<-"rate"
ratetype<-ratetypes[z]
#intervalCol<-intervalCols[z]
lineType <- lineTypes[z]
avgCol<-avgCols[z]

intervals = seq(from = 0, to = 1, by = 0.01)
smooth = FALSE 
opacity = 0.01
plot = TRUE
cex.axis = 1
cex.lab = 1.3
lwd = 2
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
    rmat$times <- (max(rmat$times) - rmat$times)[1:(length(rmat$times)-2)]

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
  offset<- 0.22

if(z == 2){
        #plot.new()
        #par(mar = mar)          
    xMin <- -110 #maxTime
        xMax <- 0
        yMin <- 0 #min(poly[[1]][, 2])
        
        if(v==4){
          yMax <- 0.65 #max(poly[[1]][, 2])
        } else {
          yMax <- 0.4+offset
        }
        if (yMin >= 0) { yMin <- 0 }
        #plot.window(xlim = xLim, ylim = c(yMin, yMax))

    plot(x = -rmat$times, y = (avg+offset),lty=2,lwd=1,xlim=c(xMin, xMax),ylim=c(yMin, yMax),type="l", col="white",main=NULL,cex.axis=0.8,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",yaxt="n",xaxt="n",bty="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
    
    abline(v=-66,lwd=2,lty=2,col="dark grey")

    if(v==length(setNames)){ 
      axis(side=1,labels=TRUE, at=c(0,-20,-40,-60,-80), cex=0.8, line=0.5) 
      mtext(side = 1, text = "time before present", line = 3, cex = 0.8, adj=0.75)
      text(x=-66,y=0.22+offset,labels="K-Pg\nboundary",cex=1,font=2,pos=4, col=grey(0.5))
    }
    if(v==1){ 
      #legend(x=-50,y=0.39+offset,legend=c("speciation","extinction","net diversification"), col=intervalCol, lty=c(5,4,1), lwd=lwd, cex=0.8, seg.len=4, bty="n") 
        #legend(x=-20,y=0.39+offset,legend=c("speciation shift","extinction shift"), fill=c("black",grey(0.5,alpha=0.5)), border="black", cex=0.8, bty="n") 
        mtext(side = 2, text = "rates", line = -6.5, cex = 0.8)
        text(x=-85,y=0.35+offset,labels="(d) Order-wide rates",cex=1.2,font=2,pos=4)
    }
    #axis(at = c(round(1.2 * xMin), axTicks(1)), cex.axis = cex.axis, side = 1)
      if(v==6){ 
        axis(at = seq(0,0.3,by=0.1)+offset, labels=c(0,NA,NA,0.3),cex.axis = cex.axis, las = 1, side = 2, line=-6.5)
      } else {
        axis(at = seq(0,0.3,by=0.1)+offset, labels=c(0,NA,NA,0.3),cex.axis = cex.axis, las = 1, side = 2, line=-8)
      }
        mtext(side = 2, text = setNames2[v], line = -5, cex = 0.75, font=2, adj=0.3)
       }

        for (i in 1:length(poly)) {
            polygon(x = -poly[[i]][, 1], y = (poly[[i]][, 2]+offset), col = transparentColor(intervalCol, opacity), border = NA)
        }
    lines(x = -rmat$times, y = (avg+offset), lwd = lwd, col = avgCol, lty=lineType)
}

# Now get the per clade CoMET data...
##
#CoMET_cols<-viridis(3,alpha=0.5)
CoMET_cols<-c("black",grey(0.7,alpha=0.6))

#v=1
CoMET_dat<-get(paste("join4runs_ALL_",setNames[v],sep=""))
output<-CoMET_dat 
type<-"speciation shift times"
criticalPP <- output[[grep(strsplit(type, " ")[[1]][1], grep("CriticalPosteriorProbabilities", names(output), value = TRUE), value = TRUE)]]

type<-"speciation shift times"
colBars<-CoMET_cols[1]

thisOutput <- output[[type]]
meanThisOutput <- colMeans(thisOutput)
meanThisOutput[(length(meanThisOutput)-1):length(meanThisOutput)]<-NA # excludes last 2 Ma from the plot
par(new=TRUE, lwd=0.1)
#barplot(rev(meanThisOutput), space = 0, col = colBars, border = "black", main = "", ylab = "", xlab = "", xaxt = "n", yaxt = "n", xlim=c(-xMin,xMax),ylim = c(0,1))
barplot(rev(meanThisOutput), space = 0, col = colBars, border = "white", main = "", ylab = "", xlab = "", xaxt = "n", yaxt = "n", xlim=c(-xMin,-xMax),ylim = c(0,0.9))

type<-"extinction shift times"
colBars<-CoMET_cols[2]

thisOutput <- output[[type]]
meanThisOutput <- colMeans(thisOutput)
meanThisOutput[(length(meanThisOutput)-1):length(meanThisOutput)]<-NA # excludes last 2 Ma from the plot
par(new=TRUE, lwd=0.1)
#barplot(rev(meanThisOutput), space = 0, col = colBars, border = "black", main = "", ylab = "", xlab = "", xaxt = "n", yaxt = "n", xlim=c(-xMin,xMax),ylim = c(0,1))
barplot(rev(meanThisOutput), space = 0, col = colBars, border = "white", main = "", ylab = "", xlab = "", xaxt = "n", yaxt = "n", xlim=c(-xMin,-xMax),ylim = c(0,0.9))

if(v==6){ 
  axis(2, at = c(0,criticalPP[1]), labels = c(0,2), line= -6.5)#, srt=90)
} else {
  axis(2, at = c(0,criticalPP[1]), labels = c(0,2), line= -8)#, srt=90)
}
#axis(4, at = c(0,criticalPP[1]), labels = c(0,2), line= -0.5, srt=90)
if(v==1){ 
  mtext(side = 2, text = "BF of\nshift", line = -8.5, cex = 0.8, adj=0.5) 
  #mtext(side = 4, text = "Bayes factor", line = 2, cex = 0.8, adj=0) 
  }

#lines(x=c(85,-3), y=c(criticalPP[1],criticalPP[1]), lty = 1, col=grey(0.5))

}

dev.off()   








# PLOT
###		
# Fig2 - Panel 2 -- LTT of pre-KPg Orders + RTT of those Order-wide rates... 
  # Now WITHOUT the rate-shifts of tree-wide analyses (TreePar and CoMET-- save those for the supplement.)
#####

pdf(file=paste("LTT_95polygon_preKPgOrdersLTT_orderWideRTT_7part_no2Ma_wBAMM-only.pdf",sep=""), onefile=TRUE, width=4, height=16)
#pdf(file=paste("LTT_95polygon_preKPgOrdersLTT_orderWideRTT_7part_no2Ma_wBAMM-CoMET-TreePar.pdf",sep=""), onefile=TRUE, width=4, height=12)

xLim<-c(-110,25)
#quartz(width=4,height=8)
firstQuad<-0.2
next5<-(1-firstQuad)/5
#RTT_per6<-next6*(3/4)
#shifts_per6<-next6*(1/4)

#layout(matrix(c(1:13), 13, 1, byrow = TRUE), widths=rep(1,13),heights=c(firstQuad,rep(c(RTT_per6,shifts_per6),6)))
layout(matrix(c(1:6), 6, 1, byrow = TRUE), widths=rep(1,6),heights=c(firstQuad,rep(next5,5)))
op <- par(#mfrow = c(3,1),
          oma = c(5,4,2,0) + 0.1, #c(bottom, left, top, right)
          mar = c(1,0,1,1) + 0, xpd=TRUE)

LegRCex<-0.9 

# part 1
   # pre/at-KPg Placentalia 
  #last5ords<-c("gold","darkgoldenrod3", grey(0.2),grey(0.5),grey(0.8))#"darkgoldenrod4","brown","coral3")
  last5ords<-c("darkgoldenrod2", grey(0), grey(0.2),grey(0.5),grey(0.8))# , ,-- "darkgoldenrod4","brown","coral3")
  last5ordsAlpha<-c()
  for(q in 1:length(last5ords)){
   asRGB<-col2rgb(last5ords[q])/255
   asAlpha<-rgb(red=asRGB[[1]],green=asRGB[[2]],blue=asRGB[[3]],alpha=0.5)
   last5ordsAlpha[q]<-asAlpha
  }

  cols<-c(viridis(5,alpha=0.5)[c(1,2,3,4)],last5ordsAlpha)
  cols2<-c(viridis(5,alpha=1)[c(1,2,3,4)],last5ords)

  x<-ltt95_Placentalia
  x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
  plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",yaxt="n",xaxt="n", bty="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
  axis(side=1,labels=FALSE, at=c(0,-20,-40,-60,-80,-100))
  axis(side=2,labels=TRUE, at=c(1,5,50,500,5000), cex=0.8, line=0.25)
  mtext(side = 2, text = "(log) lineages", line = 3, cex = 1)

  text(x=-115,y=5000,labels="C) K-Pg orders",cex=1.2,font=2,pos=4)
  abline(v=-66,lwd=2,lty=2,col="dark grey")

  xVals<-c((-ordinal$hpd95_mins[1:6]+3),(-ordinal$hpd95_maxs[7:9]-3))
  yVals<-rev(exp(seq(log(50),log(3000),length.out=9)))
  posVals<-c(rep(4,6),2,2,2)
  cexLab1<-0.85

  for (i in 1:length(ltt95_preKPg)){
   x<-ltt95_preKPg[[i]]
   x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
  par(new=TRUE)
   polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
   lines(x[,2],x[,1],lty=1, lwd=1,col=cols2[i],type="l")
   lines(x[,4],x[,1],lty=1, lwd=1,col=cols2[i],type="l")
   #text(x=0,y=nrow(ltt95_preKPg[[i]][,])-1,labels=ordNames_preKPg[i],font=2, cex=cexLab1,offset=0,adj=c(0,0))#pos=4)
   text(x=xVals[i],y=yVals[i],labels=ordNames_preKPg[i],font=1, cex=cexLab1,offset=0,pos=posVals[i])
  }

  plotCI(add=TRUE, gap=0, x=rev(-ordinal$means),y=(exp(seq(log(50),log(3000),length.out=9))), ui=rev(-ordinal$hpd95_mins),li=rev(-ordinal$hpd95_maxs), cex=1, err="x",sfrac=0, col=rev(cols2),lwd=2, pch=16)

 # part 2
  for(q in 1:length(setNames)){
  par(xpd=NA) # to prevent clipping...
  #q<-2
    # LOAD the BAMM RTT results...
    join10runs_ALL<-get(paste("rateMatrices_",setNames[q],sep=""))
    
    # LOAD the BAMM shift results...
    orderBAMMshifts<-shiftDivDat[which(shiftDivDat$orderOfShift==setNames[q]),]
    
#   # # LOAD the CoMET RTT results... is same as the shift results...
#   # output<-get(paste("join4runs_ALL_",setNames[q],sep=""))
#   # type<-"speciation shift times"
#   # criticalPP <- output[[grep(strsplit(type, " ")[[1]][1], grep("CriticalPosteriorProbabilities", names(output), value = TRUE), value = TRUE)]]
#    
#    ## LOAD the TreePar RTT results...
#    #TP_RTT_all<-list(get(paste("TP_netdiv_10trees_",setNames[q],sep=""))[,length(grid):1], get(paste("TP_speciation_10trees_",setNames[q],sep=""))[,length(grid):1],get(paste("TP_extinction_10trees_",setNames[q],sep=""))[,length(grid):1])
#    
#    # LOAD the TreePar shift histos...
#    TP_resFIN<-get(paste("TP_resFIN_",setNames[q],sep=""))
#    shiftUp<-TP_resFIN$shiftMa[which(TP_resFIN$dir=="upShift" & TP_resFIN$shiftMa < -xMin)]
#    shiftDown<-TP_resFIN$shiftMa[which(TP_resFIN$dir=="downShift" & TP_resFIN$shiftMa < -xMin)]
#    
#    magUp<-abs(TP_resFIN$mag[which(TP_resFIN$dir=="upShift" & TP_resFIN$shiftMa < -xMin)])
#    magUp[which(magUp > 10)]<-10
#
#    magDown<-abs(TP_resFIN$mag[which(TP_resFIN$dir=="downShift" & TP_resFIN$shiftMa < -xMin)])
#    magDown[which(magDown > 10)]<-10
#    
#    times<-seq(from=xMin,to=xMax,by=1)
#    histUp<-hist(-shiftUp, breaks=times, plot=FALSE)$counts
#    names(histUp)<-times[1:(length(times)-1)]
#    histDown<-hist(-shiftDown, breaks=times, plot=FALSE)$counts
#    names(histDown)<-times[1:(length(times)-1)]


# now make the plot....
# ~~~~~~
  # BAMM RTT of global mammals
library(BAMMtools)
ratetypes <- c("netdiv","speciation","extinction")
intervalCol <- cols2[q] #grey(0.5) #viridis(3) #c("blue","green","darkyellow")
#avgCol <- grey(0,alpha=1) #viridis(3) #c("blue","green","darkyellow")#rep("red",3)
#lineTypes <- c(1,5,4)
lineTypes <- c(1,1,1)
avgCols <- c(grey(0.5),"black","white")

cexRLab<-0.65
cex.lab = 1

for(z in c(2,3,1)){
rmat<-join10runs_ALL
rates <- list((rmat$lambda - rmat$mu),rmat$lambda, rmat$mu)

rate<-rates[[z]][,1:(length(rates[[z]][1,])-2)] # exclude last 2 Ma

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
    rmat$times <- (max(rmat$times) - rmat$times)[1:length(rate[1,])] # exclude last 2 Ma

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
# offsetBamm<- 0.4
  offsetBamm<- 0.25

if(z == 2){
        #plot.new()
        #par(mar = mar)          
    xMin <- -110 #maxTime
        xMax <- 25
        yMin <- 0 #min(poly[[1]][, 2])
        yMax <- 0.4+offsetBamm #max(poly[[1]][, 2])
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
    pchs<-rep(21,length(orderBAMMshifts[,1]))
    #pchs[which(shiftDivDat$Indep==1)]<-16
    alphaCol<-0.7
    ptCols<-rep(rgb(1,0,0,alpha=alphaCol),length(orderBAMMshifts[,1]))
    ptCols[which(orderBAMMshifts$avgIncDec=="downShift")]<-rgb(0,0,1,alpha=alphaCol)
    ptCols[which(orderBAMMshifts$avgIncDec=="up-or-down")]<-grey(0.5,alpha=alphaCol)

    yVar<-orderBAMMshifts$avgFactor
    xVar<-(-orderBAMMshifts$mean)
    xLow<-(-orderBAMMshifts$lower)
    xHigh<-(-orderBAMMshifts$upper)

    bammShift_yMin<-1 # -2.5
    offsetBamm2<- 0.75

unloadNamespace("BAMMtools"); unloadNamespace("gplots")

    if(length(yVar)!=0) { 
    par(new=TRUE)
    plotCI(x=xVar, y=yVar+(3*offsetBamm2),ui=xLow,li=xHigh, xlim=c(xMin, xMax), ylim=c(bammShift_yMin,6), xlab="", ylab="", sfrac=0,err="x", lwd=2,col="black",pt.bg=ptCols, scol="black",pch=pchs,font.lab=2,cex.axis=0.95,cex.lab=1, xaxt="n", yaxt="n",bty="n", cex=1.5)
    text(x=xVar, y=yVar+(3*offsetBamm2), labels=as.vector(orderBAMMshifts$ID_label), pos=1, offset=0.7)
    
    axis(side=4,line= -5,at=c(1,2,3,4)+(3*offsetBamm2),labels=c(1,2,3,4))
    mtext(side = 4, text = "Shift factor", line = -3, cex = cexRLab, adj=0+offsetBamm2)
    text(x=-115, y= 6,labels = paste(setNames[q],sep=""), cex = 1.2,font=2,pos=4)
#legend(x=3,y=6,legend=c("up","down","up/down"), pt.bg=c(rgb(1,0,0,alpha=alphaCol),rgb(0,0,1,alpha=alphaCol),grey(0.5,alpha=alphaCol)), pch=21, cex=LegRCex, pt.cex=1.5, seg.len=4, bty="n", title="BAMM\nper-lineage\nshifts") 
    } else {
    text(x=-115, y= yMax,labels = paste(setNames[q],sep=""), cex = 1.2,font=2,pos=4)
    }

  # plot CoMET SHIFTS all (histos)
#amtUP<-0.35
#par(new=TRUE, lwd=0.1, fig=c(0,1,(firstQuad+(next6*q)+amtUP),(0.15+firstQuad+(next6*q)+amtUP)),oma = c(5,4,2,0) + 0.1,  mar = c(0.3,0,1,1) + 0.1, xpd=TRUE) # mar=c(0,0,0,0), mar=c(0,0,0,0)) # fig-- c(x1, x2, y1, y2)

# type<-"extinction shift times"
# colBars<-CoMET_cols[2]
# 
# thisOutput <- output[[type]]
# meanThisOutput <- colMeans(thisOutput)
# names(meanThisOutput)<-(length(thisOutput[1,])-1):0
# meanThisOutput[(length(meanThisOutput)-1):(length(meanThisOutput))]<-NA # excludes last 2 Ma from the plot
# par(new=TRUE, lwd=0.1)
# barplot(rev(meanThisOutput), space = 0, xpd=FALSE, col = colBars, border = "white", main = "", ylab = "", xlab = "", xaxt = "n", yaxt = "n", xlim=c(-xMin,-xMax),ylim = c(0,1))
# 
# type<-"speciation shift times"
# colBars<-CoMET_cols[1]
# 
# thisOutput <- output[[type]]
# meanThisOutput <- colMeans(thisOutput)
# names(meanThisOutput)<-(length(thisOutput[1,])-1):0
# meanThisOutput[(length(meanThisOutput)-1):(length(meanThisOutput))]<-NA # excludes last 2 Ma from the plot
# meanThisOutput_up<-meanThisOutput+CoMET_UP
#   par(new=TRUE, lwd=0.1)
#   barplot(rev(meanThisOutput), space = 0, xpd=FALSE, col = colBars, border = "white", main = "", ylab = "", xlab = "", xaxt = "n", yaxt = "n", xlim=c(-xMin,-xMax),ylim = c(0,1))
# 
# #axis(2, at = c(0,criticalPP[1],criticalPP[2]), labels = c(0,2,6), line=0.25)#, line= -4)
# axis(4, at = c(0,criticalPP[1]), labels = c(0,2), line=-5) #0.25)#, line= -4)
# mtext(side = 4, text = "2 ln BF", line = -1, cex = cexRLab, adj=0) 

#legend(x=-3,y=0.4,legend=c("speciation","extinction"), pt.bg=c(rgb(1,0,0,alpha=alphaCol),rgb(0,0,1,alpha=alphaCol)), pch=22, cex=LegRCex, pt.cex=1.5, seg.len=4, bty="n", title="CoMET\ntree-wide\nshifts") 

#  # plot TreePar MAGNITUDES of shifts all (akin to the Bamm shifts...)
#par(new=TRUE, xpd=NA)
#
#plot(x= -shiftUp, y= magUp, xlim=c(xMin, xMax), ylim=c(0,20), bty="n", xlab="", ylab="", pch=21,col="black",xaxt="n", yaxt="n",cex=1.5)
#
#axis(side=4,at=seq(from=0,to=7,length.out=4),labels=TRUE, line=-5) # seq(from=0,to=max(c(histDown,histUp)),length.out=3)
#
#points(x= -shiftDown, y= magDown)
#
#
#barplot(-histUp, space = 0, xpd=NA, main="", col=rgb(1,0,0,alpha=0.5), border="white",lwd=0.5, ylim=c(0,10),yaxt="n",xaxt="n")
#barplot(add=TRUE,-histDown, space = 0, xpd=NA,  col=rgb(0,0,1,alpha=0.5), border="white", lwd=0.5, xaxt="n", yaxt="n")
#
#axis(side=4,at=seq(from=0,to=-4,length.out=5),labels=c(0,NA,0.2,NA,0.4), line=-5) # seq(from=0,to=max(c(histDown,histUp)),length.out=3)
#mtext(side=4,text="Shifts per tree", line=-3, cex=cexRLab, adj=-2)
#
#  # plot TreePar SHIFTS all (histos)
#par(new=TRUE)
#barplot(-histUp, space = 0, xpd=NA, main="", col=rgb(1,0,0,alpha=0.5), border="white",lwd=0.5, ylim=c(0,10),yaxt="n",xaxt="n")
#barplot(add=TRUE,-histDown, space = 0, xpd=NA,  col=rgb(0,0,1,alpha=0.5), border="white", lwd=0.5, xaxt="n", yaxt="n")
#
#axis(side=4,at=seq(from=0,to=-4,length.out=5),labels=c(0,NA,0.2,NA,0.4), line=-5) # seq(from=0,to=max(c(histDown,histUp)),length.out=3)
#mtext(side=4,text="Shifts per tree", line=-3, cex=cexRLab, adj=-2)

#legend(x=120,y=7,legend=c("up","down"), pt.bg=c(rgb(1,0,0,alpha=alphaCol),rgb(0,0,1,alpha=alphaCol)), pch=22, cex=LegRCex, pt.cex=1.5, seg.len=4, bty="n", title="TreePar\ntree-wide\nshifts") 

# plot(density(-TP_resFIN$shiftMa), type="n",xlim=c(xMin,xMax), yaxt="n",xaxt="n", bty="n",main="",ylab="",xlab="")
# rug(-shiftDown-1, col=grey(0.3), lwd=2,line=0.5)
# rug(-shiftUp-1, col=rgb(1,0,0,alpha=0.5), lwd=2, line=1.5)
# lines(x=c(-66,-66),y=c(-2,1),col=grey(0.5),lty=2,lwd=3)
# axis(side=1,at=seq(from=-100,to=0,by=20),labels= seq(from=-100,to=0,by=20), line=1.5)

} # end 5 order loop

    axis(side=1, at= rev(c(0,-20,-40,-60,-80,-100)), labels=rev(-c(0,-20,-40,-60,-80,-100)), cex=cex.axis, line=0.5)
    mtext(side = 1, text = "time before present (Ma)", line = 3.5, cex = cex.lab, adj=0.4)

dev.off()






##
# 3-part-- two LTT (pre and post) and then RTT global mammal
###
pdf(file=paste("LTT_95polygon_Placentals_pre-vs-postKPg_withBAMM-RTT.pdf",sep=""), onefile=TRUE, width=5, height=9)

xLim<-c(-110,25)
#quartz(width=5,height=6)
#layout(matrix(c(1:2), 1, 2, byrow = TRUE))
op <- par(mfrow = c(3,1),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)

cexLab1<-0.6
cexLab2<-0.9

# pre/at-KPg Placentalia 
cols<-viridis(length(ordNames_preKPg),alpha=0.3)
cols2<-viridis(length(ordNames_preKPg),alpha=1)

x<-ltt95_Placentalia
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",yaxt="n",xaxt="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
axis(side=1,labels=FALSE, at=c(0,-20,-40,-60,-80,-100))
axis(side=2,labels=TRUE, at=c(1,5,50,500,5000), cex=0.8)

text(x=-115,y=3000,labels="(a)",cex=1,font=2,pos=4)

abline(v=-66,lwd=2,lty=2,col="dark grey")
for (i in 1:length(ltt95_preKPg)){
	x<-ltt95_preKPg[[i]]
	x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
par(new=TRUE)
	polygon(c(x[,2], rev(x[,4])), c(x[,1], rev(x[,1])),col=cols[i], border = NA) #awesome.
	lines(x[,2],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	lines(x[,4],x[,1],lty=1, lwd=1,col=grey(0.4,alpha=0.1),type="l")
	text(x=0,y=nrow(ltt95_preKPg[[i]][,])-1,labels=ordNames_preKPg[i],font=2, cex=cexLab1,offset=0,adj=c(0,0))#pos=4)
}

plotCI(add=TRUE, gap=0, x=rev(-ordinal$means),y=(exp(seq(log(50),log(3000),length.out=9))), ui=rev(-ordinal$hpd95_mins),li=rev(-ordinal$hpd95_maxs), cex=1, err="x",sfrac=0, col=rev(cols2),lwd=2, pch=16)


# post-KPg Placentalia 
cols<-viridis(length(ordNames_postKPg),alpha=0.3)
cols2<-viridis(length(ordNames_postKPg),alpha=1)
colsSuper<-c("steelblue1","steelblue2","steelblue3","steelblue4") #rich.colors(4)

x<-ltt95_Placentalia
x[,2:4]<-x[,2:4]-max(x[,2:4]) # makes negative
plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,type="l", col="white",log="y", main=NULL,cex.axis=0.8,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",yaxt="n",xaxt="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
axis(side=1,labels=FALSE, at=c(0,-20,-40,-60,-80,-100), cex=0.8)
axis(side=2,labels=TRUE, at=c(1,5,50,500,5000), cex=0.8)
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

# 

#dev.off()

ratetypes <- c("netdiv","speciation","extinction")
intervalCols <- viridis(3) #c("blue","green","darkyellow")
avgCols <- viridis(3) #c("blue","green","darkyellow")#rep("red",3)

for(z in length(ratetypes):1){
rmat<-join10runs_ALL
#rmat<-rateMatrices[[1]]

rates <- list((rmat$lambda - rmat$mu),rmat$lambda, rmat$mu)

rate<-rates[[z]]
ratelabel<-"rate"
ratetype<-ratetypes[z]
intervalCol<-intervalCols[z]
avgCol<-avgCols[z]

intervals = seq(from = 0, to = 1, by = 0.01)
smooth = FALSE 
opacity = 0.01
plot = TRUE
cex.axis = 1
cex.lab = 1.3
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
    rmat$times <- rmat$times[which(nanCol == FALSE)]
    rmat$times <- max(rmat$times) - rmat$times

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

if(z == 3){
       	#plot.new()
        #par(mar = mar)          
 		#xMin <- 150 #maxTime
        #xMax <- 0
        yMin <- min(poly[[1]][, 2])
        yMax <- 0.3 #max(poly[[1]][, 2])
        if (yMin >= 0) { yMin <- 0 }
        #plot.window(xlim = xLim, ylim = c(yMin, yMax))

		plot(x[,3],x[,1],lty=2,lwd=1,xlim=xLim,,ylim=c(yMin, yMax),type="l", col="white",main=NULL,cex.axis=0.8,cex.lab=1,font.axis=1, font.lab=2, xlab="", ylab="",yaxt="n",xaxt="n") #xlab="Time before present (Ma)", ylab="(log) Number of lineages")
		
		axis(side=1,labels=TRUE, at=c(0,-20,-40,-60,-80,-100), cex=0.8)
		#axis(at = c(round(1.2 * xMin), axTicks(1)), cex.axis = cex.axis, side = 1)
        axis(at = seq(0,0.3,by=0.05), labels=c(0,NA,0.1,NA,0.2,NA,0.3),cex.axis = cex.axis, las = 1, side = 2)
        
        mtext(side = 1, text = "time before present", line = xline, cex = cex.lab)
        mtext(side = 2, text = ratelabel, line = yline, cex = cex.lab)

		abline(v=-66,lwd=2,lty=2,col="dark grey")
       }

        for (i in 1:length(poly)) {
            polygon(x = -poly[[i]][, 1], y = poly[[i]][, 2], col = transparentColor(intervalCol, opacity), border = NA)
        }
    lines(x = -rmat$time, y = avg, lwd = lwd, col = avgCol)
}

legend(x=-112,y=0.29,legend=c("net diversification","speciation","extinction"), col=viridis(3), lwd=lwd)

dev.off()












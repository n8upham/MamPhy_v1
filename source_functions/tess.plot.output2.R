tess.plot.output2<-function (output, fig.types = c("speciation rates", "speciation shift times", 
    "speciation Bayes factors", "extinction rates", "extinction shift times", 
    "extinction Bayes factors", "net-diversification rates", 
    "relative-extinction rates", "mass extinction times", "mass extinction Bayes factors"), 
    xlab = "million years ago", col = NULL, col.alpha = 50, treeAlpha = 0.5, xaxt = "n", 
    yaxt = "s", pch = 19, plot.tree = FALSE, ...) 
{
    validFigTypes <- c("speciation rates", "speciation shift times", 
        "speciation Bayes factors", "extinction rates", "extinction shift times", 
        "extinction Bayes factors", "net-diversification rates", 
        "relative-extinction rates", "mass extinction times", 
        "mass extinction Bayes factors")
    invalidFigTypes <- fig.types[!fig.types %in% validFigTypes]
    if (length(invalidFigTypes) > 0) {
        stop("\nThe following figure types are invalid: ", paste(invalidFigTypes, 
            collapse = ", "), ".", "\nValid options are: ", paste(validFigTypes, 
            collapse = ", "), ".")
    }
    if (is.null(col)) {
        col <- c(`speciation rates` = "#984EA3", `speciation shift times` = "#984EA3", 
            `speciation Bayes factors` = "#984EA3", `extinction rates` = "#E41A1C", 
            `extinction shift times` = "#E41A1C", `extinction Bayes factors` = "#E41A1C", 
            `net-diversification rates` = "#377EB8", `relative-extinction rates` = "#FF7F00", 
            `mass extinction times` = "#4DAF4A", `mass extinction Bayes factors` = "#4DAF4A")
    }
    else {
        names(col) <- fig.types
    }
    treeAge <- output[["tree"]]$root.age
    numIntervals <- length(output$intervals) - 1
    plotAt <- 0:numIntervals
    intervalSize <- treeAge/numIntervals
    labels <- pretty(c(0, treeAge))
    labelsAt <- numIntervals - (labels/intervalSize)
    
    for (type in fig.types) {
        if (grepl("times", type)) {
            thisOutput <- output[[type]]
            meanThisOutput <- colMeans(thisOutput)
            criticalPP<- output[[grep(strsplit(type, " ")[[1]][1], grep("CriticalPosteriorProbabilities", names(output), value = TRUE), value = TRUE)]]
            if (plot.tree) {
                plot(output$tree, show.tip.label = FALSE, edge.col = rgb(0, 
                  0, 0, treeAlpha), x.lim = c(0, treeAge))
                par(new = TRUE)
            }
            barplot(meanThisOutput, space = 0, xaxt = xaxt, col = col[type], 
                border = col[type], main = type, ylab = "posterior probability", 
                xlab = xlab, ylim = c(0, 1), ...)
            abline(h = criticalPP, lty = 2, ...)
            axis(4, at = criticalPP, labels = 2 * log(output$criticalBayesFactors), 
                las = 1, tick = FALSE, line = -0.5)
            axis(1, at = labelsAt, labels = labels)
            box()
        }
        else if (grepl("Bayes factors", type)) {
            thisOutput <- output[[type]]
            ylim <- range(c(thisOutput, -10, 10), finite = TRUE)
            if (plot.tree) {
                plot(output$tree, show.tip.label = FALSE, edge.col = rgb(0, 
                  0, 0, 0.1), x.lim = c(0, treeAge))
                par(new = TRUE)
            }
            plot(x = plotAt[-1] - diff(plotAt[1:2])/2, y = thisOutput, 
                type = "p", xaxt = xaxt, col = col[type], ylab = "Bayes factors", 
                main = type, xlab = xlab, ylim = ylim, xlim = range(plotAt), 
                pch = pch, ...)
            abline(h = 2 * log(output$criticalBayesFactors), 
                lty = 2, ...)
            axis(4, at = 2 * log(output$criticalBayesFactors), 
                las = 1, tick = FALSE, line = -0.5)
            axis(1, at = labelsAt, labels = labels)
        }
        else {
            thisOutput <- output[[type]]
            meanThisOutput <- colMeans(thisOutput)
            quantilesThisOutput <- apply(thisOutput, 2, quantile, 
                prob = c(0.025, 0.975))
            if (type %in% c("speciation rates", "extinction rates")) {
                quantilesSpeciation <- apply(output[["speciation rates"]], 
                  2, quantile, prob = c(0.025, 0.975))
                quantilesExtinction <- apply(output[["extinction rates"]], 
                  2, quantile, prob = c(0.025, 0.975))
                ylim <- c(0, max(quantilesSpeciation, quantilesExtinction))
            }
            else {
                ylim <- c(0, max(quantilesThisOutput))
            }
            if (plot.tree) {
                plot(output$tree, show.tip.label = FALSE, edge.col = rgb(0, 
                  0, 0, 0.1), x.lim = c(0, treeAge))
                par(new = TRUE)
            }
            plot(x = plotAt, y = c(meanThisOutput[1], meanThisOutput), 
                type = "l", ylim = ylim, xaxt = xaxt, col = col[type], 
                ylab = "rate", main = type, xlab = xlab, ...)
            polygon(x = c(0:ncol(quantilesThisOutput), ncol(quantilesThisOutput):0), 
                y = c(c(quantilesThisOutput[1, 1], quantilesThisOutput[1, 
                  ]), rev(c(quantilesThisOutput[2, 1], quantilesThisOutput[2, 
                  ]))), border = NA, col = paste(col[type], col.alpha, 
                  sep = ""))
            axis(1, at = labelsAt, labels = labels)
        }
    }
}
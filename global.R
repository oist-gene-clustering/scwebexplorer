##########################################
#                                        #
#      SHINY APPLICATION: UTILS          #
#                                        #
##########################################

library(shiny)
library(jsonlite)
library(scater, quietly = TRUE)
library(rPython)
library(RColorBrewer)
require(genefilter)
library(gplots)
library(NMF)
require(calibrate)
library(MAST)
library(shinyBS)
library(shinythemes)
library(shinydashboard)
library(shinyjs)
library(DT)
library(shinycssloaders)
library(shinyDND)
library(RColorBrewer)
library(ggfortify)
library(ggplot2)
library(ggbiplot)
library(plotrix)
library(shinyWidgets)

############
# DATASETS #
############

load("utils.Rdata")

##########
# LAYOUT #
##########

#_______________________#
#      Loading page     #
#_______________________#

#' Loading data animation
#'
#' Loading data animation --side effect.
#'
#' @export
load_data <- function() {
  Sys.sleep(4)
  hide("loading_page", anim=T, animType="fade")
  show("main_content")
}

##############
# PARAMETERS #
##############

ncores=detectCores()-1

## Hides error messages from the application    ##
## Can display the trace of Shiny requests      ##
## Controls the number of recursive operations  ##
## Number of cores detectCores()-1              ##
options(shiny.sanitize.errors=TRUE, 
        shiny.trace=FALSE, 
        expressions=5e5,
        mc.cores=ncores)

## Dimensions of cell rectangles ##
topside = 100
leftside = 50

## Dimension of buckets ##
n=max(sapply(scesets, function(sceset) dim(sceset)[2]))

msg="Please run MAST computation on gene and cell selection"

########################
# SANITIZING FUNCTIONS #
########################

#' Return NULL or result
#'
#' Returns NULL iff result is character(0), numeric(0), NULL, etc.
#'
#' @param e element to test
#' @return e if not character(0), numeric(0), NULL, etc. else NULL
#' 
#' @export
rn <- function(e) if (length(e) == 0) return(NULL) else return(e)

#' Convert string input from cell group selection to list
#'
#' Converts string input from cell group selection to list.
#'
#' @param v string
#' @return list
#' 
#' @export
tolist <- function(v) return(lapply(v, function(e) paste(" ", e, " ", sep="")))

####################
# STRING FUNCTIONS #
####################

#' Filter list with pattern
#'
#' Returns the list which elements matching the pattern have been deleted.
#'
#' @param pattern string pattern
#' @param ls string list
#' @return the sublist of w which elements matching the pattern have been deleted.
#' 
#' @export
filterList <- function(pattern, ls) return(rn(ls[!grepl(pattern, ls)]))

#' Filter list with patterns
#'
#' Returns the list which elements matching at least one of the patterns have been deleted.
#'
#' @param patternList list of string patterns
#' @param ls string list
#' @return the list which elements matching at least one of the patterns have been deleted.
#' 
#' @export
filterMultiList <- function(patternList, ls) {
  keep <- lapply(patternList, function(pattern) filterList(pattern, ls))
  for (keepLs in keep) ls <- intersect(ls, keepLs)
  return(ls)
}

#' Filter list with pattern -inverse
#'
#' Returns the list which elements NOT matching the pattern have been deleted.
#'
#' @param pattern pattern
#' @param ls string list
#' @return the list which elements NOT matching the pattern have been deleted.
#' 
#' @export
keepList <- function(pattern, ls) return(rn(setdiff(ls, filterList(pattern, ls))))

#' Clean list from characters
#'
#' Cleans the list.
#'
#' @param ls string list
#' @return cleaned list
#' 
#' @export
clean <- function(ls) {
  if (is.null(ls)) return(NULL)
  tmp <- strsplit(ls, " ")[[1]]
  tmp <- tmp[sapply(tmp, function(x) !(x == "" | x == "\n"))]
  tmp <- sapply(tmp, function(e) strsplit(e, "\n")[[1]])
  return(unlist(tmp))
}

##############################
# SCESET SELECTION FUNCTIONS #
##############################

#' Get gene expression data from SCESet
#'
#' Gets gene expression data from an input SCESet.
#'
#' @param sceset the SCESet associated with the gene expression matrix of the new organism
#' @return gene expression data
#' 
#' @importFrom scater counts
#' @importFrom scater fpkm
#' @importFrom scater tpm
#' @importFrom scater cpm
#' @importFrom scater exprs
#' 
#' @export
getData <- function(sceset) {
	if (!is.null(counts(sceset))) return(counts(sceset))
	if (!is.null(fpkm(sceset))) return(fpkm(sceset))
	if (!is.null(tpm(sceset))) return(tpm(sceset))
	if (!is.null(cpm(sceset))) return(cpm(sceset))
	if (!is.null(exprs(sceset))) return(exprs(sceset))
}

#' Return feature list
#'
#' Returns the list of features in the input SCESet.
#'
#' @param sceset the input SCESet
#' @return the list of features of the input SCESet.
#' 
#' @export
getFeatureList <- function(sceset) 
  if (is.null(sceset)) return(NULL) else return(rn(filterMultiList(houseKeepGenes, row.names(sceset))))

#' Return SCESet associated with selected organism
#'
#' Returns the SCESet associated with the selected organism -excluding housekeeping genes.
#'
#' @param organism string name of the organism
#' @return the SCESet.
#' 
#' @export
getSceset <- function(organism) 
  if (is.null(organism)) return(NULL) else return(scesets[[organism]][getFeatureList(scesets[[organism]]), ])

#' Return SCESet associated with selected organism
#'
#' Returns the SCESet associated with the selected organism -including housekeeping genes.
#'
#' @param organism string name of the organism
#' @return the SCESet.
#' 
#' @export
getScesetR <- function(organism) 
  if (is.null(organism)) return(NULL) else return(scesets[[organism]])

#' Return sub-SCESet restricted to selected cells and genes and sum
#'
#' Returns the sub-SCESet restricted to the selected cells and genes, and 
#' sum expression values for same cell replicates.
#'
#' @param organism input organism
#' @param cells list of selected cell names
#' @param genes list of selected gene names
#' @param hg boolean set to TRUE iff housekeeping genes should be taken into account for the sum
#' @return the corresponding sub SCESet.
#' 
#' @export
selectGroupInSceset <- function(organism, cells=NULL, genes=NULL, hg=F) {
  if (is.null(organism)) return(NULL)
  sceset <- groupedScesets[[organism]]
  if (!is.null(genes) & hg) genes <- rownames(sceset)[as.character(rownames(sceset)) %in% genes]
  if (!is.null(genes) & !hg) genes <- rownames(sceset)[as.character(rownames(sceset)) %in% filterMultiList(houseKeepGenes, genes)]
  if (is.null(genes) & !hg) genes <- getFeatureList(sceset)
  if (!is.null(genes)) sceset <- sceset[genes, ]
  if (!is.null(cells)) sceset <- sceset[, colnames(sceset)[colnames(sceset) %in% cells]]
  return(sceset)
}

#######################
# OPERATION FUNCTIONS #
#######################

#________________________#
#      Pattern Plots     #
#________________________#

#' Create Pattern plot
#'
#' Creates pattern plot on selected genes and cells.
#'
#' @param genes selected genes in feature list
#' @param cells cell group of interest
#' @param oo considered organism
#' @param allGenes feature list
#' @param thresholdType type of threshold to consider for pattern colouring
#' @param threshold the user-chosen threshold -for binary threshold
#' @param thresholdDF the data frame used for graded threshold
#' @return gene expression pattern plotting for @cellsPattern
#' 
#' @importFrom grDevices heat.colors
#' @importFrom graphics plot
#' @importFrom graphics rect
#' @importFrom graphics text
#' @importFrom graphics legend
#' @importFrom graphics axis
#' @importFrom shiny withProgress
#' @importFrom shiny incProgress
#' 
#' @export
patternPlot <- function(genes, cells, oo, allGenes, thresholdType, threshold, thresholdDF, nb=1) {
  colours <- heat.colors(11)
  genes <- sort(genes)
  lsCell <- clean(cells)
  sceset <- selectGroupInSceset(oo, lsCell, genes)
  if (is.null(sceset)) return(NULL)
  else {
      res <- getData(sceset)
      nbCells <- dim(res)[2]
      nbGenes <- length(genes)
      withProgress(message = paste("Computing pattern", nb), value = 0, {
        incProgress(1/4, detail = paste("Computation", nb))
        ansArray <- drawCell(sort(genes), lsCell, 
                             data.matrix(res), thresholdType, threshold, thresholdDF, colours)
        if (is.null(ansArray)) return(NULL)
        incProgress(2/4, detail = paste("Graphics", nb))
        pp <- plot(c(topside-110, (nbCells+1)*topside+10), c(leftside-10, (nbGenes+1)*leftside+10),
             main = paste("Expression levels of selected genes in cell group", nb),
             type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
        rr <- rect(ansArray[[1]], ansArray[[2]], ansArray[[3]], ansArray[[4]], col=ansArray[[5]])
        incProgress(3/4, detail = paste("Labels", nb))
        tt <- text(ansArray[[6]], y=ansArray[[7]], labels=genes)
        if (thresholdType == "Binary") legend("topleft", legend = c("Gene ON", "Gene OFF"), 
                                                          fill = c(colours[10], "white"))
        else legend("topleft", legend = c("0-25%", "25-50%", "50-75%", "75-100%"),
                    fill = c("white", colours[11], colours[10], colours[9]))
        incProgress(1, detail = paste("Axis", nb))
        aa <- axis(side=1, at=unique(ansArray[[6]]), labels=lsCell)
        return(list(pp, rr, tt, aa))
      })
  }
}

#' Create Celegans pattern plot
#'
#' Draws Celegans plot for gene expression pattern visualization.
#'
#' @param genes selected genes in feature list
#' @param cells cell group of interest
#' @param oo considered organism
#' @param allGenes feature list
#' @param thresholdType type of threshold to consider for pattern colouring
#' @param threshold the user-chosen threshold -for binary threshold
#' @param thresholdDF the data frame used for graded threshold
#' @return gene expression pattern plotting for @cellsPattern
#' 
#' @importFrom grDevices heat.colors
#' @importFrom graphics plot
#' @importFrom graphics rect
#' @importFrom graphics text
#' @importFrom graphics legend
#' @importFrom graphics axis
#' @importFrom plotrix draw.circle
#' @importFrom shiny withProgress
#' @importFrom shiny incProgress
#' 
#' @export
patternPlotCelegans <- function(genes, cells, oo, allGenes, thresholdType, threshold, thresholdDF, nb=1) {
  colours <- heat.colors(11)
  #Import palette from RcolorBrewer
  #celltypes <- c("red", "blue", "orange", "purple", "black", "green", "pink", "cyan")
  genes <- sort(genes)
  lsCell <- clean(cells)
  sceset <- selectGroupInSceset(oo, lsCell, genes)
  if (is.null(sceset)) return(NULL)
  else {
    res <- getData(sceset)
    nbGenes <- length(genes)
    withProgress(message = paste("Computing pattern", nb), value = 0, {
      # pp <- plot(c(topside-110, (nbCells+1)*topside+10), c(leftside-10, (nbGenes+1)*leftside+10),
      #            main = paste("Expression levels of selected genes in cell group", nb),
      #            type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
      # rr <- rect(ansArray[[1]], ansArray[[2]], ansArray[[3]], ansArray[[4]], col = ansArray[[5]])
      # tt <- text(ansArray[[6]], y = ansArray[[7]], labels = genes)
      # aa <- axis(side=1, at=unique(ansArray[[6]]), labels=lsCell)
      # return(list(pp, rr, tt, aa))
    })
  }
}

makeRectangles <- function(labels, ncol, nrow, center) {}

#' Create Ciona Pattern plot
#'
#' Draws Ciona plot for gene expression pattern visualization..
#'
#' @param genes selected genes in feature list
#' @param cells cell group of interest
#' @param oo considered organism
#' @param allGenes feature list
#' @param thresholdType type of threshold to consider for pattern colouring
#' @param threshold the user-chosen threshold -for binary threshold
#' @param thresholdDF the data frame used for graded threshold
#' @return gene expression pattern plotting for @cellsPattern
#' 
#' @importFrom grDevices heat.colors
#' @importFrom graphics plot
#' @importFrom graphics rect
#' @importFrom graphics text
#' @importFrom graphics legend
#' @importFrom graphics axis
#' @importFrom plotrix draw.circle
#' @importFrom shiny withProgress
#' @importFrom shiny incProgress
#' 
#' @export
# patternPlotCiona <- function(genes, cells, oo, allGenes, thresholdType, threshold, thresholdDF, nb=1) {
#   colours <- heat.colors(11)
#   celltypes <- c("red", "blue", "orange", "purple", "black", "green", "pink", "cyan")
#   genes <- sort(genes)
#   lsCell <- clean(cells)
#   sceset <- selectGroupInSceset(oo, lsCell, genes)
#   if (is.null(sceset)) return(NULL)
#   else {
#     res <- getData(sceset)
#     nbGenes <- length(genes)
#     withProgress(message = paste("Computing pattern", nb), value = 0, {
#       incProgress(1/4, detail = paste("Computation", nb))
#       pp <- plot(c(-1, 16), c(0, 5), main = paste("Expression levels of selected genes in cell group", nb),
#                  type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
#       rr <- {
#         draw.circle(0.8, 2.5, 1, border=celltypes[1], col="white")
#         draw.circle(5.2, 2.5, 1, border=celltypes[1], col="white")
#         draw.circle(2.1, 1.8, 1.1, border=celltypes[3], col="white")
#         draw.circle(4.1, 1.8, 1.1, border=celltypes[3], col="white")
#         draw.circle(2, 2.5, 1.2, border=celltypes[2], col="white")
#         draw.circle(4, 2.5, 1.2, border=celltypes[2], col="white")
#         draw.circle(2.3, 3.2, 1, border=celltypes[4], col="white")
#         draw.circle(3.9, 3.2, 1, border=celltypes[4], col="white")
#         
#         draw.circle(9.8, 2.8, 1, border=celltypes[5], col="white")
#         draw.circle(14.2, 2.8, 1, border=celltypes[5], col="white")
#         draw.circle(11.1, 3.3, 1.1, border=celltypes[6], col="white")
#         draw.circle(13.1, 3.3, 1.1, border=celltypes[6], col="white")
#         draw.circle(11.3, 1.8, 0.8, border=celltypes[8], col="white")
#         draw.circle(12.9, 1.8, 0.8, border=celltypes[8], col="white")
#         draw.circle(11, 2.5, 1.2, border=celltypes[7], col="white")
#         draw.circle(13, 2.5, 1.2, border=celltypes[7], col="white")
#       }
#       incProgress(2/4, detail = paste("Drawing cells...", nb))
#       tt <- text(c(3, 12, rep(-0.7, 4), rep(15.7, 4)), 
#                  c(4.5, 4.5, 1.8, 2.35, 2.65, 3.2, 1.8, 2.8, 2.5, 3.3), 
#                  label=c("Animal view", "Vegetal view", "b5.4", "b5.3", "a5.4", "a5.3", "B5.2", "A5.2", "B5.1", "A5.1"),
#                  col=c("black", "black", celltypes[4], celltypes[1], celltypes[2], celltypes[3], 
#                        celltypes[8], celltypes[5], celltypes[7], celltypes[6]))
#       
#       ansArray <- drawCell(sort(genes), lsCell,
#                            data.matrix(res), thresholdType, threshold, thresholdDF, colours)
#       if (is.null(ansArray)) return(NULL)
#       incProgress(2/4, detail = paste("Graphics", nb))
#       pp <- plot(c(topside-110, (nbCells+1)*topside+10), c(leftside-10, (nbGenes+1)*leftside+10),
#                  main = paste("Expression levels of selected genes in cell group", nb),
#                  type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
#       rr <- rect(ansArray[[1]], ansArray[[2]], ansArray[[3]], ansArray[[4]], col = ansArray[[5]])
#       incProgress(3/4, detail = paste("Labels", nb))
#       tt <- text(ansArray[[6]], y = ansArray[[7]], labels = genes)
#       if (thresholdType == "Binary") legend("topleft", legend = c("Gene ON", "Gene OFF"),
#                                             fill = c(colours[10], "white"))
#       else legend("topleft", legend = c("0-25%", "25-50%", "50-75%", "75-100%"),
#                   fill = c("white", colours[11], colours[10], colours[9]))
#       incProgress(1, detail = paste("Axis", nb))
#       aa <- axis(side=1, at=unique(ansArray[[6]]), labels=lsCell)
#       return(list(pp, rr, tt, aa))
#     })
#   }
# }

#________________________#
#      PCA Plots         #
#________________________#

#' Create PCA plot
#'
#' Creates PCA plot on cells with the given feature list -all genes by default.
#'
#' @param oo string of input organism name
#' @param main plot title
#' @param genes list of features with which PCA should be applied
#' @return resulting PCA plot
#' 
#' @importFrom ggfortify autoplot
#' @importFrom ggplot2 aes
#' @importFrom RColorBrewer scale_colour_brewer
#' @importFrom ggplot2 geom_label
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 theme
#' @importFrom stats prcomp
#' 
#' @export
customPCAPlot <- function(oo, main, genes=NULL) {
  data <- pcaPlots[[oo]][["data"]]
  df <- pcaPlots[[oo]][["df"]]
  if (!is.null(genes)) data <- prcomp(t(getData(getSceset[[oo]]))[, genes])
  gp = autoplot(data, data=df, label=F, shape=F)
  gp = gp + geom_label(aes(label=label, colour=colour)) + scale_colour_brewer(palette="Set1")
  gp = gp + ggtitle(main) + theme(plot.title=element_text(hjust=0.5))
  return(gp)
}

#' Create PCA biplot
#'
#' Creates PCA biplot -with loadings- on cells with the given feature list -all genes by default.
#'
#' @param oo string of input organism name
#' @param main plot title
#' @param genes list of features with which PCA should be applied
#' @return resulting PCA biplot with the first two loadings/eigenvectors
#' 
#' @importFrom ggbiplot ggbiplot
#' @importFrom RColorBrewer scale_colour_brewer
#' @importFrom ggplot2 geom_label
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 aes
#' @importFrom stats prcomp
#' @importFrom shiny withProgress
#' @importFrom shiny incProgress
#' @importFrom ggfortify autoplot
#' 
#' @export
customPCABiPlot <- function(oo, main, genes=NULL) {
  withProgress(message="Create biplot", value=0, {
    incProgress(1/4, detail = "Fetching data...")
    sceset <- getSceset(oo)
    if (!is.null(genes)) sceset <- sceset[genes, ]
    incProgress(2/4, detail = "Performing PCA...")
    trimmedData <- getData(sceset)[apply(getData(sceset), 1, function(x) !(var(x) == 0)), ]
    trimmedData <- log(t(trimmedData)+1)
    if (any(apply(trimmedData, 2, function(x) (var(x) == 0))) > 0) return(NULL)
    df <- pcaPlots[[oo]]$df
    df.pca <- prcomp(trimmedData, scale=T)
    incProgress(3/4, detail = "Building biplot...")
    gp = autoplot(df.pca, data=df, label=F, shape=F, loadings = TRUE, 
                  loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3)
    gp = gp + geom_label(aes(label=label, colour=colour))
    incProgress(1, detail = "Done!")
    gp <- gp + ggtitle(main) + theme(plot.title=element_text(hjust=0.5))
    return(gp)
  })
}

#________________________#
#      Cell drawing      #
#________________________#

#' Draw cell patterns 
#'
#' Draws the cell patterns, as rectangles, associated to the input list of conditions and features -the feature in one
#' cell is coloured iff the feature is considered "active" in the cell.
#'
#' @param lsGene list of features to draw
#' @param lsCell list of conditions to draw
#' @param res array containing the expression level value for each pair cell, gene
#' @param thresholdType can be "binary" or "graded" - user-selected
#' @param thresholdBinary the activation value of the selected binary threshold
#' @param threshold a dataframe containing the activation bounds for a graded visualization
#' @param colours set of colours to use for the drawing
#' @return a list containing the arrays of x coordinates, y coordinates, rectangle center coordinates and rectangle colours.
#' 
#' @export
drawCell <- function(lsGene, lsCell, res, thresholdType, thresholdBinary, threshold, colours) {
  if (length(thresholdType) == 0 | length(thresholdBinary) == 0 | length(threshold$med) == 0) return(NULL)
  if (is.null(lsGene) | is.null(lsCell) | is.null(res) 
      | is.null(thresholdType) | is.null(thresholdBinary) | is.null(threshold$med)) return(NULL)
  g <- length(lsGene)
  n <- length(lsCell)
  xleftv = NULL
  ybottomv = NULL
  xrightv = NULL
  ytopv = NULL
  colv = NULL
  ## Gene names are written at the center of the associated rectangles ##
  centersx = NULL
  centersy = NULL
  for (i in 1:n) {
    xleftv = c(xleftv, rep(i*(topside+10), g))
    xrightv = c(xrightv, rep((i+1)*topside, g))
    centersx = c(centersx, rep(((2*i+1)*topside+i*10)/2, g))
    for (j in 1:g) {
      ybottomv = c(ybottomv, j*leftside)
      ytopv = c(ytopv, (j+1)*leftside)
      centersy = c(centersy, (2*j+1)*leftside/2)
      if (thresholdType == "Binary") colv = c(colv, if (res[j,i]/threshold$total[i]*100 >= thresholdBinary) colours[11] else "white")
      else {
        col = "white"
        if (res[j,i] >= threshold$quarterBounds[i]) col = colours[11]
        if (res[j,i] >= threshold$halfBounds[i]) col = colours[10]
        if (res[j,i] >= threshold$tquarterBounds[i]) col = colours[9]
        colv = c(colv, col)
      }
    }
  }
  return(list(xleftv, ybottomv, xrightv, ytopv, colv, centersx, centersy))
}

#____________________________#
#      MAST computation      #
#____________________________#

#ref: https://bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html
#' Perform MAST computation
#'
#' Performs MAST computation on the input SCESet with the two selected cell groups and the selected feature list.
#'
#' @param organism input organism
#' @param lsGene selected feature list
#' @param lsCell1 first selected cell group
#' @param lsCell2 second selected cell group
#' @return a list containing the log-regularized counts, the associated MAST object, two arrays for fold-change, etc.,
#' the computed p-values and the resulting heatmap.
#' 
#' @importFrom MAST FromMatrix
#' @importFrom MAST zlm
#' @importFrom MAST lrTest
#' @importFrom stats relevel
#' @importFrom stats p.adjust
#' @importFrom BiocGenerics order
#' @importFrom shiny withProgress
#' @importFrom shiny incProgress
#' @importFrom SummarizedExperiment assay
#' 
#' @export
mastComputation <- function(organism, lsGene, lsCell1, lsCell2) {
  if (is.null(organism) | is.null(lsCell1) | is.null(lsCell2)) return(NULL)
  sceset <- selectGroupInSceset(organism, cells=c(lsCell1, lsCell2), genes=lsGene)
  withProgress(message = "Set: MAST computation", value = 0, {
    incProgress(1/4, detail = "Fetching data...")
    res <- getData(sceset)
    if (is.null(res)) return(NULL)
    res <- apply(as.integer(res), 2, function(x) {storage.mode(x) <- 'integer'; x})
    log_counts <- log(res+1)/log(2)
    fData <- data.frame(names=rownames(log_counts), row.names = rownames(log_counts))
    groups <- factor(c(rep("group #1", length(lsCell1)), rep("group #2", length(lsCell2))), 
                 levels = c("group #1", "group #2"))
    names(groups) <- c(lsCell1, lsCell2)
    cData = data.frame(condition=groups, row.names=names(groups))
    obj <- FromMatrix(as.matrix(log_counts), cData, fData)
    
    if (length(lsGene) < 5000) incProgress(2/4, detail = "Model expression as function of condition & number of detected genes...")
    else incProgress(2/4, detail = "Model expression as function of condition & number of detected genes... (may take ~1 min)")
    condition <- factor(colData(obj)$condition)
    condition <- relevel(condition, "group #1")
    colData(obj)$condition <- condition
    zlmCond <- zlm(~condition, obj)
    lrt <- lrTest(zlmCond, "condition")
    
    incProgress(3/4, detail = "Getting results... (may take ~1 min)")
    summaryCond <- summary(zlmCond, doLRT="conditiongroup #2") 
    summaryDt <- as.data.frame(summaryCond$datatable)
    pVals <- unlist(summaryDt[summaryDt$component == "H",4])
    names(pVals) <- unlist(summaryDt[summaryDt$component == "H",1])
    pVals <- as.data.frame(p.adjust(pVals, method = "fdr"))
    lsGene <- if (is.null(lsGene)) getFeatureList(sceset) else lsGene
    ordpv <- order(pVals)
    pVals <- data.frame(pVals[ordpv,], row.names = unlist(lsGene[ordpv], recursive = FALSE))
    colnames(pVals) <- "adjusted p-value"
    fcHurdle <- merge(summaryDt[summaryDt$contrast=="conditiongroup #2" & summaryDt$component=="H", c("primerid", "Pr(>Chisq)")],
                      summaryDt[summaryDt$contrast=="conditiongroup #2" & summaryDt$component=="logFC", 
                                c("primerid", "coef", "ci.hi", "ci.lo")], by='primerid')
    fcHurdle[,"fdr"] <- p.adjust(fcHurdle[,"Pr(>Chisq)"], method = "fdr")
    fcHurdle <- fcHurdle[order(fcHurdle$fdr), ]
    numberDEgenes <- min(10, length(lsGene))
    entrez_to_plot <- fcHurdle[1:numberDEgenes, "primerid"]
    mat_to_plot <- data.matrix(data.frame(assay(obj[entrez_to_plot, ]), row.names = entrez_to_plot))
    incProgress(1, detail = "Done!")
  })
  return(list(obj, pVals, mat_to_plot))
}

library(scater)

##############################################################################
## Please refer to README in order to add a new organism to the application ##
##############################################################################

# Step 1: create the gene expression matrix M with features as rows and conditions as columns, for organism oo
# Step 2: create the vectors associated with gene names genes, number of stage s, cell id id, 
# embryo e and replicate r (optional for replicates)
# Step 3: sceset <- createSceset(M, genes, s, id, e, r, unitcount=c("raw", "FPKM", ...))
# Step 4: updateUtils(oo, sceset)

  # path="~/Documents/project/"
  # organisms <- c("Ciona intestinalis", "Caenorhabditis elegans")
  # library(celegansdata)
  # load(paste(path, "ciona.rda", sep = ""))
  # scesets <- list(ciona, celegans_sceset)
  # names(scesets) <- organisms
  # groupedScesets <- lapply(scesets, groupSceset)
  # names(groupedScesets) <- organisms
  # houseKeepGenes <- c("^ERCC-*")
  # pcaPlots <- lapply(lapply(scesets, function(x) x[row.names(x)[!grepl(houseKeepGenes, row.names(x))], ]), 
  #                    function(x) c(data=list(prcomp(t(counts(x)))), 
  #                                  df=list(data.frame(label=data.frame(label=pData(x)[, "Cell.ID"]), 
  #                                                     colour=data.frame(colour=paste("Embryo", pData(x)$Embryo))))))
  # names(pcaPlots) <- organisms
  # nscesets=length(scesets)
  # preStages=lapply(scesets, function(sceset) unique(pData(sceset)[, "EmbryoStage"]))
  # stages <- lapply(1:nscesets, function(i) lapply(preStages[[i]], function(e) paste(e, "cell", sep = "-")))
  # names(stages) <- organisms
  # cells <- lapply(1:nscesets, function(i) lapply(preStages[[i]], 
  #                                                function(s) unique(subset(pData(scesets[[i]]), EmbryoStage == s)[, "Cell.ID"])))
  # names(cells) <- organisms
  # for (o in organisms) names(cells[[o]]) <- stages[[o]]
  # save(organisms, scesets, groupedScesets, pcaPlots, preStages, stages, cells, houseKeepGenes, file=paste(path, "utils.Rdata", sep = ""))

#_____________________#
#      Update file    #
#_____________________#

#' Update utils file
#'
#' Updates the utils.Rdata file --side effect.
#'
#' @param newOrganismName the name of the new organism to add, as it should be displayed to the user
#' @param newOrganismSceset the SCESet associated with the gene expression matrix of the new organism -see README
#' @param houseKgene a regular expression for identifying by id the house keeping genes in the dataset -should
#' be explicit and non ambiguous-
#' 
#' @importFrom scater pData
#' @importFrom scater counts
#' @importFrom stats prcomp
#' 
#' @export
updateUtils <- function(newOrganismName, newOrganismSceset, houseKgene=NULL, path="~/Documents/project/") {
    load(paste(path, "utils.Rdata", sep = ""))
    ## List of organisms whose associated SCESets are currently stored in the application ##
    organisms <- c(organisms, newOrganismName)
    ## Associated SCESets ##
    scesets <- append(scesets, newOrganismSceset)
    names(scesets) <- organisms
    ## Associated SCESets where gene expression values from same-cell replicates are summed ##
    groupedScesets <- append(groupedScesets, groupSceset(newOrganismSceset))
    names(groupedScesets) <- organisms
    ## PCA computation for each dataset ##
    s <- newOrganismSceset[row.names(newOrganismSceset)[!grepl(houseKeepGenes, row.names(newOrganismSceset))], ]
    pcaPlots <- append(pcaPlots, list(c(data=list(prcomp(t(counts(s)))), 
                                   df=list(data.frame(label=data.frame(label=pData(s)[, "Cell.ID"]), 
                                   colour=data.frame(colour=paste("Embryo", pData(s)$Embryo)))))))
    names(pcaPlots) <- organisms
    ## Housekeeping gene extensions ##
    houseKeepGenes <- c("^ERCC-*")
    if (!(is.null(houseKgene))) houseKeepGenes <- c(houseKeepGenes, houseKgene)
    nscesets=length(scesets)
    preStages <- append(preStages, unique(pData(newOrganismSceset)[, "EmbryoStage"]))
    ##Stages and cell names sorted by stage
    stages <- append(stages, lapply(preStages[[nscesets]], function(e) paste(e, "cell", sep = "-")))
    names(stages) <- organisms
    cells <- append(cells, lapply(preStages[[nscesets]], 
                                  function(s) unique(subset(pData(newOrganismSceset), EmbryoStage == s)[, "Cell.ID"])))
    names(cells) <- organisms
    for (o in organisms) names(cells[[o]]) <- stages[[o]]
    save(organisms, scesets, groupedScesets, preStages, stages, cells, houseKeepGenes, file="utils.Rdata")
    return(NULL)
}

#_______________________#
#      Computations     #
#_______________________#

#' Compute grouped SCESet
#'
#' Computes grouped SCESet, that is, where gene expression values of same-cell replicates are summed, that is, each condition has
#' its own Cell.ID label.
#'
#' @param sceset the SCESet associated with the gene expression matrix of the new organism -see README
#' @return associated grouped SCESet
#' 
#' @importFrom scater pData 
#' @importFrom scater varMetadata
#' @importFrom scater counts
#' @importFrom scater newSCESet
#' @importFrom methods new
#' 
#' @export
groupSceset <- function(sceset) {
  if (is.null(sceset)) return(NULL)
  genes <- row.names(sceset) 
  cells <- unique(pData(sceset)[, "Cell.ID"])
  tmp <- sapply(cells, function(cell) {
      s <- sceset[, pData(sceset)[, "Cell.ID"] == cell]
      if (length(colnames(s)) < 2) counts(s) else rowSums(counts(s))
    })
  phenodata <- as.data.frame(cbind(cells, pData(sceset)[sapply(cells, 
               function(cell) which(pData(sceset)[, "Cell.ID"] == cell)[1]), 2:dim(pData(sceset))[2]]))
  row.names(phenodata) <- make.unique(as.character(cells))
  colnames(phenodata) <- c("Cell", colnames(phenodata)[2:dim(pData(sceset))[2]])
  phenodata <- new("AnnotatedDataFrame", data = phenodata, varMetadata = varMetadata(sceset))
  return(newSCESet(countData = tmp, phenoData = phenodata))
}

#_______________________#
#      Create SCESet    #
#_______________________#

#' Create SCESet
#'
#' Creates well-formatted SCESet for the application.
#'
#' @param countdata gene expression matrix with genes in rows, conditions in columns
#' @param stagePerCondition integer vector of developmental stages for each column of the matrix
#' @param idPerCondition string vector of names of cell associated to each condition
#' @param embryoPerCondition string vector of embryo names from which each condition comes
#' @param replicatePerCondition integer vector of number of replicate associated to each condition
#' @param unitcount can be "raw", "RPKM", "FPKM", "TPM" or "CPM"
#' @return associated SCESet
#' 
#' @importFrom scater newSCESet
#' @importFrom methods new
#' 
#' @export
createSceset <- function(countdata, featureList, stagePerCondition, idPerCondition, 
                         embryoPerCondition, replicatePerCondition=NULL, unitcount="raw") {
  print("[1] Get file...")
  idPerCondition <- sapply(idPerCondition, function(e) gsub("\n", "_", gsub(" ", "_", e)))
  #cells <- lapply(1:dim(countdata)[2], function(i) paste(embryoPerCondition[i], "_", idPerCondition[i], "_", stagePerCondition[i], sep=""))
  cells <- colnames(countdata)
  if (!(is.null(replicatePerCondition))) cells <- lapply(cells, function(e) paste("r", replicatePerCondition[i], "_", e, sep=""))
  cc <- make.unique(sapply(cells, as.character))
  countdata <- data.frame(countdata, row.names=featureList)
  colnames(countdata) <- cc
  print("[2] Removing missing values...")
  countdata <- countdata[apply(countdata, 1, function(x) !all(is.na(x))), ]
  print("[3] Building countData...")
  print("[4] Building phenoData...")
  phenodata <- data.frame(Cell=cc, Cell.ID=idPerCondition, Embryo=embryoPerCondition, EmbryoStage=stagePerCondition, row.names=cc)
  print("[5] Done!")
  phenodata <- new("AnnotatedDataFrame", 
                   data = phenodata, 
                   varMetadata =  data.frame(labelDescription=colnames(phenodata)))
  if (unitcount == "raw") return(newSCESet(countData = countdata, phenoData = phenodata))
  if (unitcount == "TPM") return(newSCESet(tpmData = countdata, phenoData = phenodata))
  if (unitcount == "FPKM") return(newSCESet(fpkmData = countdata, phenoData = phenodata))
  if (unitcount == "RPKM") return(newSCESet(fpkmData = countdata, phenoData = phenodata))
  if (unitcount == "CPM") return(newSCESet(cpmData = countdata, phenoData = phenodata))
  print("ERROR: Wrong unit.")
  return(NULL)
}
################################################
#                                              #
#      SHINY APPLICATION: SERVER PART          #
#                                              #
################################################

##########
# SERVER #
##########

shinyServer(function(input, output, session) {
  
  load_data()
  session$onSessionEnded(stopApp)
  
  ###################
  # REACTIVE VALUES #
  ###################
  
  #__________________________________________#
  #      Reactive value initialization       #
  #__________________________________________#
  
  ## ORGANISM & FILE                                                      ##
  ## Useful values: they need to be updated when starting the application ##
                          # all features     # cell groups 1 & 2 # selected features
  values <- reactiveValues(allGenes = NULL, lsCell = NULL, lsCell1 = NULL, lsCell2 = NULL, lsGene = NULL,
                           # number of selected genes # indices of selected genes
                           noGene = 0, geneSelection = NULL, reinit = F, genedata = NULL, pcaBiplot = NULL)
  
  ## MAST RESULTS ##
  plottingMAST <- reactiveValues(obj = NULL, pvalue = NULL, heatmap = NULL)
  
  reinitMAST <- function() {plottingMAST$obj <- NULL; plottingMAST$pvalue <- NULL; plottingMAST$heatmap <- NULL}
  
  #__________________________________#
  #      Reactive value update       #
  #__________________________________#
  
  ## Organism & gene selection ##
  observe({
    if (is.null(input$organism)) return(NULL)
    values$geneSelection <- input$geneData_rows_selected
    ## Housekeeping for "max number of DLL reached" ##
    R.utils::gcDLLs()
  })
  
  getSelectedGenes <- function() {
     return(values$allGenes[values$allGenes %in% rownames(values$genedata)[input$geneData_rows_selected]])
    }
  
  ## Gene & cell selection ##
  ooo <- observe({
    withProgress(message = "Initialize selection", value=0, {
      if (is.null(input$organism)) return(NULL)
      incProgress(1/4, detail = "Checking values...")
      selectRows(dataTableProxy("geneData"), NULL)
      values$geneSelection <- NULL
      values$noGene <- 0
      reinitMAST()
      reinitPattern()
      reinitPlottingHeatmap()
      values$pcaBiplot <- NULL
      values$genedata <- NULL
      incProgress(2/4, detail = "Checking values... ...")
      values$lsCell <- NULL
      values$lsCell1 <- NULL
      values$lsCell2 <- NULL
      values$allGenes <- getFeatureList(getSceset(input$organism))
      if (is.null(input$stage1) | is.null(input$stage2)) return(NULL)
      incProgress(3/4, detail = "Retrieving same-stage cells...")
      incProgress(1, detail = "Done!")
    })
  })
  
  ## Progress message for gene selection ##
  oo <- reactive({
    withProgress(message = "Updating selected genes", value = 0, {
      incProgress(1/2, detail = "Checking values...")
      selectRows(dataTableProxy("geneData"), values$geneSelection)
      incProgress(1, detail = "Done!")
      values$noGene <- length(values$geneSelection)
    })
  })
  
  #__________________________________#
  #      Reactive value events       #
  #__________________________________#
  
  ## Gene & cell selection ##
  observeEvent(
    input$runSelect, {
      if (is.null(input$organism)) return(NULL)
      ## No cell selection pop-up event ##
      output$noCellSelectionPopup <- renderUI({
        if (is.null(input$cells1) | is.null(input$cells2)) return(NULL)
      })
      output$noGeneSelectionPopup <- renderUI({
         if (is.null(values$noGene) | values$noGene < 2) return(NULL)
      })
      values$lsCell1 <- clean(input$cells1)
      values$lsCell2 <- clean(input$cells2)
      values$lsCell <- c(values$lsCell1, values$lsCell2)
      if (is.null(values$lsCell1) | is.null(values$lsCell2) | is.null(values$geneSelection)) return(NULL)
      withProgress(message = "Computation of values", value = 0, {
        incProgress(1/3, detail = "Checking organism...")
        sceset <- getSceset(input$organism)
        incProgress(2/3, detail = "Checking features...")
        values$lsGene <- getSelectedGenes()
        incProgress(1, detail = "Done!")
      })
      
      ## MAST computation ##
      ansArray <- mastComputation(input$organism, values$lsGene, values$lsCell1, values$lsCell2)
      plottingMAST$obj <- ansArray[[1]]
      plottingMAST$pvalue <- ansArray[[2]]
      plottingMAST$heatmap <- ansArray[[3]]
        
      updateActionButton(session, "deselectAllGenes", label = "Back to gene table")
      return(NULL)
    })

  ###################
  # INPUT SELECTION #
  ###################
  
  #______________________________#
  #      Organism selection      #
  #______________________________#
  
  output$selectedOrganism <- renderUI(
    selectInput("organism", label = "Choose an organism:", choices = organisms)
  )
  
  #____________________________________________#
  #      Stage, same-stage cell selection      #
  #____________________________________________#
  
  ## Stage selection ##
  output$selectedStage <- renderUI({
    if (is.null(input$organism)) return(NULL)
    if (values$reinit) {values$reinit <- F; return(NULL)}
    selectInput("stage", label = "Developmental stage:", choices = stages[[input$organism]])
  })
  
  ## Same-stage cell selection ##
  output$stageCells <- renderUI({
    if (is.null(input$organism) | is.null(input$stage)) return(NULL)
    if (values$reinit) return(NULL)
    dragSetUI("cells", textval = tolist(cells[[input$organism]][[input$stage]]))
  })
  
  #________________________________#
  #      Cell group selection      #
  #________________________________#

  ## Cell group 1 selection ##
  output$selectedCells1 <- renderUI({
      if (is.null(input$organism)) return(NULL)
      if (values$reinit) return(NULL)
      dropUI("cells1", row_n = 1, col_n = n)
    })
  
  ## Cell group 2 selection ##
  output$selectedCells2 <- renderUI({
    if (is.null(input$organism)) return(NULL)
    if (values$reinit) return(NULL)
    dropUI("cells2", row_n = 1, col_n = n)
    })

  observeEvent(input$deselectAllCells, {values$reinit <- T})
  
  #__________________________#
  #      Gene selection      #
  #__________________________#
  
  output$selectedGenes <- renderText({
    oo()
    sprintf("No. of selected genes: %d", values$noGene)
  })
  
  output$selectionAllGenes <- renderUI({
    if (!(is.null(plottingMAST$obj))) return(NULL)
    actionButton("selectAllGenes", "Select all", icon=icon("check-square"), class="act")
    })
  
  ## Incorrect gene selection pop-up event ##
  observeEvent(
    input$selectGenes,
    {
      if (!(length(input$geneData_rows_selected) == 2)) {
        output$selectionfromPopup <- renderUI({
          showModal(modalDialog(h4("When using 'Select range' button, you should only select 
                                   two genes that define the range of genes you want to select."), title = "Warning", 
                                footer = modalButton("Go back to application"),
                                size = "m", easyClose = FALSE, fade = FALSE))
          selectRows(dataTableProxy("geneData"), NULL)
          values$geneSelection <- NULL
          values$noGene <- 0
          reinitMAST()
        })
      }
      else {
        m <- min(input$geneData_rows_selected[1], input$geneData_rows_selected[2])
        M <- max(input$geneData_rows_selected[1], input$geneData_rows_selected[2])
        values$geneSelection <- match(row.names(values$genedata)[m:M], values$allGenes)
        oo()
        reinitMAST()
      }
    }
  )
  observeEvent(
    input$deselectAllGenes,
    {
      values$genesNo <- 0
      values$lsGene <- NULL
      values$geneSelection <- NULL
      oo()
      reinitMAST()
      updateActionButton(session, "deselectAllGenes", label = "Deselect all")
    }
  )
  observeEvent(
    input$selectAllGenes,
    {
      genesNo <- length(values$allGenes)
      values$geneSelection <- 1:genesNo
      oo()
      reinitMAST()
    }
  )
  observeEvent(
    input$selectWithCutoff,
    {
      if (is.null(plottingMAST$pvalue) | is.null(input$pvalueCutoff)) return(NULL)
      pvalue <- isolate(plottingMAST$pvalue)
      geneSelection <- (1:dim(pvalue)[1])[apply(pvalue, 1, function(x) (x <= input$pvalueCutoff))]
      values$geneSelection <- geneSelection
      genesNo <- length(geneSelection)
      oo()
    }
  )
  
  #__________________________#
  #      Run computation     #
  #__________________________#
  
  output$runSelection <- renderUI({
    if (length(values$geneSelection) < 2 | is.null(clean(input$cells1)) | is.null(clean(input$cells2))) return(NULL)
    actionButton("runSelect", "Run computation on selection", icon=icon("cogs"), class="act")
  })
  
  #__________________________#
  #      Workflow progress   #
  #__________________________#
  
  output$nextbuttonDATA <- renderUI({
    if (is.null(clean(input$cells1)) | is.null(clean(input$cells2))) return(NULL)
    actionButton("continueAData", "Go to differential expression\nanalysis tab", icon("play-circle"), class = 'act')
  })
  
  observeEvent(
    input$continueAData,
    {
      updateTabItems(session, "sidebarMenu", "deAnalysis")
      updateTabsetPanel(session, "tabset1", selected = "Gene selection")
    }
  )
  
  output$nextbuttonDE <- renderUI({
    if (is.null(clean(input$cells1)) | is.null(clean(input$cells2)) | is.null(input$selectWithCutoff)) return(NULL)
    actionButton("continueADE", "See differential expression\nanalysis heatmap", icon("chevron-right"), class = 'act')
  })
  
  observeEvent(
    input$continueADE,
    {
      updateTabItems(session, "sidebarMenu", "deAnalysis")
      updateTabsetPanel(session, "tabset1", selected = "DE heatmap")
    }
  )
  
  output$nextbuttonGE <- renderUI({
    if (length(values$geneSelection) < 1 | is.null(clean(input$cells1)) | is.null(clean(input$cells2))) return(NULL)
    actionButton("continueAGE", "Go to gene expression\n pattern explorer", icon("play-circle"), class = 'act')
  })
  
  observeEvent(
    input$continueAGE,
    {
      updateTabItems(session, "sidebarMenu", "paExplorer")
    }
  )
  
  #################
  # OUTPUT UPDATE #
  #################
  
  #__________________________#
  #   Gene selection tab     #
  #__________________________#
  
  ## Data viewing ##
  output$geneData <- DT::renderDataTable({
    # Displaying values for all cells #
    if (is.null(input$organism) | is.null(cells[[input$organism]][[input$stage]])) return(NULL)
    if (is.null(plottingMAST$pvalue)) {
      withProgress(message = "Retrieving data...", value = 0, {
        incProgress(1/2, detail = "")
        values$genedata <- counts(selectGroupInSceset(input$organism, cells[[input$organism]][[input$stage]]))
        incProgress(1, detail = "Done!")
      })
    }
    # Shows results from MAST computation #
    else {
      pvalue <- plottingMAST$pvalue
      values$genedata <- as.data.frame(cbind(data.matrix(pvalue), 
                       counts(selectGroupInSceset(input$organism, cells=values$lsCell, values$lsGene))[row.names(pvalue), ]))
    }
    return(rn(values$genedata))
  },
  server = TRUE)   
  
  ## PCA visualization ##
  output$plotVisu <- renderPlot({
    if (is.null(input$organism)) return(NULL)
    return(customPCAPlot(input$organism, "PCA plot for all genes"))
  })
  
  ## Count heatmap ##
  output$plotHeatmap <- renderPlot({
    qq()
    if (is.null(plottingHeatmap$res)) return(NULL)
    if (dim(plottingHeatmap$resSceset)[2] < 2) return(NULL)
    condition <- as.data.frame(pData(plottingHeatmap$resSceset)[,"EmbryoStage"])
    colnames(condition) <- "Stage"
    if (length(unique(condition)) < 2) aheatmap(plottingHeatmap$res, 
                                                main="Most highly expressed genes (raw count)",
                                                labRow=plottingHeatmap$genes, distfun='spearman') 
    else aheatmap(plottingHeatmap$res, labRow=plottingHeatmap$genes, annCol=as.data.frame(condition), distfun='spearman')
  })
  
  output$messagegeheatmap <- renderUI({
    if (is.null(plottingHeatmap$res)) 
      return(p(class = "text-muted", "Heatmap could not be computed."))
    if (dim(plottingHeatmap$resSceset)[2] < 2) 
      return(p(class = "text-muted", "There is only one cell in the selected stage."))
    return(NULL)
  })
  
  #_________________________________________#
  #   Differential expression analysis tab  #
  #_________________________________________#
  
  ## p-value cutoff selection button ##
  output$pvalueSelection <- renderUI({
    if (is.null(plottingMAST$pvalue)) return(NULL)
    actionButton("selectWithCutoff", "Select genes having p-value < cutoff", icon("check-square"), class="act")
  })
  
  ## p-value cutoff selection slider ##
  output$pvalueCutoffSelection <- renderUI({
    if (is.null(plottingMAST$pvalue)) return(NULL)
    sliderInput("pvalueCutoff", "Choose the p-value cutoff", min = 0, max = 1, value = 0.01, step = 0.01)
  })
  
  ## MAST results ##
  output$mastData <- renderPlot({
    if (is.null(plottingMAST$heatmap)) return(NULL)
    if (sum(colSums(plottingMAST$heatmap)) == 0) return(NULL)
    cond <- colData(plottingMAST$obj)$condition
    names(cond) <- values$lsCell
    cond <- cond[colnames(plottingMAST$obj)]
    ann <- as.data.frame(cond, row.names=names(cond))
    colnames(ann) <- "Condition"
    aheatmap(plottingMAST$heatmap, 
             annCol=ann,
             main="Most DE genes (log of raw count)",
             col=rev(colorRampPalette(colors = brewer.pal(name="PiYG",n=10))(20)))
  })
  
  output$messageheatmap <- renderUI({
    if (is.null(plottingMAST$heatmap)) 
      return(p(class = "text-muted", "Heatmap not available. Have you selected cell groups and genes? Have you run MAST computation?"))
    if (sum(colSums(plottingMAST$heatmap)) == 0) 
      return(p(class = "text-muted", "Heatmap computation error."))
    return(NULL)
  })
  
  ## PCA biplot for selected genes ##
  output$biplotUI <- renderUI({
    if (length(values$geneSelection) < 1) return(NULL)
    actionButton("biplotTrigger", "Update PCA biplot", icon=icon("cogs"), class="act")
  })
  
  observeEvent(
    input$biplotTrigger,
    {
      if (is.null(input$organism) | is.null(values$geneSelection) | is.null(values$allGenes)) return(NULL)
      genes <- getSelectedGenes()
      main <- sprintf("PCA plot for genes %s ...", paste0(genes[min(3, length(genes))], sep=","))
      values$pcaBiplot <- customPCABiPlot(input$organism, main=main, genes=genes)
    }
  )
  
  output$plotDEbiplot <- renderPlot({
    if (!is.null(values$pcaBiplot)) return(values$pcaBiplot)
    if (is.null(plottingHeatmap$genes)) return(NULL)
    genes <- plottingHeatmap$genes
    return(customPCABiPlot(input$organism, 
                           main="PCA plot for 10 most highly-expressed genes", 
                           genes=genes))
  })
  
  #______________________________________________________#
  #   Average gene expression pattern visualization tab  #
  #______________________________________________________#
  
  ## Threshold values ##
  threshold <- reactiveValues(total = NULL, quarterBounds = NULL, halfBounds = NULL, tquarterBounds = NULL, med = NULL)
  thresh <- reactive({
    if (length(values$geneSelection) > 0 & !is.null(clean(input$cells1)) & !is.null(clean(input$cells2)) & !is.null(input$organism)) {
    withProgress(message="Computing gene expression value bounds", value=0, {
      incProgress(1/6, detail="Getting values...")
      res <- rn(counts(selectGroupInSceset(input$organism, cells=c(clean(input$cells1), 
                                                                   clean(input$cells2)), getSelectedGenes())))
      if (is.null(res)) return(NULL)
      res <- data.matrix(res)
      incProgress(2/6, detail="Computing bounds 1/5...")
      threshold$total <- colSums(res, na.rm = T)
      incProgress(3/6, detail="Computing bounds 2/5...")
      threshold$quarterBounds <- sapply(threshold$total, function(x) x/4)
      incProgress(4/6, detail="Computing bounds 3/5...")
      threshold$halfBounds <- sapply(threshold$total, function(x) x/2)
      incProgress(5/6, detail="Computing bounds 4/5...")
      threshold$tquarterBounds <- sapply(threshold$total, function(x) 3*x/4)
      incProgress(1, detail="Computing bounds 5/5...")
      threshold$med <- mean(colMeans(res, na.rm = T))/max(threshold$total)*100
    })
    }
    else return(NULL)
  })
  
  ## Pattern values ##
  plottingPatterns <- reactiveValues(plot1 = NULL, rect1 = NULL, text1 = NULL, axis1 = NULL,
                             plot2 = NULL, rect2 = NULL, text2 = NULL, axis2 = NULL)
  
  reinitPattern <- function() {
    plottingPatterns$plot1 <- NULL
    plottingPatterns$rect1 <- NULL
    plottingPatterns$text1 <- NULL
    plottingPatterns$axis1 <- NULL
    plottingPatterns$plot2 <- NULL
    plottingPatterns$rect2 <- NULL
    plottingPatterns$text2 <- NULL
    plottingPatterns$axis2 <- NULL
  }
  
  output$messagegepattern <- renderUI({
    if (length(values$geneSelection) < 1 | is.null(clean(input$cells1)) | is.null(clean(input$cells2))) 
      return(p(class = "text-muted", "Patterns not available. Have you selected cell groups and genes?"))
    return(NULL)
  })
  
  patt1 <- reactive({
    if (length(values$geneSelection) > 0 & !is.null(clean(input$cells1)) & !is.null(clean(input$cells2)) & !is.null(input$organism)) {
      # fct <- switch(input$organism, 
      #        "Ciona intestinalis"=patternPlotCiona,
      #        "Caenorhabditis elegans"=patternPlotCelegans,
      #        patternPlot)
      
      fct <- patternPlot
      
      ls <- do.call(fct, list(getSelectedGenes(), input$cells1, 
                              input$organism, values$allGenes, input$selectThresholdType, input$threshold, threshold), 
                    quote = FALSE, envir = parent.frame())
      plottingPatterns$plot1 <- ls[[1]]
      plottingPatterns$rect1 <- ls[[2]]
      plottingPatterns$text1 <- ls[[3]]
      plottingPatterns$axis1 <- ls[[4]]
    }
    return(NULL)
  })
  
  patt2 <- reactive({
    if (length(values$geneSelection) > 0 & !is.null(clean(input$cells1)) & !is.null(clean(input$cells2)) & !is.null(input$organism)) {
      # fct <- switch(input$organism, 
      #               "Ciona intestinalis"=patternPlotCiona,
      #               "Caenorhabditis elegans"=patternPlotCelegans,
      #               patternPlot)
      
      fct <- patternPlot
      
      ls <- do.call(fct, list(getSelectedGenes(), input$cells2, 
                              input$organism, values$allGenes, input$selectThresholdType, input$threshold, threshold, nb=2), 
                    quote = FALSE, envir = parent.frame())
      plottingPatterns$plot2 <- ls[[1]]
      plottingPatterns$rect2 <- ls[[2]]
      plottingPatterns$text2 <- ls[[3]]
      plottingPatterns$axis2 <- ls[[4]]
    }
    return(NULL)
  })
  
  ## Threshold selection ##
  output$thresholdIn <- renderUI({
    thresh()
    if (length(values$geneSelection) > 0 & !is.null(clean(input$cells1)) & !is.null(clean(input$cells2)) & !is.null(input$organism)) {
    if (input$selectThresholdType == "Binary") {
      ## If gene expression is above the threshold, then the gene is considered activated/ON ##
      sliderInput("threshold", "Choose the gene expression binary threshold (in %)",
                  min = 0, max = 100, value = threshold$med, step = 5)
    }
      else return(NULL)
    }
    else return(NULL)
  })
  
  ## Gene expression pattern visualization in cell group 1 ##
  output$plotGenePattern <- renderPlot({
    patt1()
    if (is.null(plottingPatterns$plot1)) return(NULL)
    plottingPatterns$plot1
    plottingPatterns$rect1
    plottingPatterns$text1
    plottingPatterns$axis1
  })
  
  ## Gene expression pattern visualization in cell group #2 ##
  output$plotGenePattern2 <- renderPlot({
    patt2()
    if (is.null(plottingPatterns$plot2)) return(NULL)
    plottingPatterns$plot2
    plottingPatterns$rect2
    plottingPatterns$text2
    plottingPatterns$axis2
  })
  
  ################
  # COMPUTATIONS #
  ################
  
  #_______________________________________________________#
  #   Count heatmap with the most highly expressed genes  #
  #_______________________________________________________#
  
  ## Reactive values for plot ##
  plottingHeatmap <- reactiveValues(res = NULL, genes = NULL, resSceset = NULL)
  
  ## Initialize values for plot ##
  reinitPlottingHeatmap <- function() {
    plottingHeatmap$res <- NULL 
    plottingHeatmap$resSceset <- NULL 
    plottingHeatmap$genes <- NULL 
  }
  
  ## Compute values for plot ##
  qq <- reactive({
    if (is.null(input$organism) | is.null(input$stage)) return(NULL)
    resSceset <- rn(selectGroupInSceset(input$organism, cells[[input$organism]][[input$stage]]))
    if (is.null(resSceset)) return(NULL)
    res <- rn(counts(resSceset))
    row.names(res) <- values$allGenes
    # Number of highly-expressed genes selected #
    number = min(10, dim(res)[1])
    # Selecting the number most expressed genes in data #
    withProgress(message = "Selecting genes for count heatmap...", value = 0, {
      incProgress(1/3, detail = "Get indices...")
      indices <- sort(unlist(lapply(1:dim(res)[1], function(i)sum(res[i,]))), decreasing = TRUE, index.return = TRUE)$ix[1:number]
      incProgress(2/3, detail = "Get genes...")
      genes <- row.names(res)[indices]
      incProgress(1, detail = "Done!")
      res <- rn(res[indices, ])
      if (is.null(rn(dim(res)))) return(NULL)
      res <- rn(res[apply(res, 1, function(x) !(sum(x) == 0)),])
    })
    plottingHeatmap$res <- res
    plottingHeatmap$resSceset <- rn(resSceset[indices,])
    plottingHeatmap$genes <- genes
  })
})
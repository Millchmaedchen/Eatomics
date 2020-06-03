expDesignModule_UI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      uiOutput(ns("diff.gs.collection")),
      uiOutput(ns("conditional_grouping")),
       uiOutput(ns("conditional_subselect")),
       checkboxInput(ns("ContinChoice_gsea"), "Use continuous response instead of grouping", FALSE),
      #checkboxInput("imputeFor_gsea", "Impute missing values", FALSE),
      #checkboxInput("remove_sv", "Remove surrogate variables", FALSE),
      checkboxInput(ns("includeCovariates_gsea"), "Include parameters as covariates", FALSE),
      conditionalPanel("input.includeCovariates_gsea == TRUE", ns = ns, 
                       uiOutput(ns("covariatesChoice_gsea"))),
      checkboxInput(ns("expandFilter_gsea"), "Stratification and filter",  FALSE),
      conditionalPanel("input.expandFilter_gsea == TRUE", ns = ns,
                       uiOutput(ns("filter_group_gsea"))),
      conditionalPanel("input.expandFilter == TRUE", ns = ns,
                       uiOutput(ns("filter_level_gsea"))),
      conditionalPanel("input.expandFilter_gsea == TRUE", ns = ns,
                       uiOutput(ns("selectContrast_gsea"))),
      actionButton(ns("analyze_diff_gsea"),"Analyze",class = "btn-primary")
     ), 
    mainPanel(
      fluidRow(
        column(8,
               plotlyOutput(ns("limma"),
                            height = 500
               ) %>% withSpinner()
        ),
        
        absolutePanel(
          top=5,
          right = 20,
          draggable = TRUE,
          wellPanel(
            #sidebarPanel(position = "right",
            numericInput(ns("adj.P.Val"), "Adjusted P value threshold", 0.05, min = 0, max = 1, step= 0.01),
            sliderInput(ns("logFC"), "Log Fold Change", 0, min = 0, max = 10, step = 0.1)
          )
        ),
        br(),
        column(12,
               DT::dataTableOutput(ns('up')),
               DT::dataTableOutput(ns('down')),
               br(),
               uiOutput(ns("labelColBox")),
               checkboxInput(ns("showLabels"), "Blend in PatientID", FALSE),
               plotOutput(ns("boxPlotUp")),
               plotOutput(ns("boxPlotDown")),
               htmlOutput(ns("doc1"))
        )),
      
      br(),
      #downloadButton("report", "Generate report"),
      #downloadButton("reportDataDL", "Download report data"),
      br()
    )
  )
}

# Function for module server logic
expDesignModule <- function(input, output, session, ssgsea_data_update = NULL, sessionID = NULL, gs_file_list = NULL, ClinData = NULL) {
  ns <- session$ns
  
  ClinDomit <- reactiveValues()
  proteinAbundance <- reactiveValues()
  limmaResult <- reactiveValues(gene_list = FALSE)
  reportBlocks <- reactiveValues()
  # gs_file_list_f = reactive({
  #   browser()
  #   if (is.null(gs_file_list)){NULL}
  #   else {return(gs_file_list)}
  #     })
  
  observe({
    ssgsea_data_update
    
    validate(need(!is.null(gs_file_list()), "Please perform ssGSEA on the previous tab first."))
    output$diff.gs.collection <- renderUI({
      selectInput(
        inputId = ns("diff.gs.collection_file"),
        label = strong("Choose enrichment score file"),
        choices = gs_file_list(),
        multiple = FALSE,
        selectize = TRUE
      )
    })
  })
  
  observe({input$diff.gs.collection_file
    req(input$diff.gs.collection_file)
    proteinAbundance$original = read_table2(paste(sessionID, "/", input$diff.gs.collection_file, sep = ""), 
                                            col_types = cols(Description = col_skip(), 
                                                             NA. = col_logical()), skip = 2) %>% 
      mutate(`Gene names` = Name) 
  })
  q
  ClinColClasses <- reactive({
    req(ClinData())
    ClinDomit$data = ClinData() %>% janitor::clean_names()
    df = ClinData()
    #df = ClinDomit$data
    df = lapply(df, class)
  })
  
  output$conditional_grouping <- renderUI({
    selectInput(
      inputId = ns("GR_fatcor_gsea"), 
      label = strong("Select the clinical grouping factor"),
      choices = as.list(colnames(ClinData())),
      multiple = FALSE,
      selectize = TRUE
    )
  })
  

  output$conditional_subselect <- renderUI({
    req(ClinData, 
        #ClinColClasses(), 
        input$GR_fatcor_gsea)
    if (ClinColClasses()[input$GR_fatcor_gsea]=='factor' | ClinColClasses()[input$GR_fatcor_gsea]=='logical' ){
      selectizeInput(inputId = ns("levels"),
                     label= "Select two groups to compare",
                     choices = ClinData() %>% pull(input$GR_fatcor_gsea) %>% levels(),
                     multiple = TRUE, 
                     options = list(maxItems = 2)
      )
    } else {
      d = ClinData() %>% pull(input$GR_fatcor_gsea)
      sliderInput(
        inputId = ns("num.cutoff"),
        label = "Select cutoff to divde numeric value:",
        min = min(d, na.rm = TRUE),
        max = max(d, na.rm = TRUE),
        value = colMeans(ClinData()[input$GR_fatcor_gsea], na.rm = TRUE), round = T
      )
    }
  }) 
  
  observeEvent(input$includeCovariates_gsea, {
    output$covariatesChoice_gsea<- renderUI({
      selectInput(
        inputId = ns("covariates"),
        label = strong("Select factors to include as covariates."),
        choices = as.list(colnames(ClinData())),
        multiple = TRUE
      )
    })
  })
  
  observeEvent(input$expandFilter_gsea, {
    output$filter_group_gsea <- renderUI({
      selectInput(
        inputId = ns("filter_GR_fatcor"),
        label = strong("Select a second parameter"),
        selected = 3,
        choices = as.list(colnames(ClinData())),
        multiple = FALSE,
        selectize = TRUE
      )
    })
    
    output$filter_level_gsea <- renderUI({
      req(input$filter_GR_fatcor)
      if (ClinColClasses()[input$filter_GR_fatcor]=='factor' | ClinColClasses()[input$filter_GR_fatcor]=='logical'){
        selectizeInput(inputId = ns("filter_levels"),
                       label = "Filter: Select groups to include in the analysis",
                       choices = ClinData() %>% pull(input$filter_GR_fatcor) %>% levels(),
                       multiple = TRUE
        )
      } else {
        d = ClinDomit$data %>% pull(input$filter_GR_fatcor)
        sliderInput(
          inputId = ns("filter_num.cutoff"),
          label = "Select cutoff to divide numeric value:",
          min = min(d, na.rm = TRUE),
          max = max(d, na.rm = TRUE),
          value = colMeans(ClinData()[input$filter_GR_fatcor], na.rm = TRUE), round = T
        )
      }
    })
    output$selectContrast_gsea <- renderUI({
      req(ClinDomit$mainParameter)
      selectizeInput(inputId = ns("contrastLevels"), 
                     label = "Stratify: Select the two groups you want to calculate the difference on.",
                     #choices = ClinData() %>% pull(mainParameter) %>% levels()
                     choices = ClinDomit$data %>% pull(ClinDomit$mainParameter) %>% levels(),
                     multiple = TRUE, 
                     options = list(maxItems = 2)
      )
    })
  }, ignoreNULL = FALSE, ignoreInit = TRUE)
  
  observe({
    ClinColClasses()[input$GR_fatcor_gsea] != "numeric"
    updateCheckboxInput(session, "ContinChoice_gsea", value = FALSE)
  })
  observe({
    req(input$GR_fatcor_gsea)
    input$GR_fatcor_gsea
    updateCheckboxInput(session, "expandFilter_gsea", value = FALSE)
  })
  
  ClinColClasses_2 <- reactive({
    df = ClinDomit$data
    df = lapply(df, class)
  })
  
  
  
  # set and manipulate chosen parameter to being cat from cont and form new groups from stratified setups
  observeEvent(c(input$filter_GR_fatcor, 
                 input$GR_fatcor_gsea,
                 input$ContinChoice_gsea, 
                 input$num.cutoff 
  )
  , {
    
    ClinDomit$mainParameter = janitor::make_clean_names(input$GR_fatcor_gsea)
    
    ## categorize numeric data - first parameter
    if (input$ContinChoice_gsea == FALSE & ClinColClasses_2()[ClinDomit$mainParameter]=='numeric') {
      req(input$num.cutoff)
      ClinDomit$data = ClinDomit$data %>% 
        mutate(categorizedParameter = 
                 cut(dplyr::pull(ClinDomit$data, ClinDomit$mainParameter), 
                     breaks = c(-Inf, input$num.cutoff, Inf), 
                     labels = c(paste('less than or equal to', input$num.cutoff, sep='_'), paste('greater than', input$num.cutoff,sep='_'))
                 ))  %>% 
        mutate_if(is.character, as.factor)
      #%>% 
      #dplyr::select(-c(!!mainParameter)) 
      colnames(ClinDomit$data)[colnames(ClinDomit$data) == "categorizedParameter"] = paste(input$GR_fatcor_gsea, "cat", sep = "_", collapse = "_") %>% janitor::make_clean_names() 
      ClinDomit$data = ClinDomit$data[,!duplicated(colnames(ClinDomit$data), fromLast = TRUE)]
      ClinDomit$mainParameter = paste(input$GR_fatcor_gsea,  "cat", sep = "_", collapse = "_") %>% janitor::make_clean_names() 
    }
    if (is.null(need(input$expandFilter_gsea, FALSE)) & is.null(need(input$filter_GR_fatcor, FALSE))) {
      
      ## categorize second numeric parameter
      ClinDomit$filterParameter =  janitor::make_clean_names(input$filter_GR_fatcor)
      if (ClinColClasses_2()[ClinDomit$filterParameter]=='numeric') {
        req(input$filter_num.cutoff)
        ClinDomit$data = ClinDomit$data %>% 
          mutate(filterParameter = 
                   cut(dplyr::pull(ClinDomit$data, ClinDomit$filterParameter), 
                       breaks = c(-Inf, input$filter_num.cutoff, Inf), 
                       labels = c(paste('less than or equal to', input$filter_num.cutoff, sep='_'), paste('greater than', input$filter_num.cutoff,sep='_'))
                   )) %>% 
          mutate_if(is.character, as.factor)
        # if first parameter stays continuous, the second parameter becomes a filter and needs to be saved for experimental design creation 
        if (ClinColClasses_2()[ClinDomit$mainParameter]=='numeric'){
          ## rename and save filter parameter
          colnames(ClinDomit$data)[colnames(ClinDomit$data) == "filterParameter"] = paste(ClinDomit$filterParameter) 
          # ClinDomit$filterParameter = input$filter_GR_fatcor
        } else {## unite cat first parameter and categorized second parameter, when first is cat or categorized
          ClinDomit$data = ClinDomit$data %>% 
            unite("newFactor", ClinDomit$mainParameter, filterParameter, remove = FALSE) %>% 
            mutate_if(is.character, as.factor)
          colnames(ClinDomit$data)[colnames(ClinDomit$data) == "newFactor"] = paste(ClinDomit$mainParameter, ClinDomit$filterParameter, sep = "_", collapse = "_")
          ClinDomit$data = ClinDomit$data[,!duplicated(colnames(ClinDomit$data), fromLast = TRUE)]
          ClinDomit$mainParameter = paste(ClinDomit$mainParameter, ClinDomit$filterParameter, sep = "_", collapse = "_")
        }
      } else{
        if (ClinColClasses_2()[ClinDomit$mainParameter]=='numeric'){
          ## rename and save filter parameter
          colnames(ClinDomit$data)[colnames(ClinDomit$data) == "filterParameter"] = paste(ClinDomit$filterParameter) 
          # ClinDomit$filterParameter = input$filter_GR_fatcor
        }else{
          #req(input$filter_levels)
          ## unite cat first and cat second parameter in the case of two cat parameters in the first place
          ClinDomit$data = ClinDomit$data %>% unite("newFactor", ClinDomit$mainParameter, ClinDomit$filterParameter, remove = FALSE)  %>% mutate_if(is.character, as.factor)
          colnames(ClinDomit$data)[colnames(ClinDomit$data) == "newFactor"] = paste(ClinDomit$mainParameter, ClinDomit$filterParameter, sep = "_", collapse = "_")
          ClinDomit$data = ClinDomit$data[,!duplicated(colnames(ClinDomit$data), fromLast = TRUE)]
          ClinDomit$mainParameter = paste(ClinDomit$mainParameter, ClinDomit$filterParameter, sep = "_", collapse = "_")
          
        }
        
      }
    }
  }, ignoreInit = TRUE
  #    , once = TRUE
  )
  
  observeEvent(input$analyze_diff_gsea ,{
    validate(need(proteinAbundance$original , "Please choose a file from the enrichment calculation first."))
    
    mainParameter = janitor::make_clean_names(ClinDomit$mainParameter)
    if (!is.null(ClinDomit$filterParameter)) {
      filterParameter = janitor::make_clean_names(ClinDomit$filterParameter)
    }
    #    }
    if (input$expandFilter_gsea == TRUE & input$ContinChoice_gsea == FALSE) {
      req(input$contrastLevels)
    }
    if (is.null(input$contrastLevels)){
      contrastLevels = input$levels
    } else{
      contrastLevels = input$contrastLevels
    }
    
    covariates = input$covariates#, TODO: Set up instead of filter for example
    covariates = janitor::make_clean_names(covariates)
    
    if (!is.null(ClinDomit$filterParameter) & ClinColClasses_2()[mainParameter] == "numeric") {
      ClinData = ClinDomit$data %>% 
        janitor::clean_names()  %>% 
        dplyr::select(patient_id, mainParameter, covariates, all_of(filterParameter)) %>% 
        remove_missing() %>% 
        dplyr::filter(!!sym(filterParameter) %in% !!input$filter_levels) %>% 
        select(-filterParameter)
    } else {
      ClinData = ClinDomit$data %>% 
        janitor::clean_names() %>% 
        dplyr::select(patient_id, mainParameter, covariates) %>% 
        remove_missing() 
    }
    ClinColClasses = lapply(ClinData, class)
    
    if (input$ContinChoice_gsea == FALSE & ClinColClasses[ClinDomit$mainParameter]=='numeric')
    {req(input$num.cutoff)}
    
    # Prep experimental design for first parameter being cat
    if (ClinColClasses[mainParameter]=='factor' | ClinColClasses[mainParameter]=='logical' ){
      #Reorder levels to mirror users selection of contrasts:
      if (ClinColClasses[mainParameter]=='factor') {
        ClinData[,mainParameter] = fct_relevel(dplyr::pull(ClinData, mainParameter), contrastLevels)
      }
      if (length(covariates) == 0){
        expDesign = model_matrix(ClinData, as.formula(paste("~0", mainParameter, sep = "+", collapse = "+")))       
      } else {
        expDesign = model_matrix(ClinData, as.formula(paste("~0", mainParameter, paste(covariates, sep = "+", collapse = "+"), sep = "+", collapse = "+")))       
      }
    }
    # Prep experimental design for first parameter being cont.     
    if (input$ContinChoice_gsea == TRUE){
      if (length(covariates) == 0){
        expDesign = model_matrix(ClinData, as.formula(paste("~ 1", mainParameter, sep = "+", collapse = "+")))       
      } else {
        expDesign = model_matrix(ClinData, as.formula(paste("~ 1", mainParameter, paste(covariates, sep = "+", collapse = "+"), sep = "+", collapse = "+")))       
      }
    }
    
    expDesign = data.frame(expDesign, check.names = FALSE)
    rownames(expDesign) <- ClinData$patient_id
    expDesign = matchedExpDesign(expDesign, proteinAbundance$original)
    ClinDomit$designMatrix = expDesign
    
    
    #observeEvent(input$analyzeLimma ,{
    req(ClinDomit$designMatrix)
    validate(need(proteinAbundance$original , "Please upload a proteinGroups file first (previous tab)."))
    validate(need(sum(ClinDomit$designMatrix[,1])>=3, "The experimental design does not contain three or more samples to test on."))
    
    expDesignInst = ClinDomit$designMatrix %>% janitor::clean_names()  
    if (input$ContinChoice_gsea == FALSE){
      validate(need(sum(ClinDomit$designMatrix[,2])>=3, "The experimental design does not contain three or more samples to test on."))
      validProteins =  proteinAbundance$original[, rownames(expDesignInst)]
      validProteins_1 = validProteins[,expDesignInst[,1]]
      validProteins_2 = validProteins[,expDesignInst[,2]]
      validProteinsLog = (rowSums(is.na(validProteins_1)) < length(validProteins_1)*0.5) + (rowSums(is.na(validProteins_2)) < length(validProteins_2)*0.5)
      
      ## Prepare protein abundance for limma
      #match_Protdata =  match_Protdata[, rownames(designMatrix)]
    #  if (exists(input$imputeForLimma) & input$imputeForLimma == TRUE){
    #    proteinAbundanceLimma <- proteinAbundance$imputed[as.vector(validProteinsLog > 0), c("Gene names", rownames(expDesignInst))] %>% as_tibble()
    #  } else {
        proteinAbundanceLimma <- proteinAbundance$original[as.vector(validProteinsLog > 0),c("Gene names", rownames(expDesignInst))] %>% as_tibble()
    #  }
      #proteinAbundanceLimma <- proteinAbundance$imputed[as.vector(validProteinsLog > 0),]
      proteinAbundanceLimma <- proteinAbundanceLimma %>% column_to_rownames("Gene names")
      proteinAbundanceLimma <- proteinAbundanceLimma[which(apply(proteinAbundanceLimma, 1, var, na.rm = TRUE) != 0), ]
      
      # no intercept as no reference level is assumed in this setup
      fit <- lmFit(proteinAbundanceLimma, as.matrix(expDesignInst))
      sv_contrasts <- makeContrasts(contrasts=paste(colnames(expDesignInst)[1],colnames(expDesignInst)[2], sep="-"), levels= expDesignInst)
      fit <- contrasts.fit(fit, contrasts=sv_contrasts)
    } else {
      validProteins =  proteinAbundance$original[, rownames(expDesignInst)]
      validProteinsLog = (rowSums(!is.na(validProteins)) >= 5)
      proteinAbundanceLimma <- proteinAbundance$imputed[as.vector(validProteinsLog > 0), c("Gene names", rownames(expDesignInst))]
      proteinAbundanceLimma <- proteinAbundanceLimma %>% column_to_rownames("Gene names")
      proteinAbundanceLimma <- proteinAbundanceLimma[which(apply(proteinAbundanceLimma, 1, var, na.rm = TRUE) != 0), ]
      
      fit <- lmFit(proteinAbundanceLimma, as.matrix(expDesignInst))
    }
    fit3 <- eBayes(fit, trend = TRUE) 
    if (input$ContinChoice_gsea == FALSE) {
      limmaResult$gene_list <- limma::topTable(fit3, coef=1, number = 1e+09, sort.by="P", adjust.method="BH")
    } else if (input$ContinChoice_gsea == TRUE)  {
      limmaResult$gene_list <- limma::topTable(fit3, coef=2, number = 1e+09, sort.by="P", adjust.method="BH")
    }
  }, ignoreInit = TRUE)
  
  
  #limma_input <- reactive({
  limma_input <- eventReactive(c(
    #ClinDomit$designMatrix,
    input$analyzeLimma ,
    input$adj.P.Val,
    input$logFC
  ),{
    req(limmaResult$gene_list)
    gene_list<-limmaResult$gene_list
    message("Genes in Limma: ", nrow(gene_list))
    gene_list$threshold = as.factor(abs(gene_list$logFC) > input$logFC & gene_list$adj.P.Val < input$adj.P.Val)
    
    reportBlocks$volcano_plot = 
      ggplot(data=gene_list, aes(x=logFC, y=-log10(P.Value) , colour=threshold)) +
      #ggtitle(paste(Factor,": ", names(ClinDomit$designMatrix)[2]," vs ", names(ClinDomit$designMatrix[1]), Filnames, collapse = "")) +
      geom_point(shape=20, data = gene_list, aes(text= rownames(gene_list)), size = 1, alpha = 0.4) +
      #labs(color = paste("Threshold \n adj.p < ", input$adj.P.Val, " and \n log2FC +/- ", input$logFC)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
      scale_color_tableau() +
      theme_light() +
      theme(legend.title = element_blank())
    
  })
  
  
  output$limma <- renderPlotly({
    req(limmaResult$gene_list)
    input$adj.P.Val
    input$logFC 
    input$analyzeLimma 
    gene_list <-limmaResult$gene_list
    
    ##select upregulated dataframe
    
    limmaResult$upregulated = as.factor(gene_list$logFC >= input$logFC & gene_list$adj.P.Val <= input$adj.P.Val)
    limmaResult$up<-na.omit(gene_list[limmaResult$upregulated=="TRUE",])
    limmaResult$df_up <- limmaResult$up %>% tibble::rownames_to_column("Gene.name")
    
    ##select down regulated dataframe
    
    limmaResult$downregulated= as.factor(gene_list$logFC <= (-input$logFC) & gene_list$adj.P.Val <= input$adj.P.Val)
    limmaResult$down<-na.omit(gene_list[ limmaResult$downregulated =="TRUE",])
    limmaResult$df_down<-limmaResult$down %>% tibble::rownames_to_column("Gene.name")
    
    
    reportBlocks$volcano_plot<-limma_input()
    
    #Prepare interactive plot, reactive title and legend
    
    if (!is.null(input$expandFilter_gsea) && input$expandFilter_gsea == TRUE) {
      Filnames = paste(" only ", input$filter_levels[1], collapse= "")
    } else {
      Filnames = ""
    }
    
    input$analyzeLimma 
    isolate(
      if (input$ContinChoice_gsea == FALSE){
        title_begin = paste(names(ClinDomit$designMatrix)[1], 
                            " regulated when compared to ", 
                            names(ClinDomit$designMatrix[2]), 
                            Filnames, collapse = "")
      } else {
        title_begin =  paste("Proteins regulated with regard to ", names(ClinDomit$designMatrix)[2], Filnames, collapse = "")
      })
    
    pp <- ggplotly(reportBlocks$volcano_plot, tooltip = "text") %>% 
      plotly::layout(
        legend = list(title = list(text = paste("Threshold: \n adj.p < ", 
                                                input$adj.P.Val, 
                                                " \n and \n log2FC +/- ", 
                                                input$logFC)
        )), 
        title = list(text = paste0('Volcano plot',
                                   '<br>',
                                   '<sup>',
                                   title_begin,
                                   collapse = "",
                                   '</sup>'), 
                     yref = "container", 
                     yanchor = "bottom",
                     pad = list(b = 20)
        ),
        margin = list(t = 75, b = 20)) 
    
    pp
  })
  
  output$boxPlotUp <- renderPlot({
    validate(need(input$up_rows_selected, message = "Select up- or downregulated proteins from the table to generate a display of protein abundance in the selected groups."))
    input$analyzeLimma 
    input$up_rows_selected
    input$labelColBox
    original = proteinAbundance$original %>% column_to_rownames("Gene names") %>% as.data.frame()
    isolate(if(ClinColClasses()[input$GR_fatcor_gsea]=='factor' | input$ContinChoice_gsea == FALSE) {
      experimentalDesign = ClinDomit$designMatrix[,1:2] %>% rownames_to_column(var = "PatientID") %>% gather(key= group, value = value, -PatientID ) %>% filter(value > 0) %>% dplyr::select(-value) #%>% left_join(., ClinData()[, c("PatientID", input$filter_GR_fatcor, input$labelColBox)])
      ClinData = ClinData()
      GR_fatcor = "group"
      reportBlocks$boxPlotUp = createBoxPlot( rows_selected = input$up_rows_selected, 
                                              proteinData = original, 
                                              filter_GR_fatcor = input$filter_GR_fatcor, 
                                              ClinData = ClinData(), 
                                              experimentalDesign = experimentalDesign, 
                                              GR_fatcor = GR_fatcor, 
                                              limmaResult = limmaResult$df_up, 
                                              labelColBox = input$labelColBox)
      return(reportBlocks$boxPlotUp)
    } else if(input$ContinChoice_gsea == TRUE){
      reportBlocks$scatterPlotUp = createLineScatterPlot(input$up_rows_selected, original, ClinData(), input$GR_fatcor_gsea, limmaResult$df_up, input$labelColBox)
      return(reportBlocks$scatterPlotUp)
    })
  })
  
  output$boxPlotDown <- renderPlot({
    validate(need(input$down_rows_selected, FALSE))
    original = proteinAbundance$original %>% column_to_rownames("Gene names") %>% as.data.frame()
    if(ClinColClasses()[input$GR_fatcor_gsea]=='factor' | input$ContinChoice_gsea == FALSE) {
      experimentalDesign =  ClinDomit$designMatrix[,1:2] %>% rownames_to_column(var = "PatientID") %>% gather(key= group, value = value, -PatientID ) %>% filter(value > 0) %>% dplyr::select(-value) #%>% left_join(., ClinData()[, c("PatientID", input$filter_GR_fatcor, input$labelColBox)])
      ClinData = ClinData()
      GR_fatcor = "group"
      reportBlocks$boxPlotDown = createBoxPlot( rows_selected = input$down_rows_selected, 
                                                proteinData = original, 
                                                filter_GR_fatcor = input$filter_GR_fatcor, 
                                                ClinData = ClinData, 
                                                experimentalDesign = experimentalDesign, 
                                                GR_fatcor = GR_fatcor, 
                                                limmaResult = limmaResult$df_down, 
                                                labelColBox = input$labelColBox)
      return(reportBlocks$boxPlotDown)
    } else if(input$ContinChoice_gsea == TRUE){
      reportBlocks$scatterPlotDown = createLineScatterPlot(input$down_rows_selected, original, ClinData(), input$GR_fatcor_gsea, limmaResult$df_down, input$labelColBox)
      return(reportBlocks$scatterPlotDown)
    }
  })
  
  createLineScatterPlot <- function(rows_selected, proteinData, CD_clean, GR_fatcor, limmaResult, labelColBox){
    s_up = rows_selected
    #limmaResult[s_up,]$Gene.name
    PatientID = colnames(proteinData)
    PD = as_tibble(as.data.frame(t(proteinData)))
    PD$PatientID = PatientID
    CD = CD_clean %>% dplyr::select(PatientID, GR_fatcor)
    #ClinDomit$designMatrix = ClinDomit$data
    CD = CD %>% filter(PatientID %in% rownames(ClinDomit$designMatrix))
    #ClinDomit$designMatrix
    PD = dplyr::left_join(PD, CD) %>% 
      dplyr::select(c("PatientID", 
                      #"less.than.or.equal.to_61"
                      GR_fatcor, 
                      #!! parse_expr(GR_fatcor), 
                      limmaResult[s_up,]$Gene.name)) %>% 
      gather(key = "Gene name", value = "log2 of LFQ value", c(limmaResult[s_up,]$Gene.name))
    if (is.null(labelColBox) | labelColBox == "" | labelColBox =="none") {
      groups = NULL
      labels = PD$PatientID
    } else if (!is.null(labelColBox) && labelColBox != "PatientID") {
      message(labelColBox)
      boxPlotSamples = tibble(PatientID = PD$PatientID)
      #boxPlotSamples$PatientID = gsub(pattern = "\\.", replacement = " ", x = boxPlotSamples$PatientID)
      boxPlotSamples$PatientID <- parse_factor(boxPlotSamples$PatientID, include_na = FALSE)
      coloringFactor = CD_clean[,c("PatientID", labelColBox)]
      groups = pull(right_join(coloringFactor, boxPlotSamples)[,2])
      labels = PD$PatientID
    } else {
      message("Please select a different parameter for colouring")
      groups = NULL
      labels = PD$PatientID
    }
    scatterPlot = ggplot(PD, aes(x = get(GR_fatcor), y = `log2 of LFQ value`)) + 
      geom_point(aes(x = get(GR_fatcor), y = `log2 of LFQ value`, colour = groups)) + 
      geom_smooth(method = "lm") +
      labs(x = paste(GR_fatcor)) +
      facet_wrap(~`Gene name`) +
      theme_light()
    if (input$showLabels == TRUE) {
      scatterPlot = scatterPlot + geom_text_repel(aes(label=labels), show.legend = F, size = 4)
    }
    if (is.numeric(groups)){
      scatterPlot  = scatterPlot + scale_color_continuous_tableau() 
    } else {
      scatterPlot  = scatterPlot + scale_color_tableau()
    }
    scatterPlot
  }
  
  
  createBoxPlot <- function(rows_selected, proteinData, filter_GR_fatcor, ClinData, experimentalDesign, GR_fatcor, limmaResult, labelColBox){
    s_up = rows_selected
    PatientID = colnames(proteinData)
    PD = as_tibble(as.data.frame(t(proteinData)))
    PD$PatientID = PatientID
    
    PD = dplyr::left_join(PD, experimentalDesign) %>% 
      dplyr::select(c("PatientID", 
                      GR_fatcor, 
                      limmaResult[s_up,]$Gene.name)) %>% 
      gather(key = "Gene name", value = "log2 of LFQ value", c(limmaResult[s_up,]$Gene.name))# %>% 
    PD[,GR_fatcor] = parse_factor(dplyr::pull(PD, GR_fatcor), include_na = FALSE)
    if (is.null(labelColBox) | labelColBox == "" | labelColBox =="none") {
      groups = NULL
      labels = PD$PatientID
    } else if (!is.null(labelColBox) && labelColBox != "PatientID") {
      message(labelColBox)
      boxPlotSamples = tibble(PatientID = PD$PatientID)
      boxPlotSamples$PatientID <- parse_factor(boxPlotSamples$PatientID, include_na = FALSE)
      coloringFactor = ClinData[,c("PatientID", labelColBox)]
      groups = pull(right_join(coloringFactor, boxPlotSamples)[,2, drop = FALSE])
      labels = PD$PatientID
    } else {
      message("Please select a different parameter for colouring")
      groups = NULL
      labels = PD$PatientID
    }
    
    BoxPlot = ggplot(PD, aes(x = get(GR_fatcor), y = `log2 of LFQ value`)) + 
      geom_boxplot() +
      geom_jitter(aes(colour = groups)) +
      #geom_text_repel(aes(label=labels), show.legend = F, size = 4) +
      labs(x = paste(GR_fatcor)) +
      facet_wrap(~`Gene name`) +
      theme_light()
    
    if (input$showLabels == TRUE) {
      BoxPlot = BoxPlot + geom_text_repel(aes(label=labels), show.legend = F, size = 4)
    }
    if (is.numeric(groups)){
      BoxPlot  = BoxPlot + scale_color_continuous_tableau() 
    } else {
      BoxPlot  = BoxPlot + scale_color_tableau()
    }
    BoxPlot
  }
  
  output$up <- DT::renderDataTable({
    df<-limmaResult$df_up
    DT::datatable(df[, c("Gene.name","logFC","P.Value","adj.P.Val")],
                  rownames = FALSE,
                  class = 'cell-border stripe',width= '150px',extensions ='Scroller', 
                  options = list(
                    deferRender = TRUE,
                    scrollY = 100,
                    scroller = TRUE,
                    scrollCollapse=TRUE,
                    pageLength = 100, lengthMenu = c(5,10,50,100,200)
                  )
    )
    
  })
  
  output$down <- DT::renderDataTable({
    
    df<-limmaResult$df_down
    DT::datatable(df[, c("Gene.name","logFC","P.Value","adj.P.Val")],rownames = FALSE,
                  class = 'cell-border stripe', 
                  extensions ='Scroller',
                  options = list(
                    deferRender = TRUE,
                    scrollY = 100,
                    scroller = TRUE,
                    scrollCollapse=TRUE,
                    pageLength = 100, lengthMenu = c(5, 10, 50, 100, 200)
                  )
    )
  })
  
  output$reportDataDL <- downloadHandler(
    filename = function(){
      paste('EatomicsData', input$GR_fatcor_gsea, names(ClinDomit$designMatrix)[1] , "vs.", names(ClinDomit$designMatrix)[2], Sys.Date(), '.csv', sep = '') },
    content = function(file) {
      x = list("Upregulated.Proteins" = limmaResult$df_up, 
               "Downregulated.Proteins" = limmaResult$df_down,
               "Limma.ExpDesign" = ClinDomit$designMatrix %>% rownames_to_column("PatientID"),
               "Limma.Setup.Details" = data.frame("imputed Data" = input$imputeForLimma, "eBayesTrend" = "TRUE", "Contrast" = paste(names(ClinDomit$designMatrix)[1]," regulated when compared to ", names(ClinDomit$designMatrix[2]))),
               "Upregulated.GeneSets" = gsea_regul$df_up,
               "Downregulated.GeneSets" = gsea_regul$df_down,
               "Differential.GSEA.Setup" = reportBlocks$gseaSetup, 
               "ProteinIDs_Gene_Mapping" = reportBlocks$ProteinIDMap)
      openxlsx::write.xlsx(x, file, row.names = FALSE)
      dev.off()
    }
  )
  
  
  output$doc1 <- renderUI({
    req(ClinDomit$designMatrix, limmaResult$gene_list)
    
    kable_input<-kable(ClinDomit$designMatrix)%>% kable_styling(bootstrap_options = "striped") 
    
    #kable_input<-kable(CD)%>%kable_styling() 
    
    reportBlocks$ExpSetup <- 
      HTML(markdownToHTML(fragment.only=TRUE, text=c( 
        #"* Input MaxQuant File:",protfile$name,
        #"* Input Clinicaldata File:",clinfile$name$name,
        "* Clinical grouping factor:", input$GR_fatcor_gsea,
        "* Two groups to compare:", input$levels[1],",", input$levels[2],
       # "* [Surrogate variable](http://bioconductor.org/packages/release/bioc/html/sva.html) removed:", input$remove_sv,
       # "+ Number of surrogate variables: ",surrogat$num,
      #  "* Impute missing values:",input$imputeForLimma,
        "* Apply filters:", input$expandFilter_gsea,
        " + Clinical Factor to filter: ", input$filter_GR_fatcor,
        " + One group to filter:", input$filter_levels,
        "* The samples contributing to limma are:","total",nrow(ClinDomit$designMatrix),";",input$levels[1],"(", sum(ClinDomit$designMatrix[1]),")","and",input$levels[2],"(", sum(ClinDomit$designMatrix[2]),")", 
        #kable_input,
        scroll_box(kable_input,width = "70%", height = "200px"),
        "* Total number of genes in limma are:",nrow(limmaResult$gene_list),
        "* The number of genes **upregulated** and **downregulated** in ",input$levels[1], " are ", nrow(limmaResult$up)," and " ,nrow(limmaResult$down), " respectively. ",
        "* The threshold used to highlight significant genes is [BH corrected](https://www.rdocumentation.org/packages/stats/versions/3.5.2/topics/p.adjust) adjusted P value of" , input$adj.P.Val, "and absolute log fold change of ",input$logFC 
        
      )))
    bsCollapsePanel(p("Detailed description",style = "color:#18bc9c"),
                    reportBlocks$ExpSetup)
  })
  
}
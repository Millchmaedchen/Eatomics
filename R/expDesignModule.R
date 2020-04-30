expDesignModule_UI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
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
    mainPanel(textOutput(ns("TO_Hello_user")))
  )
}

# Function for module server logic
expDesignModule <- function(input, output, session, measurementFile = NULL, ClinData = NULL) {
  ns <- session$ns
  
  ClinColClasses <- reactive({
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
    updateCheckboxInput(session, "expandFilter", value = FALSE)
  })
  
}
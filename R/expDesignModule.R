expDesignModule_UI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      uiOutput(ns("conditional_grouping_limma"))
    #   uiOutput("conditional_subselectGR_gsea"),
    #   checkboxInput("ContinChoice_gsea", "Use continuous response instead of grouping", FALSE),
    #   #checkboxInput("imputeFor_gsea", "Impute missing values", FALSE),
    #   #checkboxInput("remove_sv", "Remove surrogate variables", FALSE),
    #   checkboxInput("includeCovariates_gsea", "Include parameters as covariates", FALSE), 
    #   conditionalPanel("input.includeCovariates_gsea == TRUE",
    #                    uiOutput("covariatesChoice_gsea")),
    #   checkboxInput("expandFilter_gsea", "Stratification and filter",  FALSE),
    #   conditionalPanel("input.expandFilter_gsea == TRUE",
    #                    uiOutput("filter_group_gsea")),
    #   conditionalPanel("input.expandFilter == TRUE",
    #                    uiOutput("filter_level_gsea")),
    #   conditionalPanel("input.expandFilter_gsea == TRUE",
    #                    uiOutput("selectContrast_gsea")),
    #   actionButton(ns("analyze_diff_gsea"),"Analyze",class = "btn-primary")
     ), 
    mainPanel(textOutput(ns("TO_Hello_user")))
  )
}

# Function for module server logic
expDesignModule <- function(input, output, session, measurementFile = NULL) {
  
  output$conditional_grouping_limma<- renderUI({
    selectInput(
      inputId = "GR_fatcor", 
      label = strong("Select the clinical grouping factor"),
      choices = as.list(colnames(ClinData())),
      multiple = FALSE,
      selectize = TRUE
    )
  })
  
  output$TO_Hello_user <- renderText({
      return(paste(measurementFile))
  })
  
}
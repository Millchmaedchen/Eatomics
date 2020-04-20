## The code is developed by Milena Kraus and Mariet Stephen. 

# This is a Shiny web application (ui + server) for quantitative proteomics data analysis.
# You can run the application by clicking
# the 'Run App' button above.

# Install or load neccessary R packages 
#install.packages("pacman") # for easy install and update #+#
library(pacman) #+#

setRepositories(ind = c(1:6, 8))

#p_load("rstudioapi")

# Shiny libraries
p_load("shiny") #+#
p_load("shinythemes")
p_load('shinycssloaders')
p_load("shinyalert")
p_load('shinyFiles')
p_load("shinyWidgets")
p_load("shinyBS")

# Load data
#p_load("readxl")
#p_load("readr")
p_load(openxlsx)

# Visualization
p_load("pheatmap")
p_load("RColorBrewer")
p_load("plotly")
p_load("ggplot2")
p_load("gridExtra")
p_load("ggthemes")
#install_github("vqv/ggbiplot")
library(ggbiplot)
#p_load("ggiraph")
p_load("autoplotly")
#p_load("EnhancedVolcano")
p_load(kableExtra)
p_load("ggrepel")
p_load("gtools")


# Tidy data
p_load("tidyverse")
#p_load("dplyr")
#p_load("data.table")
#p_load(DT)
p_load("janitor")
#p_load(broom)

# Analysis logic
#p_load("sva")
p_load("imputeLCMD")
p_load("modelr")
p_load("limma")

# Report
p_load(markdown)

## no clue
#p_load("stringi")
#p_load("rgl")
#p_load("crosstalk")
#p_load("markdown")
#p_load(xtable)
#p_load("rlang")
#p_load(nlme)


# Needed for server setup
#p_load("RJDBC")
#p_load("devtools")

# Load non-reactive helper functions 
homeDir = getwd()
source(paste(homeDir, '/helpers.R', sep = ""))

# Load dependency on ssGSEA algorithm 
source(paste(homeDir, '/ssGSEA_PSEA.R', sep = ""))

ui <- fluidPage( 
  # Application title
  navbarPage("Eatomics",id="id",
             theme = shinythemes::shinytheme("flatly"),
             #tabPanel("Introduction"),
             tabPanel("Load and Prepare",
                      sidebarLayout(
                        sidebarPanel(
                          tags$div(title="Load proteinGroups.txt file",
                                   fileInput('file1',
                                             'ProteinGroups.txt',
                                             accept=c('text/csv',
                                                      'text/comma-separated-values,text/plain',
                                                      '.csv'))
                          ),
                          radioButtons("insty", "Quantification type",
                                       choices = c("LFQ" = "LFQ","iBAQ" = "iBAQ"),
                                       selected = "LFQ"),
                          uiOutput("filt"),
                          numericInput("filter_pro", "Define in how many samples a protein needs to be detected in to be included in the analysis",
                                       value = 1,
                                       min = 1,
                                       max = NA, 
                                       step = 1),
                          radioButtons("unique", "Choose method for meaningful gene names",
                                       choices = c("Summarize isoforms" = "chck_iso","Unique names for duplicate genes" = "mk_un"),
                                       selected = "mk_un"),
                          selectizeInput("imputation",
                                         "Imputation method",
                                         choices =  c("perseus-like", "QRILC", "MinDet", "knn"),
                                         selected = "perseus-like"),
                          #  p(a("Detailed information link ",
                          #       href = "https://www.rdocumentation.org/packages/MSnbase/versions/1.20.7/topics/impute-methods",
                          #       target="_blank")),

                          tags$div(title="Load the sample description file",
                                           fileInput('ClinD',
                                                     'ClinicalData.txt',
                                                     accept = c('text/csv',
                                                                'text/comma-separated-values,text/plain',
                                                                '.csv',
                                                                '.tsv'))
                          ),
                          br(),
                          actionButton("analyze","Analyze",class = "btn-primary")),
                        
                        mainPanel(tabsetPanel(id= "QCtab",type = "tabs",
                                              tabPanel(title = "PCA",
                                                       selectizeInput("PCs",
                                                                      label = "Choose the principal components for plotting",
                                                                      choices = list("PC1" = 1, "PC2" = 2, "PC3" = 3, "PC4" = 4, "PC5" = 5, "PC6" = 6, "PC7" = 7, "PC8" = 8),
                                                                      selected = c(1, 2), 
                                                                      multiple = TRUE, 
                                                                      options = list(maxItems = 2)
                                                       ),
                                                       uiOutput("labelCol"),
                                                       checkboxInput("imputeforPCA", "Use imputed data for PCA", TRUE),
                                                       plotOutput("pca_input_samples", height = 800),
                                                       downloadButton('downloadpca', 'Save')
                                              ),

                                              tabPanel(title = "Distribution overview",
                                                       plotOutput("distributionPlot",height = 800) ,
                                                       downloadButton('downloadDistributionPlot', 'Save')
                                              ),
                                              tabPanel(title = "Protein coverage",
                                                       plotOutput("numbers", height = 800),
                                                       br(),
                                                       downloadButton('downloadNumbers', 'Save')
                                              ),
                                              #tabPanel(title = "Sample coverage",
                                              #         plotOutput("coverage", height = 600),
                                              #         downloadButton('downloadCoverage', 'Save')
                                              #),
                                              
                                              # tabPanel(title = "Missing values - Heatmap",
                                              #          plotOutput("missval", height = 800),
                                              #          downloadButton('downloadMissval', 'Save')
                                              # ),
                                              # tabPanel(title = "Missing values - Quant",
                                              #          plotOutput("detect", height = 600),
                                              #          downloadButton('downloadDetect', 'Save')
                                              # ),
                                              # 
                                              # tabPanel(title = "Imputation",
                                              #          plotOutput("imputation", height = 600),
                                              #          downloadButton('downloadImputation', 'Save')
                                              # ),
                                              
                                              tabPanel(title = "Sample to sample heatmap",
                                                       selectizeInput(inputId = "distanceMetric",
                                                                      label = "Select the (dis-)similarity metric",
                                                                      choices = list("Pearson" = "Pearson", "Euclidean" = "Euclidean")
                                                       ),
                                                       plotOutput("StS_heatmap", height = 600),
                                                       downloadButton('downloadStS_heatmap', 'Save')
                                              ),
                                              tabPanel(title = "Cumulative protein intensities",
                                                       plotOutput("CumSumPlot", height = 600),
                                                       downloadButton('downloadCumSumPlot', 'Save')
                                              )
                        )   
                        )
                      )
             ),
             tabPanel("Differential Expression",
                      #actionBttn("tour_limma", icon("info"),color = "success",style = "material-circle",size = "xs"
                      #),
                      sidebarLayout(       
                        sidebarPanel(
                          uiOutput("conditional_grouping_limma"),
                          uiOutput("conditional_subselectGR_limma"),
                          checkboxInput("ContinChoice", "Use continuous response instead of grouping", FALSE),
                          checkboxInput("imputeForLimma", "Impute missing values", FALSE),
                          #checkboxInput("remove_sv", "Remove surrogate variables", FALSE),
                          checkboxInput("includeCovariates", "Include parameters as covariates", FALSE), 
                          conditionalPanel("input.includeCovariates == TRUE",
                                           uiOutput("covariatesChoice")),
                          checkboxInput("expandFilter", "Stratification and filter",  FALSE),
                          conditionalPanel("input.expandFilter == TRUE",
                                           uiOutput("filter_group_limma")),
                          conditionalPanel("input.expandFilter == TRUE",
                                           uiOutput("filter_level_limma")),
                          conditionalPanel("input.expandFilter == TRUE",
                                           uiOutput("selectContrast")),
                          actionButton("analyzeLimma","Analyze",class = "btn-primary")
                          
                        ),
                        
                        mainPanel(
                          fluidRow(
                            column(8,
                                   plotlyOutput("limma"
                                                ,height = 500
                                   ) %>% withSpinner()
                            ),
                            
                            absolutePanel(
                              top=5,
                              right = 20,
                              draggable = TRUE,
                              wellPanel(
                                #sidebarPanel(position = "right",
                                numericInput("adj.P.Val", "Adjusted P value threshold", 0.05, min = 0, max = 1, step= 0.01),
                                sliderInput("logFC", "Log Fold Change", 0, min = 0, max = 10, step = 0.1)
                              )
                            ),
                            br(),
                            column(12,
                                   DT::dataTableOutput('up'),
                                   DT::dataTableOutput('down'),
                                   br(),
                                   uiOutput("labelColBox"),
                                   checkboxInput("showLabels", "Blend in PatientID", FALSE),
                                   plotOutput("boxPlotUp"),
                                   plotOutput("boxPlotDown"),
                                   htmlOutput("doc1")
                            )),
                          
                          br(),
                          downloadButton("report", "Generate report"),
                          downloadButton("reportDataDL", "Download report data")
                        )
                      )
                      
             ),
             tabPanel("ssGSEA", 
                      # actionBttn("tour_ssGSEA", icon("info"),color = "success",style = "material-circle",size = "xs"
                      #),
                      sidebarLayout(
                        sidebarPanel(
                          tags$strong("Change the parameters & hit the analyze button  ", style="color:#18bc9c"),
                          selectInput(
                            inputId = "gs.collection", 
                            label = strong("MSigDb Gene Set Collection"),
                            choices = names(gene.set.databases)
                            #selected = "H-Hallmark50"
                          ),
                          
                          selectInput("sample.norm.type",
                                      label = "Select a normalization method",
                                      choices = c("rank", "log", "log.rank", "none"),
                                      selected = 1  
                          ),
                          numericInput("weight",
                                       label = "Select a Weight (0 to 1)",
                                       value = 0.75,
                                       min = 0, max = 1
                          ),
                          selectInput("statistic",
                                      label = "Select test statistic",
                                      choices = c("area.under.RES", "Kolmogorov-Smirnov"),
                                      selected = 1
                          ),
                          selectInput("output.score.type",
                                      label = "Select enrichment score type",
                                      choices = c("ES", "NES"),
                                      selected = 2
                          ),
                          numericInput("nperm", 
                                       label = "Enter the Number of Permutations",
                                       value = 1000
                          ),
                          numericInput("min.overlap",
                                       label = "Select the minimum overlap between gene set and data",
                                       value = 5
                          ),
                          selectInput("correl.type",
                                      label = "Select correlation type", 
                                      choices = c("rank", "z.score", "symm.rank"),
                                      selected = 1
                          ),
                          uiOutput("output.prefix"),
                          
                          tags$div(title= "Specify the type of analysis from above, then press the analyze button",
                                   actionButton("goButton", "Analyze",class = "btn-primary")         
                          )
                        ),
                        mainPanel(
                          helpText("Make sure to upload proteinGroups.txt file before running ssGSEA"),
                          bsCollapsePanel(p("Detailed description",style = "color:#18bc9c"),
                                          HTML(markdownToHTML(fragment.only=TRUE, text=c( 
                                            "* Single-sample GSEA [(ssGSEA)](http://software.broadinstitute.org/cancer/software/genepattern/modules/docs/ssGSEAProjection/4) is an extension of conventional  Gene Set Enrichment Analysis (GSEA),
                                        developed by Broad insitute<sup>1</sup>.",
                                            "* ssGSEA version used: v4",
                                            "* MSigDB version used: v6.1 ",
                                            "\n[1] Krug, K., et al., A Curated Resource for Phosphosite-specific Signature Analysis. Mol Cell
                                                Proteomics, 2019. 18(3): p. 576-593."
                                          )))),
                          #      ssGSEA to calculate separate pathway enrichment scores for each pairing of a sample and geneset.
                          #      Each ssGSEA enrichment score represents the degree to which the genes in a particular gene set are coordinately up- or down-regulated within a sample.")),
                          
                          useShinyalert()
                        )
                      )
             ),

             
             #  navbarMenu("Help",icon = icon("info-circle"),
             #             tabPanel("Quick Tour"),
             #             tabPanel("Go to vignette",
             #                     
             #                       sprintf("window.open('http://bioconductor.org/packages/%s/bioc/vignettes/iSEE/inst/doc/basic.html', '_blank')")
             #                      
             #             )
             #  )
             tabPanel("About",icon = icon("info-circle"),
                      # actionButton(
                      #   "tour_firststeps", "Click me for a quick tour",
                      #   icon("hand-o-right")
                      #   # style=.actionbutton_biocstyle
                      # ),
                      uiOutput("markdown")
             )
             
             
             #  navbarMenu("Report",icon = icon("fas fa-save"),
             #             tabPanel(downloadButton("report", "Final report",class="butt")),
             #             tags$head(tags$style(".butt{color: white !important;}")),##font color to white
             #             tabPanel(downloadButton("reportDataDL", "limma report",class="butt"))
             #  )
             
  )
)

## Server function
server <- function(input, output, session) {
  reportBlocks <- reactiveValues()
  #push the filesize upload limit
  options(shiny.maxRequestSize = 10000*1024^2, expressions = 500000)
  
  protfile <- reactiveValues()
  
  ###1 Load n Prep tab  
  
  #read .csv uploaded file  
  volumes <-c(root='./Data')
  shinyFileChoose(input, 'files', root=c(root='./Data'), filetypes=c('', 'txt')) 
  data <- reactive({
   # if (input$dataUpload == 'userFile'){
      if (is.null(input$file1)) {
        return(NULL)
      } else{
        inFile <- input$file1 
      }
    #} else if (input$dataUpload == 'serverFile') {
     # if (is.null(input$files)){
    #    return(NULL)
    #  } else {
    #    inFile <-parseFilePaths(volumes, input$files)
    #  }
      
    #}
    if (length(inFile$datapath) == "0")
      return(NULL)
    
    protfile$protfile <- inFile
    proteinGroups = read_tsv(inFile$datapath, col_types = cols(Reverse = "c", `Potential contaminant` = "c", `Only identified by site` = "c"))
    stats_proteinGroups = NULL
    stats_proteinGroups$NumFullProt = nrow(proteinGroups)
    stats_proteinGroups$NumPotCon = proteinGroups %>% 
      dplyr::filter(`Potential contaminant` == "+") %>% nrow()
    stats_proteinGroups$IdentifiedBySite  = proteinGroups %>% 
      dplyr::filter(`Only identified by site` == "+") %>% nrow()
    stats_proteinGroups$Reverse = proteinGroups %>% 
      dplyr::filter(Reverse == "+") %>% nrow()
    
    proteinGroups = proteinGroups %>% 
      dplyr::filter(is.na(`Potential contaminant`)) %>%
      dplyr::filter(is.na(`Only identified by site`)) %>%
      dplyr::filter(is.na(Reverse))
    
    stats_proteinGroups$Cleaned = proteinGroups %>% nrow()
    reportBlocks$stats_proteinGroups = stats_proteinGroups
    
    message("Loaded ", stats_proteinGroups$NumFullProt, " rows from the provided proteinGroups-file. ", 
            stats_proteinGroups$NumPotCon, " potential contaminants, ", stats_proteinGroups$IdentifiedBySite, 
            " proteins only identified by site, and ", stats_proteinGroups$Reverse, " 
            proteins only identified in the reverse database were removed from the data set. 
            The final dataset has ", stats_proteinGroups$Cleaned, " rows. " )
    proteinGroups
  })
  
  #Select columns with specific intensity prefix (LFQ/iBAQ/..)
  insty <- reactive({
    input$analyze
    data<- data()
    reportBlocks$linesFromMaxQuant = nrow(data)
      isolate(selectProteinData(data, intensityMetric = input$insty))
    
  })
  
  #Remove user defined columns/samples
  ## Filter samples UI elements
  output$filt = renderUI({
    if (!is.null(data())){
      insty<-insty()
      selectizeInput("filt",
                     "Exclude columns (samples)",
                     choices=colnames(insty),
                     multiple = TRUE,
                     selected=input$filt)
    } else {
      selectizeInput("filt",
                     "Exclude columns (samples)",
                     choices=colnames(insty))
    }
  })
  
  ## Filter user defined columns/samples
  filt <- reactive({
    input$analyze
    filt_col<-isolate(insty())
    if(is.null(input$filt)){
      filt <- isolate(filt_col)
    } else {
      filt <- isolate(filt_col[ , -which(names(filt_col) %in% input$filt)])
    }
  })
  
  # Filter proteins that were not detected in at least a user defined amount of samples
  filtpro <- reactive({
    input$analyze
    filt_pro <- isolate(filt())
    filterTH <- input$filter_pro-1
    filtpro <- isolate(filterProteins(filt_pro, filterTH = filterTH))
  })
  
  #### new universal format for the protein intensity data that should be feeded to any other requiring method is specified as
  # being a tibble
  # containing log2 transformed, normalized intensity values as specified by the user 
  # one column "Gene names" that harbours uniqe gene names and no NA values
  # the new format should be accepted by any method and also returned by any further method
  # the new format should be saved in a reactiveValues object, one containing the imputed and one the not imputed values
  
  proteinAbundance <- reactiveValues()
  
  #check for isoforms or make unique gene names and log2 transformation
  # needed parameters: 
  # input$unique
  # reactiveValue on log2-transform TRUE/FALSE
  # reactiveValue on normalization (only needed when iBAQ values are used)
  log2tansform <- reactiveVal(TRUE)
  normalizeVsn <- reactiveVal(FALSE)
  
  uniq_names <- observeEvent(c(
    input$analyze, 
    input$insty,
    input$unique, 
    input$filterTH,
    input$filt,
    normalizeVsn,
    log2tansform), ignoreNULL = TRUE,  ignoreInit = TRUE, 
    {
      req(protfile$protfile$name)
      
      filt<-filtpro()
      if (isolate(input$unique == "chck_iso" )) {
        data_unique <- checkForDuplicates(filt)
      }
      if (input$unique == "mk_un" ){
        filt$`Gene names` <- make.unique(filt$`Gene names`)
        filt = as.data.frame(filt)
        filt$`Gene names`[is.na(filt$`Gene names`)] = filt$`Majority protein IDs`[is.na(filt$`Gene names`)]
        data_unique <- filt
        reportBlocks$ProteinIDMap = data_unique %>% dplyr::select(c("Majority protein IDs", "Protein IDs", "Gene names"))
        data_unique = data_unique %>% dplyr::select(-c("Majority protein IDs", "Protein IDs"))
      }
      
      rn = data_unique %>% dplyr::select("Gene names")
      LFQ_columns = grep("Gene names", colnames(data_unique), value = FALSE, invert = TRUE )
      data_unique[data_unique == 0] <- NA
      if(log2tansform() == TRUE && protfile$protfile$name != "proteinGroups.freeze.txt") {
        data_unique = log2(data_unique[,LFQ_columns])
        if(normalizeVsn() == TRUE){
          data_unique = normalizeVSN(data_unique[,LFQ_columns])
        }
      } else if (log2tansform() == FALSE) {
        data_unique = data_unique[,LFQ_columns]
      }
      else {
        data_unique = data_unique[,LFQ_columns]
      }
      proteinAbundance$original = bind_cols(rn, data_unique)
    })
  
  impute_uniques <- observeEvent(c(proteinAbundance$original,
                                   input$imputation), ignoreNULL = TRUE, ignoreInit = TRUE,{
                                     
                                     imputed = proteinAbundance$original %>% column_to_rownames("Gene names") %>% as.data.frame()
                                     if (input$imputation == "perseus-like"){
                                       imputed[[1]] = replaceMissingFromGaussian(imputed)
                                     } else if (input$imputation == "QRILC"){
                                       imputed <- impute.QRILC(as.matrix(imputed))
                                     } else if (input$imputation == "MinDet"){
                                       imputed[[1]]<- impute.MinDet(as.matrix(imputed))
                                     } else if (input$imputation == "knn"){
                                       imputed <- impute.knn(as.matrix(imputed))
                                     }
                                     proteinAbundance$imputed = imputed[[1]] %>% as_tibble(., rownames = "Gene names")
                                   })
  
  
  
  
  #filter and plot the results 
  QCreport<-reactiveValues()
  
  
  observeEvent({
    input$analyze
  }, { 
    
    output$labelCol<- renderUI({
      conditionalPanel(condition = need(ClinData(), FALSE) , 
                       selectInput(
                         inputId = "labelCol", 
                         label = strong("Choose the clinical parameter for group colours"),
                         choices = as.list("none" = "none", colnames(ClinData())),
                         multiple = FALSE,
                         selectize = TRUE
                       )
      )
    })
    
    pca_input_samples <- reactive({
      validate(need(data(), "Please upload a proteinGroups.txt file first."))
      log2tansform(TRUE)
      if(input$insty == "iBAQ") {normalizeVsn(TRUE)}
      if (input$imputeforPCA == TRUE){
        check_table <- proteinAbundance$imputed %>% column_to_rownames("Gene names") %>% as.data.frame()
      } else {
        check_table <- proteinAbundance$original %>% column_to_rownames("Gene names") %>% as.data.frame()
      }
      check_table = t(as.matrix(check_table))
      check_table <- check_table[ , which(apply(check_table, 2, var) != 0)]
      message("calculating prcomp (again)")
      cc<-prcomp(na.omit(check_table),
                 scale. =TRUE, 
                 center = TRUE
      )
    })
    
    ## Remove? 
    PCs <- reactiveVal(input$PCs)
    
    output$pca_input_samples <- renderPlot({
      validate(
        need(length(input$PCs) == 2, "Please select exactly two principal components for plotting.")
      )
      cc <- pca_input_samples()
      if (is.null(input$labelCol) | input$labelCol == "" | input$labelCol=="none") {
        groups = NULL
        labels = rownames(cc$x)
        
      } else if (!is.null(input$labelCol) && input$labelCol != "PatientID") {
        message(input$labelCol)
        pcaSamples = tibble(PatientID = rownames(cc$x))
        pcaSamples$PatientID = gsub(pattern = "\\.", replacement = " ", x = pcaSamples$PatientID)
        pcaSamples$PatientID <- parse_factor(pcaSamples$PatientID, include_na = FALSE)
        coloringFactor = ClinData()[,c("PatientID", input$labelCol)]
        groups = pull(right_join(coloringFactor, pcaSamples)[,2])
        labels = NULL
      }  else {
        message("Please select a different parameter for colouring")
        groups = NULL
        labels = rownames(cc$x)
      }
      ellipse = TRUE
      if (is.numeric(groups)){
        ellipse = FALSE
      }
      
      biplot <- ggbiplot::ggbiplot(cc, choices = as.numeric(input$PCs), 
                         var.scale=1, 
                         obs.scale=1, 
                         var.axes=F, 
                         scale = 1, 
                         groups = groups,
                         labels = labels,
                         ellipse = ellipse
      ) 
      biplot <- biplot + 
        ggtitle("Principal component analysis") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        #scale_fill_continuous(na.value="white") + 
        theme_light()
      if (is.numeric(groups)){
        biplot  = biplot + scale_color_continuous_tableau("Red-Gold",na.value = "grey30") 
      } else {
        biplot  = biplot + scale_color_tableau()
      }
      QCreport$pca = biplot
      biplot
    })
    
    output$downloadpca <- downloadHandler(
      filename = "pca.pdf",
      # filename = function() {
      #   paste('pca-', fileName ,Sys.Date(), '.pdf', sep='') },
      content = function(file) {
        pdf(file)
        print(QCreport$pca)
        dev.off()
      })
    
    
    # #Normalization plot
    # norm_input <- reactive({
    #   plot_normalization_new(data_se_parsed(),
    #                          norm())
    # })
    # output$norm <- renderPlot({
    #   QCreport$norm<- ggsave(norm_input(),height = 800 + (length(norm_input()[[1]])*10))
    #   norm_input()
    # })
    # 
    # output$ui_norm <- renderUI({
    #   plotOutput("norm", height = 800 + (length(norm()[[1]])*10))
    # })
    # 
    # inFile <- input$file1
    # output$downloadNorm <- downloadHandler(
    #   filename = function() {
    #     paste('normalization-', fileName ,Sys.Date(), '.pdf', sep='') },
    #   content = function(file) {
    #     pdf(file)
    #     print(norm_input ())
    #     dev.off()
    #   })
    
    #Distribution plot
    
    distributionPlot_input <- reactive({
      original = proteinAbundance$original %>% column_to_rownames("Gene names") %>% as.data.frame()
      plot_distribution(original)
    })
    
    output$distributionPlot <- renderPlot({
      QCreport$distributionPlot<-distributionPlot_input()
      distributionPlot_input()
    })
    
    output$downloadDistributionPlot <- downloadHandler(
      filename = function() {
        paste('distributionPlot-', fileName ,Sys.Date(), '.pdf', sep='') },
      content = function(file) {
        pdf(file)
        print(distributionPlot_input())
        dev.off()
      })

    
    #Protein coverage plot
    numbers_input <- reactive({
      original = proteinAbundance$original %>% column_to_rownames("Gene names") %>% as.data.frame()
      plot_proteinCoverage(original)
    })
    output$numbers <- renderPlot({
      QCreport$number<-numbers_input()
      numbers_input()
    })
    output$downloadNumbers <- downloadHandler(
      filename = "ProteinNumbers.pdf",
      content = function(file) {
        pdf(file)
        print(QCreport$number)
        dev.off()
      })
    
    # #Sample Coverage plot (Seemed to not be understood by many and depends on DEP package)
    # coverage_input <- reactive({
    #   plot_coverage_new(norm())
    # })
    # output$coverage <- renderPlot({
    #   QCreport$coverage<-coverage_input()
    #   coverage_input()
    # })
    # output$downloadCoverage <- downloadHandler(
    #   filename = function() {
    #     paste('coverage-', fileName ,Sys.Date(), '.pdf', sep='') },
    #   content = function(file) {
    #     pdf(file)
    #     print(coverage_input())
    #     dev.off()
    #   })
    # 
    
    # -- Depends on DEP
    # #Missing value_heatmap
    # missval_input <- reactive({
    #   plot_missval(norm())
    # })
    # 
    # output$missval <- renderPlot({
    #   QCreport$missval<-missval_input()
    #   # utils::str(missval_input())
    #   missval_input()
    # })
    # output$downloadMissval <- downloadHandler(
    #   filename = function() {
    #     paste('missing_values_heatmap-', fileName ,Sys.Date(), '.pdf', sep='') },
    #   content = function(file) {
    #     pdf(file)
    #     print(missval_input())
    #     dev.off()
    #   })
    
    
    #Missing value _quant plot
    detect_input <- reactive({
      plot_detect(norm())
    })
    output$detect <- renderPlot({
      QCreport$detect<-detect_input()
      detect_input()
    })
    output$downloadDetect <- downloadHandler(
      filename = "missing_values_quant.pdf",
      content = function(file) {
        pdf(file)
        gridExtra::grid.arrange(QCreport$detect)
        dev.off()
      })                         
    
    
    
    #Imputation plot
    imputation_input <- reactive({
      plot_imputation_new(norm(),imp())
    })
    output$imputation <- renderPlot({
      QCreport$impute<-imputation_input()
      imputation_input()
    })
    output$downloadImputation <- downloadHandler(
      filename = "imputation.pdf",
      content = function(file) {
        pdf(file)
        print(QCreport$impute)
        dev.off()
      })
    
    #Sample-to-Sample Heatmap
    StSheatmap_input <- reactive({
      log2tansform(TRUE)
      original = proteinAbundance$original %>% column_to_rownames("Gene names") %>% as.data.frame()
      if (input$distanceMetric == "Pearson") {
        corr = TRUE
      } else {corr = FALSE}
      
      plot_StS_heatmap(original, corr = corr)
    })
    output$StS_heatmap <- renderPlot({
      QCreport$StSDistMetric = input$distanceMetric
      QCreport$StSheatmap<-StSheatmap_input()
    })
    output$downloadStS_heatmap <- downloadHandler(
      filename = "StS_heatmap.pdf",
      content = function(file) {
        pdf(file)
        print(QCreport$StSheatmap)
        print(QCreport$StSDistMetric)
        dev.off()
      })
    
    
    # Cumulative Intensities 
    CumSumPlot_input <- reactive({
      log2tansform(FALSE)
      original = proteinAbundance$original %>% column_to_rownames("Gene names") %>% as.data.frame()
      plot_CumSumIntensities(original)
    })
    output$CumSumPlot <- renderPlot({
      QCreport$cumsum<-CumSumPlot_input()
      grid.draw(CumSumPlot_input())
      
    })
    output$downloadCumSumPlot <- downloadHandler(
      filename = "CumSumPlot.pdf",
      content = function(file) {
        pdf(file)
        grid.draw(QCreport$cumsum)
        dev.off()
      }
    )
  })  
  
  
  
  ###2.Limma
  #Establish reactive Values needed several times in computations
  ClinDomit <- reactiveValues()
  protData <- reactiveValues()
  surrogat<-reactiveValues()
  clinfile <- reactiveValues()
  limmaResult <- reactiveValues(gene_list = FALSE)
  
  # Conditional UI elements
  output$labelColBox<- renderUI({
    conditionalPanel(condition = need(ClinData(), FALSE) , 
                     selectInput(
                       inputId = "labelColBox", 
                       label = strong("Choose the clinical parameter for label colours"),
                       choices = as.list("none" = "none", colnames(ClinData())),
                       multiple = FALSE,
                       selectize = TRUE
                     )
    )
  })
  
  output$conditional_grouping_limma<- renderUI({
    selectInput(
      inputId = "GR_fatcor", 
      label = strong("Select the clinical grouping factor"),
      choices = as.list(colnames(ClinData())),
      multiple = FALSE,
      selectize = TRUE
    )
  })
  
  output$conditional_subselectGR_limma <- renderUI({
    req(ClinData(), ClinColClasses(), input$GR_fatcor)
    if (ClinColClasses()[input$GR_fatcor]=='factor' | ClinColClasses()[input$GR_fatcor]=='logical' ){
      selectizeInput(inputId = "levels",
                     label= "Select two groups to compare",
                     choices = ClinData() %>% pull(input$GR_fatcor) %>% levels(),
                     multiple = TRUE, 
                     options = list(maxItems = 2)
      )
    } else {
      d = ClinData() %>% pull(input$GR_fatcor)
      sliderInput(
        inputId = "num.cutoff",
        label = "Select cutoff to divde numeric value:",
        min = min(d, na.rm = TRUE),
        max = max(d, na.rm = TRUE),
        value = colMeans(ClinData()[input$GR_fatcor], na.rm = TRUE), round = T
      )
    }
  }) 
  
  observeEvent(input$includeCovariates, {
    output$covariatesChoice<- renderUI({
      selectInput(
        inputId = "covariates",
        label = strong("Select factors to include as covariates."),
        choices = as.list(colnames(ClinData())),
        multiple = TRUE
      )
    })
  })
  
  observeEvent(input$expandFilter, {
    output$filter_group_limma<- renderUI({
      selectInput(
        inputId = "filter_GR_fatcor",
        label = strong("Select a second parameter"),
        selected = 3,
        choices = as.list(colnames(ClinData())),
        multiple = FALSE,
        selectize = TRUE
      )
    })
    
    output$filter_level_limma <- renderUI({
      req(input$filter_GR_fatcor)
      if (ClinColClasses()[input$filter_GR_fatcor]=='factor' | ClinColClasses()[input$filter_GR_fatcor]=='logical'){
        selectizeInput(inputId = "filter_levels",
                       label = "Filter: Select groups to include in the analysis",
                       choices = ClinData() %>% pull(input$filter_GR_fatcor) %>% levels(),
                       multiple = TRUE
        )
      } else {
        d = ClinDomit$data %>% pull(input$filter_GR_fatcor)
        sliderInput(
          inputId = "filter_num.cutoff",
          label = "Select cutoff to divide numeric value:",
          min = min(d, na.rm = TRUE),
          max = max(d, na.rm = TRUE),
          value = colMeans(ClinData()[input$filter_GR_fatcor], na.rm = TRUE), round = T
        )
      }
    })
    output$selectContrast <- renderUI({
      req(ClinDomit$mainParameter)
      selectizeInput(inputId = "contrastLevels", 
                     label = "Stratify: Select the two groups you want to calculate the difference on.",
                     #choices = ClinData() %>% pull(mainParameter) %>% levels()
                     choices = ClinDomit$data %>% pull(ClinDomit$mainParameter) %>% levels(),
                     multiple = TRUE, 
                     options = list(maxItems = 2)
      )
    })
  }, ignoreNULL = FALSE, ignoreInit = TRUE)
  
  observe({
    ClinColClasses()[input$GR_fatcor] != "numeric"
    updateCheckboxInput(session, "ContinChoice", value = FALSE)
  })
  observe({
    req(input$GR_fatcor)
    input$GR_fatcor
    updateCheckboxInput(session, "expandFilter", value = FALSE)
  })
  
  
  
  # Logic elements
  ## Load clinical parameters
  
  volumes2 <-c(root='./DemoData')
  shinyFileChoose(input, 'ClinDs', root=c(root='./ClinicalData'), filetypes=c('', 'txt', 'tsv'))
  
  ClinData <- reactive({

      validate(
        need(input$ClinD != "", "Please select a file for upload or choose to use the database connection.")
      )
      clinfile$name <- input$ClinD

    if (length(clinfile$name$datapath) == "0")
      return(NULL)
    ClinData = read_tsv(clinfile$name$datapath, 
                        na =c("", "NA", "N/A","0","<Null>"),
                        skip_empty_rows = TRUE, 
                        locale = locale(decimal_mark = ",") # Todo: uncomment this for US/English seperator
    ) %>% remove_empty("cols") %>% mutate_if(is.character, as.factor) %>% mutate_if(is.logical, as.factor)
    if(protfile$protfile$name == "proteinGroups.freeze.txt"){ ### TODO: Uncomment for public use!
      ClinData = ClinData %>% dplyr::rename(PatientID_DB = PatientID)
      ClinData = ClinData %>% dplyr::rename(PatientID = freezeID_log)
    }
    ClinDomit$data = ClinData %>% clean_names()
    ClinData
  })
  
  ClinColClasses <- reactive({
    df = ClinData()
    #df = ClinDomit$data
    df = lapply(df, class)
  })
  
  ClinColClasses_2 <- reactive({
    df = ClinDomit$data
    df = lapply(df, class)
  })
  
  
  
  # set and manipulate chosen parameter to being cat from cont and form new groups from stratified setups
  observeEvent(c(input$filter_GR_fatcor, 
                 input$GR_fatcor,
                 input$ContinChoice, 
                 input$num.cutoff 
                 )
               , {
                 
      ClinDomit$mainParameter = make_clean_names(input$GR_fatcor)

      ## categorize numeric data - first parameter
      if (input$ContinChoice == FALSE & ClinColClasses_2()[ClinDomit$mainParameter]=='numeric') {
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
        colnames(ClinDomit$data)[colnames(ClinDomit$data) == "categorizedParameter"] = paste(input$GR_fatcor, "cat", sep = "_", collapse = "_") %>% make_clean_names()
        ClinDomit$data = ClinDomit$data[,!duplicated(colnames(ClinDomit$data), fromLast = TRUE)]
        ClinDomit$mainParameter = paste(input$GR_fatcor,  "cat", sep = "_", collapse = "_") %>% make_clean_names()
      }
      if (is.null(need(input$expandFilter, FALSE)) & is.null(need(input$filter_GR_fatcor, FALSE))) {
        
        ## categorize second numeric parameter
        ClinDomit$filterParameter =  make_clean_names(input$filter_GR_fatcor)
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
  ## TODO: Rename Filter to "stratification"
  
  observeEvent(input$analyzeLimma ,{
    validate(need(proteinAbundance$original , "Please upload a proteinGroups file first (previous tab)."))

    mainParameter = make_clean_names(ClinDomit$mainParameter)
    if (!is.null(ClinDomit$filterParameter)) {
      filterParameter = make_clean_names(ClinDomit$filterParameter)
    }
    #    }
    if (input$expandFilter == TRUE & input$ContinChoice == FALSE) {
      req(input$contrastLevels)
    }
    if (is.null(input$contrastLevels)){
      contrastLevels = input$levels
    } else{
      contrastLevels = input$contrastLevels
    }
    
    covariates = input$covariates#, TODO: Set up instead of filter for example
    covariates = make_clean_names(covariates)
    
    if (!is.null(ClinDomit$filterParameter) & ClinColClasses_2()[mainParameter] == "numeric") {
      ClinData = ClinDomit$data %>% 
        clean_names() %>% 
        dplyr::select(patient_id, mainParameter, covariates, all_of(filterParameter)) %>% 
        remove_missing() %>% 
        dplyr::filter(!!sym(filterParameter) %in% !!input$filter_levels) %>% 
        select(-filterParameter)
    } else {
      ClinData = ClinDomit$data %>% 
        clean_names() %>% 
        dplyr::select(patient_id, mainParameter, covariates) %>% 
        remove_missing() 
    }
    ClinColClasses = lapply(ClinData, class)

    if (input$ContinChoice == FALSE & ClinColClasses[ClinDomit$mainParameter]=='numeric')
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
    if (input$ContinChoice == TRUE){
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

    expDesignInst = ClinDomit$designMatrix %>% clean_names() 
    if (input$ContinChoice == FALSE){
      validate(need(sum(ClinDomit$designMatrix[,2])>=3, "The experimental design does not contain three or more samples to test on."))
      validProteins =  proteinAbundance$original[, rownames(expDesignInst)]
      validProteins_1 = validProteins[,expDesignInst[,1]]
      validProteins_2 = validProteins[,expDesignInst[,2]]
      validProteinsLog = (rowSums(is.na(validProteins_1)) < length(validProteins_1)*0.5) + (rowSums(is.na(validProteins_2)) < length(validProteins_2)*0.5)
      
      ## Prepare protein abundance for limma
      #match_Protdata =  match_Protdata[, rownames(designMatrix)]
      if (input$imputeForLimma == TRUE){
        proteinAbundanceLimma <- proteinAbundance$imputed[as.vector(validProteinsLog > 0), c("Gene names", rownames(expDesignInst))] %>% as_tibble()
      } else {
        proteinAbundanceLimma <- proteinAbundance$original[as.vector(validProteinsLog > 0),c("Gene names", rownames(expDesignInst))] %>% as_tibble()
      }
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
    if (input$ContinChoice == FALSE) {
      limmaResult$gene_list <- limma::topTable(fit3, coef=1, number = 1e+09, sort.by="P", adjust.method="BH")
    } else if (input$ContinChoice == TRUE)  {
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

    if (!is.null(input$expandFilter) && input$expandFilter == TRUE) {
      Filnames = paste(" only ", input$filter_levels[1], collapse= "")
    } else {
      Filnames = ""
    }
    
    input$analyzeLimma 
    isolate(
    if (input$ContinChoice == FALSE){
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
    isolate(if(ClinColClasses()[input$GR_fatcor]=='factor' | input$ContinChoice == FALSE) {
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
    } else if(input$ContinChoice == TRUE){
      reportBlocks$scatterPlotUp = createLineScatterPlot(input$up_rows_selected, original, ClinData(), input$GR_fatcor, limmaResult$df_up, input$labelColBox)
      return(reportBlocks$scatterPlotUp)
    })
  })
  
  output$boxPlotDown <- renderPlot({
    validate(need(input$down_rows_selected, FALSE))
    original = proteinAbundance$original %>% column_to_rownames("Gene names") %>% as.data.frame()
    if(ClinColClasses()[input$GR_fatcor]=='factor' | input$ContinChoice == FALSE) {
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
    } else if(input$ContinChoice == TRUE){
      reportBlocks$scatterPlotDown = createLineScatterPlot(input$down_rows_selected, original, ClinData(), input$GR_fatcor, limmaResult$df_down, input$labelColBox)
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
      paste('EatomicsData', input$GR_fatcor, names(ClinDomit$designMatrix)[1] , "vs.", names(ClinDomit$designMatrix)[2], Sys.Date(), '.csv', sep = '') },
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
        "* Input MaxQuant File:",protfile$name,
        "* Input Clinicaldata File:",clinfile$name$name,
        "* Clinical grouping factor:", input$GR_fatcor,
        "* Two groups to compare:", input$levels[1],",", input$levels[2],
        "* [Surrogate variable](http://bioconductor.org/packages/release/bioc/html/sva.html) removed:", input$remove_sv,
        "+ Number of surrogate variables: ",surrogat$num,
        "* Impute missing values:",input$imputeForLimma,
        "* Apply filters:", input$expandFilter,
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
  
  
  
  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    #filename = "report.pdf",
    filename = "report.html",
    
    content = function(filename) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
    
      
      # Set up parameters to pass to Rmd document
      params <- list(
        pca = QCreport$pca,
        TSdetect = QCreport$TSdetect,
        #norm = QCreport$norm,
        number= QCreport$number,
        coverage = QCreport$coverage,
        missval = QCreport$missval,
        detect = QCreport$detect,
        impute = QCreport$impute,
        StSheatmap = QCreport$StSheatmap,
        StSheatmapDistMetric = QCreport$StSDistMetric,
        cumsum =QCreport$cumsum,
        linesFromMaxQuant = reportBlocks$linesFromMaxQuant,
        stats_proteinGroups = reportBlocks$stats_proteinGroups,
        volcano_plot = reportBlocks$volcano_plot,
        boxPlotUp = reportBlocks$boxPlotUp,
        boxPlotDown = reportBlocks$boxPlotDown,
        ExpSetup = reportBlocks$ExpSetup,
        UpRegul = limmaResult$df_up,
        DoRegul = limmaResult$df_down,
        gsea_volcano_plot=reportBlocks$gsea_volcano_plot,
        gseaSetup = reportBlocks$gseaSetup,
        gseaUpRegul = gsea_regul$up,
        gseaDoRegul = gsea_regul$down,
        separateValuesGSEA = reportBlocks$separateValuesGSEA
      )
      
      withProgress(message = 'Generating Report', {
        for (i in 1:15) {
          incProgress(1/15)
          Sys.sleep(0.25)
        }
      })
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(input = "report.Rmd", 
                      output_file = paste('EatomicsReport_', Sys.time(), '.html', sep = ''),
                       output_dir  = getwd(),
                       clean = TRUE,
                       quiet = TRUE, 
                        params = params,
                        envir = new.env(parent = globalenv())
      )
      
      
      #uploadName = paste('EatomicsReport-', Sys.Date(), Sys.time(), '.html', sep = '')
      #drive_upload("report.html", name = uploadName)
    }
  )
  
  ###3 ssGSEA tab 
  
  
  output$output.prefix <- renderUI({ 
    textInput("output.prefix",label = "Insert a Prefix for Output Files", input$gs.collection )
  })
  
  dir.create("EnrichmentScore")
  
  # run the ssGSEA analysis
  #getssgseaObj = eventReactive(input$goButton,{
  observeEvent(input$goButton, {
    
    validate(need(input$files != "", "Please upload proteomics data first (previous tab)."))
    original = proteinAbundance$original %>% column_to_rownames("Gene names") %>% as.data.frame()
    ssgsea_data = as.matrix(original)
    #ssgsea_data= imp_woNorm()
    
    
    withProgress(message = 'Calculation in progress',
                 detail = 'An alert notification will appear upon download of the file', value = 1, {
                   
                   ssgsea_obj = ssGSEA2(input.ds =ssgsea_data, gene.set.databases=gene.set.databases[input$gs.collection], sample.norm.type = input$sample.norm.type,
                                        weight = input$weight,statistic =input$statistic,output.score.type= input$output.score.type, 
                                        #nperm = input$nperm, min.overlap   = input$min.overlap ,correl.type = input$correl.type,par=F,export.signat.gct=T,param.file=T, output.prefix= input$output.prefix)           
                                        nperm = input$nperm, min.overlap   = input$min.overlap ,correl.type = input$correl.type, output.prefix= input$output.prefix)           
                 }) 
    
    shinyalert("File is downloaded", type = "success")                     
  })
  
  
  
  ##download the file
  output$download<- downloadHandler(
    filename = function() { paste('efix, .csv', sep='') },
    content = function(file) {
      write.csv(getssgseaObj(), file)
      #descs=getssgseaObj()$gs.descs.2, filename)
    }
  ) 
  
  
  
  
  ###4.Differential GSEA tab
  #TODO introduce error handling into t-test function
  #TODO enable to use imputation or not but give a hint that results will shrink dramatically if imputation is not used
  
  ClinDnamesGSEA <- reactive({
    rownames(ClinData())
  })
  
  ClinColClassesGSEA <- reactive({
    df = ClinData()
    lapply(df,class)
  })
  
  selected_gs.collection <- reactive({
    input$filterGSEA
    nc<-max(count.fields(paste(gene.set.databases[input$diff.gs.collection], sep="")))
    temp =isolate(t(read.table(paste(gene.set.databases[input$diff.gs.collection]),sep="", col.names=paste("V",1:nc,sep=""),fill=T)[,-2]))
    
    colnames(temp) = temp[1,]
    temp[-1,]
    
  })
  
  selected_gsea.score <- reactive({
    read.delim(paste("EnrichmentScore/", input$diff.gs.collection, sep = ""), stringsAsFactors=F, skip=2, row.names='Name')
  })
  
  output$conditional_groupingGSEA <- renderUI({
    selectInput(
      inputId = "GR_fatcorGSEA", 
      label = strong("Select the clinical grouping factor"),
      choices = as.list( colnames(ClinData())),
      multiple = FALSE,
      selectize = TRUE
      # selected = c(1) 
    )
  })
  
  
  output$conditional_subselectGRGSEA <- renderUI({
    req(ClinData(),ClinColClassesGSEA())
    if (ClinColClassesGSEA()[input$GR_fatcorGSEA]=='factor' | ClinColClassesGSEA()[input$GR_fatcorGSEA]=='logical'){
      selectizeInput(inputId = "levelsGSEA",
                     label= "Select two groups to compare",
                     choices = ClinData() %>% dplyr::pull(input$GR_fatcorGSEA) %>% levels(),
                     #selected = NULL,
                     multiple = TRUE, 
                     options = list(maxItems = 2)
      )
    } else {
      d = ClinData() %>% pull(input$GR_fatcorGSEA)
      sliderInput(
        inputId = "num.cutoffGSEA",
        label = "Select cutoff to divde numeric value:",
        min = min(d, na.rm = TRUE),
        max = max(d, na.rm = TRUE),
        value = colMeans(ClinData()[input$GR_fatcorGSEA], na.rm = TRUE), round = T
      )
    }
  })
  
  observeEvent(input$expandFilterGSEA, {
    output$filter_group_gsea<- renderUI({
      selectInput(
        inputId = "filter_GR_fatcorGSEA",
        label = strong("Select a group to filter on"),
        selected = 3,
        choices = as.list(colnames(ClinData())),
        multiple = FALSE,
        selectize = TRUE
      )
    })
    
    
    output$filter_level_gsea <- renderUI({
      req(input$filter_GR_fatcorGSEA)
      if (ClinColClassesGSEA()[input$filter_GR_fatcorGSEA]=='factor' | ClinColClassesGSEA()[input$filter_GR_fatcorGSEA]=='logical'){
        selectizeInput(inputId = "filter_levelsGSEA",
                       label = "Select one group to include in the analysis",
                       choices =  ClinData() %>% pull(input$filter_GR_fatcorGSEA) %>% levels(),
                       multiple = FALSE
        )
        
      } else #if(ClinColClasses()[,input$filter_GR_fatcorGSEA]=="numeric")
      {
        f = ClinData() %>% pull(input$GR_fatcorGSEA)
        sliderInput(
          inputId = "filter_num.cutoffGSEA",
          label = "Select cutoff to filter:",
          min = min(f, na.rm = TRUE),
          max = max(f, na.rm = TRUE),
          value = colMeans(ClinData()[input$filter_GR_fatcorGSEA], na.rm = TRUE), round = T
        )
      }
    })
  },ignoreNULL = FALSE, ignoreInit = TRUE)
  
  ClinDomitGSEA <- reactiveValues()
  
  observeEvent(input$filterGSEA,{
    if (ClinColClassesGSEA()[input$GR_fatcorGSEA]=='factor' | ClinColClassesGSEA()[input$GR_fatcorGSEA]=='logical'){
      design_table<-cbind(ClinData()%>% dplyr::select(PatientID), ClinData()[input$GR_fatcorGSEA] == input$levelsGSEA[1], ClinData()[input$GR_fatcorGSEA]==input$levelsGSEA[2])
      colnames(design_table)<-c("PatientID", input$levelsGSEA[1], input$levelsGSEA[2])
      ind = which(design_table[input$levelsGSEA[1]] != design_table[input$levelsGSEA[2]])
      testv = design_table[ind,]
      ClinDomitGSEA$data <- testv
    } else{
      design_table<-as.data.frame(cbind(ClinData() %>% dplyr:: select(PatientID), ClinData()[input$GR_fatcorGSEA]<= input$num.cutoffGSEA,ClinData()[input$GR_fatcorGSEA]>input$num.cutoffGSEA))
      colnames(design_table) <- make.names(c('PatientID', paste('less than or equal to',input$num.cutoffGSEA, sep='_'), paste('greater than',input$num.cutoffGSEA,sep='_')))
      ClinDomitGSEA$data <- na.omit(design_table)
    }
    
    if (input$expandFilterGSEA==TRUE){
      Ind <- NULL
      if (!is.null(input$filter_levelsGSEA)){
        filterlevels = input$filter_levelsGSEA
        for (i in 1:length(input$filter_levelsGSEA)) {
          Ind = append(Ind, ClinData()[input$filter_GR_fatcorGSEA] == filterlevels[i])
        }
      } else if (!is.null(input$filter_num.cutoffGSEA)) {
        Ind = which(ClinData()[,input$filter_GR_fatcorGSEA] <= input$num.cutoffGSEA)
      }
      if (!is.null(Ind)){
        CD = ClinDomitGSEA$data
        CD =na.omit(CD[Ind,])
        ClinDomitGSEA$data <- CD
        
      }
    }
  })
  
  ERscore_names<-reactiveValues()
  
  allsamples <- reactive({  
    test = selected_gsea.score()
    
    colnames(test) <- gsub(x = names( test),
                           pattern = "\\.",
                           replacement = " ")
    ERscore_table<-test[,-1]
    ERscore_names$data<- colnames(ERscore_table)
    ERscore_table
  })
  
  matchedSampleNamesGSEA <- reactive({
    req(ClinDomitGSEA$data)
    CD = ClinDomitGSEA$data
    l = ERscore_names$data 
    r = as.character(CD[,1])
    match_gsea<-Reduce(intersect, list(l, r))
  })
  
  clinDataMatchGSEA <-reactive({
    req(ClinDomitGSEA$data)
    CD = ClinDomitGSEA$data
    
    #match_Clindata <-CD[na.omit(matchedSampleNamesGSEA()),]
    match_Clindata<-CD %>% filter(PatientID %in% matchedSampleNamesGSEA())
    
  })
  
  # Start sifnicifance analysis
  # contrastMatrix<-eventReactive(input$filter,{
  #    validate(need(data() , "Please upload a proteinGroups file first."))
  #    validate(need(ClinData() , "Please upload a clinical information file first."))
  #    
  #    req(clinDataMatch())
  
  reducedscore<-reactiveValues()
  signGSEAResultsGroups <- eventReactive(input$filterGSEA,{
    validate(need(allsamples(), "There is no enrichment scores file available. Please execute the ssGSEA panel first."))
    validate(need(ClinData() , "Please upload a clinical information file first."))
    scoresTemp <- as.data.frame(allsamples())
    selected_samplesTmp = matchedSampleNamesGSEA()
    CD = ClinDomitGSEA$data
    
    group<-CD %>%
      gather(parameter,value,tail(names(CD), 2)) %>% spread(value,parameter) %>% dplyr::select("PatientID",`TRUE`) %>% setnames("TRUE" ,input$GR_fatcorGSEA)
    
    tmpGroupingFactor = input$GR_fatcorGSEA
    
    scoresTemp = scoresTemp[,-dim(scoresTemp)[2]]
    scoresTemp = as.data.frame(t(scoresTemp))
    scoresTemp = scoresTemp %>% rownames_to_column() 
    
    scoresTemp2 = left_join(scoresTemp, group, by = c("rowname"="PatientID"))
    
    #Remove all rows containing only NA's 
    scoresTemp2 = scoresTemp2 %>% filter_all(all_vars(!is.na(.))) %>% as.data.frame()
    
    #add_column(scoresTemp2, "new" = droplevels(scoresTemp2[,tmpGroupingFactor]))
    scoresTemp2[,length(scoresTemp2)] <- droplevels(as.factor(scoresTemp2[,length(scoresTemp2)]))
    
    reducedscore$data<-scoresTemp2
    
    customT = function( x, a ){
      dummy = bind_cols("x" = x, "a" = a)
      names(dummy) = c("x", "a")
      tidy(t.test(x= filter(dummy, a == paste(levels(a)[1]))$x, y =  filter(dummy, a == paste(levels(a)[2]))$x))
    }
    signTempNew = scoresTemp2 %>% dplyr::select(-c("rowname", tmpGroupingFactor)) %>% 
      apply(., 2, customT, as.factor(scoresTemp2[,length(scoresTemp2)]))
    
    reduced = bind_rows(signTempNew)
    reduced = add_column(reduced, "Gene Set" = names(signTempNew))
    
    reduced$FC<-foldchange( reduced$estimate1,  reduced$estimate2)
    reduced$log_ratio <-foldchange2logratio(foldchange(reduced$estimate1,reduced$estimate2),base=2)
    reduced 
    
  })
  
  # Diff <- reactive({
  #   Ind = which(signGSEAResultsGroups()$p.value<0.05)
  #   DiffNames  =   signGSEAResultsGroups()[Ind,c("Gene Set","p.value")]
  #   list(Ind, DiffNames)
  # })
  
  
  diffgsea_input<- eventReactive(input$filterGSEA,{
    #reactive({
    t.test_list<-signGSEAResultsGroups()
    message(paste("Pathways Enriched in",input$diff.gs.collection,"gene set:",nrow( t.test_list),sep = ' '))
    #gene_list$threshold = as.factor(abs(gene_list$logFC) > input$logFC & gene_list$adj.P.Val < input$adj.P.Val)
    
    t.test_list$threshold = as.factor( t.test_list$p.value < 0.05/nrow(t.test_list))
    
    CDnames = colnames(as.data.frame(clinDataMatchGSEA()))
    
    if (!is.null(input$expandFilterGSEA) && input$expandFilterGSEA == TRUE) {
      Filnames = paste(" only ", input$filter_levelsGSEA[1], collapse= "")
    } else {
      Filnames = ""
    }
    Factor = input$GR_fatcorGSEA
    
    ggplot(data= t.test_list, aes(x=log_ratio, y=-log10(p.value), colour=threshold)) +
      ggtitle(paste(Factor,": ", CDnames[2]," vs ", CDnames[3], Filnames, collapse = "")) +
      geom_point(shape=1, data = t.test_list, aes(text= `Gene Set`), size = 1, alpha = 0.4) +
      xlim(-1,1) + #labs(color = "Adjusted p-value < 0.05") + 
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) + scale_color_tableau() +
      theme_light()+ theme(legend.title = element_blank())
    
  })
  
  
  gsea_regul <- reactiveValues()
  
  output$diffgsea <- renderPlotly( {
    t.test_list<-signGSEAResultsGroups()
    t.test_list$P.adj = t.test_list$p.value*nrow(t.test_list)
    ##select upregulated dataframe
    gsea_up = input$up_gsea_selected
    t.test_list$upregulated = as.factor(t.test_list$log_ratio> 0 & t.test_list$p.value < 0.05/nrow(t.test_list))
    
    gsea_regul$up<-na.omit( t.test_list[ t.test_list$upregulated =="TRUE",])
    gsea_regul$df_up<- gsea_regul$up %>% dplyr::rename("Geneset (upregulated)"= "Gene Set" )
    
    ##select down regulated dataframe
    gsea_down = input$down_gsea_selected
    t.test_list$downregulated = as.factor(t.test_list$log_ratio < 0 & t.test_list$p.value < 0.05/nrow(t.test_list))
    
    gsea_regul$down<-na.omit( t.test_list[ t.test_list$downregulated =="TRUE",])
    gsea_regul$df_down<- gsea_regul$down %>% dplyr::rename("Geneset (downregulated)"= "Gene Set" )
    reportBlocks$gsea_volcano_plot<-diffgsea_input()
    
    legendtitle <- list(yref='paper',xref="paper",y=1,x=1.1, text=paste("Threshold:\nadj.p < 0.05"), showarrow=F)
    pp<-ggplotly(reportBlocks$gsea_volcano_plot, tooltip = "text")%>% layout(legend=list(y=0.89, yanchor="top"), annotations=legendtitle ) 
    
    pp
    if (length(gsea_up)) {
      pp<- reportBlocks$gsea_volcano_plot + geom_point(data = gsea_regul$df_up[gsea_up,], color='red')
      
    }
    else if (length(gsea_down)) {
      pp<- reportBlocks$gsea_volcano_plot + geom_point(data = gsea_regul$df_down[gsea_down,], color='green')
    }
    pp
    
  })
  
  output$pathway_up <- DT::renderDataTable({
    df<-gsea_regul$df_up
    
    DT::datatable(df[, c("Geneset (upregulated)","log_ratio","P.adj")],rownames = FALSE,
                  class = 'cell-border stripe',width= '150px',extensions ='Scroller', options = list(
                    deferRender = TRUE,
                    scrollY = 100,
                    scroller = TRUE,
                    scrollCollapse=TRUE,
                    pageLength = 100, lengthMenu = c(5,10,50,100,200)
                  )
    )
    
  })
  
  output$pathway_down <- DT::renderDataTable({
    
    df<-gsea_regul$df_down
    print(colnames(df))
    DT::datatable(df[, c("Geneset (downregulated)","log_ratio","P.adj")],rownames = FALSE,
                  class = 'cell-border stripe', 
                  extensions ='Scroller',
                  options = list(
                    deferRender = TRUE,
                    scrollY = 100,
                    scroller = TRUE,
                    scrollCollapse=TRUE,
                    pageLength = 100, lengthMenu = c(5,10,50,100,200)
                  )
    )
  })
  
  output$gsea_doc <- renderUI({
    CD = as.data.frame(clinDataMatchGSEA())
    CC= CD[, -1]
    
    kable_input<-kable(CD)%>%kable_styling(bootstrap_options = "striped") 
    # kable_input<-knitr::kable(CD)
    
    reportBlocks$separateValuesGSEA <- list(
      "diff.gs.collection" = input$diff.gs.collection, 
      "gseafile" = gseafile$name$name,
      "GR_fatcorGSEA" = input$GR_fatcorGSEA,
      "groupGSEA_1" = input$levelsGSEA[1],
      "groupGSEA_2" = input$levelsGSEA[2],
      "expandFilterGSEA" = input$expandFilterGSEA,
      "filter_GR_fatcorGSEA" = input$filter_GR_fatcorGSEA,
      "filter_levelsGSEA"= input$filter_levelsGSEA
    )
    
    reportBlocks$gseaSetup <- HTML(markdownToHTML(fragment.only=TRUE, text=c( 
      "* Enriched Pathway Geneset:",input$diff.gs.collection,
      "* Input Clinicaldata File:",gseafile$name$name,
      "* Clinical grouping factor:", input$GR_fatcorGSEA,
      "* Two groups compared:", input$levelsGSEA[1],",", input$levelsGSEA[2],
      #"* Impute missing values:",input$imputeForLimma,
      "* Apply filters:", input$expandFilterGSEA,
      " + Clinical Factor to filter: ", input$filter_GR_fatcorGSEA,
      " + One group to filter:", input$filter_levelsGSEA,
      "* The samples contributing to differential pathway analysis:","total",nrow(CC),";",input$levelsGSEA[1],"(",length(which(CC[1]=="TRUE")),")","and",input$levelsGSEA[2],"(",length(which(CC[2]=="TRUE")),")", 
      #kable_input,
      scroll_box(kable_input,width = "70%", height = "200px"),
      "* Total number of gene sets the comparison was performed on:",nrow(signGSEAResultsGroups()),
      "* The number of pathways **upregulated** and **downregulated** in ",input$levelsGSEA[1], " are ", nrow(gsea_regul$up)," and " ,nrow(gsea_regul$down), " ,respectively. ",
      "* The threshold used to highlight significant genes is [bonferroni corrected](https://www.rdocumentation.org/packages/stats/versions/3.5.2/topics/p.adjust) P value of 0.05" 
      #input$adj.P.Val, "and absolute log fold change of ",input$logFC 
      
    )))
    bsCollapsePanel(p("Detailed description",style = "color:#18bc9c"),
                    reportBlocks$gseaSetup)
  })
  
  
  
  
  ###5. About tabpanel
  
  output$markdown <- renderUI({
    includeHTML("markdown_docs/About.html")
  })
  
  ####tour guide
  
  
  # observeEvent(input$tour_firststeps, {
  #   #if (input$id == "ssGSEA") {
  #   if(is.null(tour)) {
  #     tour
  #   }
  #   introjs(session, options=list(steps=tour))
  #   #}
  # })
  # 
  # observeEvent(input$tour_limma, {
  #   if(is.null(tour_limma)) {
  #     tour_limma
  #   }
  #   introjs(session, options=list(steps=tour_limma))
  #   #}
  # })
  
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
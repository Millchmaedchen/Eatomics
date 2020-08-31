## The code is written by Milena Kraus and Mariet Stephen. 
# This is a Shiny application  for quantitative proteomics data analysis.
# You can run the application by clicking
# the 'Run App' button above.

# Set Repositories for correct installation of all dependencies/packages
setRepositories(ind = c(1:6, 8))

# Install or load neccessary R packages 
list_of_packages = c("shiny","shinythemes", 'shinycssloaders', "shinyalert", "shinyjs", 'shinyFiles', "shinyWidgets", "shinyBS", "openxlsx", "pheatmap", "RColorBrewer", "plotly", "ggplot2",
                     "gridExtra", "ggthemes", "autoplotly", "kableExtra", "ggrepel", "gtools", "verification", "tidyverse", "broom", "imputeLCMD", "modelr", "limma", "markdown", "processx", "pkgload")

lapply(list_of_packages, 
       function(x) if(!require(x,character.only = TRUE)) install.packages(x, dependencies = TRUE))

#install.packages("https://cran.r-project.org/src/contrib/Archive/janitor/janitor_1.2.1.tar.gz", repos=NULL, type='source')

# Load non-reactive helper functions 
homeDir = getwd()
source(paste(homeDir, '/helpers.R', sep = ""))

# Load dependency on ssGSEA algorithm 
source(paste(homeDir, '/ssGSEA_PSEA.R', sep = ""))

# Load experimental design module
source(paste(homeDir, '/expDesignModule.R', sep = ""))

# Load download module
source(paste(homeDir, '/plot_download_module.R', sep = ""))

# Load plot color module
source(paste(homeDir, '/configurePlotCol_module.R', sep = ""))

# Load available gene sets from Data 
gene.set.databases = list.files(path = paste(homeDir, "/../Data/GeneSetDBs/", sep = ""), pattern = ".gmt", full.names = TRUE)
names(gene.set.databases) <- list.files(path = paste(homeDir, "/../Data/GeneSetDBs/", sep = ""), pattern = ".gmt")

#Change for multi-user implementation
#sessionID = paste( "EnrichmentScores_User", idmaker(1), sep = "")
#dir.create(sessionID) # flush directory when session ends
sessionID = paste(homeDir,"/../Data/EnrichmentScores_User", sep = "")


ui <- shiny::fluidPage( 
  useShinyjs(),
  # Application title
  shiny::navbarPage("Eatomics",id="id",
             theme = shinythemes::shinytheme("flatly"),
             #tabPanel("Introduction"),
             shiny::tabPanel("Load and Prepare",
                      shiny::sidebarLayout(
                        shiny::sidebarPanel(
                          tags$div(class = "header", checked = NA,
                                   tags$b("Load proteinGroups.txt file")
                          ),
                          br(),
                                  #shinyFilesButton(id = 'files', 
                                  #                 label='Select demo proteinGroups file', 
                                  #                 title = 'Select demo proteinGroups file', 
                                  #                 multiple=FALSE)
                                   fileInput('file1',
                                             'ProteinGroups.txt',
                                             accept=c('text/csv',
                                                      'text/comma-separated-values,text/plain',
                                                      '.csv'))
                          ,
                          br(),
                          shiny::radioButtons("insty", "Quantification type",
                                       choices = c("LFQ" = "LFQ","iBAQ" = "iBAQ"),
                                       selected = "LFQ"),
                          shiny::uiOutput("filt"),
                          shiny::numericInput("filter_pro", "Define in how many samples a protein needs to be detected in to be included in the analysis",
                                       value = 1,
                                       min = 1,
                                       max = NA, 
                                       step = 1),
                          shiny::radioButtons("unique", "Choose method for meaningful gene names",
                                       choices = c("Summarize isoforms" = "chck_iso","Unique names for duplicate genes" = "mk_un"),
                                       selected = "mk_un"),
                          shiny::selectizeInput("imputation",
                                         "Imputation method",
                                         choices =  c("perseus-like", "QRILC", "MinDet", "knn"),
                                         selected = "perseus-like"),
                          #  p(a("Detailed information link ",
                          #       href = "https://www.rdocumentation.org/packages/MSnbase/versions/1.20.7/topics/impute-methods",
                          #       target="_blank")),
                          tags$div(class = "header", checked = NA,
                                   tags$b("Load the sample description file")
                          ),
                          br(),

                          tags$div(title="Load the sample description file",
                                           fileInput('ClinD',
                                                     'ClinicalData.txt',
                                                     accept = c('text/csv',
                                                                'text/comma-separated-values,text/plain',
                                                                '.csv',
                                                                '.tsv'))
                                   #shinyFilesButton(id = 'demo_clin_data', 
                                  #                  label='Select demo clinical data', 
                                  #                  title = 'Select demo clinical data', 
                                  #                  multiple=FALSE)
                                   
                          ),
                          br(),
                          shinyjs::disabled(shiny::actionButton("analyze","Analyze",class = "btn-primary"))
                          ),
                        
                        shiny::mainPanel(
                          shiny::tabsetPanel(id= "QCtab",type = "tabs",
                                              shiny::tabPanel(title = "PCA",
                                                      shiny::selectizeInput("PCs",
                                                                      label = "Choose the principal components for plotting",
                                                                      choices = list("PC1" = "PC1", "PC2" = "PC2", "PC3" = "PC3", "PC4" = "PC4", "PC5" = "PC5", "PC6" = "PC6", "PC7" = "PC7", "PC8" = "PC8"),
                                                                      selected = c("PC1", "PC1"), 
                                                                      multiple = TRUE, 
                                                                      options = list(maxItems = 2)
                                                       ),
                                                       shiny::uiOutput("labelCol"),
                                                       shiny::checkboxInput("imputeforPCA", "Use imputed data for PCA", TRUE),
                                                       shiny::plotOutput("pca_input_samples", height = 800),
                                                      br(),
                                                       shiny::downloadButton('downloadpca', 'Save'), 
                                                      br(),
                                                       downloadObjUI(id = "pca")  
                                              ),

                                             shiny::tabPanel(title = "Distribution overview",
                                                             configurePlotColUI(id = "distr"),
                                                             shiny::plotOutput("distributionPlot",height = 800) ,
                                                             br(),
                                                             shiny::downloadButton('downloadDistributionPlot', 'Save'),
                                                             br(),
                                                             downloadObjUI(id = "distr")  
                                             ),
                                              shiny::tabPanel(title = "Protein coverage",
                                                              configurePlotColUI(id = "coverage"),
                                                       shiny::plotOutput("numbers", height = 800),
                                                       br(),
                                                       shiny::downloadButton('downloadNumbers', 'Save'),
                                                       br(),
                                                       downloadObjUI(id = "coverage")  
                                              ),
                                              shiny::tabPanel(title = "Sample to sample heatmap",
                                                       shiny::selectizeInput(inputId = "distanceMetric",
                                                                      label = "Select the (dis-)similarity metric",
                                                                      choices = list("Pearson" = "Pearson", "Euclidean" = "Euclidean")
                                                       ),
                                                       shiny::plotOutput("StS_heatmap", height = 600),
                                                       br(),
                                                       shiny::downloadButton('downloadStS_heatmap', 'Save'),
                                                       br(),
                                                       downloadObjUI(id = "sts_hm")  
                                              ),
                                             shiny::tabPanel(title = "Missing value densities",
                                                             shiny::plotOutput("misVal", height = 600),
                                                             br(),
                                                             shiny::downloadButton('downloadmisVal', 'Save'),
                                                             br(),
                                                             downloadObjUI(id = "misVal")  
                                             ),
                                              shiny::tabPanel(title = "Cumulative protein intensities",
                                                       shiny::plotOutput("CumSumPlot", height = 600),
                                                       br(),
                                                       DT::dataTableOutput('cumsum_table'),
                                                       br(),
                                                       shiny::downloadButton('downloadCumSumPlot', 'Save'), 
                                                       br(),
                                                       downloadObjUI(id = "cumsum_int")  
                                              )
                        ),
                        br(),
                        shiny::uiOutput("Load_info")
                        )
                      ),            
                      hr(),
                      div(class = "footer", style="text-align:center", "Impressum",
                          shiny::includeHTML(paste(homeDir, "/../Impressum.html", sep = ""))
                      )
             ),
             shiny::tabPanel(title = "Differential Abundance",
                      shiny::sidebarLayout(       
                        shiny::sidebarPanel(
                          shiny::uiOutput("conditional_grouping_limma"),
                          shiny::uiOutput("conditional_subselectGR_limma"),
                          shiny::checkboxInput("ContinChoice", "Use continuous response instead of grouping", FALSE),
                          shiny::checkboxInput("imputeForLimma", "Impute missing values", FALSE),
                          #checkboxInput("remove_sv", "Remove surrogate variables", FALSE),
                          shiny::checkboxInput("includeCovariates", "Include parameters as covariates", FALSE), 
                          shiny::conditionalPanel("input.includeCovariates == true",
                                           shiny::uiOutput("covariatesChoice")),
                          shiny::checkboxInput("expandFilter", "Stratification and filter",  FALSE),
                          shiny::conditionalPanel("input.expandFilter == true",
                                           shiny::uiOutput("filter_group_limma")),
                          shiny::conditionalPanel("input.expandFilter == true",
                                                  shiny::uiOutput("categorizeSecNum")),
                          shiny::conditionalPanel("input.expandFilter == true",
                                           shiny::uiOutput("filter_level_limma")),
                          shiny::conditionalPanel("input.expandFilter == true",
                                           shiny::uiOutput("selectContrast")),
                          shiny::actionButton("analyzeLimma", "Analyze",class = "btn-primary"),
                          shiny::textOutput('analyzeAlerts')
                          
                        ),
                        
                        shiny::mainPanel(
                          shiny::fluidRow(
                            column(8,
                                   plotly::plotlyOutput("limma"
                                                ,height = 500
                                   ) %>% shinycssloaders::withSpinner()
                            ),
                            
                            shiny::absolutePanel(
                              top=5,
                              right = 20,
                              draggable = TRUE,
                              shiny::wellPanel(
                                #sidebarPanel(position = "right",
                                shiny::numericInput("adj.P.Val", "Adjusted P value threshold", 0.05, min = 0, max = 1, step= 0.01),
                                shiny::sliderInput("logFC", "Log Fold Change", 0, min = 0, max = 10, step = 0.1)
                              )
                            ),
                            br(),
                            column(12,
                                   DT::dataTableOutput('up'),
                                   DT::dataTableOutput('down'),
                                   br(),
                                   shiny::uiOutput("labelColBox"),
                                   shiny::checkboxInput("showLabels", "Blend in PatientID", FALSE),
                                   shiny::plotOutput("boxPlotUp"),
                                   shiny::plotOutput("boxPlotDown"),
                                   shiny::htmlOutput("doc1")
                            )),
                          
                          br(),
                          shiny::downloadButton("report", "Generate report"),
                          shiny::downloadButton("reportDataDL", "Download report data"),
                          br()
                        )
                      ),            
                      hr(),
                      div(class = "footer", style="text-align:center", "Impressum",
                          shiny::includeHTML(paste(homeDir, "/../Impressum.html", sep = ""))
                      )
                      
             ),
             shiny::tabPanel("ssGSEA", 
                      shiny::sidebarLayout(
                        shiny::sidebarPanel(
                          tags$strong("Change the parameters & hit the analyze button  ", style="color:#18bc9c"),
                          shiny::selectInput(
                            inputId = "gs.collection", 
                            label = strong("Gene Set Collection"),
                            choices = names(gene.set.databases)
                          ),
                          
                          shiny::selectInput("sample.norm.type",
                                      label = "Select a normalization method",
                                      choices = c("rank", "none"),
                                      selected = 1  
                          ),
                          shiny::numericInput("weight",
                                       label = "Select a Weight (0 to 1)",
                                       value = 0.75,
                                       min = 0, max = 1
                          ),
                          shiny::selectInput("statistic",
                                      label = "Select test statistic",
                                      choices = c("area.under.RES", "Kolmogorov-Smirnov"),
                                      selected = 1
                          ),
                          shiny::selectInput("output.score.type",
                                      label = "Select enrichment score type",
                                      choices = c("ES", "NES"),
                                      selected = 2
                          ),
                          shiny::numericInput("nperm", 
                                       label = "Enter the Number of Permutations",
                                       value = 1000
                          ),
                          shiny::numericInput("min.overlap",
                                       label = "Select the minimum overlap between gene set and data",
                                       value = 10
                          ),
                          shiny::selectInput("correl.type",
                                      label = "Select correlation type", 
                                      choices = c("rank", "z.score", "symm.rank"),
                                      selected = 1
                          ),
                          shiny::uiOutput("output.prefix"),
                          
                          tags$div(title= "Specify the type of analysis from above, then press the analyze button",
                                   actionButton("goButton", "Analyze",class = "btn-primary")         
                          )
                        ),
                        shiny::mainPanel(
                          shiny::helpText("Make sure to upload proteinGroups.txt file before running ssGSEA"),
                          shinyBS::bsCollapsePanel(p("Detailed description",style = "color:#18bc9c"),
                                          shiny::HTML(markdown::markdownToHTML(fragment.only=TRUE, text=c( 
                                            "* Single-sample GSEA [(ssGSEA)](http://software.broadinstitute.org/cancer/software/genepattern/modules/docs/ssGSEAProjection/4) is an extension of conventional  Gene Set Enrichment Analysis (GSEA),
                                        developed by Broad insitute<sup>1</sup>.",
                                            "* ssGSEA version used: v4",
                                            "* MSigDB version used: v6.1 ",
                                            "\n[1] Krug, K., et al., A Curated Resource for Phosphosite-specific Signature Analysis. Mol Cell
                                                Proteomics, 2019. 18(3): p. 576-593."
                                          )))),
                          #      ssGSEA to calculate separate pathway enrichment scores for each pairing of a sample and geneset.
                          #      Each ssGSEA enrichment score represents the degree to which the genes in a particular gene set are coordinately up- or down-regulated within a sample.")),
                          
                          shinyalert::useShinyalert()
                        )
                      ),            
                      hr(),
                      div(class = "footer", style="text-align:center", "Impressum",
                          shiny::includeHTML(paste(homeDir, "/../Impressum.html", sep = ""))
                      )
             ), 
             
             shiny::tabPanel("Differential Enrichment",
                      #uiOutput("diff.gs.collection"),
                      expDesignModule_UI(id = "gsea"),            
                      hr(),
                      div(class = "footer", style="text-align:center", "Impressum",
                          shiny::includeHTML(paste(homeDir, "/../Impressum.html", sep = ""))
                      )
                 ),
            shiny::tabPanel("Help",icon = icon("info-circle"),
                      shiny::includeHTML(paste(homeDir, "/../Vignette/Test2.html", sep = "")),            
                      hr(),
                      div(class = "footer", style="text-align:center", "Impressum",
                          shiny::includeHTML(paste(homeDir, "/../Impressum.html", sep = ""))
                      )
             )
  )
)

## Server function
server <- function(input, output, session) {
  reportBlocks <- shiny::reactiveValues()
  session$onSessionEnded(function() { unlink(paste(sessionID, "/*.gct", sep = "")) })
  #push the filesize upload limit
  options(shiny.maxRequestSize = 10000*1024^2, expressions = 500000)
  
  protfile <- shiny::reactiveValues()
  
  ###1 Load n Prep tab  
  
  #read .csv uploaded file  
  volumes <-c(root = paste(homeDir, '/../Data', sep = ""))
  shinyFiles::shinyFileChoose(input, 'files', root=c(root='./../Data'), filetypes=c('', 'txt')) 
  shinyFiles::shinyFileChoose(input, 'demo_clin_data', root=c(root='./../Data'), filetypes=c('', 'txt')) 
  
  data <- shiny::reactive({
   # if (input$dataUpload == 'userFile'){
       if (is.null(input$file1)) {
         return(NULL)
       } else{
         inFile <- input$file1 
       }
    #} else if (input$dataUpload == 'serverFile') {
    #if (is.null(input$files)){
    #   return(NULL)
    # } else {
    #   inFile <-parseFilePaths(volumes, input$files)
     #}

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
    stats_proteinGroups$IntensityMetric = input$insty
    stats_proteinGroups$ImputationMethod = input$imputation
    reportBlocks$stats_proteinGroups = stats_proteinGroups
    
    message("Loaded ", stats_proteinGroups$NumFullProt, " rows from the provided proteinGroups-file. ", 
            stats_proteinGroups$NumPotCon, " potential contaminants, ", stats_proteinGroups$IdentifiedBySite, 
            " proteins only identified by site, and ", stats_proteinGroups$Reverse, " 
            proteins only identified in the reverse database were removed from the data set. 
            The final dataset has ", stats_proteinGroups$Cleaned, " rows. " )
    proteinGroups
  })
  
  #Select columns with specific intensity prefix (LFQ/iBAQ/..)
  insty <- shiny::reactive({
    input$analyze
    data <- data()
    reportBlocks$linesFromMaxQuant = nrow(data)
      shiny::isolate(selectProteinData(data, intensityMetric = input$insty))
    
  })
  
  output$Load_info = shiny::renderUI({
    input$analyze
    validate(need(reportBlocks$stats_proteinGroups, FALSE))
    shinyBS::bsCollapsePanel(p("Detailed description",style = "color:#18bc9c"),
                             shiny::HTML(markdown::markdownToHTML(fragment.only=TRUE, text = paste(
                               "Loaded ", reportBlocks$stats_proteinGroups$NumFullProt, " rows from the provided proteinGroups-file. ",
                               reportBlocks$stats_proteinGroups$NumPotCon, " potential contaminants, ", reportBlocks$stats_proteinGroups$IdentifiedBySite,
                               " proteins only identified by site, and ", reportBlocks$stats_proteinGroups$Reverse, "
                               proteins only identified in the reverse database were removed from the data set.
                               The final dataset has ", reportBlocks$stats_proteinGroups$Cleaned, " rows. "
                             ))))
  })
  
  #Remove user defined columns/samples
  ## Filter samples UI elements
  output$filt = shiny::renderUI({
    if (!is.null(data())){
    req(data())
      insty<-insty()
      shiny::selectizeInput("filt",
                     "Exclude columns (samples)",
                     choices=colnames(insty),
                     multiple = TRUE,
                     selected=input$filt)
    } else {
      shiny::selectizeInput("filt",
                     "Exclude columns (samples)",
                     choices=colnames(insty))
    }
  })
  
  ## Filter user defined columns/samples
  filt <- shiny::reactive({
    input$analyze
    filt_col<-shiny::isolate(insty())
    if(is.null(input$filt)){
      filt <- shiny::isolate(filt_col)
    } else {
      filt <- shiny::isolate(filt_col[ , -which(names(filt_col) %in% input$filt)])
    }
  })
  
  # Filter proteins that were not detected in at least a user defined amount of samples
  filtpro <- shiny::reactive({
    input$analyze
    filt_pro <- shiny::isolate(filt())
    filterTH <- input$filter_pro-1
    filtpro <- shiny::isolate(filterProteins(filt_pro, filterTH = filterTH))
  })
  
  #### new universal format for the protein intensity data that should be feeded to any other requiring method is specified as
  # being a tibble
  # containing log2 transformed, normalized intensity values as specified by the user 
  # one column "Gene names" that harbours uniqe gene names and no NA values
  # the new format should be accepted by any method and also returned by any further method
  # the new format should be saved in a reactiveValues object, one containing the imputed and one the not imputed values
  
  proteinAbundance <- shiny::reactiveValues()
  
  #check for isoforms or make unique gene names and log2 transformation
  # needed parameters: 
  # input$unique
  # reactiveValue on log2-transform TRUE/FALSE
  # reactiveValue on normalization (only needed when iBAQ values are used)
  log2tansform <- reactiveVal(TRUE)
  normalizeVsn <- reactiveVal(FALSE)
  
  uniq_names <- shiny::observeEvent(c(
    input$analyze, 
    input$insty,
    input$unique, 
    input$filterTH,
    input$filt,
    normalizeVsn(),
    log2tansform()), ignoreNULL = TRUE,  ignoreInit = TRUE, 
    {
      shiny::req(protfile$protfile$name)
      filt<-filtpro()
      if (shiny::isolate(input$unique == "chck_iso" )) {
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
          data_unique = limma::normalizeVSN(data_unique[,LFQ_columns])
        }
      } else if (log2tansform() == FALSE) {
        data_unique = data_unique[,LFQ_columns]
      }
      else {
        data_unique = data_unique[,LFQ_columns]
      }
      proteinAbundance$original = dplyr::bind_cols(rn, data_unique)
    })
  
  impute_uniques <- shiny::observeEvent(c(proteinAbundance$original,
                                   input$imputation), ignoreNULL = TRUE, ignoreInit = TRUE,{
                                     
                                     imputed = proteinAbundance$original %>% tibble::column_to_rownames("Gene names") %>% as.data.frame()
                                     if (input$imputation == "perseus-like"){
                                       imputed[[1]] = replaceMissingFromGaussian(imputed)
                                     } else if (input$imputation == "QRILC"){
                                       imputed <-  imputeLCMD::impute.QRILC(as.matrix(imputed))
                                     } else if (input$imputation == "MinDet"){
                                       imputed[[1]]<-  imputeLCMD::impute.MinDet(as.matrix(imputed))
                                     } else if (input$imputation == "knn"){
                                       imputed <- impute::impute.knn(as.matrix(imputed))
                                     }
                                     proteinAbundance$imputed = imputed[[1]] %>% tibble::as_tibble(., rownames = "Gene names")
                                   })
  
  
  
  
  #filter and plot the results 
  QCreport<- shiny::reactiveValues()
  
  
  shiny::observeEvent({
    input$ClinD
  }, { 
    
    output$labelCol<- shiny::renderUI({
      shiny::conditionalPanel(condition = need(ClinData(), FALSE) , 
                       shiny::selectInput(
                         inputId = "labelCol", 
                         label = strong("Choose the clinical parameter for group colours"),
                         choices = as.list("none" = "none", colnames(ClinData())),
                         multiple = FALSE,
                         selectize = TRUE
                       )
      )
    })
    
    
    
    pca_input_samples <- shiny::reactive({
      shiny::validate(need(data(), "Please upload a proteinGroups.txt file first."))
      log2tansform(TRUE)
      if(input$insty == "iBAQ") {normalizeVsn(TRUE)}
      if (input$imputeforPCA == TRUE){
        check_table <- proteinAbundance$imputed %>% tibble::column_to_rownames("Gene names") %>% as.data.frame()
      } else {
        check_table <- proteinAbundance$original %>% tibble::column_to_rownames("Gene names") %>% as.data.frame()
      }
      check_table = t(as.matrix(check_table))
      check_table <- check_table[ , which(apply(check_table, 2, var) != 0)]
      cc<-stats::prcomp(na.omit(check_table),
                 scale. =TRUE, 
                 center = TRUE
      )
    })
    
    ## Remove? 
    PCs <- shiny::reactiveVal(input$PCs)
    
    output$pca_input_samples <- shiny::renderPlot({
      shiny::validate(
        need(length(input$PCs) == 2, "Please select two principal components for plotting.")
      )
      cc <- pca_input_samples()
      dummy = as.data.frame(cc$x) %>% rownames_to_column("PatientID")
      if (is.null(input$labelCol) | input$labelCol == "" | input$labelCol=="none" ) {
        groups = NULL
        labels = "PatientID"
        
      } else if (!is.null(input$labelCol) && input$labelCol != "PatientID") {
        dummy = dummy %>% left_join(., dplyr::select(ClinData(), PatientID, groups = input$labelCol))
        labels = "PatientID"
      }  else {
        message("Please select a different parameter for colouring")
        groups = NULL
        labels = "PatientID"
      }
      ellipse = TRUE
      if (is.numeric(dummy$groups)){
        ellipse = FALSE
      }

      biplot = ggplot(dummy, aes_string(x = input$PCs[1], y = input$PCs[2])) + 
        geom_point(aes(color = groups)) 
      if (ellipse == TRUE) {
        biplot = biplot + 
          stat_ellipse(aes(group = groups, color = groups))
      }
      
      biplot <- biplot + 
        ggplot2::ggtitle("Principal component analysis") +
        ggplot2::theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        ggrepel::geom_text_repel(aes_string(label=labels), show.legend = F, size = 4) +
        #scale_fill_continuous(na.value="white") + 
        ggplot2::theme_light()
      if (is.numeric(dummy$groups)){
        biplot  = biplot + ggthemes::scale_color_continuous_tableau("Red-Gold",na.value = "grey30") 
      } else {
        biplot  = biplot + ggthemes::scale_color_tableau()
      }
      QCreport$pca = biplot
      biplot
    })

    
    output$downloadpca <- shiny::downloadHandler(
      filename = "pca.pdf",
      content = function(file) {
        plot_parameters = callModule(downloadObj, id = "pca", title = "Principal component analysis", filename = "pca.pdf")
        QCreport$pca = add_plot_info(QCreport$pca, plot_parameters)
        grDevices::pdf(file)
        print(QCreport$pca)
        grDevices::dev.off()
    })

    #Distribution plot
    
    labelCol = callModule(configurePlotCol, id = "distr", ClinData = ClinData())
    
    output$distributionPlot <- shiny::renderPlot({
      log2tansform(TRUE)
      ClinData = ClinData()
      dummy = tibble()
        if (!is.null(labelCol()) && labelCol() != "PatientID") {
        dummy = dplyr::select(ClinData, PatientID, groups = labelCol())
      }  else {
       # message("Please select a different parameter for colouring")
        dummy = dplyr::select(ClinData, PatientID) %>% bind_cols(groups = rep("constant", nrow(.)))
      }
      original = proteinAbundance$original %>% as_tibble() #%>% tibble::column_to_rownames("Gene names") #%>% as.data.frame()
      QCreport$distributionPlot = plot_distribution(original, dummy)
      QCreport$distributionPlot
    })
    
    output$downloadDistributionPlot <- shiny::downloadHandler(
      filename = "DistributionPlot.pdf",
      content = function(file) {
        plot_parameters = callModule(downloadObj, id = "distr", title = "Distribution of proetin intensities", filename = "DistributionPlot.pdf")
        QCreport$distributionPlot = add_plot_info(QCreport$distributionPlot, plot_parameters)
        grDevices::pdf(file)
        print(QCreport$distributionPlot)
        grDevices::dev.off()
      })

    
    #Protein coverage plot
    labelCol_cov = callModule(configurePlotCol, id = "coverage", ClinData = ClinData())
    
    output$numbers <- renderPlot({      
      original = proteinAbundance$original %>% as_tibble()# column_to_rownames("Gene names") %>% as.data.frame()
      ClinData = ClinData()
      dummy = tibble()
      
      if (!is.null(labelCol_cov()) && labelCol_cov() != "PatientID") {
        dummy = dplyr::select(ClinData, PatientID, groups = labelCol_cov())
      }  else {
        dummy = dplyr::select(ClinData, PatientID) %>% bind_cols(groups = rep("constant", nrow(.)))
      }
      
      QCreport$number <- plot_proteinCoverage(original, dummy)
      QCreport$number
    })
    
    output$downloadNumbers <- downloadHandler(
      filename = "ProteinNumbers.pdf",
      content = function(file) {
        plot_parameters = callModule(downloadObj, id = "coverage", title = "Measured proteins per sample", filename = "ProteinNumbers.pdf")
        QCreport$number = add_plot_info(QCreport$number, plot_parameters)
        grDevices::pdf(file)
        print(QCreport$number)
        grDevices::dev.off()
      })

    #Sample-to-Sample Heatmap
    StSheatmap_input <- shiny::reactive({
      log2tansform(TRUE)
      original = proteinAbundance$original %>% tibble::column_to_rownames("Gene names") %>% as.data.frame()
      if (input$distanceMetric == "Pearson") {
        corr = TRUE
      } else {corr = FALSE}
      
      plot_StS_heatmap(original, corr = corr)
    })
    
    output$StS_heatmap <- shiny::renderPlot({
      QCreport$StSDistMetric = input$distanceMetric
      QCreport$StSheatmap <- StSheatmap_input()
      StSheatmap_input()
    })
    output$downloadStS_heatmap <- shiny::downloadHandler(
      filename = "StS_heatmap.pdf",
      content = function(file) {
        plot_parameters = callModule(downloadObj, id = "sts_hm", title = "Sample to sample heatmap", filename = "StS_heatmap.pdf")
        QCreport$StSheatmap = add_plot_info(QCreport$StSheatmap, plot_parameters)
        grDevices::pdf(file)
        print(QCreport$StSheatmap)
        print(QCreport$StSDistMetric)
        grDevices::dev.off()
      })
    
    # Missing values density
    output$misVal = shiny::renderPlot({
      log2tansform(TRUE)
      QCreport$misVal = plot_misValDens(proteinAbundance)
      QCreport$misVal
    })
    
    output$downloadNumbers <- downloadHandler(
      filename = "MisValDensity.pdf",
      content = function(file) {
        plot_parameters = callModule(downloadObj, id = "misVal", title = "Distribution of Intensity for origanal and missing values", filename = "MisValDensity.pdf")
        QCreport$number = add_plot_info(QCreport$misVal, plot_parameters)
        grDevices::pdf(file)
        print(QCreport$misVal)
        grDevices::dev.off()
      })
    
    
    
    # Cumulative Intensities 
    output$CumSumPlot <- shiny::renderPlot({
      input$analyze
      log2tansform(FALSE)
      original = proteinAbundance$original %>% tibble::column_to_rownames("Gene names") %>% as.data.frame()
      QCreport$cumsum = CumSumIntensities(original)[1]
      QCreport$cumsum_table = CumSumIntensities(original)[[2]]

     # grid::grid.draw(CumSumPlot_input())
      QCreport$cumsum
    })
    
    output$cumsum_table <- DT::renderDataTable({
      input$analyze
      original = proteinAbundance$original %>% tibble::column_to_rownames("Gene names") %>% as.data.frame()
        df <- CumSumIntensities(original)[[2]]

       df = DT::datatable(df,
                      rownames = TRUE,
                      class = 'cell-border stripe',
                      width= '150px',
                      extensions ='Scroller', 
                      options = list(
                        deferRender = TRUE,
                        scrollY = 100,
                        scroller = TRUE,
                        scrollCollapse = TRUE,
                        pageLength = 100, 
                        lengthMenu = c(5,10,50,100,200)
                      )) 

    })
    
    output$downloadCumSumPlot <- downloadHandler(
      filename = "CumSumPlot.pdf",
      content = function(file) {
        plot_parameters = callModule(downloadObj, id = "cumsum_int", title = "Cumulative protein intensities", filename = "CumSumPlot.pdf")
        original = proteinAbundance$original %>% tibble::column_to_rownames("Gene names") %>% as.data.frame()
        QCreport$cumsum = add_plot_info(CumSumIntensities(original)[[1]], plot_parameters)
        grDevices::pdf(file)
        print(QCreport$cumsum)
        print(CumSumIntensities(original)[[2]])
        grDevices::dev.off()
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
  output$labelColBox<- shiny::renderUI({
    shiny::conditionalPanel(condition = need(ClinData(), FALSE) , 
                     shiny::selectInput(
                       inputId = "labelColBox", 
                       label = strong("Choose the clinical parameter for label colours"),
                       choices = as.list("none" = "none", colnames(ClinData())),
                       multiple = FALSE,
                       selectize = TRUE
                     )
    )
  })
  
  output$conditional_grouping_limma<- shiny::renderUI({
    shiny::selectInput(
      inputId = "GR_fatcor", 
      label = strong("Select the clinical grouping factor"),
      choices = as.list(colnames(ClinData())),
      multiple = FALSE,
      selectize = TRUE
    )
  })
  
  output$conditional_subselectGR_limma <- shiny::renderUI({
    req(ClinData(), ClinColClasses(), input$GR_fatcor)
    if (ClinColClasses()[input$GR_fatcor]=='factor' | ClinColClasses()[input$GR_fatcor]=='logical' ){
      shiny::selectizeInput(inputId = "levels",
                     label= "Select two groups to compare",
                     choices = ClinData() %>% dplyr::pull(input$GR_fatcor) %>% levels(),
                     multiple = TRUE, 
                     options = list(maxItems = 2)
      )
    } else {
      d = ClinData() %>% dplyr::pull(input$GR_fatcor)
      shiny::sliderInput(
        inputId = "num.cutoff",
        label = "Select cutoff to divde numeric value:",
        min = min(d, na.rm = TRUE),
        max = max(d, na.rm = TRUE),
        value = colMeans(ClinData()[input$GR_fatcor], na.rm = TRUE), round = T
      )
    }
  }) 
  
  shiny::observeEvent(input$includeCovariates, {
    output$covariatesChoice<- shiny::renderUI({
      shiny::selectInput(
        inputId = "covariates",
        label = strong("Select factors to include as covariates."),
        choices = as.list(colnames(ClinData())),
        multiple = TRUE
      )
    })
  })
  
  shiny::observeEvent(input$expandFilter, {
    
    output$filter_group_limma<- shiny::renderUI({
      shiny::selectInput(
        inputId = "filter_GR_fatcor",
        label = strong("Select a second parameter"),
        selected = 3,
        choices = as.list(colnames(ClinData())),
        multiple = FALSE,
        selectize = TRUE
      )
    })
    
   output$categorizeSecNum <- shiny::renderUI({
     if (ClinColClasses()[input$filter_GR_fatcor]=='numeric') {
       d = ClinData() %>% dplyr::pull(input$filter_GR_fatcor)
       shiny::sliderInput(
         inputId = "filter_num.cutoff",
         label = "Select cutoff to divide numeric value:",
         min = min(d, na.rm = TRUE),
         max = max(d, na.rm = TRUE),
         value = colMeans(ClinData()[input$filter_GR_fatcor], na.rm = TRUE), 
         round = T
       )
     } else {NULL}
   }) 

    output$filter_level_limma <- shiny::renderUI({
      req(input$filter_GR_fatcor)
      req(ClinDomit$filterParameter)
      if (#(ClinColClasses()[input$filter_GR_fatcor]=='factor' | ClinColClasses()[input$filter_GR_fatcor]=='logical') & 
        input$ContinChoice == TRUE){
     # if (ClinColClasses()[input$filter_GR_fatcor]=='factor' | ClinColClasses()[input$filter_GR_fatcor]=='logical' | ClinColClasses_2()[ClinDomit$filterParameter]=='factor'){
        shiny::selectizeInput(inputId = "filter_levels",
                       label = "Filter: Select groups to include in the analysis",
                       choices = ClinData() %>% dplyr::pull(input$filter_GR_fatcor) %>% levels(),
                       multiple = TRUE
        )
      } #else {
        #d = ClinData() %>% dplyr::pull(input$filter_GR_fatcor)
        #shiny::sliderInput(
        #  inputId = "filter_num.cutoff",
        #  label = "Select cutoff to divide numeric value:",
        #  min = min(d, na.rm = TRUE),
        #  max = max(d, na.rm = TRUE),
        #  value = colMeans(ClinData()[input$filter_GR_fatcor], na.rm = TRUE), 
        #  round = T
        #)
      #}
    })
    output$selectContrast <- shiny::renderUI({
      req(ClinDomit$mainParameter)
      if (input$ContinChoice == FALSE){
        shiny::selectizeInput(inputId = "contrastLevels", 
                              label = "Stratify: Select the two groups you want to calculate the difference on.",
                              #choices = ClinData() %>% pull(mainParameter) %>% levels()
                              choices = ClinDomit$data %>% dplyr::pull(ClinDomit$mainParameter) %>% levels(),
                              multiple = TRUE, 
                              options = list(maxItems = 2)
        )
      }

    })
  }, ignoreNULL = FALSE, ignoreInit = TRUE)
  
  shiny::observe({
    shiny::req(input$GR_fatcor)
    shiny::req(input$ContinChoice)
    input$GR_fatcor
    if(ClinColClasses()[input$GR_fatcor] != "numeric") {
      shiny::showNotification("Please make sure that you have selected a continuous variable.")
      shiny::updateCheckboxInput(session, "ContinChoice", value = FALSE)
    }
  })
  
  shiny::observe({
    shiny::req(input$GR_fatcor)
    shiny::req(input$ContinChoice)
    #shiny::need(input$ContinChoice)
    input$GR_fatcor
    shiny::updateCheckboxInput(session, "expandFilter", value = FALSE)
  })
  
  
  
  # Logic elements
  ## Load clinical parameters
  
  volumes2 <-c(root = paste(homeDir, '/../Data', sep = ""))
  shinyFiles::shinyFileChoose(input, 'ClinDs', root=c(root='./../Data'), filetypes=c('', 'txt', 'tsv'))
  
  ClinData <- shiny::reactive({
    
    if (is.null(input$ClinD)) {
      return(NULL)
    } else{
      clinfile$name <- input$ClinD
    }

      #shiny::validate(
      #  need(input$ClinD != "", "Please provide a sample description file for upload on the previous tab.")
        #need(input$demo_clin_data != "", "Please select a file for upload or choose to use the database connection.")
    #  )
      #clinfile$name <- input$ClinD
    #clinfile$name <-parseFilePaths(volumes2, input$demo_clin_data)

    if (length(clinfile$name$datapath) == "0"){
      disable("analyze")
      return(NULL)
    }

      
    ClinData = readr::read_tsv(clinfile$name$datapath, 
    #ClinData = readr::read_tsv(input$demo_clin_data$datapath, 
                        na =c("", "NA", "N/A","0","<Null>"),
    # locale = locale(decimal_mark = ","), # Todo: uncomment this for US/English seperator
                        skip_empty_rows = TRUE 

    ) %>% janitor::remove_empty("cols") %>% dplyr::mutate_if(is.character, as.factor) %>% dplyr::mutate_if(is.logical, as.factor)
  #   if(protfile$protfile$name == "proteinGroups.freeze.txt"){ ### TODO: Uncomment for public use!
  #     ClinData = ClinData %>% dplyr::rename(PatientID_DB = PatientID)
  #     ClinData = ClinData %>% dplyr::rename(PatientID = freezeID_log)
   #  }
    ClinDomit$data = ClinData %>% janitor::clean_names()
    shinyjs::enable("analyze") 
    ClinData
  })
  
  ClinColClasses <- shiny::reactive({
    df = ClinData()
    #df = ClinDomit$data
    df = lapply(df, class)
  })
  
  ClinColClasses_2 <- shiny::reactive({
    df = ClinDomit$data
    df = lapply(df, class)
  })
  
  
  
  # set and manipulate chosen parameter to being cat from cont and form new groups from stratified setups
  shiny::observeEvent(c(input$filter_GR_fatcor, 
                 input$GR_fatcor,
                 input$ContinChoice, 
                 input$num.cutoff,
                 input$filter_num.cutoff
                 )
               , {
                 limmaResult$gene_list = NULL #refresh volcano plot 
                 ClinDomit$filterParameter = NULL
                 ClinDomit$mainParameter = NULL
      ClinDomit$mainParameter = janitor::make_clean_names(input$GR_fatcor)

      ## categorize numeric data - first parameter
      if (input$ContinChoice == FALSE & ClinColClasses_2()[ClinDomit$mainParameter]=='numeric') {
        shiny::req(input$num.cutoff)
        ClinDomit$data = ClinDomit$data %>% 
          dplyr::mutate(categorizedParameter = 
                   cut(dplyr::pull(ClinDomit$data, ClinDomit$mainParameter), 
                       breaks = c(-Inf, input$num.cutoff, Inf), 
                       labels = c(paste('less than or equal to', input$num.cutoff, sep='_'), paste('greater than', input$num.cutoff,sep='_'))
                   ))  %>% 
          dplyr::mutate_if(is.character, as.factor)
        #%>% 
          #dplyr::select(-c(!!mainParameter)) 
        colnames(ClinDomit$data)[colnames(ClinDomit$data) == "categorizedParameter"] = paste(input$GR_fatcor, "cat", sep = "_", collapse = "_") %>% janitor::make_clean_names()
        ClinDomit$data = ClinDomit$data[,!duplicated(colnames(ClinDomit$data), fromLast = TRUE)]
        ClinDomit$mainParameter = paste(input$GR_fatcor,  "cat", sep = "_", collapse = "_") %>% janitor::make_clean_names()
      }
      if (is.null(shiny::need(input$expandFilter, FALSE)) & is.null(shiny::need(input$filter_GR_fatcor, FALSE))) {
        
        ## categorize second numeric parameter
        ClinDomit$filterParameter =  janitor::make_clean_names(input$filter_GR_fatcor)
        if (ClinColClasses()[input$filter_GR_fatcor]=='numeric') {
          req(input$filter_num.cutoff)
          ClinDomit$data = ClinDomit$data %>% 
            dplyr::mutate(filterParameter = 
                     cut(dplyr::pull(ClinDomit$data, ClinDomit$filterParameter), 
                         breaks = c(-Inf, input$filter_num.cutoff, Inf), 
                         labels = c(paste('less than or equal to', input$filter_num.cutoff, sep='_'), paste('greater than', input$filter_num.cutoff,sep='_'))
                     )) %>% 
            dplyr::mutate_if(is.character, as.factor)
          # if first parameter stays continuous, the second parameter becomes a filter and needs to be saved for experimental design creation 
          if (ClinColClasses_2()[ClinDomit$mainParameter]=='numeric' & input$ContinChoice == TRUE){
            ## rename and save filter parameter
            colnames(ClinDomit$data)[colnames(ClinDomit$data) == "filterParameter"] = paste(ClinDomit$filterParameter, "cat", sep = "_", collapse = "_") %>% janitor::make_clean_names() 
            ClinDomit$data = ClinDomit$data[,!duplicated(colnames(ClinDomit$data), fromLast = TRUE)]
            ClinDomit$filterParameter = paste(ClinDomit$filterParameter, "cat", sep = "_", collapse = "_") %>% janitor::make_clean_names() 
           # ClinDomit$filterParameter = input$filter_GR_fatcor
          } else {## unite cat first parameter and categorized second parameter, when first is cat or categorized
            ClinDomit$data = ClinDomit$data %>% 
              tidyr::unite("newFactor", ClinDomit$mainParameter, filterParameter, remove = FALSE) %>% 
              dplyr::mutate_if(is.character, as.factor)
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
            ClinDomit$data = ClinDomit$data %>% tidyr::unite("newFactor", ClinDomit$mainParameter, ClinDomit$filterParameter, remove = FALSE)  %>% dplyr::mutate_if(is.character, as.factor)
            colnames(ClinDomit$data)[colnames(ClinDomit$data) == "newFactor"] = paste(ClinDomit$mainParameter, ClinDomit$filterParameter, sep = "_", collapse = "_")
            ClinDomit$data = ClinDomit$data[,!duplicated(colnames(ClinDomit$data), fromLast = TRUE)]
            ClinDomit$mainParameter = paste(ClinDomit$mainParameter, ClinDomit$filterParameter, sep = "_", collapse = "_")
            
          }
          
        }
      }
    }, ignoreInit = TRUE
#    , once = TRUE
    )
  
  analyzeAlerts <- reactiveValues("somelist" = c(FALSE, NULL))
  
  output$analyzeAlerts <- renderText({
    input$analyzeLimma
    if(analyzeAlerts$somelist[1] == TRUE) {
      analyzeAlerts$somelist[2]
    } else {""}
  })

  observeEvent(input$analyzeLimma ,{
    limmaResult$gene_list = NULL #refresh volcano plot
    reportBlocks$volcano_plot = NULL 
    
   if (!is.null(need(proteinAbundance$original, "TRUE"))) { # check if protein abundance is loaded
      analyzeAlerts$somelist = c(TRUE, "Please upload a proteinGroups file first (previous tab).")  
    } else {analyzeAlerts$somelist = c(FALSE, NULL)}
    shiny::validate(need(proteinAbundance$original , "Validate statement"))
    
    mainParameter = janitor::make_clean_names(ClinDomit$mainParameter)
    if (!is.null(ClinDomit$filterParameter)) {
      filterParameter = janitor::make_clean_names(ClinDomit$filterParameter)
    }
    #    }
    if (input$expandFilter == TRUE & input$ContinChoice == FALSE) {
      shiny::req(input$contrastLevels)
    }
    if ((!is.null(need(input$levels, "TRUE")) | length(input$levels) != 2) & ClinColClasses()[input$GR_fatcor]!='numeric') {
      analyzeAlerts$somelist = c(TRUE, "Please select two groups to compare.")  
      shiny::validate(need(input$levels, "Validate statement"))
      shiny::validate(need(length(input$levels) == 2, "Validate statement." ))
    }
    else {analyzeAlerts$somelist = c(FALSE, NULL)}
    
    if (is.null(input$contrastLevels)){
      contrastLevels = input$levels
    } else{
      contrastLevels = input$contrastLevels
    }
    reportBlocks$contrastLevels = contrastLevels
    
    covariates = input$covariates
    reportBlocks$covariates = covariates
    covariates = janitor::make_clean_names(covariates)

    
    if (!is.null(ClinDomit$filterParameter) & ClinColClasses_2()[mainParameter] == "numeric") {
      ClinData = ClinDomit$data %>% 
        janitor::clean_names() %>% 
        dplyr::select(patient_id, mainParameter, covariates, !!sym(filterParameter)) %>% 
        ggplot2::remove_missing() %>% 
        dplyr::filter(!!sym(filterParameter) %in% !!input$filter_levels) %>% 
        dplyr::select(-filterParameter)
    } else {
      ClinData = ClinDomit$data %>% 
        janitor::clean_names() %>% 
        dplyr::select(patient_id, mainParameter, covariates) %>% 
        ggplot2::remove_missing() 
    }
    ClinColClasses = lapply(ClinData, class)

    if (input$ContinChoice == FALSE & ClinColClasses[ClinDomit$mainParameter]=='numeric')
      {shiny::req(input$num.cutoff)}
    
    # Prep experimental design for first parameter being cat
    if (ClinColClasses[mainParameter]=='factor' | ClinColClasses[mainParameter]=='logical' ){
      #Reorder levels to mirror users selection of contrasts:
      if (ClinColClasses[mainParameter]=='factor') {
        ClinData[,mainParameter] = forcats::fct_relevel(pull(ClinData, mainParameter), contrastLevels)
      }
      if (length(covariates) == 0){
        expDesign = modelr::model_matrix(ClinData, as.formula(paste("~0", mainParameter, sep = "+", collapse = "+")))       
      } else {
        expDesign = modelr::model_matrix(ClinData, as.formula(paste("~0", mainParameter, paste(covariates, sep = "+", collapse = "+"), sep = "+", collapse = "+")))       
      }
    }
    # Prep experimental design for first parameter being cont.     
    if (input$ContinChoice == TRUE){
      if (length(covariates) == 0){
        expDesign = modelr::model_matrix(ClinData, as.formula(paste("~ 1", mainParameter, sep = "+", collapse = "+")))       
      } else {
        expDesign = modelr::model_matrix(ClinData, as.formula(paste("~ 1", mainParameter, paste(covariates, sep = "+", collapse = "+"), sep = "+", collapse = "+")))       
      }
    }
    
    expDesign = data.frame(expDesign, check.names = FALSE)
    rownames(expDesign) <- ClinData$patient_id
    expDesign = matchedExpDesign(expDesign, proteinAbundance$original)
    ClinDomit$designMatrix = expDesign
  
    #even <- input$n %% 2 == 0
    #shinyFeedback::feedbackWarning("n", !even, "Please select an even number")
    
    shiny::req(ClinDomit$designMatrix)

    if (!is.null(need(sum(ClinDomit$designMatrix[,1])>=3, "TRUE"))) {
      analyzeAlerts$somelist = c(TRUE, "The experimental design does not contain three or more samples to test on.")  
    }
    shiny::validate(need(sum(ClinDomit$designMatrix[,1])>=3, "Validate statement"))

    expDesignInst = ClinDomit$designMatrix %>% janitor::clean_names()
    if (input$ContinChoice == FALSE){
      if (!is.null(need(sum(ClinDomit$designMatrix[,2])>=3, "TRUE"))) {
        analyzeAlerts$somelist = c(TRUE, "The experimental design does not contain three or more samples to test on.")  
      }
      shiny::validate(need(sum(ClinDomit$designMatrix[,2])>=3, "Validate statement"))
      
      validProteins =  proteinAbundance$original[, rownames(expDesignInst)]
      validProteins_1 = validProteins[,expDesignInst[,1]]
      validProteins_2 = validProteins[,expDesignInst[,2]]
      validProteinsLog = (rowSums(is.na(validProteins_1)) < length(validProteins_1)*0.5) + (rowSums(is.na(validProteins_2)) < length(validProteins_2)*0.5)
      
      ## Prepare protein abundance for limma
      #match_Protdata =  match_Protdata[, rownames(designMatrix)]
      if (input$imputeForLimma == TRUE){
        proteinAbundanceLimma <- proteinAbundance$imputed[as.vector(validProteinsLog > 0), c("Gene names", rownames(expDesignInst))] %>% tibble::as_tibble()
      } else {
        proteinAbundanceLimma <- proteinAbundance$original[as.vector(validProteinsLog > 0),c("Gene names", rownames(expDesignInst))] %>% tibble::as_tibble()
      }
      #proteinAbundanceLimma <- proteinAbundance$imputed[as.vector(validProteinsLog > 0),]
      proteinAbundanceLimma <- proteinAbundanceLimma %>% tibble::column_to_rownames("Gene names")
      proteinAbundanceLimma <- proteinAbundanceLimma[which(apply(proteinAbundanceLimma, 1, var, na.rm = TRUE) != 0), ]
      
      # no intercept as no reference level is assumed in this setup
      fit <- limma::lmFit(proteinAbundanceLimma, as.matrix(expDesignInst))
      sv_contrasts <- limma::makeContrasts(contrasts=paste(colnames(expDesignInst)[1],colnames(expDesignInst)[2], sep="-"), levels= expDesignInst)
      fit <- limma::contrasts.fit(fit, contrasts=sv_contrasts)
    } else {
      validProteins =  proteinAbundance$original[, rownames(expDesignInst)]
      validProteinsLog = (rowSums(!is.na(validProteins)) >= 5)
      if (input$imputeForLimma == TRUE){
        proteinAbundanceLimma <- proteinAbundance$imputed[as.vector(validProteinsLog > 0), c("Gene names", rownames(expDesignInst))] %>% tibble::as_tibble()
      } else {
        proteinAbundanceLimma <- proteinAbundance$original[as.vector(validProteinsLog > 0),c("Gene names", rownames(expDesignInst))] %>% tibble::as_tibble()
      }
     # proteinAbundanceLimma <- proteinAbundance$imputed[as.vector(validProteinsLog > 0), c("Gene names", rownames(expDesignInst))]
      proteinAbundanceLimma <- proteinAbundanceLimma %>% tibble::column_to_rownames("Gene names")
      proteinAbundanceLimma <- proteinAbundanceLimma[which(apply(proteinAbundanceLimma, 1, var, na.rm = TRUE) != 0), ]
      
      fit <- limma::lmFit(proteinAbundanceLimma, as.matrix(expDesignInst))
    }
    fit3 <- limma::eBayes(fit, trend = TRUE) 
    if (input$ContinChoice == FALSE) {
      limmaResult$gene_list <- limma::topTable(fit3, coef=1, number = 1e+09, sort.by="P", adjust.method="BH")
    } else if (input$ContinChoice == TRUE)  {
      limmaResult$gene_list <- limma::topTable(fit3, coef=2, number = 1e+09, sort.by="P", adjust.method="BH")
    }
  }, ignoreInit = TRUE)
  
  
  #limma_input <- reactive({
  limma_input <- shiny::eventReactive(c(
    limmaResult$gene_list,
    input$analyzeLimma ,
    input$adj.P.Val,
    input$logFC
  ),{
    reportBlocks$volcano_plot = NULL
    shiny::req(limmaResult$gene_list)
    gene_list <- limmaResult$gene_list
    message("Genes in Limma: ", nrow(gene_list))
    gene_list$threshold = as.factor(abs(gene_list$logFC) > input$logFC & gene_list$adj.P.Val < input$adj.P.Val)
    
    reportBlocks$volcano_plot = 
      ggplot2::ggplot(data=gene_list,aes(x=logFC, y=-log10(P.Value) , colour=threshold)) +
      ggtitle("this should be a title") +
    #  ggtitle(paste(ClinDomit$mainParameter,": ", names(ClinDomit$designMatrix)[2]," vs ", names(ClinDomit$designMatrix[1]), collapse = "")) +
       ggplot2::geom_point(shape=20, data = gene_list, aes(text= rownames(gene_list)), size = 1, alpha = 0.4) +
      #labs(color = paste("Threshold \n adj.p < ", input$adj.P.Val, " and \n log2FC +/- ", input$logFC)) +
      ggplot2::theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
      ggthemes::scale_color_tableau() +
      ggplot2::theme_light() +
      #ggplot2::theme(legend.title = element_blank())
      guides(color=guide_legend(title="Significant"))

  })
  
  
  output$limma <- plotly::renderPlotly({
    shiny::req(limmaResult$gene_list)
    input$adj.P.Val
    input$logFC 
    input$analyzeLimma 
    gene_list <-limmaResult$gene_list
    
    ##select upregulated dataframe
    
    limmaResult$upregulated = as.factor(gene_list$logFC >= input$logFC & gene_list$adj.P.Val <= input$adj.P.Val)
    limmaResult$up <- na.omit(gene_list[limmaResult$upregulated=="TRUE",])
    limmaResult$df_up <- limmaResult$up %>% 
      tibble::rownames_to_column("Gene names") %>% 
      dplyr::left_join(., reportBlocks$ProteinIDMap, by = ("Gene.name" = "Gene names"))

    
    ##select down regulated dataframe
    
    limmaResult$downregulated= as.factor(gene_list$logFC <= (-input$logFC) & gene_list$adj.P.Val <= input$adj.P.Val)
    limmaResult$down <- na.omit(gene_list[ limmaResult$downregulated =="TRUE",])
    limmaResult$df_down<-limmaResult$down %>% 
      tibble::rownames_to_column("Gene names") %>% 
      dplyr::left_join(., reportBlocks$ProteinIDMap, by = ("Gene.name" = "Gene names"))
    
    temp_volcano = limma_input()
    #reportBlocks$volcano_plot <- limma_input()
    
    #Prepare interactive plot, reactive title and legend
   

    if (!is.null(input$expandFilter) && input$expandFilter == TRUE) {
      Filnames = paste(" only ", input$filter_levels[1], collapse= "")
    } else {
      Filnames = ""
    }
    
    input$analyzeLimma 
    shiny::isolate(
    if (input$ContinChoice == FALSE){
      title_begin = paste("Proteins with a change in abundance in ", names(ClinDomit$designMatrix)[1], 
                                           " when compared to \n", 
                                           names(ClinDomit$designMatrix[2]), 
                                           Filnames, collapse = "")
    } else {
      title_begin =  paste("Proteins abundance with regard to ", names(ClinDomit$designMatrix)[2], Filnames, collapse = "")
    }
    )
    reportBlocks$volcano_plot = temp_volcano + ggtitle(title_begin)

    pp <- plotly::ggplotly(temp_volcano, tooltip = "text") %>% 
      plotly::layout(title = paste0('Volcano plot',
                                      '<br>',
                                      '<sup>',
                                      title_begin,
                                      collapse = "",
                                       '</sup>'
                                     ),
                     legend = list(title =  paste0("Threshold: \n adj.p < ", 
                                                               input$adj.P.Val, 
                                                               " \n and \n log2FC +/- ", 
                                                               input$logFC)#,
                                                #side = "left", 
                                               # bordercolor = '#444', 
                                               # bgcolor = "blue"
                                )
                     )
    pp
  })
  
  output$boxPlotUp <- shiny::renderPlot({
    shiny::validate(need(input$up_rows_selected, message = "Select up- or downregulated proteins from the table to generate a display of protein abundance in the selected groups."))
    input$analyzeLimma 
    input$up_rows_selected
    input$labelColBox
    original = proteinAbundance$original %>% tibble::column_to_rownames("Gene names") %>% as.data.frame()
    isolate(if(ClinColClasses()[input$GR_fatcor]=='factor' | input$ContinChoice == FALSE) {
      experimentalDesign = ClinDomit$designMatrix[,1:2] %>% tibble::rownames_to_column(var = "PatientID") %>% tidyr::gather(key= group, value = value, -PatientID ) %>% dplyr::filter(value > 0) %>% dplyr::select(-value) #%>% left_join(., ClinData()[, c("PatientID", input$filter_GR_fatcor, input$labelColBox)])
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
    shiny::validate(need(input$down_rows_selected, FALSE))
    original = proteinAbundance$original %>% tibble::column_to_rownames("Gene names") %>% as.data.frame()
    if(ClinColClasses()[input$GR_fatcor]=='factor' | input$ContinChoice == FALSE) {
      experimentalDesign =  ClinDomit$designMatrix[,1:2] %>% tibble::rownames_to_column(var = "PatientID") %>% tidyr::gather(key= group, value = value, -PatientID ) %>% dplyr::filter(value > 0) %>% dplyr::select(-value) #%>% left_join(., ClinData()[, c("PatientID", input$filter_GR_fatcor, input$labelColBox)])
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
    PD = tibble::as_tibble(as.data.frame(t(proteinData)))
    PD$PatientID = PatientID
    CD = CD_clean %>% dplyr::select(PatientID, GR_fatcor)
    #ClinDomit$designMatrix = ClinDomit$data
    CD = CD %>% dplyr::filter(PatientID %in% rownames(ClinDomit$designMatrix))
    #ClinDomit$designMatrix
    PD = dplyr::left_join(PD, CD) %>% 
      dplyr::select(c("PatientID", 
                      #"less.than.or.equal.to_61"
                      GR_fatcor, 
                      #!! parse_expr(GR_fatcor), 
                      limmaResult[s_up,]$`Gene names`)) %>% 
      tidyr::gather(key = "Gene name", value = "log2 of LFQ value", c(limmaResult[s_up,]$`Gene names`))
    if (is.null(labelColBox) | labelColBox == "" | labelColBox =="none") {
      groups = NULL
      labels = PD$PatientID
    } else if (!is.null(labelColBox) && labelColBox != "PatientID") {
      message(labelColBox)
      boxPlotSamples = tibble::tibble(PatientID = PD$PatientID)
      #boxPlotSamples$PatientID = gsub(pattern = "\\.", replacement = " ", x = boxPlotSamples$PatientID)
      boxPlotSamples$PatientID <- readr::parse_factor(boxPlotSamples$PatientID, include_na = FALSE)
      coloringFactor = CD_clean[,c("PatientID", labelColBox)]
      groups = dplyr::pull(dplyr::right_join(coloringFactor, boxPlotSamples)[,2])
      labels = PD$PatientID
    } else {
      message("Please select a different parameter for colouring")
      groups = NULL
      labels = PD$PatientID
    }
    scatterPlot = ggplot2::ggplot(PD, aes(x = get(GR_fatcor), y = `log2 of LFQ value`)) + 
      ggplot2::geom_point(aes(x = get(GR_fatcor), y = `log2 of LFQ value`, colour = groups)) + 
      ggplot2::geom_smooth(method = "lm") +
      ggplot2::labs(x = paste(GR_fatcor)) +
      ggplot2::facet_wrap(~`Gene name`) +
      ggplot2::theme_light()
    if (input$showLabels == TRUE) {
      scatterPlot = scatterPlot + ggrepel::geom_text_repel(aes(label=labels), show.legend = F, size = 4)
    }
    if (is.numeric(groups)){
      scatterPlot  = scatterPlot + ggthemes::scale_color_continuous_tableau() 
    } else {
      scatterPlot  = scatterPlot + ggthemes::scale_color_tableau()
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
                      limmaResult[s_up,]$`Gene names`)) %>% 
      tidyr::gather(key = "Gene name", value = "log2 of LFQ value", c(limmaResult[s_up,]$`Gene names`))# %>% 
    PD[,GR_fatcor] = readr::parse_factor(dplyr::pull(PD, GR_fatcor), include_na = FALSE)
    if (is.null(labelColBox) | labelColBox == "" | labelColBox =="none") {
      groups = NULL
      labels = PD$PatientID
    } else if (!is.null(labelColBox) && labelColBox != "PatientID") {
      message(labelColBox)
      boxPlotSamples = tibble::tibble(PatientID = PD$PatientID)
      boxPlotSamples$PatientID <- readr::parse_factor(boxPlotSamples$PatientID, include_na = FALSE)
      coloringFactor = ClinData[,c("PatientID", labelColBox)]
      groups = dplyr::pull(dplyr::right_join(coloringFactor, boxPlotSamples)[,2, drop = FALSE])
      labels = PD$PatientID
    } else {
      message("Please select a different parameter for colouring")
      groups = NULL
      labels = PD$PatientID
    }
    
    BoxPlot = ggplot2::ggplot(PD, aes(x = get(GR_fatcor), y = `log2 of LFQ value`)) + 
      ggplot2::geom_boxplot() +
      ggplot2::geom_jitter(aes(colour = groups)) +
      #geom_text_repel(aes(label=labels), show.legend = F, size = 4) +
      ggplot2::labs(x = paste(GR_fatcor)) +
      ggplot2::facet_wrap(~`Gene name`) +
      ggplot2::theme_light()
    
    if (input$showLabels == TRUE) {
      BoxPlot = BoxPlot + ggrepel::geom_text_repel(aes(label=labels), show.legend = F, size = 4)
    }
    if (is.numeric(groups)){
      BoxPlot  = BoxPlot + ggthemes::scale_color_continuous_tableau() 
    } else {
      BoxPlot  = BoxPlot + ggthemes::scale_color_tableau()
    }
    BoxPlot
  }
  
  output$up <- DT::renderDataTable({
    df <- limmaResult$df_up
    shiny::validate(need(limmaResult$df_up, message = "Please configure an experiment first."))
    DT::datatable(df[, c("Gene names", "logFC","P.Value","adj.P.Val", "Majority protein IDs")],
                  rownames = FALSE,
                  class = 'cell-border stripe',
                  width= '150px',
                  extensions ='Scroller', 
                  options = list(
                    deferRender = TRUE,
                    scrollY = 100,
                    scroller = TRUE,
                    scrollCollapse = TRUE,
                    pageLength = 100, 
                    lengthMenu = c(5,10,50,100,200)
                  )) %>% 
      DT::formatRound(columns = 2, digits = 2) %>% 
      DT::formatRound(columns = c(3:4), digits = 5)
    
  })
  
  output$down <- DT::renderDataTable({
    
    df<-limmaResult$df_down
    shiny::validate(need(limmaResult$df_down, message = "Please configure an experiment first."))
    DT::datatable(df[, c("Gene names", "logFC", "adj.P.Val", "P.Value", "Majority protein IDs")],
                  rownames = FALSE,
                  class = 'cell-border stripe', 
                  extensions ='Scroller',
                  options = list(
                    deferRender = TRUE,
                    scrollY = 100,
                    scroller = TRUE,
                    scrollCollapse=TRUE,
                    pageLength = 100, 
                    lengthMenu = c(5, 10, 50, 100, 200)
                  )) %>% 
                    DT::formatRound(columns = 2, digits = 2) %>% 
                    DT::formatRound(columns = c(3:4), digits = 5)

  })
  
  output$reportDataDL <- shiny::downloadHandler(
    filename = function(){
      paste('EatomicsData', input$GR_fatcor, names(ClinDomit$designMatrix)[1] , "vs.", names(ClinDomit$designMatrix)[2], Sys.Date(), '.csv', sep = '') },
    content = function(file) {
      x = list("Upregulated.Proteins" = limmaResult$df_up, 
               "Downregulated.Proteins" = limmaResult$df_down,
               "Limma.ExpDesign" = ClinDomit$designMatrix %>% tibble::rownames_to_column("PatientID"),
               "Limma.Setup.Details" = data.frame("imputed Data" = input$imputeForLimma, "eBayesTrend" = "TRUE", "Contrast" = paste(names(ClinDomit$designMatrix)[1]," regulated when compared to ", names(ClinDomit$designMatrix[2]))),
               "ProteinIDs.Gene.Mapping" = reportBlocks$ProteinIDMap,
               "Raw.data.incl.imputed.values" = proteinAbundance$imputed)
      openxlsx::write.xlsx(x, file, row.names = FALSE)
      grDevices::dev.off()
    }
  )
  
  
  output$doc1 <- shiny::renderUI({
    req(ClinDomit$designMatrix, limmaResult$gene_list)
    
    kable_input<-knitr::kable(ClinDomit$designMatrix) %>% kableExtra::kable_styling(bootstrap_options = "striped") 
    
    reportBlocks$ExpSetup <- 
      shiny::HTML(markdown::markdownToHTML(fragment.only=TRUE, text=c( 
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
         kableExtra::scroll_box(kable_input,width = "70%", height = "200px"),
        "* Total number of genes in limma are:",nrow(limmaResult$gene_list),
        "* The number of genes **upregulated** and **downregulated** in ",input$levels[1], " are ", nrow(limmaResult$up)," and " ,nrow(limmaResult$down), " respectively. ",
        "* The threshold used to highlight significant genes is [BH corrected](https://www.rdocumentation.org/packages/stats/versions/3.5.2/topics/p.adjust) adjusted P value of" , input$adj.P.Val, "and absolute log fold change of ",input$logFC 
        
      )))
    shinyBS::bsCollapsePanel(p("Detailed description", style = "color:#18bc9c"),
                    reportBlocks$ExpSetup)
  })
  
  output$report <- shiny::downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = paste('EatomicsReport-', Sys.Date(), '.html', sep = ''),
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
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
              distributionPlot = QCreport$distributionPlot,
              linesFromMaxQuant = reportBlocks$linesFromMaxQuant,
              stats_proteinGroups = reportBlocks$stats_proteinGroups,
              volcano_plot = reportBlocks$volcano_plot,
              boxPlotUp = reportBlocks$boxPlotUp,
              boxPlotDown = reportBlocks$boxPlotDown,
              ExpSetup = reportBlocks$ExpSetup,
              UpRegul = limmaResult$df_up,
              DoRegul = limmaResult$df_down
              )
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  ) 
  

  
  ###3 ssGSEA tab 
  
  ssgsea_data = shiny::reactiveValues(file = NULL)
  output$output.prefix <- shiny::renderUI({ 
    shiny::textInput("output.prefix",label = "Insert a Prefix for Output Files", input$gs.collection )

  })

  # run the ssGSEA analysis
  #getssgseaObj = eventReactive(input$goButton,{


  shiny::observeEvent(input$goButton, {
    
    shiny::validate(need(proteinAbundance$original, "Please upload proteomics data first (previous tab)."))
    original = proteinAbundance$original %>% tibble::column_to_rownames("Gene names") %>% as.data.frame()
    ssgsea_data$data = as.matrix(original)
        ssgsea_data$gs_collection = input$gs.collection
        ssgsea_data$sample.norm.type = input$sample.norm.type
        ssgsea_data$weight = input$weight
        ssgsea_data$statistic = input$statistic
        ssgsea_data$output.score.type = input$output.score.type
        ssgsea_data$nperm = input$nperm
        ssgsea_data$min.overlap = input$min.overlap
        ssgsea_data$correl.type = input$correl.type
        ssgsea_data$output.prefix = input$output.prefix
        ssgsea_data$linesFromMaxQuant = reportBlocks$linesFromMaxQuant
        ssgsea_data$stats_proteinGroups = reportBlocks$stats_proteinGroups

        shiny::withProgress(message = 'Calculation in progress',
                            detail = 'An alert notification will appear upon download of the file', value = 1, {
                              
                              ssgsea_obj = ssGSEA2(input.ds = ssgsea_data$data, gene.set.databases = gene.set.databases[input$gs.collection], sample.norm.type = input$sample.norm.type,
                                                   weight = input$weight, statistic = input$statistic, output.score.type = input$output.score.type, 
                                                   nperm = input$nperm, min.overlap = input$min.overlap, correl.type = input$correl.type, output.prefix = input$output.prefix, directory = sessionID)           
                            }) 

    ssgsea_data$prefix = ssgsea_obj
    shinyalert::shinyalert("Enrichment scores are ready - proceed to the next tabpanel.", type = "success", showConfirmButton = TRUE, timer = 5000)                     

  })
  
  
  
  ##download the file
  output$download<- shiny::downloadHandler(
    filename = function() { paste('efix, .csv', sep='') },
    content = function(file) {
      utils::write.csv(getssgseaObj(), file)
      #descs=getssgseaObj()$gs.descs.2, filename)
    }
  ) 
  
  
  
  
  ###4.Differential GSEA tab

  gs_file_list = observeEvent(ssgsea_data$prefix,{
    
    if(is.null(ssgsea_data$prefix)){
      ssgsea_data$file = NULL}
    else{
      ssgsea_data$file = list.files(path = paste(sessionID, "/", sep = ""))
    }
  })
  
  callModule(module = expDesignModule, id = "gsea", 
             ssgsea_data_update = reactive(ssgsea_data$prefix),
             gs_file_list = reactive(ssgsea_data$file),
             sessionID = sessionID,
             ClinData = ClinData, 
             reportData = ssgsea_data
             #ClinData = reactive(ClinDomit$data)
  )
  
  
  
  ###5. About tabpanel
  
 # output$markdown <- shiny::renderUI({
#    shiny::includeHTML(paste(homeDir, "/../Vignette/Help.html", sep = ""))
#  })
  
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
shiny::shinyApp(ui = ui, server = server)

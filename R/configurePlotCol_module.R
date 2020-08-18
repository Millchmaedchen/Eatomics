configurePlotColUI <- function(id, ...) {
  ns <- NS(id)
  shiny::uiOutput(ns("labelCol_out"))
}

configurePlotCol <- function(input, output, session,
                        ClinData) {
  ns <- session$ns
  output$labelCol_out <- shiny::renderUI({
    isolate(
      shiny::conditionalPanel(condition = need(ClinData, FALSE) , 
                              shiny::selectInput(
                                inputId = ns("labelCol_d"), 
                                label = strong("Choose the clinical parameter for group colours"),
                                choices = as.list("none" = "none", colnames(ClinData)),
                                multiple = FALSE,
                                selectize = TRUE
                              )
    )

    )
  })
  dummy = input$labelCol_d
  return(reactive({input$labelCol_d}))
}


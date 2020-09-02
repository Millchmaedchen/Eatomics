downloadObjUI <- function(id, ...) {
  ns <- NS(id)
  tagList(
    #br(),
    shinyBS::bsCollapsePanel(p("Customize figure text for download", style = "color:#18bc9c"),

  shiny::textInput(inputId = ns("title"), label = "Type a custom the plot title for download", value = "" ),
  shiny::textInput(inputId = ns("subtitle"), label = "Type a custom the plot subtitle for download", value = "" ),
  shiny::textInput(inputId = ns("caption"), label = "Type a custom the plot caption for download", value = "" )
    )
  )
}

downloadObj <- function(input, output, session,
                        title,
                        filename) {
  plot_parameters = list("title" = title, "subtitle" = "", "caption" = "", "filename" = filename)


  ns <- session$ns
   if(input$title != "") {
     plot_parameters$title = input$title
   }
  if(input$subtitle != "") {
    plot_parameters$subtitle = input$subtitle
  }
  if(input$caption != "") {
    plot_parameters$caption = input$caption
  }

  plot_parameters$filename = ifelse(input$filename != "", paste(input$filename, ".pdf"), filename)
  return(plot_parameters)
  }
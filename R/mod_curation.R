#' curation UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#'
mod_curation_ui <- function(id){
  ns <- NS(id)
  tagList(
    actionButton(ns("AB_curate_data"), label = "Curate data", class = "btn btn-large btn-danger")
  )
}

#' curation Server Functions
#'
#' @noRd
#'
mod_curation_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    # quick test
    observeEvent(input$AB_curate_data, {
      print("Hi there, I want to curate some data!")
    })

  })
}

## To be copied in the UI
# mod_curation_ui("curation_ui_1")

## To be copied in the server
# mod_curation_server("curation_ui_1")

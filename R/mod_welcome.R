#' welcome UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_welcome_ui <- function(id){
  ns <- NS(id)
  tagList(
    wellPanel(
      id = "about",
      includeMarkdown(app_sys("app/www/welcome.Rmd"))
    )
  )
}

#' welcome Server Functions
#'
#' @noRd
mod_welcome_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

  })
}

## To be copied in the UI
# mod_welcome_ui("welcome_ui_1")

## To be copied in the server
# mod_welcome_server("welcome_ui_1")

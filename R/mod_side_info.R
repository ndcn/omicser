#' side_info UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_side_info_ui <- function(id){
  ns <- NS(id)
  tagList(
    includeMarkdown(app_sys("app/www/info.Rmd"))
  )
}

#' side_info Server Functions
#'
#' @noRd
mod_side_info_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

  })
}

## To be copied in the UI
# mod_side_info_ui("side_info_ui_1")

## To be copied in the server
# mod_side_info_server("side_info_ui_1")

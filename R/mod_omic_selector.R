#' omic_selector UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_omic_selector_ui <- function(id){
  ns <- NS(id)
  tagList(
 
  )
}
    
#' omic_selector Server Functions
#'
#' @noRd 
mod_omic_selector_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_omic_selector_ui("omic_selector_ui_1")
    
## To be copied in the server
# mod_omic_selector_server("omic_selector_ui_1")

#' pg_table UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_pg_table_ui <- function(id){
  ns <- NS(id)
  tagList(
 
  )
}
    
#' pg_table Server Functions
#'
#' @noRd 
mod_pg_table_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_pg_table_ui("pg_table_ui_1")
    
## To be copied in the server
# mod_pg_table_server("pg_table_ui_1")

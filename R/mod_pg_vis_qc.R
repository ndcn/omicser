#' pg_vis_qc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_pg_vis_qc_ui <- function(id){
  ns <- NS(id)
  tagList(

  )
}

#' pg_vis_qc Server Functions
#'
#' @noRd
mod_pg_vis_qc_server <- function(id,rv_data, rv_selections){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

  })
}

## To be copied in the UI
# mod_pg_vis_qc_ui("pg_vis_qc_ui_1")

## To be copied in the server
# mod_pg_vis_qc_server("pg_vis_qc_ui_1")

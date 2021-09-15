#' pg_tables_tab UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_pg_tables_tab_ui <- function(id){
  ns <- NS(id)
  tagList(



    tabsetPanel(
      type = 'pills',  #'hidden' and a radio might work best
      id = 'tab',

      tabPanel(
        title = "Omic Meta", value = 'var',
        mod_pg_table_ui(ns("pg_table_var"))
      ),

      tabPanel(
        title = "Sample Meta",value = 'obs',
        mod_pg_table_ui(ns("pg_table_obs"))
      ),

      tabPanel(
        title = "diff expr", value='de',
        mod_pg_table_ui(ns("pg_table_de"))
      ),

      # tabPanel(
      #   title = "Data Matrix",value = 'X',
      #   mod_pg_table_ui(ns("pg_table_X"))
      # )

    ) #tabsetpanel


  )
}

#' pg_tables_tab Server Functions
#'
#' @noRd
mod_pg_tables_tab_server <- function(id,rv_in, p){
  moduleServer( id, function(input, output, session){
    ns <- session$ns


    dt_obs <- reactive({
      req(rv_in$ad)
        rv_in$ad$obs
    })

    dt_var <- reactive({
      req(rv_in$ad)
      rv_in$ad$var
    })

    dt_de <- reactive({
      req(rv_in)
      rv_in$de
    })


    mod_pg_table_server("pg_table_obs",dt_obs)
    mod_pg_table_server("pg_table_var",dt_var)
    mod_pg_table_server("pg_table_de",dt_de)


  })
}

## To be copied in the UI
# mod_pg_tables_tab_ui("pg_tables_tab_ui_1")

## To be copied in the server
# mod_pg_tables_tab_server("pg_tables_tab_ui_1")

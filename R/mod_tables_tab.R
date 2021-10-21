#' tables_tab UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_tables_tab_ui <- function(id){
  ns <- NS(id)
  tagList(



    tabsetPanel(
      type = 'pills',  #'hidden' and a radio might work best
      id = 'tab',

      tabPanel(
        title = "Feature Meta", value = 'var',
        mod_table_ui(ns("pg_table_var"))
      ),

      tabPanel(
        title = "Sample Meta",value = 'obs',
        mod_table_ui(ns("pg_table_obs"))
      ),

      tabPanel(
        title = "Diff. Expr.", value='de',
        mod_table_ui(ns("pg_table_de"))
      )

      # tabPanel(
      #   title = "Data Matrix",value = 'X',
      #   mod_pg_table_ui(ns("pg_table_X"))
      # )

    ) #tabsetpanel


  )
}

#' tables_tab Server Functions
#'
#' @param id shiny internal
#' @param rv_data reactive data
#' @param rv_selections side selector reactives
#'
#' @noRd
mod_tables_tab_server <- function(id,rv_data, rv_selections){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    dt_obs <- reactive({
      req(rv_data$anndata)
        rv_data$anndata$obs
    })

    dt_var <- reactive({
      req(rv_data$anndata)
      rv_data$anndata$var
    })

    dt_de <- reactive({
      req(rv_data)
      rv_data$de
    })


    mod_table_server("pg_table_obs",dt_obs)
    mod_table_server("pg_table_var",dt_var)
    mod_table_server("pg_table_de",dt_de)


  })
}

## To be copied in the UI
# mod_tables_tab_ui("tables_tab_ui_1")

## To be copied in the server
# mod_tables_tab_server("tables_tab_ui_1")

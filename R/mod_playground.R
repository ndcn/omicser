#' playground UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_playground_ui <- function(id){
  ns <- NS(id)
  tagList(
    # TODO: change these to dynamically render or not based on "selector"
    # # could be "side selector" or some radio/dropdown choices...
    #

    # shinyWidgets::verticalTabsetPanel(
    #   id = 'tab',
    #   shinyWidgets::verticalTabPanel(
    #     title = "Table", value='table',
    #     mod_pg_table_ui(id=ns("pg_table_ui_2"))
    #   ),
    #   # ingest tab
    #   shinyWidgets::verticalTabPanel(
    #     title = "Values", value = 'raw',
    #     mod_pg_vis_raw_ui(id=ns("pg_vis_raw_ui_1"))
    #   ),
    #   # table tab
    #   shinyWidgets::verticalTabPanel(
    #     title = "Comparative",value = 'comp',
    #     mod_pg_vis_comp_ui(id=ns("pg_vis_comp_ui_1"))
    #
    #   )
    # ) #verticalTabsetPanel
    tabsetPanel(
      type = 'pills',  #'hidden' and a radio might work best
      id = 'tab',
      tabPanel(
        title = "Table", value='table',
        mod_pg_table_ui(id=ns("pg_table_ui_2"))
      ),
      # ingest tab
      tabPanel(
        title = "Quantities", value = 'raw',
        mod_pg_vis_raw_ui(id=ns("pg_vis_raw_ui_1"))
      ),
      # table tab
      tabPanel(
        title = "Comparisons",value = 'comp',
        mod_pg_vis_comp_ui(id=ns("pg_vis_comp_ui_1"))

      )
    ) #tabsetpanel


  )
}

#' playground Server Functions
#'
#' @noRd
mod_playground_server <- function(id ,rv_in, p) {
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    mod_pg_table_server("pg_table_ui_2",rv_in, p)
    mod_pg_vis_raw_server("pg_vis_raw_ui_1",rv_in, p)
    mod_pg_vis_comp_server("pg_vis_comp_ui_1",rv_in, p)


  })
}

## To be copied in the UI
# mod_playground_ui("playground_ui_1")

## To be copied in the server
# mod_playground_server("playground_ui_1")


#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # Your application server logic



  ############################ +
  ## Module 4 : "Ingest" Data
  ##
  ##
  ############################ +
  {
    # Call module "ingest" - returns reactive data values
    rv_data <- mod_ingestor_server("ingestor_ui_1")

  }

  ############################ +
  ## Module 1 : SIDESELECT   ####
  ##
  ##
  ############################ +
  {
    rv_selections <- mod_side_selector_server("side_selector_ui_1", rv_data = rv_data)

  }

  ############################ +
  ## Module 2 :  info for alternative sidebar
  ##
  ############################ +
  {
    mod_side_info_server("side_info_ui_1")
  }

  ############################ +
  ## Module 3,: Welcome
  ##
  ############################ +
  {
    mod_welcome_server("welcome_ui_1")
  }


  ############################ +
  ## Module 5 : table
  ##
  ##
  ############################ +
  {

    mod_tables_tab_server("tables_tab_ui_1", rv_data = rv_data, rv_selections = rv_selections)

  }
  # selected_db <- reactive_selection_function(db,params)
  # mod_playground_server("playground_ui_1",db = data_table)


  ############################ +
  ## Module 6 : playground   ####
  ##     id call = "playground_ui_1"
  ##     calls:
  ##      6a: mod_pg_vis_comp
  ##      6b: mod_pg_vis_raw
  ##      6c: mod_pg_table?
  ##      6d: mod_pc_qc
  ##
  ##
  ##     ###+
  ############################ +
  {
    # mod_playground_server("playground_ui_1", rvs = rv, p = vis_params  )
    mod_playground_server("playground_ui_1", rv_data = rv_data, rv_selections = rv_selections)
  }

  ############################ +
  ## Module 7 : Help   ####
  ##
  ##
  ############################ +
  {
    mod_help_server("help_ui_1")
  }

  # ############################ +
  # ## Module 8 : Export   ####
  # ##
  # ##
  # ############################ +
  # {
  #   mod_export_server("export_ui_1", rv_in = rv, p = p_vis)
  # }

}


#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # Your application server logic

  # config is in the ./inst/db_info/ directory
  # for now we have NOT added these config options to the golem-config.yml
  #CONFIG <- configr::read.config( "./omxr_options.yml" )
  #


  CONFIG <- omicser::get_config()

  DB_NAMES <- CONFIG$database_names
  CONDA_ENV <- CONFIG$conda_environment
  DB_ROOT_PATH <- CONFIG$db_root_path
  #CODA_EXE <- CONFIG$conda_exe
  #TODO: curate CONDA_EXE?  Maybe set up an environment variable?
  #      is there a tool to read from sys env?

  #required = TRUE, or explicitly set the RETICULATE_PYTHON environment variable.

  reticulate::use_condaenv(
      required = TRUE,
      condaenv = CONDA_ENV,
      conda = "auto") #TODO: change this to CONFIG$conda_exe

  ############################ +
  ## Module 4 : "Ingest" Data
  ##
  ##
  ############################ +
  {
    # Call module "ingest" - returns reactive data values
    rv_data <- mod_ingestor_server("ingestor_ui_1", db_names=DB_NAMES, db_root_path=DB_ROOT_PATH)

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
    # output$my_datatable_0 <- DT::renderDataTable({
    #   DT::datatable(rv$data_table)
    # })
    #   Should we pass all the reactive values, or just the datatable?
    #   currently returns a subtable of the datatable
    # filtered_db <- mod_table_server("table_ui_1", dt = rv$data_table, p = vis_params)
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

  # ############################ +
  # ## Module 7 : Export   ####
  # ##
  # ##
  # ############################ +
  # {
  #   mod_export_server("export_ui_1", rv_in = rv, p = p_vis)
  # }

}

# hardwire for now


CONFIG <- configr::read.config( "./omxr_options.yml" )
print(getwd())

DATASET_NAMES <- CONFIG$dataset_names
CONDA_ENV <- CONFIG$conda_environment
DS_ROOT_PATH <- CONFIG$ds_root_path


#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # Your application server logic
  print(getwd())

  reticulate::use_condaenv(
                          required = TRUE,
                          condaenv =  CONDA_ENV
                          )


  ############################ +
  ## Module 4 : "Ingest" Data
  ##
  ##
  ############################ +
  {
    # Call module "ingest"
    rv <- mod_ingestor_server("ingestor_ui_1")

    # rv <- mod_ingest_server("ingest_ui_1")
    # print(isolate(rv$data_table))
    # Ingest "LOAD" button triggers sharing the rv structure:

    #   rv$data_table - datatable
    #   rv$database_name  - internal name of db
    #   rv$gene_names - full list of genes
    #   rv$var0 (column name)
    #   rv$factors0  - factors in that column
    #   rv$use_var1 - logical flag to idnciate if var1/factor 1 are used
    #   rv$var1  (column name)
    #   rv$factors1  - factors in that column
    #   rv$trigger
  }




  ############################ +
  ## Module 1 : SIDESELECT   ####
  ##
  ##
  ############################ +
  {
    p_vis <- mod_side_selector_server("side_selector_ui_1", rv_in = rv)

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
    mod_pg_table_server("pg_table_ui_1")
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
    mod_playground_server("playground_ui_1", rv_in = rv, p = p_vis)
  }

  ############################ +
  ## Module 7 : Export   ####
  ##
  ##
  ############################ +
  {
    mod_export_server("export_ui_1", rv_in = rv, p = p_vis)
  }

  ############################ +
  ## Module 8 : Additional Info   ####
  ##
  ##
  ############################ +
  {
    mod_additional_info_server("additional_info_ui_1")
  }
}


#' ingestor UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_ingestor_ui <- function(id) {
  ns <- NS(id)


  tagList(
    fluidRow(
      column(
        width=12,
        selectizeInput(
          ns("SI_database"), "Database",
          "",multiple=FALSE, options = list(placeholder = " database first")
          )

        )
      ),

    fluidRow(
      column(
        width=3,
        shinyjs::disabled(
          actionButton(ns("AB_ingest_load"), label = "Load Database")
        )
      )
    ),
    fluidRow(
      column(
        width=3,
        textOutput(ns("ui_datatype"))
      )
    ),
    mod_additional_info_ui(id = ns("additional_info_ui_ingest")),

    fluidRow(
      column(
        width=6,
        uiOutput(ns("ui_omics"))
      )
    ),

    fluidRow(
      helpText(
        HTML("<i>observations</i>  <b> experimental</b>.")
      )
      ),
    fluidRow(
      column(
          width=8,
          textOutput(ns("ui_obs_exp"))
        )
    )
  ) #tagList
}


#' ingestor Server Functions
#'
#' @noRd
mod_ingestor_server <- function(id,DB_NAMES, DB_ROOT_PATH) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # # TODO: deltet this or make it dynamic?
    # database_names =  golem::get_golem_options( "database_names" )
    # ds_root_path =  golem::get_golem_options( "ds_root_path" )
    # global from omxr_options.yml (app_server.R)
    # CONFIG <- configr::read.config( "./omxr_options.yml" )
    #
    # DB_NAMES <- CONFIG$database_names
    # CONDA_ENV <- CONFIG$conda_environment
    # DB_ROOT_PATH <- CONFIG$ds_root_path
    database_names <- DB_NAMES

    db_name = reactiveValues(name=NULL)
    observe({
      db_name$name <- to_return$db_meta$name
    })
print("making additinal _info")
    mod_additional_info_server("additional_info_ui_ingest",
                               db_name = db_name,
                               DB_ROOT_PATH = DB_ROOT_PATH)

    ## rv_data (to_return) REACTIVE VALUES  ===================
    to_return <- reactiveValues(
      # these values hold the database contents (only reactive because we can choose)
      database_name = NULL,
      omics_type = NULL, # Transcript-omics, Prote-omics, Lipid-omics, Metabol-omics, misc-
      measure = NULL, #normalized count, concentration, composotion, etc.

      # should this be packed into ad$uns?  written into .yml
      # provides context and  plot labeling choices...
      # TODO: makefunctions like ShinyProt "protein_conversions.R"
      db_meta = list( name = NULL,
                      omics_type = NULL,
                      measurment = NULL,
                      organism = NULL,
                      publication = NULL,
                      etc = NULL
                      ),

      #  Everything is packed into an anndata object
      ad = NULL,

      # omics key feature i.e. genes, proteins, lipids
      omics = NULL, #the omics columnname...

      config = NULL,
      default = NULL,
      meta = NULL,
      de = NULL,
      trigger = 0
    )

    updateSelectizeInput(session, "SI_database", choices = database_names, selected = DB_NAMES[1], server=TRUE)

## load dataset (observeEvent) ===================
    observeEvent(input$SI_database, {
      req(input$SI_database)
     print("in observeEvent(input$SI_database, ")
      db_name <- (input$SI_database)
      # ad <- anndata::read_h5ad(filename=paste0("data-raw/",db_name,"/omxr_data.h5ad"))
      # diff_exp = readRDS(file = paste0("data-raw/",db_name,"/diff_expr_table.rds"))
      #globals: DS_ROOT_PATH
      #
      ad <- anndata::read_h5ad(filename=file.path(DB_ROOT_PATH,db_name,"db_data.h5ad"))
      diff_exp = readRDS(file = file.path(DB_ROOT_PATH,db_name,"db_de_table.rds"))
      conf_def <- gen_config_table(ad, db_name, DB_ROOT_PATH)

      omics <- ad$var_names
      names(omics) <- omics

      to_return$de <- diff_exp

      # Rico: why is this here?
      # to_return$omics_type <- to_return$db_meta$omics_type
      to_return$ad <- ad
      to_return$omics <- omics
      # to_return$meta <- omicmeta  # this might be too redundant

      to_return$config <- conf_def$conf
      to_return$default <- conf_def$def

      to_return$db_meta$name<-db_name
      to_return$db_meta$omics_type<-conf_def$def$db_meta$omic_type
      to_return$db_meta$measurment<-conf_def$def$db_meta$measurement
      to_return$db_meta$organism<-conf_def$def$db_meta$organizm
      to_return$db_meta$publication<-conf_def$def$db_meta$pub
      to_return$db_meta$etc<-conf_def$def$db_meta$annotation_database

    })

    ## observes ===================
    print("pre database_name")

    observe({
      if ( !is.null(to_return$database_name) ) {
        shinyjs::enable("AB_ingest_load")
      } else {
        shinyjs::disable("AB_ingest_load")
        print(" no database loaded... .")
      }
    })
    print("post database_name")

    # keep these up to date for side_select..
    # # Should this just be done in     observeEvent(input$SI_database,?
    observe({
      #to_return$database_name <- input$SI_database
      to_return$database_name  <- names(which(database_names==input$SI_database))
    })
    # load button :: send the reactive object back to the app...
    observeEvent(input$AB_ingest_load, {
      # all other return values set with SI_database
      to_return$trigger <- to_return$trigger + 1
    })

    ## render info (TODO:) ===================

# TODO: clean up these renderPrint...
    output$ui_obs_exp <- renderPrint({
      obs_choices <- isolate(to_return$ad$obs_keys())
      if (is.null(obs_choices)) {
        print("(obs_exp) no datbase loaded")
      } else {
        print(paste0("observ exp: ", paste(obs_choices[1:min(10,length(obs_choices))], collapse = ","),"... ") )
      }
    })



    # show the factors that have been loaded
    # # TODO: render a nicely formatted version of the meta data and all this stuff...
    output$ui_datatype <- renderText({
      print(paste0("db type - ", to_return$db_meta$omics_type , "-omics"))
    })

    output$ui_omics <- renderPrint({
      if (is.null(to_return$omics)) {
        print("no datbase loaded")
      } else {
        omics <- isolate(to_return$omics)
        print(paste0("Current database omics: ", paste( names(omics)[1:10], collapse = ",")) )
      }
    })



    return(to_return)
  })
}

## To be copied in the UI
# mod_ingestor_ui("ingestor_ui_1")

## To be copied in the server
# mod_ingestor_server("ingestor_ui_1")

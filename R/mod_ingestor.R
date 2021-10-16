
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
    h4("Load databases for browsing and exploration"),

    fluidRow(
        selectizeInput(ns("SI_database"),
                            "Database",
                            "", multiple=FALSE, options = list(placeholder = " database first")
                            ),
        shinyjs::disabled(
               actionButton(ns("AB_ingest_load"), label = "Load Database")
                      )

      ),

    hr(style = "border-top: 1px dashed grey;"),
    fluidRow(
      column(width = 2,
             offset = 0,
             "DATABASE",br(),
             "info"
      ),
      column(width = 10,
             offset = 0,
             uiOutput(ns("ui_data_meta")),
      )
    ),
    mod_additional_info_ui(id = ns("additional_info_ui_ingest"))
  ) #tagList
}


#' ingestor Server Functions
#'
#' @noRd
mod_ingestor_server <- function(id,db_names, db_root_path) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # # TODO: deltet this or make it dynamic?
    # database_names =  golem::get_golem_options( "database_names" )
    # ds_root_path =  golem::get_golem_options( "ds_root_path" )
    # global from omxr_options.yml (app_server.R)
    # CONFIG <- configr::read.config( "./omxr_options.yml" )
    #
    # db_names <- CONFIG$database_names
    # CONDA_ENV <- CONFIG$conda_environment
    # db_root_path <- CONFIG$ds_root_path
    database_names <- db_names

    db_name = reactiveValues(name=NULL,
                             dir=NULL)
    observe({
      db_name$name <- to_return$db_meta$name
    })

    mod_additional_info_server("additional_info_ui_ingest",
                               db_name = db_name,
                               db_root_path = db_root_path)

    to_return <- reactiveValues(
      # these values hold the database contents (only reactive because we can choose)
      database_name = NULL,
      omics_type = NULL, # Transcript-omics, Prote-omics, Lipid-omics, Metabol-omics, misc-

      # TODO:
      omic_units = NULL, #normalized count, concentration, composotion, etc.

      #  Everything is packed into an anndata object
      anndata = NULL,
      de = NULL,

      # omics key feature i.e. genes, proteins, lipids
      omics = NULL, #the omics columnname...

      config = NULL,
      default = NULL,


      trigger = 0
    )

    updateSelectizeInput(session, "SI_database", choices = database_names, selected = db_names[1], server=TRUE)


    loaded_data <- reactiveValues(
        anndata = NULL,
        de = NULL, # Tra
        conf_def = NULL
    )

    ## load dataset (observeEvent) ===================
    observeEvent(input$SI_database, {
      req(input$SI_database)

      print("making loaded_data rv")
      db_name$dir <- input$SI_database
      db_name$name <- names(which(database_names==input$SI_database))
      loaded_data$de <- readRDS(file = file.path(db_root_path,db_name$dir,"db_de_table.rds"))
      ad <- anndata::read_h5ad(filename=file.path(db_root_path,db_name$dir,"db_data.h5ad"))
      loaded_data$anndata <- ad
      loaded_data$conf_def <- gen_config_table(ad, db_name$dir, db_root_path)
    })
    # is there overhead to store these in a reactive conduit?

    ## render info (TODO:) ===================

    # will this work?  or is the "observing" affecting a reactive with this render?
    output$ui_data_meta <- renderUI({
      req(loaded_data$conf_def,
          db_name$name)
      #omics_type <- loaded_data$conf_def$def$db_meta$omic_type
      #type_text <- HTML(paste("omics type: <b>", omics_type, " </b> x "))
      out_text <- HTML(paste("<b>omics type: </b> ,<i>",
                             loaded_data$conf_def$def$db_meta$omic_type,
                             " </i> from the <i>", loaded_data$conf_def$def$db_meta$lab),"</i>")

      ret_tags <-  tagList(
        h4("Database Info >> "),
        #"some more text",
        out_text,
        br()
      )
      return( ret_tags )
    })



    ## observes ===================

    observe({
      if ( !is.null( db_name$name ) ) {
        shinyjs::enable("AB_ingest_load")
      } else {
        shinyjs::disable("AB_ingest_load")
        print(" no database loaded... .")
      }
    })

    # keep these up to date for side_select..
    # # Should this just be done in     observeEvent(input$SI_database,?

    # load button :: send the reactive object over to the side selector.
    observeEvent(input$AB_ingest_load, {
      # all other return values set with SI_database

      to_return$de <- loaded_data$de
      to_return$anndata <- loaded_data$anndata

      to_return$config <- loaded_data$conf_def$conf
      to_return$default <- loaded_data$conf_def$def

      to_return$db_meta$name<-db_name$name
      to_return$db_meta$db_dir<-db_name$dir

      to_return$db_meta$omics_type<-loaded_data$conf_def$def$db_meta$omic_type
      to_return$db_meta$measurment<-loaded_data$cconf_def$def$db_meta$measurement
      to_return$db_meta$organism<-loaded_data$cconf_def$def$db_meta$organizm
      to_return$db_meta$publication<-loaded_data$cconf_def$def$db_meta$pub
      to_return$db_meta$etc<-loaded_data$cconf_def$def$db_meta$annotation_database

      to_return$trigger <- to_return$trigger + 1

    })




    return(to_return)
  })
}

## To be copied in the UI
# mod_ingestor_ui("ingestor_ui_1")

## To be copied in the server
# mod_ingestor_server("ingestor_ui_1")

require(anndata)


# TODO: pack this into a .rda .rds to load for dynamic updates
dataset_names <- c(
  "Domenico DIA" = "domenico_stem_cell",
  "Vilas Microglia" = "vilas_microglia",
  "Vilas Microglia (seu)" = "vilas_microglia_seu",
  "Yassene Lipid Concentrations" ="yassene_A_conc",
  "Yassene Lipid Compositions" ="yassene_A_compos",
  "Oscar Microglia" ="oscar_microglia"
)

dataset_type <- c(
  "Transcriptomics" = "transcript",
  "Proteomics" = "prote",
  "Lipidomics" = "lipid",
  "Metabolomics" = "metabol",
  "Other" = "X-"
)

observation_type <- c(
  "Raw Quantity" = "raw",
  "Comparitive Quantity" = "comp",
  "Grouping" = "group"
)
#require(anndata)

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
          ns("SI_dataset"), "Dataset",
          choices = dataset_names,
          select = NULL #"VilasB"  # should i have a default??
          )
        )
      ),
    fluidRow(
      column(
        width=3,
        textOutput(ns("ui_datatype"))
      )
    ),
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
    ),

    fluidRow(
      column(
        width=3,
        shinyjs::disabled(
          actionButton(ns("AB_ingest_load"), label = "Load Database")
        )
      )
    )
  ) #tagList
}


#' ingestor Server Functions
#'
#' @noRd
mod_ingestor_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns


    ############################ +
    ## initiate reactive database structure
    ##
    ############################ +
    to_return <- reactiveValues(
      # these values hold the database contents (only reactive because we can choose)
      database_name = NULL,
      omics_type = NULL, # Transcript-omics, Prote-omics, Lipid-omics, Metabol-omics, misc-

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


    ############################ +
    ## load dataset
    ############################ +
    observeEvent(input$SI_dataset, {
      if (!is.null(input$SI_dataset)) { # unnesscasary defensive?
        ds_name <- (input$SI_dataset)
        ds_label <- names(which(dataset_names==ds_name))
          if (ds_name == "vilas_microglia") { # Vilas transc
            ad <- anndata::read_h5ad(filename="data-raw/vilas_microglia/core_data_plus_de.h5ad")
            #ad <- anndata::read_h5ad(filename="data-raw/vilas_microglia/core_data.h5ad")
            diff_exp = readRDS(file ="data-raw/vilas_microglia/de_table.rds")
            to_return$de <- diff_exp

            omicconf <- omicser::vilas_microglia_conf
            omicdef <- omicser::vilas_microglia_def
            omics <- omicser::vilas_microglia_omics
            omicmeta <- omicser::vilas_microglia_meta



            to_return$database_name <- ds_label
            to_return$omics_type <- "transcript"
            to_return$ad <- ad
            to_return$omics <- omics
            to_return$meta <- omicmeta  # this might be too redundant


            to_return$config <- omicconf
            to_return$default <- omicdef

          } else if (ds_name == "vilas_microglia_seu") { # Vilas microglia
            ad <- anndata::read_h5ad(filename="data-raw/vilas_microglia_seu/core_data_plus_de.h5ad")
            #ad <- anndata::read_h5ad(filename="data-raw/vilas_microglia_seu/core_data.h5ad")
            diff_exp = readRDS(file ="data-raw/vilas_microglia_seu/de_table.rds")
            to_return$de <- diff_exp

            omicconf <- omicser::vilas_microglia_seu_conf
            omicdef <- omicser::vilas_microglia_seu_def
            omics <- omicser::vilas_microglia_seu_omics
            omicmeta <- omicser::vilas_microglia_seu_meta

            to_return$database_name <- ds_label
            to_return$omics_type <- "transcript"

            to_return$ad <- ad
            # set a default to avoid funny business in side_select
            to_return$omics <- omics
            to_return$meta <- omicmeta  # this might be too redundant

            to_return$config <- omicconf
            to_return$default <- omicdef

          } else if (ds_name == "yassene_A_conc") { # Vilas transc
            ad <- anndata::read_h5ad(filename="data-raw/yassene_A_conc/core_data_plus_de.h5ad")
            diff_exp = readRDS(file ="data-raw/yassene_A_conc/de_table.rds")
            to_return$de <- diff_exp

            omicconf <- omicser::yassene_A_conc_conf
            omicdef <- omicser::yassene_A_conc_def
            omics <- omicser::yassene_A_conc_omics
            omicmeta <- omicser::yassene_A_conc_meta

            to_return$database_name <- ds_label
            to_return$omics_type <- "lipid"

            to_return$ad <- ad
            to_return$omics <- omics
            to_return$meta <- omicmeta  # this might be too redundant

            to_return$config <- omicconf
            to_return$default <- omicdef

          } else if (ds_name == "yassene_A_compos") { # Vilas transc

            ad <- anndata::read_h5ad(filename="data-raw/yassene_A_compos/core_data_plus_de.h5ad")
            #ad <- anndata::read_h5ad(filename="data-raw/yassene_A_compos/core_data.h5ad")
            diff_exp = readRDS(file ="data-raw/yassene_A_compos/de_table.rds")
            to_return$de <- diff_exp


            omicconf <- omicser::yassene_A_compos_conf
            omicdef <- omicser::yassene_A_compos_def
            omics <- omicser::yassene_A_compos_omics
            omicmeta <- omicser::yassene_A_compos_meta

            to_return$database_name <- ds_label
            to_return$omics_type <- "lipid"

            to_return$ad <- ad
            to_return$omics <- omics
            to_return$meta <- omicmeta  # this might be too redundant

            to_return$config <- omicconf
            to_return$default <- omicdef
          } else if (ds_name == "domenico_stem_cell") { # "Domenico Tsc Long"

            ad <- anndata::read_h5ad(filename="data-raw/domenico_stem_cell/core_data_plus_de.h5ad")
            #ad <- anndata::read_h5ad(filename="data-raw/domenico_stem_cell/core_data.h5ad")
            diff_exp = readRDS(file ="data-raw/domenico_stem_cell/de_table.rds")
            to_return$de <- diff_exp

            omicconf <- omicser::domenico_stem_cell_conf
            omicdef <- omicser::domenico_stem_cell_def
            omics <- omicser::domenico_stem_cell_omics
            omicmeta <- omicser::domenico_stem_cell_meta

            to_return$database_name <- ds_label
            to_return$omics_type <- "prote"

            to_return$ad <- ad
            to_return$omics <- omics
            to_return$meta <- omicmeta  # this might be too redundant


            to_return$config <- omicconf
            to_return$default <- omicdef

          } else if (ds_name == "oscar_microglia") { # "Oscar Microglia"


            ad <- anndata::read_h5ad(filename="data-raw/oscar_microglia/core_data_plus_de.h5ad")
            #ad <- anndata::read_h5ad(filename="data-raw/oscar_microglia/core_data.h5ad")
            diff_exp = readRDS(file ="data-raw/oscar_microglia/de_table.rds")
            to_return$de <- diff_exp

            omicconf <- omicser::oscar_microglia_conf
            omicdef <- omicser::oscar_microglia_def
            omics <- omicser::oscar_microglia_omics
            omicmeta <- omicser::oscar_microglia_meta

            to_return$database_name <- ds_label
            to_return$omics_type <- "transcript"

            to_return$ad <- ad
            to_return$omics <- omics
            to_return$meta <- omicmeta  # this might be too redundant


            to_return$config <- omicconf
            to_return$default <- omicdef

          # } else if (ds_name == "OscarB") { # Vilas transc
          #
        # } else if (ds_name == "novel") {
        #   # TODO:  code here for upload module
        #   # upload VAR, OBS, X, meta with choosers for organizing code...
        #   to_return$obs <- NULL
        #   print("NO DATA LOADED")
        #   to_return$omics_type <- "NA"
        #   # leave everythin NULL
        } else { # TODO: code the ingests
          print("NO DATA LOADED")
          to_return$database_name <- NULL
          #to_return$omics_type <- "NA"
          # leave everythin NULL

        }
      }
    })


    observe({
      if ( !is.null(to_return$database_name) ) {
        shinyjs::enable("AB_ingest_load")
      } else {
        shinyjs::disable("AB_ingest_load")
        print(" no database loaded... try Vilas Trans or Domenico or Yassene.")
      }
    })


    output$ui_obs_exp <- renderPrint({
      obs_choices <- isolate(to_return$ad$obs_keys())
      if (is.null(obs_choices)) {
        print("(obs_exp) no datbase loaded")
      } else {
        print(paste0("observ exp: ", paste(obs_choices[1:min(10,length(obs_choices))], collapse = ","),"... ") )
      }
    })

    output$ui_omics <- renderPrint({
      if (is.null(to_return$omics)) {
        print("no datbase loaded")
      } else {
        omics <- isolate(to_return$omics)
        print(paste0("Current database omics: ", paste( names(omics)[1:10], collapse = ",")) )
      }
    })



    # keep these up to date for side_select..
    # # Should this just be done in     observeEvent(input$SI_dataset,?
    observe({
      #to_return$database_name <- input$SI_dataset
      to_return$database_name  <- names(which(dataset_names==input$SI_dataset))

    })

    # show the factors that have been loaded
    output$ui_datatype <- renderText({
      print(paste0("db type - ", to_return$omics_type , "-omics"))
    })


    # load button :: send the reactive object back to the app...
    observeEvent(input$AB_ingest_load, {
      # all other return values set with SI_dataset
      to_return$trigger <- to_return$trigger + 1
    })


    return(to_return)
  })
}

## To be copied in the UI
# mod_ingestor_ui("ingestor_ui_1")

## To be copied in the server
# mod_ingestor_server("ingestor_ui_1")

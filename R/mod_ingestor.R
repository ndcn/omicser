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

  # hardwire for now
  dataset_names <- c(
    "Skeletal Muscle" = "TMT",
    "Muscle Stem Cell" = "DIA",
    "Domenico Tsc Long" = "DomenicoA",
    "Vilas Transcript" = "VilasA",
    "Vilas Alt" = "VilasB",
    "Oscar 1" = "OscarA",
    "Upload New data" = "novel"
  )

  dataset_type <- c(
    "Transcriptomics" = "transcript",
    "Proteomics" = "prote",
    "Lipidomics" = "lipid",
    "Metabolomics" = "metabol",
    "Other" = "X-"
  )

  tagList(
    fluidRow(
      col_12( # full width
        selectizeInput(
          ns("SI_dataset"), "Dataset",
          choices = dataset_names,
          select = "DomenicoA"  # should i have a default??
        ),

        # Var -> Genes/ Proteins
        selectizeInput(
          ns("SI_var_name"), "Choose Primary Variable", "",
          multiple = FALSE, options = list(placeholder = "choose dataset first")
        ),
        shinyjs::disabled(checkboxInput(ns("CB_aux_vars"), "More Variable Annotations?", value = FALSE)),
        shinyjs::disabled(
          selectizeInput(ns("SI_aux_vars"), "Choose Variable Annotations", "",
            multiple = TRUE, options = list(placeholder = "choose dataset first")
          )
        ),

        # Obs -> Cells
        selectizeInput(
          ns("SI_obs_exp"), "Choose Experimental Observable", "",
          multiple = FALSE, options = list(placeholder = "choose dataset first")
        ),
        shinyjs::disabled(checkboxInput(ns("CB_aux_obs"), "use Additional Observables?", value = FALSE)),

        shinyjs::disabled(
          selectizeInput(ns("SI_aux_obs"), "Choose Additional Observables", "",
            multiple = TRUE, options = list(placeholder = "choose dataset first")
          )
        ),
        shinyjs::disabled(
          actionButton(ns("AB_ingest_load"), label = "(Re) load !")
        ),
        # shinyjs::disabled(
        #   selectizeInput(ns("SI_datatype"), "Data Type",
        #     choices = dataset_type, select = "X-"
        #   )
        # ),
        textOutput(ns("ui_datatype"))

      )
    )
  )
}

# ############################ +
# ## UI for interactive upload/choosing...
# ## i.e.  if "novel"
# TODO:  make a module for this..
# ############################ +


# # LOAD VAR
# selectizeInput(
#   ns("SI_load_var_names"),"Choose Primary Variable", "",
#   multiple=FALSE, options = list(placeholder = "choose dataset first")
# ),
#
# selectizeInput(
#   ns("SI_load_var_annots"),"Choose Variable Annotations", "",
#   multiple=TRUE, options = list(placeholder = "choose dataset first")
# ),
#
# # LOAD OBS
# selectizeInput(
#   ns("SI_load_obs_names"),"Choose Primary Observable", "",
#   multiple=FALSE, options = list(placeholder = "choose dataset first")
# ),
#
# selectizeInput(
#   ns("SI_load_obs_annots"),"Choose Observable Annotations", "",
#   multiple=TRUE, options = list(placeholder = "choose dataset first")
# ),
#
#
# # LOAD data-matrix
# if (input$data_type == 'transcriptomics'){
#   #LOAD RAW DATA TABLE
# }
#



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
      X = NULL,
      var = NULL,
      obs = NULL,
      database_name = NULL,
      omics_type = NULL, # Transcript-omics, Prote-omics, Lipid-omics, Metabol-omics, misc-

      var_names = NULL, # eg.  var_names: Genes, Proteins, Lipids
      var_annotations = NULL, # colnames of variable annotations.  e.g. gene/protein-families, ontologies, lipid class

      obs_names = NULL, # name of observations e.g. cell ID
      obs_annotations = NULL,
      # do i need this here?
      uns = list(), # NULL if not used, otherwise list of unstructured annotations...
      uns_keys = NULL,
      varm = list(), # NULL if not used, otherwise list of unstructured annotations...
      varm_keys = NULL,
      obsm = list(), # NULL if not used, otherwise list of unstructured annotations...
      obsm_keys = NULL,

      # these reflect the choices made in the ingestor
      # omics key feature i.e. genes, proteins, lipids
      #  primary & auxilarry
      #
      omics_feature = NULL,
      aux_features = NULL,
      # observables:  experiemntal factor... i.e. patient/control, old/young,
      #   plust aucuilarry observations. e.g. QC- batch, sex, etc  or inferred: cluster, label, cell-type
      exp_factor = NULL,
      aux_factors = NULL, # e.g.

      trigger = 0
    )

    ############################ +
    ## Update selectInput according to dataset
    ############################ +
    # observe({
    #   if (!is.null(input$SI_dataset)) {
    # load the data every time this changes...
    observeEvent(input$SI_dataset, {

      if (!is.null(input$SI_dataset)) { # unnesscasary defensive?

        ds_name <- (input$SI_dataset)

        if (ds_name == "DomenicoA") { # Skeletal Muscle

          X <- omicser::Domenico_A_X
          obs <- omicser::Domenico_A_obs
          var_ <- omicser::Domenico_A_var
          obsm <- omicser::Domenico_A_obsm
          varm <- omicser::Domenico_A_varm
          uns <- omicser::Domenico_A_uns

          to_return$database_name <- ds_name
          to_return$omics_type <- "prote"

          to_return$var <- var_
          to_return$obs <- obs
          to_return$X <- X

          to_return$var_names <-row.names(var_)
          to_return$var_annotations <- colnames(var_)

          to_return$obs_names <- row.names(obs)
          to_return$obs_annotations <- colnames(obs)

          to_return$uns <- uns
          to_return$varm <- varm
          to_return$obsm <- obsm

          # set a default to avoid funny business in side_select
          to_return$omics_feature <- to_return$var_annotations[1]

        } else if (ds_name == "VilasA") { # Vilas transc

          X <- omicser::Vilas_A_X
          obs <- omicser::Vilas_A_obs
          var_ <- omicser::Vilas_A_var
          obsm <- omicser::Vilas_A_obsm
          varm <- omicser::Vilas_A_varm
          uns <- list()

          to_return$database_name <- ds_name
          to_return$omics_type <- "transcript"

          to_return$var <- var_
          to_return$obs <- obs
          to_return$X <- X

          to_return$var_names <-row.names(var_) # the variable IDs
          to_return$var_annotations <- colnames(var_)

          to_return$obs_names <- row.names(obs)  # the sample IDs
          to_return$obs_annotations <- colnames(obs)

          to_return$uns <- uns
          to_return$varm <- varm
          to_return$obsm <- obsm

          # set a default to avoid funny business in side_select
          to_return$omics_feature <- to_return$var_annotations[1]

          # } else if (ds_name == "DIA") { # Muscle Stem Cell
          #
          #   obs <- omicser::domenico_DIA_obs
          #   #vars <- omicser::domenico_DIA_vars #NULL
          #   # X <- omicser::domenico_DIA_X  # NULL...don't load it if we don't want to use it
          #   vars <- obs$gene_name
          #   # force var_names to be "gene_names"... should this be reactive??
          #   to_return$var_names <- reactive({
          #     unique(obs$gene_name)
          #   })
          #   to_return$database_name <- ds_name
          #
          #   to_return$omics_type <- "prote"
          #
          #   to_return$var <- vars
          #   to_return$obs <- obs
          #   to_return$X <- NULL
          #   to_return$uns_meta <- NULL
          #
          #
          #   var_choices <- unique(vars)
          #   obs_choices <- colnames(obs)
          #
          # } else if (ds_name == "TMT") { # Skeletal Muscle
          #
          #   obs <- omicser::domenico_TMT_obs # should this be reactive?
          #
          #
          #   to_return$var_names <- reactive({
          #     unique(obs$gene_name)
          #   })
          #   to_return$database_name <- ds_name
          #   to_return$omics_type <- "prote"
          #
          #   to_return$vars <- vars
          #   to_return$obs <- obs
          #   to_return$X <- NULL
          #
          #   to_return$var_annotations <- colnames(vars)
          #   to_return$obs_names <- colnames(obs)
          #   to_return$obs_annotations <- colnames(obs)
          #   to_return$uns_meta <- NULL
          #
          #
          #
        #   # } else if (ds_name == "VilasB") { #Vilas Alt
        # } else if (ds_name == "OscarA") { # from Oscar
        #   X <- omicser::oscar_microglia1_X # should this be reactive?
        #   obs <- omicser::oscar_microglia1_obs # should this be reactive?
        #   vars <- omicser::oscar_microglia1_vars # should this be reactive?
        #
        #   data_table <- obs
        #
        #   # need to add some summary / differential columns to obs.
        #
        #   to_return$var_names <- reactive({
        #     unique(obs$gene_name)
        #   })
        #   to_return$database_name <- ds_name
        #   to_return$omics_type <- "transcript"
        #
        #   to_return$vars <- vars
        #   to_return$obs <- obs
        #   to_return$X <- NULL
        #
        #   to_return$var_annotations <- colnames(vars)
        #   to_return$obs_names <- colnames(obs)
        #   to_return$obs_annotations <- colnames(obs)
        #   to_return$uns_meta <- NULL
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
        #TODO: make var_names and aux_vars (and obs_exp & aux_obs) mutually exclusive
        var_choices <- isolate(to_return$var_annotations)
        updateSelectizeInput(session, "SI_var_name", choices = var_choices, selected = var_choices[1], server = TRUE)
        # TODO: update so we cant choose both at the same time
        #
        obs_choices <- isolate(to_return$obs_annotations)
        updateSelectizeInput(session, "SI_obs_exp", choices = obs_choices, selected = obs_choices[1], server = TRUE)

        shinyjs::enable("CB_aux_obs")
        shinyjs::enable("CB_aux_vars")

      } else {
        # disable check boxes
        shinyjs::disable("CB_aux_obs")
        shinyjs::disable("CB_aux_vars")
        print(" no database loaded... try Vilas Trans or Domenico ...")
      }
    })

    # keep these up to date for side_select... can maybe remove teh defaults
    observe({
      to_return$database_name <- input$SI_dataset
      to_return$omics_feature <- input$SI_var_name
    })


    # Update selectInput according to dataset
    observe({
      if (input$CB_aux_vars) {
        print("enabled aux variables ")
        shinyjs::enable("SI_aux_vars")
        var_choices <- isolate(to_return$var_annotations)
        updateSelectizeInput(session, "SI_aux_vars", choices = var_choices, selected = var_choices[1], server = TRUE)

      } else {
        print("disabled  aux variables ")
        shinyjs::disable("SI_aux_vars")
      }
    })

    observe({
      if (input$CB_aux_obs) {
        print("enabled aux observations ")
        shinyjs::enable("SI_aux_obs")
        obs_choices <- isolate(to_return$obs_annotations)
        updateSelectizeInput(session, "SI_aux_obs", choices = obs_choices, selected = obs_choices[2], server = TRUE)

      } else {
        print("disabled  aux observations ")
        shinyjs::disable("SI_aux_obs")
      }
    })

    # Update selectInput SI_var_name
    observe({
      if (input$SI_var_name != "") { #|| ( !is.na(input$SI_var0) )
        shinyjs::enable("AB_ingest_load")
      } else {
        print(' SI_var_name == "" ')
        shinyjs::disable("AB_ingest_load")
      }
    })


    # observe({
    #   if (input$SI_obs_exp == "") { #|| !(is.na(input$SI_var1))
    #     to_return$factors1 <- unique(factor(to_return$data_table[[input$SI_obs_exp]]))
    #
    #   } else {
    #     print("skipped SI_obs_exp")
    #   }
    # })
    # show the factors that have been loaded
    output$ui_datatype <- renderText({
      print(paste0("db type - ", to_return$omics_type , "-omics"))
    })


    # (Re)load button :: send the reactive object back to the app...
    observeEvent(input$AB_ingest_load, {

      if (input$CB_aux_vars) {
        to_return$aux_features <- input$SI_aux_vars
      } else {
        to_return$aux_features <- NULL
      }

      to_return$exp_factor <- input$SI_obs_exp
      if (input$CB_aux_obs) {
        to_return$aux_factors <- input$SI_aux_obs
      } else {
        to_return$aux_factors <- NULL
      }

      to_return$trigger <- to_return$trigger + 1
    })

     # NULL if not used, otherwise list of unstructured annotations...

    return(to_return)
  })
}

## To be copied in the UI
# mod_ingestor_ui("ingestor_ui_1")

## To be copied in the server
# mod_ingestor_server("ingestor_ui_1")

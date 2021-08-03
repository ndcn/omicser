# hardwire for now

# TODO: pack this into a .rda .rds to load for dynamic updates
dataset_names <- c(
  "Skeletal Muscle" = "TMT",
  "Muscle Stem Cell" = "DIA",
  "Domenico Tsc Long" = "DomenicoA",
  "Vilas Transcript" = "VilasA",
  "Vilas Microglia" = "VilasB",
  "Oscar Toy" = "OscarA",
  "Oscar Microglia" = "OscarB",
  "Upload New data" = "novel",
  "Yassene Lipidomics" ="YasseneA"
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
require(anndata)

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
      col_6( # full width
        selectizeInput(
          ns("SI_dataset"), "Dataset",
          choices = dataset_names,
          select = "VilasB"  # should i have a default??
          )
      )),
    fluidRow(
      col_6(
        # Var -> Genes/ Proteins
        selectizeInput(
          ns("SI_var_name"), "Choose -Omic Variable", "",
          multiple = FALSE, options = list(placeholder = "choose dataset first")
        )),
      col_4(
        ofset = 1,
        actionButton(ns("AB_subset_tog"), "Toggle subset observations"), #TODO: change text when toggled
        shinyjs::disabled(selectizeInput(ns("SI_subset"), "Obs information to subset:", "",
                       multiple = FALSE, options = list(placeholder = "choose dataset first"))),

        uiOutput(ns("ui_subset")),
        shinyjs::disabled(actionButton(ns("CB_sub_all"), "Select all groups", class = "btn btn-primary")),
        shinyjs::disabled(actionButton(ns("CB_sub_none"), "Deselect all groups", class = "btn btn-primary") )
      )
    ),
    fluidRow(
      col_4(
        shinyjs::disabled(checkboxInput(ns("CB_aux_groups"), "Omic Groupings?", value = FALSE))),
      col_4(
        shinyjs::disabled(checkboxInput(ns("CB_aux_vars"), "Aux Variable Annotations?", value = FALSE))),
    ),
    fluidRow(
      col_4(
        shinyjs::disabled(selectizeInput(ns("SI_aux_groups"), "Choose omic groupings", "",
                                  multiple = TRUE, options = list(placeholder = "choose dataset first"))
                          )),
      col_4(
        shinyjs::disabled(
          selectizeInput(ns("SI_aux_vars"), "Choose Variable Annotations", "",
            multiple = TRUE, options = list(placeholder = "choose dataset first")
          )
        ))

        ),
    fluidRow(
        col_8(

        # Obs -> Cells
        selectizeInput(
          ns("SI_obs_exp"), "Choose Experimental Variable", "",
          multiple = FALSE, options = list(placeholder = "choose dataset first")
          )
        )
    ),
    fluidRow(
      col_4(
        shinyjs::disabled( checkboxInput(ns("CB_obs_raw"), "Raw measures?",value=FALSE))
      ),
      col_4(
        shinyjs::disabled( checkboxInput(ns("CB_obs_comp"), "comparative measure?",value= FALSE))
      ),
      col_4(
        shinyjs::disabled( checkboxInput(ns("CB_obs_group"), "grouping var?",value=FALSE))
        )
    ),
    fluidRow(
      col_4(
        shinyjs::disabled(
          selectizeInput(ns("SI_obs_raw"), "Raw", "",
            multiple = FALSE, options = list(placeholder = "")
          )
        )
      ),
      col_4(
        shinyjs::disabled(
          selectizeInput(ns("SI_obs_comp"), "comparative", "",
                         multiple = FALSE, options = list(placeholder = "")
          )
        )),
      col_4(
        shinyjs::disabled(
          selectizeInput(ns("SI_obs_group"), "Grouping", "",
                         multiple = FALSE, options = list(placeholder = "")
          )
        )
        )
      ),
    fluidRow(
        # col_3(
          shinyjs::disabled(
            actionButton(ns("AB_ingest_load"), label = "(Re) load !")
          ),
          # ),
        # col_1(),
        # col_6(
        # shinyjs::disabled(
        #   selectizeInput(ns("SI_datatype"), "Data Type",
        #     choices = dataset_type, select = "X-"
        #   )
        # ),
        textOutput(ns("ui_datatype"))
        # )
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
                #
                # X = NULL,
                # var = NULL,
                # obs = NULL,
                # var_names = NULL, # eg.  var_names: Genes, Proteins, Lipids
                # var_annotations = NULL, # colnames of variable annotations.  e.g. gene/protein-families, ontologies, lipid class
                #
                # obs_names = NULL, # name of observations e.g. cell ID
                # obs_annotations = NULL,
                # # do i need this here?
                # uns = list(), # NULL if not used, otherwise list of unstructured annotations...
                # uns_keys = NULL,
                # varm = list(), # NULL if not used, otherwise list of unstructured annotations...
                # varm_keys = NULL,
                # obsm = list(), # NULL if not used, otherwise list of unstructured annotations...
                # obsm_keys = NULL,
      database_name = NULL,
      omics_type = NULL, # Transcript-omics, Prote-omics, Lipid-omics, Metabol-omics, misc-

      # these reflect the choices made in the ingestor
      # omics key feature i.e. genes, proteins, lipids
      #  primary & auxilarry

      #  Everything is packed into an anndata object
      ad = NULL,
      omics_feature = NULL, #the omics columnname...
      aux_features = NULL,
      # observables:  experiemntal factor... i.e. patient/control, old/young,
      #   plust aucuilarry observations. e.g. QC- batch, sex, etc  or inferred: cluster, label, cell-type
      exp_factor = NULL,
      aux_group = NULL,
      aux_comp = NULL,
      aux_raw = NULL,

      defs = NULL,
      confs = NULL,
      meta = NULL,
      trigger = 0
    )

    omic_rv <- reactiveValues(
      config = NULL,
      default = NULL)

    ############################ +
    ## Update selectInput according to dataset
    ############################ +
    # observe({
    #   if (!is.null(input$SI_dataset)) {
    # load the data every time this changes...
    observeEvent(input$SI_dataset, {
      if (!is.null(input$SI_dataset)) { # unnesscasary defensive?

        ds_name <- (input$SI_dataset)

        # if (ds_name == "DomenicoA") { # Skeletal Muscle
        #
        #   X <- omicser::Domenico_A_X
        #   obs <- omicser::Domenico_A_obs
        #   var_ <- omicser::Domenico_A_var
        #   obsm <- omicser::Domenico_A_obsm
        #   varm <- omicser::Domenico_A_varm
        #   uns <- omicser::Domenico_A_uns
        #
        #   to_return$database_name <- ds_name
        #   to_return$omics_type <- "prote"
        #
        #   to_return$var <- var_
        #   to_return$obs <- obs
        #   to_return$X <- X
        #
        #   to_return$var_names <-row.names(var_)
        #   to_return$var_annotations <- colnames(var_)
        #
        #   to_return$obs_names <- row.names(obs)
        #   to_return$obs_annotations <- colnames(obs)
        #
        #   to_return$uns <- uns
        #   to_return$varm <- varm
        #   to_return$obsm <- obsm
        #
        #   # set a default to avoid funny business in side_select
        #   to_return$omics_feature <- to_return$var_annotations[1]

        # } else
          if (ds_name == "VilasA") { # Vilas transc

          # X <- omicser::Vilas_A_X
          # obs <- omicser::Vilas_A_obs
          # var_ <- omicser::Vilas_A_var
          # obsm <- omicser::Vilas_A_obsm
          # varm <- omicser::Vilas_A_varm
          # uns <- list()


          ad <- anndata::read_h5ad(filename="data-raw/vilas_A.h5ad")


          omicconf <- omicser::vilas_A_conf
          omicdef <- omicser::vilas_A_def
          omics <- omicser::vilas_A_omics
          omicmeta <- omicser::vilas_A_meta


          to_return$database_name <- ds_name
          to_return$omics_type <- "transcript"

          to_return$ad <- ad
          # set a default to avoid funny business in side_select
          to_return$omics_feature <- omics
          to_return$meta <- omicmeta  # this might be too redundant
          # to_return$defaults <- omicdef
          # to_return$configs <- omicconf

          omicconf$meta <- as.data.table(omicconf$meta)
          omicconf$mat <- as.data.table(omicconf$mat)
          omic_rv$config <- omicconf

          omic_rv$default <- omicdef

          } else if (ds_name == "VilasB") { # Vilas microglia


            ad <- anndata::read_h5ad(filename="data-raw/vilas_B.h5ad")


            omicconf <- omicser::vilas_B_conf
            omicdef <- omicser::vilas_B_def
            omics <- omicser::vilas_B_omics
            omicmeta <- omicser::vilas_B_meta

            to_return$database_name <- ds_name
            to_return$omics_type <- "transcript"

            to_return$ad <- ad
            # set a default to avoid funny business in side_select
            to_return$omics_feature <- omics
            to_return$meta <- omicmeta  # this might be too redundant
            # to_return$defaults <- omicdef
            # to_return$configs <- omicconf

            omic_rv$config <- omicconf

            omic_rv$default <- omicdef
            # raw_obs <-
            # group_obs <-
            # comp_obs <-
          } else if (ds_name == "YasseneA") { # Vilas transc


            # X <- omicser::yassene_A_X
            # obs <- omicser::yassene_A_obs
            # var_ <- omicser::yassene_A_var
            # obsm <- omicser::yassene_A_obsm
            # varm <- omicser::yassene_A_varm
            # uns <- omicser::yassene_A_varm
            # layers <- omicser::yassene_A_layers

            ad <- anndata::read_h5ad(filename="data-raw/yassene_A.h5ad")


            omicconf <- omicser::yassene_A_conf
            omicdef <- omicser::yassene_A_def
            omics <- omicser::yassene_A_omics
            omicmeta <- omicser::yassene_A_meta

            to_return$database_name <- ds_name
            to_return$omics_type <- "lipid"

            to_return$ad <- ad
            # set a default to avoid funny business in side_select
            to_return$omics_feature <- omics
            to_return$meta <- omicmeta  # this might be too redundant
            # to_return$defaults <- omicdef
            # to_return$configs <- omicconf

            omic_rv$default <- omicdef
            # raw_obs <-
            # group_obs <-
            # comp_obs <-
          } else if (ds_name == "DomenicoA") { # "Domenico Tsc Long"


            # X <- omicser::domenico_A_X
            # obs <- omicser::domenico_A_obs
            # var_ <- omicser::domenico_A_var
            # obsm <- omicser::domenico_A_obsm
            # varm <- omicser::domenico_A_varm
            # uns <- omicser::domenico_A_varm
            # layers <- omicser::domenico_A_layers

            ad <- anndata::read_h5ad(filename="data-raw/domenico_A.h5ad")
            #ad <- anndata::read_h5ad(filename="data-raw/domenico_A.h5ad",backed = 'r')
            # If 'r', load ~anndata.AnnData in backed mode instead of fully
            # loading it into memory (memory mode). If you want to modify
            # backed attributes of the AnnData object, you need to choose 'r+'.

            omicconf <- omicser::domenico_A_conf
            omicdef <- omicser::domenico_A_def
            omics <- omicser::domenico_A_omics
            omicmeta <- omicser::domenico_A_meta

            to_return$database_name <- ds_name
            to_return$omics_type <- "prote"

            to_return$ad <- ad
            # set a default to avoid funny business in side_select
            to_return$omics_feature <- omics
            to_return$meta <- omicmeta  # this might be too redundant
            # to_return$defaults <- omicdef
            # to_return$configs <- omicconf


            omic_rv$config <- omicconf

            omic_rv$default <- omicdef
            # raw_obs <-
            # group_obs <-
            # comp_obs <-


          # } else if (ds_name == "OscarB") { # Vilas transc
          #
          #   X <- omicser::oscar_microglia_X
          #   obs <- omicser::oscar_microglia_obs
          #   var_ <- omicser::oscar_microglia_var
          #   obsm <- omicser::oscar_microglia_obsm
          #   varm <- omicser::oscar_microglia_varm
          #   uns <- list()
          #
          #   to_return$database_name <- ds_name
          #   to_return$omics_type <- "transcript"
          #
          #   to_return$var <- var_
          #   to_return$obs <- obs
          #   to_return$X <- X
          #
          #   to_return$var_names <-row.names(var_) # the variable IDs
          #   to_return$var_annotations <- colnames(var_)
          #
          #   to_return$obs_names <- row.names(obs)  # the sample IDs
          #   to_return$obs_annotations <- colnames(obs)
          #
          #   to_return$uns <- uns
          #   to_return$varm <- varm
          #   to_return$obsm <- obsm
          #
          #   # set a default to avoid funny business in side_select
          #   to_return$omics_feature <- to_return$var_annotations[1]
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
        #TODO: make var_names and aux_vars (and obs_exp & aux_obs) mutually exclusive
        var_choices <- isolate(to_return$ad$var_keys())
        var_choices <- var_choices[1]  # short circuit for now...
        updateSelectizeInput(session, "SI_var_name", choices = var_choices, selected = var_choices[1], server = TRUE)
        # TODO: update so we cant choose both at the same time
        #
        obs_choices <- isolate(to_return$ad$obs_keys())
        updateSelectizeInput(session, "SI_obs_exp", choices = obs_choices, selected = obs_choices[1], server = TRUE)
        print(paste0("experimental observations: ",obs_choices))
        #shinyjs::enable("CB_aux_obs")

        shinyjs::enable("CB_aux_vars")
        shinyjs::enable("CB_obs_raw")
        shinyjs::enable("CB_obs_comp")
        shinyjs::enable("CB_obs_group")

      } else {
        # disable check boxes
        # shinyjs::disable("CB_aux_obs")
        # shinyjs::disable("CB_aux_vars")

        shinyjs::disable("CB_aux_vars")
        shinyjs::disable("CB_obs_raw")
        shinyjs::disable("CB_obs_comp")
        shinyjs::disable("CB_obs_group")

        print(" no database loaded... try Vilas Trans or Domenico or Oscar Microglia...")
      }
    })


    # Update selectInput according to dataset
    observe({
      req(omic_rv$config)
      browser()
      if (input$AB_subset_tog) {
        print("enabled subset ")
        shinyjs::enable("SI_subset")
        updateSelectizeInput(session, "SI_subset","Obs information to subset:",
                             choices = omic_rv$config$meta[grp == TRUE]$UI,
                             selected = omic_rv$default$grp1,  server = TRUE)

        shinyjs::enable("CB_sub_all")
        shinyjs::enable("CB_sub_none")
      } else {
        print("disabled subset ")
        shinyjs::disable("SI_subset")
        shinyjs::disable("CB_sub_all")
        shinyjs::disable("CB_sub_none")
      }
    })

    output$ui_subset <- renderUI({
      sub = strsplit(omic_rv$config$meta[UI == input$SI_subset]$fID, "\\|")[[1]]
      checkboxGroupInput("CB_sub_inner1", "Select which groups to show", inline = TRUE,
                         choices = sub, selected = sub)
    })
    observeEvent(input$CB_sub_all, {
      sub = strsplit(omic_rv$config$meta[UI == input$SI_subset]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(session, inputId = "CB_sub_inner1", label = "Select which groups to show",
                               choices = sub, selected = NULL, inline = TRUE) # WARNING.  make sure subset is null is checked (length(0?))
    })
    observeEvent(input$CB_sub_none, {
      sub = strsplit(omic_rv$config$meta[UI == input$SI_subset]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(session, inputId = "CB_sub_inner1", label = "Select which groups to show",
                               choices = sub, selected = sub, inline = TRUE)
    })

    # keep these up to date for side_select... can maybe remove teh defaults
    observe({
      to_return$database_name <- input$SI_dataset
      #to_return$omics_feature <- input$SI_var_name  #send the loaded value...
    })


    # Update selectInput according to dataset
    observe({

      if (input$CB_aux_vars) {
        print("enabled aux variables ")
        shinyjs::enable("SI_aux_vars")
        var_choices <- isolate(to_return$ad$var_keys())
        updateSelectizeInput(session, "SI_aux_vars", choices = var_choices, selected = var_choices[1], server = TRUE)
        print(var_choices[1:max(10,length(var_choices))])

      } else {
        print("disabled  aux variables ")
        shinyjs::disable("SI_aux_vars")
      }
    })

    observe({
      if (input$CB_obs_raw) {
        print("enabled raw observations ")
        obs_choices <- isolate(to_return$ad$varm_keys())
        pattern_raw <- "^(frac_|mean_).+"
        obs_choices <- obs_choices[grepl(x=obs_choices,pattern=pattern_raw)]
        if (length(obs_choices[grepl(x=obs_choices,pattern=pattern_raw)])>0) {
          shinyjs::enable("SI_obs_raw")
          updateSelectizeInput(session, "SI_obs_raw", choices = obs_choices, selected = obs_choices[2], server = TRUE)
          print(obs_choices)
        } else {
          print("no raw values to consider ")
        }

      } else {
        print("disabled  raw observations ")
        shinyjs::disable("SI_obs_raw")
        # change back to placeholder??
        #updateSelectizeInput(session, "SI_obs_raw", choices = obs_choices, selected = obs_choices[2], server = TRUE)
      }
    })


    observe({
      if (input$CB_obs_comp) {
        print("enable comparative observations ")

        obs_choices <- isolate(to_return$ad$varm_keys())
        pattern_raw <- "^.*(Ratio|FC|Q).+"
        obs_choices <- obs_choices[grepl(x=obs_choices,pattern=pattern_raw)]

        if (length(obs_choices[grepl(x=obs_choices,pattern=pattern_raw)])>0) {
          print(obs_choices)
          shinyjs::enable("SI_obs_comp")
          updateSelectizeInput(session, "SI_obs_comp", choices = obs_choices, selected = obs_choices[1], server = TRUE)
        } else {
          print("no comparitives to consider ")
        }
      } else {
        print("disabled comparative observations ")
        shinyjs::disable("SI_obs_comp")
        # change back to placeholder??
        #updateSelectizeInput(session, "SI_obs_raw", choices = obs_choices, selected = obs_choices[2], server = TRUE)
      }
    })

    observe({
      if (input$CB_obs_group) {
        print("enable grouping observations ")

        shinyjs::enable("SI_obs_group")
        obs_choices <- isolate(colnames(to_return$obs))
        pattern_raw <- "^.*(cluster|batch|oup|grp|ond|abel|lbl).+"
        obs_choices <- obs_choices[grepl(x=obs_choices,pattern=pattern_raw)]

        if (length(obs_choices[grepl(x=obs_choices,pattern=pattern_raw)])>0) {
          print(obs_choices)
          shinyjs::enable("SI_obs_group")
          updateSelectizeInput(session, "SI_obs_group", choices = obs_choices, selected = obs_choices[1], server = TRUE)
        } else {
          print("no groupings to consider ")
        }

      } else {
        print("disabled grouping observations ")
        shinyjs::disable("SI_obs_group")
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


      if (input$CB_obs_raw) {
        to_return$aux_raw <- input$SI_obs_raw
      } else {
        to_return$aux_raw <- NULL
      }

      if (input$CB_obs_comp) {
        to_return$aux_comp <- input$SI_obs_comp
      } else {
        to_return$aux_comp <- NULL
      }

      if (input$CB_obs_group) {
        to_return$aux_group <- input$SI_obs_group
      } else {
        to_return$aux_group <- NULL
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

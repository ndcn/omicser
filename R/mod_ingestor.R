# hardwire for now

# TODO: pack this into a .rda .rds to load for dynamic updates
dataset_names <- c(
  "Domenico Tsc Long" = "DomenicoA",
  "Vilas Transcript" = "VilasA",
  "Vilas Microglia" = "VilasB",
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
      col_3(
        textOutput(ns("ui_datatype"))
      )
    ),
    fluidRow(
      col_6(
        # Var -> Genes/ Proteins
        uiOutput(ns("ui_omics"))

        # selectizeInput(
        #   ns("SI_var_name"), "Choose -Omic Variable", "",
        #   multiple = FALSE, options = list(placeholder = "choose dataset first")
        # )
        ),
    #   col_4(
    #     offset = 1,
    #     actionButton(ns("AB_subset_tog"), "Toggle subset observations"), #TODO: change text when toggled
    #     shinyjs::disabled(selectizeInput(ns("SI_subset"), "Obs information to subset:", "",
    #                    multiple = FALSE, options = list(placeholder = "choose dataset first"))),
    #
    #     uiOutput(ns("ui_subset")),
    #     shinyjs::disabled(actionButton(ns("CB_sub_all"), "Select all groups", class = "btn btn-primary")),
    #     shinyjs::disabled(actionButton(ns("CB_sub_none"), "Deselect all groups", class = "btn btn-primary") )
    #   )
    ),
    # fluidRow(
    #   col_4(
    #     shinyjs::disabled(checkboxInput(ns("CB_aux_groups"), "Omic Groupings?", value = FALSE))
    #     )
    #   # ,
    #   # col_4(offset = 1,
    #   #   shinyjs::disabled(checkboxInput(ns("CB_aux_vars"), "Aux Variable Annotations?", value = FALSE))),
    # ),
    # fluidRow(
    #   col_4(
    #     shinyjs::disabled(selectizeInput(ns("SI_aux_groups"), "Choose omic groupings", "",
    #                               multiple = TRUE, options = list(placeholder = "choose dataset first"))
    #                       )),
    #   col_4(
    #     shinyjs::disabled(
    #       selectizeInput(ns("SI_aux_vars"), "Choose Variable Annotations", "",
    #         multiple = TRUE, options = list(placeholder = "choose dataset first")
    #       )
    #     ))
    #
    #     ),
    fluidRow(
      helpText(HTML("<i>observations</i>  <b> experimental</b>.")
      )
      ),
    fluidRow(
        col_8(
          textOutput(ns("ui_obs_exp"))

        # # Obs -> Cells
        # selectizeInput(
        #   ns("SI_obs_exp"), "Choose Experimental Variable", "",
        #   multiple = TRUE, options = list(placeholder = "choose dataset first")
        #   )
        )
    ),
    # fluidRow(
    #   col_4(
    #     shinyjs::disabled( checkboxInput(ns("CB_obs_raw"), "Raw measures?",value=FALSE))
    #   ),
    #   col_4(
    #     shinyjs::disabled( checkboxInput(ns("CB_obs_comp"), "comparative measure?",value= FALSE))
    #   )
    #   # ,
    #   # col_4(
    #   #   shinyjs::disabled( checkboxInput(ns("CB_obs_group"), "grouping var?",value=FALSE))
    #   #   )
    # ),
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
        ))
      ),
    fluidRow(
      col_3(
        shinyjs::disabled(
          actionButton(ns("AB_ingest_load"), label = "(Re) load !")
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


      # these reflect the choices made in the ingestor
      # omics key feature i.e. genes, proteins, lipids
      #  primary & auxilarry
      omics_feature = NULL, #the omics columnname...
      aux_features = NULL,
      # observables:  experiemntal factor... i.e. patient/control, old/young,
      #   plust aucuilarry observations. e.g. QC- batch, sex, etc  or inferred: cluster, label, cell-type
      exp_factor = NULL,
      aux_group = NULL,
      aux_comp = NULL,
      aux_raw = NULL,

      config = NULL,
      default = NULL,
      meta = NULL,
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
          if (ds_name == "VilasA") { # Vilas transc

            ad <- anndata::read_h5ad(filename="data-raw/vilas_A.h5ad")

            omicconf <- omicser::vilas_A_conf
            omicdef <- omicser::vilas_A_def
            omics <- omicser::vilas_A_omics
            omicmeta <- omicser::vilas_A_meta

            to_return$database_name <- ds_name
            to_return$omics_type <- "transcript"
            to_return$ad <- ad
            to_return$omics_feature <- omics
            to_return$meta <- omicmeta  # this might be too redundant

            omicconf$meta <- as.data.table(omicconf$meta)
            omicconf$mat <- as.data.table(omicconf$mat)
            to_return$config <- omicconf
            to_return$default <- omicdef

          } else if (ds_name == "VilasB") { # Vilas microglia

            ad <- anndata::read_h5ad(filename="data-raw/vilas_B.h5ad")

            omicconf <- omicser::vilas_B_conf
            omicdef <- omicser::vilas_B_def
            omics <- omicser::vilas_B_omics
            omicsomicmeta <- omicser::vilas_B_meta

            to_return$database_name <- ds_name
            to_return$omics_type <- "transcript"

            to_return$ad <- ad
            # set a default to avoid funny business in side_select
            to_return$omics_feature <- omics
            to_return$meta <- omicmeta  # this might be too redundant

            to_return$config <- omicconf
            to_return$default <- omicdef

          } else if (ds_name == "YasseneA") { # Vilas transc

            ad <- anndata::read_h5ad(filename="data-raw/yassene_A.h5ad")


            omicconf <- omicser::yassene_A_conf
            omicdef <- omicser::yassene_A_def
            omics <- omicser::yassene_A_omics
            omicmeta <- omicser::yassene_A_meta

            to_return$database_name <- ds_name
            to_return$omics_type <- "lipid"

            to_return$ad <- ad
            to_return$omics_feature <- omics
            to_return$meta <- omicmeta  # this might be too redundant

            to_return$config <- omicconf
            to_return$default <- omicdef
          } else if (ds_name == "DomenicoA") { # "Domenico Tsc Long"

            ad <- anndata::read_h5ad(filename="data-raw/domenico_A.h5ad")

            omicconf <- omicser::domenico_A_conf
            omicdef <- omicser::domenico_A_def
            omics <- omicser::domenico_A_omics
            omicmeta <- omicser::domenico_A_meta

            to_return$database_name <- ds_name
            to_return$omics_type <- "prote"

            to_return$ad <- ad
            to_return$omics_feature <- omics
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
        # var_choices <- isolate(to_return$ad$var_keys())
        # var_choices <- var_choices[1]  # short circuit for now...
        # updateSelectizeInput(session, "SI_var_name", choices = var_choices, selected = var_choices[1], server = TRUE)
        shinyjs::enable("SI_obs_raw")
        #shinyjs::enable("CB_aux_vars")

        #shinyjs::enable("CB_obs_raw")
        #shinyjs::enable("CB_obs_comp")
        #shinyjs::enable("CB_obs_group")

      } else {
        # disable check boxes
        shinyjs::disable("SI_obs_raw")
        #shinyjs::disable("CB_aux_vars")
        # shinyjs::disable("CB_obs_raw")
        # shinyjs::disable("CB_obs_comp")
        #shinyjs::disable("CB_obs_group")

        print(" no database loaded... try Vilas Trans or Domenico or Yassene.")
      }
    })
    output$ui_obs_exp <- renderPrint({
      obs_choices <- isolate(to_return$ad$obs_keys())

      if (is.null(obs_choices)) {
        print("(obs_exp) no datbase loaded")
      } else {
        print(paste0("observ exp: ", paste(obs_choices, collapse = ",")) )
      }
    })

    output$ui_omics <- renderPrint({
      if (is.null(to_return$omics_feature)) {
        print("no datbase loaded")
      } else {
        omics <- isolate(to_return$omics_feature)
        print(paste0("Current database is: ", paste( names(omics)[1:10], collapse = ",")) )
      }
    })

    # # Update selectInput according to dataset
    # observe({
    #   req(omic_rv$config)
    #
    #   if (input$AB_subset_tog) {
    #     print("enabled subset ")
    #     shinyjs::enable("SI_subset")
    #     updateSelectizeInput(session, "SI_subset","Obs information to subset:",
    #                          choices = omic_rv$config$meta[grp == TRUE]$UI,
    #                          selected = omic_rv$default$grp1,  server = TRUE)
    #
    #     shinyjs::enable("CB_sub_all")
    #     shinyjs::enable("CB_sub_none")
    #   } else {
    #     print("disabled subset ")
    #     shinyjs::disable("SI_subset")
    #     shinyjs::disable("CB_sub_all")
    #     shinyjs::disable("CB_sub_none")
    #   }
    # })
#
#     output$ui_subset <- renderUI({
#       sub = strsplit(omic_rv$config$meta[UI == input$SI_subset]$fID, "\\|")[[1]]
#       checkboxGroupInput("CB_sub_inner1", "Select which groups to show", inline = TRUE,
#                          choices = sub, selected = sub)
#     })
#     observeEvent(input$CB_sub_all, {
#       sub = strsplit(omic_rv$config$meta[UI == input$SI_subset]$fID, "\\|")[[1]]
#       updateCheckboxGroupInput(session, inputId = "CB_sub_inner1", label = "Select which groups to show",
#                                choices = sub, selected = NULL, inline = TRUE) # WARNING.  make sure subset is null is checked (length(0?))
#     })
#     observeEvent(input$CB_sub_none, {
#       sub = strsplit(omic_rv$config$meta[UI == input$SI_subset]$fID, "\\|")[[1]]
#       updateCheckboxGroupInput(session, inputId = "CB_sub_inner1", label = "Select which groups to show",
#                                choices = sub, selected = sub, inline = TRUE)
#     })

    # keep these up to date for side_select..
    observe({
      to_return$database_name <- input$SI_dataset
    })

    # show the factors that have been loaded
    output$ui_datatype <- renderText({
      print(paste0("db type - ", to_return$omics_type , "-omics"))
    })


    # # Update selectInput according to dataset
    # observe({
    #
    #   if (input$CB_aux_vars) {
    #     print("enabled aux variables ")
    #     shinyjs::enable("SI_aux_vars")
    #     var_choices <- isolate(to_return$ad$var_keys())
    #     updateSelectizeInput(session, "SI_aux_vars", choices = var_choices, selected = var_choices[1], server = TRUE)
    #     print(var_choices[1:max(10,length(var_choices))])
    #
    #   } else {
    #     print("disabled  aux variables ")
    #     shinyjs::disable("SI_aux_vars")
    #   }
    # })

    observe({
      req(to_return$config)
      raw_choices = to_return$config$mat[observ == TRUE & ID!="raw"]$fID
      if (length(raw_choices)>0) {
        updateSelectizeInput(session, "SI_obs_raw","raw measures:",
                                                    choices = raw_choices,
                                                    selected = raw_choices[1],  server = TRUE)
      } else {
        print("disabled  raw observations ")
        shinyjs::disable("SI_obs_raw")
        # change back to placeholder??
        updateSelectizeInput(session, "SI_obs_raw", "Raw ", "", options = list(placeholder = ""))
      }

    })

    observe({
      req(to_return$config)

      comp_choices = to_return$config$mat[comp == TRUE ]$fID
      if (length(comp_choices)>0) {
        updateSelectizeInput(session, "SI_obs_comp","comparatives:",
                             choices = comp_choices,
                             selected = comp_choices[1],  server = TRUE)
      } else {
        print("disabled  comparative observations ")
        shinyjs::disable("SI_obs_comp")
        updateSelectizeInput(session, "SI_obs_comp", "comparative", "",
                        options = list(placeholder = "") )
      }


    })



    # Update selectInput SI_obs_exp
    observe({
      req(input$SI_obs_raw)
      if (input$SI_obs_raw != "") { #|| ( !is.na(input$SI_var0) )
        shinyjs::enable("AB_ingest_load")
      } else {
        print(' SI_obs_exp == <empty> ')
        shinyjs::disable("AB_ingest_load")
      }
    })



    # (Re)load button :: send the reactive object back to the app...
    observeEvent(input$AB_ingest_load, {

      to_return$exp_factor <- input$SI_obs_exp
      # if (input$CB_aux_vars) {
      #   to_return$aux_features <- input$SI_aux_vars
      # } else {
      #   to_return$aux_features <- NULL
      # }
      to_return$aux_raw <- input$SI_obs_raw
      to_return$aux_comp <- input$SI_obs_comp

      # if (input$CB_obs_raw) {
      #   to_return$aux_raw <- input$SI_obs_raw
      # } else {
      #   to_return$aux_raw <- NULL
      # }
      #
      # if (input$CB_obs_comp) {
      #   to_return$aux_comp <- input$SI_obs_comp
      # } else {
      #   to_return$aux_comp <- NULL
      # }

      # if (input$CB_obs_group) {
      #   to_return$aux_group <- input$SI_obs_group
      # } else {
      #   to_return$aux_group <- NULL
      # }

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

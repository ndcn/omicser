
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
    column(width = 3,
           offset = 0,
           selectizeInput(ns("SI_database"),
                            "Choose Database",
                            "", multiple=FALSE, options = list(placeholder = " database first")
                            )
           ),
        column(width = 6,
               offset = 0,
               shinyjs::disabled(
               actionButton(ns("AB_ingest_load"), label = "Load Data", class = "btn btn-large btn-danger")
                      ),
               "loaded data |------> ",
               br(),
               uiOutput(ns("ui_data_meta")),
        )
      ),

    hr(style = "border-top: 1px dashed grey;"),
    fluidRow(


      column(width = 4,
             offset = 0,

             h4("Select the relavent meta-table columns:"),
             hr(style = "border-top: 1px dashed grey;"),

             h4("sample meta (obs)"),
             fluidRow(
               column(6,
                      selectizeInput(
                               ns("SI_exp_fact"), "experimental factors (grouping)", "",
                               multiple=TRUE, options = list(placeholder = "Choose sample factors")
                             )),
               column(6,
                      selectizeInput(
                           ns("SI_exp_annot"), "sample annotation (heatmap)", "",
                           multiple=TRUE, options = list(placeholder = "Choose sample-meta annotations (heatmap)")
                         )
               )
             ),
             hr(style = "border-top: 1px dashed grey;"),
             h4("feature meta (var)"),
             fluidRow(
               column(6,
               selectizeInput(
                 ns("SI_omic_feat"), "omic feature groups", "",
                 multiple=TRUE, options = list(placeholder = "Choose omic feature groups")
               )),
               column(6,
                 selectizeInput(
                 ns("SI_feat_annot"), "feature annotation (heatmap)", "",
                 multiple=TRUE, options = list(placeholder = "Choose omic feature annotatins (heatmap)")
               ))
              ),
               hr(style = "border-top: 1px dashed grey;"),

             fluidRow(
               column(6,
                      selectizeInput(
                             ns("SI_feat_deets"), "feature details", "",
                             multiple=TRUE, options = list(placeholder = "Choose omic feature meta-info")
                           )
               ),


               column(6,
                      HTML("<b>feature filtering variable</b>"),br(),
                      checkboxInput(ns("CB_fano_fact"), "calculate? (relative variance)", value = TRUE),

                    shinyjs::disabled(
                           selectizeInput(
                             ns("SI_feature_filter"), "filter variable", "",
                             multiple=FALSE, options = list(placeholder = "Choose omic feature for filtering `highly variable`")
                           ))
               )
             )
          ),

          column(width = 8,
                 offset = 0,
                 mod_additional_info_ui(id = ns("additional_info_ui_ingest"), title = "Databse Information")
          )
     ) #fluidrow
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


    db_name = reactiveValues(name=NULL,
                             dir=NULL)
    # observe({
    #   db_name$name <- to_return$db_meta$name
    # })

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
      #  should these just reference a top-level reactive??
      anndata = NULL,
      de = NULL,
      config = NULL,
      default = NULL,


      trigger = 0,

      shaddow_defs = NULL,
      fano_factor = NULL
    )

    updateSelectizeInput(session, "SI_database", choices = db_names, selected = db_names[1], server=TRUE)


    idx_of_disp <- reactive({
        matrixStats::colVars(raw,na.rm = TRUE)/colMeans(raw,na.rm = TRUE)
    })


  shaddow_defs <- reactiveValues(
    # holding the interactive selected factors and features which supercede what is from the defaults/configs
    expr_annot = NULL,
    exp_fact = NULL,
    exp_annot = NULL,
    omic_feat = NULL,
    feat_annot = NULL,
    feat_deets = NULL,
    feature_filter = NULL

   )

  temp_rv <- reactiveValues(
    # holding the interactive selected factors and features which supercede what is from the defaults/configs
    anndata = NULL,
    de = NULL,
    config = NULL,
    default = NULL,
    fano_factor = NULL

  )



    ## load dataset (observeEvent) ===================
    #observeEvent(input$SI_database, {
    observe({
        req(input$SI_database)

      db_name$dir <- input$SI_database
      db_name$name <- names(which(db_names==input$SI_database))
      # zero out the data_meta until we "load"
      output$ui_data_meta <- renderUI({h4("Press `LOAD` to activate data") })


      # load data
      temp_rv$de <- readRDS(file = file.path(db_root_path,db_name$dir,"db_de_table.rds"))
      anndata <- anndata::read_h5ad(filename=file.path(db_root_path,db_name$dir,"db_data.h5ad"))
      temp_rv$anndata <- anndata
      conf_def <- gen_config_table(anndata, db_name$dir, db_root_path)
      #update reactives...
      temp_rv$config <- conf_def$conf
      temp_rv$default <- conf_def$def

      # computed index of dispersion
      # TODO: fix this HACK
      # computer no-matter what,
      temp_rv$fano_factor <- (matrixStats::colVars(anndata$X,na.rm = TRUE)/(colMeans(anndata$X,na.rm = TRUE)+10^-40))

      obs_choices <- anndata$obs_keys()

      group_obs <- temp_rv$config[grp == TRUE & field=="obs"]$UI # <- choices_x
      def_grp_o <- temp_rv$default$obs_subset

      #will still be disabled if we are calculating...


      updateSelectInput(inputId = "SI_exp_fact", #label="Choose sample factors"
                        choices = obs_choices, #multiple=TRUE
                        selected = group_obs)



      updateSelectInput(inputId = "SI_exp_annot", #label="Choose sample-meta annotations (heatmap)"
                        choices = obs_choices, #multiple=TRUE
                        selected = group_obs)


      var_choices <- anndata$var_keys()
      group_var <- temp_rv$config[grp == TRUE & field=="var"]$UI # <- choices_x
      def_grp_v <- temp_rv$default$var_subset


      updateSelectInput(inputId = "SI_omic_feat", #label="Choose omic feature groups"
                        choices = var_choices, #multiple=TRUE
                        selected = def_grp_v)

      updateSelectInput(inputId = "SI_feat_annot", #label= "Choose omic feature annotatins (heatmap)"
                        choices = var_choices, #multiple=TRUE
                        selected = group_var)

      updateSelectInput(inputId = "SI_feat_deets", #label= "Choose omic feature meta-info"
                        choices = var_choices, #multiple=TRUE
                        selected = var_choices)

      updateSelectInput(inputId = "SI_feature_filter",# label="Choose omic feature for filtering `highly variable`"
                        choices = var_choices, #multiple=FALSE
                        selected = "")

    })



    # is there overhead to store these in a reactive conduit?

    ## render info (TODO:) ===================



    ## observes ===================

    observe({
      if ( !is.null( db_name$name ) ) {
        shinyjs::enable("AB_ingest_load")
      } else {
        shinyjs::disable("AB_ingest_load")
        print(" no database loaded... .")
      }
    })








    observe({
      if ( input$CB_fano_fact ) {
        shinyjs::disable("SI_feature_filter")
        # compute and set ? value??

      } else {
        shinyjs::enable("SI_feature_filter")
        message(" choose.. .")
      }
    })




    # keep these up to date for side_select..
    # # Should this just be done in     observeEvent(input$SI_database,?

    # load button :: send the reactive object over to the side selector.
    observeEvent(input$AB_ingest_load, {
      # all other return values set with SI_database
      # TODO: FIXE CONF_DEF TABLE with the new selections...
      to_return$de <- temp_rv$de
      to_return$anndata <- temp_rv$anndata
      to_return$config <- temp_rv$config
      to_return$default <- temp_rv$default


      to_return$db_meta$name<-db_name$name
      to_return$db_meta$db_dir<-db_name$dir

      to_return$db_meta$omics_type<-temp_rv$default$db_meta$omic_type
      to_return$db_meta$measurment<-temp_rv$default$db_meta$measurement
      to_return$db_meta$organism<-temp_rv$default$db_meta$organizm
      to_return$db_meta$publication<-temp_rv$default$db_meta$pub
      to_return$db_meta$etc<-temp_rv$default$db_meta$annotation_database

      to_return$trigger <- to_return$trigger + 1



      shaddow_defs$exp_fact <- input$SI_exp_fact
      shaddow_defs$exp_annot <- input$SI_exp_annot
      shaddow_defs$omic_feat <- input$SI_omic_feat
      shaddow_defs$feat_annot <- input$SI_feat_annot
      shaddow_defs$feat_deets  <- input$SI_feat_deets
      shaddow_defs$feature_filter  <- if (input$CB_fano_fact) "fano factor" else input$SI_feature_filter


      to_return$shaddow_defs <- shaddow_defs
      to_return$fano_factor <- temp_rv$fano_factor


      # will this work?  or is the "observing" affecting a reactive with this render?
      output$ui_data_meta <- renderUI({
        out_text <- HTML(paste("from the <i>", temp_rv$default$db_meta$lab),"</i> lab")
        ret_tags <-  tagList(
          h4(to_return$db_meta$omics_type),
          #"some more text",
          out_text,
          br()
        )
        return( ret_tags )
      })




    })




    return(to_return)
  })
}

## To be copied in the UI
# mod_ingestor_ui("ingestor_ui_1")

## To be copied in the server
# mod_ingestor_server("ingestor_ui_1")

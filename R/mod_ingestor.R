
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
      column(width = 2,
             offset = 0,
             shinyjs::disabled(
             actionButton(ns("AB_ingest_load"), label = "Load Data", class = "btn btn-large btn-danger")
                    ),
             br(),
             "   loaded data |------> "
             ),
      column(width = 3,
             offset = 0,
             uiOutput(ns("ui_data_meta"))
            )
    ),
    fluidRow(
      column(width = 3,
             offset = 0,
             actionButton(ns("AB_add_database_load"), label = "Add Database", class = "btn btn-large btn-danger")
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
                      checkboxInput(ns("CB_var_to_mean_ratio"), "calculate VMR? (idx if disp.)", value = TRUE),

                    shinyjs::disabled(
                           selectizeInput(
                             ns("SI_feature_filter"), "filter variable", "",
                             multiple=FALSE, options = list(placeholder = "Choose omic feature for filtering `highly variable`")
                           )),

                    br(),br(),
                    shinyjs::disabled(
                      actionButton(ns("AB_update_defaults"), label = "Update Defaults", class = "btn btn-large btn-danger")
                    )
               )
             )
          ),

          column(width = 8,
                 offset = 0,
                 mod_additional_info_ui(id = ns("additional_info_ui_ingest"), title = "Database Information")
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


    # # config is in the ./inst/db_info/ directory
    # CONFIG <- omicser::get_config()
    # DB_NAMES <- CONFIG$database_names
    # DB_ROOT <- CONFIG$db_root_path
    # INSTALL_TYPE <- CONFIG$install #production, dev,


    DB_NAMES <- golem::get_golem_options("database_names")
    DB_ROOT <- golem::get_golem_options("db_root")
    INSTALL_TYPE <- golem::get_golem_options("install_type") #production, dev,

    # potential / modal values
    modal_db_info = reactiveValues(db_root = DB_ROOT, #"UNDEFINED",
                                   db_names = DB_NAMES,  #list(),   ## named list
                                   finished_db_root=NULL,
                                   finished_db_name=NULL,
                                   config_ready=FALSE,
                                   meta="UNDEFINED",  #meta info ()
                                   root_trig=0,
                                   db_trig=0)  ## This prevents the values$to_print output fro

    # vetted values
    db <- reactiveValues(name=NULL, #current
                              dir=NULL, #current
                              root=NULL, #db_root
                              list=NULL, #all databases
                              regenerate=TRUE) #always regenerate on initial load...

    #instead of shinyFiles i could use find .hd5 and strip of the two parents in the path...
    shinyFiles::shinyDirChoose(input, "database_root_path", roots = c(shinyFiles::getVolumes()()), session = session, allowDirCreate = FALSE)
    mod_additional_info_server("additional_info_ui_ingest",
                               db_info = db)

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
      VMR = NULL
    )

    idx_of_disp <- reactive({
      matrixStats::colVars(raw,na.rm = TRUE)/colMeans(raw,na.rm = TRUE)
    })

    #TODO:  remove shaddow_defs and unify with temp_rv / to_return
    shaddow_defs <- reactiveValues(
      # holding the interactive selected factors and features which supercede what is from the defaults/configs
      exp_fact = NULL,
      exp_annot = NULL,
      omic_feat = NULL,
      feat_annot = NULL,
      feat_deets = NULL,
      feature_filter = NULL  #TODO:  change this to feature_filter_name
    )

    temp_rv <- reactiveValues(
      # holding the interactive selected factors and features which supercede what is from the defaults/configs
      anndata = NULL,
      de = NULL,
      config = NULL,
      default = NULL,
      VMR = NULL  #TODO:  change this to feature_filter_values
    )


    # modal comlete!  get select database ready...
    observe({
      trig1 <- isolate(modal_db_info$root_trig )
      if ( modal_db_info$db_root=="UNDEFINED") {
        modal_db_info$root_trig <- trig1 + 1
        message(paste0("trigger1: ",trig1))
      }
    })

    observe({
      req(db$root)
      trig2 <- isolate(modal_db_info$db_trig )
      if ( !any(modal_db_info$db_names != "UNDEFINED") ) {
          modal_db_info$db_trig <- trig2 + 1
          message(paste0("trigger2: ",trig2))
        } else {
          db$list <- isolate(modal_db_info$db_names)
        }

    })

    # event if either button press or "not a database name that is not UNDEFINED"
    virtual_add_database_button <- reactive({
      paste(modal_db_info$db_trig, input$AB_add_database_load)
    })

    # db add MODAL
    observeEvent(
      virtual_add_database_button(),{
        message("add database modal: ")

        in_db_dirs <- list.dirs(path = db$root, recursive = FALSE)
        in_db_dirs <- basename(in_db_dirs) #just keep the directory names not the full paths
        #names(in_db_dirs) <- basename(in_db_dirs)

        req(in_db_dirs)
        showModal(modalDialog(
          title="Add Database",
          size = "m",
          easyClose = FALSE,
          fade = TRUE,
          footer = tagList(
            modalButton("Cancel"),
            shinyjs::disabled(actionButton(ns("ab_save_config"), label= "Update & Close", icon = NULL, width = NULL))
          ),
          #uiOutput(ns("ui_ingest_modal")),
          p("databases"),
          #need 1. a directory selector (button modal?)
          #     2. select and name databases ()
          #     3. text field... annotate (upload .doc or write text)
          #
          #     list.dirs(path = ".",recursive=FALSE)
          selectizeInput(
            ns("database_folder"), "Database Folder",  choices = in_db_dirs, selected = in_db_dirs[1],
            multiple=FALSE
            ),

          textInput(ns("txt_database_name"),label = "Database Name", placeholder="what to call it"),
          actionButton(ns("ab_add_database_name"), label= "Add Database!", icon = NULL, width = NULL),
          textInput(ns("txt_study_meta"),label = "Database context", placeholder=" "),
          actionButton(ns("ab_study_meta"), label= "Add Info", icon = NULL, width = NULL)
        ))
        modal_db_info$finished_db_name <- FALSE # assert

    },ignoreInit = TRUE
    )

    # Modal 1
    observeEvent(modal_db_info$root_trig,{
      message("db_root modal: ")
      if (modal_db_info$db_root != "UNDEFINED" ){ #triggered on init...
        db$root <- modal_db_info$db_root
        modal_db_info$finished_db_root <- TRUE
      } else {
        showModal(modalDialog(
          title="Configure Database Root",
          size = "m",
          easyClose = FALSE,
          fade = TRUE,
          footer = tagList(
            #modalButton("Cancel"),
            actionButton(ns("ab_save_db_root"), label= "Save DB root", icon = NULL, width = NULL)
          ),
          p("Choose the path where database folders are located (`db root`)"),
          shinyFiles::shinyDirButton(id=ns("database_root_path"), label="Input directory", title="DB root path"),
          verbatimTextOutput(ns("txt_database_root_path")),
        ))
        modal_db_info$finished_db_root <- FALSE #ASSERT
      }
    })

    observe({
      selected_root <- shinyFiles::parseDirPath(roots = c(shinyFiles::getVolumes()()), input$database_root_path)
      # restrictions = system.file(package = "base")
      output$txt_database_root_path <- renderPrint({
        if (is.integer(input$database_root_path)) {
          cat("No directory has been selected (shinyDirChoose)")
        } else {
          #as.character(input$database_root_path)
          as.character(selected_root)
        }
      })

      if (length(selected_root)>0){
        if (dir.exists(paths = selected_root)){
          modal_db_info$db_root <- selected_root
          }
        } else {
          message("no db chosen yet..")
      }

    })

    observeEvent(input$ab_save_db_root, {
      # update config and close modal
      if ( modal_db_info$db_root != "UNDEFINED"){
        out_conf <- omicser::get_config()
        out_conf$db_root_path <- modal_db_info$db_root
        omicser::write_config(out_conf)
        removeModal(session = session)
        modal_db_info$finished_db_root <- TRUE
        db$regenerate <- TRUE
        db$root <- out_conf$db_root_path
      } else {
        message("urkk! [save db_root]  still undefined")

      }
    })

    # enable / disable Update&close button
    # if configuration is done, it is ok to enable button
    observeEvent(input$AB_add_database_load, {
      if (modal_db_info$config_ready) {
        shinyjs::enable("ab_save_config")
      } else {
        shinyjs::disable("ab_save_config")
        message(" add named database")
      }
    })

    observeEvent(input$ab_save_config, {
      # update config and close modal
      # should be impossible for db_root to be UNDEFINED, but safe
      if ( modal_db_info$db_root!="UNDEFINED" && any(modal_db_info$db_names != "UNDEFINED") ){
        out_conf <- omicser::get_config()
        out_conf$database_names <- modal_db_info$db_names
        out_conf$db_root_path <- modal_db_info$db_root # do I need to re-assert this??
        omicser::write_config(out_conf)
        modal_db_info$finished_db_name <- TRUE
        removeModal(session = session)
        db$regenerate <- TRUE
        db$list <- isolate(modal_db_info$db_names)
      } else {
        message("eek! [save config] something still undefined")
      }

    })

    observeEvent(input$ab_add_database_name, {
      req(input$database_folder)
      new_db <- input$database_folder

      new_db_nm <- input$txt_database_name
      if (length(new_db_nm)>0){
        names(new_db) <-new_db_nm
      } else {
        names(new_db) <- new_db
      }

      curr_db_names <- isolate(modal_db_info$db_names)
      if (curr_db_names == "UNDEFINED"){
        curr_db_names <- list()
      }
      if (!(new_db %in% curr_db_names) ){
        modal_db_info$db_names <- c(curr_db_names,new_db)
      } else { #the added one was redundant... and db$list is not null
        # just update the name
        names(curr_db_names) <- names(new_db)
        modal_db_info$db_names <- curr_db_names
      }

    })


    observeEvent(input$ab_study_meta, {
      # TODO: interact with the db_meta.yml
      #if (input$txt_study_meta ) is a .txt .doc etc file then load and render it
      #
      #else render the text as parameter in an .Rmd?
      # add

      # render the text somewhere...
      # nput$txt_database_name
      db$meta <- input$txt_study_meta
      print(input$txt_study_meta)
    })


    # database modal observer...
    observe({
      if ( any(modal_db_info$db_names != "UNDEFINED") ) { # fallback has "UNDEFINED"
        modal_db_info$config_ready <- TRUE
      } else {
        modal_db_info$config_ready <- FALSE
      }

    })



    observe({
      req(db$list) # leave this null until its ready.
      updateSelectizeInput(session, "SI_database", choices = db$list, selected = db$list[1], server=TRUE)
    })


    # observe({
    #   db$name <- to_return$db_meta$name
    # })


    ## load dataset (observeEvent) ===================
    #observeEvent(input$SI_database, {
    observe({
        req(input$SI_database)
      db$dir <- input$SI_database
      db$name <- names(which(db$list==input$SI_database))
      # zero out the data_meta until we "load"
      output$ui_data_meta <- renderUI({h4("Press `LOAD` to activate data") })
      # load data
      temp_rv$de <- readRDS(file = file.path(db$root,db$dir,"db_de_table.rds"))
      anndata <- anndata::read_h5ad(filename=file.path(db$root,db$dir,"db_data.h5ad"))
      temp_rv$anndata <- anndata
      conf_def <- gen_config_table(anndata, db$dir, db$root, regenerate=db$regenerate)
      db$regenerate <- FALSE
      #update reactives...
      temp_rv$config <- conf_def$conf
      temp_rv$default <- conf_def$def
      # computed index of dispersion
      # TODO: fix this HACK, change to VMR (index of dispersion)
      # computer no-matter what,
      if (min(anndata$X)<=0) {
        tmp_mu <- abs(colMeans(anndata$X,na.rm = TRUE))
        tmp_vmr <- matrixStats::colVars(anndata$X,na.rm = TRUE)/tmp_mu
        # set vmr to zero when mean is zero
        tmp_vmr[tmp_mu==0] <- 0
      } else {
        tmp_X <- log1p(anndata$X)
        tmp_mu <- colMeans(tmp_X,na.rm = TRUE)
        tmp_vmr <- matrixStats::colVars(tmp_X,na.rm = TRUE)-tmp_mu

      }

      temp_rv$VMR <- tmp_vmr  #actuall logVMR

      obs_choices <- anndata$obs_keys()

      #group_obs <- temp_rv$config[grp == TRUE & field=="obs"]$UI # <- choices_x
      #def_grp_o <- temp_rv$config[grp == TRUE & field=="obs" & default==1]$ID
      #will still be disabled if we are calculating...
      updateSelectInput(inputId = "SI_exp_fact", #label="Choose sample factors"
                        choices = obs_choices, #multiple=TRUE
                        selected = temp_rv$default$group_obs)

      updateSelectInput(inputId = "SI_exp_annot", #label="Choose sample-meta annotations (heatmap)"
                        choices = obs_choices, #multiple=TRUE
                        selected = temp_rv$default$obs_annots)

      var_choices <- anndata$var_keys()
      #group_var <- temp_rv$config[grp == TRUE & field=="var"]$UI # <- choices_x
      #group_var <- temp_rv$config[grp == TRUE & field=="var"]$UI # <- choices_x

      updateSelectInput(inputId = "SI_omic_feat", #label="Choose omic feature groups"
                        choices = var_choices, #multiple=TRUE
                        selected = temp_rv$default$group_var)

      updateSelectInput(inputId = "SI_feat_annot", #label= "Choose omic feature annotatins (heatmap)"
                        choices = var_choices, #multiple=TRUE
                        selected = temp_rv$default$var_annots)

      updateSelectInput(inputId = "SI_feat_deets", #label= "Choose omic feature meta-info"
                        choices = var_choices, #multiple=TRUE
                        selected = temp_rv$default$feature_details)

      updateSelectInput(inputId = "SI_feature_filter",# label="Choose omic feature for filtering `highly variable`"
                        choices = var_choices, #multiple=FALSE
                        selected = temp_rv$default$filter_feature[1])

      shinyjs::enable("AB_update_defaults")

    })

    ## observes ===================
    observe({
      if ( !is.null( db$name ) ) {
        shinyjs::enable("AB_ingest_load")
      } else {
        shinyjs::disable("AB_ingest_load")
        message(" no database loaded yet... .")
      }
    })

    observe({
      if ( input$CB_var_to_mean_ratio ) {
        shinyjs::disable("SI_feature_filter")
        # compute and set ? value??

      } else {
        shinyjs::enable("SI_feature_filter")
        message(" choose.. .")
      }
    })

    # update defaults: write to configuration
    observeEvent(input$AB_update_defaults, {
      # all other return values set with SI_database
      #conf_list_old <- configr::read.config( file.path(db$root,db$dir,"db_config.yml" ) )

      conf_list_old <- omicser::get_db_conf(db$dir, db_root = db$root)
      conf_list_old$group_obs <- input$SI_exp_fact
      conf_list_old$obs_annots <- input$SI_exp_annot
      conf_list_old$group_var <- input$SI_omic_feat
      conf_list_old$var_annots <- input$SI_feat_annot
      conf_list_old$feature_details <- input$SI_feat_deets
      conf_list_old$filter_feature <- input$SI_feature_filter

      # use a modal?  or simply send everything from the ingest UI as update?
      omicser::write_db_conf(conf_list_old,db$dir, db_root =  db$root)

      # set flag to regenerate the configuration gen_config_table
      db$regenerate <- TRUE

      }
    )

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


      to_return$db_meta$name<-db$name
      to_return$db_meta$db_dir<-db$dir



      to_return$db_meta$omics_type<-temp_rv$default$omic_type
      to_return$db_meta$meta_info<-temp_rv$default$meta_info

      to_return$trigger <- to_return$trigger + 1

      shaddow_defs$exp_fact <- input$SI_exp_fact
      shaddow_defs$exp_annot <- input$SI_exp_annot
      shaddow_defs$omic_feat <- input$SI_omic_feat
      shaddow_defs$feat_annot <- input$SI_feat_annot
      shaddow_defs$feat_deets  <- input$SI_feat_deets
      shaddow_defs$feature_filter  <- if (input$CB_var_to_mean_ratio) "VMR" else input$SI_feature_filter

      to_return$shaddow_defs <- shaddow_defs
      to_return$VMR <- temp_rv$VMR


      # will this work?  or is the "observing" affecting a reactive with this render?
      output$ui_data_meta <- renderUI({
        out_text <- HTML(paste("from the <i>", temp_rv$default$lab),"</i> lab")
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



  }

)}

## To be copied in the UI
# mod_ingestor_ui("ingestor_ui_1")

## To be copied in the server
# mod_ingestor_server("ingestor_ui_1")

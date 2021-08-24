

#' side_selector UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_side_selector_ui <- function(id){
  ns <- NS(id)

  quant_plot_types <- c("Volcano Plot"="volcano", "Dot Plot"="dotplot","Heat Map"="hmap","Box/violin Plot"="box")
  comp_plot_types <- c("Volcano Plot"="volcano", "Dot Plot"="dotplot","Heat Map"="hmap","Box/violin Plot"="box")

  selector_tags <- tagList(
    # helpText(
    #    HTML("Load a Dataset first, and then selected the variables and quantiteis below.")
    #    ),
    #  hr(style = "border-top: 1px solid #000000;"),
    # #TODO:  change these to render prper HTML

      htmlOutput( ns( "ui_curr_database" )),
      uiOutput( ns( "ui_DIV_warn" )),
      htmlOutput( ns( "ui_db_type" )),

      fluidRow(
        hr(style = "border-top: 1px solid #000000;"),
        column(
          width=2,
          "choose:"
        ),
        column(
          width = 5,
          offset = 0,
          shinyjs::disabled(radioButtons( ns("RB_obs_var"), "variable type",
                                          choices = c("obs","var","X"),
                                          selected = "X" , inline = TRUE)
          )
        )
      ),# options = list(placeholder = ""))
      fluidRow(
      column(
        width=2,
        "Choose:",
        "variables",
        "X-axis"
      ),
      column(
        width=5,
        offset=0,
        #h5("Sample- information select"),
        # "Choose cells/samples, grouping/subsetting variables, observables, and omic features ",
        #"to visualize in the playground",
        #br(),br(),
        shinyjs::disabled(
          selectizeInput(ns("SI_obs_x"), "observ- (obs)", choices=NULL) # options = list(placeholder = ""))
        )
        #actionButton(ns("AB_tobble_obs_var"), "Toggle obs/var"), #TODO: change text when toggled
      ),
      column(
        width=5,
        offset=0,
        shinyjs::disabled(
          selectizeInput(ns("SI_var_x"), "omic- (var)", choices=NULL)
        )
      )
    ),
    fluidRow(
      hr(style = "border-top: 1px solid #000000;"),
      column(
        width=2,
        "Choose:",
        "measure",
        "Y-axis"
      ),
     column(
        width = 5,
        offset = 0,
        shinyjs::disabled(
          selectizeInput(ns("SI_obs_y"), "quant (obs):", choices=NULL)
        )
      ),
    column(
      width = 5,
      offset = 0,
      shinyjs::disabled(
        selectizeInput(ns("SI_var_y"), "quant (var):", choices=NULL)
       )
      )

      ),

  hr(style = "border-top: 1px solid #000000;"),
  fluidRow(
    column(
      width=2,
      "group-by:"
    ),
    column(
      width = 5,
      offset = 0,
      shinyjs::disabled(
        selectizeInput(ns("SI_obs_groupA"), "groupby (obs):", choices=NULL)
      )
    ),
    column(
      width = 5,
      offset = 0,
      shinyjs::disabled(
        selectizeInput(ns("SI_var_groupA"), "groupby (var):", choices=NULL)
      )
    )

  ),
  hr(style = "border-top: 1px solid #000000;"),
  ## TODO:  make this entire row a dynamic uiOutput

    fluidRow(
      column(
        width = 5,
        shinyjs::disabled(
          selectizeInput(ns("SI_obs_subset"),
                         "Obs | subset:",
                         "",
                         multiple = FALSE,
                         options = list(placeholder = "choose dataset first"))
        ),
        br(),
        shinyjs::disabled(
          actionButton(ns("CB_obs_sub_none"), "<none>", class = "btn btn-primary")
          ),
        shinyjs::disabled(
          actionButton(ns("CB_obs_sub_all"), "<all>", class = "btn btn-primary"))
        ),
      column(width = 7,
        offset = 0,
        uiOutput(ns("ui_obs_subset"))
        )
      ),

## TODO:  fix var subsetting to only load if we have "grouping"
  # fluidRow(
  #   column(
  #     width = 4,
  #    #actionButton(ns("AB_subset_tog"), "Toggle subset observations"), #TODO: change text when toggled
  #    shinyjs::disabled(
  #      selectizeInput(ns("SI_var_subset"),
  #                     "Var | subset:",
  #                     "",
  #                     multiple = FALSE,
  #                     options = list(placeholder = "choose dataset first"))
  #     ),
  #     br(),
  #     shinyjs::disabled(
  #       actionButton(ns("CB_var_sub_none"), "<none>", class = "btn btn-primary")
  #     ),
  #     shinyjs::disabled(
  #       actionButton(ns("CB_var_sub_all"), "<all>", class = "btn btn-primary"))
  #   ),
  #   column(width = 8,
  #          offset = 0,
  #          uiOutput(ns("ui_var_subset"))
  #   )
  # ),
  hr(style = "border-top: 1px solid #000000;"),

    mod_omic_selector_ui(ns("omic_selector_ui_1")),

    fluidRow(
      uiOutput(ns("ui_text_warn"), width = "100%"),
        ), #fluidRow

    fluidRow(
      column(
        width = 10,
        offset = 1,
        uiOutput(ns("ui_groupB"))
      )
    )

  ) #taglist

  return(selector_tags)

}

#' side_selector Server Functions
#'
#' @noRd
mod_side_selector_server <- function(id, rv_in){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    ### Reactive expressions ============================================
    out_params <- reactiveValues(
      #omics_names = NULL,
      omics_list = NULL,
      # aggregate obs
      feat_grp = NULL,
      feat_subset = NULL,
      feat_subsel = NULL,
      feat_x = NULL,
      feat_y = NULL,
      observ_grpA = NULL,
      observ_subsetA = NULL,
      observ_subselA = NULL,
      observ_grpB = NULL,
      # observ_subsetB = NULL,
      # observ_subselB = NULL,
      observ_x = NULL,
      observ_y = NULL,
      # observ_y_comp = NULL,
      # until this is set we don't actually plot ayhting.
      measure_type = NULL, #"raw" or "comp"
      # TODO:  make this conditional on what kind of data is loaded.. and just a single plot type
      raw_plot_type = NULL,
      comp_plot_type = NULL,
      obs_type = NULL
    )


    all_omics <- reactive( names(rv_in$omics) )
    def_omics <- reactive( rv_in$default$omics)
    new_db_trig <- reactive( rv_in$trigger )
    omics_list <- mod_omic_selector_server("omic_selector_ui_1", all_omics ,def_omics,rv_in$trigger)


    ### Outputs =========================================================

    output$ui_curr_database <- renderUI({
      if (is.null(rv_in$database_name)) {
        out_text <- "No data loaded"
        } else {
          out_text <- paste("Current dataset: ", rv_in$database_name)
        }
      out_text <- h4(out_text)
      return(out_text)
    })

    # Warning if no data loaded
    output$ui_DIV_warn <- renderUI( {
      if (is.null(rv_in$database_name)) {
        div(
          tags$br(),
          span(class = "warn", "No dataset loaded")
         )
        }
      })

    # show the factors that have been loaded
    output$ui_db_type <- renderUI({
      req(rv_in$omics_type)
      out_text <- paste("<h6>Data type: <i>", rv_in$omics_type, "</i>-omics</h6>")
      out_text <- HTML(out_text)
      return(out_text)
      })

    observeEvent(
      (rv_in$trigger>0),  #why does this happen twice?
      {
        req(rv_in$config,  # set when database is chosen
            rv_in$default)

        ### INITIALIZE ALL UI ELEMENTS ++++++++++++++++++++++++

        rv_in$trigger <- 0

        # update full omics list
        all_omics <- rv_in$omics


        # OBSERVATIONS / VARIABLES (x-axis) -------------------------------
        # GROUP BY -------------------------------
        # SUBSET -------------------------------
        choices_obs_x <- rv_in$config[grp == TRUE & field=="obs"]$UI
        choices_var_x = rv_in$config[grp == TRUE & field=="var"]$UI


        if (length(choices_obs_x)>0){
          freezeReactiveValue(input, "SI_obs_x")
          updateSelectizeInput(session, "SI_obs_x", server = TRUE,
                               choices = choices_obs_x,
                               selected = rv_in$default$obs_x[1])
          freezeReactiveValue(input, "SI_obs_groupA")
          updateSelectizeInput(session, "SI_obs_groupA", server = TRUE,
                               choices = choices_obs_x,
                               selected = rv_in$default$obs_x[1])
          freezeReactiveValue(input, "SI_obs_subset")
          updateSelectizeInput(session, "SI_obs_subset",
                               choices =  choices_obs_x,
                               selected = rv_in$default$grp2,  server = TRUE)

          enable_flagO <- TRUE
        } else {
          enable_flagO <- FALSE
        }

        if (length(choices_var_x)>0){
          freezeReactiveValue(input, "SI_var_x")
          updateSelectizeInput(session, "SI_var_x", server = TRUE,
                               choices = choices_var_x,
                               selected = rv_in$default$var_x[1])
          # same defaults/choices as SI_obs_x
          freezeReactiveValue(input, "SI_var_groupA")
          updateSelectizeInput(session, "SI_var_groupA", server = TRUE,
                               choices = choices_var_x,
                               selected = rv_in$default$var_x[1])
          freezeReactiveValue(input, "SI_var_subset")
          updateSelectizeInput(session, "SI_var_subset",
                               choices = choices_var_x,
                               selected = rv_in$default$var_subset,  server = TRUE)

          enable_flagV <- TRUE
        } else {
          enable_flagV <- FALSE
        }

        # QUANTITIES / MNEASURES (y-axis) -------------------------------
        choices_var_y = rv_in$config[measure == TRUE & field=="var"]$UI
        choices_obs_y <- rv_in$config[measure == TRUE & field=="obs"]$UI


        if (length(choices_obs_y)>0){
          freezeReactiveValue(input, "SI_obs_y")
          updateSelectizeInput(session, "SI_obs_y",
                             choices = choices_obs_y,
                             selected = choices_obs_y[1],  server = TRUE) #rv_in$default$obs_y[1])
          enable_flagO <- enable_flagO
        } else {
          enable_flagO <- FALSE
        }

        if (length(choices_var_y)>0){
          freezeReactiveValue(input, "SI_var_y")
          updateSelectizeInput(session, "SI_var_y",
                               choices = choices_var_y,
                               selected = choices_var_y[1],  server = TRUE)
          enable_flagV <- enable_flagV
        } else {
          enable_flagV <- FALSE
        }


        #TODO:  make the choices *named* paste0(rv_in$config$fID, fUI)
        # dat_source = rv_in$config$mat[fID == rv_in$aux_raw]$ID
        # paste(raw_choices, dat_sourcce)


        if (enable_flagO) {
          shinyjs::enable("SI_obs_subset")
          shinyjs::enable("SI_obs_groupA")
          shinyjs::enable("CB_obs_sub_all")
          shinyjs::enable("CB_obs_sub_none")
          shinyjs::enable("SI_obs_x")
          shinyjs::enable("SI_obs_y")
          if (enable_flagV) {
            shinyjs::enable("RB_obs_var")
          }

        } else {
          if (enable_flagV) {
            # udpate cb to select var
            freezeReactiveValue(input, "RB_obs_var")
            updateRadioButtons(session, "RB_obs_var", "variable type",
                               choices = c("obs","var"),
                               selected = "var" , inline = TRUE)
            }
          }

        if (enable_flagV) {
          shinyjs::enable("SI_var_x")
          shinyjs::enable("SI_var_y")
          shinyjs::enable("SI_var_subset")
          shinyjs::enable("SI_var_groupA")
          shinyjs::enable("CB_var_sub_all")
          shinyjs::enable("CB_var_sub_none")
        } # disabled by default


      },
    ignoreNULL = TRUE,
    ignoreInit = TRUE
    ) #observe event



    observe({

      if (input$RB_obs_var == "obs") {
        shinyjs::enable("SI_obs_y")
        shinyjs::enable("SI_obs_x")
        print("disabled obs  ")
        shinyjs::disable("SI_var_y")
        shinyjs::disable("SI_var_x")

        # choices_obs_y = rv_in$config[measure == TRUE & field=="obs"]$UI
        # choices_obs_x = rv_in$config[grp == TRUE & field=="obs"]$UI
        #
        # freezeReactiveValue(input, "SI_obs_y")
        # updateSelectizeInput(session, "SI_obs_y",
        #                      choices = choices_obs_y,
        #                      selected = choices_obs_y[1],  server = TRUE) #rv_in$default$obs_y[1])
        #
        # freezeReactiveValue(input, "SI_obs_x")
        # updateSelectizeInput(session, "SI_obs_x", server = TRUE,
        #                      choices = choices_obs_x,
        #                      selected = rv_in$default$obs_x[1])
      } else { #(input$RB_obs_var == "obs")

        shinyjs::enable("SI_var_x")
        shinyjs::enable("SI_var_y")
        print("disabled var  ")
        shinyjs::disable("SI_obs_y")
        shinyjs::disable("SI_obs_x")

        # choices_var_y = rv_in$config[measure == TRUE & field=="var"]$UI
        # choices_var_x = rv_in$config[grp == TRUE & field=="var"]$UI
        # freezeReactiveValue(input, "SI_var_y")
        # updateSelectizeInput(session, "SI_var_y",
        #                      choices = choices_var_y,
        #                      selected = choices_var_y[1],  server = TRUE)
        #
        # freezeReactiveValue(input, "SI_var_x")
        # updateSelectizeInput(session, "SI_var_x", server = TRUE,
        #                      choices = choices_var_x,
        #                      selected = rv_in$default$var_x[1])
       # TODP: try and set the choices to "var"
        # updateSelectizeInput(session, "SI_obs_x", server = TRUE,
        #                      choices = rv_in$config[grp == TRUE]$UI,
        #                      selected = rv_in$default$obs_x[1])
        # return_value$obs_type <- "var"


      }

    })

    # observe({
    #   req(rv_in$config)
    #   updateSelectizeInput(session, "SI_obs_x", server = TRUE,
    #                        choices = rv_in$config[grp == TRUE]$UI,
    #                        selected = rv_in$default$obs_x[1])
    # })
    # observe({
    #   req(rv_in$config)
    #   raw_choices = rv_in$config[measure == TRUE]$UI
    #   #TODO:  make the choices *named* paste0(rv_in$config$fID, fUI)
    #   # dat_source = rv_in$config$mat[fID == rv_in$aux_raw]$ID
    #   # paste(raw_choices, dat_sourcce)
    #
    #   if (length(raw_choices)>0) {
    #     freezeReactiveValue(input, "SI_obs_y")
    #     updateSelectizeInput(session, "SI_obs_y","Raw vals:",
    #                          choices = raw_choices,
    #                          selected = raw_choices[1],  server = TRUE)
    #   } else {
    #     print("disabled  raw observations ")
    #     shinyjs::disable("SI_obs_y")
    #     # change back to placeholder??
    #     freezeReactiveValue(input, "SI_obs_y")
    #     updateSelectizeInput(session, "SI_obs_y", "Raw vals: ", "", options = list(placeholder = ""))
    #   }
    #
    # })
    #


    #
    # conditionalPanel(
    #   condition = "input.Van_a2togL % 2 == 1",
    #   selectInput("Van_a2sub1", "Cell information to subset:",
    #               choices = Van_conf[grp == TRUE]$UI,
    #               selected = Van_def$grp1),
    #   uiOutput("Van_a2sub1.ui"),
    #   actionButton("Van_a2sub1all", "Select all groups", class = "btn btn-primary"),
    #   actionButton("Van_a2sub1non", "Deselect all groups", class = "btn btn-primary")
    # )


    output$ui_obs_subset <- renderUI({
      req(input$SI_obs_subset)

        if (input$SI_obs_subset == ""){
          print("loading...SI_obs_subset.")
          return()
        }
        cfg <- isolate(rv_in$config)

        subs <- strsplit(cfg[UI == input$SI_obs_subset]$fID, "\\|")[[1]]

        subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
        if (all(!is.na(subs2))){
          names(subs2) <- seq_along(subs2)
          subs <- subs[as.numeric(names(sort(subs2)))]
        }

        subs_label <- paste0("select ",isolate(input$SI_obs_subset),"s: ")

        checkboxGroupInput( ns("CB_obs_subsel"),
                           label = subs_label,
                           inline = TRUE,
                           choices = subs,
                           selected = subs)

    })


    output$ui_var_subset <- renderUI({
      req(input$SI_var_subset) # do i need this?

        if (input$SI_var_subset == ""){
          print("loading...ui_var_subset.")
          return()
        }
        cfg <- isolate(rv_in$config)

        subs <- strsplit(cfg[UI == input$SI_var_subset]$fID, "\\|")[[1]]
        subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
        if (all(!is.na(subs2))){
          names(subs2) <- seq_along(subs2)
          subs <- subs[as.numeric(names(sort(subs2)))]
        }

        subs_label <- paste0("select ",isolate(input$SI_var_subset),"s: ")

        checkboxGroupInput( ns("CB_var_subsel"),
                            label = subs_label,
                            inline = TRUE,
                            choices = subs,
                            selected = subs)

    })


    #
    # output$ui_groupB <- renderUI({
    #
    #   cfg <- isolate(rv_in$config)
    #
    #   subs <- strsplit(cfg[UI == input$SI_obs_subset]$fID, "\\|")[[1]]
    #   subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
    #   names(subs2) <- seq_along(subs2)
    #   subs <- subs[as.numeric(names(sort(subs2)))]
    #
    #   subs_label <- paste0("select ",isolate(input$SI_obs_subset),"s: ")
    #
    #   checkboxGroupInput( ns("CB_obs_subsel"),
    #                       label = subs_label,
    #                       inline = TRUE,
    #                       choices = subs,
    #                       selected = subs)
    # })
    #

    ### observes =========================================================
    observeEvent(input$CB_obs_sub_all, {
      #req(input$CB_obs_subsel) this makes the button useless if everything is deselected
      req(rv_in$meta)
      cfg <- isolate(rv_in$config)

      subs <- strsplit(cfg[UI == input$SI_obs_subset]$fID, "\\|")[[1]]
      subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
      if (all(!is.na(subs2))){
        names(subs2) <- seq_along(subs2)
        subs <- subs[as.numeric(names(sort(subs2)))]
      }

      subs_label <- paste0("select ",isolate(input$SI_obs_subset),"s: ")

      freezeReactiveValue(input, "CB_obs_subsel")
      updateCheckboxGroupInput(session,
                               inputId = "CB_obs_subsel",
                               label = subs_label,
                               choices = subs,
                               selected = subs,
                               inline = TRUE) # WARNING.  make sure subset is null is checked (length(0?))
    })


    observeEvent(input$CB_obs_sub_none, {
      req(input$CB_obs_subsel) #short circuits if already de-selected
      req(rv_in$meta)
      cfg <- isolate(rv_in$config)

      curr_sel <- isolate(input$CB_obs_subsel)
      subs <- strsplit(cfg[UI == input$SI_obs_subset]$fID, "\\|")[[1]]
      subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
      if (all(!is.na(subs2))){
        names(subs2) <- seq_along(subs2)
        subs <- subs[as.numeric(names(sort(subs2)))]
      }

      subs_label <- paste0("select ",isolate(input$SI_obs_subset),"s: ")

      freezeReactiveValue(input, "CB_obs_subsel")
      updateCheckboxGroupInput(session,
                               inputId = "CB_obs_subsel",
                               label = subs_label,
                               choices = subs,
                               selected = character(0), # character(0), #, "" , NULL
                               inline = TRUE) # WARNING.  make sure subset is null is checked (length(0?))
    })


    observeEvent(input$CB_var_sub_all, {
      #req(input$CB_obs_subsel) this makes the button useless if everything is deselected
      req(rv_in$meta)
      if (input$SI_var_subset == ""){
        return()
      }

      cfg <- isolate(rv_in$config)
      subs <- strsplit(cfg[UI == input$SI_var_subset]$fID, "\\|")[[1]]
      subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
      if (all(!is.na(subs2))){
        names(subs2) <- seq_along(subs2)
        subs <- subs[as.numeric(names(sort(subs2)))]
      }

      subs_label <- paste0("select ",isolate(input$SI_var_subset),"s: ")



      freezeReactiveValue(input, "CB_var_subsel")
      updateCheckboxGroupInput(session,
                               inputId = "CB_var_subsel",
                               label = subs_label,
                               choices = subs,
                               selected = subs,
                               inline = TRUE) # WARNING.  make sure subset is null is checked (length(0?))
    })


    observeEvent(input$CB_var_sub_none, {
      req(input$CB_var_subsel) #short circuits if already de-selected
      req(rv_in$meta)
      cfg <- isolate(rv_in$config)


      curr_sel <- isolate(input$CB_var_subsel)
      subs <- strsplit(cfg[UI == input$SI_var_subset]$fID, "\\|")[[1]]
      subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
      if (all(!is.na(subs2))){
        names(subs2) <- seq_along(subs2)
        subs <- subs[as.numeric(names(sort(subs2)))]
      }

      subs_label <- paste0("select ",isolate(input$SI_var_subset),"s: ")

      freezeReactiveValue(input, "CB_var_subsel")
      updateCheckboxGroupInput(session,
                               inputId = "CB_var_subsel",
                               label = subs_label,
                               choices = subs,
                               selected = character(0), # character(0), #, "" , NULL
                               inline = TRUE) # WARNING.  make sure subset is null is checked (length(0?))
    })


   #  observe({
   #
   #    req(rv_in$config)
   #    comp_choices = rv_in$config[diff_exp == TRUE ]$UI
   #
   #
   #    if (length(comp_choices)>0) {
   #      freezeReactiveValue(input, "SI_obs_comp")
   #      updateSelectizeInput(session, "SI_obs_comp","Comparatives:",
   #                           choices = comp_choices,
   #                           selected = comp_choices[1],  server = TRUE)
   #    } else {
   #      print("disabled  comparative observations ")
   #      shinyjs::disable("SI_obs_comp")
   #      freezeReactiveValue(input, "SI_obs_comp")
   #      updateSelectizeInput(session, "SI_obs_comp", "Comparatives:", "",
   #                      options = list(placeholder = "") )
   #    }
   #
   #
   # })


    # raw2_choices = rv_in$config$mat[observ == TRUE & ID!="raw"]$fIDloc
    # comp2_choices = rv_in$config$mat[comp == TRUE ]$fIDloc
    #
    # names(comp2_choices) <- comp2_choices
    # names(raw2_choices) <- raw2_choices
    #
    # if (length(raw2_choices)<1){
    #
    #   if (length(comp2_choices)<1){
    #     #no Y!!!
    #     sel2 <- charachter(0)
    #     sel2_choices <- raw2_choices #NULL
    #
    #   } else {
    #     sel2 <- comp2_choices[1]
    #     sel2_choices <- comp2_choices
    #
    #   }
    # } else {
    #   if (length(comp2_choices)<1){
    #     sel2 <- raw2_choices[1]
    #     sel2_choices <- raw2_choices
    #   }  else {
    #     sel2_choices <- list(raw=raw2_choices , comp=comp2_choices)
    #     sel2 <- raw2_choices[1]
    #   }
    # }

    #TODO:  make the choices *named* paste0(rv_in$config$fID, fUI)
    # dat_source = rv_in$config$mat[fID == rv_in$aux_raw]$ID
    # paste(raw_choices, dat_sourcce)
    #{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }"

    #
      # updateSelectizeInput(session, "SI_y_info", server = TRUE,
      #                      choices = sel2_choices,
      #                      selected = sel2)
    # })
    # TODO:  wrap this in a module..
    #
    #


    observeEvent( omics_list$viz_now, {
      # route the chosen data type out...
      req(input$CB_obs_subsel,
          input$SI_obs_subset)
      if (omics_list$viz_now) {

        out_params$observ_grpA = input$SI_obs_groupA
        out_params$observ_subsetA = input$SI_obs_subset
        out_params$observ_subselA = input$CB_obs_subsel

        out_params$feat_grp = input$SI_var_groupA
        out_params$feat_subset = input$SI_var_subset
        out_params$feat_subsel = input$CB_var_subsel

        # # aggregate var
        # out_params$feat_grp <- input$SI_feat_grp
        # out_params$feat_subsel <- input$SI_feat_subsel

        out_params$observ_x <- input$SI_obs_x
        out_params$observ_y <- input$SI_obs_y
        # out_params$observ_y_rcomp <- input$SI_obs_comp
        #
        out_params$obs_type <- input$RB_obs_var

        out_params[["omics_list"]] <- omics_list  # value & viz_now

      }

    })


  ### RETURN =========================================================
  return(out_params)

  })
}

## To be copied in the UI
# mod_side_selector_ui("side_selector_ui_1")

## To be copied in the server
# mod_side_selector_server("side_selector_ui_1")

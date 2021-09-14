

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
        width=10,
        offset=0,
        h5("Sample/Omics- information select"),
        # "Choose cells/samples, grouping/subsetting variables, observables, and omic features ",
        "to visualize in the playground",
        br(),br(),
        selectizeInput(ns("SI_x_info"), "Annotation information (X-axis):", choices=NULL), # options = list(placeholder = ""))
        #actionButton(ns("AB_tobble_obs_var"), "Toggle obs/var"), #TODO: change text when toggled
        shinyjs::disabled(radioButtons( ns("RB_x_obs_var"), "variable type",
                               choices = c("obs","var"),
                               selected = "obs" ))
        )
      ),
    fluidRow(
      column(
        width=3,
        offset=0,
        br(),
        "",
        h5("Measures"),
        h5("(Y-Axis)")
        ),
      column(
        width = 8,
        offset = 0,
        selectizeInput(ns("SI_obs_raw"), "Raw vals", choices=NULL)
        )

      ),

    fluidRow(
      column(width = 5,
        offset = 1,
        #actionButton(ns("AB_subset_tog"), "Toggle subset observations"), #TODO: change text when toggled
        shinyjs::disabled(selectizeInput(ns("SI_groupA"), "Obs information to subset:", "",
                                         multiple = FALSE, options = list(placeholder = "choose dataset first"))),
      )
    ),
    fluidRow(
      column(width = 10,
        offset = 1,
        uiOutput(ns("ui_subsetA"))
        )
      ),
    fluidRow(
      column(
        width = 5,
        offset = 1,
        shinyjs::disabled(actionButton(ns("CB_sub_all"), "Select All", class = "btn btn-primary"))),
      column(
        width = 5,
        offset=0,
        shinyjs::disabled(actionButton(ns("CB_sub_none"), "Select None", class = "btn btn-primary") )
        )
      ),


    fluidRow(
      column(width = 10,
             offset = 1,
             uiOutput(ns("ui_groupB"))
      )
    ),


    mod_omic_selector_ui(ns("omic_selector_ui_1")),

    fluidRow(
      uiOutput(ns("ui_text_warn"), width = "100%"),
        ) #fluidRow


  ) #taglist

  return(selector_tags)

}

#' side_selector Server Functions
#'
#' @noRd
mod_side_selector_server <- function(id, rv_in){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    ############################ +
    ## reactive data structures
    ## out_params:  generated here and shared out
    ## rv_in:  input values shared by modules
    ############################ +
    out_params <- reactiveValues(

      #omics_names = NULL,
      omics_list = NULL,

      # aggregate obs
      feat_grp = NULL,
      feat_subsel = NULL,

      observ_grpA = NULL,
      observ_subselA = NULL,
      observ_grpB = NULL,
      observ_subselB = NULL,
      observ_x = NULL,
      observ_y_raw = NULL,
      # observ_y_comp = NULL,

      # until this is set we don't actually plot ayhting.
      measure_type = NULL, #"raw" or "comp"
      # TODO:  make this conditional on what kind of data is loaded.. and just a single plot type
      raw_plot_type = NULL,
      comp_plot_type = NULL,

      obs_type = NULL
    )

    # omics_for_selector <- reactive(rv_in$omics)
    # omics_list <- mod_omic_selector_server("omic_selector_ui_1", omics_for_selector)
    all_omics <- reactive( rv_in$omics )
    # pass_trig <- reactive( rv_in$trigger)
    # omics_list <- mod_omic_selector_server("omic_selector_ui_1", all_omics , pass_trig)
    #
    omics_list <- mod_omic_selector_server("omic_selector_ui_1", all_omics )

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



    ############################+
    ## update factors if data just loaded
    ############################+

    observeEvent(
      (rv_in$trigger>0),  #why does this happen twice?
      {
        req(rv_in$config)  # set when database is chosen
        req(rv_in$default)  # set when database is chosen

        #rv_in$trigger  <-  0
        shinyjs::enable("SI_groupA")
        freezeReactiveValue(input, "SI_groupA")
        updateSelectizeInput(session, "SI_groupA","Obs information to subset:",
                             choices = rv_in$config[grp == TRUE]$UI,
                             selected = rv_in$default$grp2,  server = TRUE)

        shinyjs::enable("CB_sub_all")
        shinyjs::enable("CB_sub_none")
        shinyjs::enable("SI_obs_raw")
        shinyjs::enable("RB_x_obs_var")

        # update omics_out
        all_omics <- rv_in$omics

        raw_choices = rv_in$config[measure == TRUE]$UI
        #TODO:  make the choices *named* paste0(rv_in$config$fID, fUI)
        # dat_source = rv_in$config$mat[fID == rv_in$aux_raw]$ID
        # paste(raw_choices, dat_sourcce)

        if (length(raw_choices)>0) {
          freezeReactiveValue(input, "SI_obs_raw")
          updateSelectizeInput(session, "SI_obs_raw","Raw vals:",
                               choices = raw_choices,
                               selected = raw_choices[1],  server = TRUE)
        } else {  #this should be impossible!
          print("disabled raw observations ")
          shinyjs::disable("SI_obs_raw")
          # change back to placeholder??
          freezeReactiveValue(input, "SI_obs_raw")
          updateSelectizeInput(session, "SI_obs_raw", "Raw vals: ", "", options = list(placeholder = ""))
        }

        updateSelectizeInput(session, "SI_x_info", server = TRUE,
                             choices = rv_in$config[grp == TRUE]$UI,
                             selected = rv_in$default$obs_x[1])

      },
    ignoreNULL = TRUE,
    ignoreInit = TRUE
    ) #observe event

    observe({
      if (input$RB_x_obs_var == "obs") {
        updateSelectizeInput(session, "SI_x_info", server = TRUE,
                             choices = rv_in$config[grp == TRUE]$UI,
                             selected = rv_in$default$obs_x[1])
        out_params$obs_type <- input$RB_x_obs_var
      } else {
       # TODP: try and set the choices to "var"
        print("haven't implimented plotting omics annots yet")
        # updateSelectizeInput(session, "SI_x_info", server = TRUE,
        #                      choices = rv_in$config[grp == TRUE]$UI,
        #                      selected = rv_in$default$obs_x[1])
        # return_value$obs_type <- "var"

      }

    })

    # observe({
    #   req(rv_in$config)
    #   updateSelectizeInput(session, "SI_x_info", server = TRUE,
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
    #     freezeReactiveValue(input, "SI_obs_raw")
    #     updateSelectizeInput(session, "SI_obs_raw","Raw vals:",
    #                          choices = raw_choices,
    #                          selected = raw_choices[1],  server = TRUE)
    #   } else {
    #     print("disabled  raw observations ")
    #     shinyjs::disable("SI_obs_raw")
    #     # change back to placeholder??
    #     freezeReactiveValue(input, "SI_obs_raw")
    #     updateSelectizeInput(session, "SI_obs_raw", "Raw vals: ", "", options = list(placeholder = ""))
    #   }
    #
    # })
    #


    output$ui_subsetA <- renderUI({
      #req(input$SI_groupA) # do i need this?
      if (input$SI_groupA == ""){
        print("loading....")
        return()
      }
      cfg <- isolate(rv_in$config)

      subs <- strsplit(cfg[UI == input$SI_groupA]$fID, "\\|")[[1]]
      subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
      names(subs2) <- seq_along(subs2)
      subs <- subs[as.numeric(names(sort(subs2)))]

      subs_label <- paste0("select ",isolate(input$SI_groupA),"s: ")

      checkboxGroupInput( ns("CB_sub_innerA"),
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
    #   subs <- strsplit(cfg[UI == input$SI_groupA]$fID, "\\|")[[1]]
    #   subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
    #   names(subs2) <- seq_along(subs2)
    #   subs <- subs[as.numeric(names(sort(subs2)))]
    #
    #   subs_label <- paste0("select ",isolate(input$SI_groupA),"s: ")
    #
    #   checkboxGroupInput( ns("CB_sub_innerA"),
    #                       label = subs_label,
    #                       inline = TRUE,
    #                       choices = subs,
    #                       selected = subs)
    # })
    #

    observeEvent(input$CB_sub_all, {
      #req(input$CB_sub_innerA) this makes the button useless if everything is deselected
      req(rv_in$meta)
      cfg <- isolate(rv_in$config)

      subs <- strsplit(cfg[UI == input$SI_groupA]$fID, "\\|")[[1]]
      subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
      names(subs2) <- seq_along(subs2)
      subs[as.numeric(names(sort(subs2)))]

      subs_label <- paste0("select ",isolate(input$SI_groupA),"s: ")



      freezeReactiveValue(input, "CB_sub_innerA")
      updateCheckboxGroupInput(session,
                               inputId = "CB_sub_innerA",
                               label = subs_label,
                               choices = subs,
                               selected = subs,
                               inline = TRUE) # WARNING.  make sure subset is null is checked (length(0?))
    })


    observeEvent(input$CB_sub_none, {
      req(input$CB_sub_innerA) #short circuits if already de-selected
      req(rv_in$meta)
      cfg <- isolate(rv_in$config)

      curr_sel <- isolate(input$CB_sub_innerA)
      subs <- strsplit(cfg[UI == input$SI_groupA]$fID, "\\|")[[1]]
      subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
      names(subs2) <- seq_along(subs2)
      subs[as.numeric(names(sort(subs2)))]

      subs_label <- paste0("select ",isolate(input$SI_groupA),"s: ")

      freezeReactiveValue(input, "CB_sub_innerA")
      updateCheckboxGroupInput(session,
                               inputId = "CB_sub_innerA",
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
      req(input$CB_sub_innerA,
          input$SI_groupA)
      out_params$observ_grpA = input$SI_groupA
      out_params$observ_subselA = input$CB_sub_innerA

      # # aggregate var
      # out_params$feat_grp <- input$SI_feat_grp
      # out_params$feat_subsel <- input$SI_feat_subsel

      out_params$observ_x <- input$SI_x_info
      out_params$observ_y_raw <- input$SI_obs_raw
      # out_params$observ_y_rcomp <- input$SI_obs_comp
      #
      out_params$obs_type <- input$RB_x_obs_var


    })

#
#   observe(input$observ_y_raw)({
#
#     # TODO
#     # set whether we are updating "raw" or "comp" data for visualization based on selected Y
#
#     out_params$measure_type <- c("raw")
#
#   })


    # keep these updated
    observe({
      #out_params[["exp_fact"]] <- input$SI_exp_fact_select #NULL
      out_params[["omics_list"]] <-    omics_list  # value & viz_now
      # out_params[["raw_plot_type"]] <-    input$RB_raw_plot_type
      # out_params[["comp_plot_type"]] <-    input$RB_comp_plot_type
    })

  return(out_params)

  })
}

## To be copied in the UI
# mod_side_selector_ui("side_selector_ui_1")

## To be copied in the server
# mod_side_selector_server("side_selector_ui_1")

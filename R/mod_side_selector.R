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
    helpText(
       HTML("Choose a Dataset. Change visualization by selecting muscle-type and condition. Proteins can be selected directly from the <i>Visualize</i> tab,  from <i>Browse Table</i> tab, or by typing the protein name into the <i>Selected genes</i> window below.")
       ),
     hr(style = "border-top: 1px solid #000000;"),
    #TODO:  change these to render prper HTML

    fluidRow(
        htmlOutput( ns( "ui_curr_database" ))
      ),
    fluidRow(
      uiOutput( ns( "ui_DIV_warn" ))
      ),
    fluidRow(
      htmlOutput( ns( "ui_db_type" ))
      ), #fluidRow 1a

    fluidRow(

      column(
        width=10,
        offset=0,
        h5("Sample/Omics- information"),
        "Choose cells/samples, grouping/subsetting variables, observables, and omic features ",
        "to visualize in the playground",
        br(),br(),
        selectizeInput(ns("SI_x_info"), "Annotation information (X-axis):", choices=NULL), # options = list(placeholder = ""))

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

      column(
        width = 8,
        offset=3,
        selectizeInput(ns("SI_obs_comp"), "Comparatives", choices=NULL)
      )

    ),


    # fluidRow(
    #       col_8(
    #         selectizeInput(ns("SI_exp_fact_select"), "choose experimental factor: ", "",
    #                   options = list(placeholder = "load database first"))
    #         )
    # ), #fluidRow 1b
    #########################################
    # subset
    #########################################
    fluidRow(
      column(width = 5,
        offset = 1,
        #actionButton(ns("AB_subset_tog"), "Toggle subset observations"), #TODO: change text when toggled
        shinyjs::disabled(selectizeInput(ns("SI_subset"), "Obs information to subset:", "",
                                         multiple = FALSE, options = list(placeholder = "choose dataset first"))),
      )
    ),
    fluidRow(
      column(width = 10,
        offset = 1,
        uiOutput(ns("ui_subset"))
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



    # fluidRow(
    #   # col_3(offset=1,
    #   #       shinyjs::hidden(textOutput(ns("ui_group_factor"))),
    #   #       shinyjs::disabled(selectizeInput(ns("SI_aux_group_select"),  "group values", "",
    #   #                                        multiple=TRUE,
    #   #                                        options = list(placeholder = "load db/choose aux factor "))),
    #   # ),
    #   column(5,offset=1,
    #         shinyjs::hidden(textOutput(ns("ui_raw_factor"))),
    #         shinyjs::disabled(selectizeInput(ns("SI_aux_raw_select"),  "raw measurments", "",
    #                                          multiple=TRUE,
    #                                          options = list(placeholder = ""))),
    #   ),
    #   column(5,offset=0,
    #         shinyjs::hidden(textOutput(ns("ui_comp_factor"))),
    #         shinyjs::disabled(selectizeInput(ns("SI_aux_comp_select"),  "comparative measures", "",
    #                                          multiple=TRUE,
    #                                          options = list(placeholder = "")))
    #   )
    # ),
#
#     fluidRow(
#       # col_3(offset=1,
#       #       shinyjs::hidden(textOutput(ns("ui_group_factor1")))
#       #     ),
#       column(
#             width = 3,
#             offset=4,
#             shinyjs::hidden(textOutput(ns("ui_raw_factor1")))
#       ),
#       column(
#             width = 3,
#             offset=1,
#             shinyjs::hidden(textOutput(ns("ui_comp_factor1")))
#       )
#     ), #fluidRow 1c

    #


#optCrt="
    mod_omic_selector_ui(ns("omic_selector_ui_1")),

    fluidRow(
      uiOutput(ns("ui_text_warn"), width = "100%"),
        ) #fluidRow

# # TODO:  make this conditional on what kind of data is loaded.. and just a single plot type
#     fluidRow(
#       column(
#            width=6,
#            style="border-right: 2px solid black",
#            h4("quantaties"),
#            fluidRow(
#              radioButtons(ns("RB_raw_plot_type"), "Plot type",
#                            choices = quant_plot_types,
#                            selected = quant_plot_types[3], inline = TRUE )
#              )
#       ),
#       column(
#         width=6,
#         h4("comparisons"),
#         fluidRow(
#           radioButtons(ns("RB_comp_plot_type"), "Plot type",
#                        choices = comp_plot_types,
#                        selected = comp_plot_types[1], inline = TRUE )
#           )
#       )
#     ) #fluidRow 4
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

      observ_grp = NULL,
      observ_subsel = NULL,

      #
      observ_x = NULL,
      observ_y_raw = NULL,
      observ_y_comp = NULL,

      # aggregate var
      feat_grp = NULL,
      feat_subsel = NULL,

      # until this is set we don't actually plot ayhting.
      measure_type = NULL, #"raw" or "comp"
      # TODO:  make this conditional on what kind of data is loaded.. and just a single plot type
      raw_plot_type = NULL,
      comp_plot_type = NULL

      # out_params$aux_raw <- input$SI_obs_raw
      # out_params$aux_comp <- input$SI_obs_comp
    )


    # omics_for_selector <- reactive(rv_in$omics_feature)
    # omics_list <- mod_omic_selector_server("omic_selector_ui_1", omics_for_selector)
    omics_out <- reactive( rv_in$omics_feature )

    omics_list <- mod_omic_selector_server("omic_selector_ui_1", omics_out )


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

        shinyjs::enable("SI_subset")
        freezeReactiveValue(input, "SI_subset")
        updateSelectizeInput(session, "SI_subset","Obs information to subset:",
                             choices = rv_in$config$meta[grp == TRUE]$UI,
                             selected = rv_in$default$grp2,  server = TRUE)


        shinyjs::enable("CB_sub_all")
        shinyjs::enable("CB_sub_none")

        shinyjs::enable("SI_obs_raw")
        shinyjs::enable("SI_obs_comp")

        # update omics_out
        omics_out <- rv_in$omics_feature




        # options = list(
        #   maxOptions = length(Van_conf[is.na(fID)]$UI) + 3,
        #   create = TRUE, persist = TRUE, render = I(optCrt)))

        # # Update selectInput according to dataset
        # if (!is.null(rv_in$config)) {
        #   print("enabled subset ")
        #   shinyjs::enable("SI_subset")
        #   freezeReactiveValue(input, "SI_subset")
        #   updateSelectizeInput(session, "SI_subset","Obs information to subset:",
        #                        choices = rv_in$config$meta[grp == TRUE]$UI,
        #                        selected = rv_in$default$grp2,  server = TRUE)
        #
        #   shinyjs::enable("CB_sub_all")
        #   shinyjs::enable("CB_sub_none")
        # } else {
        #   print("disabled subset ")
        #   shinyjs::disable("SI_subset")
        #   shinyjs::disable("CB_sub_all")
        #   shinyjs::disable("CB_sub_none")
        # }


      },
    ignoreNULL = TRUE,
    ignoreInit = TRUE
    ) #observe event




    output$ui_subset <- renderUI({
      #req(input$SI_subset) # do i need this?
      if (input$SI_subset == ""){
        print("loading....")
        return()
      }
      cfg <- isolate(rv_in$config$meta)
      subs <- strsplit(cfg[UI == input$SI_subset]$fID, "\\|")[[1]]
      subs <- sort(subs)
      subs_label <- paste0("select ",isolate(input$SI_subset),"s: ")

      checkboxGroupInput( ns("CB_sub_inner1"),
                         label = subs_label,
                         inline = TRUE,
                         choices = subs,
                         selected = subs)
    })


    observeEvent(input$CB_sub_all, {
      #req(input$CB_sub_inner1) this makes the button useless if everything is deselected
      req(rv_in$meta)
      cfg <- isolate(rv_in$config$meta)

      subs <- strsplit(cfg[UI == input$SI_subset]$fID, "\\|")[[1]]
      subs <- sort(subs)
      subs_label <- paste0("select ",isolate(input$SI_subset),"s: ")

      freezeReactiveValue(input, "CB_sub_inner1")
      updateCheckboxGroupInput(session,
                               inputId = "CB_sub_inner1",
                               label = subs_label,
                               choices = subs,
                               selected = subs,
                               inline = TRUE) # WARNING.  make sure subset is null is checked (length(0?))
    })


    observeEvent(input$CB_sub_none, {
      req(input$CB_sub_inner1) #short circuits if already de-selected
      req(rv_in$meta)
      cfg <- isolate(rv_in$config$meta)

      curr_sel <- isolate(input$CB_sub_inner1)
      subs <- strsplit(cfg[UI == input$SI_subset]$fID, "\\|")[[1]]
      subs <- sort(subs)
      subs_label <- paste0("select ",isolate(input$SI_subset),"s: ")

      freezeReactiveValue(input, "CB_sub_inner1")
      updateCheckboxGroupInput(session,
                               inputId = "CB_sub_inner1",
                               label = subs_label,
                               choices = subs,
                               selected = character(0), # character(0), #, "" , NULL
                               inline = TRUE) # WARNING.  make sure subset is null is checked (length(0?))
    })


    observe({
      req(rv_in$config)
      raw_choices = rv_in$config$mat[observ == TRUE & ID!="raw"]$fIDloc
      #TODO:  make the choices *named* paste0(rv_in$config$fID, fUI)
      # dat_source = rv_in$config$mat[fID == rv_in$aux_raw]$ID
      # paste(raw_choices, dat_sourcce)

      if (length(raw_choices)>0) {
        freezeReactiveValue(input, "SI_obs_raw")
        updateSelectizeInput(session, "SI_obs_raw","Raw vals:",
                                                    choices = raw_choices,
                                                    selected = raw_choices[1],  server = TRUE)
      } else {
        print("disabled  raw observations ")
        shinyjs::disable("SI_obs_raw")
        # change back to placeholder??
        freezeReactiveValue(input, "SI_obs_raw")
        updateSelectizeInput(session, "SI_obs_raw", "Raw vals: ", "", options = list(placeholder = ""))
      }

    })

    observe({

      req(rv_in$config)
      comp_choices = rv_in$config$mat[comp == TRUE ]$fIDloc
      if (length(comp_choices)>0) {
        freezeReactiveValue(input, "SI_obs_comp")
        updateSelectizeInput(session, "SI_obs_comp","Comparatives:",
                             choices = comp_choices,
                             selected = comp_choices[1],  server = TRUE)
      } else {
        print("disabled  comparative observations ")
        shinyjs::disable("SI_obs_comp")
        freezeReactiveValue(input, "SI_obs_comp")
        updateSelectizeInput(session, "SI_obs_comp", "Comparatives:", "",
                        options = list(placeholder = "") )
      }


   })




    observe({
        req(rv_in$config)
        updateSelectizeInput(session, "SI_x_info", server = TRUE,
                           choices = rv_in$config$meta[grp == TRUE]$UI,
                           selected = rv_in$default$grp2)
    })

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
      req(input$CB_sub_inner1,
          input$SI_subset)
      out_params$observ_grp = input$SI_subset
      out_params$observ_subsel = input$CB_sub_inner1

      # # aggregate var
      # out_params$feat_grp <- input$SI_feat_grp
      # out_params$feat_subsel <- input$SI_feat_subsel

      out_params$observ_x <- input$SI_x_info
      out_params$observ_y_raw <- input$SI_obs_raw
      out_params$observ_y_rcomp <- input$SI_obs_comp


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
      out_params[["raw_plot_type"]] <-    input$RB_raw_plot_type
      out_params[["comp_plot_type"]] <-    input$RB_comp_plot_type
    })

  return(out_params)

  })
}

## To be copied in the UI
# mod_side_selector_ui("side_selector_ui_1")

## To be copied in the server
# mod_side_selector_server("side_selector_ui_1")

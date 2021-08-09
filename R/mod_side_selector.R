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
    fluidRow(textOutput(ns("ui_curr_database"))),
    fluidRow(uiOutput(ns("ui_DIV_warn"))),
    fluidRow(textOutput(ns("ui_exp_fact_name"))), #fluidRow 1a

    # fluidRow(
    #       col_8(
    #         selectizeInput(ns("SI_exp_fact_select"), "choose experimental factor: ", "",
    #                   options = list(placeholder = "load database first"))
    #         )
    # ), #fluidRow 1b
    fluidRow(
      column(width = 6,
        offset = 1,
        #actionButton(ns("AB_subset_tog"), "Toggle subset observations"), #TODO: change text when toggled
        shinyjs::disabled(selectizeInput(ns("SI_subset"), "Obs information to subset:", "",
                                         multiple = FALSE, options = list(placeholder = "choose dataset first"))),
      )
      #,  # not implimented ... need to generate var annotations groups
      # column(width = 4,
      #        offset = 1,
      #        #actionButton(ns("AB_subset_tog"), "Toggle subset observations"), #TODO: change text when toggled
      #        shinyjs::disabled(selectizeInput(ns("SI_var_subset"), "omic subset?:", "",
      #                                         multiple = FALSE, options = list(placeholder = "choose dataset first"))),
      # )
      ),
    fluidRow(
      column(width = 10,
        offset = 1,
        uiOutput(ns("ui_subset"))
        )
      ),
    fluidRow(
      column(width = 5,shinyjs::disabled(actionButton(ns("CB_sub_all"), "Select all groups", class = "btn btn-primary"))),
      column(width = 5,
        offset=1,
        shinyjs::disabled(actionButton(ns("CB_sub_none"), "Deselect all groups", class = "btn btn-primary") )
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

    # fluidRow(
    #   selectizeInput(
    #     ns("SI_omics_select"), "Select omic features", "",
    #     multiple=TRUE, options = list(placeholder = "Choose omic feature (i.e. genes,proteins,lipids...)")
    #   )
    # ), #fluidRow 2
    #
    # fluidRow(
    #   col_4(offset = 0, style='padding-left:0px; padding-right:1px',
    #          actionButton(ns("AB_omics_submit"),"Submit",class="hidableSubmit")
    #   ),
    #   col_4( offset = 2, style='padding-left:0px; padding-right:1px',
    #          actionButton(ns("AB_omics_reset"), "Clear",class="hidableClear")
    #   )
    #   ),

    mod_omic_selector_ui(ns("omic_selector_ui_1")),

    fluidRow(
      uiOutput(ns("ui_text_warn"), width = "100%"),
        ), #fluidRow
    fluidRow(
      column(
           width=6,
           style="border-right: 2px solid black",
           h4("quantaties"),
           fluidRow(
             radioButtons(ns("RB_raw_plot_type"), "Plot type",
                           choices = quant_plot_types,
                           selected = quant_plot_types[3], inline = TRUE )
             )
      ),
      column(
        width=6,
        h4("comparisons"),
        fluidRow(
          radioButtons(ns("RB_comp_plot_type"), "Plot type",
                       choices = comp_plot_types,
                       selected = comp_plot_types[1], inline = TRUE )
          )
      )
    ) #fluidRow 4
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
    out_params <-reactiveValues(

      #omics_names = NULL,
      omics_list = NULL,
      raw_plot_type = NULL,
      comp_plot_type = NULL,
      observ_grp = NULL,
      observ_subsel = NULL,

      # get ride of these
      exp_fact = NULL,
      aux_fact = NULL,
      comp_dat = NULL,
      obs_dat = NULL

    )

    # omics_for_selector <- reactive(rv_in$omics_feature)
    # omics_list <- mod_omic_selector_server("omic_selector_ui_1", omics_for_selector)
    #
    omics_out <- reactive( rv_in$omics_feature )

    omics_list <- mod_omic_selector_server("omic_selector_ui_1", omics_out )


    output$ui_curr_database <- renderPrint({
      if (is.null(rv_in$database_name)) {
        print("no datbase loaded")
        } else {
        print(paste("Current database is: ", rv_in$database_name))
        }
      })

    # Warning if no data loaded
    output$ui_DIV_warn <- renderUI( {
      if (is.null(rv_in$database_name)) {
        div(
          tags$br(),
          span(class = "warn", "No dataset loaded"),
          print("warning:::!")
          ) }
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

        # update omics_out
        omics_out <- rv_in$omics_feature

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


    # show the factors that have been loaded
    output$ui_exp_fact_name <- renderText({
      req(input$SI_subset)
      print(paste("current experimental factor: ", input$SI_subset))
    })


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


    # TODO:  wrap this in a module..
    #
    #

    # ############################+
    # ## gene selection:
    # ############################+
    # # keep track of which proteins have been selected
    # omics_list <- reactiveValues(value=character(0), viz_now=FALSE)
    #
    # # Error handling: check if too many omics selected
    # # NOT within_limit  --> print error statement
    # # within_limit --> Update the protein selection
    # # THIS GETS RUN ON EVERY CLICK!!
    # max_omic_feats <- 100
    #
    # observe({
    #   if ( length(unique(omics_list$value ) ) >= max_omic_feats ) {
    #     # will this work?  or is the "observing" affecting a reactive with this render?
    #     output$ui_text_warn <- renderUI({
    #       tags$div(class = "warning", checked = NA,
    #                HTML(
    #                  paste('
    #                   </head>
    #                   <style>
    #                   .warning div {text-align: left; padding: 50px 70px 100px;}
    #                   </style>
    #                   </head>
    #
    #                   <body>
    #                   <div>
    #                   Omic-selection limit reached, choose fewer items for
    #                   faster computation ! <br> Continue by pressing "Clear".
    #                   </div>
    #                   </body>')))
    #     })
    #
    #
    #   } else {
    #     #omics_choice_list <- rv_in$var[[rv_in$omics_feature]]
    #     omics_choice_list <- isolate(rv_in$omics_feature)
    #     #if ( !is.null(omics_choice_list$var) & !is.null(rv_in$omics_feature) ) {
    #     if ( is.null(omics_choice_list)  ) {
    #       omics_choice_list <- "" #rownames(rv_in$var)
    #     } else {
    #       omics_choice_list <- names(omics_choice_list)
    #     }
    #
    #     omics_choices <- isolate(omics_list$value)
    #     freezeReactiveValue(input, "SI_omics_select")
    #     updateSelectizeInput(session, "SI_omics_select",
    #                          choices = isolate(omics_choice_list),
    #                          selected = omics_choices, server=TRUE)
    #     # # DEBUG
    #     # print("   ::updated omics list::  ")
    #     # print(isolate(omics_list$value))
    #   }
    #   #print(paste0("current number of selected omics:  ",length(unique(isolate(omics_list$value)))))
    #
    # })
    #
    #
    # ############################+
    # ## "reset" and "submit" simply sets the viz_now flag
    # ############################+
    # ############################
    # observe({ # turn on if the "placeholder" is gone
    #   shinyjs::toggleState("AB_omics_submit", !all(input$SI_omics_select == "Choose omic feature (i.e. genes,proteins,lipids...)"))
    # })
    #
    # observeEvent(input$AB_omics_reset, {
    #   omics_list$viz_now = FALSE
    #   omics_list$value <- character(0)
    #   # get rid of error message when resetting selection
    #   output$ui_text_warn <- renderUI({ })
    # })
    #
    # observeEvent(input$AB_omics_submit, {
    #   if(length(unique(omics_list$value)) < max_omic_feats ) { #defensive
    #     omics_list$value <- input$SI_omics_select #include direct selection from protein-box when pressing submit
    #     omics_list$viz_now = TRUE
    #   } # else ?? print warning??
    # })


    observeEvent( omics_list$viz_now, {
      # route the chosen data type out...
      req(input$CB_sub_inner1,
          input$SI_subset)
      out_params$observ_grp = input$SI_subset
      out_params$observ_subsel = input$CB_sub_inner1

    })
      #
      # if (!is.null(rv_in$aux_raw)){
      #   # could also extract from the name instead of lookkup
      #   #       subs <- strsplit(rv_in$aux_raw, "\\|")[[1]]
      #   # rv_in$ad[[subs[1]]][[subs[2]]]
      #   dat_source = rv_in$config$mat[fIDloc == rv_in$aux_raw]$ID
      #   dat_key = rv_in$config$mat[fIDloc == rv_in$aux_raw]$fID
      #   #dd <- isolate(rv_in$ad[[rv_in$config$mat[fID == rv_in$aux_raw]$ID]][[rv_in$aux_raw]])
      #   if (dat_source == "obs") {
      #     aux_raw_names <- isolate(rv_in$ad$obs[[dat_key]])
      #   } else if (dat_source == "var") {
      #     aux_raw_names <- isolate(rv_in$ad$var[[dat_key]])
      #   } else {
      #     aux_raw_names <- NULL
      #   }
      #   out_params$obs_dat = aux_raw_names
      # }
      #
      # # comparatives are packed in obsm and varm
      # if (!is.null(rv_in$aux_comp)){
      #   dat_source = rv_in$config$mat[fIDloc == rv_in$aux_comp]$ID
      #   dat_key = rv_in$config$mat[fIDloc == rv_in$aux_comp]$fID
      #   if (dat_source == "obsm") {
      #     aux_comp_names <- isolate(rv_in$ad$uns[[rv_in$aux_comp]])
      #   } else if (dat_source == "varm") {
      #     aux_comp_names <- isolate(rv_in$ad$uns[[rv_in$aux_comp]])
      #   } else {
      #     aux_comp_names <- NULL
      #
      #   }
      #   out_params$comp_dat <- aux_comp_names
      # }
      #
      #
      #
      #


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

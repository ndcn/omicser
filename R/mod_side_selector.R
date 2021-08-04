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

  plot_types <- c("Volcano Plot"="volcano", "Dot Plot"="dotplot","Heat Map"="hmap","Box Plot"="box")

  selector_tags <- tagList(
    helpText(HTML("Choose a Dataset.
                          Change visualization by selecting muscle-type and condition.
                          Proteins can be selected directly from the <i>Visualize</i> tab, from <i>Browse Table</i> tab, or by typing the protein name into the <i>Selected genes</i> window below.")
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
        col_6(
        ofset = 1,
        #actionButton(ns("AB_subset_tog"), "Toggle subset observations"), #TODO: change text when toggled
        shinyjs::disabled(selectizeInput(ns("SI_subset"), "Obs information to subset:", "",
                                         multiple = FALSE, options = list(placeholder = "choose dataset first"))),

        uiOutput(ns("ui_subset")),
        shinyjs::disabled(actionButton(ns("CB_sub_all"), "Select all groups", class = "btn btn-primary")),
        shinyjs::disabled(actionButton(ns("CB_sub_none"), "Deselect all groups", class = "btn btn-primary") )
      )
    ),

    fluidRow(
      # col_3(offset=1,
      #       shinyjs::hidden(textOutput(ns("ui_group_factor"))),
      #       shinyjs::disabled(selectizeInput(ns("SI_aux_group_select"),  "group values", "",
      #                                        multiple=TRUE,
      #                                        options = list(placeholder = "load db/choose aux factor "))),
      # ),
      col_3(offset=4,
            shinyjs::hidden(textOutput(ns("ui_raw_factor"))),
            shinyjs::disabled(selectizeInput(ns("SI_aux_raw_select"),  "raw measurments", "",
                                             multiple=TRUE,
                                             options = list(placeholder = "load db/choose aux factor "))),
      ),
      col_3(offset=0,
            shinyjs::hidden(textOutput(ns("ui_comp_factor"))),
            shinyjs::disabled(selectizeInput(ns("SI_aux_comp_select"),  "comparative measures", "",
                                             multiple=TRUE,
                                             options = list(placeholder = "load db/choose aux factor ")))
      )
    ),

    fluidRow(
      # col_3(offset=1,
      #       shinyjs::hidden(textOutput(ns("ui_group_factor1")))
      #     ),
      col_3(offset=4,
            shinyjs::hidden(textOutput(ns("ui_raw_factor1")))
      ),
      col_3(offset=1,
            shinyjs::hidden(textOutput(ns("ui_comp_factor1")))
      )
    ), #fluidRow 1c

    fluidRow(
      selectizeInput(
        ns("SI_omics_select"), "Select omic features", "",
        multiple=TRUE, options = list(placeholder = "Choose omic feature (i.e. genes,proteins,lipids...)")
      )
    ), #fluidRow 2

    fluidRow(
      col_4(offset = 0, style='padding-left:0px; padding-right:1px',
             actionButton(ns("AB_omics_submit"),"Submit",class="hidableSubmit")
      ),
      col_4( offset = 2, style='padding-left:0px; padding-right:1px',
             actionButton(ns("AB_omics_reset"), "Clear",class="hidableClear")
      )
      ),
    fluidRow(
      uiOutput("ui_text_warn", width = "100%"),
    ), #fluidRow 3

    fluidRow(
      radioButtons(ns("RB_plot_type"), "Plot type",
                   choices = plot_types,
                   selected = plot_types[3], inline = TRUE )
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
      exp_fact = NULL,
      aux_fact = NULL,
      omics_names = NULL,
      omics_list = NULL,
      plot_type = NULL,
      comp_dat = NULL,
      obs_dat = NULL
    )


    output$ui_curr_database <- renderPrint({
      if (is.null(rv_in$database_name)) {
        print("no datbase loaded")
        } else {
        print(paste("Current database is: ", rv_in$database_name))
        }
      })

    # Warning if no data loaded
    output$ui_DIV_warn <- renderUI( {
      if (!is.null(rv_in$database_name)) {
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
      rv_in$trigger,
      {
        req(rv_in$ad)
        print("TRIGGERED -------->")
        browser()
        exp_fact_name <- isolate(rv_in$exp_factor)
        # show the factors that have been loaded
        output$ui_exp_fact_name <- renderText({
          req(rv_in$exp_factor)
          print(paste("current experimental factor: ", exp_fact_name))

        })
        # #print(isolate(rv_in$factors0))
        # exp_factor_vals <- isolate(rv_in$ad$obs[[rv_in$exp_factor]])
        # updateSelectizeInput(session, "SI_exp_fact_select", paste0(exp_fact_name," (exp fact):"), # wrap in ns?
        #                      choices = exp_factor_vals,
        #                      selected = exp_factor_vals[1], server=TRUE)
        #
        # reset trigger for next time (if we don't reset, will everything just be responsive?)
        rv_in$trigger  <-  0

        # Update selectInput according to dataset

        if (!is.null(rv_in$config)) {
          print("enabled subset ")
          shinyjs::enable("SI_subset")
          updateSelectizeInput(session, "SI_subset","Obs information to subset:",
                               choices = rv_in$config$meta[grp == TRUE]$UI,
                               selected = rv_in$default$grp1,  server = TRUE)

          shinyjs::enable("CB_sub_all")
          shinyjs::enable("CB_sub_none")
        } else {
          print("disabled subset ")
          shinyjs::disable("SI_subset")
          shinyjs::disable("CB_sub_all")
          shinyjs::disable("CB_sub_none")
        }


      },
    ignoreNULL = TRUE,
    ignoreInit = TRUE
    ) #observe event



    output$ui_subset <- renderUI({
      if (input$SI_subset == ""){
        print("loading....")
      } else {
      sub = strsplit(rv_in$config$meta[UI == input$SI_subset]$fID, "\\|")[[1]]
      checkboxGroupInput("CB_sub_inner1", "Select which groups to show", inline = TRUE,
                         choices = sub, selected = sub)
      }
    })
    observeEvent(input$CB_sub_all, {
      sub = strsplit(rv_in$config$meta[UI == input$SI_subset]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(session, inputId = "CB_sub_inner1", label = "Select which groups to show",
                               choices = sub, selected = NULL, inline = TRUE) # WARNING.  make sure subset is null is checked (length(0?))
    })
    observeEvent(input$CB_sub_none, {
      sub = strsplit(rv_in$config$meta[UI == input$SI_subset]$fID, "\\|")[[1]]
      updateCheckboxGroupInput(session, inputId = "CB_sub_inner1", label = "Select which groups to show",
                               choices = sub, selected = sub, inline = TRUE)
    })

    #
    ############################+
    ## gene selection:
    ############################+
    # keep track of which proteins have been selected
    omics_list <- reactiveValues(value=character(0), viz_now=FALSE)

    # Error handling: check if too many proteins selected
    # NOT within_limit  --> print error statement
    # within_limit --> Update the protein selection
    # THIS GETS RUN ON EVERY CLICK!!@
    max_omic_feats <- 100

    observe({
      if ( length(unique(omics_list$value ) ) >= max_omic_feats ) {
        # will this work?  or is the "observing" affecting a reactive with this render?
        output$ui_text_warn <- renderUI({
          tags$div(class = "warning", checked = NA,
                   HTML(
                     paste('
                      </head>
                      <style>
                      .warning div {text-align: left; padding: 50px 70px 100px;}
                      </style>
                      </head>

                      <body>
                      <div>
                      Omic-selection limit reached, choose fewer items for
                      faster computation ! <br> Continue by pressing "Clear".
                      </div>
                      </body>')))
        })


      } else {
        #omics_choice_list <- rv_in$var[[rv_in$omics_feature]]
        omics_choice_list <- isolate(rv_in$omics_feature)
        #if ( !is.null(omics_choice_list$var) & !is.null(rv_in$omics_feature) ) {
        if ( is.null(omics_choice_list)  ) {
          omics_choice_list <- "" #rownames(rv_in$var)
        } else {
          omics_choice_list <- names(omics_choice_list)
        }

        omics_choices <- isolate(omics_list$value)
        updateSelectizeInput(session, "SI_omics_select",
                             choices = isolate(omics_choice_list),
                             selected = omics_choices, server=TRUE)
        # # DEBUG
        # print("   ::updated omics list::  ")
        # print(isolate(omics_list$value))
      }
      #print(paste0("current number of selected omics:  ",length(unique(isolate(omics_list$value)))))

    })


    ############################+
    ## "reset" and "submit" simply sets the viz_now flag
    ############################+
    ############################
    observe({ # turn on if the "placeholder" is gone
      shinyjs::toggleState("AB_omics_submit", !all(input$SI_omics_select == "Choose omic feature (i.e. genes,proteins,lipids...)"))
    })

    observeEvent(input$AB_omics_reset, {
      omics_list$viz_now = FALSE
      omics_list$value <- character(0)
      # get rid of error message when resetting selection
      output$ui_text_warn <- renderUI({ })
    })

    observeEvent(input$AB_omics_submit, {
      if(length(unique(omics_list$value)) < max_omic_feats ) { #defensive
        omics_list$value <- input$SI_omics_select #include direct selection from protein-box when pressing submit
        omics_list$viz_now = TRUE
      }

      # route the chosen data type out...
      if (!is.null(rv_in$aux_raw)){
        dat_source = rv_in$config$mat[fID == rv_in$aux_raw]$ID
        #dd <- isolate(rv_in$ad[[rv_in$config$mat[fID == rv_in$aux_raw]$ID]][[rv_in$aux_raw]])
        if (dat_source == "obs") {
          aux_raw_names <- isolate(rv_in$ad$obs[[rv_in$aux_raw]])
        } else if (dat_source == "var") {
          aux_raw_names <- isolate(rv_in$ad$var[[rv_in$aux_raw]])
        }
        out_params$obs_dat = aux_raw_names
      }

      if (!is.null(rv_in$aux_comp)){
        dat_source = rv_in$config$mat[fID == rv_in$aux_raw]$ID
        if (dat_source == "obsm") {
          aux_comp_names <- isolate(rv_in$ad$uns[[rv_in$aux_raw]])
        } else if (dat_source == "varm") {
          aux_comp_names <- isolate(rv_in$ad$uns[[rv_in$aux_raw]])
        }

        out_params$comp_dat = aux_comp_names
      }
    })

    observe({
      out_params[["exp_fact"]] <- input$SI_exp_fact_select
      out_params[["aux_fact"]] <-  input$SI_aux_fact_select
      out_params[["omics_names"]] <-    input$SI_omics_select
      out_params[["omics_list"]] <-    omics_list  # value & viz_now
      out_params[["plot_type"]] <-    input$RB_plot_type
    })

  return(out_params)

  })
}

## To be copied in the UI
# mod_side_selector_ui("side_selector_ui_1")

## To be copied in the server
# mod_side_selector_server("side_selector_ui_1")

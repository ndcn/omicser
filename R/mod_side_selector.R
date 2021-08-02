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

    fluidRow(
          col_8(
            selectizeInput(ns("SI_exp_fact_select"), "choose experimental factor: ", "",
                      options = list(placeholder = "load database first"))
            )
    ), #fluidRow 1b
    fluidRow(
      col_3(offset=1,
            shinyjs::hidden(textOutput(ns("ui_group_factor"))),
            shinyjs::disabled(selectizeInput(ns("SI_aux_group_select"),  "group values", "",
                                             multiple=TRUE,
                                             options = list(placeholder = "load db/choose aux factor "))),
      ),
      col_3(offset=0,
            shinyjs::hidden(textOutput(ns("ui_raw_factor"))),
            shinyjs::disabled(selectizeInput(ns("SI_aux_raw_select"),  "raw measurments", "",
                                             multiple=TRUE,
                                             options = list(placeholder = "load db/choose aux factor "))),
      ),
      col_3(offset=0,
            shinyjs::hidden(textOutput(ns("ui_comp_factor"))),
            shinyjs::disabled(selectizeInput(ns("SI_aux_comp_select"),  "comparative1 measures", "",
                                             multiple=TRUE,
                                             options = list(placeholder = "load db/choose aux factor ")))
      )
    ),

    fluidRow(
      col_3(offset=1,
            shinyjs::hidden(textOutput(ns("ui_group_factor1")))
          ),
      col_3(offset=1,
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
      plot_type = NULL
    )


    ############################ +
    ## rv_in:  input values shared by modules
    ############################ +
   #      database_name = NULL,
   #      omics_type    = NULL,   # Transcript-omics, Prote-omics, Lipid-omics, Metabol-omics, misc-
   #
   #      ad = NULL,   # this is the DATA  packed into an anndata object
   #
   #                              # these reflect the choices made in the ingestor
   #                              #    omics key feature i.e. genes, proteins, lipids
   #                              #  primary & auxilarry
   #
   #      omics_feature = NULL,
   #      aux_features = NULL,
   #      exp_factor = NULL,        # observables:  experiemntal factor... i.e. patient/control, old/young,
   #
   #
   #0      exp_factor = NULL,
  #  aux_group = NULL,
  #  aux_comp = NULL,
  #  aux_raw = NULL,


   #      defaults = NULL,
   #      configs = NULL,
   #      meta = NULL,
   #
   #      trigger = 0
   # #



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
      rv_in$trigger,
      {
        req(rv_in$ad)
        print("TRIGGERED -------->")

        exp_fact_name <- isolate(rv_in$exp_factor)
        # show the factors that have been loaded
        output$ui_exp_fact_name <- renderText({
          req(rv_in$exp_factor)
          print(paste("current experimental factor: ", exp_fact_name))

        })
        #print(isolate(rv_in$factors0))
        exp_factor_vals <- isolate(rv_in$ad$obs[[rv_in$exp_factor]])
        updateSelectizeInput(session, "SI_exp_fact_select", paste0(exp_fact_name," (exp fact):"), # wrap in ns?
                             choices = exp_factor_vals,
                             selected = exp_factor_vals[1], server=TRUE)
        # reset trigger for next time (if we don't reset, will everything just be responsive?)
        rv_in$trigger  <-  0


        if (!is.null(rv_in$aux_group)){
          aux_group_names <-  unique(isolate(rv_in$obs[[rv_in$aux_group]]))
          # this should also be backed into conf or def
          output$ui_group_factor <- renderText({
            req(rv_in$aux_group)
            print(paste("grouping factors: ", paste0( aux_group_names,collapse="; ")))
          })
          updateSelectizeInput(session, "SI_aux_group_select",
                               choices = aux_group_names,
                               selected = aux_group_names[1], server=TRUE)
          shinyjs::enable("SI_aux_group_select")
        }

        if (!is.null(rv_in$aux_raw)){
          aux_raw_names <- isolate(rv_in$ad$uns[[rv_in$aux_raw]])
          output$ui_raw_factor <- renderText({
            req(rv_in$aux_raw)
            print(paste("raw measures: ", paste0( aux_raw_names,collapse="; ")))
          })
          updateSelectizeInput(session, "SI_aux_raw_select",
                               choices = aux_raw_names,
                               selected = aux_raw_names[1], server=TRUE)
          shinyjs::enable("SI_aux_raw_select")

        }

        if (!is.null(rv_in$aux_comp)){
          aux_comp_names <- isolate(rv_in$ad$uns[[rv_in$aux_comp]])
          output$ui_comp_factor <- renderText({
            req(rv_in$aux_comp)
            print(paste("comparative measures: ", paste0( aux_comp_names,collapse="; ")))
          })
          updateSelectizeInput(session, "SI_aux_comp_select",
                               choices = aux_comp_names,
                               selected = aux_comp_names[1], server=TRUE)
          shinyjs::enable("SI_aux_comp_select")
          }

      # # reset omics list
      # print("reset omics_list!")
      # omics_list$viz_now = FALSE
      # omics_list$value <- character(0)
      # # get rid of error message when resetting selection
      # output$textWarn <- renderUI({ })
      # output$textWarn_comp <- renderUI({ })
      #
      # do i need to wrap this in a reactive?
      },
    ignoreNULL = TRUE,
    ignoreInit = TRUE
    ) #observe event


    observe({
      req(rv_in$aux_group)
      req(input$SI_aux_group_select)
      if (input$SI_aux_group_select[1] != "") { #|| ( !is.na(input$SI_var0) )
        shinyjs::show("ui_group_factor1")
      } else {
        shinyjs::hide("ui_group_factor1")
      }
    })
    output$ui_group_factor1 <- renderText({
      req(rv_in$aux_group)
      aux_factor1_name <- isolate(input$SI_aux_group_select)
      print( paste0("grouping factors selected: ", paste( aux_factor1_name,collapse="; ") ))
    })


    observe({
      req(rv_in$aux_raw)
      req(input$SI_aux_raw_select)
      if (input$SI_aux_raw_select[1] != "") { #|| ( !is.na(input$SI_var0) )
        shinyjs::show("ui_raw_factor1")
      } else {
        shinyjs::hide("ui_raw_factor1")
      }
    })
    output$ui_raw_factor1 <- renderText({
      req(rv_in$aux_raw)
      aux_raw1_name <- isolate(input$SI_aux_raw_select)
      print( paste0("raw measures selected: ", paste( aux_raw1_name,collapse="; ") ))
    })

    observe({
      req(rv_in$aux_comp)
      req(input$SI_aux_comp_select)
      if (input$SI_aux_comp_select[1] != "") { #|| ( !is.na(input$SI_var0) )
        shinyjs::show("ui_comp_factor1")
      } else {
        shinyjs::hide("ui_comp_factor1")
      }
    })
    output$ui_comp_factor1 <- renderText({
      req(rv_in$aux_comp)
      aux_comp1_name <- isolate(input$SI_aux_comp_select)
      print( paste0("comp measures selected: ", paste( aux_comp1_name,collapse="; ") ))
    })


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

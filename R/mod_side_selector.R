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
          col_8(selectizeInput(ns("SI_exp_fact_select"), "choose experimental factor: ", "",
                      options = list(placeholder = "load database first")))
            ), #fluidRow 1b

    fluidRow(textOutput(ns("ui_aux_fact_names"))),

    fluidRow(
      col_8(selectizeInput(ns("SI_aux_fact_select"),  "choose aux factors: ", "",
                    options = list(placeholder = "load database first")))
            ),
    fluidRow(
      col_8(offset=2,
            shinyjs::hidden(textOutput(ns("ui_aux_factor1"))),
            shinyjs::disabled(selectizeInput(ns("SI_aux_factor1_select"),  "aux factor values", "",
                                        multiple=TRUE,
                                        options = list(placeholder = "load db/choose aux factor ")))
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
   #      var_names = NULL,       # eg.  var_names: Genes, Proteins, Lipids
   #      var_annotations = NULL, # colnames of variable annotations.  e.g. gene/protein-families, ontologies, lipid class
   #
   #      obs_names = NULL,       # name of observations e.g. cell ID
   #      obs_annotations = NULL,
   #
   #      X   = NULL,
   #      var = NULL,
   #      obs = NULL,
   #
   #      uns = list(),
   #      uns_keys = NULL,
   #      varm = list(),
   #      varm_keys = NULL,
   #      obsm = list(),
   #      obsm_keys = NULL,
   #
   #                              # these reflect the choices made in the ingestor
   #                              #    omics key feature i.e. genes, proteins, lipids
   #                              #  primary & auxilarry
   #
   #      omics_feature = NULL,
   #      aux_features = NULL,
   #      exp_factor = NULL,        # observables:  experiemntal factor... i.e. patient/control, old/young,
   #      aux_factors = NULL,        # e.g. #   plust aucuilarry observations. e.g. QC- batch, sex, etc  or inferred: cluster, label, cell-type
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

    observeEvent(rv_in$trigger, {
      req(rv_in$trigger>0)
      print("TRIGGERED -------->")

      exp_fact_name <- isolate(rv_in$exp_factor)
      # show the factors that have been loaded
      output$ui_exp_fact_name <- renderText({
        req(rv_in$exp_factor)
        print(paste("current experimental factor: ", exp_fact_name))

      })
      #print(isolate(rv_in$factors0))
      exp_factor_vals <- isolate(rv_in$obs[[rv_in$exp_factor]])
      updateSelectizeInput(session, "SI_exp_fact_select", paste0(exp_fact_name," (exp fact):"), # wrap in ns?
                           choices = exp_factor_vals,
                           selected = exp_factor_vals[1], server=TRUE)
      # reset trigger for next time (if we don't reset, will everything just be responsive?)
      rv_in$trigger  <-  0


      if (!is.null(rv_in$aux_factors)){
        aux_fact_names <-  isolate(rv_in$aux_factors)
        output$ui_aux_fact_names <- renderText({
          req(rv_in$aux_factors)
          print(paste("auxillary factors: ", paste0( aux_fact_names,collapse="; ")))
        })
        updateSelectizeInput(session, "SI_aux_fact_select",
                             choices = aux_fact_names,
                             selected = aux_fact_names[1], server=TRUE)
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

    }) #observe event

    observe({
      if (input$SI_aux_fact_select != "") { #|| ( !is.na(input$SI_var0) )
        shinyjs::enable("SI_aux_factor1_select")

        aux_facts1 <- isolate(rv_in$obs[[input$SI_aux_fact_select]])
        aux_factor1_name <- isolate(input$SI_aux_fact_select)

        updateSelectizeInput(session, "SI_aux_factor1_select", paste0(aux_factor1_name," (aux factor)"),
                             choices = aux_facts1,
                             selected = aux_facts1[1], server=TRUE)

        shinyjs::show("ui_aux_factor1")
      } else {
        print("skipped SI_aux_factor1_select")
        shinyjs::disable("SI_aux_factor1_select")
        updateSelectizeInput(session, "SI_aux_factor1_select",
                             "aux factor values", "",
                             options = list(placeholder = "load db/choose aux factor "),
                             server=TRUE)
        shinyjs::hide("ui_aux_factor1")

      }
    })

    output$ui_aux_factor1 <- renderText({
      req(rv_in$aux_factors)
      aux_factor1_name <- isolate(input$SI_aux_fact_select)
      print( paste0(" auxillary factors selected: ", paste( aux_factor1_name,collapse="; ") ))
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
        omics_choice_list <- rv_in$var[[rv_in$omics_feature]]

        #if ( !is.null(omics_choice_list$var) & !is.null(rv_in$omics_feature) ) {
        if ( is.null(omics_choice_list)  ) {
          omics_choice_list <- "" #rownames(rv_in$var)
        }

        omics_choices <- isolate(omics_list$value)
        updateSelectizeInput(session, "SI_omics_select",
                             choices = isolate(omics_choice_list),
                             selected = omics_choices, server=TRUE)
        # DEBUG
        print("   ::updated omics list::  ")
        print(isolate(omics_list$value))
      }
      print(paste0("current number of selected omics:  ",length(unique(isolate(omics_list$value)))))

    })



    ############################+
    ## "reset" and "submit" simply sets the viz_now flag
    ############################+
    ############################
    observe({
      # turn on if the "placeholder" is gone
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
      } else {

      } # else
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

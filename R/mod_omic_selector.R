#' omic_selector UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_omic_selector_ui <- function(id){
  ns <- NS(id)
  omic_sel_tags <- tagList(

    fluidRow(
      selectizeInput(
        ns("SI_omics_select"), "Select omic features", "",
        multiple=TRUE, options = list(placeholder = "Choose omic feature (i.e. genes,proteins,lipids...)")
      )
    ), #fluidRow 2

    fluidRow(
      column(width = 4 ,offset = 0,
             style='padding-left:0px; padding-right:1px',
             actionButton(ns("AB_omics_submit"),"Submit",class="hidableSubmit")
      ),
      column(width = 4, offset = 1,
             style='padding-left:0px; padding-right:1px',
             actionButton(ns("AB_omics_reset"), "Clear",class="hidableClear")
      )
    ),
    fluidRow(
      uiOutput(ns("ui_text_warn"), width = "100%"),
    ) #fluidRow 3

  )

  return(omic_sel_tags)
}

#' omic_selector Server Functions
#'
#' @noRd
mod_omic_selector_server <- function(id, omics_in){
  moduleServer( id, function(input, output, session){
    ns <- session$ns


    # keep track of which proteins have been selected
    omics_list <- reactiveValues(value=character(0), viz_now=FALSE)

    # Error handling: check if too many omics selected
    # NOT within_limit  --> print error statement
    # within_limit --> Update the omics selection
    # THIS GETS RUN ON EVERY CLICK!!
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
        #assert that the choices are updated...
        l_omics_in <- omics_in()
        omics_choice_list <- isolate(l_omics_in)
        if ( is.null(omics_choice_list)  ) {
          omics_choice_list <- "" #rownames(rv_in$var)
        } else {
          omics_choice_list <- names(omics_choice_list)
        }

        omics_choices <- isolate(omics_list$value)

        omics_choices <- omics_choices[omics_choices %in% omics_choice_list] # will this fix stale list?

        freezeReactiveValue(input, "SI_omics_select")
        updateSelectizeInput(session, "SI_omics_select",
                             choices = omics_choice_list,
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

    ############################
    observe({ # turn on if the "placeholder" is gone
      shinyjs::toggleState("AB_omics_submit", !all(input$SI_omics_select == "Choose omic feature (i.e. genes,proteins,lipids...)"))
    })

    # TODO:  force the list to reset when the database is re-loaded...
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
      } # else ?? print warning??
    })


  return(omics_list)

  })

}

## To be copied in the UI
# mod_omic_selector_ui("omic_selector_ui_1")

## To be copied in the server
# mod_omic_selector_server("omic_selector_ui_1", rv_in)

# start_omics <- c("SNAP25,ENO2,SLC17A7,SLC17A6,GAD1,GAD2,SLC32A1,LAMP5,SST,CHODL,PVALB,VIP,CUX2,RORB,RBP4,GJA1,FGFR3,GFAP,OLIG1,OPALIN,PDGFRA,AIF1,TYROBP,NOSTRIN")


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

      column(width = 2, offset = 0,
             #style='padding-left:0px; padding-right:1px',
             actionButton(ns("AB_omics_reset"), "Clear", class = "btn btn-primary" ) #class="hidableClear"
      ),
      column(width = 2, offset = 0,
             #style='padding-left:0px; padding-right:1px',
             actionButton(ns("AB_omics_def"), "Default", class = "btn btn-primary" ) #class="hidableDefault"
      ),
      column(width = 2, offset = 0,style="border-right: 2px solid black"),
      column(width = 2 ,offset = 4,
             style='padding-left:0px; padding-right:1px',
             actionButton(ns("AB_omics_submit"),"Submit",class = "btn btn-primary" ) #class="hidableSubmit"
      )
    ),
    fluidRow(
      uiOutput(ns("ui_text_warn"), width = "100%"),
    ), #fluidRow 3
    # COPY/PASTE HACK JAVASCRIPT BELOW
    tags$script(
      HTML(
        ' console.log("page will load now");
          document.addEventListener("DOMContentLoaded", function(){
            console.log("page loaded");
            document.addEventListener("copy", (event) => {
              console.log("coppying from item:", event.target);
              const anchorNode = document.getSelection().anchorNode
              if (anchorNode instanceof HTMLElement && anchorNode.classList.contains("selectize-input")) {
                const items = Array.from(anchorNode.getElementsByClassName("item active"))
                const selectedItemsAsString = items.map(i => i.innerText).join(", ")
                console.log("coppied content:", selectedItemsAsString);
                event.clipboardData.setData("text/plain", selectedItemsAsString)
                event.preventDefault()
              }
            })
          });'
      )
    )

  )

  return(omic_sel_tags)
}

#' omic_selector Server Functions
#'
#' @param id id for module
#' @param all_omics complete list
#' @param def_omics default (target omics)
#' @param new_db_trig signal that we have a new database loaded
#'
#' @noRd
mod_omic_selector_server <- function(id, all_omics, def_omics, new_db_trig) {
  moduleServer( id, function(input, output, session) {
    ns <- session$ns



    observeEvent(
      new_db_trig(),  # new database loaded
      {
        req(all_omics,
            def_omics)  # set when database is chosen
        if (new_db_trig()>0) {
          omics_list$viz_now = FALSE
          omics_list$value <- isolate(def_omics())
          # get rid of error message when resetting selection
          output$ui_text_warn <- renderUI({ })
        }
      },
      ignoreNULL = TRUE,
      ignoreInit = TRUE
    ) #observe



    # keep track of which proteins have been selected start with default for
    # TODO:  check if we need to "isolate"
    omics_list <- reactiveValues(value=isolate(def_omics()), viz_now=FALSE)

    max_omic_feats <- 100


    observe({
      #toggle the vis_now flag if we are changing features
      #omics_list$viz_now = FALSE

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
        l_all_omics <- all_omics()
        omics_choice_list <- isolate(l_all_omics)
        if ( is.null(omics_choice_list)  ) {
          omics_choice_list <- "" #rownames(rv_in$var)
        } else {
          if (!is.null(names(omics_choice_list))) {
          omics_choice_list <- names(omics_choice_list)
        }}

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

    observeEvent(input$AB_omics_reset, {
      omics_list$viz_now <- FALSE
      omics_list$value <- character(0)
    })

    observeEvent(input$AB_omics_def, {
      omics_list$viz_now <- FALSE
      omics_list$value <- isolate(def_omics())
    })

    observeEvent(input$AB_omics_submit, {
      if(length(unique(omics_list$value)) < max_omic_feats ) { #defensive
        omics_list$value <- input$SI_omics_select #include direct selection from protein-box when pressing submit
        omics_list$viz_now <- TRUE

        output$ui_text_warn <- renderUI({
         p("updated plot")
        })
      }
    })

    observe({
      req(omics_list$viz_now)
      if(omics_list$viz_now == FALSE) { #defensive
        output$ui_text_warn <- renderUI({ })
      } # else ?? print warning??
    })

    # # dOES THIS PATTERN HELP US UPDATE PROPERLY???
    # reset_omics <- reactiveValues(reset = NULL)
    # observeEvent(input$AB_omics_reset, {
    #   reset_omics$reset <- TRUE
    # })
    # observeEvent(all_omics, {
    #   reset_omics$reset <- TRUE
    # })
    # #    observeEvent(reset_omics$reset, {
    # omics_list$viz_now = FALSE
    # omics_list$value <- character(0)
    # # get rid of error message when resetting selection
    # output$ui_text_warn <- renderUI({ })


  return(omics_list)

  })

}

## To be copied in the UI
# mod_omic_selector_ui("omic_selector_ui_1")

## To be copied in the server
# mod_omic_selector_server("omic_selector_ui_1", rv_in)

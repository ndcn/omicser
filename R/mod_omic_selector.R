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


    h4("Target *-Omics"),
    fluidRow(
         selectizeInput(
                        ns("SI_omics_select"), "Select (<100)", "",
                        multiple=TRUE, options = list(placeholder = "Choose omic feature (i.e. genes,proteins,lipids...)")
                      )

    ), #fluidRow 2

    fluidRow(
      column(width = 2, offset = 0,
             #style='padding-left:0px; padding-right:1px',
             actionButton(ns("AB_omics_reset"), "Clear", class = "btn btn-primary" ) #class="hidableClear"
      ),
      column(width = 2, offset = 0,style="border-right: 2px solid black",
             #style='padding-left:0px; padding-right:1px',
             actionButton(ns("AB_omics_def"), "Default", class = "btn btn-primary" ) #class="hidableDefault"
      ),
      column(width = 6, offset = 0,
              textInput(
                inputId = ns("text_omic_add"),
                label = "omics to add", width = "100%",
                value = ""
                #value = "Gene1,Gene2"
                )
      ),
      column(
        2,
        style='padding-left:0px; padding-right:1px',
        actionButton(ns("AB_text_omic_add"),"add ⤴",class = "btn btn-primary" )
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
#' @param active_omics complete list
#' @param def_omics default (target omics)
#' @param new_db_trig signal that we have a new database loaded
#'
#' @noRd
mod_omic_selector_server <- function(id, active_omics, def_omics){ #, new_db_trig) {
  moduleServer( id, function(input, output, session) {
    ns <- session$ns



    # Define return rvs: `selected_omics`  -----------------
    selected_omics <- reactiveValues(
                              target_omics=NULL, #isolate(def_omics()),
                              all_omics= NULL, #isolate(active_omics()),
                              freeze = NULL
    )


    # OBSERVES -----------------
    observe({
      #req(active_omics())  # set when database is chosen ... this is essentially a reset...
      selected_omics$all_omics <- active_omics()
      selected_omics$freeze <- 0 #reset ffreeze?
      output$ui_text_warn <- renderUI({ })
    })


    observe({
      #req(def_omics() )  # set when database is chosen
      selected_omics$target_omics <- def_omics()
    })

    # # REACTIVES    -----------------
    # plus_omics <- eventReactive( input$AB_text_omic_add, {
    #   print("plus_omics()")
    #   return(  strsplit(input$text_omic_add, ",| |;|, ")[[1]] )
    # })

    # REACTIVES    -----------------
    plus_omics <- reactive({
      input$AB_text_omic_add  #reactive on the button rather than any text input
      return( strsplit(isolate(input$text_omic_add), ",| |;|, ")[[1]] )

    })



    # observeEvent(
    #   new_db_trig(),  # new database loaded
    #   {
    #     req(active_omics(),
    #         def_omics)  # set when database is chosen
    #
    #     if (new_db_trig()>0) {
    #       selected_omics$target_omics <- def_omics()
    #       selected_omics$all_omics <- active_omics()
    #       selected_omics$freeze <- 0
    #
    #       # get rid of error message when resetting selection
    #       output$ui_text_warn <- renderUI({ })
    #     }
    #   },
    #   ignoreNULL = TRUE,
    #   ignoreInit = TRUE
    # ) #observe


    # # text input add to list...
    # omic_add_input <- reactive({
    #   strsplit(input$text_omic_add, ",| |;|, ")[[1]]
    # })


    max_omic_feats <- 100

    observe({
      #toggle the vis_now flag if we are changing features
      req(selected_omics$all_omics)# don't trip

      omics_choices <- selected_omics$target_omics #don't update if the target_omics change (avoid loop)

      if ( length( unique( omics_choices ) ) >= max_omic_feats ) {
        print("TOOO MANY FEATURES...")
        print(paste0("n omics choices. 2A.. ", length(omics_choices)) )


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

        omics_choice_list <- active_omics()

        # if ( is.null(omics_choice_list)  ) {
        #   omics_choice_list <- "" #rownames(rv_in$var)
        # } else {
        #   if (!is.null(names(omics_choice_list))) {
        #   omics_choice_list <- names(omics_choice_list)
        # }}

        #omics_choices <- isolate(selected_omics$target_omics) #don't update if the target_omics change (avoid loop)
        omics_choices <- unique( append( omics_choices,plus_omics() ) )
        omics_choices <- omics_choices[omics_choices %in% omics_choice_list] # will this fix stale list?

        freezeReactiveValue(input, "SI_omics_select")
        updateSelectizeInput(session, "SI_omics_select",
                             choices = omics_choice_list,
                             selected = omics_choices, server=TRUE)

      }

    })



    ############################+
    ## "reset" and "submit" simply sets the viz_now flag
    ############################+
    ############################

    observeEvent(input$AB_omics_reset, {
      selected_omics$target_omics <- character(0)
      #selected_omics$all_omics <- active_omics()
      print(length(selected_omics$target_omics))

      })

    observeEvent(input$AB_omics_def, {
      selected_omics$target_omics <- def_omics()
      print(length(selected_omics$target_omics))
      #selected_omics$all_omics <- active_omics() #unnesscessary...
    })

    observeEvent(
      selected_omics$freeze, {  # new database loaded
      # selected_omics$target_omics <- isolate(input$SI_omics_select) #include direct selection from protein-box when pressing submit
      # selected_omics$all_omics <- isolate(active_omics())
      selected_omics$target_omics <- input$SI_omics_select #include direct selection from protein-box when pressing submit

      })


  return(selected_omics)

  })

}

## To be copied in the UI
# mod_omic_selector_ui("omic_selector_ui_1")

## To be copied in the server
# mod_omic_selector_server("omic_selector_ui_1", rv_in)

#' subset_selector UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_subset_selector_ui <- function(id){
  ns <- NS(id)

  tags <- uiOutput(ns("ui_subset"))

  return(tags)
}

#' subset_selector Server Functions
#'
#' @noRd
mod_subset_selector_server <- function(id,rv_conf,meta_field){
  moduleServer( id, function(input, output, session){
    ns <- session$ns


    ### Reactive expressions ============================================
    ret_rv_sub <- reactiveValues(
      set = NULL,
      select = NULL
    )


    ## dynamic subset UI group
    output$ui_subset <- renderUI({
      cfg <- rv_conf()

      req(cfg)

      #subs_label <- paste0("select ",isolate(input$SI_subset),"s: ")
      choices <- cfg[grp == TRUE & field==meta_field]$UI

      ret_tags <- tagList(
        fluidRow(
          column(
            width = 4,
            offset = 0,
            selectizeInput(ns("SI_subset"),
                           "subset by:",
                           choices =  choices,
                           selected = choices[1])
          ),
          column(width = 5,
                 offset = 0,
                 uiOutput(ns("ui_subset_sel")),
          ),
          column(width = 1,
                 offset = 0,
                 actionButton(ns("CB_sub_none"), "NONE", class = "btn btn-primary" ),
                 br(),
                 actionButton(ns("CB_sub_all"), "ALL", class = "btn btn-primary")
          )

        ) #fluidrow
      ) #taglist

      return(ret_tags)
    })


    ## dynamic subset selector checkboxes
    output$ui_subset_sel <- renderUI({
      cfg <- rv_conf()

      req(input$SI_subset,
          cfg)

      # check to see if we have the input$SI_subset in our config
      # if not return andwait for the SI_subset to be updated
      if ( !(input$SI_subset %in% cfg$UI) ) {
        print("input$SI_subset no yet updated")
        return(NULL)
      }

      subs <- strsplit(cfg[UI == input$SI_subset]$fID, "\\|")[[1]]
      # sorting hack
      subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
      if (all(!is.na(subs2))){
        names(subs2) <- seq_along(subs2)
        subs <- subs[as.numeric(names(sort(subs2)))]
      }
      subs_label <- paste0("select ",isolate(input$SI_subset),"s: ")

      ret_tags <- tagList(
        checkboxGroupInput( ns("CB_subsel"),
                            label = subs_label,
                            inline = TRUE,
                            choices = subs,
                            selected = subs)
      ) #taglist

      return(ret_tags)

    })


    ### observe s =========================================================
    observeEvent(input$CB_sub_all, {
      cfg <- rv_conf()
      req(cfg)

      subs <- strsplit(cfg[UI == input$SI_subset]$fID, "\\|")[[1]]
      subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
      if (all(!is.na(subs2))){
        names(subs2) <- seq_along(subs2)
        subs <- subs[as.numeric(names(sort(subs2)))]
      }
      subs_label <- paste0("select ",isolate(input$SI_subset),"s: ")

      freezeReactiveValue(input, "CB_subsel")
      updateCheckboxGroupInput(session,
                               inputId = "CB_subsel",
                               label = subs_label,
                               choices = subs,
                               selected = subs,
                               inline = TRUE) # WARNING.  make sure subset is null is checked (length(0?))
    })


    observeEvent(input$CB_sub_none, {
      cfg <- rv_conf()
      req(cfg)

      subs <- strsplit(cfg[UI == input$SI_subset]$fID, "\\|")[[1]]
      subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
      if (all(!is.na(subs2))){
        names(subs2) <- seq_along(subs2)
        subs <- subs[as.numeric(names(sort(subs2)))]
      }

      subs_label <- paste0("select ",isolate(input$SI_subset),"s: ")

      freezeReactiveValue(input, "CB_subsel")
      updateCheckboxGroupInput(session,
                               inputId = "CB_subsel",
                               label = subs_label,
                               choices = subs,
                               selected = character(0), # character(0), #, "" , NULL
                               inline = TRUE) # WARNING.  make sure subset is null is checked (length(0?))
    })


    observe({
            ret_rv_sub$set <- input$SI_subset
            ret_rv_sub$select <- input$CB_subsel
            })

    ### RETURN =========================================================
    return(ret_rv_sub)

  })
}

## To be copied in the UI
# mod_subset_selector_ui("subset_selector_ui_1")

## To be copied in the server
# mod_subset_selector_server("subset_selector_ui_1")

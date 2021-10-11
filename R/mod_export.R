#' export UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_export_ui <- function(id){
  ns <- NS(id)
  tagList(

    # textOutput(ns("var0")),
    # textOutput(ns("var1")),
    # textOutput(ns("target_omics")),

    verbatimTextOutput(ns("var0"),placeholder = TRUE),
    verbatimTextOutput(ns("var1"),placeholder = TRUE),
    verbatimTextOutput(ns("target_omics"),placeholder = TRUE),

    verbatimTextOutput(ns("verb"),placeholder = TRUE),
    verbatimTextOutput(ns("text1"),placeholder = TRUE),
    verbatimTextOutput(ns("names"),placeholder = TRUE),


    verbatimTextOutput(ns("rv_1"),placeholder = TRUE),
    verbatimTextOutput(ns("rv_2"),placeholder = TRUE)

    ## some dummies to plot our recative values to test passing...

  )
}

#' export Server Functions
#'
#' @noRd
mod_export_server <- function(id, rv_in, p){
  moduleServer( id, function(input, output, session){
    ns <- session$ns


    #mod_vis_server("vis_ui_1",db=db,p=p)

    # exp_fact = NULL,
    # aux_fact = NULL,
    # omics_names = NULL,
    # target_omics = NULL,
    # plot_type = NULL,
    # comp_dat = NULL,
    # obs_dat = NULL


    #output$text1     <- renderText({ vecFun()})
    output$var0      <- renderText({ paste0(isolate(p[["comp_dat"]]),collapse = ", ")  })

    output$var1      <- renderText({ paste0(isolate(p[["obs_dat"]]),collapse = ", ")   })

    output$target_omics <- renderText({ paste0(isolate(p$target_omics$value),collapse = ", ")  })


    output$verb       <- renderText({ p[["plot_type"]] })



    output$text1<- renderPrint({ list((p$comp_dat),
                                      (p$obs_dat),
                                      (p$omics_names),
                                      (p$target_omics$value),
                                      (p$target_omics$viz_now),
                                      (p$plot_type)) })

    output$names0 <- renderText({ p[["omics_names"]] })

    output$rv_1 <- renderText({

      omics_choice_list <- isolate(rv_in$omics_feature)
      #if ( !is.null(omics_choice_list$var) & !is.null(rv_in$omics_feature) ) {
      if ( is.null(omics_choice_list)  ) {
        omics_choice_list <- "" #rownames(rv_in$var)
      } else {
        omics_choice_list <- names(omics_choice_list)[1:25]
        }
      paste0(omics_choice_list,collapse = ", ")
      })

    output$rv_2 <-  renderText({
      paste0(isolate( rv_in$ad$var[[ rv_in$ad$var_keys()[1] ]] )[1:20],collapse = ", ")
      })

  })
}

## To be copied in the UI
# mod_export_ui("export_ui_1")

## To be copied in the server
# mod_export_server("export_ui_1")

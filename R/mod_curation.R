#' curation UI Function
#'
#' @description A shiny Module for the curation of new data sets. This module
#'     doesn't contain any real functionality for curation, but is only the
#'     starting point for curation.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @author Rico Derks
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#'
mod_curation_ui <- function(id){
  ns <- NS(id)
  tagList(
    actionButton(ns("AB_curate_data"), label = "Curate data", class = "btn btn-large btn-danger")
  )
}

#' curation Server Functions
#'
#' @noRd
#'
mod_curation_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    # which omic is selected
    omic_selected <- reactiveVal()

    # show the modal
    observeEvent(input$AB_curate_data, {
      # show a popup for the curation part
      showModal(
        modalDialog(
          title = "Curation",
          easyClose = FALSE,
          fade = TRUE,
          size = "l",
          footer = tagList(
            modalButton("Cancel"),
          ),
          radioButtons(
            inputId = ns("rb_select_omic"),
            label = "Select omics for curation :",
            choices = c("Lipidomics" = "lipid",
                        "Proteomics" = "prot",
                        "Transcriptomics" = "trans"),
            selected = "lipid"
          ),
          # show the curation stuff for each omics
          uiOutput(outputId = ns("ui_omic_params"))
        )
      )
    })

    # get which omic parameter is selected
    observeEvent(input$rb_select_omic, {
      req(input$rb_select_omic)

      omic_selected(input$rb_select_omic)
    })


    # for which omics to show parameters
    output$ui_omic_params <- renderUI({
      req(omic_selected())

      # get the selected omic
      selected_omic <- omic_selected()

      tagList(
        hr(),
        # show the module for the selected omic
        switch(
          selected_omic,
          "lipid" = mod_curation_lipid_ui(id = ns("curation_lipid_ui"))
        )
      )
    })

  })
}

#' additional_info UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_additional_info_ui <- function(id){
  ns <- NS(id)
  tagList(
    wellPanel(
      id = ns("about"),
      htmlOutput(outputId = ns("additional_info_md")),
      uiOutput(outputId = ns("additional_info_md2"))
    )
  )
}

#' additional_info Server Functions
#'
#' @param id shiny internal
#' @param rv_data main reactive value (just getting the database directory)
#' @param DB_ROOT_PATH where do our databases live
#'
#' @noRd
mod_additional_info_server <- function(id,rv_data, DB_ROOT_PATH){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    # show the factors that have been loaded
    output$additional_info_md <- renderUI({
      req(rv_data$db_meta$name)

      md_path <- file.path(DB_ROOT_PATH,rv_data$db_meta$name,"additional_info.Rmd")

      out_htm <- HTML(markdown::markdownToHTML(md_path))

      return(out_htm)
    })
    # show the factors that have been loaded
    output$additional_info_md2 <- renderUI({
      req(rv_data$db_meta$name)

      md_path <- file.path(DB_ROOT_PATH,rv_data$db_meta$name,"additional_info.Rmd")

      out_htm <- HTML(markdown::markdownToHTML(md_path))

      return(out_htm)
    })

  })
}

## To be copied in the UI
# mod_additional_info_ui("additional_info_ui_1")

## To be copied in the server
# mod_additional_info_server("additional_info_ui_1")

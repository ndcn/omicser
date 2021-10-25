#' help UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_help_ui <- function(id){
  ns <- NS(id)
  tagList(
    h4("WORK IN PROGRESS"),

    HTML('

              <body>
              <hr>
              <div>

              MORE INFO: <br>
              <a href="https://ndcn.github.io/omicser" target="_blank">Documentaion Site</a>
<br>
               <a href="https://github.com/ndcn/omicser" target="_blank">ndcn/omicser@github</a>

              </div>
              </body>
              '
    )
 )
}

#' help Server Functions
#'
#' @noRd
mod_help_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

  })
}

## To be copied in the UI
# mod_help_ui("help_ui_1")

## To be copied in the server
# mod_help_server("help_ui_1")

#' Convert a Word (docx) to a Markdown (md) file
#'
#' This function converts a word document to a markdown file.
#' If output path is NULL (as default), the function will save this to a temporary path.
#'
#' @param file The path to your .docx file.
#' @param output The path to save your new markdown file
#'
#' @return The same file as a Markdown document.
#' @importFrom rmarkdown pandoc_convert
#' @examples
#'
#' #ADD LATER
#'
word_to_markdown <- function(file, output = NULL) {

  if (base::is.null(output)) {
    outputFile <- base::tempfile(fileext = ".md")
  } else {
    outputFile <- output
  }

  pandoc_convert(input = file, to = "markdown", output = outputFile)

  return(outputFile)

}

#' Display a document in Shiny.
#'
#' This function will take a Markdown document (.md) or Word document (.docx) and
#' return the UI for a Shiny app to display it. This should be embedded within a
#' UI code (see examples for details).
#'
#'@param file The file (path) containing a document to display.
#'
#'@return html tags of the rendered doc
#'
#' @examples
#'
#' \dontrun{
#' ui <- shiny::fluidPage(
#' shiny::mainPanel(
#' display_document('path-to-your-file')
#' )
#' )
#'}
display_document <- function(file) {

  if (base::grepl(".docx", file)) {
    #TODO: test .doc and .docx rendering works
    markdown_file <- word_to_markdown(file, output = NULL)
    base::on.exit(base::unlink(markdown_file))
  } else if (base::grepl(".md", file)) { # NOTE: grepl returns true for .Rmd and .md
    markdown_file <- file
  } else {
    stop("Invalid file type. Please add a Markdown .md/.Rmd or Word .docx file.")
  }

  shiny::includeMarkdown(markdown_file)

}







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
mod_additional_info_server <- function(id,db_name, DB_ROOT_PATH){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    # show the factors that have been loaded
    output$additional_info_md <- renderUI({
      req(db_name$name)

      md_path <- file.path(DB_ROOT_PATH,db_name$name,"additional_info.Rmd")
      #TODO:  allow a .docx file instead which will need to be rendered and then output
      #
      if (!file.exists(md_path)) { #LOAD DEFAULT MESSAGE
        md_path <- system.file("app/www/additional_info.Rmd",package='omicser')
      }

      # TODO: check for .docx OR .md
      out_htm <- display_document(md_path)

      #out_htm <- shiny::includeMarkdown(md_path)


      return(out_htm)
    })

  })
}

## To be copied in the UI
# mod_additional_info_ui("additional_info_ui_1")

## To be copied in the server
# mod_additional_info_server("additional_info_ui_1")

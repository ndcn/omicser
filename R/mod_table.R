#' pg_table UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @importFrom DT DTOutput
mod_table_ui <- function(id){
  ns <- NS(id)
  tagList(
    column(width =12,
           DTOutput(ns("table_ui_1"),
                    width = "100%",
                    height = "800px"
           )
    )
  )
}

#' pg_table Server Functions
#'
#' @noRd
#' @importFrom DT datatable renderDT
mod_table_server <- function(id,dt_in) {
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    output$table_ui_1 <- renderDT({
      # in case we didn't clean up factors during curation ...
      #dat <- dt_in() %>% dplyr::mutate(dplyr::across(where(is.factor), as.character ))
      dat <- dt_in()
      dat[] <- lapply(dat, function(x) if (is.factor(x)) as.character(x) else {x})

      # only text fields might be too long....
      select_char_target <-  which(sapply(dat, is.character))  #filter these ones..
      names(select_char_target) <- NULL
      select_num_target <-  which(sapply(dat, is.numeric))  #filter these ones..
      names(select_num_target) <- NULL

      #dynamically display the shortened version
      render_char_js <- JS(
        "function(data, type, row, meta) {",
        "return type === 'display' && data.length > 20 ?",
        "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
        "}")


      render_num_js = JS(
        "function(data, type, row, meta){",
        "return type === 'display' ?",
        "data.toExponential(2) : data;",
        "}")


      # hover "tips" to get the full text field
      callback_js <- JS(
        "table.on('mouseover', 'td', function(){",
        "  var index = table.cell(this).index();",
        "  Shiny.setInputValue('cell', index, {priority: 'event'});",
        "});"
      )

      datatable(dat,
              class   = 'row-border order-column stripe hover display',
              rownames = TRUE, # need to add +1 to column indexes
              extensions = "Buttons",
              selection = 'single',
              fillContainer = TRUE,
              options = list(
                scrollX = TRUE,
                autoWidth  = FALSE,
                columnDefs = list(
                    list(
                      targets = select_char_target,
                      render = render_char_js),
                    list(
                      targets = select_num_target,
                      render = render_num_js)
                  ),
                pageLength = 25,
                lengthMenu = c(25, 50, 100,200,500),
                dom = 'Bfrtip',
                buttons = list('pageLength', 'csv', 'print'),
                processing = TRUE, # processing & serverSide must be FALSE in order for scientific to work proper
                serverSide = TRUE
              ),
              escape  = FALSE,
              callback = callback_js)
  },
  server = TRUE  # make this explicit
  )



  })
}

## To be copied in the UI
# mod_table_ui("table_ui_1")

## To be copied in the server
# mod_table_server("table_ui_1")

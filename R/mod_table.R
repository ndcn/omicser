#' table UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_table_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      shinydashboardPlus::box(
        uiOutput("table_header"),
        DT::dataTableOutput("table_over_info")
      )
    ),
    fluidRow(
      shinydashboardPlus::box(id = ns("omic_bp_bx"),
          title = "Metabolite Distribution",
          solidHeader = TRUE,
          status = "primary",
          plotlyOutput(ns("omic_bp")),
          width = 8),
      # Info on metab selected in table
      shinydashboardPlus::box(id = ns("omic_info_bx"),
          title = "Metabolite Details",
          solidHeader = TRUE,
          status = "primary",
          htmlOutput(ns("omic_info")),
          width = 4)
    )
    )

}

#' table Server Functions
#'
#' @param id
#' @param rv_data
#' @param prm
#'
#' @noRd
mod_table_server <- function(id, rv_data, rv_selections) {
  moduleServer( id, function(input, output, session){
    ns <- session$ns


    output$table_over_info <- DT::renderDataTable({

      dat <- rv_data$ad$obs
      # re-index based on our side-selector proteins... and highlight them
      # after reindexing dataset, highlight rows of already selected proteins
      proxy %>% selectRows(row.names(dat)[dat$gene.name %in% unique(c(isolate(input$protein_sel), isolate(input$protein_sel_table)))])
      target <- switch(input$data_sel, 'TMT' = c(1,7), 'DIA' = c(1,5))
      order_col <- switch(input$data_sel, 'TMT' = 7, 'DIA' = 5)

      DT::datatable(dat,
                    class   = 'row-border order-column stripe hover',
                    rownames = FALSE,
                    extensions = "Buttons",
                    selection = list(mode = 'multiple', target = 'row'),
                    options = list(
                      scrollX = TRUE,
                      autoWidth  = FALSE,
                      initComplete = JS(
                        "function(settings, json) {",
                        "$(this.api().table().header()).css({'background-color': '#008b42', opacity: .7, 'color': '#ffffff'});",
                        "}"),
                      columnDefs = list(list(
                        targets = target,
                        render = JS(
                          "function(data, type, row, meta) {",
                          "return type === 'display' && data.length > 10 ?",
                          "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
                          "}"))
                      ),
                      LengthMenu = c(10, 30, 50),
                      pageLength = 10,
                      ordering = TRUE,
                      order = list(c(1, 'asc'), c(order_col, 'desc')),
                      orderClasses = TRUE,
                      dom = 'Bfrtip',
                      buttons = list('pageLength', 'csv', 'print'),
                      processing = FALSE, # processing & serverSide must be FALSE in order for scientific to work proper
                      serverSide = FALSE),
                    escape  = FALSE)
    })

    output$table_header <- renderUI({
      datname <- switch(input$data_sel, 'TMT' = 'Skeletal Muscle', 'DIA' = 'Muscle Stem Cell')
      tagList(tags$h3(paste(datname, "Table"))
      )
    })

    values <- reactiveValues()

    observeEvent(input$table_rows_selected, {
      values$metab_by_subject <- metab_by_subject
    })

    # Prepare data for boxplots
    metab_by_subject <- reactive({
      validate(
        need(
          input$table_rows_selected,
          "Select a metabolite to view plots"
        )
      )
      selected <- input$table_rows_selected
      metab_counts <- data.frame(
        colnames(dt_in)[1:454],
        unlist(c(dt_in[dt_in$BIOCHEMICAL == metab_meta[selected,1],1:454]))
      )
      colnames(metab_counts) <- c("id", "count")
      #metab_counts$count <- log10(metab_counts$count)

      merge(subject_status, metab_counts, by = "id")
    })

    values <- reactiveValues()

    observeEvent(input$table_rows_selected, {
      values$metab_by_subject <- metab_by_subject
    })

    # Reading distribution plot ouput
    output$omic_bp <- renderPlotly({
      req(values$metab_by_subject)
      # Plot by CACO status
      if (box_var() == "Status") {
        plot_ly(
          data = metab_by_subject()[metab_by_subject()$CACO %in% c("CO", "AD", "TREM2", "ADAD"),],
          y = ~count,
          color = ~CACO,
          colors = "Dark2",
          type = "box",
          boxpoints = "all",
          pointpos = 0,
          text = paste("Transformed Reading: ",
                       format(metab_by_subject()$count[metab_by_subject()$CACO %in% c("CO", "AD", "TREM2", "ADAD")], digits = 3)),
          hoverinfo = list("median", "text")
        ) %>%
          layout(
            title = metab_meta[input$table_rows_selected,1],
            yaxis = list(title = "Transformed Reading",
                         hoverformat = ".2f"),
            showlegend = FALSE
          ) %>%
          config(displayModeBar = FALSE)
      }

      # Plot by Sex
      else if (box_var() == "Sex") {
        plot_ly(
          data = metab_by_subject(),
          y = ~count,
          color = ~SEX,
          colors = "Dark2",
          type = "box",
          boxpoints = "all",
          pointpos = 0,
          text = paste("Log10 Reading: ",
                       format(metab_by_subject()$count, digits = 3)),
          hoverinfo = list("median", "text")
        ) %>%
          layout(
            title = metab_meta[input$table_rows_selected,1],
            yaxis = list(title = "Log10 Reading"),
            showlegend = FALSE
          ) %>%
          config(displayModeBar = FALSE)
      }
    })


    # Information on metabolite selected on table
    output$omic_bp <- renderText({

      metab_name <- as.character(metab_meta[input$table_rows_selected,1])

      paste(
        "Metabolite:", metab_name, "</br>",
        "Super Pathway:", metab_meta$`Super Pathway`[metab_meta$Metabolite == metab_name], "</br>",
        "Sub Pathway:", metab_meta$`Sub Pathway`[metab_meta$Metabolite == metab_name],
        # Link to HMDB, if ID is available
        ifelse(
          !is.na(metab_meta$`HMDB ID`[metab_meta$Metabolite == metab_name]),
          paste(
            "</br>HMDB ID: ",
            metab_meta$`HMDB ID`[metab_meta$Metabolite == metab_name]
          ),
          ""
        ),
        # Link to PubChem, if ID is available
        ifelse(
          !is.na(metab_meta$`PubChem ID`[metab_meta$Metabolite == metab_name]),
          paste("</br>PubChem ID: ",
                metab_meta$`PubChem ID`[metab_meta$Metabolite == metab_name]
          ),
          ""
        )
      )
    })

  })
}

## To be copied in the UI
# mod_table_ui("pg_table_ui_1")

## To be copied in the server
# mod_table_server("pg_table_ui_1")

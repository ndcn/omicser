#' table UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#' @importFrom shinydashboardPlus box
#' @importFrom shiny NS tagList
#' @importFrom DT DTOutput
mod_pg_table_ui <- function(id){
  ns <- NS(id)
  tagList(
    uiOutput(ns("table_header")),

    fluidRow(
      # box(id = ns("table_over_bx"),
      #     title = "Omic Data Table",
      #     status = "primary",
      #     width = 12,
      #     height = "1800px",
      #     solidHeader = TRUE,
      #     collapsible = FALSE,
          # uiOutput(ns("table_header")),
      column(width =12,
             DTOutput(ns("table_over_info"),
                   width = "100%",
                   height = "800px"
                   )
        )
    ),
    fluidRow(
      box(id = ns("omic_bp_bx"),
          title = "Expression Distribution",
          status = "warning",
          solidHeader = TRUE,
          plotlyOutput(ns("omic_bp")),
          width = 8),
      # Info on metab selected in table
      box(id = ns("omic_info_bx"),
          title = "Feature Meta-info",
          status = "warning",
          solidHeader = TRUE,
          htmlOutput(ns("omic_info")),
          width = 4)
    )
  )

}

#' table Server Functions
#'
#' @param id shiny internal
#' @param rv_data reactive data
#' @param rv_selections side selector reactives
#' @importFrom DT JS datatable renderDT dataTableProxy selectRows
#' @noRd
mod_pg_table_server <- function(id, rv_data, rv_selections, active_layer_data) {
  moduleServer( id, function(input, output, session){
    ns <- session$ns


    #create an omics-annotation table that includes differential expressions
    omic_table <- reactiveValues(
      dt = NULL
    )

   #proxy = dataTableProxy('table_over_info')
    observe({
      req(rv_data$anndata)

      rv_selections$GO #forces update on "GO"

      # create the table to render
      omics <- rv_selections$selected_omics$all_omics
      omics_idx <- which(rv_data$anndata$var_names %in% omics)

      join_data <- rv_data$anndata$var %>%
                        dplyr::mutate(dplyr::across(where(is.factor), as.character ))
      wide_data <- rv_data$de %>%
                        tidyr::pivot_wider(id_cols = names,
                                            names_from = c(group,test_type),
                                            values_from = c(logfoldchanges, scores ,pvals,pvals_adj,versus),
                                            names_repair = "check_unique"
                                          )
      comp_data <- join_data %>%
                        dplyr::inner_join(wide_data,by = c("feature_name"="names") )
      #comp_data <- as.data.frame(comp_data1)
      rownames(comp_data) <- comp_data$feature_name

      # SET THE DATA
      omic_table$dt <- comp_data

      #reset hover if we have a new dataset
      table_sel_values$omic_expr_values = NULL
      table_sel_values$omic_name = NULL
    })


    # Prepare data for boxplots
    omic_expr_values <- reactive({
      validate(
        need(
          input$table_over_info_rows_selected,
          "Select an omic to view plots"
        )
      )
      selected <- input$table_over_info_rows_selected
      omic_name <- omic_table$dt$feature_name[selected]
      data_vec <- active_layer_data$data[,omic_name]
      grouping_var <- rv_selections$observ_group_by
      # have NOT subsetted samples
      omic_counts <- data.frame( rv_data$anndata$obs_names,
                                 data_vec,
                                 rv_data$anndata$obs[[grouping_var]] )
      colnames(omic_counts) <- c("id", "val","grp")
      return(omic_counts)
    })


    # force waiting for a selection row
    table_sel_values <- reactiveValues(
      omic_expr_values = NULL,
      omic_name =NULL
    )

    observeEvent(input$table_over_info_rows_selected, {
      table_sel_values$omic_expr_values <- omic_expr_values()
      table_sel_values$omic_name <- omic_table$dt$feature_name[input$table_over_info_rows_selected]
    })

    output$table_over_info <- renderDT({
      req(rv_data$de,
          omic_table$dt)
      dt <- omic_table$dt

      # only text fields might be too long....
      select_char_target <-  which(sapply(dt, is.character))  #filter these ones..
      names(select_char_target) <- NULL
      select_num_target <-  which(sapply(dt, is.numeric))  #filter these ones..
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

      init_js <- JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#008b42', opacity: .7, 'color': '#ffffff'});",
        "}")


      callback_js <- JS(
        "table.on('mouseover', 'td', function(){",
        "  var index = table.cell(this).index();",
        "  Shiny.setInputValue('cell', index, {priority: 'event'});",
        "});"
      )

      #  TODO: enable this for list of columns that have "pvals" or "logfoldchange"
      # row_callback_js = JS(
      #   "function(row, data) {",
      #   "for (i = 8; i < data.length; i++) {",
      #   "$('td:eq('+i+')', row).html(data[i].toExponential(2));",
      #   "}",
      #   "}")

      datatable(dt,
              class   = 'row-border order-column stripe hover display',
              rownames = TRUE, # need to add +1 to column indexes
              extensions = "Buttons",
              selection = list(mode = 'single', selected = NULL, target = 'row', selectable = NULL), #'single',
              fillContainer = TRUE,
              options = list(
                scrollX = TRUE,
                autoWidth  = FALSE,
                initComplete = init_js,
                columnDefs = list(
                  list(
                    targets = select_char_target,
                    render = render_char_js),
                  list(
                    targets = select_num_target,
                    render = render_num_js)
                ),
                pageLength = 50,
                scrollY = "800px",
                # paging = FALSE,
                lengthMenu = c(25, 50, 100,200,500),
                # ordering = TRUE,
                # order = list(c(order_col, 'desc')),
                # orderClasses = TRUE,
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

    output$table_header <- renderUI({
      req(rv_data$db_meta$omics_type)
      datname <- rv_data$db_meta$omics_type
      tagList(tags$h3(paste(datname, "Table"))
      )
    })




    # Reading distribution plot ouput
    output$omic_bp <- renderPlotly({
      req(table_sel_values$omic_expr_values)

      # TODO:  make this into a boxplot function
      bx_data <- table_sel_values$omic_expr_values
      plot_ly(
        data = bx_data,
        y = ~val,
        color = ~grp,
        type = "box",
        boxpoints = "all",
        pointpos = 0,
        text = paste("Vals: ",
                     format(table_sel_values$omic_expr_values$val, digits = 3)),
        hoverinfo = list("median", "text"),
        colors = "Dark2"
      ) %>%
        plotly::layout(
          title = table_sel_values$omic_name,
          yaxis = list(title = "value",
                       hoverformat = ".2f"),
          showlegend = FALSE,
          xaxis = list(title = rv_selections$observ_group_by,
                       hoverformat = ".2f")
        ) %>%
        plotly::config(displayModeBar = FALSE)
    })



    # Information on metabolite selected on table
    output$omic_info <- renderText({
      req(table_sel_values$omic_name)

      print("in renderText: omic_info")

      omic_name <- table_sel_values$omic_name # omic_table$dt$feature_name[input$table_over_info_rows_selected ]

      omic_details <- rv_data$shaddow_defs$feat_deets

      annots <- rv_data$anndata$var[omic_name,]

      deet <- paste("<b>", rv_data$db_meta$omics_type, " </b> features </br>" )

      for (detail in omic_details) {
        value <- annots[[detail]]
        if (is.numeric(value) ) {
          value <- format(value, digits=3,scientific=TRUE)
        }
        deet <- paste(deet, "<b>", detail, "</b>:", value, "</br> ")
      }

      deet

      })


  })
}

## To be copied in the UI
# mod_pg_table_ui("pg_pg_table_ui_1")

## To be copied in the server
# mod_pg_table_server("pg_pg_table_ui_1")

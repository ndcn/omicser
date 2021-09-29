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
          title = "Omic Distribution",
          status = "warning",
          solidHeader = TRUE,
          plotlyOutput(ns("omic_bp")),
          width = 8),
      # Info on metab selected in table
      box(id = ns("omic_info_bx"),
          title = "Omic Details",
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
mod_pg_table_server <- function(id, rv_data, rv_selections) {
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    #create an omics-annotation table that includes differential expressions
    my_table <- reactive({
      req(rv_data$de)

      # # not sure why, but processing this as data.table causes the reactive
      # # return to jam pu the DT::datatable
      # # TODO:  find bug and remove dplyr
      # join_data <- as.data.table(rv_data$ad$var)[]
      # fix_cols <- names(join_data)[sapply(join_data, is.factor)]
      # join_data <- join_data[,(fix_cols):= lapply(.SD, as.character), .SDcols = fix_cols]
      # wide_data <- dcast(as.data.table(rv_data$de),
      #                    names ~ group+test_type,
      #                     value.var = c("logfoldchanges", "scores" ,"pvals","pvals_adj","versus"))
      # #
      # comp_data <- join_data[wide_data,on=c("omics_name"="names"),nomatch=0]


      join_data <- rv_data$ad$var %>% dplyr::mutate(dplyr::across(where(is.factor), as.character ))

      wide_data <- rv_data$de %>% tidyr::pivot_wider(
                    id_cols = names,
                    names_from = c(group,test_type),
                    values_from = c(logfoldchanges, scores ,pvals,pvals_adj,versus),
                    names_repair = "check_unique"
                  )
      comp_data <- join_data %>%
        dplyr::inner_join(wide_data,by = c("omics_name"="names") )

      #comp_data <- as.data.frame(comp_data1)
      rownames(comp_data) <- comp_data$omics_name

      return(comp_data)
    })

   #proxy = dataTableProxy('table_over_info')


    output$table_over_info <- renderDT({
      req(rv_data$omics)

      dat <- my_table()

      # this is a stub for collecting for teh side selector...
      #data.table(dat)
      # # re-index based on our side-selector proteins... and highlight them
      # # after reindexing dataset, highlight rows of already selected proteins
      # proxy %>% selectRows(row.names(dat)[dat$omics_name %in% isolate(rv_data$omics)])
      #order_col <- c(which(colnames(dat)=="omics_name"), which(colnames(data)==rv_selections$plot_var_x))

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

      datatable(dat,
              class   = 'row-border order-column stripe hover display',
              rownames = TRUE, # need to add +1 to column indexes
              extensions = "Buttons",
              selection = 'single',
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
      req(rv_data$omics_type)
      datname <- rv_data$omics_type
      tagList(tags$h3(paste(datname, "Table"))
      )
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

      omic_name <- my_table()$omics_name[selected]
      # we have the active layer in a reative in the pg: e.t.
      #data_vec <- active_layer_data$data[,omic_name]

      data_vec <- rv_data$ad$X[,omic_name]
      grouping_var <- rv_selections$plot_x
      omic_counts <- data.frame( rv_data$ad$obs_names,
                                   data_vec,
                                   rv_data$ad$obs[[grouping_var]] )
      colnames(omic_counts) <- c("id", "val","grp")
      return(omic_counts)

    })


    # force waiting for a selection row
    table_values <- reactiveValues()
    observeEvent(input$table_over_info_rows_selected, {
      table_values$omic_expr_values <- omic_expr_values()
      table_values$omic_name <- my_table()$omics_name[input$table_over_info_rows_selected]

    })

    # Reading distribution plot ouput
    output$omic_bp <- renderPlotly({
      req(table_values$omic_expr_values)
      # Plot by CACO status
      # TODO:  make this into a boxplot function
      #omic_name <- my_table()$omics_name[input$table_over_info_rows_selected]
      bx_data <- table_values$omic_expr_values
      plot_ly(
        data = bx_data, #table_values$omic_expr_values, #omic_expr_values(),
        y = ~val,
        color = ~grp,
        type = "box",
        boxpoints = "all",
        pointpos = 0,
        text = paste("Vals: ",
                     format(table_values$omic_expr_values$val, digits = 3)),
        hoverinfo = list("median", "text"),
        colors = "Dark2"
      ) %>%
        plotly::layout(
          title = table_values$omic_name,#my_table()$omics_name[input$table_over_info_rows_selected],
          yaxis = list(title = "value",
                       hoverformat = ".2f"),
          showlegend = FALSE
        ) %>%
        plotly::config(displayModeBar = FALSE)
    })



    # Information on metabolite selected on table
    output$omic_info <- renderText({
      req(table_values$omic_name)
      omic_name <- table_values$omic_name # my_table()$omics_name[input$table_over_info_rows_selected ]
      omic_details <- rv_data$default$omic_details
      annots <- rv_data$ad$var[omic_name,]

      deet <- paste("<b> Name </b>: ", omic_name, "</br>")

      for (detail in omic_details) {
        value <- annots[[detail]]
        if (is.numeric(value) ) {
          value <- format(value, digits=3,scientific=TRUE)
        }
        deet <- paste(deet, "<b>", detail, "</b>:", value, "</br> ")
      }

      deet

      })

        # "Super Pathway:", metab_meta$`Super Pathway`[metab_meta$Metabolite == metab_name], "</br>",
        # "Sub Pathway:", metab_meta$`Sub Pathway`[metab_meta$Metabolite == metab_name],
        # # Link to HMDB, if ID is available
        # ifelse(
        #   !is.na(metab_meta$`HMDB ID`[metab_meta$Metabolite == metab_name]),
        #   paste(
        #     "</br>HMDB ID: ",
        #     metab_meta$`HMDB ID`[metab_meta$Metabolite == metab_name]
        #   ),
        #   ""
        # ),
        # # Link to PubChem, if ID is available
        # ifelse(
        #   !is.na(metab_meta$`PubChem ID`[metab_meta$Metabolite == metab_name]),
        #   paste("</br>PubChem ID: ",
        #         metab_meta$`PubChem ID`[metab_meta$Metabolite == metab_name]
        #   ),
        #   ""
        #)

  })
}

## To be copied in the UI
# mod_pg_table_ui("pg_pg_table_ui_1")

## To be copied in the server
# mod_pg_table_server("pg_pg_table_ui_1")

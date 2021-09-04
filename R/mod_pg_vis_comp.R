
#' pg_vis_comp UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @importFrom plotly plotlyOutput
#' @importFrom shinydashboardPlus box
mod_pg_vis_comp_ui <- function(id){
  ns <- NS(id)
  tagList(


    HTML("Volcano plot"),
    h4("Relative expression of omics between two conditions"),
    "In this tab, users can visualise comparative measures of omic-expression in our data",
    "(e.g. logFC vs p-value).",
    br(),br(),
    fluidRow(
      column(
        3, style="border-right: 2px solid black",
        br(),

        # this dynamicall makes all of our test/group selectors
          uiOutput(outputId = ns("UI_comp_group_selection"))

      ), # End of column (6 space)
      column(9,
             #uiOutput(ns("UI_viz_output")),
               plotlyOutput(outputId = ns("volcano_plot"),height = "800px")
             )
      ),
      fluidRow(
        box(id = ns("omic_boxplot"),
            title = "Omic Distribution",
            solidHeader = TRUE,
            status = "primary",
            plotlyOutput(outputId = ns("test_boxplot")),
            width = 8),
        # Info on metab selected on volcano plot
        box(id = ns("omic_info"),
            title = "Omic Details",
            solidHeader = TRUE,
            status = "primary",
            htmlOutput(ns("UI_meta_box")),
            #uiOutput(outputId = ns("UI_meta_box")), #htmlOutput(ns("metabinfo")),
            width = 4)
      )

  ) #taglist

}
#' pg_vis_comp Server Functions
#'
#' @noRd
#' @importFrom plotly renderPlotly plotlyOutput plot_ly add_markers event_data
mod_pg_vis_comp_server <- function(id,rv_in, p){
  moduleServer( id, function(input, output, session){
    ns <- session$ns


    data_point_sz = .65
    plot_size = "Small"
    font_size = "Small"
    color_schemes = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple")
    color_scheme = color_schemes[2]


    #### compare samples
    #### TODO: change name to diff_exp_filter or something since we are no longer testing here
    ####    change to use dtdplyr for speed... and remove magrittr dependency
    filtered_de <- reactive({
      req(rv_in$de,
          input$SI_comp_fact,
          input$RB_select_test)

      diff_exp <- isolate(rv_in$de)
      # TODO: change to data.frame
      # filter according to current "cases"
      de <- diff_exp %>%

        dplyr::filter(test_type == input$RB_select_test &
                        versus == input$SI_comp_fact) %>%

        dplyr::mutate(f=1,
                      significant = pvals_adj < 0.01,
                      point_color =  dplyr::recode(as.character(.data$significant), "TRUE" = "#FF7F00", "FALSE" = "#1F78B4"))

      return(de)
    })


    # create some ui output
    output$UI_comp_group_selection <- renderUI({
      req(rv_in$config,
          rv_in$default,
          rv_in$de)

      cfg <- isolate(rv_in$config)
      test_choices <- levels(factor(isolate(rv_in$de$test_type)))


      to_return <-  tagList(
        radioButtons(inputId = ns("RB_select_test"),
                     label = "sig test:",
                     choices = test_choices,
                     selected = rv_in$default$test[1]),
        selectInput(inputId = ns("SI_comp_fact"),
                    label = "compare: ",
                    choices =  strsplit(cfg[ID=="diff_exp_comps"]$fID, "\\|")[[1]],
                    selected = rv_in$default$comp_fact),

        shinyjs::disabled(selectInput(inputId = ns("SI_color_grp"),
                    label = "Color by:",
                    choices = rv_in$config[grp == TRUE]$UI,
                    selected = rv_in$default$color_grp)
        ),
        checkboxInput(inputId = ns("CB_drop_bottom"),
                               label = "drop non-sigs?:",
                               value = TRUE)

      )

      return(to_return)
    })



    observeEvent(input$SI_comp_type, {
      req(rv_in$config,
          rv_in$default)
      # get the groups
      cfg <- isolate(rv_in$config)
      subs <- strsplit(cfg[UI == "diff_exp_comps"]$fID, "\\|")[[1]]
      subs <- sort(subs)

      #hack to sort the names with numbers
      subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
      if (all(!is.na(subs2))){
        names(subs2) <- seq_along(subs2)
        subs <- subs[as.numeric(names(sort(subs2)))]
      }

      #subs_label <- paste0("select ",isolate(input$SI_subset),"s: ")

      updateSelectInput(inputId = "SI_comp_fact",
                        label = "compare :",
                        choices = subs,
                        selected = rv_in$default$comp_fact[1])


    })


    output$volcano_plot <- renderPlotly({
      req(filtered_de, #input$SI_comp_type,
          input$SI_comp_fact)

      colorby_group <- input$SI_color_grp
      de <- filtered_de()


      # title_text <- paste0(input$SI_comp_type, " || ", input$SI_comp_fact)
      if( dim(de)[1]>0 ) {
        # pg_volcano_ly(in_data = filtered_de(),
        #              pvalue_adjust = input$test_cor_pvalue,
        #              title = paste0(input$SI_comp_type, " vs ", input$test_group2))
        title_text <- paste0("comp: ", input$SI_comp_fact)


        # TODO:  make this a function and move to fct_pg_vis_comp
        # volc <- pg_volc_ly(de=de_local , title = title_text)
        # return(volc)
        volc <- plot_ly(
          x = de$logfoldchange,
          y = -log10(de$pvals),
          name = "FDR > 0.05",
          type = "scatter",
          showlegend = FALSE,
          mode = "markers",
          # Hovertext
          text = paste(de$names,
                       "</br></br>Beta: ",
                       format( de$logfoldchanges, digits = 3, scientific = TRUE),
                       " (score: ",
                       format( de$scores, digits = 3, scientific = TRUE),
                       "</br>Q-value: ",
                       format(de$pvals_adj, digits = 3, scientific = TRUE)),
          hoverinfo = "text",
          color = ~I(de$point_color) )

        volc <- volc %>%
          # Adding markers for a custom legend.  Technically,
          # the entire volcano plot trace is colored blue,
          # but we need a legend to indicate the meaning of the orange points,
          # so we add traces with orange and blue and relabel.
          # It's hacky but it works better for animation and plotly_click purposes.

          # Blue/not significant
          plotly::add_markers(x= 0.8, y = 6.5, color = I("#1F78B4"), showlegend = FALSE, hoverinfo = "skip") %>%
          plotly::add_annotations(x=0.8, y=6.5, xref = "x", yref = "y", text = "FDR > 0.01",
                                  xanchor = 'left', showarrow = F, xshift = 10) %>%
          # Orange/significant
          plotly::add_markers(x= 0.8, y = 7, color = I("#FF7F00"), showlegend = FALSE, hoverinfo = "skip") %>%
          plotly::add_annotations(x=0.8, y=7, xref = "x", yref = "y", text = "FDR < 0.01",
                                  xanchor = 'left', showarrow = F, xshift = 10) %>%

          plotly::layout(
            title = title_text,
            xaxis = list(title = "Effect (logFC)", range = c(-4, 4)),
            yaxis = list(title = "-log10 p-value", range = c(-0.1, 10.25))
          ) %>%
          # Disable the legend click since our traces do not correspond to the
          # actual legend labels
          htmlwidgets::onRender("function(el,x){el.on('plotly_legendclick', function(){ return false; })}") %>%
          plotly::config(displayModeBar = FALSE)

      }

      return(volc)

    })

    ### Set up data for boxplots
    omic_by_subject <- reactive({
      selected <- event_data(event = "plotly_click")$pointNumber + 1
      scaled_counts <- scale(rv_in$ad$X[,selected])
      grouping_var <- filtered_de()$obs_name[1]
      omic_counts <- data.frame( rv_in$ad$obs_names,
                                 scaled_counts,
                                 rv_in$ad$obs[[grouping_var]] )
      colnames(omic_counts) <- c("id", "val","grp")
      #metab_counts$count <- log10(metab_counts$count)
      return(omic_counts)
    })


    output$test_boxplot <- renderPlotly({
      validate(
        need(
          event_data(event = "plotly_click") & event_data(event = "plotly_click")$curveNumber == 0,
          "Select a metabolite on the volcano plot above to view distribution among statuses."
        )
      )
      # TODO:  make this into a boxplot function
      plot_ly(
        data = omic_by_subject(),
        y = ~val,
        color = ~grp,
        type = "box",
        boxpoints = "all",
        pointpos = 0,
        text = paste("Transformed Reading: ",
                     format(omic_by_subject()$val, digits = 3)),
        hoverinfo = list("median", "text"),
        colors = "Dark2"
      ) %>%
        plotly::layout(
          title = filtered_de()[event_data(event = "plotly_click")$pointNumber + 1,1],
          yaxis = list(title = "Transformed Reading",
                       hoverformat = ".2f"),
          showlegend = FALSE
        ) %>%
        plotly::config(displayModeBar = FALSE)
    })

    # Information on omics selected on volcano plot
    # TODO:  collect the meta-info and render to HTML for visualization...
    output$UI_meta_box <- renderText({
      validate(
        need(
          event_data(event = "plotly_click") & event_data(event = "plotly_click")$curveNumber == 0,
          "Select an omic on the volcano plot above to view distribution among statuses."
        )
      )
      # PUT rv_in$ad$obs... or something here
      # om_name <- volcano_data()[event_data(event = "plotly_click")$pointNumber + 1, 1]
      metab_name = "dummy"
      paste(
        "<b>Metabolite</b>: ", metab_name, "</br>"
        # "<b>Super Pathway</b>: ", metab_meta$`Super Pathway`[metab_meta$Metabolite == metab_name], "</br>",
        # "<b>Sub Pathway:</b> ", metab_meta$`Sub Pathway`[metab_meta$Metabolite == metab_name],
        # # Link to HMDB, if ID is available
        # ifelse(
        #   !is.na(metab_meta$`HMDB ID`[metab_meta$Metabolite == metab_name]),
        #   paste(
        #     "</br><b>HMDB ID:</b> ",
        #     metab_meta$`HMDB ID`[metab_meta$Metabolite == metab_name]
        #   ),
        #   ""
        # ),
        # # Link to PubChem, if ID is available
        # ifelse(
        #   !is.na(metab_meta$`PubChem ID`[metab_meta$Metabolite == metab_name]),
        #   paste("</br><b>PubChem ID:</b> ",
        #         metab_meta$`PubChem ID`[metab_meta$Metabolite == metab_name]
        #   ),
        #   ""
        # ),
        # ifelse(
        #   !is.na(metab_name),
        #   paste(
        #     # "</br><b>ADAD vs CO effect (q):</b> &nbsp&nbsp&nbsp",
        #     # format(volcano_data_adadvsco$effects[volcano_data_adadvsco$metabs == metab_name], digits = 3, scientific = TRUE),
        #     # " (", format(volcano_data_adadvsco$pvals_adj[volcano_data_adadvsco$metabs == metab_name], digits = 3, scientific = TRUE), ")", "</br>",
        #     # "<b>AD vs CO effect (q)</b>: &nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp",
        #     # format(volcano_data_cavsco$effects[volcano_data_cavsco$metabs == metab_name], digits = 3, scientific = TRUE),
        #     # " (", format(volcano_data_cavsco$pvals_adj[volcano_data_cavsco$metabs == metab_name], digits = 3, scientific = TRUE), ")", "</br>",
        #     # "<b>TREM2 vs CO effect (q)</b>: &nbsp",
        #     # format(volcano_data_trem2vsco$effects[volcano_data_trem2vsco$metabs == metab_name], digits = 3, scientific = TRUE),
        #     # " (", format(volcano_data_trem2vsco$pvals_adj[volcano_data_trem2vsco$metabs == metab_name], digits = 3, scientific = TRUE), ")", "</br>",
        #     # # "<b>ADAD vs AD effect (q)</b>: &nbsp&nbsp&nbsp",
        #     # # format(volcano_data_adadvsca$effects[volcano_data_adadvsca$metabs == metab_name], digits = 3, scientific = TRUE),
        #     # # " (", format(volcano_data_adadvsca$pvals_adj[volcano_data_adadvsca$metabs == metab_name], digits = 3, scientific = TRUE), ")",
        #     sep = ""),
        #   ""
        # )
      )
    })


  })
}



## To be copied in the UI
# mod_pg_vis_comp_ui("pg_vis_comp_ui_1")

## To be copied in the server
# mod_pg_vis_comp_server("pg_vis_comp_ui_1")

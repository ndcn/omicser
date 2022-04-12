
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
mod_pg_diff_expr_ui <- function(id){
  ns <- NS(id)
  tagList(


    HTML("Volcano plot"),
    h4("Relative expression of -omic features"),
    "In this tab, users can visualise the expression of one -omic feature in one experimental condition relative to another condition (or baseline set of conditions).",
    "(i.e. log2(Fold Change) vs -log10(p-value).",
    br(),br(),
    fluidRow(
      column(
        2, style="border-right: 2px solid black",
        br(),

        # this dynamicall makes all of our test/group selectors
        uiOutput(outputId = ns("UI_comp_group_selection")),
        # select which p-value you want
        radioButtons(inputId = ns("RB_select_pval"),
                     label = "Select p-value:",
                     choices = c("raw p-value" = "pvals",
                                 "adjusted p-value" = "pvals_adj"),
                     selected = "pvals")

      ), # End of column (6 space)
      column(10,
             #uiOutput(ns("UI_viz_output")),
             plotlyOutput(outputId = ns("volcano_plot"),height = "800px")
      )
    ),
    fluidRow(
      box(id = ns("omic_boxplot"),
          title = "Expression Distribution",
          solidHeader = TRUE,
          status = "primary",
          plotlyOutput(outputId = ns("test_boxplot")),
          width = 8),
      # Info on metab selected on volcano plot
      box(id = ns("omic_info"),
          title = "Feature Meta-info",
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
#' @param id shiny internal
#' @param rv_data reactive from ingest
#' @param rv_selections reactives from side select
#' @param active_layer_data current data matrix
#'
#' @noRd
#'
#' @importFrom plotly renderPlotly plotlyOutput plot_ly add_markers event_data
#' @importFrom utils head
#'
mod_pg_diff_expr_server <- function(id,rv_data, rv_selections, active_layer_data){
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
      req(rv_data$de,
          input$SI_comp_fact,
          input$RB_select_test)

      # add coding for target genes to be red...
      diff_exp <- as.data.table(rv_data$de)
      # # TODO: change to data.frame
      # # filter according to current "cases"
      #filter
      de <- diff_exp[test_type == input$RB_select_test & versus == input$SI_comp_fact][
        , f := 1][
          , significant :=lapply(.SD, function(x) x < 0.01), .SDcols = input$RB_select_pval][
            , point_color := ifelse(significant, "#FF7F00", "#1F78B4")][
              , neglogpval := lapply(.SD, function(x) -log10(x)), .SDcols = input$RB_select_pval]

      selected_ <- de$names %in% rv_selections$selected_omics$target_omics
      de$point_color[selected_] <- "#FF6347"

      # de <- diff_exp %>%
      #   dplyr::filter(test_type == input$RB_select_test &
      #                    versus == input$SI_comp_fact) %>%
      #   dplyr::mutate(f=1,
      #                 significant = pvals_adj < 0.01,
      #                 point_color =  dplyr::recode(as.character(.data$significant), "TRUE" = "#FF7F00", "FALSE" = "#1F78B4"))

      return(de)
    })


    # create some ui output
    output$UI_comp_group_selection <- renderUI({
      req(rv_data$config,
          rv_data$default,
          rv_data$de)

      cfg <- rv_data$config
      test_choices <- levels(factor(rv_data$de$test_type))
      #comp_types <- levels(factor(isolate(rv_data$de$comp_type)))

      to_return <-  tagList(
        radioButtons(inputId = ns("RB_select_test"),
                     label = "sig test:",
                     choices = test_choices,
                     selected = rv_data$default$test[1]),
        # TODO: choose comp_type first and then have a simplified comparison set
        # selectInput(inputId = ns("SI_comp_type"),
        #             label = "compare: ",
        #             choices =  comp_types,
        #             selected = comp_types[1],

        selectInput(inputId = ns("SI_comp_fact"),
                    label = "compare: ",
                    choices =  strsplit(cfg[ID=="diff_exp_comps"]$fID, "\\|")[[1]],
                    selected = rv_data$default$comp_fact)

      )

      return(to_return)
    })



    observeEvent(input$SI_comp_type, {
      req(rv_data$config,
          rv_data$default)
      # get the groups
      cfg <- rv_data$config
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
                        selected = rv_data$default$comp_fact[1])


    })


    # make a waiter to signal when the volcano is rendering
    vwtr <- waiter::Waiter$new(
      id = ns("volcano_plot"),
      html = waiter::spin_loaders(id=12,color="red"),
      color = waiter::transparent(.5)  #waiter::transparent(.5)
    )


    output$volcano_plot <- renderPlotly({
      req(filtered_de, #input$SI_comp_type,
          input$SI_comp_fact)

      vwtr$show()

      de <- filtered_de()
      y_leg <- max(de$neglogpval)
      x_leg <- max(de$logfoldchanges)

      # title_text <- paste0(input$SI_comp_type, " || ", input$SI_comp_fact)
      if( dim(de)[1]>0 ) {
        # pg_volcano_ly(in_data = filtered_de(),
        #              pvalue_adjust = input$test_cor_pvalue,
        #              title = paste0(input$SI_comp_type, " vs ", input$test_group2))
        title_text <- paste0("differential expression - ", input$SI_comp_fact)

        # group - the comparison {names}V{reference}
        # names - what are we comparing?
        # obs_name - name of the meta data variable
        # test_type - what statistic are we using - reference - the denomenator. or the condition we are comparing expressions values to
        # comp_type - grpVref or grpVrest. rest is all other conditions
        # logfoldchanges - log2(name/reference)
        # scores - statistic score
        # pvals - pvalues from the stats test. e.g. t-test
        # pvals_adj - adjusted pvalue (Q)
        # versus - label which we will choose in the browser
        ### ADDED in filtered_de()
        # significant
        # point_color

        #TODO: color the target omics a different color

        # TODO:  make this a function and move to fct_pg_vis_comp
        # volc <- pg_volc_ly(de=de_local , title = title_text)
        # return(volc)
        volc <- plot_ly(
          x = de$logfoldchanges,
          y = de$neglogpval,
          name = "FDR < 0.01",
          type = "scatter",
          showlegend = FALSE,
          mode = "markers",
          alpha = 0.5,
          # Hovertext
          text = paste(de$names,
                       "</br></br>log2FC: ",
                       format( de$logfoldchanges, digits = 3, scientific = TRUE),
                       " (score: ",
                       format( de$scores, digits = 3, scientific = TRUE),
                       "</br>Q-value: ",
                       format(de[ , input$RB_select_pval, with = FALSE], digits = 3, scientific = TRUE)),
          hoverinfo = "text",
          color = ~I(de$point_color) )

        # make the y-label
        y_label_vol_plot <- switch(input$RB_select_pval,
                                   "pvals" = "-log10 p-value",
                                   "pvals_adj" = "-log10 adjusted p-value")

        volc <- volc %>%
          # Adding markers for a custom legend.  Technically,
          # the entire volcano plot trace is colored blue,
          # but we need a legend to indicate the meaning of the orange points,
          # so we add traces with orange and blue and relabel.
          # It's hacky but it works better for animation and plotly_click purposes.

          # Blue/not significant
          plotly::add_markers(x= 0.75*x_leg ,
                              y = 0.80*y_leg,
                              color = I("#1F78B4"), showlegend = FALSE, hoverinfo = "skip") %>%
          plotly::add_annotations(x= 0.75*x_leg ,
                                  y = 0.80*y_leg,
                                  xref = "x", yref = "y", text = "FDR > 0.01",
                                  xanchor = 'left', showarrow = F, xshift = 10) %>%
          # Orange/significant
          plotly::add_markers(x= 0.75*x_leg ,
                              y = 0.75*y_leg,
                              color = I("#FF7F00"), showlegend = FALSE, hoverinfo = "skip") %>%
          plotly::add_annotations(x= 0.75*x_leg ,
                                  y = 0.75*y_leg,
                                  xref = "x", yref = "y", text = "FDR < 0.01",
                                  xanchor = 'left', showarrow = F, xshift = 10) %>%
          # red target features
          plotly::add_markers(x= 0.75*x_leg,
                              y = 0.85*y_leg,
                              color = I("#FF6347"), showlegend = FALSE, hoverinfo = "skip") %>%
          plotly::add_annotations(x= 0.75*x_leg ,
                                  y = 0.85*y_leg,
                                  xref = "x", yref = "y", text = "target feature",
                                  xanchor = 'left', showarrow = F, xshift = 10) %>%

          # plotly::add_trace(
          #   x = de[point_color == "#FF6347"]$logfoldchanges,
          #   y = de[point_color == "#FF6347"]$neglogpval,
          #   marker = list(
          #     color = "#FF6347",
          #     size = 10,
          #     line = list(
          #       color = 'rgb(231, 99, 250)',
          #       width = 1
          #     )
        #   ),
        #   showlegend = FALSE) %>%

        plotly::layout(
          title = list(text = title_text,
                       x = 0),
          xaxis = list(zeroline = FALSE,
                       title = "log2FC"), #, range = c(-8, 8)),
          yaxis = list(title = y_label_vol_plot), #, range = c(-0.1, 60.25))
          shapes = list(vline(-1),
                        vline(1),
                        hline(-log10(0.01))),
          legend = list(orientation = "h",x = 0.1, y = 0.9)) %>%


          # plotly::layout(
          #   title = title_text,
          #   xaxis = list(title = "Effect (logFC)", range = c(-8, 8)),
          #   yaxis = list(title = "-log10 p-value", range = c(-0.1, 60.25))
          # ) %>%
          # Disable the legend click since our traces do not correspond to the
          # actual legend labels
          htmlwidgets::onRender("function(el,x){el.on('plotly_legendclick', function(){ return false; })}") %>%
          plotly::config(displayModeBar = TRUE,
                         toImageButtonOptions = list(format= 'png', # one of png, svg, jpeg, webp
                                                  filename= 'diff_expr_volcano',
                                                  height= 600,
                                                  width= 700,
                                                  scale= 1 )) %>%

          plotly::event_register('plotly_click')


      } else {
        volc <- NULL
      }

      return(volc)

    })

    ### Set up data for boxplots
    ### # TODO: rename e.g. omic_expr_values (in mod_table)
    omic_by_subject <- reactive({
      grouping_var <- filtered_de()$obs_name[1]

      selected <- event_data(event = "plotly_click")$pointNumber + 1
      omic_name <- filtered_de()$names[selected]

      data_vec <- active_layer_data$data[,omic_name]
      # check if there is data present
      validate(
        need(
          !is.null(data_vec),
          "Please click **Plot now!** to load the data for the boxplot!"
        )
      )

      # to be sure that this is a character and not a factor
      grouping_var <- as.character(filtered_de()$obs_name[1])
      omic_counts <- data.frame( rv_data$anndata$obs_names,
                                 data_vec,
                                 rv_data$anndata$obs[[grouping_var]] )
      colnames(omic_counts) <- c("id", "val","grp")

      return(omic_counts)
    })


    output$test_boxplot <- renderPlotly({
      validate(
        need(
          event_data(event = "plotly_click") & event_data(event = "plotly_click")$curveNumber == 0,
          "Select an *-omic feature on the volcano plot above to view distribution among experimental factors."
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
        text = paste("Vals: ",
                     format(omic_by_subject()$val, digits = 3)),
        hoverinfo = list("median", "text"),
        colors = "Dark2",
        source = "B"
      ) %>%
        plotly::layout(
          title = filtered_de()$names[event_data(event = "plotly_click")$pointNumber + 1],
          yaxis = list( title = "value",
                        hoverformat = ".2f" ),
          xaxis = list( title = filtered_de()$obs_name[1] ),
          showlegend = FALSE
        ) %>%
        #plotly::config(displayModeBar = FALSE)
        plotly::config(displayModeBar = TRUE,
                       toImageButtonOptions = list(format= 'png', # one of png, svg, jpeg, webp
                                                 filename= 'expr_box',
                                                 height= 300,
                                                 width= 600,
                                                 scale= 1 ))
    })

    # Information on omics selected on volcano plot
    # TODO:  collect the meta-info and render to HTML for visualization...
    output$UI_meta_box <- renderText({
      validate(
        need(
          event_data(event = "plotly_click") & event_data(event = "plotly_click")$curveNumber == 0,
          "Select an *-omic feature on the volcano plot above to see related meta-info"
        )
      )



      selected <- event_data(event = "plotly_click")$pointNumber + 1
      omic_name <- filtered_de()$names[selected]

      # TODO: replace.
      #omic_details <- rv_data$default$omic_details
      omic_details <- rv_data$shaddow_defs$feat_deets

      annots <- rv_data$anndata$var[omic_name,]

      deet <- paste("<b>", rv_data$db_meta$omics_type, " </b> features </br>")

      for (detail in omic_details) {
        value <- annots[[detail]]
        if (is.numeric(value) ) {
          value <- format(value, digits=3,scientific=TRUE)
        }
        deet <- paste(deet, "<b>", detail, "</b>:", value, "</br> ")
      }

      deet

    })

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



  })
}



## To be copied in the UI
# mod_pg_diff_expr_ui("pg_diff_expr_ui_1")

## To be copied in the server
# mod_pg_diff_expr_server("pg_diff_expr_ui_1")

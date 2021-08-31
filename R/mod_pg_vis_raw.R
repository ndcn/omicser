
#' pg_vis_raw UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_pg_vis_raw_ui <- function(id){
  ns <- NS(id)
  tagList(
    HTML("Violinplot / Boxplot"),
    h4("Cell information / omic expression violin plot / box plot"),
    "Here we visualise the expression or continuous cell information ",
    br(),br(),
    fluidRow(
      column(
        3, style = "border-right: 2px solid black",
        # selectInput("Van_c1inp1", "Cell information (X-axis):",
        #             choices = Van_conf[grp == TRUE]$UI,
        #             selected = Van_def$grp1) ,
        # # %>%
        #   helper(type = "inline", size = "m", fade = TRUE,
        #          title = "Cell information to group cells by",
        #          content = c("Select categorical cell information to group cells by",
        #                      "- Single cells are grouped by this categorical covariate",
        #                      "- Plotted as the X-axis of the violin plot / box plot")),
        # selectInput("Van_c1inp2", "measures", choices=NULL),# %>%
          # helper(type = "inline", size = "m", fade = TRUE,
          #        title = "Cell Info / Gene to plot",
          #        content = c("Select cell info / gene to plot on Y-axis",
          #                    "- Can be continuous cell information (e.g. nUMIs / scores)",
          #                    "- Can also be gene expression")),
        radioButtons(ns("RB_dist_plot_type"), "Plot type:",
                     choices = c("violin", "boxplot"),
                     selected = "violin", inline = TRUE),
        checkboxInput(ns("CB_show_data_points"), "Show data points", value = FALSE),
        # actionButton("Van_c1togL", "Toggle to subset cells"),
        # conditionalPanel(
        #   condition = "input.Van_c1togL % 2 == 1",
        #   selectInput("Van_c1sub1", "Cell information to subset:",
        #               choices = Van_conf[grp == TRUE]$UI,
        #               selected = Van_def$grp1),
        #   uiOutput("Van_c1sub1.ui"),
        #   actionButton("Van_c1sub1all", "Select all groups", class = "btn btn-primary"),
        #   actionButton("Van_c1sub1non", "Deselect all groups", class = "btn btn-primary")
        # ), br(),
        br()
        # actionButton("Van_c1tog", "Toggle graphics controls"),
        # conditionalPanel(
        #   condition = "input.Van_c1tog % 2 == 1",
        #   sliderInput("Van_c1siz", "Data point size:",
        #               min = 0, max = 4, value = 1.25, step = 0.25),
        #   radioButtons("Van_c1psz", "Plot size:",
        #                choices = c("Small", "Medium", "Large"),
        #                selected = "Medium", inline = TRUE),
        #   radioButtons("Van_c1fsz", "Font size:",
        #                choices = c("Small", "Medium", "Large"),
        #                selected = "Medium", inline = TRUE))
      ), # End of column (6 space)
      column(9,
             uiOutput(ns("UI_viz_output")#,
             # downloadButton("Van_c1oup.pdf", "Download PDF"),
             # downloadButton("Van_c1oup.png", "Download PNG"), br(),
             # div(style="display:inline-block",
             #     numericInput("Van_c1oup.h", "PDF / PNG height:", width = "138px",
             #                  min = 4, max = 20, value = 8, step = 0.5)),
             # div(style="display:inline-block",
             #     numericInput("Van_c1oup.w", "PDF / PNG width:", width = "138px",
             #                  min = 4, max = 20, value = 10, step = 0.5))
      )  # End of column (6 space)
    )    # End of fluidRow (4 space)

   ),

   HTML("Bubbleplot / Heatmap"),
   h4("Gene expression bubbleplot / heatmap"),
   "In this tab, users can visualise the expression patterns of ",
   "multiple omics grouped by categorical cell information (e.g. library / cluster).", br(),
   "The normalised expression values are group averaged and scaled/thresholded (?)).",
   br(),br(),
   fluidRow(
     column(
       3, style="border-right: 2px solid black",
       # textAreaInput("Van_d1inp", HTML("List of gene names <br />
       #                                    (Max 50 genes, separated <br />
       #                                     by , or ; or newline):"),
       #               height = "200px",
       #               value = paste0(Van_def$genes, collapse = ", ")) %>%
       #   helper(type = "inline", size = "m", fade = TRUE,
       #          title = "List of genes to plot on bubbleplot / heatmap",
       #          content = c("Input genes to plot",
       #                      "- Maximum 50 genes (due to ploting space limitations)",
       #                      "- Genes should be separated by comma, semicolon or newline")),
       # selectInput("Van_d1grp", "Group by:",
       #             choices = Van_conf[grp == TRUE]$UI,
       #             selected = Van_conf[grp == TRUE]$UI[1]), #%>%
       #   # helper(type = "inline", size = "m", fade = TRUE,
       #   #        title = "Cell information to group cells by",
       #   #        content = c("Select categorical cell information to group cells by",
       #   #                    "- Single cells are grouped by this categorical covariate",
       #   #                    "- Plotted as the X-axis of the bubbleplot / heatmap")),
       #   #
       radioButtons(ns("RB_heat_plot_type"), "Plot type:",
                    choices = c("Bubbleplot", "Heatmap"),
                    selected = "Bubbleplot", inline = TRUE),
       checkboxInput(ns("CB_scale"), "Scale omic expression", value = TRUE),
       checkboxInput(ns("CB_cluster_rows"), "Cluster rows (omics)", value = TRUE),
       checkboxInput(ns("CB_cluster_cols"), "Cluster columns (samples)", value = FALSE),
       br(),
     #   actionButton("Van_d1togL", "Toggle to subset cells"),
     #   conditionalPanel(
     #     condition = "input.Van_d1togL % 2 == 1",
     #     selectInput("Van_d1sub1", "Cell information to subset:",
     #                 choices = Van_conf[grp == TRUE]$UI,
     #                 selected = Van_def$grp1),
     #     uiOutput("Van_d1sub1.ui"),
     #     actionButton("Van_d1sub1all", "Select all groups", class = "btn btn-primary"),
     #     actionButton("Van_d1sub1non", "Deselect all groups", class = "btn btn-primary")
     #   ), br(), br(),
     #   actionButton("Van_d1tog", "Toggle graphics controls"),
     #   conditionalPanel(
     #     condition = "input.Van_d1tog % 2 == 1",
     #     radioButtons("Van_d1cols", "Colour scheme:",
     #                  choices = c("White-Red", "Blue-Yellow-Red",
     #                              "Yellow-Green-Purple"),
     #                  selected = "Blue-Yellow-Red"),
     #     radioButtons("Van_d1psz", "Plot size:",
     #                  choices = c("Small", "Medium", "Large"),
     #                  selected = "Medium", inline = TRUE),
     #     radioButtons("Van_d1fsz", "Font size:",
     #                  choices = c("Small", "Medium", "Large"),
     #                  selected = "Medium", inline = TRUE))
     ), # End of column (6 space)
     column(9,
            h4(htmlOutput("HTML_header")),
            uiOutput(ns("UI_heatmap"))
            # downloadButton("Van_d1oup.pdf", "Download PDF"),
            # downloadButton("Van_d1oup.png", "Download PNG"), br(),
            # div(style="display:inline-block",
            #     numericInput("Van_d1oup.h", "PDF / PNG height:", width = "138px",
            #                  min = 4, max = 20, value = 10, step = 0.5)),
            # div(style="display:inline-block",
            #     numericInput("Van_d1oup.w", "PDF / PNG width:", width = "138px",
            #                  min = 4, max = 20, value = 10, step = 0.5))
     )  # End of column (6 space)
   )    # End of fluidRow (4 space)

  )  #taglist
}

#' pg_vis_raw Server Functions
#'
#' @noRd
mod_pg_vis_raw_server <- function(id, rv_in, p, heat_data, box_data){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    # # omics_list
    #   feat_grp
    #   feat_subsel
    #   observ_grpA
    #   observ_subselA
    #   observ_grpB
    #   observ_subselB
    #   observ_x
    #   observ_y
    #   measure_type
    #   raw_plot_type
    #   comp_plot_type
    #   obs_type
    # #

    data_point_sz = .65
    plot_size = "Small"
    font_size = "Small"
    color_schemes = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple")
    color_scheme = color_schemes[2]

    #go_signal <- reactive( p$omics_list$viz )

    output$plot_box_out <- renderPlot({
      #req(rv_in$ad)
      req( p$omics_list)

      # box_data <- list(
      #   x_name = in_fact,
      #   y_name = dat_key,
      #   dat_source = dat_source,
      #   data = bx_data,
      # )
      plot_type <- input$RB_dist_plot_type
      show_data_points <- input$CB_show_data_points

      notch <- TRUE
      data_point_sz = .65
      plot_size = "Small"
      font_size = "Small"
      grp <- ifelse(p$group_action=="none",FALSE,TRUE)

      vplot <- violin_box(box_data$data, box_data$x_name,box_data$y_name,
                          box_data$colors,
                          plot_type,
                          show_data_points,
                          data_point_sz,
                          font_size,
                          notch,
                          grp)


      # #TODO: make these interactive? (e.g. shinycell)
      # data_point_sz = .65
      # plot_size = "Small"
      # font_size = "Small"
      # color_schemes = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple")
      # color_scheme = color_schemes[2]
      # in_quant <- dat_key #(maybe) just observ_y
      # pg_violin_box(in_conf, in_meta, in_fact, in_quant,
      #               in_grp, in_subset, in_subsel,
      #               in_data, in_omic, plot_type,
      #               show_data_points, data_point_sz, font_size )
      return(vplot)
    })

    output$UI_viz_output <- renderUI({
      #if (p$omics_list$viz_now) {
        plotOutput(ns("plot_box_out"), height = pList2[plot_size])
      #}
    })


    output$HTML_header <- renderUI({
      print("renderUI HTML header")

      omic_list = p$omics_list$value
      if(nrow(omic_list) > 50){
        HTML("More than 50 input genes! Please reduce the gene list!")
      } else {
        oup = paste0(nrow(omic_list[present == TRUE]), " genes OK and will be plotted")
        if(nrow(omic_list[present == FALSE]) > 0){
          oup = paste0(oup, "<br/>",
                       nrow(omic_list[present == FALSE]), " genes not found (",
                       paste0(omic_list[present == FALSE]$omic, collapse = ", "), ")")
        }
        HTML(oup)
      }
    })


    # plot_heatmap_out renderPlot---------------------------------
    output$plot_heatmap_out <- renderPlot({

      # hm_data <- list(
      #   x_name = X_fact,
      #   y_name = group_by$y,
      #   x_source = x_is,
      #   type = dat_loc,
      #   data = hm_data,
      #   mat = X_data,
      #   meta = tmp_meta,
      #   ready = TRUE
      # )

      print("in plot_heatmap_out")


      input_data <- heat_data$data
      #
      x_names <- heat_data$x_names
      y_names <- heat_data$y_names

browser()
      plot_type <- input$RB_heat_plot_type
      in_do_scale <- input$CB_scale
      in_clust_row <- input$CB_cluster_rows
      in_clust_col <- input$CB_cluster_cols

      plot_size = "Small"
      color_schemes = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple")
      color_scheme = color_schemes[2]

      grp <- ifelse(p$group_action=="none",FALSE,TRUE)

      hmap <- bubble_heatmap(input_data, x_names, y_names, plot_type,
                                 in_do_scale, in_clust_row, in_clust_col,
                                 color_scheme, plot_size, grp, save = FALSE)
#
#       #TODO: make these interactive? (e.g. shinycell)
#       plot_size = "Small"
#       color_schemes = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple")
#       color_scheme = color_schemes[2]
#       pg_bubble_heatmap(in_conf, in_meta, in_omics, in_fact, plot_type,
#                         in_grp, in_subset,in_subsel,
#                         in_data, all_omics,
#                         input$CB_scale, input$CB_cluster_rows, input$CB_cluster_cols,
#                         color_scheme, plot_size)

    return(hmap)
    })

    output$UI_heatmap <- renderUI({
      #if (p$omics_list$viz_now) {
        print("renderUI UI heatmap")

        plotOutput(ns("plot_heatmap_out"), height = pList3[plot_size])

    })

    # output$Van_d1oup.pdf <- downloadHandler(
    #   filename = function() { paste0("Van_",input$Van_d1plt,"_",input$Van_d1grp,".pdf") },
    #   content = function(file) { ggsave(
    #     file, device = "pdf", height = input$Van_d1oup.h, width = input$Van_d1oup.w,
    #     plot = pg_bubble_heatmap(Van_conf, Van_meta, input$Van_d1inp, input$Van_d1grp, input$Van_d1plt,
    #                       input$Van_d1sub1, input$Van_d1sub2, "Van_gexpr.h5", Van_gene,
    #                       input$Van_d1scl, input$Van_d1row, input$Van_d1col,
    #                       input$Van_d1cols, input$Van_d1fsz, save = TRUE) )
    #   })
    # output$Van_d1oup.png <- downloadHandler(
    #   filename = function() { paste0("Van_",input$Van_d1plt,"_",input$Van_d1grp,".png") },
    #   content = function(file) { ggsave(
    #     file, device = "png", height = input$Van_d1oup.h, width = input$Van_d1oup.w,
    #     plot = pg_bubble_heatmap(Van_conf, Van_meta, input$Van_d1inp, input$Van_d1grp, input$Van_d1plt,
    #                       input$Van_d1sub1, input$Van_d1sub2, "Van_gexpr.h5", Van_gene,
    #                       input$Van_d1scl, input$Van_d1row, input$Van_d1col,
    #                       input$Van_d1cols, input$Van_d1fsz, save = TRUE) )
    #   })


  })
}

## To be copied in the UI
# mod_pg_vis_raw_ui("pg_vis_raw_ui_1")

## To be copied in the server
# mod_pg_vis_raw_server("pg_vis_raw_ui_1")

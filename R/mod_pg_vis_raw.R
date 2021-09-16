
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
    h4("Cell information violin plot / box plot"),
    "Here we visualise the marginal expression values ",
    br(),br(),
    fluidRow(
      column(
        3, style = "border-right: 2px solid black",
        radioButtons(ns("RB_dist_plot_type"), "plot type:",
                     choices = c("violin", "boxplot"),
                     selected = "violin", inline = TRUE),
        checkboxInput(ns("CB_show_data_points"), "Show data points", value = FALSE),
               br()
        # actionButton("Van_c1tog", "Toggle graphics controls"),
        # conditionalpanel(
        #   condition = "input.Van_c1tog % 2 == 1",
        #   sliderInput("Van_c1siz", "Data point size:",
        #               min = 0, max = 4, value = 1.25, step = 0.25),
        #   radioButtons("Van_c1psz", "plot size:",
        #                choices = c("Small", "Medium", "Large"),
        #                selected = "Medium", inline = TRUE),
        #   radioButtons("Van_c1fsz", "Font size:",
        #                choices = c("Small", "Medium", "Large"),
        #                selected = "Medium", inline = TRUE))
      ), # End of column (6 space)
      column(9,
             uiOutput(ns("UI_box_output")#,
             # downloadButton("Van_c1oup.pdf", "Download pDF"),
             # downloadButton("Van_c1oup.png", "Download pNG"), br(),
             # div(style="display:inline-block",
             #     numericInput("Van_c1oup.h", "pDF / pNG height:", width = "138px",
             #                  min = 4, max = 20, value = 8, step = 0.5)),
             # div(style="display:inline-block",
             #     numericInput("Van_c1oup.w", "pDF / pNG width:", width = "138px",
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
       radioButtons(ns("RB_heat_plot_type"), "plot type:",
                    choices = c("Bubbleplot", "Heatmap"),
                    selected = "Bubbleplot", inline = TRUE),
       checkboxInput(ns("CB_scale"), "Scale omic expression", value = TRUE),
       checkboxInput(ns("CB_cluster_rows"), "Cluster rows (omics)", value = TRUE),
       checkboxInput(ns("CB_cluster_cols"), "Cluster columns (samples)", value = FALSE),
       br(),

     ), # End of column (6 space)
     column(9,
            h4(htmlOutput("HTML_header")),
            uiOutput(ns("UI_heatmap"))
            # downloadButton("Van_d1oup.pdf", "Download pDF"),
            # downloadButton("Van_d1oup.png", "Download pNG"), br(),
            # div(style="display:inline-block",
            #     numericInput("Van_d1oup.h", "pDF / pNG height:", width = "138px",
            #                  min = 4, max = 20, value = 10, step = 0.5)),
            # div(style="display:inline-block",
            #     numericInput("Van_d1oup.w", "pDF / pNG width:", width = "138px",
            #                  min = 4, max = 20, value = 10, step = 0.5))
     )  # End of column (6 space)
   ),    # End of fluidRow (4 space)
   HTML("Violinplot / Boxplot"),
   h4("omic annotaion violin plot / box plot"),
   "variable annotations (omic category marginals)",
   br(),br(),
   fluidRow(
     column(
       3, style = "border-right: 2px solid black",
       radioButtons(ns("RB_dist_plot_type2"), "plot type:",
                    choices = c("violin", "boxplot"),
                    selected = "violin", inline = TRUE),
       checkboxInput(ns("CB_show_data_points2"), "Show data points", value = FALSE),

       br()

     ), # End of column (6 space)
     column(9,
            uiOutput(ns("UI_var_box_output")#,
                     # downloadButton("Van_c1oup.pdf", "Download pDF"),
                     # downloadButton("Van_c1oup.png", "Download pNG"), br(),
                     # div(style="display:inline-block",
                     #     numericInput("Van_c1oup.h", "pDF / pNG height:", width = "138px",
                     #                  min = 4, max = 20, value = 8, step = 0.5)),
                     # div(style="display:inline-block",
                     #     numericInput("Van_c1oup.w", "pDF / pNG width:", width = "138px",
                     #                  min = 4, max = 20, value = 10, step = 0.5))
            )  # End of column (6 space)
     )    # End of fluidRow (4 space)

   )



  )  #taglist
}

#' pg_vis_raw Server Functions
#'
#' @param id
#' @param rv_data
#' @param rv_selections
#' @param heat_data
#' @param box_data
#' @param varbox_data
#'
#' @noRd
mod_pg_vis_raw_server <- function(id, rv_data, rv_selections, heat_data, box_data, varbox_data){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

   # TODO: enable interactive setting of these values
    data_point_sz = .65
    plot_size = "Small"
    font_size = "Small"
    color_schemes = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-purple")
    color_scheme = color_schemes[2]

    #go_signal <- reactive( rv_selections$omics_list$viz )
    output$plot_box_out <- renderPlot({
      #req(rv_data$ad)
      req( rv_selections$omics_list,
           box_data$data)

      plot_type <- input$RB_dist_plot_type
      show_data_points <- input$CB_show_data_points

      notch <- TRUE
      data_point_sz = .65
      plot_size = "Small"
      font_size = "Small"
      grp <- ifelse(rv_selections$group_action=="none",FALSE,TRUE)

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
      # color_schemes = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-purple")
      # color_scheme = color_schemes[2]
      # in_quant <- dat_key #(maybe) just observ_y
      # pg_violin_box(in_conf, in_meta, in_fact, in_quant,
      #               in_grp, in_subset, in_subsel,
      #               in_data, in_omic, plot_type,
      #               show_data_points, data_point_sz, font_size )
      return(vplot)
    })


    output$plot_varbox_out <- renderPlot({
      #req(rv_data$ad)
      req(rv_selections$omics_list,
          varbox_data$data)

      plot_type <- input$RB_dist_plot_type2
      show_data_points <- input$CB_show_data_points2

      notch <- TRUE
      data_point_sz = .65
      plot_size = "Small"
      font_size = "Small"
      grp <- ifelse(rv_selections$group_action=="none",FALSE,TRUE)

      vplot <- violin_box(varbox_data$data, varbox_data$x_name, varbox_data$y_name,
                          varbox_data$colors,
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
      # color_schemes = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-purple")
      # color_scheme = color_schemes[2]
      # in_quant <- dat_key #(maybe) just observ_y
      # pg_violin_box(in_conf, in_meta, in_fact, in_quant,
      #               in_grp, in_subset, in_subsel,
      #               in_data, in_omic, plot_type,
      #               show_data_points, data_point_sz, font_size )
      return(vplot)
    })

    output$UI_box_output <- renderUI({
      #if (rv_selections$omics_list$viz_now) {
        plotOutput(ns("plot_box_out"), height = pList2[plot_size])
      #}
    })

    output$UI_var_box_output <- renderUI({
      #if (rv_selections$omics_list$viz_now) {
      plotOutput(ns("plot_varbox_out"), height = pList2[plot_size])
      #}
    })


    # plot_heatmap_out renderPlot---------------------------------
    output$plot_heatmap_out <- renderPlot({
      req(heat_data$data)

      input_data <- heat_data$data
      #
      x_names <- heat_data$x_names
      y_names <- heat_data$y_names

      plot_type <- input$RB_heat_plot_type
      in_do_scale <- input$CB_scale
      in_clust_row <- input$CB_cluster_rows
      in_clust_col <- input$CB_cluster_cols

      plot_size = "Small"
      color_schemes = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-purple")
      color_scheme = color_schemes[2]

      grp <- ifelse(rv_selections$group_action=="none",FALSE,TRUE)

      hmap <- bubble_heatmap(input_data, x_names, y_names, plot_type,
                                 in_do_scale, in_clust_row, in_clust_col,
                                 color_scheme, plot_size, grp, save = FALSE)
#


    return(hmap)
    })

    output$UI_heatmap <- renderUI({
      #if (rv_selections$omics_list$viz_now) {
        plotOutput(ns("plot_heatmap_out"), height = pList3[plot_size])

    })


    output$HTML_header <- renderUI({

      omic_list = rv_selections$omics_list$value
      if(nrow(omic_list) > 50){
        HTML("More than 50 input omics! please reduce the omic list!")
      } else {
        oup = paste0(nrow(omic_list[present == TRUE]), " omic OK and will be plotted")
        if(nrow(omic_list[present == FALSE]) > 0){
          oup = paste0(oup, "<br/>",
                       nrow(omic_list[present == FALSE]), " omic not found (",
                       paste0(omic_list[present == FALSE]$omic, collapse = ", "), ")")
        }
        HTML(oup)
      }
    })
    # TODO: enable save/export
    # output$Van_d1oup.pdf <- downloadHandler(
    #   filename = function() { paste0("Van_",input$Van_d1plt,"_",input$Van_d1grp,".pdf") },
    #   content = function(file) { ggsave(
    #     file, device = "pdf", height = input$Van_d1oup.h, width = input$Van_d1oup.w,
    #     plot = bubble_heatmap(Van_conf, Van_meta, input$Van_d1inp, input$Van_d1grp, input$Van_d1plt,
    #                       input$Van_d1sub1, input$Van_d1sub2, "Van_gexpr.h5", Van_gene,
    #                       input$Van_d1scl, input$Van_d1row, input$Van_d1col,
    #                       input$Van_d1cols, input$Van_d1fsz, save = TRUE) )
    #   })
    # output$Van_d1oup.png <- downloadHandler(
    #   filename = function() { paste0("Van_",input$Van_d1plt,"_",input$Van_d1grp,".png") },
    #   content = function(file) { ggsave(
    #     file, device = "png", height = input$Van_d1oup.h, width = input$Van_d1oup.w,
    #     plot = bubble_heatmap(Van_conf, Van_meta, input$Van_d1inp, input$Van_d1grp, input$Van_d1plt,
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

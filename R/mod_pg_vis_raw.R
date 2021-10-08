
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
  tag_list <- tagList(
   HTML("Bubbleplot / Heatmap"),
   h4("Gene expression bubbleplot / heatmap"),
   "In this tab, users can visualise the expression patterns of ",
   "multiple omics grouped by categorical cell information (e.g. library / cluster).", br(),
   "The normalised expression values are group averaged and scaled/thresholded (?)).",
   br(),br(),
   fluidRow(
     hr(style = "border-top: 1px solid #000000;"),
     column(
       3, style="border-right: 2px solid black",
       radioButtons(ns("RB_heat_plot_type"), "plot type:",
                    choices = c("Bubbleplot", "Heatmap"),
                    selected = "Heatmap", inline = TRUE),
       checkboxInput(ns("CB_scale"), "Scale omic expression", value = TRUE),
       br()

     ), # End of column (6 space)
     column(
       3, style="border-right: 2px solid black",
       checkboxInput(ns("CB_cluster_rows"), "Cluster rows (omics)", value = TRUE),
       checkboxInput(ns("CB_cluster_cols"), "Cluster columns (samples)", value = FALSE),
       br()
     ),
     column(
       3,
       checkboxInput(ns("CB_show_grp_rows"), "Show row groupings (omics)", value = FALSE),
       checkboxInput(ns("CB_show_grp_cols"), "Show col groupings (samples)", value = FALSE),
       checkboxInput(ns("CB_group_agg"), "aggregate groupings? (samples)", value = FALSE),

       br()
     )
   ),
   fluidRow(
     hr(style = "border-top: 1px dashed grey;"),
            h4(htmlOutput("HTML_header")),
            uiOutput(ns("UI_heatmap_agg"))
            # downloadButton("Van_d1oup.pdf", "Download pDF"),
            # downloadButton("Van_d1oup.png", "Download pNG"), br(),
            # div(style="display:inline-block",
            #     numericInput("Van_d1oup.h", "pDF / pNG height:", width = "138px",
            #                  min = 4, max = 20, value = 10, step = 0.5)),
            # div(style="display:inline-block",
            #     numericInput("Van_d1oup.w", "pDF / pNG width:", width = "138px",
            #                  min = 4, max = 20, value = 10, step = 0.5))
     #)  # End of column (6 space)
   ),    # End of fluidRow (4 space)
   fluidRow(
     hr(style = "border-top: 1px dashed grey;"),
     uiOutput(ns("UI_heatmap_all"))
     # downloadButton("Van_d1oup.pdf", "Download pDF"),
     # downloadButton("Van_d1oup.png", "Download pNG"), br(),
     # div(style="display:inline-block",
     #     numericInput("Van_d1oup.h", "pDF / pNG height:", width = "138px",
     #                  min = 4, max = 20, value = 10, step = 0.5)),
     # div(style="display:inline-block",
     #     numericInput("Van_d1oup.w", "pDF / pNG width:", width = "138px",
     #                  min = 4, max = 20, value = 10, step = 0.5))
     #)  # End of column (6 space)
   )    # End of fluidRow (4 space)
   )

  return(tag_list)
}

#' pg_vis_raw Server Functions
#'
#' @param id shiny internal
#' @param rv_data reactives from ingest
#' @param rv_selections reactives from side selector
#' @param heat_data heatmap data
#' @param box_data boxplot data
#' @param varbox_data var boxplot data
#'
#' @noRd
mod_pg_vis_raw_server <- function(id, rv_data, rv_selections, heat_data){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

   # TODO: enable interactive setting of these values
    data_point_sz = .65
    plot_size = "Small"
    font_size = "Small"
    color_schemes = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-purple")
    color_scheme = color_schemes[2]



    # plot_heatmap_out renderPlot---------------------------------
    output$plot_heatmap_agg_out <- renderPlot({
      req(heat_data$data)

      input_data <- heat_data$data
      #
      x_names <- heat_data$x_names
      y_names <- heat_data$y_names

      plot_type <- input$RB_heat_plot_type
      in_do_scale <- input$CB_scale
      in_clust_row <- input$CB_cluster_rows
      in_clust_col <- input$CB_cluster_cols

      in_agg_grp <- input$CB_group_agg

      plot_size = "Small"
      color_schemes = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-purple")
      color_scheme = color_schemes[2]

      grp <- ifelse(rv_selections$group_action=="none",FALSE,TRUE)
      in_hdata <- isolate(input_data)

      # Aggregate:  exp(mean) + proportion -> log(val)
      # disabling expm1 and log1p because we will be putting "NORMALIZED" quantities in...
      if (max(abs(in_hdata$val),na.rm=TRUE) < 700.0 ){
        fwd_scale <- expm1
        inv_scale <- log1p
      } else {
        fwd_scale <- function( input ){ input }
        inv_scale <- function( input ){ input }
      }

      # in_hdata:
      #   X_ID
      #   sub
      #   Y_nm
      #   value

      unit_name = "scaled 'X'"
      # remove NA
      in_hdata <- in_hdata[!is.na(val)]
      in_hdata$val = fwd_scale(in_hdata$val)

      # AGGREGATE
      in_hdata = in_hdata[, .(val = mean(val), prop = sum(val>0) / length(X_ID)),
                      by = c("Y_nm", "X_nm")]
      # in_hdata = in_hdata[, .(val = mean(val,na.rm=TRUE), prop = sum(val>0,na.rm=TRUE) / length(sample_ID)),
      #                 by = c("omic", "grp_by")]

      in_hdata$val = inv_scale(in_hdata$val) # do we need this??? we are already normalized...

      # remove zeros

      # Scale if required
      if(in_do_scale){
        in_hdata[, val:= scale(val), keyby = "Y_nm"]
      }

      # hclust row/col if necessary
      gg_mat = dcast.data.table(in_hdata, Y_nm~X_nm, value.var = "val")
      tmp = gg_mat$Y_nm
      gg_mat = as.matrix(gg_mat[, -1])
      rownames(gg_mat) = tmp

      ht <- ComplexHeatmap::Heatmap(gg_mat,
                                        cluster_rows = in_clust_row,
                                        cluster_columns = in_clust_col,
                                        #column_split = grp_x,
                                        #row_split = grp_y,
                                        name = unit_name,
                                        column_title = "factor",
                                        row_title = "omics",
                                        show_parent_dend_line = FALSE)

      return(ht)
    })



    # plot_heatmap_out renderPlot---------------------------------
    output$plot_heatmap_all_out <- renderPlot({
      req(heat_data$data)

      input_data <- heat_data$data
      #
      x_names <- heat_data$x_names
      y_names <- heat_data$y_names

      plot_type <- input$RB_heat_plot_type
      in_do_scale <- input$CB_scale
      in_clust_row <- input$CB_cluster_rows
      in_clust_col <- input$CB_cluster_cols

      in_agg_grp <- input$CB_group_agg
      unit_name = "scaled 'X'"


      plot_size = "Small"
      color_schemes = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-purple")
      color_scheme = color_schemes[2]

      grp <- ifelse(rv_selections$group_action=="none",FALSE,TRUE)

      # hmap <- bubble_heatmap(input_data, x_names, y_names, plot_type,
      #                            in_do_scale, in_clust_row, in_clust_col,
      #                            color_scheme, plot_size, grp, save = FALSE)

      # hmap <- simple_bubble_heatmap(input_data, x_names, y_names, plot_type,
      #                             in_do_scale, in_clust_row, in_clust_col,
      #                             color_scheme, plot_size, grp, grp,save = FALSE)
      #
      #
      # simple_bubble_heatmap <- function(in_hdata, x_names, y_names, plot_type,
      #                                   in_do_scale, in_clust_row, in_clust_col,
      #                                   color_scheme, plot_size, grp_x, grp_y, save = FALSE)

      in_data <- isolate(input_data)

      # Scale if required
      if(in_do_scale){
        in_data[, val:= scale(val), keyby = "Y_nm"]
      }

      rawgg_mat = dcast.data.table(in_data, Y_nm~X_ID, value.var = "val")
      tmp = rawgg_mat$Y_nm
      rawgg_mat = as.matrix(rawgg_mat[, -1])
      rownames(rawgg_mat) = tmp

      grp_x = unique(in_data$X_ID)
      tmp <- input_data[group_y==in_data$group_y[1]]
      grp_x <- tmp$group_x


      #ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:6), labels = LETTERS[1:5]))
      ha = ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_block(labels = levels(grp_x)))
      # need to fix the raster here...
      # suppress the "colnames" when we are plotting all the samples... all the genes

      ht <- ComplexHeatmap::Heatmap(rawgg_mat,
                                      cluster_rows = in_clust_row,
                                      cluster_columns = in_clust_col,
                                      column_split = grp_x,
                                      #row_split = grp_y,
                                      top_annotation = ha,
                                      name = unit_name,
                                      column_title = "factor",
                                      row_title = "omics",
                                      show_parent_dend_line = FALSE,
                                      use_raster = FALSE)
      #top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:4))),


      return(ht)
    })


    output$UI_heatmap_agg <- renderUI({
      #if (rv_selections$omics_list$viz_now) {
        plotOutput(ns("plot_heatmap_agg_out"), height = pList3[plot_size])

    })

    output$UI_heatmap_all <- renderUI({
      #if (rv_selections$omics_list$viz_now) {
      plotOutput(ns("plot_heatmap_all_out"), height = pList3[plot_size])

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

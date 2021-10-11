#TODO:  change name of file and module to viz_expression or something
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
       hr(style = "border-top: 1px dashed grey;"),
       checkboxInput(ns("CB_scale"), "Scale values?", value = TRUE),
       br()

     ), # End of column (6 space)
     column(
       3, style="border-right: 2px solid black",
       checkboxInput(ns("CB_target_omics"), "target -omics only?", value = TRUE),
       br(),br()
     ),
     column(
       3, style="border-right: 2px solid black",
       checkboxInput(ns("CB_cluster_rows"), "Cluster rows (omics)", value = FALSE),
       checkboxInput(ns("CB_cluster_cols"), "Cluster columns (samples)", value = TRUE),
       br()
     ),
     column(
       3,
       checkboxInput(ns("CB_show_grp_rows"), "Show row groupings (omics)", value = FALSE),
       checkboxInput(ns("CB_show_grp_cols"), "Show col groupings (samples)", value = FALSE),
       #checkboxInput(ns("CB_group_agg"), "aggregate groupings? (samples)", value = FALSE),
       br()
     )
   ),
   fluidRow(
     hr(style = "border-top: 1px dashed grey;"),
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
     h4(htmlOutput(ns("HTML_header"))),
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

      #TODO:
      # -limit levels of clustering columns
      # -
      #

      input_data <- heat_data$data
      #
      plot_type <- input$RB_heat_plot_type
      in_clust_row <- input$CB_cluster_rows
      in_clust_col <- input$CB_cluster_cols

      in_agg_grp <- input$CB_group_agg

      plot_size = "Small"
      color_schemes = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-purple")
      color_scheme = color_schemes[2]

      dmat <- t(heat_data$mat)


      #grouping <- ifelse(rv_selections$group_action=="none",FALSE,TRUE)
      in_hdata <- isolate(input_data)

      # Aggregate:  exp(mean) + proportion -> log(val)
      # disabling expm1 and log1p because we will be putting "NORMALIZED" quantities in...
     test_v <- max( abs(in_hdata$val) , na.rm=TRUE)
      if ( test_v < 700.0 ){
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
      if(input$CB_scale){
        in_hdata[, val:= scale(val), keyby = "Y_nm"]
      }

      # hclust row/col if necessary
      gg_mat = dcast.data.table(in_hdata, Y_nm~X_nm, value.var = "val")

      gg_mat2 = dcast.data.table(in_hdata, Y_nm~X_nm, value.var = "prop")

      tmp = gg_mat$Y_nm
      gg_mat = as.matrix(gg_mat[, -1])
      rownames(gg_mat) = tmp

      if (dim(gg_mat)[2] < 3){
        in_clust_col <- FALSE
      }

      if (input$CB_target_omics){
        omics <- heat_data$selected_omics$target_omics
        in_hdata <- in_hdata[in_hdata$Y_nm %in% omics,]
        #gg_mat2 <- gg_mat[rownames(gg_mat) %in% omics,]
        gg_mat <- gg_mat[omics,]
        #filter
      }

      if (plot_type=="Bubbleplot"){
        grp <- heat_data$x_group

        ht <- bubble_heatmap2(in_hdata, gg_mat, plot_type,
                                   in_do_scale, in_clust_row, in_clust_col,
                                   color_scheme, plot_size, grp, save = FALSE)

      } else {
            ht <- ComplexHeatmap::Heatmap(gg_mat,
                                              cluster_rows = in_clust_row,
                                              cluster_columns = in_clust_col,
                                              #column_split = grp_x,
                                              #row_split = grp_y,
                                              name = unit_name,
                                              column_title = "factor",
                                              row_title = "omics",
                                              show_parent_dend_line = FALSE)

    }

      return(ht)

    })



    # plot_heatmap_out renderPlot---------------------------------
    output$plot_heatmap_all_out <- renderPlot({
      req(heat_data$data)

      unit_name = "scaled 'X'"

      plot_type <- input$RB_heat_plot_type
      in_do_scale <- input$CB_scale
      in_clust_row <- FALSE # input$CB_cluster_rows
      in_clust_col <- input$CB_cluster_cols

      in_agg_grp <- input$CB_group_agg

      plot_size = "Small"
      color_schemes = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-purple")
      color_scheme = color_schemes[2]

      dmat <- t(heat_data$mat)

      in_data <- isolate(heat_data$data)

      # Scale if required
      if (in_do_scale) {
        in_data[, val:= scale(val), keyby = "Y_nm"]
      }

      rawgg_mat = dcast.data.table(in_data, Y_nm~X_ID, value.var = "val")
      tmp = rawgg_mat$Y_nm
      rawgg_mat = as.matrix(rawgg_mat[, -1])
      rownames(rawgg_mat) = tmp

      grp_x = unique(in_data$X_ID)
      #tmp <- input_data[group_y==in_data$group_y[1]]

      grp_x <- heat_data$meta[[heat_data$x_group]]

      target_omics <- heat_data$selected_omics$target_omics
      omics_at <- which(rownames(dmat) %in% target_omics)
      #ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:6), labels = LETTERS[1:5]))
      ha = ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_block(labels = levels(grp_x)))

      ha2 = ComplexHeatmap::rowAnnotation(foo = ComplexHeatmap::anno_mark(at = omics_at,
                                          labels = target_omics))


      # need to fix the raster here...
      # suppress the "colnames" when we are plotting all the samples... all the genes
      #TODO:
      # -limit levels of clustering columns
      # -
      ht <- ComplexHeatmap::Heatmap(dmat,
                                      cluster_rows = in_clust_row,
                                      cluster_columns = in_clust_col,
                                      column_split = grp_x,
                                      #row_split = grp_y,
                                      top_annotation = ha,
                                      right_annotation = ha2,
                                      row_names_side = "left",
                                      row_names_gp = grid::gpar(fontsize = 4),
                                      name = unit_name,
                                      column_title = "factor",
                                      row_title = "omics",
                                      show_parent_dend_line = FALSE,
                                      use_raster = FALSE)
      #top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:4))),


      return(ht)
    })


    output$UI_heatmap_agg <- renderUI({
      #if (rv_selections$target_omics$viz_now) {
      req(heat_data$data)
      out_text <- paste("hm dim: ", paste(dim(heat_data$data),collapse = " x "))

      to_return <-  tagList(
        h4("grouped average  expression"),
        "(Aggregated)",
        out_text,
        br(),br(),
        plotOutput(ns("plot_heatmap_agg_out"), height = "1200px")
      )

      return(to_return)


    })# Panel sizes
    # pList = c("400px", "600px", "800px")
    # names(pList) = c("Small", "Medium", "Large")
    # pList2 = c("500px", "700px", "900px")
    # names(pList2) = c("Small", "Medium", "Large")
    # pList3 = c("600px", "800px", "1000px")
    # names(pList3) = c("Small", "Medium", "Large")
    # sList = c(18,24,30)
    # names(sList) = c("Small", "Medium", "Large")
    # lList = c(5,6,7)
    # names(lList) = c("Small", "Medium", "Large")
    output$UI_heatmap_all <- renderUI({
      #if (rv_selections$target_omics$viz_now) {
      req(heat_data$data)

      to_return <-  tagList(
        h4("full heatmap"),
        #"some more text",
        br(),br(),
        plotOutput(ns("plot_heatmap_all_out"), height = "1200px")
      )

    })

    output$HTML_header <- renderUI({
      omic_list = rv_selections$target_omics$value
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



  })
}

## To be copied in the UI
# mod_pg_vis_raw_ui("pg_vis_raw_ui_1")

## To be copied in the server
# mod_pg_vis_raw_server("pg_vis_raw_ui_1")

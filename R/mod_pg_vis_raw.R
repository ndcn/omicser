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
     hr(style = "border-top: 1px dashed grey;"),
     column(
       2, style="border-right: 2px solid black",

       checkboxInput(ns("CB_aggregate_conditions"), "aggregate conditions?", value = FALSE),
       checkboxInput(ns("CB_scale"), "Scale values?", value = TRUE),

       hr(style = "border-top: 1px dashed grey;"),
       checkboxInput(ns("CB_cluster_rows"), "Cluster rows (omics)", value = TRUE),
       checkboxInput(ns("CB_cluster_cols"), "Cluster columns (samples)", value = TRUE),

       br(),
       checkboxInput(ns("CB_show_grp_rows"), "Show row groupings (omics)", value = FALSE),
       checkboxInput(ns("CB_show_grp_cols"), "Show col groupings (samples)", value = FALSE),

       # Input: Decimal interval with step value ----
       hr(style = "border-top: 1px dashed grey;"),
       checkboxInput(ns("CB_show_all_feats"), "show all features?", value = TRUE),
       "warning! this will be slow",

       shinyjs::disabled(
         radioButtons(ns("RB_heat_plot_type"), "plot type:",
                      choices = c("Bubbleplot", "Heatmap"),
                      selected = "Heatmap", inline = TRUE)
       )
       # this selector for aggregated heatmap
       # A
       #uiOutput(outputId = ns("UI_comp_group_selection"))

     ), # End of column (6 space)
     column(10,
            uiOutput(ns("UI_heatmap_header")),
            plotOutput(ns("plot_heatmap_out"), height = "1200px")

            # downloadButton("Van_d1oup.pdf", "Download pDF"),
            # downloadButton("Van_d1oup.png", "Download pNG"), br(),
            # div(style="display:inline-block",
            #     numericInput("Van_d1oup.h", "pDF / pNG height:", width = "138px",
            #                  min = 4, max = 20, value = 10, step = 0.5)),
            # div(style="display:inline-block",
            #     numericInput("Van_d1oup.w", "pDF / pNG width:", width = "138px",
            #                  min = 4, max = 20, value = 10, step = 0.5))
     #)  # End of column (6 space)
     )
   # ),    # End of fluidRow (4 space)
   # fluidRow(
   #   hr(style = "border-top: 1px dashed grey;"),
   #   h4(htmlOutput(ns("HTML_header"))),
   #   uiOutput(ns("UI_heatmap_all"))
   #   # downloadButton("Van_d1oup.pdf", "Download pDF"),
   #   # downloadButton("Van_d1oup.png", "Download pNG"), br(),
   #   # div(style="display:inline-block",
   #   #     numericInput("Van_d1oup.h", "pDF / pNG height:", width = "138px",
   #   #                  min = 4, max = 20, value = 10, step = 0.5)),
   #   # div(style="display:inline-block",
   #   #     numericInput("Van_d1oup.w", "pDF / pNG width:", width = "138px",
   #   #                  min = 4, max = 20, value = 10, step = 0.5))
   #   #)  # End of column (6 space)
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
mod_pg_vis_raw_server <- function(id, rv_data, rv_selections, heat_data){ #}, agg_heat){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

   # TODO: enable interactive setting of these values
    data_point_sz = .65
    plot_size = "Small"
    font_size = "Small"
    color_schemes = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-purple")
    color_scheme = color_schemes[2]

    wtr <- waiter::Waiter$new(
      id = ns("plot_heatmap_out"),
      html = waiter::spin_loaders(id=12,color="red"),
      color = waiter::transparent(.5)  #waiter::transparent(.5)
    )

    observeEvent(input$plot_heatmap_out_waiter_hidden, {
      message(paste0("hide: ",input$plot_heatmap_out_waiter_hidden))
    })

    observeEvent(input$plot_heatmap_out_waiter_shown, {
      message(paste0("hide: ",input$plot_heatmap_out_waiter_shown))
    })

    observe({
      req(input$SLD_nfeats)
      if ( input$CB_show_all_feats ) {
        shinyjs::enable("SLD_nfeats")
      } else {
        shinyjs::disable("SLD_nfeats")
      }
    })



    output$UI_heatmap_header <- renderUI({
      #if (rv_selections$target_omics$viz_now) {
      req(heat_data$mat,
          heat_data$selected_omics$target_omics)

      out_text1 <- paste("hm dim: ",
                         paste(dim(heat_data$mat),collapse = " x "))
      out_text2 <-  paste("n target features: ",
                          length(heat_data$selected_omics$target_omics))

      to_return <-  tagList(
        h4("*-omic expression"),
        out_text1,
        ",", out_text2,
        br(),
      )
      return(to_return)
    })



    # plot_heatmap_out renderPlot---------------------------------
    output$plot_heatmap_out <- renderPlot( {
      req(heat_data$mat)
      message("starting: plot_heatmap_out packer")

      #TODO:  turn this list of "annotations" into heathmap annotations...
      wtr$show()

      in_mat <- heat_data$mat  #imported as samples X features

      # Scale if required
      units_label <- 'omic\nexpr'
      if(input$CB_scale) {
        in_mat <- scale(in_mat) # scale operates on columns... and for now we have omics in the columns
        units_label <-  'omic\nexpr\n(Z-score)'
        if ( dim(in_mat)[1]<5)   {
          message("warning!!! small number of groups... could be indeterminate...")
        }
      }


      # need to figure out the types of annotations here.?..
      top_annotations <- heat_data$x_annot

      right_annotations <- heat_data$y_annot


      # aggregate if required
      if( input$CB_aggregate_conditions ) {
        in_hdata <- as.data.table(in_mat)
        in_hdata$grp <- as.character(heat_data$x_names)

        agg_hdata <-   in_hdata[, lapply(.SD, mean), by = grp]
        tmp <- agg_hdata$grp
        agg_mat <- as.matrix(agg_hdata[,grp:=NULL])
        rownames(agg_mat) <- tmp

        in_mat <- agg_mat

        x_aggregated <- TRUE
        x_grp <- heat_data$x_names # NULL  #TESTING

        x_title <- paste0("samples X",heat_data$x_group)

        # what to do with top_annotations when grouping?  distributions of each variable?
      } else {
        x_aggregated <- FALSE
        x_grp <- heat_data$x_names

        x_title <- "samples" #heat_data$x_group




      }

      in_mat <- t(in_mat)
      #  depricated bubbleplot
      # plot_type <- input$RB_heat_plot_type

      omics <- heat_data$selected_omics$target_omics

      # if we are trying to show everything we need to wubset...
      if (input$CB_show_all_feats) {

        omics_at <- which(rownames(in_mat) %in% omics)
        show_row_names <- FALSE #replace with annotation
        omics_title <- "all (sampled) features"

        # if (input$CB_cluster_rows){
        #   cluster_rows <- cluster_within_group(t(in_mat), mygroup[samp_omics_at])
        # } else {
        #   cluster_rows <- FALSE
        # }

      } else {
        in_mat <- in_mat[omics,]
        show_row_names = TRUE
        ha2 <- NULL
        omics_at <- which(rownames(in_mat) %in% omics)

        # if (input$CB_cluster_rows){
        #   cluster_rows <- hclust(dist(in_mat))
        # } else {
        #   cluster_rows <- FALSE
        # }

        show_row_names <- TRUE
        omics_title <- "target features"

        # subset side annotations...?
        right_annotations <- right_annotations[omics_at,]
      }

      cluster_rows <- if (input$CB_cluster_rows) hclust(dist(in_mat)) else FALSE


      cluster_columns <- input$CB_cluster_cols
      if (dim(in_mat)[2] < 3){
        cluster_columns <- FALSE
      }

      row_split <- NULL
      show_column_names <- TRUE

      hm <- make_cx_heatmap(in_mat,
                          cluster_rows, row_split, show_row_names, omics_title,
                          cluster_columns, column_split,show_column_names, x_aggregated, x_title,x_grp,x_names,
                          units_label, omics, omics_at, top_annotations, right_annotations)


      return(hm)

    })




#
#
#       # zero NAs?
#       #
#       # Log? (assume all the counts tables are normalized prior to aggregation...)
#       #
#       # Scale if required
#       units_label <- 'omic\nexpr'
#
#       if ( input$CB_scale ) {
#         in_mat <- scale(in_mat) # scale operates on columns... and for now we have omics in the columns
#         units_label <-  'omic\nexpr\n(Z-score)'
#       }
#
#
#       # START FIX HERE
#       in_mat <- t(in_mat)
#
#       x_title <- "samples" #heat_data$x_group
#
#       grp_x <- as.character(heat_data$x_names) #make sure its not a factor
#       grp_x <- heat_data$x_names #make sure its not a factor
#
#       omics <- heat_data$selected_omics$target_omics
#       omics_at <- which(rownames(in_mat) %in% omics)
#
#       omics_title <- "*-omics"
#
#       if (input$CB_cluster_cols) {
#         ha <- NULL
#         grp_x <- NULL
#         #ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:6), labels = LETTERS[1:5]))
#       } else {
#         ha <- ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_block(labels = levels(grp_x)))
#       }
#
#       ha2 <- ComplexHeatmap::rowAnnotation(foo = ComplexHeatmap::anno_mark(at = omics_at,
#                                           labels = omics))
#       show_row_names = FALSE
#
#
#       # need to fix the raster here...
#       # suppress the "colnames" when we are plotting all the samples... all the genes
#       #TODO:
#       # -limit levels of clustering columns
#       # -
#       message("plot_heatmap_all_out: CALLING ComplexHeatmap::Heatmap")
#       ht <- ComplexHeatmap::Heatmap(in_mat,
#                                     cluster_rows = in_clust_row,
#                                     cluster_columns = in_clust_col,
#                                     column_split = grp_x,
#                                     #row_split = grp_y,
#                                     top_annotation = ha,
#                                     show_row_names = show_row_names,
#                                     right_annotation = ha2,
#                                     row_names_side = "right",
#                                     #row_names_side = "left",
#                                     row_names_gp = grid::gpar(fontsize = 7),
#                                     name = units_label,
#                                     column_title = x_title,
#                                     row_title = omics_title,
#                                     show_parent_dend_line = TRUE,
#                                     use_raster = FALSE) #,
#                                     #raster_device = "png")
#       message("plot_heatmap_all_out: FINISHED ComplexHeatmap::Heatmap")

      #top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:4))),
#
#       # parameters for the colour-bar that represents gradient of expression
#       heatmap_legend_param = list(
#         color_bar = 'continuous',
#         legend_direction = 'vertical',
#         legend_width = unit(8, 'cm'),
#         legend_height = unit(5.0, 'cm'),
#         title_position = 'topcenter',
#         title_gp=gpar(fontsize = 12, fontface = 'bold'),
#         labels_gp=gpar(fontsize = 12, fontface = 'bold')),

      # # row (gene) parameters
      # cluster_rows = TRUE,
      # show_row_dend = TRUE,
      # #row_title = 'Statistically significant genes',
      # row_title_side = 'left',
      # row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
      # row_title_rot = 90,
      # show_row_names = FALSE,
      # row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
      # row_names_side = 'left',
      # row_dend_width = unit(25,'mm'),

      # # column (sample) parameters
      # cluster_columns = TRUE,
      # show_column_dend = TRUE,
      # column_title = '',
      # column_title_side = 'bottom',
      # column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
      # column_title_rot = 0,
      # show_column_names = FALSE,
      # column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
      # column_names_max_height = unit(10, 'cm'),
      # column_dend_height = unit(25,'mm'),
#
#       return(ht)
#     })
# Panel sizes
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
      req(heat_data$mat)
      out_text <- paste("hm dim: ", paste(dim(heat_data$mat),collapse = " x "))


      to_return <-  tagList(
        h4("full heatmap"),
        #"some more text",
        out_text,
        br(),br(),
        plotOutput(outputId = ns("plot_heatmap_all_out"),
                   height = "1200px")
      )


      return(to_return)
    })

    output$HTML_header <- renderUI({
      omic_list <- heat_data$selected_omics$all_omics
      if( length(omic_list) > 100){
        HTML("More than 100 input omics! please reduce the omic list!")
      } else {
        HTML(paste0(length(omic_list), " omics OK and will be plotted"))
      }
    })



  })
}

## To be copied in the UI
# mod_pg_vis_raw_ui("pg_vis_raw_ui_1")

## To be copied in the server
# mod_pg_vis_raw_server("pg_vis_raw_ui_1")

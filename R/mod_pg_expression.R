#TODO:  change name of file and module to viz_expression or something
#' pg_vis_raw UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList downloadButton
mod_pg_expression_ui <- function(id){
  ns <- NS(id)
  tag_list <- tagList(
   HTML("Annotated Heatmap"),
   h4("-omic expression heatmap"),
   "In this tab, users can visualise the expression patterns of the data.",
   #"Multiple -omic features grouped by categorical cell information (e.g. library / cluster).", br(),
   #"The normalised expression values can be aggregated over conditions.",
   br(),br(),

   fluidRow(
     hr(style = "border-top: 1px dashed grey;"),
     column(
       2, style="border-right: 2px solid black",

       checkboxInput(ns("CB_aggregate_conditions"), "aggregate conditions?", value = TRUE),
       checkboxInput(ns("CB_scale"), "Scale values?", value = TRUE),

       hr(style = "border-top: 1px dashed grey;"),
       checkboxInput(ns("CB_cluster_rows"), "Cluster features (rows)", value = TRUE),
       checkboxInput(ns("CB_cluster_cols"), "Cluster samples (columns)", value = TRUE),

       br(),
       checkboxInput(ns("CB_show_grp_rows"), "Show feature groupings (rows)", value = FALSE),
       checkboxInput(ns("CB_show_grp_cols"), "Show sample groupings (columns)", value = FALSE),

       # Input: Decimal interval with step value ----
       hr(style = "border-top: 1px dashed grey;"),
       ">DISABLED<",
       shinyjs::disabled(checkboxInput(ns("CB_show_all_feats"), "show all features?", value = TRUE)),
       "(coming soon)",

       #"warning! this will be slow",
       # shinyjs::disabled(
       #   radioButtons(ns("RB_heat_plot_type"), "plot type:",
       #                choices = c("Bubbleplot", "Heatmap"),
       #                selected = "Heatmap", inline = TRUE)
       # )
       # this selector for aggregated heatmap
       # A
       #uiOutput(outputId = ns("UI_comp_group_selection"))

       downloadButton(ns("DB_heatmap_pdf"), "Download PDF"),
       downloadButton(ns("DB_heatmap_png"), "Download PNG"), br(),
       div(style="display:inline-block",
           numericInput(ns("heatmap_h"), "PDF / PNG height:", width = "138px",
                        min = 4, max = 20, value = 10, step = 0.5)),
       div(style="display:inline-block",
           numericInput(ns("heatmap_w"), "PDF / PNG width:", width = "138px",
                        min = 4, max = 20, value = 10, step = 0.5))

     ), # End of column (6 space)
     column(10,
            uiOutput(ns("UI_heatmap_header")),
            plotOutput(ns("plot_heatmap_out"), height = "1200px")


     )  # End of column (6 space)

   )#fluid_row
   )#tag_list

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
mod_pg_expression_server <- function(id, rv_data, rv_selections, heat_data){ #}, agg_heat){
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
      html = waiter::spin_loaders(id=8,color="black"),
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

    heatmap_initialized = reactiveVal(FALSE)


    # REACTIVEVALUES =================================
    reactive_heatmap <- reactive({

      req(heat_data$mat)
      message("starting: plot_heatmap_out packer")

      #TODO:  turn this list of "annotations" into heathmap annotations...
      wtr$show()

      in_mat <- heat_data$mat  #imported as samples X features
      # Scale if required
      units_label <- 'expression'
      if(input$CB_scale) {
        in_mat <- scale(in_mat) # scale operates on columns... and for now we have omics in the columns
        units_label <-  'expression\n(Z-score)'
        if ( dim(in_mat)[1]<5)   {
          message("warning!!! small number of groups... could be indeterminate...")
        }
      }

      # need to figure out the types of annotations here.?..
      top_annotations <- heat_data$samp_annot

      # aggregate if required
      if( input$CB_aggregate_conditions ) {
        in_hdata <- as.data.table(in_mat)
        in_hdata$grp <- as.character(heat_data$samp_grp)

        agg_hdata <-   in_hdata[, lapply(.SD, mean), by = grp]
        tmp <- agg_hdata$grp
        agg_mat <- as.matrix(agg_hdata[,grp:=NULL])
        rownames(agg_mat) <- tmp

        in_mat <- agg_mat

        samp_aggregated <- TRUE
        samp_grp <- tmp # NULL  #TESTING
        samp_grp_nm <- heat_data$samp_grp_nm

        samp_title <- paste0(" avg. over `",samp_grp_nm,"`")


        # what to do with top_annotations when grouping?  distributions of each variable?
      } else {
        samp_aggregated <- FALSE
        samp_grp <- heat_data$samp_grp
        samp_grp_nm <- heat_data$samp_grp_nm

        samp_title <- "samples" #heat_data$samp_group

      }

      # makd sure that samp_grp_nm is in annotations...
      if (!any(colnames(top_annotations)==samp_grp_nm)) {
        top_annotations[[samp_grp_nm]] <- as.character(heat_data$samp_grp)
      }

      in_mat <- t(in_mat) #transpose so its features X samples

      right_annotations <- heat_data$feat_annot

      omics <- heat_data$selected_omics$target_omics

      feat_grp <- heat_data$feat_grp
      feat_grp_nm <- heat_data$feat_grp_nm

      cluster_feats <- input$CB_cluster_rows
      cluster_samps <- input$CB_cluster_cols

      samp_split <- input$CB_show_grp_cols
      # if we are trying to show everything we need to wubset...
      if (input$CB_show_all_feats) {
        omics_at <- which(rownames(in_mat) %in% omics)
        show_row_names <- FALSE #replace with annotation
        feats_title <- "all (sampled) features"

        # if (input$CB_cluster_rows){
        #   cluster_rows <- cluster_within_group(t(in_mat), mygroup[samp_omics_at])
        # } else {
        #   cluster_rows <- FALSE
        # }

      } else {
        in_mat <- in_mat[omics,]
        show_row_names = TRUE
        omics_at <- which(rownames(in_mat) %in% omics)

        # if (input$CB_cluster_rows){
        #   cluster_rows <- hclust(dist(in_mat))
        # } else {
        #   cluster_rows <- FALSE
        # }
        show_row_names <- TRUE
        feats_title <- "target features"

        # subset side annotations...?
        right_annotations <- right_annotations[omics_at,]
        feat_grp <- feat_grp[omics_at]

      }

      if (input$CB_show_grp_rows) {
        feat_grp <- feat_grp

      } else {
        feat_grp <- NULL
      }

      hm <- make_cx_heatmap(in_mat,
                            cluster_samps, samp_grp, samp_grp_nm, samp_title, samp_aggregated, samp_split,
                            cluster_feats, feat_grp, feat_grp_nm, feats_title,
                            units_label,
                            omics, omics_at,
                            top_annotations, right_annotations)

      # tryCatch({
      #   dev.null()
      #   if(inherits(hm, "Heatmap")) {
      #     hm = draw(hm + NULL)
      #     message("inherited hm + NULL")
      #   } else {
      #     if(!hm@layout$initialized) {
      #       hm = draw(hm)
      #       message("!initialized ")
      #     }
      #   }
      # }, finally = ComplexHeatmap:::dev.off2() )


      return( hm )

      })


    # plot_heatmap_out renderPlot---------------------------------
    output$plot_heatmap_out <- renderPlot( {
      hm <- reactive_heatmap()
      message(class(hm))

      ComplexHeatmap::draw(hm, merge_legend = TRUE)

    })



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


    # output$heatmap_download_button = downloadHandler(
    #
    #   filename = function() {
    #     format = as.numeric(input$heatmap_download_format)
    #     fm = c("png", "pdf", "svg")[format]
    #     paste0(heatmap_id, "_heatmap",fm)
    #   },
    #   content = function(file) {
    #     format = as.numeric(input$heatmap_download_format)
    #     fm = c("png", "pdf", "svg")[format]
    #     dev = list(png, pdf, svglite::svglite)[[format]]
    #
    #     showNotification(paste0("Download heatmap in ",fm), duration = 2, type = "message")
    #     message(paste0(Sys.time(), " Download heatmap in ",fm))
    #
    #     temp = tempfile()
    #     width = input$heatmap_w
    #     height = input$heatmap_h
    #
    #     if(fm == "png") {
    #       dev(temp, width = width*2, height = height*2, res = 72*2)
    #     } else if(fm == "pdf") {
    #       dev(temp, width = width/100*4/3, height = height/100*4/3)
    #     } else {
    #       dev(temp, width = width, height = height)
    #     }
    #     if(heatmap_initialized()) {
    #       draw(ht_list())
    #     } else {
    #       grid.newpage()
    #       grid.text("No heatmap is available.")
    #     }
    #     dev.off()
    #
    #     file.copy(temp, file)
    #   }
    # )


    output$DB_heatmap_pdf <- downloadHandler(
              filename = function() {
                #TODO: make the filename something sensible wrt database name...
                #  e.g.  db_name_layer_omics_type.pdf?
                fnm <- paste0(Sys.time(),"heatmap_",input$heatmap_h,"_",input$heatmap_w,".pdf")
                message(fnm)
                return( fnm )
               },
              content = function(file) {

                showNotification("Download heatmap in PDF", duration = 2, type = "message")
                message(paste0(Sys.time(), " Download heatmap as PDF") )
                temp = tempfile()
                width = input$heatmap_w
                height = input$heatmap_h

                pdf(temp, width = width*4/3, height = height*4/3)

                hm <- reactive_heatmap()
                print(class(hm))
                #draw(hm)
                ComplexHeatmap::draw(hm, merge_legend = TRUE)

                dev.off()

                file.copy(temp, file)
              })

    #
    # content = function(file) { ggplot2::ggsave(file,
    #                          device = "pdf",
    #                          height = input$heatmap_h,
    #                          width = input$heatmap_w,
    #                          plot = reactive_heatmap() )}
    #
    output$DB_heatmap_png <- downloadHandler(
                filename = function() {
                  fnm <- paste0(Sys.time(),"heatmap_",input$heatmap_h,"_",input$heatmap_w,".png")
                  message(fnm)
                  return(fnm)
                  },
                content = function(file) {
                  showNotification("Download heatmap as PNG", duration = 2, type = "message")
                  message(paste0(Sys.time(), " Download heatmap as PNG "))

                  hm <- reactive_heatmap()
                  print(class(hm))

                  temp = tempfile()
                  res <- 72*2
                  width = input$heatmap_w*res
                  height = input$heatmap_h*res

                  png(temp, width = width*4/3, height = height*4/3, res = res)

                  #draw(hm)
                  ComplexHeatmap::draw(hm, merge_legend = TRUE)
                  dev.off()

                  file.copy(temp, file)
                  }
            )
   }
)}


## To be copied in the UI
# mod_pg_expression_ui("pg_vis_raw_ui_1")

## To be copied in the server
# mod_pg_expression_server("pg_vis_raw_ui_1")

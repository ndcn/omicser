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
        # radioButtons(ns("RB_plot_type"), "Plot type:",
        #              choices = c("volcano", "other(TBD"),
        #              selected = "volcano", inline = TRUE),
        #checkboxInput(ns("CB_show_data_points"), "Show data points", value = FALSE),

         # actionButton("Van_c1togL", "Toggle to subset cells"),
        #
        # conditionalPanel(
        #   condition = "input.Van_c1togL % 2 == 1",
        #   selectInput("Van_c1sub1", "Cell information to subset:",
        #               choices = Van_conf[grp == TRUE]$UI,
        #               selected = Van_def$grp1),
        #   uiOutput("Van_c1sub1.ui"),
        #   actionButton("Van_c1sub1all", "Select all groups", class = "btn btn-primary"),
        #   actionButton("Van_c1sub1non", "Deselect all groups", class = "btn btn-primary")
        # ), br(),
        br(),
        uiOutput("UI_test_select"),

        # radioButtons(inputId = ns("select_test_normalization"),
        #              label = "Normalization:",
        #              choices = c("Raw data" = "raw",
        #                          "Total area normalization" = "tot_area"),
        #              selected = "tot_area"),
        # radioButtons(inputId = ns("select_test_transformation"),
        #              label = "Transformation:",
        #              choices = c("No transformation" = "none",
        #                          "Log10 transformation" = "log10"),
        #              selected = "none"),
        # checkboxInput(inputId = ns("test_cor_pvalue"),
        #               label = "Show corrected p-value",
        #               value = FALSE),
        # this dynamicall makes all of our test/group selectors
        uiOutput(outputId = ns("comp_group_selection"))
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
             #uiOutput(ns("UI_viz_output")),
              plotlyOutput(outputId = ns("volcano_plot"),height = "800px"),
              plotlyOutput(outputId = ns("test_boxplot"))


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


 ) #taglist

}

#' pg_vis_comp Server Functions
#'
#' @noRd
#' @importFrom plotly renderPlotly plotlyOutput plot_ly add_markers event_data
#'
mod_pg_vis_comp_server <- function(id,rv_in, p){
  moduleServer( id, function(input, output, session){
    ns <- session$ns


    data_point_sz = .65
    plot_size = "Small"
    font_size = "Small"
    color_schemes = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple")
    color_scheme = color_schemes[2]


    #### compare samples
    test_result <- reactive({
      req(
          rv_in$de,
          input$SI_comp_type,
          input$SI_comp_fact,
          input$SI_color_grp,
          input$RB_select_test,
          input$select_test_normalization,
          input$select_test_transformation)


      diff_exp <- isolate(rv_in$de)

      # filter according to current "cases"
      # TODO: fixe test_type

      de <- diff_exp %>%
        dplyr::filter(test_type == input$RB_select_test &
                        group == input$SI_comp_fact &
                        condition == input$SI_comp_type) %>%
        dplyr::mutate(f=1,
                      significant = pvals_adj < 0.01,
                      point_color =  dplyr::recode(as.character(.data$significant), "TRUE" = "#FF7F00", "FALSE" = "#1F78B4"))

      return(de)
    })


    output$UI_test_select <- renderUI({
      req(rv_in$de,
          rv_in$defaults)

      radioButtons(inputId = ns("RB_select_test"),
                   label = "Test:",
                   choices = factors( isolate(rv_in$de$select_test) ),
                   selected = rv_in$default$test[1])


    })
      # need to make in_data a tibble? based on rv_in$ad$X
        # prep_test_data <- lipid_data %>%
        # ## 1. filter out bad data
        # filter(.data$keep == TRUE) %>%
        #
        # ## 2. choos the rows with the two groups being compared.
        # rename(my_group_info = !!sym(group)) %>%
        # filter(.data$my_group_info == group1_name |
        #          .data$my_group_info == group2_name) %>%
        #
        # ## 3. keep minimal number of columns - ID, omics-name, omics-class(agg), sampel name, group info, area,
        # select(.data$my_id, .data$ShortLipidName, .data$LipidClass, .data$sample_name, .data$my_group_info, .data$area) %>%
        #
        # ## 4a. Normalizaiton..
        # # total area normalisation
        # group_by(.data$sample_name) %>%
        # mutate(norm_area = .data$area / sum(.data$area)) %>%
        # ungroup() %>%
        # # select which normalization to use
        # mutate(value = case_when(
        #   normalization == "raw" ~ .data$area,
        #   normalization == "tot_area" ~ .data$norm_area
        # )) %>%
        #
        # ## 4b. scaling
        # # do transformations and select which transformation to keep
        # mutate(value = case_when(
        #   transformation == "none" ~ .data$value,
        #   transformation == "log10" ~ log10(.data$value + 1) # the +1 is correct for any zero's
        # )) %>%
        #
        # # 5. column clean up.
        # # remove the 2 area columns
        # select(-.data$area, -.data$norm_area) %>%
        #
        # ## 6. make nexted "test_data" column
        # nest(test_data = c(.data$sample_name, .data$my_group_info, .data$value)) %>%
        #
        # ## 7. calculate fold count + log fold cont.
        # mutate(fc = map_dbl(.x = .data$test_data,
        #                     .f = ~ mean(.x$value[.x$my_group_info == group1_name]) / mean(.x$value[.x$my_group_info == group2_name])),
        #        fc_log2 = log2(.data$fc))
        # ## 8. test for difference
        # result <- switch(test,
        #       "ttest" = do_ttest(lipid_data = prep_test_data),
        #       "mwtest" = do_mwtest(lipid_data = prep_test_data))

      # colorby_group <- input$SI_comp_type

      # browser()
      #
      # pg_volc_ly <- function(in_data, pvalue_adjust = FALSE, title = "")
      #
      #
      # # results_test <- do_stat_test(in_meta = isolate(rv_in$meta),
      # #                              in_data = isolate(rv_in$ad$X),
      # #                              group = input$test_select_group,
      # #                              group1_name = input$SI_comp_type,
      # #                              group2_name = input$test_group2,
      # #                              colorby_group = colorby_group,
      # #                              normalization = input$select_test_normalization,
      # #                              transformation = input$select_test_transformation,
      # #                              test = input$select_test)
      #
      # return(results_test)


    # create some ui output
    output$comp_group_selection <- renderUI({
      req(rv_in$config,
          rv_in$default)


      cfg <- isolate(rv_in$config)

   # subs_label <- paste0("select ",isolate(input$SI_subset),"s: ")
      #

      # checkboxGroupInput( ns("CB_sub_inner1"),
      #                     label = subs_label,
      #                     inline = TRUE,
      #                     choices = subs,
      #                     selected = subs)


      to_return <-  tagList(
        selectInput(inputId = ns("SI_comp_type"),
                    label = "Comp type",
                    choices =  strsplit(cfg[ID=="diff_exp_conds"]$fID, "\\|")[[1]],
                    selected = rv_in$default$comp_type),
        # uiOutput(outputId = ns("test_vs_groups")),
        selectInput(inputId = ns("SI_comp_fact"),
                    label = "compare :",
                    choices = NULL ),

        selectInput(inputId = ns("SI_color_grp"),
                    label = "Color by:",
                    choices = rv_in$config[grp == TRUE]$UI,
                    selected = rv_in$default$color_grp)
      )

      return(to_return)
    })

    # output$test_vs_groups <- renderUI({
    #   tagList(
    #     selectInput(inputId = ns("test_group1"),
    #                 label = "group 1:",
    #                 choices = "none"),
    #     HTML("<center><b>vs.</b></center>"),
    #     selectInput(inputId = ns("test_group2"),
    #                 label = "group 2:",
    #                 choices = "none")
    #   )
    # })



    observeEvent(input$SI_comp_type, {
      req(rv_in$config,
          rv_in$default)

      browser()

      # get the groups
      cfg <- isolate(rv_in$config)
      subs <- strsplit(cfg[UI == "diff_exp_groups"]$fID, "\\|")[[1]]
      subs <- sort(subs)

      subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
      names(subs2) <- seq_along(subs2)
      subs2 <- subs2[!is.na(subs2)]
      #subs_label <- paste0("select ",isolate(input$SI_subset),"s: ")
      if (length(subs2)>0) {
        subs <- subs[as.numeric(names(sort(subs2)))]
      }

      updateSelectInput(inputId = "SI_comp_fact",
                        label = "compare :",
                        choices = subs,
                        selected = rv_in$default$comp_fact[1])

      # updateSelectInput(inputId = "test_group2",
      #                   label = "Group 2:",
      #                   choices = subs,
      #                   selected = subs[2])

    })


    output$volcano_plot <- renderPlotly({
      req(test_result,
          input$SI_comp_type,
          input$SI_comp_fact)
      browser()

      colorby_group <- input$SI_comp_type
      de_local <- test_result()

      if( dim(de_local)[1]>0 ) {
        # pg_volcano_ly(in_data = test_result(),
        #              pvalue_adjust = input$test_cor_pvalue,
        #              title = paste0(input$SI_comp_type, " vs ", input$test_group2))
        pg_volc_ly(de=de_local , title = paste0(input$SI_comp_type, " || ", input$SI_comp_fact) )
      }
    })


    # output$test_boxplot <- renderPlotly({
    #   req(test_result)
    #
    #   if(!is.null(test_result())) {
    #     # capture the click event
    #     # this contains a column with the shortlipidname (column name: customdata)
    #     my_data <- event_data(event = "plotly_click",
    #                           source = "volcano_plot_click")
    #
    #     if(!is.null(my_data)) {
    #       # restructure data
    #       plot_data <- test_result() %>%
    #         filter(.data$ShortLipidName %in% my_data$customdata) %>%
    #         select(.data$test_data) %>%
    #         unnest(.data$test_data)
    #
    #       # show the boxplot
    #       box_plot(lipid_data = plot_data,
    #                title = paste0("Lipid: ", my_data$customdata))
    #     }
    #   }
    #
    #
    # })


#     #output$plot_volcano_out <- renderPlot({
#     output$plot_volcano_out <- renderPlotly({
# browser()
#       print("renderPlot volcano")
#       req(p$omics_list,
#           p$observ_grp,
#           p$observ_subsel,
#           p$observ_y_raw)
#
#
#       in_conf <- rv_in$config$meta
#       in_meta <- rv_in$meta
#
#
#       in_fact <- p$observ_grp
#
#
#       dat_source = rv_in$config$mat[fIDloc == p$observ_y_raw]$ID
#       dat_key = rv_in$config$mat[fIDloc == p$observ_y_raw]$fID
#       #  probably add a "matrix" column to config...
#       # maybe use the data_source to indicate if we want raw or layers... short circuit for now
#       in_data <- isolate(rv_in$ad$X)  # do we need to isolate it??
#
#
#       in_quant <- "X" #dat_key #(maybe) just observ_y_raw
#
#       # these are the "groups" to show on the x axis
#       in_group <- in_fact
#       #in_subset1 <- p$observ_subsel
#       in_subset1 <- p$observ_grp
#       in_subset2 <- p$observ_subsel
#
#       # these are the groups to show on thye y axis
#       all_omics <- rv_in$ad$var_names
#       names(all_omics) <- all_omics
#
#       all_obs <- rv_in$ad$obs_names
#       names(all_obs) <- all_obs
#
#
#       #
#       #
#       # pg_volcano_ly(in_conf, in_meta, p$omics_list$value, in_group, input$RB_heat_plot_type,
#       #                   p$observ_grp, p$observ_subsel, in_data, all_omics,
#       #                   input$CB_scale, input$CB_cluster_rows, input$CB_cluster_cols,
#       #                   color_scheme, plot_size)
#       #
#       # pg_bubble_heatmap(in_conf, in_meta, p$omics_list$value, in_group, input$RB_heat_plot_type,
#       #                   p$observ_grp, p$observ_subsel, in_data, all_omics,
#       #                   input$CB_scale, input$CB_cluster_rows, input$CB_cluster_cols,
#       #                   color_scheme, plot_size)
#       return(NULL)
#       })
#
#
    # output$UI_heatmap <- renderPlotly({
    #   print("renderUI UI volcano")
    #
    #   #plotOutput(ns("plot_volcano_out"), height = pList3[plot_size])
    #   plotlyOutput(ns("plot_volcano_out"), height = pList3[plot_size])
    #
    #
    # })





    # output$plot_box_out <- renderPlot({
    #   #req(rv_in$ad)
    #   req(p$omics_list,
    #       p$observ_grp,
    #       p$observ_subsel,
    #       p$observ_y_raw)
    #
    #   # observ_y_raw = NULL,
    #   # observ_y_comp = NULL,
    #
    #   in_conf <- rv_in$config$meta
    #   in_meta <- rv_in$meta
    #
    #
    #   in_fact <- p$observ_grp
    #   # in_fact <- rv_in$config$meta[[p$observ_grp]]
    #   #
    #   # in_fact <- in_meta[[p$observ_grp]]
    #
    #   #rv_in$ad$obs[[p$observ_grp]] %in% p$observ_subsel
    #
    #   # could also extract from the name instead of lookkup
    #   #       subs <- strsplit(rv_in$observ_y_raw, "\\|")[[1]]
    #   # rv_in$ad[[subs[1]]][[subs[2]]]
    #   dat_source = rv_in$config$mat[fIDloc == p$observ_y_raw]$ID
    #   dat_key = rv_in$config$mat[fIDloc == p$observ_y_raw]$fID
    #   #dd <- isolate(rv_in$ad[[rv_in$config$mat[fID == rv_in$observ_y_raw]$ID]][[rv_in$observ_y_raw]])
    #   if (dat_source == "obs") {
    #     in_data <- isolate(rv_in$ad$obs[[dat_key]])
    #   } else if (dat_source == "var") {
    #     in_data <- isolate(rv_in$ad$var[[dat_key]])
    #     names(in_data) <- ad$var_names
    #   } else {
    #     in_data <- NULL
    #     print('boxplots only for obs and var ?@!?')
    #     #TODO:  spawn warning box
    #     return(NULL)
    #   }
    #   in_quant <- dat_key #(maybe) just observ_y_raw
    #
    #   pg_violin_box(in_conf, in_meta, in_fact, in_quant,
    #                 p$observ_grp, p$observ_subsel,
    #                 in_data, p$omics_list$value, input$RB_plot_type, input$CB_show_data_points,
    #                 data_point_sz, font_size)
    #
    # })
    #
    # output$UI_viz_output <- renderUI({
    #   plotOutput(ns("plot_box_out"), height = pList2[plot_size])
    # })


    ###
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



  })
}

## To be copied in the UI
# mod_pg_vis_comp_ui("pg_vis_comp_ui_1")

## To be copied in the server
# mod_pg_vis_comp_server("pg_vis_comp_ui_1")

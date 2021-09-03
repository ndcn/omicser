#TODO: clean up unused functions / code here
#

#' @title Do a statitical test
#'
#' @description Do a statitical test, t-test or Mann-Whitney U test
#'
#' @param in_data tibble in tidy format
#' @param group what column contains the group info
#' @param group1_name is the name of group 1
#' @param group2_name is the name of group 2
#' @param normalization what normalization to use, none (raw data) or total area normalization
#' @param transformation what transformation to use
#' @param test do a t-test or Mann-Whitney U test
#'
#' @return a tibble ready for statistical testing
#'
#' @import tidyselect
#' @importFrom dplyr select filter mutate rename group_by ungroup
#' @importFrom tidyr pivot_wider everything nest
#' @importFrom purrr map_dbl
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @author Rico Derks
#'
do_stat_test <- function(in_meta, in_data, group, group1_name, group2_name, colorby_group, normalization = c("raw", "tot_area"), transformation = c("none", "log10"),
                         test) {

  if(group1_name != "none" &
     group2_name != "none" &
     group1_name != group2_name &
     group1_name %in% in_meta[[group]] &
     group2_name %in% in_meta[[group]]) {

    if (length(which(colnames(in_meta)=="sampleID"))){
      in_meta <- in_meta[,-1]
    }

    l_id <-rownames(in_data)
    l_omic <-colnames(in_data)

    # meta <- tibble::as_tibble(in_meta) %>%
    #   rename(my_group_info = !!rlang::sym(group),
    #          colorby_group = !!rlang::sym(colorby_group))  %>%
    #   mutate(my_group_info = as.character(my_group_info),
    #          colorby_group = as.character(colorby_group))
    #
    meta <- tibble::as_tibble(in_meta) %>%
            mutate(my_group_info = as.character(!!rlang::sym(group)),
                   colorby_group = as.character(!!rlang::sym(colorby_group)))

    test_data <- tibble::as_tibble(in_data)
    test_data$sampleID <- rownames(in_data)
    #this takes a while.. we should do this immediately and once when the database is loaded

    prep_test_data <- dplyr::inner_join(meta,test_data, on=sampleID) %>%
      dplyr::filter(.data$my_group_info == group1_name | .data$my_group_info == group2_name) #%>%

    prep_test_data <- prep_test_data %>%
      tidyr::pivot_longer(
        cols = l_omic,
        names_to = "omic",
        values_to = "value",
      ) %>%

      select(.data$sampleID, .data$omic, .data$value, .data$my_group_info , .data$colorby_group)

      # do transformations and select which transformation to keep
    prep_test_data <- prep_test_data %>%
      mutate(value = case_when(
        transformation == "none" ~ .data$value,
        transformation == "log10" ~ log10(.data$value + 1) # the +1 is correct for any zero's
      ))

    prep_test_data <- prep_test_data %>%
      nest(test_data = c(.data$sampleID, .data$my_group_info, .data$value))

    prep_test_data2 <- prep_test_data %>%
      mutate(fc = map_dbl(.x = .data$test_data,
                          .f = ~ mean(.x$value[.x$my_group_info == group1_name]) / mean(.x$value[.x$my_group_info == group2_name])),
             fc_log2 = log2(.data$fc))


    result <- switch(test,
                     "ttest" = do_ttest(in_data = prep_test_data),
                     "mwtest" = do_mwtest(in_data = prep_test_data))
  } else {
    result <- NULL
  }

  return(result)
}

# Function - Volcano Plot --------------------------------------------------
#' @title Create volcano plot
#'
#' @description Create volcano plot.
#'
#' @param in_data tibble with all the lipid data and test data
#' @param pvalue_adjust show the corrected p value, default is FALSE
#' @param title title of the plot
#'
#' @return plotly object
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate case_when
#' @importFrom rlang .data
#' @importFrom plotly plot_ly add_markers layout event_register
#' @importFrom grDevices rainbow
#'
#' @author Andy Henrie modified from Rico Derks
pg_volcano_ly <- function(in_data, pvalue_adjust = FALSE, title = "") {



  # pg_volcano_ly <- function(in_conf, in_meta, in_omics, in_fact, plot_type,
  #                           in_grp, in_subsel, in_data, all_omics, in_do_scale, in_clust_row, in_clust_col,
  #                           color_scheme, plot_size,
  #                           save = FALSE) {

    #
    # (in_conf, in_meta, p$omics_list$value, in_group, input$RB_heat_plot_type,
    #                           p$observ_grp, p$observ_subsel, in_data, all_omics,
    #                           input$CB_scale, input$CB_cluster_rows, input$CB_cluster_cols,
    #                           color_scheme, plot_size)
  # create y-axis title
  y_title <- ifelse(pvalue_adjust == FALSE,
                    "-log10(p value)",
                    "-log10(cor. p value)")

  # create the plot
  p <- in_data %>%
    mutate(show_p = case_when(
      pvalue_adjust == FALSE ~ .data$p_log10,
      pvalue_adjust == TRUE ~ .data$p_log10_adj
    )) %>%

    plot_ly(x = ~fc_log2,
            y = ~show_p,
            text = ~omic,
            colors = rainbow(n = 100),
            customdata = in_data$omic,
            source = "volcano_plot_click") %>%

    add_markers(color = ~colorby_group, size = 3) %>%

    layout(xaxis = list(zeroline = FALSE,
                        title = "log2(fold change)"),
           yaxis = list(title = y_title),
           shapes = list(vline(-1),
                         vline(1),
                         hline(-log10(0.05))),
           legend = list(orientation = "h"),
           title = list(text = title,
                        x = 0)) %>%

    event_register(event = "plotly_click")

  return(p)
}

vline <- function(x = 0, color = "blue") {
  list(
    type = "line",
    y0 = 0,
    y1 = 1,
    yref = "paper",
    x0 = x,
    x1 = x,
    line = list(color = color,
                width = 1,
                dash = "dash")
  )
}

hline <- function(y = 0, color = "blue") {
  list(
    type = "line",
    x0 = 0,
    x1 = 1,
    xref = "paper",
    y0 = y,
    y1 = y,
    line = list(color = color,
                width = 1,
                dash = "dash")
  )
}


#' @title Do a t-test on all lipids
#'
#' @description Do a t-test on all lipids.
#'
#' @param in_data tibble in tidy format, already nested
#'
#' @details A t-test will be done for each lipid. Also the p-value will be corrected
#'     for multiple testing.
#'
#' @return a tibble with the results and the following columns
#'     model_test contains the model information for each lipid (nested)
#'     pvalue contains the uncorrected pvalue
#'     pvalue_cor contains the corrected pvalue
#'     fc the fold change estimate1/estimate2
#'
#' @import tidyselect
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom purrr map map_dbl
#' @importFrom broom tidy
#' @importFrom stats p.adjust
#'
#' @author Rico Derks
#'
do_ttest <- function(in_data) {
  results <- in_data %>%
    mutate(model_test = map(.x = .data$test_data,
                            .f = ~ broom::tidy(t.test(formula = value ~ my_group_info,
                                                      data = .x)) %>%
                              # calculate the fold change
                              mutate(fc = estimate1 / estimate2)),
           pvalue = map_dbl(.x = .data$model_test,
                            .f = ~ .x$p.value),
           pvalue_adj = p.adjust(.data$pvalue,
                                 method = "BH"),
           p_log10 = -log10(.data$pvalue),
           p_log10_adj = -log10(.data$pvalue_adj)
           )

  return(results)
}

#' @title Do a Mann-Whitney U test on all lipids
#'
#' @description Do a Mann-Whitney U test on all lipids
#'
#' @param in_data tibble in tidy format, already nested
#'
#' @details A Mann-Whitney U test will be done for each lipid. Also the p-value will be corrected
#'     for multiple testing.
#'
#' @return a tibble with the results and the following columns
#'     model_test contains the model information for each lipid (nested)
#'     pvalue contains the uncorrected pvalue
#'     pvalue_cor contains the corrected pvalue
#'     fc the fold change estimate1/estimate2
#'
#' @import tidyselect
#' @importFrom dplyr mutate group_by summarise pull
#' @importFrom tidyr pivot_wider
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom purrr map map_dbl
#' @importFrom broom tidy
#' @importFrom stats p.adjust
#'
#' @author Rico Derks
#'
do_mwtest <- function(in_data) {
  results <- in_data %>%
    mutate(model_test = map(.x = .data$test_data,
                            .f = ~ broom::tidy(wilcox.test(formula = value ~ my_group_info,
                                                           data = .x))),
           pvalue = map_dbl(.x = .data$model_test,
                            .f = ~ .x$p.value),
           pvalue_adj = p.adjust(.data$pvalue,
                                 method = "BH"),
           p_log10 = -log10(.data$pvalue),
           p_log10_adj = -log10(.data$pvalue_adj))

  return(results)
}

#' @title Create box plot
#'
#' @description Create box plot.
#'
#' @param in_data tibble with all the lipid data and test data
#' @param title title of the plot
#'
#' @return plotly object
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select
#' @importFrom tidyr unnest
#' @importFrom rlang .data
#' @importFrom plotly plot_ly layout hide_legend
#'
#' @author Rico Derks
#'
box_plot <- function(in_data, title = "") {
  # create the plot
  plt <- in_data %>%
    plot_ly(x = ~my_group_info,
            y = ~value,
            text = ~sample_name,
            color = ~my_group_info,
            type = "box",
            boxpoints = "all",
            jitter = 0.4,
            pointpos = 0) %>%
    layout(title = list(text = title, x = 0),
           yaxis = list(title = "Value",exponentformat = "E"),
           xaxis = list(title = "Group")) %>%
    hide_legend()

  return(plt)
}


#
#   #(in_data,var1, val1, group=NULL, genes, ...) {
#
#   # assumes:  adj.P.Val, logFC, gene_name
#
#   data <- df %>%  dplyr::filter( df[,which(colnames(df)==var1)]==val1)
#
#   color_map <-  c(Down='skyblue3', Not_Sign='gray', Up='darkorange')
#
#   # Add column Color with labels for coloring the data
#   # data <- Sol_all
#   stroke_map <- c(mod='rgb(50,205,50)') # defines the color for selected proteins
#
#   # create new column "sel" that encodes the key for selected proteins
#   data$sel <- "none"
#   # include "mod" key only for those genes that are selected by now
#   data$sel[data$gene_name %in% genes ] <- "mod"
#
#   data <- data %>%
#     mutate(Color = ifelse(logFC > 0, "Up", "Down")) %>%
#     mutate(Color = ifelse(abs(adj.P.Val) < 0.05, Color, "Not_Sign"))
#
#   plot1 <- plot_ly(data,
#                    x = ~logFC,
#                    y = ~(-log10(adj.P.Val)),
#                    color = ~Color,
#                    colors = color_map,
#                    stroke = ~sel,
#                    strokes = stroke_map,
#                    type="scatter",
#                    ids = ~gene_name,
#                    hoverinfo= "text",
#                    text = ~paste('Gene: ', gene_name,
#                                  '<br>logFC: ', round(logFC,3),
#                                  '<br>FDR: ', round(adj.P.Val,10)
#                    ),
#                    textposition = "middle left",
#                    mode= "markers",
#                    marker=list(symbol = "0",
#                                size=5,
#                                bgcolor = "#e5e5e5",
#                                opacity = 0.5,
#                                line = list(width = 2)
#                    ),
#                    source = "volcano_plot",
#                    customdata = ~gene_name
#                 ) %>% plotly::layout(showlegend = TRUE,
#                          legend = list(legend = 1, itemsizing = "constant"),
#                          title = list(text=~paste(group, "in", val1), x = 0),
#                          yaxis = list(title = 'FDR (Ajusted p-value)'),
#                          xaxis = list(title = 'Log(2) Fold Change', zeroline=FALSE) #showgrid = FALSE
#                   )
#
# }
#


#' #' pg_volcano
#' #' @title Create volcano plot
#' #'
#' #' @description Create volcano plot.
#' #'
#' #' @param in_data tibble with all the lipid data and test data
#' #' @param pvalue_adjust show the corrected p value, default is FALSE
#' #' @param title title of the plot
#' #'
#' #' @return plotly object
#' #'
#' #' @importFrom magrittr %>%
#' #' @importFrom dplyr mutate case_when
#' #' @importFrom rlang .data
#' #' @importFrom plotly plot_ly add_markers layout event_register
#' #' @importFrom grDevices rainbow
#' #'
#' #' @author Andy Henrie (derived from Rico Dirks lipidomics)
#'
#' pg_volcano <- function(in_conf, in_meta, in_omics, in_fact, plot_type,
#'                               in_grp, in_subsel, in_data, all_omics, in_do_scale, in_clust_row, in_clust_col,
#'                               color_scheme, plot_size,
#'                               save = FALSE){
#'   if(is.null(in_grp)){in_grp = in_conf$UI[1]}
#'   # Identify genes that are in our dataset
#'   omic_list = get_omic_list(in_omics, all_omics)
#'   omic_list = omic_list[present == TRUE]
#'   shiny::validate(need(nrow(omic_list) <= 50, "More than 50 genes to plot! Please reduce the gene list!"))
#'   shiny::validate(need(nrow(omic_list) > 1, "Please input at least 2 genes to plot!"))
#'
#'   # Prepare ggData
#'   #
#'
#'
#' volcano_plot <- function(in_data, pvalue_adjust = FALSE, title = "") {
#'   # create y-axis title
#'
#'   y_title <- ifelse(pvalue_adjust == FALSE,
#'                     "-log10(p value)",
#'                     "-log10(cor. p value)")
#'
#'   # create the plot
#'   p <- in_data %>%
#'     mutate(show_p = case_when(
#'       pvalue_adjust == FALSE ~ .data$p_log10,
#'       pvalue_adjust == TRUE ~ .data$p_log10_adj
#'     )) %>%
#'     plot_ly(x = ~fc_log2,
#'             y = ~show_p,
#'             text = ~ShortLipidName,
#'             colors = rainbow(n = 100),
#'             customdata = in_data$ShortLipidName,
#'             source = "volcano_plot_click") %>%
#'     add_markers(color = ~LipidClass,
#'                 size = 3) %>%
#'     layout(xaxis = list(zeroline = FALSE,
#'                         title = "log2(fold change)"),
#'            yaxis = list(title = y_title),
#'            shapes = list(vline(-1),
#'                          vline(1),
#'                          hline(-log10(0.05))),
#'            legend = list(orientation = "h"),
#'            title = list(text = title,
#'                         x = 0)) %>%
#'     event_register(event = "plotly_click")
#'
#'   return(p)
#' }
#'
#' vline <- function(x = 0, color = "blue") {
#'   list(
#'     type = "line",
#'     y0 = 0,
#'     y1 = 1,
#'     yref = "paper",
#'     x0 = x,
#'     x1 = x,
#'     line = list(color = color,
#'                 width = 1,
#'                 dash = "dash")
#'   )
#' }
#'
#' hline <- function(y = 0, color = "blue") {
#'   list(
#'     type = "line",
#'     x0 = 0,
#'     x1 = 1,
#'     xref = "paper",
#'     y0 = y,
#'     y1 = y,
#'     line = list(color = color,
#'                 width = 1,
#'                 dash = "dash")
#'   )
#' }


pg_volc_ly <- function(de, title = "") {


  # create the plot
  # p <- in_data %>%

    plt <- plot_ly(
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

    plt <- plt %>%
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
        title = title,
        xaxis = list(title = "Effect (logFC)", range = c(-4, 4)),
        yaxis = list(title = "-log10 p-value", range = c(-0.1, 10.25))
      ) %>%
      # Disable the legend click since our traces do not correspond to the
      # actual legend labels
      htmlwidgets::onRender("function(el,x){el.on('plotly_legendclick', function(){ return false; })}") %>%
      plotly::config(displayModeBar = FALSE)


    return(plt)

#  return(p)
}



# Function - Heat Map DIA -------------------------------------------------
heat_DIA <- function(..., proteins, dimension_x) {
  # if (identical(proteins, character(0))) {
  #     print("Empty")
  #     return(NULL)
  # } else {
  #proteins <- c("Ppt1", "Cul4b", "Camsap1")
  # filter data based on selected proteins (MsC-Data and Transcriptomics data)
  allDat2 <- transcriptomics_long %>% filter(Genes %in% proteins)
  allDat2 <- rename(allDat2,
                    Comparison..group1.group2. = label,
                    gene.name = Genes)
  allDat2$datatype <- "transcriptome"
  allDat <- filter_prot(proteins, "DIA") %>%
              subset(select = c("Comparison..group1.group2.", "logFC",  "adj.P.Val", "UniProtIds", "gene.name", "ProteinDescriptions"))
  allDat$datatype <- "proteome"
  # # Reorder Variable-levels for visualization
  # allDat$Comparison..group1.group2. <- as.character(allDat$Comparison..group1.group2.)
  # allDat$Comparison..group1.group2. <- factor(allDat$Comparison..group1.group2., levels = c("o / y", "g / y", "g / o"))

  mult_fact <- length(unique(allDat$gene.name))
  rwb <- colorRampPalette(colors = c('skyblue3', 'white', 'darkorange'))

  #proteins <- c("Ppt1", "Cul4b", "Camsap1")
  plot <- allDat %>% ggplot( aes(x = Comparison..group1.group2.,
                                 y = gene.name,
                                 size = -log10(adj.P.Val),
                                 colour = logFC)) +
    geom_point(data = allDat2) +  # do not change order of the arguments geom_point and scale_x… since this would change the order of variables on x axis
    geom_point() +
    facet_wrap(.~datatype,
               scales="free_x",
               ncol = 2, dir = "h",
               drop=TRUE, strip.position = "bottom",
               labeller = label_wrap_gen(width = 5, multi_line = TRUE)) +
    xlab("") + #needed for subsequent scale_x_discrete-command, as too long labels need to be split
    scale_x_discrete("Condition",
                     drop=TRUE, na.translate = FALSE,
                     labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n")) +
    scale_colour_gradientn(colours = rwb(100),
                           limits=c(-1.5, 1.5),
                           oob = scales::squish) +
    scale_size(name = "-log10(adj.P.Val)",
               breaks = c(1.3,2,4,6)) +
    theme(panel.grid.major.y = element_line(colour="grey", size=0.1, inherit.blank = FALSE),
          panel.background = element_blank(),
          panel.spacing = unit(0, "lines"),
          axis.ticks=element_blank(),
          strip.text.x = element_text(size=12),
          strip.background = element_blank(),
          axis.text = element_text(size = 12, inherit.blank = FALSE),
          axis.title = element_text(size = 17, inherit.blank = FALSE),
          axis.text.x.bottom = element_text(margin = margin(t = 15, r = 15, b = 15, l = 15, unit = "pt"))) +
    labs(x="Condition", y="Proteins")
  #annotate("text", x = c(2, 6), y = c(0.6,0.6), label = c("proteome", "transcriptome") , color="gray30", size=4)

  if(length(proteins) <= 8) {
    plot <- plot + theme(legend.position = "bottom",
                         legend.direction = "horizontal",
                         legend.text = element_text(angle = 90, vjust = 0.9))
  }
  if(!dimension_x){
    plot <- plot + theme(axis.text.x = element_text(angle = 90))
  }
  return(plot)
  # }
}


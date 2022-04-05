#TODO: clean up unused functions / code here
#

# NOT CALLED/WORKING
pg_volc_ly <- function(de, title = "") {


  # create the plot
  # p <- in_data %>%

    plt <- plot_ly( x = de$logfoldchanges,
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



# NOT CALLED/WORKING
curr_volc_ly <- function(de, title = "") {

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
  # point_color (mapped to significant < threshold)

  #TODO: color the target omics a different color

  # TODO:  make this a function and move to fct_pg_vis_comp
  # volc <- pg_volc_ly(de=de_local , title = title_text)
  # return(volc)
  volc <- plot_ly(
    x = de$logfoldchanges,
    y = -log10(de$pvals),
    name = "FDR > 0.01",
    type = "scatter",
    showlegend = FALSE,
    mode = "markers",
    # Hovertext
    text = paste(de$names,
                 "</br></br>log2FC: ",
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
    plotly::add_markers(x=3.0, y = 9.5, color = I("#1F78B4"), showlegend = FALSE, hoverinfo = "skip") %>%
    plotly::add_annotations(x=3.0, y=9.5, xref = "x", yref = "y", text = "FDR > 0.01",
                            xanchor = 'left', showarrow = F, xshift = 10) %>%
    # Orange/significant
    plotly::add_markers(x=3.0, y = 9.0, color = I("#FF7F00"), showlegend = FALSE, hoverinfo = "skip") %>%
    plotly::add_annotations(x=3.0, y=9.0, xref = "x", yref = "y", text = "FDR < 0.01",
                            xanchor = 'left', showarrow = F, xshift = 10) %>%
#
#
#     %>%
#     add_markers(x
#                 y
#                 color = ~LipidClass,
#                 size = 3) %>%

    plotly::layout(
      title = list(text = title,
                     x = 0),
      xaxis = list(zeroline = FALSE,
                   title = "log2FC"), #, range = c(-8, 8)),
      yaxis = list(title = "-log10 p-val"), #, range = c(-0.1, 60.25))
      shapes = list(vline(-1),
                    vline(1),
                    hline(-log10(0.05))),
      legend = list(orientation = "h") ) %>%
    # Disable the legend click since our traces do not correspond to the
    # actual legend labels
    htmlwidgets::onRender("function(el,x){el.on('plotly_legendclick', function(){ return false; })}") %>%
    plotly::config(displayModeBar = FALSE) %>%
    plotly::event_register('plotly_click')

  return(volc)

}


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
#@importFrom grDevices rainbow
#'
#'
#' @author Rico Derks
#'
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
# point_color (mapped to significant < threshold)

volcano_plot <- function(in_data, pvalue_adjust = FALSE, title = "") {
  # create y-axis title
  y_title <- ifelse(pvalue_adjust == FALSE,
                    "-log10(p value)",
                    "-log10(cor. p value)")

  # create the plot
  p <- in_data %>%
    mutate(show_p = case_when(
      pvalue_adjust == FALSE ~ .data$pvals,
      pvalue_adjust == TRUE ~ .data$pvals_adj
    )) %>%
    plot_ly(x = ~logfoldchanage,
            y = ~show_p,
            text = ~names,
            alpha = 0.3,
            #colors = rainbow(n = 100),
            customdata = in_data$names,
            source = "volcano_plot_click") %>%
    add_markers(color = ~LipidClass,
                size = 3) %>%
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


#' @title create plotly volcano plot from ggplot2
#'
#' @param in_data data
#' @param pvalue_adjust adjust p-value
#' @param title title of the plot
#'
#' @return plotly object
#'
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal scale_color_manual geom_vline geom_hline
#' @importFrom rlang .data
#'
#' @noRd
#'
volc_ggplotly <- function(in_data, pvalue_adjust = FALSE, title = "") {

  # plot adding up all layers we have seen so far
  ggplot(data = in_data,
         aes(x = .data$log2FoldChange,
             y = -log10(.data$pvalue),
             col = .data$diffexpressed,
             label = .data$delabel)) +
    geom_point() +
    theme_minimal() +
    scale_color_manual(values = c("blue", "black", "red")) +
    geom_vline(xintercept = c(-0.6, 0.6),
               col = "red") +
    geom_hline(yintercept = -log10(0.05),
               col = "red")

}


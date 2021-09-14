#TODO: clean up unused functions / code here
#

# NOT CALLED/WORKING
pg_volc_ly <- function(de, title = "") {


  # create the plot
  # p <- in_data %>%

    plt <- plot_ly(
                  x = de$logfoldchanges,
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


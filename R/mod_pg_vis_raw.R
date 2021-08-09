# Function to extract legend
g_legend <- function(a.gplot){
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

# Get omic list
get_omic_list <- function(inp, inp_omics){
  omic_list = data.table(omic = inp,present=TRUE)# unique(trimws(strsplit(inp, ",|;|")[[1]])),

  omic_list[!omic %in% (inp_omics)]$present = FALSE
  return(omic_list)
}

# Plot gene expression bubbleplot / heatmap
scBubbHeat <- function(inpConf, inpMeta, inp, inpGrp, inpPlt,
                       inpsub1, inpsub2, in_data, inpGene, inpScl, inpRow, inpCol,
                       inpcols, inpfsz, save = FALSE){
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]}
  # Identify genes that are in our dataset
  geneList = get_omic_list(inp, inpGene)
  geneList = geneList[present == TRUE]
  shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!"))
  shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!"))

  # Prepare ggData
  #h5file <- H5File$new(inpH5, mode = "r")
  #h5data <- h5file[["grp"]][["data"]]
  ggData = data.table()
  for(iGene in geneList$omic){
    tmp = inpMeta[, c("sampleID", inpConf[UI == inpsub1]$ID), with = FALSE]
    colnames(tmp) = c("sampleID", "sub")
    tmp$grpBy = inpMeta[[inpConf[UI == inpGrp]$ID]]
    tmp$geneName = iGene
    tmp$val = in_data[,inpGene[iGene]] #h5data$read(args = list(inpGene[iGene], quote(expr=)))
    ggData = rbindlist(list(ggData, tmp))
  }
 # h5file$close_all()


  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){
    ggData = ggData[sub %in% inpsub2]
  }
  shiny::validate(need(uniqueN(ggData$grpBy) > 1, "Only 1 group present, unable to plot!"))

  # Aggregate:  exp(mean) + proportion -> log(val)
  # disabling expm1 and log1p because we will be putting "NORMALIZED" quantities in...

  if (max(abs(ggData$val),na.rm=TRUE) < 700.0 ){
    fwd_scale <- expm1
    inv_scale <- log1p
  } else {
    fwd_scale <- function( input ){ input }
    inv_scale <- function( input ){ input }
  }

  # remove NA
  ggData <- ggData[!is.na(val)]
  ggData$val = fwd_scale(ggData$val)
  ggData = ggData[, .(val = mean(val), prop = sum(val>0) / length(sampleID)),
                  by = c("geneName", "grpBy")]
  # ggData = ggData[, .(val = mean(val,na.rm=TRUE), prop = sum(val>0,na.rm=TRUE) / length(sampleID)),
  #                 by = c("geneName", "grpBy")]

  ggData$val = inv_scale(ggData$val) # do we need this??? we are already normalized...

  # remove zeros

  # Scale if required
  colRange = range(ggData$val)
  if(inpScl){
    ggData[, val:= scale(val), keyby = "geneName"]
    colRange = c(-max(abs(range(ggData$val))), max(abs(range(ggData$val))))
  }

  # hclust row/col if necessary
  ggMat = dcast.data.table(ggData, geneName~grpBy, value.var = "val")
  tmp = ggMat$geneName
  ggMat = as.matrix(ggMat[, -1])
  rownames(ggMat) = tmp
  if(inpRow){
    hcRow = ggdendro::dendro_data(as.dendrogram(hclust(dist(ggMat))))
    ggRow = ggplot2::ggplot() + ggplot2::coord_flip() +
      ggplot2::geom_segment(data = hcRow$segments, ggplot2::aes(x=x,y=y,xend=xend,yend=yend)) +
      ggplot2::scale_y_continuous(breaks = rep(0, uniqueN(ggData$grpBy)),
                         labels = unique(ggData$grpBy), expand = c(0, 0)) +
      ggplot2::scale_x_continuous(breaks = seq_along(hcRow$labels$label),
                         labels = hcRow$labels$label, expand = c(0, 0.5)) +
      sctheme(base_size = sList[inpfsz]) +
      ggplot2::theme(axis.title = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(color="white", angle = 45, hjust = 1))
    ggData$geneName = factor(ggData$geneName, levels = hcRow$labels$label)
  } else {
    ggData$geneName = factor(ggData$geneName, levels = rev(geneList$omic))
  }
  if(inpCol){
    hcCol = ggdendro::dendro_data(as.dendrogram(hclust(dist(t(ggMat)))))
    ggCol = ggplot2::ggplot() +
      ggplot2::geom_segment(data = hcCol$segments, ggplot2::aes(x=x,y=y,xend=xend,yend=yend)) +
      ggplot2::scale_x_continuous(breaks = seq_along(hcCol$labels$label),
                         labels = hcCol$labels$label, expand = c(0.05, 0)) +
      ggplot2::scale_y_continuous(breaks = rep(0, uniqueN(ggData$geneName)),
                         labels = unique(ggData$geneName), expand=c(0,0)) +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      ggplot2::theme(axis.title = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(color = "white"))
    ggData$grpBy = factor(ggData$grpBy, levels = hcCol$labels$label)
  }

  # Actual plot according to plottype
  if(inpPlt == "Bubbleplot"){
    # Bubbleplot
    ggOut = ggplot2::ggplot(ggData, ggplot2::aes(grpBy, geneName, color = val, size = prop)) +
      ggplot2::geom_point() +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      ggplot2::scale_x_discrete(expand = c(0.05, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0.5)) +
      ggplot2::scale_size_continuous("proportion", range = c(0, 8),
                            limits = c(0, 1), breaks = c(0.00,0.25,0.50,0.75,1.00)) +
      ggplot2::scale_color_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) +
      ggplot2::guides(color = ggplot2::guide_colorbar(barwidth = 15)) +
      ggplot2::theme(axis.title = ggplot2::element_blank(), legend.box = "vertical")
  } else {
    # Heatmap
    ggOut = ggplot2::ggplot(ggData, ggplot2::aes(grpBy, geneName, fill = val)) +
      ggplot2::geom_tile() +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      ggplot2::scale_x_discrete(expand = c(0.05, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0.5)) +
      ggplot2::scale_fill_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) +
      ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 15)) +
      ggplot2::theme(axis.title = ggplot2::element_blank())
  }

  # Final tidy
  ggLeg = g_legend(ggOut)
  ggOut = ggOut + ggplot2::theme(legend.position = "none")
  if(!save){
    if(inpRow & inpCol){ggOut =
      gridExtra::grid.arrange(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),
                   layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))
    } else if(inpRow){ggOut =
      gridExtra::grid.arrange(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),
                   layout_matrix = rbind(c(1,3),c(2,NA)))
    } else if(inpCol){ggOut =
      gridExtra::grid.arrange(ggOut, ggLeg, ggCol, heights = c(1,7,2),
                   layout_matrix = rbind(c(3),c(1),c(2)))
    } else {ggOut =
      gridExtra::grid.arrange(ggOut, ggLeg, heights = c(7,2),
                   layout_matrix = rbind(c(1),c(2)))
    }
  } else {
    if(inpRow & inpCol){ggOut =
      gridExtra::arrangeGrob(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),
                  layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))
    } else if(inpRow){ggOut =
      gridExtra::arrangeGrob(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),
                  layout_matrix = rbind(c(1,3),c(2,NA)))
    } else if(inpCol){ggOut =
      gridExtra::arrangeGrob(ggOut, ggLeg, ggCol, heights = c(1,7,2),
                  layout_matrix = rbind(c(3),c(1),c(2)))
    } else {ggOut =
      gridExtra::arrangeGrob(ggOut, ggLeg, heights = c(7,2),
                  layout_matrix = rbind(c(1),c(2)))
    }
  }
  return(ggOut)
}




#  plot_type = c("violin", "boxplot")
pg_violin_box <- function(in_conf, in_meta, in_fact, in_quant,
                          in_grp, in_subsel,
                          in_data, in_omic,
                          plot_type, show_data_points,
                          data_point_sz, font_size){
  data_point_sz = 1.25
  plot_size = "Medium"
  font_size = "Medium"

  if(is.null(in_grp)){in_grp = in_conf$UI[1]}
  # Prepare gg_data
  gg_data = in_meta[, c(in_conf[UI == in_fact]$ID, in_conf[UI == in_grp]$ID),
                    with = FALSE]
  colnames(gg_data) = c("X", "sub")

  # gg_data is three column subset (X,sub,val)
  # Load in either cell meta or gene expr
  if(in_quant %in% in_conf$UI){ #i.e. its in "obs"
    gg_data$val = in_meta[[in_conf[UI == in_quant]$ID]]
  } else {
    # need to unpack from other ad components.
    print("no _val_ in gg_data")
    #subset omics

  }
  # } else {
  #   h5file <- H5File$new(inpH5, mode = "r")
  #   h5data <- h5file[["grp"]][["data"]]
  #   gg_data$val = h5data$read(args = list(in_omic[in_quant], quote(expr=)))
  #   gg_data[val < 0]$val = 0
  #   set.seed(42)
  #   tmpNoise = rnorm(length(gg_data$val)) * diff(range(gg_data$val)) / 1000
  #   gg_data$val = gg_data$val + tmpNoise
  #   h5file$close_all()
  # }


  if(length(in_subsel) != 0 & length(in_subsel) != nlevels(gg_data$sub)){
    gg_data = gg_data[sub %in% in_subsel]
  }

  # Do factoring # IS THIS NESCESSARY?
  gg_col = strsplit(in_conf[UI == in_fact]$fCL, "\\|")[[1]]
  names(gg_col) = levels(gg_data$X)
  gg_lvl = levels(gg_data$X)[levels(gg_data$X) %in% unique(gg_data$X)]
  gg_data$X = factor(gg_data$X, levels = gg_lvl)
  gg_col = gg_col[gg_lvl]

  # Actual ggplot
  if(plot_type == "violin"){
    gg_out = ggplot2::ggplot(gg_data, ggplot2::aes(X, val, fill = X)) + ggplot2::geom_violin(scale = "width")
  } else {
    gg_out = ggplot2::ggplot(gg_data, ggplot2::aes(X, val, fill = X)) + ggplot2::geom_boxplot()
  }
  if(show_data_points){
    gg_out = gg_out + ggplot2::geom_jitter(size = data_point_sz, shape = 16)
  }
  gg_out = gg_out + ggplot2::xlab(in_fact) + ggplot2::ylab(in_quant) +
    sctheme(base_size = sList[font_size], Xang = 45, XjusH = 1) +
    ggplot2::scale_fill_manual("", values = gg_col) +
    ggplot2::theme(legend.position = "none")

  return(gg_out)
}

# Plot theme
sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5){
  oupTheme = ggplot2::theme(
    text =             ggplot2::element_text(size = base_size, family = "Helvetica"),
    panel.background = ggplot2::element_rect(fill = "white", colour = NA),
    axis.line =   ggplot2::element_line(colour = "black"),
    axis.ticks =  ggplot2::element_line(colour = "black", size = base_size / 20),
    axis.title =  ggplot2::element_text(face = "bold"),
    axis.text =   ggplot2::element_text(size = base_size),
    axis.text.x = ggplot2::element_text(angle = Xang, hjust = XjusH),
    legend.position = "bottom",
    legend.key =      ggplot2::element_rect(colour = NA, fill = NA)
  )
  if(!XYval){
    oupTheme = oupTheme + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
  }
  return(oupTheme)
}

### Useful stuff
# Colour palette
cList = list(c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84",
               "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000"),
             c("#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF",
               "#FEE090","#FDAE61","#F46D43","#D73027")[c(1,1:9,9)],
             c("#FDE725","#AADC32","#5DC863","#27AD81","#21908C",
               "#2C728E","#3B528B","#472D7B","#440154"))
names(cList) = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple")

# Panel sizes
pList = c("400px", "600px", "800px")
names(pList) = c("Small", "Medium", "Large")
pList2 = c("500px", "700px", "900px")
names(pList2) = c("Small", "Medium", "Large")
pList3 = c("600px", "800px", "1000px")
names(pList3) = c("Small", "Medium", "Large")
sList = c(18,24,30)
names(sList) = c("Small", "Medium", "Large")
lList = c(5,6,7)
names(lList) = c("Small", "Medium", "Large")

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
    h4("Cell information / gene expression violin plot / box plot"),
    "In this tab, users can visualise the gene expression or continuous cell information ",
    "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).",
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
        radioButtons(ns("RB_plot_type"), "Plot type:",
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
   "In this tab, users can visualise the gene expression patterns of ",
   "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(),
   "The normalised expression are averaged, log-transformed and then plotted.",
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
mod_pg_vis_raw_server <- function(id,rv_in, p){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    # p
    # omics_list = NULL,
    # raw_plot_type = NULL,
    # comp_plot_type = NULL,
    # observ_grp = NULL,
    # observ_subsel = NULL,

    data_point_sz = 1.25
    plot_size = "Medium"
    font_size = "Medium"
    color_schemes = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple")
    color_scheme = color_schemes[2]
    # observe({
    #   req(p$observ_grp,
    #       p$observ_subsel,
    #       rv_in$aux_raw)
    #   browser()
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
    #   if (!is.null(rv_in$aux_raw)){
    #     # could also extract from the name instead of lookkup
    #     #       subs <- strsplit(rv_in$aux_raw, "\\|")[[1]]
    #     # rv_in$ad[[subs[1]]][[subs[2]]]
    #     dat_source = rv_in$config$mat[fIDloc == rv_in$aux_raw]$ID
    #     dat_key = rv_in$config$mat[fIDloc == rv_in$aux_raw]$fID
    #     #dd <- isolate(rv_in$ad[[rv_in$config$mat[fID == rv_in$aux_raw]$ID]][[rv_in$aux_raw]])
    #     if (dat_source == "obs") {
    #       in_data <- isolate(rv_in$ad$obs[[dat_key]])
    #     } else if (dat_source == "var") {
    #       in_data <- isolate(rv_in$ad$var[[dat_key]])
    #     } else {
    #       in_data <- NULL
    #     }
    #     in_quant <- dat_key #(maybe) just aux_raw
    #   }
    #
    # })


    output$plot_box_out <- renderPlot({
      #req(rv_in$ad)
      req(p$omics_list,
          p$observ_grp,
          p$observ_subsel,
          rv_in$aux_raw)

      in_conf <- rv_in$config$meta
      in_meta <- rv_in$meta


      in_fact <- p$observ_grp
      # in_fact <- rv_in$config$meta[[p$observ_grp]]
      #
      # in_fact <- in_meta[[p$observ_grp]]

      #rv_in$ad$obs[[p$observ_grp]] %in% p$observ_subsel

        # could also extract from the name instead of lookkup
        #       subs <- strsplit(rv_in$aux_raw, "\\|")[[1]]
        # rv_in$ad[[subs[1]]][[subs[2]]]
        dat_source = rv_in$config$mat[fIDloc == rv_in$aux_raw]$ID
        dat_key = rv_in$config$mat[fIDloc == rv_in$aux_raw]$fID
        #dd <- isolate(rv_in$ad[[rv_in$config$mat[fID == rv_in$aux_raw]$ID]][[rv_in$aux_raw]])
        if (dat_source == "obs") {
          in_data <- isolate(rv_in$ad$obs[[dat_key]])
        } else if (dat_source == "var") {
          in_data <- isolate(rv_in$ad$var[[dat_key]])
          names(in_data) <- ad$var_names
        } else {
          in_data <- NULL
          print('boxplots only for obs and var ?@!?')
          #TODO:  spawn warning box
          return(NULL)
        }
        in_quant <- dat_key #(maybe) just aux_raw

      pg_violin_box(in_conf, in_meta, in_fact, in_quant,
                    p$observ_grp, p$observ_subsel,
                    in_data, p$omics_list$value, input$RB_plot_type, input$CB_show_data_points,
                    data_point_sz, font_size)

    })

    output$UI_viz_output <- renderUI({
      plotOutput(ns("plot_box_out"), height = pList2[plot_size])
    })

    ### Plots for tab c1
    # output$Van_c1sub1.ui <- renderUI({
    #   sub = strsplit(Van_conf[UI == input$Van_c1sub1]$fID, "\\|")[[1]]
    #   checkboxGroupInput("Van_c1sub2", "Select which cells to show", inline = TRUE,
    #                      choices = sub, selected = sub)
    # })
    # observeEvent(input$Van_c1sub1non, {
    #   sub = strsplit(Van_conf[UI == input$Van_c1sub1]$fID, "\\|")[[1]]
    #   updateCheckboxGroupInput(session, inputId = "Van_c1sub2", label = "Select which cells to show",
    #                            choices = sub, selected = NULL, inline = TRUE)
    # })
    # observeEvent(input$Van_c1sub1all, {
    #   sub = strsplit(Van_conf[UI == input$Van_c1sub1]$fID, "\\|")[[1]]
    #   updateCheckboxGroupInput(session, inputId = "Van_c1sub2", label = "Select which cells to show",
    #                            choices = sub, selected = sub, inline = TRUE)
    # })


    # output$Van_c1oup <- renderPlot({
    #   pg_violin_box(in_conf, in_meta, input$Van_c1inp1, input$Van_c1inp2,
    #                 input$Van_c1sub1, input$Van_c1sub2,
    #                 "Van_gexpr.h5", Van_gene, input$Van_c1typ, input$Van_c1pts,
    #                 data_point_sz, font_size)
    # })
    #
    # output$Van_c1oup.ui <- renderUI({
    #   plotOutput("Van_c1oup", height = pList2[plot_size])
    # })
    # output$Van_c1oup.pdf <- downloadHandler(
    #   filename = function() { paste0("Van_",input$Van_c1typ,"_",input$Van_c1inp1,"_",
    #                                  input$Van_c1inp2,".pdf") },
    #   content = function(file) { ggsave(
    #     file, device = "pdf", height = input$Van_c1oup.h, width = input$Van_c1oup.w, useDingbats = FALSE,
    #     plot = scVioBox(Van_conf, Van_meta, input$Van_c1inp1, input$Van_c1inp2,
    #                     input$Van_c1sub1, input$Van_c1sub2,
    #                     "Van_gexpr.h5", Van_gene, input$Van_c1typ, input$Van_c1pts,
    #                     input$Van_c1siz, input$Van_c1fsz) )
    #   })
    # output$Van_c1oup.png <- downloadHandler(
    #   filename = function() { paste0("Van_",input$Van_c1typ,"_",input$Van_c1inp1,"_",
    #                                  input$Van_c1inp2,".png") },
    #   content = function(file) { ggsave(
    #     file, device = "png", height = input$Van_c1oup.h, width = input$Van_c1oup.w,
    #     plot = scVioBox(Van_conf, Van_meta, input$Van_c1inp1, input$Van_c1inp2,
    #                     input$Van_c1sub1, input$Van_c1sub2,
    #                     "Van_gexpr.h5", Van_gene, input$Van_c1typ, input$Van_c1pts,
    #                     input$Van_c1siz, input$Van_c1fsz) )
    #   })



    ### Plots for tab d1
    # output$Van_d1sub1.ui <- renderUI({
    #   sub = strsplit(Van_conf[UI == input$Van_d1sub1]$fID, "\\|")[[1]]
    #   checkboxGroupInput("Van_d1sub2", "Select which cells to show", inline = TRUE,
    #                      choices = sub, selected = sub)
    # })
    # observeEvent(input$Van_d1sub1non, {
    #   sub = strsplit(Van_conf[UI == input$Van_d1sub1]$fID, "\\|")[[1]]
    #   updateCheckboxGroupInput(session, inputId = "Van_d1sub2", label = "Select which cells to show",
    #                            choices = sub, selected = NULL, inline = TRUE)
    # })
    # observeEvent(input$Van_d1sub1all, {
    #   sub = strsplit(Van_conf[UI == input$Van_d1sub1]$fID, "\\|")[[1]]
    #   updateCheckboxGroupInput(session, inputId = "Van_d1sub2", label = "Select which cells to show",
    #                            choices = sub, selected = sub, inline = TRUE)
    # })
    output$HTML_header <- renderUI({
      print("renderUI HTML header")

      browser()
      gene_list = p$omics_list$value
      if(nrow(gene_list) > 50){
        HTML("More than 50 input genes! Please reduce the gene list!")
      } else {
        oup = paste0(nrow(gene_list[present == TRUE]), " genes OK and will be plotted")
        if(nrow(geneList[present == FALSE]) > 0){
          oup = paste0(oup, "<br/>",
                       nrow(geneList[present == FALSE]), " genes not found (",
                       paste0(geneList[present == FALSE]$omic, collapse = ", "), ")")
        }
        HTML(oup)
      }
    })

    output$plot_heatmap_out <- renderPlot({

      print("renderPlot heatmap")
      req(p$omics_list,
          p$observ_grp,
          p$observ_subsel,
          rv_in$aux_raw)

      in_conf <- rv_in$config$meta
      in_meta <- rv_in$meta


      in_fact <- p$observ_grp



      dat_source = rv_in$config$mat[fIDloc == rv_in$aux_raw]$ID
      dat_key = rv_in$config$mat[fIDloc == rv_in$aux_raw]$fID
      #  probably add a "matrix" column to config...
      # maybe use the data_source to indicate if we want raw or layers... short circuit for now
      in_data <- isolate(rv_in$ad$X)  # do we need to isolate it??

      # if (dat_source == "obs") {
      #   in_data <- isolate(rv_in$ad$obs[[dat_key]])
      # } else if (dat_source == "var") {
      #   in_data <- isolate(rv_in$ad$var[[dat_key]])
      #   names(in_data) <- rv_in$ad$var_names
      # } else {
      #   in_data <- NULL
      #   print('boxplots only for obs and var ?@!?')
      #   #TODO:  spawn warning box
      #   return(NULL)
      # }

      in_quant <- "X" #dat_key #(maybe) just aux_raw

      # these are the "groups" to show on the x axis
      in_group <- in_fact
      #in_subset1 <- p$observ_subsel
      in_subset1 <- p$observ_grp
      in_subset2 <- p$observ_subsel

      # these are the groups to show on thye y axis
      all_omics <- rv_in$ad$var_names
      names(all_omics) <- all_omics

      all_obs <- rv_in$ad$obs_names
      names(all_obs) <- all_obs

      #
      # in_group  input$Van_d1grp
      # in_subset1  p$observ_grp, input$Van_d1sub1 (cell)
      # in_subset2 p$observ_subsel, input$Van_d1sub2 (cell)
      #
      browser()
      scBubbHeat(in_conf, in_meta, p$omics_list$value, in_group, input$RB_heat_plot_type,
                 p$observ_grp, p$observ_subsel, in_data, all_omics,
                 input$CB_scale, input$CB_cluster_rows, input$CB_cluster_cols,
                 color_scheme, plot_size)
    })
    output$UI_heatmap <- renderUI({
      print("renderUI UI heatmap")

      plotOutput(ns("plot_heatmap_out"), height = pList3[plot_size])
    })

    # output$Van_d1oup.pdf <- downloadHandler(
    #   filename = function() { paste0("Van_",input$Van_d1plt,"_",input$Van_d1grp,".pdf") },
    #   content = function(file) { ggsave(
    #     file, device = "pdf", height = input$Van_d1oup.h, width = input$Van_d1oup.w,
    #     plot = scBubbHeat(Van_conf, Van_meta, input$Van_d1inp, input$Van_d1grp, input$Van_d1plt,
    #                       input$Van_d1sub1, input$Van_d1sub2, "Van_gexpr.h5", Van_gene,
    #                       input$Van_d1scl, input$Van_d1row, input$Van_d1col,
    #                       input$Van_d1cols, input$Van_d1fsz, save = TRUE) )
    #   })
    # output$Van_d1oup.png <- downloadHandler(
    #   filename = function() { paste0("Van_",input$Van_d1plt,"_",input$Van_d1grp,".png") },
    #   content = function(file) { ggsave(
    #     file, device = "png", height = input$Van_d1oup.h, width = input$Van_d1oup.w,
    #     plot = scBubbHeat(Van_conf, Van_meta, input$Van_d1inp, input$Van_d1grp, input$Van_d1plt,
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

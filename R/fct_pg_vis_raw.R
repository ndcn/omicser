# Plot gene expression bubbleplot / heatmap
#' pg_bubble_heatmap
#'#' pg_bubble_heatmap
#'
#' @description playground bubbleplot / heatmap function on a ggplot backbone
#'
#' @param in_conf config
#' @param in_meta meta
#' @param in_omics chosen omics
#' @param in_fact x-axis variable
#' @param plot_type bubble or heat
#' @param in_grp  group for subsetting
#' @param in_subsel subselections for grouping
#' @param in_data data_matrix (ie..  X)
#' @param all_omics master list of omics
#' @param in_do_scale Boolean "to scale"
#' @param in_clust_row cluster row, boolean
#' @param in_clust_col cluster columns, boolean
#' @param color_scheme color maps
#' @param plot_size small,med or large
#' @param save prepping for export (or render)
#'
#' @return ggplot object
#' @export
#'
#' @examples  TODO
#' @noRd
pg_bubble_heatmap <- function(in_conf, in_meta, in_omics, in_fact, plot_type,
                              in_grp, in_subsel, in_data, all_omics, in_do_scale, in_clust_row, in_clust_col,
                              color_scheme, plot_size,
                              save = FALSE){
  if(is.null(in_grp)){in_grp = in_conf$UI[1]}
  # Identify genes that are in our dataset
  omic_list = get_omic_list(in_omics, all_omics)
  omic_list = omic_list[present == TRUE]
  shiny::validate(need(nrow(omic_list) <= 50, "More than 50 genes to plot! Please reduce the gene list!"))
  shiny::validate(need(nrow(omic_list) > 1, "Please input at least 2 genes to plot!"))

  # Prepare ggData
  #h5file <- H5File$new(inpH5, mode = "r")
  #h5data <- h5file[["grp"]][["data"]]
  ggData = data.table()
  for(gene_i in omic_list$omic){
    tmp = in_meta[, c("sampleID", in_conf[UI == in_grp]$ID), with = FALSE]
    colnames(tmp) = c("sampleID", "sub")
    tmp$grpBy = in_meta[[in_conf[UI == in_fact]$ID]]
    tmp$geneName = gene_i
    tmp$val = in_data[,all_omics[gene_i]] #h5data$read(args = list(all_omics[gene_i], quote(expr=)))
    ggData = rbindlist(list(ggData, tmp))
  }
  # h5file$close_all()


  if(length(in_subsel) != 0 & length(in_subsel) != nlevels(ggData$sub)){
    ggData = ggData[sub %in% in_subsel]
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
  if(in_do_scale){
    ggData[, val:= scale(val), keyby = "geneName"]
    colRange = c(-max(abs(range(ggData$val))), max(abs(range(ggData$val))))
  }

  # hclust row/col if necessary
  ggMat = dcast.data.table(ggData, geneName~grpBy, value.var = "val")
  tmp = ggMat$geneName
  ggMat = as.matrix(ggMat[, -1])
  rownames(ggMat) = tmp
  if(in_clust_row){
    ####
    ####
    hcRow = ggdendro::dendro_data(as.dendrogram(hclust(dist(ggMat))))
    ggRow = ggplot2::ggplot() + ggplot2::coord_flip() +
      ggplot2::geom_segment(data = hcRow$segments, ggplot2::aes(x=x,y=y,xend=xend,yend=yend)) +
      ggplot2::scale_y_continuous(breaks = rep(0, uniqueN(ggData$grpBy)),
                                  labels = unique(ggData$grpBy), expand = c(0, 0)) +
      ggplot2::scale_x_continuous(breaks = seq_along(hcRow$labels$label),
                                  labels = hcRow$labels$label, expand = c(0, 0.5)) +
      sctheme(base_size = sList[plot_size]) +
      ggplot2::theme(axis.title = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(color="white", angle = 45, hjust = 1))
    ggData$geneName = factor(ggData$geneName, levels = hcRow$labels$label)
  } else {
    ggData$geneName = factor(ggData$geneName, levels = rev(omic_list$omic))
  }
  if(in_clust_col){
    hcCol = ggdendro::dendro_data(as.dendrogram(hclust(dist(t(ggMat)))))
    ggCol = ggplot2::ggplot() +
      ggplot2::geom_segment(data = hcCol$segments, ggplot2::aes(x=x,y=y,xend=xend,yend=yend)) +
      ggplot2::scale_x_continuous(breaks = seq_along(hcCol$labels$label),
                                  labels = hcCol$labels$label, expand = c(0.05, 0)) +
      ggplot2::scale_y_continuous(breaks = rep(0, uniqueN(ggData$geneName)),
                                  labels = unique(ggData$geneName), expand=c(0,0)) +
      sctheme(base_size = sList[plot_size], Xang = 45, XjusH = 1) +
      ggplot2::theme(axis.title = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(color = "white"))
    ggData$grpBy = factor(ggData$grpBy, levels = hcCol$labels$label)
  }

  # Actual plot according to plottype
  if(plot_type == "Bubbleplot"){
    # Bubbleplot
    ggOut = ggplot2::ggplot(ggData, ggplot2::aes(grpBy, geneName, color = val, size = prop)) +
      ggplot2::geom_point() +
      sctheme(base_size = sList[plot_size], Xang = 45, XjusH = 1) +
      ggplot2::scale_x_discrete(expand = c(0.05, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0.5)) +
      ggplot2::scale_size_continuous("proportion", range = c(0, 8),
                                     limits = c(0, 1), breaks = c(0.00,0.25,0.50,0.75,1.00)) +
      ggplot2::scale_color_gradientn("expression", limits = colRange, colours = cList[[color_scheme]]) +
      ggplot2::guides(color = ggplot2::guide_colorbar(barwidth = 15)) +
      ggplot2::theme(axis.title = ggplot2::element_blank(), legend.box = "vertical")
  } else {
    # Heatmap
    ggOut = ggplot2::ggplot(ggData, ggplot2::aes(grpBy, geneName, fill = val)) +
      ggplot2::geom_tile() +
      sctheme(base_size = sList[plot_size], Xang = 45, XjusH = 1) +
      ggplot2::scale_x_discrete(expand = c(0.05, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0.5)) +
      ggplot2::scale_fill_gradientn("expression", limits = colRange, colours = cList[[color_scheme]]) +
      ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 15)) +
      ggplot2::theme(axis.title = ggplot2::element_blank())
  }

  # Final tidy
  ggLeg = g_legend(ggOut)
  ggOut = ggOut + ggplot2::theme(legend.position = "none")
  if(!save){
    if(in_clust_row & in_clust_col){ggOut =
      gridExtra::grid.arrange(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),
                              layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))
    } else if(in_clust_row){ggOut =
      gridExtra::grid.arrange(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),
                              layout_matrix = rbind(c(1,3),c(2,NA)))
    } else if(in_clust_col){ggOut =
      gridExtra::grid.arrange(ggOut, ggLeg, ggCol, heights = c(1,7,2),
                              layout_matrix = rbind(c(3),c(1),c(2)))
    } else {ggOut =
      gridExtra::grid.arrange(ggOut, ggLeg, heights = c(7,2),
                              layout_matrix = rbind(c(1),c(2)))
    }
  } else {
    if(in_clust_row & in_clust_col){ggOut =
      gridExtra::arrangeGrob(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),
                             layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))
    } else if(in_clust_row){ggOut =
      gridExtra::arrangeGrob(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),
                             layout_matrix = rbind(c(1,3),c(2,NA)))
    } else if(in_clust_col){ggOut =
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
#' pg_violin_box
#' pg_violin_box
#'#' pg_violin_box
#'
#' @description playground bubbleplot / heatmap function on a ggplot backbone
#
#' @param in_conf config
#' @param in_meta meta
#' @param in_fact x-axis grouping
#' @param in_quant y-axis quantity to visualize
#' @param in_grp grouping variable
#' @param in_subsel subset of group
#' @param in_data  raw_data (in_quant values)
#' @param in_omic omic for subset (Not implimented)
#' @param plot_type violin or box
#' @param show_data_points boolean
#' @param data_point_sz size - smaller if too many data points is better
#' @param font_size font size
#'
#' @return rendered ggplot
#' @export
#'
#' @examples TBD
pg_violin_box <- function(in_conf, in_meta, in_fact, in_quant,
                          in_grp, in_subsel,
                          in_data, in_omic,
                          plot_type, show_data_points,
                          data_point_sz, font_size){


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
    gg_out = gg_out + ggplot2::geom_jitter(size = data_point_sz, shape = 16, alpha=0.25) #shape="."
  }
  gg_out = gg_out + ggplot2::xlab(in_fact) + ggplot2::ylab(in_quant) +
    sctheme(base_size = sList[font_size], Xang = 45, XjusH = 1) +
    ggplot2::scale_fill_manual("", values = gg_col) +
    ggplot2::theme(legend.position = "none")

  return(gg_out)
}

# Function to extract legend
g_legend <- function(a.gplot){
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

# Get omic list
get_omic_list <- function(omic, full_omics){
  omic_list = data.table(omic = omic,present=TRUE)# unique(trimws(strsplit(inp, ",|;|")[[1]])),

  omic_list[!omic %in% (full_omics)]$present = FALSE
  return(omic_list)
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

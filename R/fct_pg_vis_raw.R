#TODO: clean up unused functions / code here


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
                              in_grp, in_subset,in_subsel,
                              in_data, all_omics,
                              in_do_scale, in_clust_row, in_clust_col,
                              color_scheme, plot_size, save = FALSE){


  if(is.null(in_grp)){in_grp = in_conf$UI[1]}
  # Identify genes that are in our dataset
  omic_list = get_omic_list(in_omics, all_omics)
  omic_list = omic_list[present == TRUE]
  shiny::validate(need(nrow(omic_list) <= 50, "More than 50 genes to plot! Please reduce the gene list!"))
  shiny::validate(need(nrow(omic_list) > 1, "Please input at least 2 genes to plot!"))

  # Prepare gg_data
  #h5file <- H5File$new(inpH5, mode = "r")
  #h5data <- h5file[["grp"]][["hm_dat"]]
  gg_data = data.table()
  for(gene_i in omic_list$omic){
    tmp = in_meta[, c("sample_ID", in_conf[UI == in_subset]$ID), with = FALSE]
    colnames(tmp) = c("sample_ID", "sub")
    tmp$grp_by = in_meta[[in_conf[UI == in_fact]$ID]]
    tmp$omic_nm = gene_i
    tmp$val = in_data[,all_omics[gene_i]] #h5data$read(args = list(all_omics[gene_i], quote(expr=)))
    gg_data = rbindlist(list(gg_data, tmp))
  }
  # h5file$close_all()


  if(length(in_subsel) != 0 & length(in_subsel) != nlevels(gg_data$sub)){
    gg_data = gg_data[sub %in% in_subsel]
  }
  shiny::validate(need(uniqueN(gg_data$grp_by) > 1, "Only 1 group present, unable to plot!"))

  # Aggregate:  exp(mean) + proportion -> log(val)
  # disabling expm1 and log1p because we will be putting "NORMALIZED" quantities in...

  if (max(abs(gg_data$val),na.rm=TRUE) < 700.0 ){
    fwd_scale <- expm1
    inv_scale <- log1p
  } else {
    fwd_scale <- function( input ){ input }
    inv_scale <- function( input ){ input }
  }

  # remove NA
  gg_data <- gg_data[!is.na(val)]
  gg_data$val = fwd_scale(gg_data$val)
  gg_data = gg_data[, .(val = mean(val), prop = sum(val>0) / length(sample_ID)),
                  by = c("omic_nm", "grp_by")]
  # gg_data = gg_data[, .(val = mean(val,na.rm=TRUE), prop = sum(val>0,na.rm=TRUE) / length(sample_ID)),
  #                 by = c("omic_nm", "grp_by")]

  gg_data$val = inv_scale(gg_data$val) # do we need this??? we are already normalized...

  # remove zeros

  # Scale if required
  colRange = range(gg_data$val)
  if(in_do_scale){
    gg_data[, val:= scale(val), keyby = "omic_nm"]
    colRange = c(-max(abs(range(gg_data$val))), max(abs(range(gg_data$val))))
  }

  # hclust row/col if necessary
  gg_mat = dcast.data.table(gg_data, omic_nm~grp_by, value.var = "val")
  tmp = gg_mat$omic_nm
  gg_mat = as.matrix(gg_mat[, -1])
  rownames(gg_mat) = tmp
  if(in_clust_row){
    ####
    ####
    hc_row = ggdendro::dendro_data(as.dendrogram(hclust(dist(gg_mat))))
    gg_row = ggplot2::ggplot() + ggplot2::coord_flip() +
      ggplot2::geom_segment(data = hc_row$segments, ggplot2::aes(x=x,y=y,xend=xend,yend=yend)) +
      ggplot2::scale_y_continuous(breaks = rep(0, uniqueN(gg_data$grp_by)),
                                  labels = unique(gg_data$grp_by), expand = c(0, 0)) +
      ggplot2::scale_x_continuous(breaks = seq_along(hc_row$labels$label),
                                  labels = hc_row$labels$label, expand = c(0, 0.5)) +
      sctheme(base_size = sList[plot_size]) +
      ggplot2::theme(axis.title = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(color="white", angle = 45, hjust = 1))
    gg_data$omic_nm = factor(gg_data$omic_nm, levels = hc_row$labels$label)
  } else {
    gg_data$omic_nm = factor(gg_data$omic_nm, levels = rev(omic_list$omic))
  }
  if(in_clust_col){
    hc_col = ggdendro::dendro_data(as.dendrogram(hclust(dist(t(gg_mat)))))
    gg_col = ggplot2::ggplot() +
      ggplot2::geom_segment(data = hc_col$segments, ggplot2::aes(x=x,y=y,xend=xend,yend=yend)) +
      ggplot2::scale_x_continuous(breaks = seq_along(hc_col$labels$label),
                                  labels = hc_col$labels$label, expand = c(0.05, 0)) +
      ggplot2::scale_y_continuous(breaks = rep(0, uniqueN(gg_data$omic_nm)),
                                  labels = unique(gg_data$omic_nm), expand=c(0,0)) +
      sctheme(base_size = sList[plot_size], Xang = 45, XjusH = 1) +
      ggplot2::theme(axis.title = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(color = "white"))
    gg_data$grp_by = factor(gg_data$grp_by, levels = hc_col$labels$label)
  }

  # Actual plot according to plottype
  if(plot_type == "Bubbleplot"){
    # Bubbleplot
    gg_out = ggplot2::ggplot(gg_data, ggplot2::aes(grp_by, omic_nm, color = val, size = prop)) +
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
    gg_out = ggplot2::ggplot(gg_data, ggplot2::aes(grp_by, omic_nm, fill = val)) +
      ggplot2::geom_tile() +
      sctheme(base_size = sList[plot_size], Xang = 45, XjusH = 1) +
      ggplot2::scale_x_discrete(expand = c(0.05, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0.5)) +
      ggplot2::scale_fill_gradientn("expression", limits = colRange, colours = cList[[color_scheme]]) +
      ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 15)) +
      ggplot2::theme(axis.title = ggplot2::element_blank())
  }

  # Final tidy
  gg_leg = g_legend(gg_out)
  gg_out = gg_out + ggplot2::theme(legend.position = "none")
  if(!save){
    if(in_clust_row & in_clust_col){gg_out =
      gridExtra::grid.arrange(gg_out, gg_leg, gg_col, gg_row, widths = c(7,1), heights = c(1,7,2),
                              layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))
    } else if(in_clust_row){gg_out =
      gridExtra::grid.arrange(gg_out, gg_leg, gg_row, widths = c(7,1), heights = c(7,2),
                              layout_matrix = rbind(c(1,3),c(2,NA)))
    } else if(in_clust_col){gg_out =
      gridExtra::grid.arrange(gg_out, gg_leg, gg_col, heights = c(1,7,2),
                              layout_matrix = rbind(c(3),c(1),c(2)))
    } else {gg_out =
      gridExtra::grid.arrange(gg_out, gg_leg, heights = c(7,2),
                              layout_matrix = rbind(c(1),c(2)))
    }
  } else {
    if(in_clust_row & in_clust_col){gg_out =
      gridExtra::arrangeGrob(gg_out, gg_leg, gg_col, gg_row, widths = c(7,1), heights = c(1,7,2),
                             layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))
    } else if(in_clust_row){gg_out =
      gridExtra::arrangeGrob(gg_out, gg_leg, gg_row, widths = c(7,1), heights = c(7,2),
                             layout_matrix = rbind(c(1,3),c(2,NA)))
    } else if(in_clust_col){gg_out =
      gridExtra::arrangeGrob(gg_out, gg_leg, gg_col, heights = c(1,7,2),
                             layout_matrix = rbind(c(3),c(1),c(2)))
    } else {gg_out =
      gridExtra::arrangeGrob(gg_out, gg_leg, heights = c(7,2),
                             layout_matrix = rbind(c(1),c(2)))
    }
  }
  return(gg_out)
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
                          in_grp, in_subset, in_subsel,
                          in_data, in_omic,
                          plot_type, show_data_points,
                          data_point_sz, font_size){
  #in_meta is either obs or var... depending on which dimension we are visualizing...

  if(is.null(in_subset)) {
    in_subset = in_conf[grp==TRUE]$UI[1]
    }

  # Prepare gg_data
  gg_data <- in_meta[, c(in_conf[ID == in_fact]$ID, in_conf[UI == in_subset]$ID),  with = FALSE]
  #gg_data$X <- rownames(gg_data) #we don't have teh sample ID col anynore
  #
  colnames(gg_data) = c("X", "sub")
  if (is.numeric(in_data)) {  #obs or var
    gg_data$val <- in_data # we already got our vector...
  } else if (endsWith(class(in_data)[[1]] ,"Matrix")) { # there has to be a better way
    # need to append matrix onto
    gg_data$val <- as.data.frame(Matrix::rowMeans(in_data[row.names(in_meta),],na.rm=TRUE))
    #add noise?
      gg_data[val < 0]$val = 0
      set.seed(42)
      tmpNoise = rnorm(length(gg_data$val)) * diff(range(gg_data$val)) / 1000
      gg_data$val = gg_data$val + tmpNoise
  } else {
     print("in_data seems wrong")
   }
# TODO: fix subsetting logic... we can only subset the same grouping for metadata...
#  OR we need to subset the raw-matrix and THEN group...
  if(length(in_subsel) != 0 & length(in_subsel) != nlevels(gg_data$sub)){
      gg_data = gg_data[gg_data$sub %in% in_subsel]
    }


  # map color onto the plot from our config
  gg_col = strsplit(in_conf[UI == in_fact]$fCL, "\\|")[[1]]  #pre-computed colors..
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




# violin_box ==========================
violin_box <- function(vb_dat, x_name,y_name,grp_colors,
                          plot_type,
                          show_data_points,
                          data_point_sz,
                          font_size,
                          notch,
                          grp) {
  #in_meta is either obs or var... depending on which dimension we are visualizing...

  if (grp){
    gg_out = ggplot2::ggplot(vb_dat, ggplot2::aes(x=X, y=val, fill = grp))
  } else {
    gg_out = ggplot2::ggplot(vb_dat, ggplot2::aes(x=X, y=val, fill = X))
  }

  # Actual ggplot
  if(plot_type == "violin"){
    gg_out = gg_out + ggplot2::geom_violin(scale = "width")
  } else {
    gg_out = gg_out + ggplot2::geom_boxplot()
  }

  if(show_data_points){
    gg_out = gg_out + ggplot2::geom_jitter(size = data_point_sz, shape = 16, alpha=0.25) #shape="."
  }

  gg_out = gg_out + ggplot2::xlab(x_name) + ggplot2::ylab(y_name) +
    sctheme(base_size = sList[font_size], Xang = 45, XjusH = 1) +
    ggplot2::theme(legend.position = "none")

  if (!grp) {
    gg_out = gg_out + ggplot2::scale_fill_manual("", values = grp_colors)
  }

  return(gg_out)
}

# bubble_heatmap ==========================
bubble_heatmap <- function(heat_data, x_names, y_names, plot_type,
                              in_do_scale, in_clust_row, in_clust_col,
                              color_scheme, plot_size, grp, save = FALSE){

  hm_dat <- isolate(heat_data)

  # Aggregate:  exp(mean) + proportion -> log(val)
  # disabling expm1 and log1p because we will be putting "NORMALIZED" quantities in...

  if (max(abs(hm_dat$val),na.rm=TRUE) < 700.0 ){
    fwd_scale <- expm1
    inv_scale <- log1p
  } else {
    fwd_scale <- function( input ){ input }
    inv_scale <- function( input ){ input }
  }

# hm_dat:
#   X_ID
#   sub
#   Y_nm
#   value

  # remove NA
  hm_dat <- hm_dat[!is.na(val)]
  hm_dat$val = fwd_scale(hm_dat$val)
  hm_dat = hm_dat[, .(val = mean(val), prop = sum(val>0) / length(X_ID)),
                    by = c("Y_nm", "X_nm")]
  # hm_dat = hm_dat[, .(val = mean(val,na.rm=TRUE), prop = sum(val>0,na.rm=TRUE) / length(sample_ID)),
  #                 by = c("omic", "grp_by")]

  hm_dat$val = inv_scale(hm_dat$val) # do we need this??? we are already normalized...

  # remove zeros

  # Scale if required
  color_range = range(hm_dat$val)
  if(in_do_scale){
    hm_dat[, val:= scale(val), keyby = "Y_nm"]
    color_range = c(-max(abs(range(hm_dat$val))), max(abs(range(hm_dat$val))))
  }

  # hclust row/col if necessary
  gg_mat = dcast.data.table(hm_dat, Y_nm~X_nm, value.var = "val")
  tmp = gg_mat$Y_nm
  gg_mat = as.matrix(gg_mat[, -1])
  rownames(gg_mat) = tmp
  if(in_clust_row){
    ####
    ####
    hc_row = ggdendro::dendro_data(as.dendrogram(hclust(dist(gg_mat))))
    gg_row = ggplot2::ggplot() + ggplot2::coord_flip() +
      ggplot2::geom_segment(data = hc_row$segments, ggplot2::aes(x=x,y=y,xend=xend,yend=yend)) +
      ggplot2::scale_y_continuous(breaks = rep(0, uniqueN(hm_dat$X_nm)),
                                  labels = unique(hm_dat$X_nm), expand = c(0, 0)) +
      ggplot2::scale_x_continuous(breaks = seq_along(hc_row$labels$label),
                                  labels = hc_row$labels$label, expand = c(0, 0.5)) +
      sctheme(base_size = sList[plot_size]) +
      ggplot2::theme(axis.title = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(color="white", angle = 45, hjust = 1))

    hm_dat$Y_nm = factor(hm_dat$Y_nm, levels = hc_row$labels$label)
  } else {

    hm_dat$Y_nm = factor(hm_dat$Y_nm, levels = rev(y_names))

  }
  if(in_clust_col){
    hc_col = ggdendro::dendro_data(as.dendrogram(hclust(dist(t(gg_mat)))))
    gg_color = ggplot2::ggplot() +
      ggplot2::geom_segment(data = hc_col$segments, ggplot2::aes(x=x,y=y,xend=xend,yend=yend)) +
      ggplot2::scale_x_continuous(breaks = seq_along(hc_col$labels$label),
                                  labels = hc_col$labels$label, expand = c(0.05, 0)) +
      ggplot2::scale_y_continuous(breaks = rep(0, uniqueN(hm_dat$omic)),
                                  labels = unique(hm_dat$omic), expand=c(0,0)) +
      sctheme(base_size = sList[plot_size], Xang = 45, XjusH = 1) +
      ggplot2::theme(axis.title = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(color = "white"))
    hm_dat$X_nm = factor(hm_dat$X_nm, levels = hc_col$labels$label)
  }

  # Actual plot according to plottype
  if(plot_type == "Bubbleplot"){
    # Bubbleplot
    gg_out = ggplot2::ggplot(hm_dat, ggplot2::aes(X_nm, Y_nm, color = val, size = prop)) +
      ggplot2::geom_point() +
      sctheme(base_size = sList[plot_size], Xang = 45, XjusH = 1) +
      ggplot2::scale_x_discrete(expand = c(0.05, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0.5)) +
      ggplot2::scale_size_continuous("proportion", range = c(0, 8),
                                     limits = c(0, 1), breaks = c(0.00,0.25,0.50,0.75,1.00)) +
      ggplot2::scale_color_gradientn("expression", limits = color_range, colours = cList[[color_scheme]]) +
      ggplot2::guides(color = ggplot2::guide_colorbar(barwidth = 15)) +
      ggplot2::theme(axis.title = ggplot2::element_blank(), legend.box = "vertical")
  } else {
    # Heatmap
    gg_out = ggplot2::ggplot(hm_dat, ggplot2::aes(X_nm, Y_nm, fill = val)) +
      ggplot2::geom_tile() +
      sctheme(base_size = sList[plot_size], Xang = 45, XjusH = 1) +
      ggplot2::scale_x_discrete(expand = c(0.05, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0.5)) +
      ggplot2::scale_fill_gradientn("expression", limits = color_range, colours = cList[[color_scheme]]) +
      ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 15)) +
      ggplot2::theme(axis.title = ggplot2::element_blank())
  }

  # Final tidy
  gg_leg = g_legend(gg_out)
  gg_out = gg_out + ggplot2::theme(legend.position = "none")
  if(!save){
    if(in_clust_row & in_clust_col){gg_out =
      gridExtra::grid.arrange(gg_out, gg_leg, gg_color, gg_row, widths = c(7,1), heights = c(1,7,2),
                              layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))
    } else if(in_clust_row){gg_out =
      gridExtra::grid.arrange(gg_out, gg_leg, gg_row, widths = c(7,1), heights = c(7,2),
                              layout_matrix = rbind(c(1,3),c(2,NA)))
    } else if(in_clust_col){gg_out =
      gridExtra::grid.arrange(gg_out, gg_leg, gg_color, heights = c(1,7,2),
                              layout_matrix = rbind(c(3),c(1),c(2)))
    } else {gg_out =
      gridExtra::grid.arrange(gg_out, gg_leg, heights = c(7,2),
                              layout_matrix = rbind(c(1),c(2)))
    }
  } else {
    if(in_clust_row & in_clust_col){gg_out =
      gridExtra::arrangeGrob(gg_out, gg_leg, gg_color, gg_row, widths = c(7,1), heights = c(1,7,2),
                             layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))
    } else if(in_clust_row){gg_out =
      gridExtra::arrangeGrob(gg_out, gg_leg, gg_row, widths = c(7,1), heights = c(7,2),
                             layout_matrix = rbind(c(1,3),c(2,NA)))
    } else if(in_clust_col){gg_out =
      gridExtra::arrangeGrob(gg_out, gg_leg, gg_color, heights = c(1,7,2),
                             layout_matrix = rbind(c(3),c(1),c(2)))
    } else {gg_out =
      gridExtra::arrangeGrob(gg_out, gg_leg, heights = c(7,2),
                             layout_matrix = rbind(c(1),c(2)))
    }
  }
  return(gg_out)
}



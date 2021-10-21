#TODO: clean up unused functions / code here



# Function to extract legend
g_legend <- function(a.gplot){
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

# # Get omic list
# get_omic_list <- function(omic, full_omics){
#   omic_list = data.table(omic = omic,present=TRUE)# unique(trimws(strsplit(inp, ",|;|")[[1]])),
#
#   omic_list[!omic %in% (full_omics)]$present = FALSE
#   return(omic_list)
# }


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



# bubble_heatmap ==========================
bubble_heatmap2 <- function(hm_dat, gg_mat, plot_type,
                           in_do_scale, in_clust_row, in_clust_col,
                           color_scheme, plot_size, grp, save = FALSE){




  # Scale if required
  color_range = c(-max(abs(range(gg_mat))), max(abs(range(gg_mat))))

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

# simple_bubble_heatmap ==========================
simple_bubble_heatmap <- function(in_hdata, x_names, y_names, plot_type,
                           in_do_scale, in_clust_row, in_clust_col,
                           color_scheme, plot_size, grp_x, grp_y, save = FALSE){

  hm_dat <- isolate(in_hdata) #isolate doesn't mean anything here...
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

unit_name = "scaled 'X'"
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

  ht <- ComplexHeatmap::Heatmap(gg_mat,
                cluster_rows = in_clust_row,
                cluster_columns = in_clust_col,
                column_split = grp_x,
                row_split = grp_y,
                name = unit_name,
                show_parent_dend_line = FALSE)


  return(ht)

}



make_cx_heatmap = function(in_mat,
                           cluster_rows, row_split, show_row_names, omics_title,
                           cluster_columns, column_split,show_column_names, x_aggregated, x_title,x_grp,x_names,
                           units_label, omics, omics_at, top_annotations, right_annotations){

    if (x_aggregated)
      grp_x <- colnames(in_mat) #make sure its not a factor
      #x_names are what went into x_grp
    else {
      grp_x <- as.character(x_grp) #make sure its not a factor
    }


    #top_annotations and ha below should work together (and x_names)

    if (cluster_columns) {
      column_split = NULL
      ha <- ComplexHeatmap::HeatmapAnnotation(groups = ComplexHeatmap::anno_block(labels = levels(grp_x)))
    } else {
      column_split = grp_x
      ha <- NULL
    }



    #right_annotations and ha2 below should work together (and x_names)
    ha2 <- ComplexHeatmap::rowAnnotation(feats = ComplexHeatmap::anno_mark(at = omics_at, labels = omics))


    ht <- ComplexHeatmap::Heatmap(in_mat,
                                  cluster_rows = cluster_rows,
                                  cluster_columns = cluster_columns,
                                  top_annotation = ha,
                                  show_row_names = show_row_names,
                                  show_column_names = show_column_names,
                                  row_names_side = "right",
                                  row_names_gp = grid::gpar(fontsize = 7),
                                  column_split = column_split,
                                  row_split = row_split,
                                  # border = border,
                                  name = units_label,
                                  column_title = x_title,
                                  row_title = omics_title) + ha2

    # ,
    #                               raster_quality = 1,
    #                               #show_parent_dend_line = TRUE,
    #                               use_raster = FALSE)
    # # ,
    # #                               raster_device = "png")

    message("plot_heatmap_out: FINISHED ComplexHeatmap::Heatmap")

    ht <- ComplexHeatmap::draw(ht, merge_legend = TRUE)
}
#
#
#   env$row_index = which(l)
#
#
#   ht = draw(ht,


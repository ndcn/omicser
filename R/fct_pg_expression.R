#TODO: clean up unused functions / code here


col_unif = list(c("white", "orange"),
                 c("white", "purple"),
                 c("black", "orange"),
                 c("black", "purple"))

col_norm = list(c("green", "white", "red"),
                 c("purple", "white", "orange"),
                 c("blue", "white", "red"),
                 c("orange", "white", "pink")
                 )

col_cats10 = list("Set3","Paired","Pastel1") #Ncols==10,11,12
col_cats9 = list("Set1","Pastel1") #Ncols==9,10
col_cats8 = list("Accent","Dark2","Pastel2","Set2") #1-8

col_cats <- c(col_cats10,col_cats9,col_cats8)

get_my_cols2 <- function(top_annotations){
  # if aggregated don't show top annotations...
  max_levels <- 12

  rpt_cats <- 1 #index to non-repeating colormaps.
  rpt_unif <- 1
  rpt_norm <- 1
  # A- "observations"  (ad_in$obs) pack in the observations.
  top_colors <- list()
  annot_colnms <- colnames(top_annotations)

  for (annot_i in annot_colnms) {
    # Additional pre-processing for categorical metadata
    meta_i <- top_annotations[[annot_i]]
    n_levels <- length(unique(meta_i))
    if (n_levels > 2 & n_levels <= max_levels) {
      if (is.factor(meta_i)) {
        col_i <- structure(brewer.pal(n_levels, col_cats[[rpt_cats]]),
                           names = levels(meta_i))
        rpt_cats <- rpt_cats %% 4 + 1
        top_colors[[annot_i]] <- col_i
      } else if (is.numeric(meta_i)) {
        mx <- max(meta_i)
        mn <- min(meta_i)
        if (mn < 0) {
          mx <- round(max(abs(meta_i)))
          col_i <- circlize::colorRamp2(c(-mx, 0, mx), col_norm[[rpt_norm]])
          rpt_norm <- rpt_norm %% 4 + 1

        } else {
          col_i <- circlize::colorRamp2(c(mn, mx), col_unif[[rpt_unif]])
          rpt_unif <- rpt_unif %% 4 + 1
        }
        top_colors[[annot_i]] <- col_i

      } else { #charachter
        col_i <- structure(brewer.pal(n_levels, col_cats[[rpt_cats]]),
                           names = levels(factor(meta_i)))
        rpt_cats <- rpt_cats %% 4 + 1
        top_colors[[annot_i]] <- col_i

      }

    } else {
      if (is.numeric(meta_i)) {
        mx <- max(meta_i)
        mn <- min(meta_i)
        if (mn < 0) {
          mx <- round(max(abs(meta_i)))
          col_i <- circlize::colorRamp2(c(-mx, 0, mx), col_norm[[rpt_norm]])
          rpt_norm <- rpt_norm %% 4 + 1

        } else {
          col_i <- circlize::colorRamp2(c(mn, mx), col_unif[[rpt_unif]])
          rpt_unif <- rpt_unif %% 4 + 1
        }
        top_colors[[annot_i]] <- col_i
      }
    }
  }

  return(top_colors)
}

make_cx_heatmap = function(in_mat,
                           cluster_samps, samp_grp, samp_grp_nm, samp_title,samp_aggregated,samp_split,
                           cluster_feats, feat_grp, feat_grp_nm, feats_title,
                           units_label,
                           omics, omics_at,
                           top_annotations, right_annotations) {

  # TODO: create top annotations that show the distribution within each group
  if (samp_aggregated) {
      grp_columns <- colnames(in_mat) #make sure its not a factor
      # make a new top_annotation
      # TODO:  make other top_annotaions
      #top_annotations <- top_annotations[,samp_grp_nm, with=FALSE]
      if (length(unique(grp_columns))<=10) {
        top_annotations <-data.frame(as.character(grp_columns))
      } else {
        top_annotations <-data.frame(grp_columns)
      }
      colnames(top_annotations) <- samp_grp_nm

    } else {
      grp_columns <- as.character(samp_grp) #make sure its not a factor
    }

  # if feat_grp is NULL don;t split
  show_row_names <- FALSE
  if (cluster_feats) {
    if (!is.null(feat_grp)){
      #split first then group...
      #hc <- hclust(dist(t(in_mat)))
      k <- length(unique(feat_grp))
      #group = cutree(hc, k = k)
      cluster_feats <- ComplexHeatmap::cluster_within_group(t(in_mat), feat_grp)
      row_split <- k
    } else {
       row_split <- NULL
       cluster_feats <- hclust(dist(in_mat))
    }
   } else {
     cluster_feats <- FALSE
     row_split <- feat_grp
   }

  #samp_names are what went into samp_grp
  show_column_names <- TRUE
  if (!cluster_samps || dim(in_mat)[2]<3) {
    cluster_samps <- FALSE
    if (samp_split){
      column_split <- samp_grp
    } else {
      column_split <- NULL
    }
    #hblock <-  ComplexHeatmap::anno_block(labels = levels(samp_grp))
  } else {
    # When `cluster_columns` is a dendrogram, `column_split` can only be a single
    # number.
    if (samp_split){

      #hc <- hclust(dist(t(in_mat)))
      k <- length(unique(samp_grp))
      #group = cutree(hc, k = k)
      cluster_samps <- ComplexHeatmap::cluster_within_group(in_mat, samp_grp)
      column_split <- k

    } else {
      column_split <- NULL
      cluster_samps <- hclust(dist(t(in_mat)))

    }

    #hblock <- NULL
  }


   # COLOR STUFF
    mx <- round(max(abs(in_mat)))
    mat_colors <- circlize::colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))


    #top_annotations and ha below should work together (and samp_names)

    # if isTruthy(cluster_samps) {
    #   # don't split and add block labels if we are clustering...
    #   column_split = NULL
    #   hblock <- NULL
    # } else {
    #   column_split = grp_columns
    #   hblock <- ComplexHeatmap::HeatmapAnnotation(groups = ComplexHeatmap::anno_block(labels = levels(grp_columns)))
    # }

    if (is.null(top_annotations)) {
      t_annots <-  NULL #ComplexHeatmap::HeatmapAnnotation(group=hblock)
    } else {
      if (samp_aggregated) {
        # make a new top_annotation
        # TODO:
        #top_annotations <- top_annotations[,samp_grp_nm, with=FALSE]
        top_annotations <-data.frame(grp_columns)
        colnames(top_annotations) <- samp_grp_nm
      }
      top_colors <- get_my_cols2(top_annotations)
      t_annots <- ComplexHeatmap::HeatmapAnnotation(df = top_annotations, col = top_colors)
    }

    # annots <- colnames(top_annotations)
    # t_annots <- list()
    # for (annot_i in annots){
    #   # figure out what type (is.char, is.numeric... etc)
    #
    #   curr_annot <- top_annotations[[annot_i]]
    #   if (is.factor(curr_annot)){
    #     an <- ComplexHeatmap::anno_summary()
    #   } else if (is.numeric(curr_annot)) {
    #
    #   } else {
    #     message(paste0("skipping ", annot_i," (unknown type)"))
    #   }
    #
    #   # add to list.
    #   #t_annots <- list()

    #}


    #right_annotations and ha2 below should work together (and samp_names)
    #
    rt_colors <- get_my_cols2(right_annotations)

    r_annots <- ComplexHeatmap::rowAnnotation(df = right_annotations, col=rt_colors,
                                              feats = ComplexHeatmap::anno_mark(at = omics_at, labels = omics) )

    # annots <- colnames(right_annotations)
    # r_annots <- list()
    # for (annot_i in annots){
    #   # figure out what type (is.char, is.numeric... etc)
    #
    #         feats = ComplexHeatmap::anno_mark(at = omics_at, labels = omics)
    #
    #   # create the thing...s
    #   # add to list.
    #
    #   r_annots <- append(r_annots, curr_annot)
    #
    # }
    # # add the names to our list.
    # names(r_annots) <- annots
#
#     # pack into a rowAnnotation
#     h_ra <- ComplexHeatmap::rowAnnotation( r_annots )

    ht <- ComplexHeatmap::Heatmap(in_mat,
                                  col = mat_colors,
                                  cluster_rows = cluster_feats,
                                  cluster_columns = cluster_samps,
                                  top_annotation = t_annots,
                                  show_row_names = show_row_names,
                                  show_column_names = show_column_names,
                                  row_names_side = "right",
                                  row_names_gp = grid::gpar(fontsize = 7),
                                  column_split = column_split,
                                  row_split = row_split,
                                  # border = border,
                                  name = units_label,
                                  column_title = samp_title,
                                  row_title = feats_title,
                                  raster_quality = 2,
                                  #show_parent_dend_line = TRUE,
                                  #use_raster = TRUE,
                                  #raster_by_magick = TRUE
                                  ) + r_annots

    # # ,
    # #                               raster_device = "png")

    message("plot_heatmap_out: FINISHED ComplexHeatmap::Heatmap")
    return(ht)
}



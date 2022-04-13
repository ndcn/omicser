


#' gen_config_table
#'
#' @description function to read the data object and create the config files needed by the UI
#' @param ad_in  anndata object input
#' @param db_name where is the database name
#' @param db_root_path where is out database located
#' @param regenerate regenerate the tables if TRUE, use saved if FALSE
#'
#' @return list containint config table, defaults list and omics vector
#'
#' @noRd
#' @import data.table
#' @import RColorBrewer
#'
gen_config_table <- function(ad_in, db_name, db_root_path, regenerate = FALSE) {
  # load or generate configs..

  config_files <- c(file.path(db_root_path,db_name,"omxr_conf.rds" ),
                    file.path(db_root_path,db_name,"omxr_def.rds" ))

  conf_yml <- file.path(db_root_path,db_name,"db_config.yml" )

  # check if we have it or are forcing
  if (!regenerate) {
    # do we have the files
    if (all(file.exists( config_files ))) {
      # are the saved files newer than the yaml?
      if (difftime(file.info(config_files[1])$ctime, file.info(conf_yml)$ctime, units = "days") < 0) {
        # read the rds files
        omxr_conf <- readRDS(file = config_files[1])
        omxr_def <- readRDS(file = config_files[2])

        return(list(conf  = omxr_conf,
                     def   = omxr_def))
      }
    }
  }


    # should this max levels be according to colors??  i.e. 7, or 9?
    max_levels <- 12 # ceiling for considering somethign a factor


    # Get defaults / last saved...
    #conf_list <- configr::read.config( file.path(db_root_path,db_name,"db_config.yml" ) )
    conf_list <- omicser::get_db_conf(db_name, db_root = db_root_path)

    # depricate db_meta... force using a text file / .Rmd and force other things into the db_config.yml

    # PREPROCESS --------------------------------
    samples <- ad_in$obs_names
    obs_meta <- ad_in$obs

    features <- ad_in$var_names
    var_meta <- ad_in$var

    X_dims <- dim(ad_in$X)
    meta_names <- ad_in$obs_keys()
    # default list of features for subsetting/choosing


    # A observation meta ($obs)  ----------------------------------------
    obs_meta <- data.table::as.data.table(ad_in$obs)
      # val_type <- c("character","numeric","integer")
      # col_type <- c("factor","value","annotation")

    # Make sure that all of our categorical variables are factors...
    # "categorical" is defined as anything with less than `max_levels`
    for (i_meta in colnames(obs_meta) ) {
      #TODO: lapply?
      levels <- sort(unique(obs_meta[[i_meta]]))
      nlevels <- length(levels)
      if (nlevels <= max_levels) {
        obs_meta[[i_meta]] <-  factor(obs_meta[[i_meta]], levels = levels)
        if ( typeof(i_meta) == "double" ) print("Warning less than `max_levels` double type")
      }
    }
    # include everything
    meta_to_include <- colnames(obs_meta) #TODO: make this a parameter
    # Start making config data.table
    omxr_conf <- data.table()
    # A- "observations"  (ad_in$obs) pack in the observations.
    for (i_meta in meta_to_include) {

      levels <- sort(unique(obs_meta[[i_meta]]))
      nlevels <- length(levels)
      if (nlevels <= max_levels) {
        obs_meta[[i_meta]] <-  factor(obs_meta[[i_meta]], levels = levels)
        if ( typeof(i_meta) == "double" ) print("Warning less than `max_levels` double type")
      }

      tmp_conf <- data.table(
        ID = i_meta, UI = i_meta, fID = NA, fUI = NA, fCL = NA, field = "obs",
        default = 0, grp = FALSE,diff_exp = FALSE, dimred = FALSE
      )

      # Additional pre-processing for categorical metadata
      n_levels <- nlevels(obs_meta[[i_meta]])
      if ( n_levels !=0 ) {  # && (i_meta %in% conf_list$group_obs)
        if ( n_levels >= 2 ) { # grouping variable...
          tmp_conf$fID <- paste0(levels(obs_meta[[i_meta]]), collapse = "|")
          tmp_conf$fUI <- tmp_conf$fID
          #tmp_conf$fCL <- paste0(colorRampPalette(brewer.pal(12, "Paired"))(n_levels), collapse = "|")
          #tmp_conf$grp <- TRUE
        } else if ( n_levels == 1 ) {
          tmp_conf$fID <- levels(obs_meta[[i_meta]])
          tmp_conf$fUI <- tmp_conf$fID
          #tmp_conf$fCL <- "black"
          print('we should drop constants and put the value in uns')
        }
      }
      omxr_conf <- rbindlist(list(omxr_conf, tmp_conf))
    }


    # B variable annotations (var-meta; ad_in$var)  ----------------------------------------
    var_meta <- data.table::as.data.table(ad_in$var)
    # val_type <- c("character","numeric","integer")
    # col_type <- c("factor","value","annotation")

    var_to_include <- ad_in$var_keys()# the first one should be "features
    for (i_var in var_to_include) {
      tmp_conf <- data.table(
        ID = i_var, UI = i_var, fID = NA, fUI = NA, fCL = NA, field = "var",
        default = 0, grp = FALSE, diff_exp = FALSE, dimred = FALSE
      )

      var_vect <-ad_in$var[[i_var]]
      var_vect <- factor(var_vect)
      n_levels <- nlevels(var_vect)

      if (n_levels !=0 &
          n_levels <= max_levels ) {
        if (n_levels >= 2 ) {
          #check to see if its in the list
          tmp_conf$fID <- paste0(levels(var_vect), collapse = "|")
          tmp_conf$fUI <- tmp_conf$fID
          #tmp_conf$fCL <- paste0(colorRampPalette(brewer.pal(12, "Paired"))(n_levels), collapse = "|")
          #tmp_conf$grp <- TRUE
        } else if (n_levels == 1) {
          tmp_conf$fID <- levels(var_vect)
          tmp_conf$fUI <- tmp_conf$fID
          #tmp_conf$fCL <- "black"
        }
        #TODO: test for comp measure...
      }
      omxr_conf <- rbindlist(list(omxr_conf, tmp_conf))
    }

    # C differential expression   ----------------------------------------
    diffs <- conf_list$diffs
    # diffs <- list(diff_exp_comps =  levels(factor(diff_exp$group)),
    #               diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
    #               diff_exp_tests =  levels(factor(diff_exp$test_type)))
    de_to_include <- names(diffs) # the first one should be "features
    for (i_diff in de_to_include) {
      tmp_conf <- data.table(
        ID = i_diff, UI = i_diff, fID = NA, fUI = NA, fCL = NA, field = "de",
        default = 0, grp = FALSE, diff_exp = TRUE, dimred = FALSE
      )
      tmp_list <-diffs[[i_diff]]
      tmp_list <- factor(tmp_list)
      n_levels <- nlevels(tmp_list)
      if (n_levels <= max_levels) {
        if (n_levels >= 2) {
          tmp_conf$fID <- paste0(levels(tmp_list), collapse = "|")
          tmp_conf$fUI <- tmp_conf$fID
          #tmp_conf$fCL <- paste0(colorRampPalette(brewer.pal(12, "Paired"))(n_levels), collapse = "|")
          #tmp_conf$grp <- TRUE
        } else if (n_levels == 1) {
          tmp_conf$fID <- levels(tmp_list)
          tmp_conf$fUI <- tmp_conf$fID
          #tmp_conf$fCL <- "black"
        }
        #TODO: test for comp measure...
        omxr_conf <- rbindlist(list(omxr_conf, tmp_conf))
      }
    }


    # TODO:  first pass infer dimreds and layers from anndata table
    # D- dimension reductions (ad_in$varm ) ---------------------
    dimreds <- conf_list$dimreds
    dr_to_include <- names(dimreds) # the first one should be "features
    for (i_dr in dr_to_include) {
      tmp_list <-dimreds[[i_dr]]
      tmp_conf <- data.table(
        ID = i_dr, UI = i_dr, fID = NA, fUI = NA, fCL = NA, field = "dr",
        default = 0, grp = FALSE, diff_exp = FALSE, dimred = TRUE
      )

      # Additional pre-processing for categorical metadata
      tmp_list <- factor(tmp_list)
      n_levels <- nlevels(tmp_list)
      if (n_levels <= max_levels) {
        if (n_levels >= 2) {
          tmp_conf$fID <- paste0(levels(tmp_list), collapse = "|")
          tmp_conf$fUI <- tmp_conf$fID
          #tmp_conf$grp <- TRUE
        } else if (n_levels == 1) {
          tmp_conf$fID <- levels(tmp_list)
          tmp_conf$fUI <- tmp_conf$fID
        }
        #TODO: test for comp measure...
        omxr_conf <- rbindlist(list(omxr_conf, tmp_conf))
      }
    }

    # E- matrix layers (ad_in$layers ) ---------------------
    lr_to_include <- conf_list$layer_values
    lr_names <- if (length(conf_list$layer_names)==length(lr_to_include)) conf_list$layer_names else conf_list$layer_values

    for (i in 1:length(lr_to_include)) {
      i_lr <- lr_to_include[i]
      i_nm <- lr_names[i]

      tmp_conf <- data.table(
        ID = i_lr, UI = i_lr, fID = 0, fUI = i_nm, fCL = NA,field = "layer",
        default = 0, grp = FALSE,  diff_exp = FALSE, dimred = TRUE
      )

      #TODO: test for comp measure...
      omxr_conf <- rbindlist(list(omxr_conf, tmp_conf))
    }


    # DEFAULTS ----------------------------------------
    # obs to subset # default selection is all (if multi) or first (if only 1)
    for (def_i in conf_list$group_obs) {
      omxr_conf[ UI==def_i & field=="obs"]$grp <- TRUE
    }

    i <- 0
    for (def_i in conf_list$default_obs) {
      i <- i + 1
      omxr_conf[ UI==def_i & field=="obs"]$default <- i
    }


    for (def_i in conf_list$group_var) {
      omxr_conf[ UI==def_i & field=="var"]$grp <- TRUE
    }

    i <- 0
    for (def_i in conf_list$default_var){
      i <- i + 1
      omxr_conf[ UI==def_i & field=="var"]$default <- i
    }

    #TODO: add color information here
    #     actually maps for GRP
    #     function for annotations..
    #
    #

    omxr_def <- list()
    #copy things over...
    omxr_def$group_obs <- conf_list$group_obs
    omxr_def$group_var <- conf_list$group_var

    omxr_def$obs_annots <- conf_list$obs_annots
    omxr_def$var_annots <- conf_list$var_annots
    omxr_def$target_features <- conf_list$target_features # first 20
    omxr_def$feature_details <- conf_list$feature_deets
    omxr_def$filter_feature  <- conf_list$filter_feature

    # #meta info
    # meta_info
    #   annotation_database
    #   publication
    #   method
    #   organism
    #   lab
    #   source
    #   title
    #   measurment
    #   pub
    #   url
    #   date

    omxr_def$omic_type <- conf_list$omic_type
    omxr_def$meta_info <- conf_list$meta_info
    # depricate these?
    omxr_def$aggregate_by_default <- conf_list$aggregate_by_default
    omxr_def$organism <- conf_list$organism
    omxr_def$lab <- conf_list$lab
    omxr_def$annotation_database <- conf_list$annotation_database

    # WRITE FILES  ----------------------------------------

    saveRDS(omxr_conf, file = file.path(db_root_path,db_name,"omxr_conf.rds") )
    saveRDS(omxr_def, file = file.path(db_root_path,db_name,"omxr_def.rds") )
    #saveRDS(omics_, file = file.path(db_root_path,db_name,"omxr_omics.rds") )


  out_vals <- list(conf  = omxr_conf,
                   def   = omxr_def)

  return(out_vals)

}



# DEPRICATED COLOR FUNCTIONS BELOW
# #
# col_unif = list(c("white", "orange"),
#                 c("white", "purple"),
#                 c("black", "orange"),
#                 c("black", "purple"))
#
# col_norm = list(c("green", "white", "red"),
#                 c("purple", "white", "orange"),
#                 c("blue", "white", "red"),
#                 c("orange", "white", "pink")
# )
#
# col_cats10 = list("Set3","Paired","Pastel1") #Ncols==10,11,12
# col_cats9 = list("Set1","Pastel1") #Ncols==9,10
# col_cats8 = list("Accent","Dark2","Pastel2","Set2") #1-8
# col_unif = list(c("white", "orange"),
#                 c("white", "purple"),
#                 c("black", "orange"),
#                 c("black", "purple"))
#
# col_norm = list(c("green", "white", "red"),
#                 c("purple", "white", "orange"),
#                 c("blue", "white", "red"),
#                 c("orange", "white", "pink")
# )
#
# col_cats10 = list("Set3","Paired","Pastel1") #Ncols==10,11,12
# col_cats9 = list("Set1","Pastel1") #Ncols==9,10
# col_cats8 = list("Accent","Dark2","Pastel2","Set2") #1-8
#
# get_my_cols <- function(top_annotations){
#   # if aggregated don't show top annotations...
#   max_levels <- 12
#
#   rpt_cats <- 1 #index to non-repeating colormaps.
#   rpt_unif <- 1
#   rpt_norm <- 1
#   # A- "observations"  (ad_in$obs) pack in the observations.
#   top_colors <- list()
#   annot_colnms <- colnames(top_annotations)
#   for (annot_i in annot_colnms) {
#     # Additional pre-processing for categorical metadata
#     meta_i <- top_annotations[[annot_i]]
#     n_levels <- length(unique(meta_i))
#     if (n_levels > 2 & n_levels <= max_levels) {
#       if (is.factor(meta_i)) {
#         col_i <- structure(brewer.pal(n_levels, col_cats[[rpt_cats]]),
#                            names = levels(meta_i))
#         rpt_cats <- rpt_cats %% 4 + 1
#         top_colors[[annot_i]] <- col_i
#       } else if (is.numeric(meta_i)) {
#         mx <- max(meta_i)
#         mn <- min(meta_i)
#         if (mn < 0) {
#           mx <- round(max(abs(meta_i)))
#           col_i <- circlize::colorRamp2(c(-mx, 0, mx), col_norm[[rpt_norm]])
#           rpt_norm <- rpt_norm %% 4 + 1
#
#         } else {
#           col_i <- circlize::colorRamp2(c(mn, mx), col_unif[[rpt_unif]])
#           rpt_unif <- rpt_unif %% 4 + 1
#         }
#         top_colors[[annot_i]] <- col_i
#
#       } else { #charachter
#         col_i <- structure(brewer.pal(n_levels, col_cats[[rpt_cats]]),
#                            names = levels(factor(meta_i)))
#         rpt_cats <- rpt_cats %% 4 + 1
#         top_colors[[annot_i]] <- col_i
#
#       }
#
#     } else {
#       if (is.numeric(meta_i)) {
#         mx <- max(meta_i)
#         mn <- min(meta_i)
#         if (mn < 0) {
#           mx <- round(max(abs(in_mat)))
#           col_i <- circlize::colorRamp2(c(-mx, 0, mx), col_norm[[rpt_norm]])
#           rpt_norm <- rpt_norm %% 4 + 1
#         } else {
#           col_i <- circlize::colorRamp2(c(mn, mx), col_unif[[rpt_unif]])
#           rpt_unif <- rpt_unif %% 4 + 1
#         }
#         top_colors[[annot_i]] <- col_i
#       }
#     }
#   }
#
#   return(top_colors)
# }
#



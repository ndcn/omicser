
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
gen_config_table <- function(ad_in, db_name, db_root_path, regenerate = FALSE) {
  # load or generate configs..

  config_files <- c(file.path(db_root_path,db_name,"omxr_conf.rds" ),
                    file.path(db_root_path,db_name,"omxr_def.rds" ))
  # DISABLED THIS... we can just read from the anndata
  # ,
  #                   file.path(db_root_path,db_name,"omxr_omics.rds" ))

  # check if we have it or are forcing
  if ( all(file.exists( config_files )) && !regenerate) {
    omxr_conf <- readRDS(file = config_files[1])
    omxr_def <- readRDS(file = config_files[2])
  } else {
    max_levels <- 25 # ceiling for considering somethign a factor
    max_levels_ui = 25  # ceiling for UI levels


  # Get defaults / last saved...

    conf_list <- configr::read.config( file.path(db_root_path,db_name,"db_config.yml" ) )
    db_meta <- configr::read.config( file.path(db_root_path,db_name,"db_meta.yml" ) )


    # PREPROCESS --------------------------------
    samples <- ad_in$obs_names
    obs_meta <- ad_in$obs

    features <- ad_in$var_names
    var_meta <- ad_in$var

    X_dims <- dim(ad_in$X)
    meta_names <- ad_in$obs_keys()
    # default list of features for subsetting/choosing


    # A observation meta ($obs)  ----------------------------------------
    obs_meta <- as.data.table(ad_in$obs)
    for (i_meta in colnames(obs_meta)) {
      obs_meta[[i_meta]] = unlist(obs_meta[[i_meta]]) # unlist and refactor
      if (is.character(obs_meta[[i_meta]])) {
        levels = sort(unique(obs_meta[[i_meta]]))
        if (length(levels) < max_levels) {
          obs_meta[[i_meta]] = factor(obs_meta[[i_meta]], levels = levels)
        }
      }
    }

    # include everything
    meta_to_include <- colnames(obs_meta)
    # Start making config data.table
    omxr_conf <- data.table()

    # A- "observations"  (ad_in$obs) pack in the observations.
    for (i_meta in meta_to_include) {
      tmp_conf <- data.table(
        ID = i_meta, UI = i_meta, fID = NA, fUI = NA,
        fCL = NA, Qobs = TRUE, field = "obs",
        default = 0, grp = FALSE, measure= FALSE,
        diff_exp = FALSE, dimred = FALSE
      )
      # Additional pre-processing for categorical metadata
      n_levels <- nlevels(obs_meta[[i_meta]])

      if (n_levels !=0 &
          n_levels <= max_levels &
          n_levels < X_dims[1]*.6 ) {
        if (n_levels >= 2) {
          tmp_conf$fID <- paste0(levels(obs_meta[[i_meta]]), collapse = "|")
          tmp_conf$fUI <- tmp_conf$fID
          tmp_conf$fCL <- paste0(colorRampPalette(brewer.pal(12, "Paired"))(n_levels), collapse = "|")
          tmp_conf$grp <- TRUE
        } else if (n_levels == 1) {
          tmp_conf$fID <- levels(obs_meta[[i_meta]])
          tmp_conf$fUI <- tmp_conf$fID
          tmp_conf$fCL <- "black"
          print('we should drop constants and put the value in uns')
        } #else 0
      } # just append the "blank"
      #TODO: test for comp measure...
      omxr_conf <- rbindlist(list(omxr_conf, tmp_conf))
    }

    # B variable annotations (var-meta; ad_in$var)  ----------------------------------------
    var_to_include <- ad_in$var_keys()# the first one should be "features
    for (i_var in var_to_include) {
      tmp_conf <- data.table(
        ID = i_var, UI = i_var, fID = NA, fUI = NA,
        fCL = NA, Qobs = FALSE, field = "var",
        default = 0, grp = FALSE, measure= FALSE,
        diff_exp = FALSE, dimred = FALSE
      )
      var_vect <-ad_in$var[[i_var]]
      var_vect <- factor(var_vect)
      n_levels <- nlevels(var_vect)

      if (n_levels !=0 &
          n_levels <= max_levels &
          n_levels < length(var_vect)*.6) {
        if (n_levels >= 2 &
            (i_var %in% conf_list$x_var) ) {
          #check to see if its in the list
          tmp_conf$fID <- paste0(levels(var_vect), collapse = "|")
          tmp_conf$fUI <- tmp_conf$fID
          tmp_conf$fCL <- paste0(colorRampPalette(brewer.pal(12, "Paired"))(n_levels), collapse = "|")
          tmp_conf$grp <- TRUE
        } else if (n_levels == 1) {
          tmp_conf$fID <- levels(var_vect)
          tmp_conf$fUI <- tmp_conf$fID
          tmp_conf$fCL <- "black"
        }
        #TODO: test for comp measure...
      }
      omxr_conf <- rbindlist(list(omxr_conf, tmp_conf))
    }

    # C differential expression   ----------------------------------------
    diffs <- conf_list$diffs
    # diffs <- list(diff_exp_groups =  levels(factor(diff_exp$group)),
    #               diff_exp_comp_type =  levels(factor(diff_exp$comp_type)),
    #               diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
    #               diff_exp_tests =  levels(factor(diff_exp$test_type)))
    de_to_include <- names(diffs) # the first one should be "features
    for (i_diff in de_to_include) {
      tmp_conf <- data.table(
        ID = i_diff, UI = i_diff, fID = NA, fUI = NA,
        fCL = NA, Qobs = FALSE, field = "de",
        default = 0, grp = FALSE, measure= FALSE,
        diff_exp = TRUE, dimred = FALSE
      )
      tmp_list <-diffs[[i_diff]]
      tmp_list <- factor(tmp_list)
      n_levels <- nlevels(tmp_list)
      if (n_levels <= max_levels) {
        if (n_levels >= 2) {
          tmp_conf$fID <- paste0(levels(tmp_list), collapse = "|")
          tmp_conf$fUI <- tmp_conf$fID
          tmp_conf$fCL <- paste0(colorRampPalette(brewer.pal(12, "Paired"))(n_levels), collapse = "|")
          tmp_conf$grp <- TRUE
        } else if (n_levels == 1) {
          tmp_conf$fID <- levels(tmp_list)
          tmp_conf$fUI <- tmp_conf$fID
          tmp_conf$fCL <- "black"
        }
        #TODO: test for comp measure...
        omxr_conf <- rbindlist(list(omxr_conf, tmp_conf))
      }
    }


    # D- dimension reductions (ad_in$varm ) ---------------------
    dimreds <- conf_list$dimreds
    dr_to_include <- names(dimreds) # the first one should be "features
    for (i_dr in dr_to_include) {
      tmp_list <-dimreds[[i_dr]]
      tmp_conf <- data.table(
        ID = i_dr, UI = i_dr, fID = NA, fUI = NA,
        fCL = NA, Qobs = FALSE, field = "dr",
        default = 0, grp = FALSE, measure= FALSE,
        diff_exp = FALSE, dimred = TRUE
      )

      # Additional pre-processing for categorical metadata
      tmp_list <- factor(tmp_list)
      n_levels <- nlevels(tmp_list)
      if (n_levels <= max_levels) {
        if (n_levels >= 2) {
          tmp_conf$fID <- paste0(levels(tmp_list), collapse = "|")
          tmp_conf$fUI <- tmp_conf$fID
          tmp_conf$fCL <- paste0(colorRampPalette(brewer.pal(12, "Paired"))(n_levels),
                                 collapse = "|"
          )
          tmp_conf$grp <- TRUE
        } else if (n_levels == 1) {
          tmp_conf$fID <- levels(tmp_list)
          tmp_conf$fUI <- tmp_conf$fID
          tmp_conf$fCL <- "black"
        }
        #TODO: test for comp measure...
        omxr_conf <- rbindlist(list(omxr_conf, tmp_conf))
      }
    }

    # E- matrix layers (ad_in$layers ) ---------------------
    lr_to_include <- conf_list$layers
    for (i_lr in lr_to_include) {
      tmp_conf <- data.table(
        ID = i_lr, UI = i_lr, fID = NA, fUI = NA,
        fCL = NA, Qobs = FALSE, field = "layer",
        default = 0, grp = FALSE, measure= FALSE,
        diff_exp = FALSE, dimred = TRUE
      )
      tmp_conf$fID <- 0
      tmp_conf$fUI <- tmp_conf$fID
      tmp_conf$fCL <- "black"
        #TODO: test for comp measure...
        omxr_conf <- rbindlist(list(omxr_conf, tmp_conf))
    }

    # obs to subset # default selection is all (if multi) or first (if only 1)
    # features
    # vars subset

    # OMIC MAPPING (disabled) ----------------------------------------
    # just make a sorted list
    omic_feature_mapping <- features
    names(omic_feature_mapping) <- features #

    features_ <- seq(X_dims[2])
    names(omic_feature_mapping) <- NULL
    names(features_) <- omic_feature_mapping

    # sort alphabetically and by length
    features_ <- features_[order(names(features_))]
    features_ <- features_[order(nchar(names(features_)))]


    for (def_i in 1:length(conf_list$default_factors)){
      # need errorcheck?
      omxr_conf[ UI==conf_list$default_factors[def_i] ]$default <- def_i
    }


    all_measures <- c(conf_list$y_obs, conf_list$y_var)
    for (def_i in 1:length(all_measures)){
      # need errorcheck?
      omxr_conf[ UI==all_measures[def_i] ]$measure <- TRUE
    }


  # DEFAULTS ----------------------------------------
    omxr_def <- list()
    omxr_def$omics <- conf_list$target_omics # first 20

    omxr_def$obs_x <- omxr_conf[default == 1]$UI #

    omxr_def$obs_y <- omxr_conf[measure == TRUE & field == "obs"]$UI[1] #

    # in case we don't have more than 1 default we need to check...
    #
    if (length(conf_list$default_factors)>1){
      def_n <- 2
    } else {
      def_n <- 1
    }

    omxr_def$obs_subset <- omxr_conf[default == def_n]$UI #
    omxr_def$color_grp <- omxr_conf[default == def_n]$UI # color by...


    omxr_def$var_x <- omxr_conf[grp==TRUE & field=="var"]$UI[1] #
    omxr_def$var_y <- omxr_conf[measure == TRUE & field == "var"]$UI[1] #
    omxr_def$var_subset <- omxr_conf[grp==TRUE & field=="var"]$UI[1]  # Us

    # TODO: add new parameter defaults here.
    omxr_def$scale <- NA # raw / area / z
    omxr_def$trans <- NA # log10 or no
    omxr_def$p_val <- NA # corrected?

    omxr_def$comp_name <- diffs$diff_exp_obs_name[1]  # x vs y / x vs rest,
    omxr_def$comp_type <- diffs$diff_exp_comp_type[1]  # x vs y / x vs rest,
    omxr_def$comp_fact <- diffs$diff_exp_comps[1] # choice
    omxr_def$test <- diffs$diff_exp_tests[1]  # t-test, xx

    omxr_def$db_meta <- db_meta
    omxr_def$omic_details <- conf_list$omic_details

    # WRITE FILES  ----------------------------------------

    saveRDS(omxr_conf, file = file.path(db_root_path,db_name,"omxr_conf.rds") )
    saveRDS(omxr_def, file = file.path(db_root_path,db_name,"omxr_def.rds") )
    #saveRDS(omics_, file = file.path(db_root_path,db_name,"omxr_omics.rds") )

  }

  out_vals <- list(conf  = omxr_conf,
                   def   = omxr_def)
                  #    ,
                   #omics = omics_)
  return(out_vals)

}





#' Title
#'
#' @param ad_in
#' @param measures
#' @param diffs
#' @param dimreds
#' @param default_factors
#' @param db_prefix
#' @param db_dir
#'
#' @return
#' @export
#'
#' @examples
#' @import data.table RColorBrewer
create_config_table <- function(ad_in,
                                measures,
                                diffs,
                                dimreds,
                                default_factors,
                                db_prefix = "test1",
                                db_dir = "data-raw") {
  max_levels <- 50 # ceiling for considering somethign a factor
  legend_cols = 4
  max_levels_ui = 50  # ceiling for UI levels

  conf_list <- configr::read.config( file.path("data-raw",DB_NAME,"config.yml" ) )


  # PREPROCESS --------------------------------
  samples <- ad_in$obs_names
  obs_meta <- ad_in$obs

  omics <- ad_in$var_names
  var_meta <- ad_in$var

  # # just in case...
  # if (!is.null(obs_meta$sampleID)){ #if we already have sampleID
  #   # call sampleID sampleIDog
  #   obs_meta <- dplyr::rename(obs_meta,sampleID_0=sampleID)
  # }
  X_dims <- dim(ad_in$X)


  meta_names <- ad_in$obs_keys()







  def_omics = omics[1:20]  # first 20

  obs_meta = data.table(sampleID = samples)  # redundant but makes naming explicit..
  obs_meta = cbind(obs_meta,ad$obs)
  colnames(obs_meta) = c("sampleID", meta_names)

  for (i_meta in colnames(obs_meta)[-1]) {
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

  # A- "observations"  (ad$obs) pack in the observations.
  for (i_meta in meta_to_include) {
    tmp_conf <- data.table(
      ID = i_meta, UI = i_meta, fID = NA, fUI = NA,
      fCL = NA, fRow = NA, Qobs = TRUE, field = "obs",
      default = 0, grp = FALSE, measure= FALSE,
      diff_exp = FALSE, dimred = FALSE
    )
    # Additional pre-processing for categorical metadata
    n_levels <- nlevels(obs_meta[[i_meta]])
    print(levels(obs_meta[[i_meta]]))
    if (n_levels <= max_levels |
        n_levels < X_dims[1]*.6) {
      if (n_levels >= 2) {
        tmp_conf$fID <- paste0(levels(obs_meta[[i_meta]]), collapse = "|")
        tmp_conf$fUI <- tmp_conf$fID
        tmp_conf$fCL <- paste0(colorRampPalette(brewer.pal(12, "Paired"))(n_levels), collapse = "|")
        tmp_conf$fRow <- ceiling(n_levels / legend_cols)
        tmp_conf$grp <- TRUE
      } else if (n_levels == 1) {
        tmp_conf$fID <- levels(obs_meta[[i_meta]])
        tmp_conf$fUI <- tmp_conf$fID
        tmp_conf$fCL <- "black"
        tmp_conf$fRow <- 1
      }
      #TODO: test for comp measure...
      omxr_conf <- rbindlist(list(omxr_conf, tmp_conf))
    }
  }

  # B- "variable" annotations  (ad$var)
  var_to_include <- ad_in$var_keys()[-1] # the first one should be "omics
  for (i_var in var_to_include) {
    tmp_conf <- data.table(
      ID = i_var, UI = i_var, fID = NA, fUI = NA,
      fCL = NA, fRow = NA, Qobs = FALSE, field = "var",
      default = 0, grp = FALSE, measure= FALSE,
      diff_exp = FALSE, dimred = FALSE
    )
    var_vect <-ad_in$var[[i_var]]
    var_vect <- factor(var_vect)
    n_levels <- nlevels(var_vect)
    if (n_levels <= max_levels |
        n_levels < length(var_vect)*.6) {
      if (n_levels >= 2) {
        tmp_conf$fID <- paste0(levels(var_vect), collapse = "|")
        tmp_conf$fUI <- tmp_conf$fID
        tmp_conf$fCL <- paste0(colorRampPalette(brewer.pal(12, "Paired"))(n_levels), collapse = "|")
        tmp_conf$fRow <- ceiling(n_levels / legend_cols)
        tmp_conf$grp <- TRUE
      } else if (n_levels == 1) {
        tmp_conf$fID <- levels(var_vect)
        tmp_conf$fUI <- tmp_conf$fID
        tmp_conf$fCL <- "black"
        tmp_conf$fRow <- 1
      }
      #TODO: test for comp measure...
      omxr_conf <- rbindlist(list(omxr_conf, tmp_conf))
    }
  }

  # diffs <- list(diff_exp_groups =  levels(factor(diff_exp$group)),
  #               diff_exp_comp_type =  levels(factor(diff_exp$comp_type)),
  #               diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
  #               diff_exp_tests =  levels(factor(diff_exp$test_type)))
  # C- "diff_exp" annotations  (ad$var)
  de_to_include <- names(diffs) # the first one should be "omics
  for (i_diff in de_to_include) {
    tmp_conf <- data.table(
      ID = i_diff, UI = i_diff, fID = NA, fUI = NA,
      fCL = NA, fRow = NA, Qobs = FALSE, field = "de",
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
        tmp_conf$fRow <- ceiling(n_levels / legend_cols)
        tmp_conf$grp <- TRUE
      } else if (n_levels == 1) {
        tmp_conf$fID <- levels(tmp_list)
        tmp_conf$fUI <- tmp_conf$fID
        tmp_conf$fCL <- "black"
        tmp_conf$fRow <- 1
      }
      #TODO: test for comp measure...
      omxr_conf <- rbindlist(list(omxr_conf, tmp_conf))
    }
  }

  # D- dimension reductions (ad$varm )
  dr_to_include <- names(dimreds) # the first one should be "omics
  for (i_dr in dr_to_include) {
    tmp_list <-dimreds[[i_dr]]
    tmp_conf <- data.table(
      ID = i_dr, UI = i_dr, fID = NA, fUI = NA,
      fCL = NA, fRow = NA, Qobs = FALSE, field = "dr",
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
        tmp_conf$fRow <- ceiling(n_levels / legend_cols)
        tmp_conf$grp <- TRUE

      } else if (n_levels == 1) {
        tmp_conf$fID <- levels(tmp_list)
        tmp_conf$fUI <- tmp_conf$fID
        tmp_conf$fCL <- "black"
        tmp_conf$fRow <- 1
      }
      #TODO: test for comp measure...
      omxr_conf <- rbindlist(list(omxr_conf, tmp_conf))
    }
  }

  # obs to subset # default selection is all (if multi) or first (if only 1)
  # omics
  # vars subset


  # # force sampleID to match X column names
  # if (!isTRUE(all.equal(obs_meta$sampleID, X_colnm))) {
  #   obs_meta$sampleID <- factor(obs_meta$sampleID, levels = X_colnm)
  #   obs_meta <- obs_meta[order(sampleID)]
  #   obs_meta$sampleID <- as.character(obs_meta$sampleID)
  # }

  omic_mapping <- omics
  names(omic_mapping) <- omics #

  omics_ <- seq(X_dims[2])
  names(omic_mapping) <- NULL
  names(omics_) <- omic_mapping

  # sort alphabetically and by length
  omics_ <- omics_[order(names(omics_))]
  omics_ <- omics_[order(nchar(names(omics_)))]


  for (def_i in 1:length(default_factors)){
    # need errorcheck?
    omxr_conf[ UI==default_factors[def_i] ]$default <- def_i
  }

  all_measures <- c(measures$obs, measures$var)
  for (def_i in 1:length(all_measures)){
    # need errorcheck?
    omxr_conf[ UI==all_measures[def_i] ]$measure <- TRUE
  }


  # DEFAULTS ----------------------------------------
  omxr_def <- list()
  omxr_def$omics <- def_omics #

  omxr_def$obs_x <- omxr_conf[default == 1]$UI #

  omxr_def$obs_y <- omxr_conf[measure == TRUE & field == "obs"]$UI[1] #
  omxr_def$obs_subset <- omxr_conf[default == 2]$UI #

  omxr_def$var_x <- omxr_conf[grp==TRUE & field=="var"]$UI[1] #
  omxr_def$var_y <- omxr_conf[measure == TRUE & field == "var"]$UI[1] #
  omxr_def$var_subset <- omxr_conf[grp==TRUE & field=="var"]$UI[1]  # Us

  omxr_def$color_grp <- omxr_conf[default == 2]$UI # color by...

  omxr_def$scale <- NA # raw / area / z
  omxr_def$trans <- NA # log10 or no
  omxr_def$p_val <- NA # corrected?

  omxr_def$comp_name <- diffs$diff_exp_obs_name[1]  # x vs y / x vs rest,
  omxr_def$comp_type <- diffs$diff_exp_comp_type[1]  # x vs y / x vs rest,
  omxr_def$comp_fact <- diffs$diff_exp_groups[1] # choice
  omxr_def$test <- diffs$diff_exp_tests[1]  # t-test, xx

  #
  #

  saveRDS(omxr_conf, file = paste0(db_dir, "/", db_prefix, "_conf.rds"))
  saveRDS(omxr_def, file = paste0(db_dir, "/", db_prefix, "_def.rds"))

  # make these??
  saveRDS(obs_meta, file = paste0(db_dir, "/", db_prefix, "_meta.rds"))
  saveRDS(omics_, file = paste0(db_dir, "/", db_prefix, "_omics.rds"))

  out_vals <- list(conf=omxr_conf,def=omxr_def)
  return(out_vals)

}

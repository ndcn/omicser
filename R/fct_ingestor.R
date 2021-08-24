#' ingestor
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#'
#' @noRd

omxr_pack_anndata <- function (data_in){

  #tools::file_path_sans_ext(data_in)
  if ( class(data_in)[1] == "list" ) {
    # multple files in c("object", "data_mat","obs_meta","var_annot","omics","sample_ID",etc")
    # data_mat - data matrix.  need to assert that this matrix object has omics and sample+ID names

    if (dim(data_in$data_mat)[2] != dim(data_in$var_annot)[1]) {
      X <-Matrix::t(data_in$data_mat)
    } else {
      X <- data_in$data_mat
    }
    if (is.null(dimnames(X)[1])) { rownames(X) <- data_in$sample_ID }
    if (is.null(dimnames(X)[2])) { colnames(X) <- data_in$omics  }

    if(is.null(rownames(obs)) ){
      rownames(obs) <- data_in$sample_ID
    }
    # obs_meta - ensure that we are factored and sample_ID is first column
    obs <- data_in$obs_meta
    id_col <- colnames(obs)[1]
    if (all(obs[[id_col]] == data_in$sample_ID)) {
      obs <- obs %>% dplyr::rename(sample_ID=all_of(id_col))
    } else {
      obs <- obs %>% dplyr::mutate(sample_ID=data_in$sample_ID)
    }

    if(is.null(rownames(obs)) ){
      rownames(obs) <- data_in$sample_ID
    }

    # var_annot - ensure that "omics" is first column
    var_ <- data_in$var_annot
    omics_col <- colnames(var_)[1]

    if ( all(var_[[omics_col]] == data_in$omics) ) {
      var_ <- var_ %>% dplyr::rename(omics_name=all_of(omics_col))
    } else {
      var_ <- var_ %>% dplyr::mutate(omics_name=data_in$omics)
    }
    if(is.null(rownames(var_)) ){
      rownames(var_) <- data_in$omics
    }

    # etc goes into an uns entry
    #
    if ( class(data_in$uns)=="list" ){
      uns <- data_in$uns
    } else {
      uns <- list(etc=data_in$uns)
    }

    ad <- anndata::AnnData(
      X = X,
      obs = obs,
      var = var_,
      uns = uns
    )

  } else if (tolower(tools::file_ext(data_in)) == "rds") {

    data_in = readRDS( file = data_in )


    if(class(data_in)[1] == "Seurat"){
      # how stereotyped is this pattern?  check for Oscar...
      ad <- sceasy::convertFormat(data_in, from="seurat", to="anndata",
                                  outFile = NULL,
                                  assay = 'SCT',
                                  main_layer = 'data',
                                  transfer_layers = c('data', 'counts', 'scale.data')
      )

      raw <- sceasy::convertFormat(data_in, from="seurat", to="anndata",
                                   outFile = NULL,
                                   assay = 'RNA',
                                   main_layer = 'counts',
                                   transfer_layers = NULL)


      ad$raw <- raw


    } else if (class(data_in)[1] == "SingleCellExperiment") {
      print("SingleCellExperiment not enabled")
      ad <- NULL

    } else if ("data.frame" %in% class(data_in)) {
      # could _everything be in a dataframe???
      # yes... lipidomic... strip off first two columns?
      print("enable this for lipidomics? not enabled")
      ad <- NULL

    }

  } else if (tolower(tools::file_ext(data_in)) == "h5ad") {
    ad <- anndata::read_h5ad(data_in)


  } else if (tolower(tools::file_ext(data_in)) == "loom"){
    print("loom loading not enabled")
    ad <- NULL
  }


  return(ad)
}



#' ingestor
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#'
#' @noRd
setup_database <- function(database_name, data_in, db_meta , re_pack=TRUE){
  #LOAD & PACK into ANNDATA
  ##if data_in contains filenames they must be the full path (i.e. RAW_DIR inlcuded)

  DB_NAME <- database_name


  require(reticulate)
  reticulate::use_condaenv(required = TRUE, condaenv = 'omxr')
  require(anndata)

  DB_DIR = file.path("data-raw",DB_NAME)
  if (!dir.exists(DB_DIR)) {
    dir.create(DB_DIR)
  }

  # # db_meta e.g.
  #   organism <- ""
  #   lab <- "Menon"
  #   annotation_database <- "NA"
  # check to see what level of data we were given in our data_list
  #
  if (length(data_in) & names(data_in[1])=="object") {
    data_in <- unlist(data_in)
  }

  # check if we have it or are forcing
  if ( !file.exists(paste0(DB_DIR,"/core_data.h5ad"))
       | re_pack ) {
    #create it
    # sub-functions to deal with what kind of data we have...
    ad <- omxr_pack_anndata(data_in)

  } else {
    #load it
    ad <- anndata::read_h5ad(file.path(DB_DIR,"core_data.h5ad"))
  }

  return(ad)
}




require("data.table")
require("RColorBrewer")

gen_config_table <- function(ad_in, ds_name) {
  max_levels <- 50 # ceiling for considering somethign a factor
  legend_cols = 4
  max_levels_ui = 50  # ceiling for UI levels

  conf_list <- configr::read.config( file.path("data-raw",ds_name,"config.yml" ) )


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


  # default list of omics for subsetting/choosing
  def_omics = omics[1:20]  # first 20

  #obs_meta = data.table(sampleID = samples)  # redundant but makes naming explicit..
  #obs_meta = cbind(obs_meta,ad$obs)
  #colnames(obs_meta) = c("sampleID", meta_names)
  # for (i_meta in colnames(obs_meta)[-1]) {


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

    if (n_levels !=0 &
        n_levels <= max_levels &
        n_levels < X_dims[1]*.6 ) {
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
        print('we should drop constants and put the value in uns')
      } #else 0
    } # just append the "blank"
    #TODO: test for comp measure...
    omxr_conf <- rbindlist(list(omxr_conf, tmp_conf))
  }

  # B- "variable" annotations  (ad$var)
  var_to_include <- ad_in$var_keys()# the first one should be "omics
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

    if (n_levels !=0 &
        n_levels <= max_levels &
        n_levels < length(var_vect)*.6) {
      if (n_levels >= 2 &
          (i_var %in% conf_list$x_var) ) {
        #check to see if its in the list
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

  diffs <- conf_list$diffs
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


  dimreds <- conf_list$dimreds
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


  for (def_i in 1:length(conf_list$default_factors)){
    # need errorcheck?
    omxr_conf[ UI==default_factors[def_i] ]$default <- def_i
  }

  all_measures <- c(conf_list$y_obs, conf_list$y_var)
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
  omxr_def$comp_fact <- diffs$diff_exp_comps[1] # choice
  omxr_def$test <- diffs$diff_exp_tests[1]  # t-test, xx

  omxr_def$db_meta <- conf_list$db_meta

  out_vals <- list(conf  = omxr_conf,
                   def   = omxr_def,
                   omics = omics_)
  return(out_vals)

}


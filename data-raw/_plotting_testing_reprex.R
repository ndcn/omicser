require(reticulate)
reticulate::use_condaenv(required = TRUE, condaenv = 'omxr')
require(anndata)

source('data-raw/compute_de_table.R')

source('data-raw/create_config_table.R')

source('R/fct_ingestor.R')
# gen_config_table

DB_NAME <- "vilas_microglia_sceasy"
# ad <- read_h5ad(filename=file.path("data-raw",DB_NAME,"normalized_data_with_de.h5ad"))
# diff_exp <- readRDS( file = file.path("data-raw",DB_NAME, "diff_expr_table.rds"))
#
# conf_list_out <- configr::read.config( file.path("data-raw",DB_NAME,"config.yml" ) )


dataset_names <- c(
  "Domenico DIA" = "domenico_stem_cell",
  "Vilas Microglia" = "vilas_microglia",
  "Vilas Microglia (seu)" = "vilas_microglia_seu",
  "Vilas Microglia (sceasy)" = "vilas_microglia_sceasy",
  "Yassene Lipid Concentrations" ="yassene_A_conc",
  "Yassene Lipid Compositions" ="yassene_A_compos",
  "Oscar Microglia" ="oscar_microglia"
)

dataset_type <- c(
  "Transcriptomics" = "transcript",
  "Proteomics" = "prote",
  "Lipidomics" = "lipid",
  "Metabolomics" = "metabol",
  "Other" = "X-"
)

# INGESTOR --------
to_return <- list()

ds_name <- DB_NAME
ds_label <- names(dataset_names[match(ds_name,dataset_names)])

  ad <- anndata::read_h5ad(filename=paste0("data-raw/",ds_name,"/omxr_data.h5ad"))
  diff_exp = readRDS(file = paste0("data-raw/",ds_name,"/diff_expr_table.rds"))

  conf_def <- gen_config_table(ad, ds_name)


  omics <- ad$var_names
  names(omics) <- omics
  # omicmeta <- ad$obs
  #omics_sorted <- conf_def$omics


  to_return$de <- diff_exp

  to_return$database_name <- ds_label
  to_return$omics_type <- to_return$db_meta$omics_type
  to_return$ad <- ad
  to_return$omics <- omics
  # to_return$meta <- omicmeta  # this might be too redundant

  to_return$config <- conf_def$conf
  to_return$default <- conf_def$def

  to_return$db_meta$name<-ds_name
  to_return$db_meta$omics_type<-conf_def$def$db_meta$omic_type
  to_return$db_meta$measurment<-conf_def$def$db_meta$measurement
  to_return$db_meta$organism<-conf_def$def$db_meta$organizm
  to_return$db_meta$publication<-conf_def$def$db_meta$pub
  to_return$db_meta$etc<-conf_def$def$db_meta$annotation_database



  rv_in <- to_return

# SIDE SELECTOR ----------------
  d_source = "obs"

  out_params <- list()
  out_params$data_source <- d_source #rv_in$config[UI == d_source]$field
  out_params$plot_x <-  to_return$default$obs_x[1]
  out_params$plot_y <-  to_return$default$obs_y[1]

  out_params$observ_subset <-  rv_in$default$obs_x # rv_in$config[grp == TRUE & field=="obs"]$UI
  out_params$observ_subsel <- strsplit(rv_in$config[UI == out_params$observ_subset]$fID, "\\|")[[1]]

  # # Disabled
  # out_params$feat_subset = input$SI_var_subset
  # out_params$feat_subsel = input$CB_var_subsel
  agg_obs <- rv_in$config[grp == TRUE & field=="obs"]$UI # subset for aggregating?
  agg_var <- rv_in$config[grp == TRUE & field=="var"]$UI # subset for aggregating
  group_obs <- rv_in$config[grp == TRUE & field=="obs"]$UI # subset for aggregating?
  group_var <- rv_in$config[grp == TRUE & field=="var"]$UI # subset for aggregating

  # group (plotting)
  out_params$feat_group_by <- group_var[1] #input$SI_group_var
  out_params$observ_group_by <- group_obs[1] # input$SI_group_obs

  # aggregate collapsing data matrix
  #         # group (plotting)
  out_params$feat_agg <- agg_var #input$SI_agg_obs
  out_params$observ_agg <- agg_obs# input$SI_agg_var


  out_params$group_action <- "group by" #input$RB_agg_or_grp

  def_omics <- rv_in$default$omics
  omics_list <- list(value=def_omics,
                     viz_now=FALSE)

  out_params[["omics_list"]] <- omics_list  # value & viz_now


  # X
  out_params$data_source <- "X" #rv_in$config[UI == to_return$default$obs_x[1]]$field
  out_params$plot_x <-  to_return$default$obs_x[1]
  out_params$plot_y <-  to_return$default$obs_y[1]
  out_params$group_action <- "none" #input$RB_agg_or_grp



p <- out_params
# PLAYGROUND  ===========================



############# playground reactives
  ret_vals <- list(
    x = NULL,
    y = NULL,
    type = NULL,
    dat = NULL
  )

  # set feat_subset to NA (not enabled) default to "omics list"
  p$feat_subset <- "highly_variable"
  p$feat_subsel <- levels(factor(rv_in$ad$var[[ p$feat_subset ]]))
  p$feat_subsel  <- p$feat_subsel[2] #just choose "highly_variable == TRUE
  p$feat_subset <- NA
  p$feat_subsel <- NA

  p$omics_list$viz_now <- FALSE  # reset it here...

  in_conf <- rv_in$config
  in_meta <- as.data.table(rv_in$ad$obs) # can probably just acces ad$obs directly since we don' tneed it to be a data_table?
  row.names(in_meta) <- rv_in$ad$obs_names

  dat_source <- p$data_source
  #in_data <- isolate(rv_in$ad[[dat_source]][[dat_key]])
  dat_source <- "X"


  # HEATMAP DATA ======================
  # ---------- : prep data_matrix for heat_map/bubble
  # TODO: make a oaram in side_selector
  x_is <- "obs" #haven't implimented visualizing by var yet...
  X_fact <- p$plot_x #in_fact <- p$observ_grpA
  dat_loc <- p$plot_y
  dat_loc ="X"

  if (dat_loc=="X") {
    X_data <- (rv_in$ad$X) #isolate
  } else if (dat_loc == "raw") {
    X_data <- (rv_in$ad$raw$X)#isolate
  } else if (dat_loc == "layers") { #must be a layer
    if (dat_loc %in% rv_in$ad$layers_keys()) {
      X_data <- (rv_in$ad$layers[[dat_loc]])#isolate
    } else {
      print("data not found")
      return(ret_vals)
    }
  }
  # makesure i have rownames right for X_data?
  in_subset <- p$observ_subset
  in_subsel <- p$observ_subsel[1:2]
  # subset by
  group_by <- list( x = p$observ_group_by,
                    y = p$feat_group_by) # NOTE: could be individual "omics" from selector


  ## AGGREGATE along x-axis
  if (x_is == "var"){
    # transpose data
    # subset obs (samples)
    if (!is.na( in_subset ) ) {
      dat_js_set <- in_meta[[ in_subset ]]
      X_data <-  X_data[dat_js_set %in% in_subsel ,]
    } else {
      dat_js_set <- rv_in$ad$obs_names
    }


    dat_js <- p$observ_subsel #if NA or NULL check? (choose all...)
    dat_j_grp <-p$observ_subset

    X_data <- Matrix::t(X_data)
    # switch these around
    in_subset <- p$feat_subset
    in_subsel <- p$feat_subsel
    group_by <- list( y = p$observ_group_by,
                      x = p$feat_group_by)

    tmp_meta <- as.data.table(rv_in$ad$var) # can probably just acces ad$obs directly since we don' tneed it to be a data_table?
    row.names(tmp_meta) <- rv_in$ad$var_names
    # DO SUBSET

    X_ID <- in_conf[UI == X_fact]$ID
    group_ys <- rv_in$ad$obs[[ group_by$y ]] [which( dat_js_set %in% dat_js )]
    names(group_ys) <- dat_js


  } else { #x_is <- "obs"
    # subset var (omics)
    if (!is.na( p$feat_subset ) ) {  #maybe don't need this check?
      dat_js <- p$feat_subsel #index by number
      dat_js_set <- rv_in$ad$var[[ p$feat_subset ]]
      omic_js <- rv_in$ad$var_names[ dat_js_set %in% dat_js ]
      #X_data <- X_data[,dat_js_set %in% p$feat_subsel ]
    } else { #subset to omics_list
      dat_js <- omics_list$value

      dat_js_set <- rv_in$ad$var_names
      omic_js <-  dat_js  # dat_js_set[rv_in$ad$var_names %in% dat_js]

      #convert to indices...
      # this already indexes into X_data rows by name
      #X_data <- X_data[,omics_list$value] # value?
    }
    tmp_meta <- as.data.table(rv_in$ad$obs) # can probably just acces ad$obs directly since we don' tneed it to be a data_table?
    row.names(tmp_meta) <- rv_in$ad$obs_names

    group_ys <- rv_in$ad$var[[ group_by$y ]] [which( dat_js_set %in% dat_js )]
    names(group_ys) <- dat_js

    X_ID = "sample_ID"
  }
  # create the table
  hm_data = data.table()
  #TODO:  we already subsetted, so don't need to loop?
  # loop over sample_IDs / subset category
  for(omic_j in omic_js){
    tmp <- tmp_meta[, c( in_conf[UI == X_ID]$ID, in_conf[UI == in_subset]$ID, in_conf[UI == group_by$x]$ID), with = FALSE]
    colnames(tmp) <- c("X_ID", "sub","group_x")
    #tmp$X_ID = tmp_meta[[in_conf[UI == X_fact]$ID]]

    tmp$Y_nm <- omic_j

    tmp$group_y <- group_ys[omic_j]
    #tmp$val = X_data[,which(dat_js_set==dat_j)] #h5data$read(args = list(all_omics[gene_i], quote(expr=)))
    tmp$val = X_data[,omic_j ] #h5data$read(args = list(all_omics[gene_i], quote(expr=)))

    hm_data = rbindlist(list(hm_data, tmp))
  }


  # use X_data for in_data below??
  # Already ahve tmp_meta, in)fact, etc...
  # BOX/VIOLIN DATA ======================
  if (dat_source != "X") {
    if (dat_source == "obs") {
      in_fact <- p$plot_x #in_fact <- p$observ_grpA
      dat_key <- p$plot_y
      in_data <- isolate(rv_in$ad$obs[[dat_key]])
      names(in_data) <- rv_in$ad$obs_names
      group_by <- p$observ_group_by
      in_subset <- p$observ_subset
      in_subsel <- p$observ_subsel

      tmp_meta <- as.data.table(isolate(rv_in$ad$obs))
      row.names(tmp_meta) <- rv_in$ad$obs_names

      x_is <- "obs"


    } else if (dat_source == "var") {
      in_fact <- p$plot_x #in_fact <- p$observ_grpA
      dat_key <- p$plot_y
      in_data <- isolate(rv_in$ad$var[[dat_key]])
      names(in_data) <- isolate(rv_in$ad$var_names)
      group_by <- p$feat_group_by # NOTE: could be individual "omics" from selector
      in_subset <- p$feat_subset
      in_subsel <- p$feat_subsel

      tmp_meta <- as.data.table(isolate(rv_in$ad$var))
      row.names(tmp_meta) <- rv_in$ad$var_names

      x_is <- "var"

    }

    if (is.null(in_subset)) {
      in_subset <- in_fact
      in_subsel <- levels(factor(tmp_meta[[in_fact]]))
    }
  # now subset in_data
    # #SUBSET HERE


  } else {
    # aggregate X  into data vector

    ##TODO: is there a cleaner way to do these means than double coersion?
    ##
    in_data <- as.data.frame(Matrix::rowMeans(X_data[row.names(tmp_meta),],na.rm=TRUE))
    # force to positive???
    #add noise?
    #bx_data[val < 0]$val = 0
    set.seed(42)
    viz_noise <- rnorm(length(in_data)) * diff(range(in_data)) / 1000 # part per thousand  "visualization jitter"
    in_data <- in_data + viz_noise

    in_fact <- X_fact
  }


  # do we need different ID and UI??
  bx_data <- tmp_meta[, c(in_conf[ID == in_fact]$ID,  in_conf[UI == in_subset]$ID),  with = FALSE]
  colnames(bx_data) = c("X", "sub")
  bx_data$val <- in_data



  # SUBSET
  if(length(in_subsel) != 0 & length(in_subsel) != nlevels(bx_data$sub)){
    bx_data <-  bx_data[bx_data$sub %in% in_subsel,]
  }

  if(length(in_subsel) != 0 & length(in_subsel) != nlevels(hm_data$sub)){
    hm_data <-  hm_data[hm_data$sub %in% in_subsel,]
  }





    # #  (dat_source == "X")
    # # use the "obs" fact a
    # # HOW ARE WE KEEPING TRACK OF WHETHER X/Y ARE OBS OR FEAT?
    # # # FOR NOW THEY ARE OBS... NO TRANSPOSE NESCESSARY
    # in_fact <- p$plot_x #in_fact <- p$observ_grpA
    # dat_key <- p$plot_y
    #
    # if (dat_key=="X") {
    #   in_data <- isolate(rv_in$ad$X)
    # } else if (dat_key == "raw") {
    #   in_data <- isolate(rv_in$ad$raw$X)
    # } else if (dat_key == "layers") { #must be a layer
    #   if (dat_key %in% rv_in$ad$layers_keys()) {
    #     in_data <- isolate(rv_in$ad$layers[[dat_key]])
    #   } else {
    #     print("data not found")
    #     return(ret_vals)
    #   }
    #
    # }
    #
    # in_subset <- p$observ_subset
    # in_subsel <- p$observ_subsel
    # # subset by
    # group_by <- list( x = p$observ_group_by,
    #                   y = p$feat_group_by) # NOTE: could be individual "omics" from selector
    #
    # # TODO: make a oaram in side_selector
    # x_is <- "obs" #haven't implimented visualizing by var yet...
    #
    # if (x_is == "var"){
    #   # transpose data
    #   # subset obs (samples)
    #   if (!is.na( in_subset ) ) {
    #     in_data <-  in_data[in_meta[[ in_subset ]] %in% in_subsel ,]
    #
    #   }
    #   in_data <- Matrix::t(in_data)
    #   # switch these around
    #   in_subset <- p$feat_subset
    #   in_subsel <- p$feat_subsel
    #   group_by <- list( y = p$observ_group_by,
    #                     x = p$feat_group_by)
    #
    # } else {
    #   # subset var (omics)
    #   if (!is.na( p$feat_subset ) ) {  #maybe don't need this check?
    #     in_data <- in_data[,rv_in$var[[ p$feat_subset ]] %in% p$feat_subset ]
    #
    #   } else { #subset to omics_list
    #     in_data <- in_data[,omics_list$value]
    #   }
    #
    # }
    #
    #
    # # Aggregate ?
    #
    # # TODO: is X BEING PLOTTED VS var or obs?  need to transpose and aggregate
    # #         FOR NOW ALWASY "OBS"
    # #      SUBSET before AGGREGATE on the aggregating dimension...
    #
    # # #SUBSET HERE ?  currently does this in gg_data assignment
    #



  # TODO: fix subsetting logic... we can only subset the same grouping for metadata...
  #  OR we need to subset the raw-matrix and THEN group...

  # SUBSET LATER?

  if(length(in_subsel) != 0 & length(in_subsel) != nlevels(bx_data$sub)){
    bx_data <-  bx_data[bx_data$sub %in% in_subsel,]
  }









in_conf <- rv_in$config
in_meta <- rv_in$meta


in_fact <- p$plot_x
dat_key= p$plot_y
dat_source <-p$data_source
#in_data <- isolate(rv_in$ad[[dat_source]][[dat_key]])
if (dat_source == "obs") {
  in_data <- (rv_in$ad$obs[[dat_key]])
  names(in_data) <- rv_in$ad$obs_names
} else if (dat_source == "var") {
  in_data <- (rv_in$ad$var[[dat_key]])
  names(in_data) <- (rv_in$ad$var_names)
} else { #  (dat_source == "obs")
  in_data <- rv_in$ad$X
  print('boxplots only for obs and var ?@!?')
  #TODO:  spawn warning box
  return(NULL)
}


in_meta <- as.data.table(rv_in$ad$obs)
row.names(in_meta) <- rv_in$ad$obs_names

in_grp <- rv_in$default$obs_x[1] #p$observ_grpA,
in_subsel <- rv_in$default$grp2 # p$observ_subselA,
in_omic <- omics[1:20] #p$omics_list$value,
plot_type <- "Bubbleplot" #input$RB_plot_type
show_data_points <- TRUE


data_point_sz = .65
plot_size = "Small"
font_size = "Small"
color_schemes = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple")
color_scheme = color_schemes[2]

in_quant <- dat_key #(maybe) just observ_y
pg_violin_box(in_conf, in_meta, in_fact, in_quant,
              in_grp, in_subsel,
              in_data, in_omic,
              plot_type, show_data_points,
              data_point_sz, font_size)






# ------------------------



in_conf <- rv_in$config
#in_meta <- rv_in$meta
in_meta <- as.data.table(rv_in$ad$obs)

in_fact <- p$observ_grpA
#  probably add a "matrix" column to config...
# maybe use the data_source to indicate if we want raw or layers... short circuit for now
in_data <- (rv_in$ad$X)  # do we need to isolate it??

in_quant <- "X" #dat_key #(maybe) just observ_y
# these are the "groups" to show on the x axis
in_group <- in_fact

# these are the groups to show on thye y axis
all_omics <- rv_in$ad$var_names
names(all_omics) <- all_omics

all_obs <- rv_in$ad$obs_names
names(all_obs) <- all_obs

in_grp <- rv_in$default$obs_x[1] #p$observ_grpA,
in_subsel <- rv_in$default$grp2 # p$observ_subselA,
in_omics <- omics[1:20] #p$omics_list$value,
plot_type <- "Bubbleplot" #input$RB_plot_type
in_do_scale =TRUE
in_clust_row =TRUE
in_clust_col =TRUE

pg_bubble_heatmap(in_conf, in_meta, in_omics, in_group, plot_type,
                  in_grp, in_subsel, in_data, all_omics,
                  in_do_scale, in_clust_row, in_clust_col,
                  color_scheme, plot_size)

# pg_bubble_heatmap <- function(in_conf, in_meta, in_omics, in_fact, plot_type,
#                               in_grp, in_subsel, in_data, all_omics, in_do_scale, in_clust_row, in_clust_col,
#                               color_scheme, plot_size,
#                               save = FALSE)

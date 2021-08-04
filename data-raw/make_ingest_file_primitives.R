
# require(reticulate)
# require(hdf5r)
require("data.table")
require("RColorBrewer")


#' Generate data files required for omicser app
#'
#' derived/modified from ShinyCell https://github.com/SGDDNB/ShinyCell
#'
#' Generate data files for easy ingest to internals of Omicser shiny app. Five files will be generated,
#' namely (i) the  config \code{prefix_conf.rds}, (ii) the -omics feature
#' mapping object config \code{prefix_gene.rds}, (iii) the single-cell gene
#' expression \code{prefix_gexpr.h5} for transcriptomic datasets, (iv) the metadata
#' \code{prefix_meta.rds} and (v) the defaults for the Shiny app
#' \code{prefix_def.rds}.
#'
#'
#' A prefix is specified in each file to allow for multiple datasets to be storied in a single
#' repositor/instance of the Omicser app. N
#'
#' all input params are primitives of anndata internals.
#' @param X - count matrix
#' @param obs - metadata on observations (cells)
#' @param var_ - omics metadata (e.g. genes).  not called "var" to prevent problems
#' @param obsm - e.g. dimred
#' @param varm
#' @param uns
#'  es
#' @param X,obs,var_,obsm=NA,varm=NA,uns=NA,
#'
#' @param gex.assay assay in single-cell data object to use for plotting
#'   gene expression, which must match one of the following:
#'   \itemize{
#'     \item{Seurat objects}: "RNA" or "integrated" assay,
#'       default is "RNA"
#'     \item{SCE objects}: "logcounts" or "normcounts" or "counts",
#'       default is "logcounts"
#'     \item{h5ad files}: "X" or any assay in "layers",
#'       default is "X"
#'     \item{loom files}: "matrix" or any assay in "layers",
#'       default is "matrix"
#'   }
#' @param gex.slot slot in single-cell assay to plot. This is only used
#'   for Seurat objects (v3+). Default is to use the "data" slot
#'
#' @param gene_mapping specifies whether to convert human / mouse Ensembl gene
#'   IDs (e.g. ENSG000xxx / ENSMUSG000xxx) into "user-friendly" gene symbols.
#'   Set this to \code{TRUE} if you are using Ensembl gene IDs. Default is
#'   \code{FALSE} which is not to perform any conversion. Alternatively, users
#'   can supply a named vector where \code{names(gene_mapping)} correspond
#'   to the actual gene identifiers in the gene expression matrix and
#'   \code{gene_mapping} correspond to new identifiers to map to
#'
#' @param db_prefix specify file prefix
#' @param db_dir specify directory to place files (e.g. ingest / data-raw)
#'
#' @param default_omic specify secondary default omics to show
#' @param default_omic character vector specifying default omics_ to
#'   show in bubbleplot / heatmap
#' @param default_dimred character vector specifying the two default dimension
#'   reductions. Default is to use UMAP if not TSNE embeddings
#' @param chunk_size number of genes written to h5file at any one time. Lower
#'   this number to reduce memory consumption. Should not be less than 10
#'
#' @return data files required for ingest
#'  ie.  "XXconf.rds"
#'        "XXXmeta.rds
#'       "XXXomics.rds"
#'       "XXXdef.rds"
#'
#' @author Andy Henrie (from og code https://github.com/SGDDNB/ShinyCell)
#'
#' @import data.table hdf5r reticulate hdf5r
#'
#' @examples
#' make_ingest_files(seu, ui_conf,
#'   gex.assay = "RNA", gex.slot = "data",
#'   db_prefix = "sc1", db_dir = "shinyApp/",
#'   default.gene1 = "GATA3", default.gene2 = "DNMT3L",
#'   default.multigene = c(
#'     "ANPEP", "NANOG", "ZIC2", "NLGN4X", "DNMT3L",
#'     "DPPA5", "SLC7A2", "GATA3", "KRT19"
#'   ),
#'   default.dimred = c("UMAP_1", "UMAP_2")
#' )
#' @export
make_ingest_file_primitives <- function(X,obs,var_,obsm,varm,uns, layers,
                                        observables, comparables,dimreds,
                                        default_omic = NA, default_dimred = NA, meta_to_include = NA,
                                        gex.assay = NA, gex.slot = c("data", "scale.data", "counts"),
                                        gene_mapping = FALSE, db_prefix = "test1", db_dir = "data-raw",
                                        chunk_size = 500,  legend_cols = 4,
                                        max_levels_ui = 50){
  max_levels <- 100
  ##################################
  ## PREPROCESS
  ###################################
  ### Preprocessing and checks

  X_matdim <- rev(dim(X))
  X_rownm <- rownames(var_)
  X_colnm <- rownames(obs)
  def_omics = X_rownm[1:10]

  obs_meta = data.table(sampleID = X_colnm)  # redundant but makes naming explicit..
  obs_meta = cbind(obs_meta, obs)
  #colnames(obs_meta) = c("sampleID", colnames(obs))

  for (i_meta in colnames(obs_meta)[-1]) {
    obs_meta[[i_meta]] = unlist(obs_meta[[i_meta]]) # unlist and refactor
    if (is.character(obs_meta[[i_meta]])) {
      levels = sort(unique(obs_meta[[i_meta]]))
      if (length(levels) < max_levels) {
        obs_meta[[i_meta]] = factor(obs_meta[[i_meta]], levels = levels)
      }
    }
  }


  ############

  obsm_exist <-  !is.null(obsm)
  varm_exist <-  !is.null(varm)
  uns_exist <-  !is.null(uns)


  if (isTRUE(obsm_exist)) {
    obsm_keys <- names(obsm)
    keys <- obsm_keys
    i_v <- 'obsm'
    m_conf <- data.table()
    for (key_i in keys) {
      tmp_conf <- data.table(
        ID = i_v, UI = i_v, fID = NA, fUI = NA,
        fCL = NA, fRow = NA, default = 0,  observ=FALSE,comp=FALSE,dimred=FALSE
      )

      if (is.data.frame(obsm[[key_i]])) {
        i_levels <- levels(obsm[[key_i]])
      } else {
        i_levels <- 1:dim(obsm[[key_i]])[2] #maybe should check uns?
      }
      n_levels <- length(i_levels)
      if (n_levels <= max_levels) {
        if (n_levels >= 2) {
          tmp_conf$fID <- paste0(i_levels, collapse = "|")
          tmp_conf$fUI <- tmp_conf$fID
          tmp_conf$fCL <- paste0(colorRampPalette(brewer.pal(12, "Paired"))(n_levels),
                                 collapse = "|"
          )
          tmp_conf$fRow <- ceiling(n_levels / legend_cols)
          tmp_conf$grp <- TRUE
        } else if (n_levels == 1) {
          tmp_conf$fID <- NA
          tmp_conf$fUI <- tmp_conf$fID
          tmp_conf$fCL <- "black"
          tmp_conf$fRow <- 1
        }
      }
      #TODO: test for comp measure...
      if (!is.na(dimreds$obsm)){
        if (key_i %in% dimreds$obsm){
          tmp_conf$dimred <- TRUE
        }
      }
      if (!is.na(observables$obsm[1])){
        if (key_i %in% comparables$obsm){
          tmp_conf$observ <- TRUE
        }
      }
      if (!is.na(comparables$obsm[1])){
        if (key_i %in% comparables$obsm){
          tmp_conf$comp <- TRUE
        }
      }
      m_conf <- rbindlist(list(m_conf, tmp_conf))

    }
    obsm_conf <- m_conf
  } else {
    print('no matrix obs')
    obsm_keys <- NA
    obsm_conf <- data.table()

  }


  if (isTRUE(varm_exist)) {
    varm_keys <- names(varm)

    keys <- varm_keys
    i_v <- 'varm'
    m_conf <- data.table()
    for (key_i in keys) {
      tmp_conf <- data.table(
        ID = i_v, UI = i_v, fID = NA, fUI = NA,
        fCL = NA, fRow = NA, default = 0,  observ=FALSE,comp=FALSE,dimred=FALSE
      )

      if (is.data.frame(varm[[key_i]])) {
        i_levels <- levels(varm[[key_i]])
      } else {
        i_levels <- 1:dim(varm[[key_i]])[2]
      }
      n_levels <- length(i_levels)
      if (n_levels <= max_levels) {
        if (n_levels >= 2) {
          tmp_conf$fID <- paste0(i_levels, collapse = "|")
          tmp_conf$fUI <- tmp_conf$fID
          tmp_conf$fCL <- paste0(colorRampPalette(brewer.pal(12, "Paired"))(n_levels),
                                 collapse = "|"
          )
          tmp_conf$fRow <- ceiling(n_levels / legend_cols)
          tmp_conf$grp <- TRUE
        } else if (n_levels == 1) {
          tmp_conf$fID <- NA
          tmp_conf$fUI <- tmp_conf$fID
          tmp_conf$fCL <- "black"
          tmp_conf$fRow <- 1
        }
      }
        #TODO: test for comp measure...
      if (!is.na(dimreds$varm[1])){
        if (key_i %in% dimreds$varm){
          tmp_conf$dimred <- TRUE
        }
      }

      if (!is.na(comparables$varm[1])){
        if (key_i %in% comparables$varm){
          tmp_conf$comp <- TRUE
        }
      }
    m_conf <- rbindlist(list(m_conf, tmp_conf))
    }
    varm_conf <- m_conf
  } else {
    print('no matrix vars')
    varm_keys <- NA
    varm_conf <- data.table()

  }

  if (uns_exist) {
    uns_keys <- names(uns)
    keys <- uns_keys
    i_v <- 'uns'
    # Start making config data.table
    m_conf <- data.table()
    tmp_conf <- data.table(
      ID = i_v, UI = i_v, fID = NA, fUI = NA,
      fCL = NA, fRow = NA, default = 0,  observ=FALSE,comp=FALSE,dimred=FALSE
    )
    # Additional preprocessing for categorical metadata
    n_levels <- length(keys)
    if (n_levels <= max_levels) {
      if (n_levels >= 2) {
        tmp_conf$fID <- paste0(keys, collapse = "|")
        tmp_conf$fUI <- tmp_conf$fID
        tmp_conf$fCL <- paste0(colorRampPalette(brewer.pal(12, "Paired"))(n_levels),collapse = "|" )
        tmp_conf$fRow <- ceiling(n_levels / legend_cols)
      } else if (n_levels == 1) {
        tmp_conf$fID <- keys
        tmp_conf$fUI <- tmp_conf$fID
        tmp_conf$fCL <- "black"
        tmp_conf$fRow <- 1
      }

    }
    uns_conf <- rbindlist(list(m_conf, tmp_conf))

    # # add the varm keys to the uns
    # if (varm_exist) {
    #   uns_keys <- names(varm)
    #   for (i_uns in uns_keys) {
    #     uns[[i_uns]] <- colnames(varm[[i_uns]])
    #   }
    # }
    # if (obsm_exist) {
    #   uns_keys <- names(obsm)
    #   for (i_uns in uns_keys) {
    #     uns[[i_uns]] <- colnames(obsm[[i_uns]])
    #   }
    # }
  } else {
    uns_keys <- NA
    print('no unstructured meta')
    uns_conf <- data.table()
    }




  # now add the rest of the observables...
  if (!is.na(observables$var[1])){

    keys <- observables$var

    i_v <- 'var'
    m_conf <- data.table()
    for (key_i in keys) {
      tmp_conf <- data.table(
        ID = i_v, UI = i_v, fID = NA, fUI = NA,
        fCL = NA, fRow = NA, default = 0,  observ=FALSE,comp=FALSE,dimred=FALSE
      )
      tmp_conf$fID <- key_i
      tmp_conf$fUI <- tmp_conf$fID

      tmp_conf$observ <- TRUE
      m_conf <- rbindlist(list(m_conf, tmp_conf))
    }
    var_conf <- m_conf
  } else {var_conf <- data.table() }

  # now add the rest of the observables...
  if (!is.na(observables$obs[1])){

    keys <- observables$obs

    i_v <- 'obs'
    m_conf <- data.table()
    for (key_i in keys) {
      tmp_conf <- data.table(
        ID = i_v, UI = i_v, fID = NA, fUI = NA,
        fCL = NA, fRow = NA, default = 0,  observ=FALSE,comp=FALSE,dimred=FALSE
      )
      tmp_conf$fID <- key_i
      tmp_conf$fUI <- tmp_conf$fID

      tmp_conf$observ <- TRUE
      m_conf <- rbindlist(list(m_conf, tmp_conf))
    }
    obs_conf <- m_conf
  } else { obs_conf <- data.table() }

  # now add the rest of the observables...
  if (!is.na(observables$layers[1])){

    keys <- observables$layers

    i_v <- 'layer'
    m_conf <- data.table()
    for (key_i in keys) {
      tmp_conf <- data.table(
        ID = i_v, UI = i_v, fID = NA, fUI = NA,
        fCL = NA, fRow = NA, default = 0,  observ=FALSE,comp=FALSE,dimred=FALSE
      )
      tmp_conf$fID <- key_i
      tmp_conf$fUI <- tmp_conf$fID

      tmp_conf$observ <- TRUE
      m_conf <- rbindlist(list(m_conf, tmp_conf))
    }
    layer_conf <- m_conf
  } else { layer_conf <- data.table() }

  # now add the rest of the observables...
  if (!is.na(observables$raw[1])){

    keys <- observables$raw

    i_v <- 'raw'
    m_conf <- data.table()
    for (key_i in keys) {
      tmp_conf <- data.table(
        ID = i_v, UI = i_v, fID = NA, fUI = NA,
        fCL = NA, fRow = NA, default = 0,  observ=FALSE,comp=FALSE,dimred=FALSE
      )
      tmp_conf$fID <- key_i
      tmp_conf$fUI <- tmp_conf$fID

      tmp_conf$observ <- TRUE
      m_conf <- rbindlist(list(m_conf, tmp_conf))
    }
    raw_conf <- m_conf
  } else { raw_conf <- data.table() }

  m_conf <- rbindlist(list(obsm_conf, varm_conf,uns_conf, obs_conf, var_conf, layer_conf, raw_conf), fill= TRUE)
  dimred_exist <- any(m_conf$dimred)


  # Checks and get list of metadata to include
  if (is.na(meta_to_include[1])) {
    meta_to_include <- colnames(obs_meta)
  }
  if (length(meta_to_include) < 2) {
    stop("At least 2 metadata is required!")
  }

  # Start making config data.table
  ui_conf <- data.table()
  for (i_meta in meta_to_include) {
    tmp_conf <- data.table(
      ID = i_meta, UI = i_meta, fID = NA, fUI = NA,
      fCL = NA, fRow = NA, default = 0, grp = FALSE
    )

    # Additional preprocessing for categorical metadata
    n_levels <- nlevels(obs_meta[[i_meta]])
    print(levels(obs_meta[[i_meta]]))
    if (n_levels <= max_levels) {
      if (n_levels >= 2) {
        tmp_conf$fID <- paste0(levels(obs_meta[[i_meta]]), collapse = "|")
        tmp_conf$fUI <- tmp_conf$fID
        tmp_conf$fCL <- paste0(colorRampPalette(brewer.pal(12, "Paired"))(n_levels),
                               collapse = "|"
        )
        tmp_conf$fRow <- ceiling(n_levels / legend_cols)
        tmp_conf$grp <- TRUE
      } else if (n_levels == 1) {
        tmp_conf$fID <- levels(obs_meta[[i_meta]])
        tmp_conf$fUI <- tmp_conf$fID
        tmp_conf$fCL <- "black"
        tmp_conf$fRow <- 1
      }
      #TODO: test for comp measure...
      ui_conf <- rbindlist(list(ui_conf, tmp_conf))
    }
  }

  # this isn't actually doing anything... leave for now
  # TODO:  remove $default
  # Set defaults
  def1 <- grep("ident|library", ui_conf$ID, ignore.case = TRUE)[1]
  def2 <- grep("clust", ui_conf$ID, ignore.case = TRUE)
  def2 <- setdiff(def2, def1)[1]
  if (is.na(def1)) {
    def1 <- setdiff(c(1, 2), def2)[1]
  }
  if (is.na(def2)) {
    def2 <- setdiff(c(1, 2), def1)[1]
  }

  ui_conf[def1]$default <- 1
  ui_conf[def2]$default <- 2

  # STOP if there is no single multi-level covariate
  if (nrow(ui_conf[grp == TRUE]) == 0) {
    stop(paste0(
      "did not detect any multi-group cell metadata \n",
      "       e.g. library / cluster. Has any analysis been performed?"
    ))
  }

  ##################################
  # Perform dimred_exist if specified (also map defGenes)
  ###################################
  if (gene_mapping[1] == TRUE) {
    if (sum(grepl("^ENSG000", X_rownm)) >= sum(grepl("^ENMUSG000", X_rownm))) {
      tmp1 = fread(system.file("extdata", "geneMapHS.txt.gz",
                               package = "omicser"
      ))
    } else {
      tmp1 = fread(system.file("extdata", "geneMapMM.txt.gz",
                               package = "omicser"
      ))
    }
    gene_mapping = tmp1$geneName
    names(gene_mapping)= tmp1$geneID
  }
  # TODO:  fix the if else logic here... .or maybe gen_mapping[1] could be set to FALSE above?
  # Check if dimred_exist is partial or not
  if (gene_mapping[1] == FALSE) {
    gene_mapping <- X_rownm
    names(gene_mapping) <- X_rownm # Basically no mapping
  } else {
    if (!all(X_rownm %in% names(gene_mapping))) {
      # warning("Mapping for some gene identifiers are not provided!")
      tmp1 <- X_rownm[X_rownm %in% names(gene_mapping)]
      tmp1 <- gene_mapping[tmp1]
      tmp2 <- X_rownm[!X_rownm %in% names(gene_mapping)]
      names(tmp2) <- tmp2
      gene_mapping <- c(tmp1, tmp2)
    }
    gene_mapping <- gene_mapping[X_rownm]
  }
  def_omics <- gene_mapping[def_omics]


  # Check default.gene1 / default.gene2 / default.multigene
  if (is.na(default_omic[1])) {
    default_omic <- def_omics
  } else if (all(default_omic %in% gene_mapping)) {
    default_omic <- default_omic
  } else if (all(default_omic %in% names(gene_mapping))) {
    default_omic <- gene_mapping[default_omic]
  } else {
    warning(paste0(
      "default_omic doesn't exist in gene expression, ",
      "using defaults..."
    ))
    default_omic <- def_omics
  }


  ##################################
  # Make XXXmeta.rds and XXXconf.rds (updated with dimred info)
  ###################################
  ui_conf$dimred <- FALSE
  obs_meta <- obs_meta[, c("sampleID", as.character(ui_conf$ID)), with = FALSE]
  # Factor metadata again
  for (i in as.character(ui_conf[!is.na(fID)]$ID)) {
    obs_meta[[i]] <- factor(obs_meta[[i]],
                            levels = strsplit(ui_conf[ID == i]$fID, "\\|")[[1]]
    )
    levels(obs_meta[[i]]) <- strsplit(ui_conf[ID == i]$fUI, "\\|")[[1]]
    ui_conf[ID == i]$fID <- ui_conf[ID == i]$fUI
  }

  if (isTRUE(dimred_exist)) {
    # TODO:  code extracting this from obsm "pca", "tsne", "umap" etc
    # Extract dimred and append to both XXXmeta.rds and XXXconf.rds...
    for (idr in colnames(obsm)) { #py_to_r(inp_H5$obsm_keys())) {
      dr_mat = obsm[[idr]] #py_to_r(inp_H5$obsm)[[idr]]
      #dr_mat = py_to_r(inp_H5$obsm[idr])
      tmp_name = gsub("pca", "pc", gsub("X_", "", idr))
      tmp_name = paste0(tmp_name, "_", 1:ncol(dr_mat))
      colnames(dr_mat) <- tmp_name
      if (ncol(dr_mat) > 5) {
        dr_mat = dr_mat[, 1:5]
      } # Take first 5 components only
      dr_mat = as.data.table(dr_mat)
      obs_meta = cbind(obs_meta, dr_mat)

      # Update ui_conf accordingly
      tmp = data.table(
        ID = colnames(dr_mat), UI = colnames(dr_mat),
        fID = NA, fUI = NA, fCL = NA, fRow = NA,
        default = 0, grp = FALSE, dimred = TRUE
      )
      tmp$UI = gsub("_", "", tmp$UI) # now underscores allowed
      ui_conf = rbindlist(list(ui_conf, tmp))
    }
  } else {
    print(" no dimension reductions available.  compute them?")
  }


  ##################################
  # Make XXXgexpr.h5
  ###################################
  if (!dir.exists(db_dir)) {
    dir.create(db_dir)
  }


  #scGEX <- t(X)


  if (!isTRUE(all.equal(obs_meta$sampleID, X_colnm))) {
    obs_meta$sampleID <- factor(obs_meta$sampleID, levels = X_colnm)
    obs_meta <- obs_meta[order(sampleID)]
    obs_meta$sampleID <- as.character(obs_meta$sampleID)
  }

  ###################################
  # Make XXXgenes.rds
  ####################################
  om1_omics <- seq(X_matdim[1])
  names(gene_mapping) <- NULL
  names(om1_omics) <- gene_mapping

  om1_omics <- om1_omics[order(names(om1_omics))]
  om1_omics <- om1_omics[order(nchar(names(om1_omics)))]


  # ###################################
  # # Make XXXdef.rds (list of defaults)
  # #####################################
  #
  # if ( isTRUE(dimred_exist) ) {
  #   # wwe have dimred available so make the defaults
  #   if (all(default_dimred %in% ui_conf[dimred == TRUE]$ID)) {
  #     default_dimred[1] <- ui_conf[ID == default_dimred[1]]$UI
  #     default_dimred[2] <- ui_conf[ID == default_dimred[2]]$UI
  #   } else if (all(default_dimred %in% ui_conf[dimred == TRUE]$UI)) {
  #     default_dimred <- default_dimred # Nothing happens
  #   } else {
  #     warn <- TRUE
  #     if (is.na(default_dimred[1])) {
  #       default_dimred <- "umap"
  #       warn <- FALSE
  #     }
  #     # Try to guess... and give a warning
  #     guess <- gsub("[0-9]", "", default_dimred[1])
  #     if (length(grep(guess, ui_conf[dimred == TRUE]$UI, ignore.case = TRUE)) >= 2) {
  #       default_dimred <- ui_conf[dimred == TRUE]$UI[
  #         grep(guess, ui_conf[dimred == TRUE]$UI, ignore.case = TRUE)[1:2]
  #       ]
  #     } else {
  #       nDR <- length(ui_conf[dimred == TRUE]$UI)
  #       default_dimred <- ui_conf[dimred == TRUE]$UI[(nDR - 1):nDR]
  #     }
  #     if (warn) {
  #       warning(paste0(
  #         "default_dimred not found, switching to ",
  #         default_dimred[1], " and ", default_dimred[1]
  #       ))
  #     } # Warn if user-supplied default_dimred is not found
  #   }
  #
  # } else { #dimred_exist
  #   # no dimension reduction (dimred_exist==FALSE) so assert default_dimred is NA
  #   default_dimred <- NA
  # }

  # Note that we stored the display name here
  # # make sure that our 1 and 2 are likely to be sensible things to visualize....
  #
  om1_def <- list()
  om1_def$meta1 <- ui_conf[default == 1]$UI # Use display name
  om1_def$meta2 <- ui_conf[default == 2]$UI # Use display name
  om1_def$omics1 <- default_omic[1] # Actual == Display name
  om1_def$omics2 <- default_omic[2] # Actual == Display name
  om1_def$omics <- default_omic # Actual == Display name
  om1_def$dimred <- default_dimred # Use display name
  tmp <- nrow(ui_conf[default != 0 & grp == TRUE])
  if (tmp == 2) {
    om1_def$grp1 <- om1_def$meta1
    om1_def$grp2 <- om1_def$meta2
  } else if (tmp == 1) {
    om1_def$grp1 <- ui_conf[default != 0 & grp == TRUE]$UI
    if (nrow(ui_conf[default == 0 & grp == TRUE]) == 0) {
      om1_def$grp2 <- om1_def$grp1
    } else {
      om1_def$grp2 <- ui_conf[default == 0 & grp == TRUE]$UI[1]
    }
  } else {
    om1_def$grp1 <- ui_conf[default == 0 & grp == TRUE]$UI[1]
    if (nrow(ui_conf[default == 0 & grp == TRUE]) < 2) {
      om1_def$grp2 <- om1_def$grp1
    } else {
      om1_def$grp2 <- ui_conf[default == 0 & grp == TRUE]$UI[2]
    }
  }

  # TODO: use the m_conf default/grp to set things

  # remove fUI and default from conf...
  ui_conf <- ui_conf[, -c("fUI", "default"), with = FALSE]
  m_conf <- m_conf[, -c("fUI", "default"), with = FALSE]

  ### Saving objects
  ###   TODO: sort by om1_omics FIRST
  ###
  ###
  conf <- list(meta=ui_conf,mat=m_conf)
  saveRDS(conf, file = paste0(db_dir, "/", db_prefix, "conf.rds"))
  saveRDS(obs_meta, file = paste0(db_dir, "/", db_prefix, "meta.rds"))
  saveRDS(om1_omics, file = paste0(db_dir, "/", db_prefix, "omics.rds"))
  saveRDS(om1_def, file = paste0(db_dir, "/", db_prefix, "def.rds"))


  ###   sort X by omics
  saveRDS(X, file = paste0(db_dir, "/", db_prefix, "X.rds"))

  saveRDS(obs, file = paste0(db_dir, "/", db_prefix, "obs.rds"))

  saveRDS(var_, file = paste0(db_dir, "/", db_prefix, "var.rds"))
  saveRDS(layers, file = paste0(db_dir, "/", db_prefix, "layers.rds"))
  saveRDS(obsm, file = paste0(db_dir, "/", db_prefix, "obsm.rds"))
  saveRDS(varm, file = paste0(db_dir, "/", db_prefix, "varm.rds"))
  saveRDS(uns, file = paste0(db_dir, "/", db_prefix, "uns.rds"))

}

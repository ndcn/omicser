#' playground UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_playground_ui <- function(id){
  ns <- NS(id)
  tagList(
    # TODO: change these to dynamically render or not based on "selector"
    # # could be "side selector" or some radio/dropdown choices...

    tabsetPanel(
      type = 'pills',  #'hidden' and a radio might work best
      id = 'tab',
      # summary stats tab
      tabPanel(
        title = "Expression", value = 'raw',
        mod_pg_vis_raw_ui(id=ns("pg_vis_raw_ui_1"))
      ),
      # volcano tab
      tabPanel(
        title = "Diff. Expr.",value = 'comp',
        mod_pg_vis_comp_ui(id=ns("pg_vis_comp_ui_1"))

      ),
      # table tab
      tabPanel(
        title = "Table", value='table',
        mod_pg_table_ui(id=ns("pg_pg_table_ui_1"))

      ),
      # QC tab
      tabPanel(
        title = "QC",value = 'qc',
        mod_pg_vis_qc_ui(id=ns("pg_vis_qc_ui_1"))

      )
      ) #tabsetpanel


  )
}

#' playground Server Functions
#'
#' @noRd
mod_playground_server <- function(id ,rv_data, rv_selections) {
  moduleServer( id, function(input, output, session){
    ns <- session$ns

# MODULES =================================
    mod_pg_table_server("pg_pg_table_ui_1",rv_data, rv_selections, active_layer_data)

    mod_pg_vis_raw_server("pg_vis_raw_ui_1",rv_data, rv_selections, heat_data)#, agg_heat) #,box_data,varbox_data)
    mod_pg_vis_comp_server("pg_vis_comp_ui_1",rv_data, rv_selections, active_layer_data)
    mod_pg_vis_qc_server("pg_vis_qc_ui_1",rv_data, rv_selections)

# REACTIVEVALUES =================================
    heat_data <- reactiveValues(
      x_names = NULL,
      y_names = NULL,
      x_group = NULL,
      y_group = NULL,
      x_source =NULL,
      subsel = NULL,
      subset = NULL,
      type = NULL,
      data = NULL,
      mat = NULL,
      meta = NULL,
      ready = FALSE
    )

    # TODO:  get units/label for dat_loc
    # grab the right matrix for the heatmap and for the volcano plot distribution)
    active_layer_data <- reactiveValues(
      layer = NULL,
      data = NULL
    )





# OBSERVES =================================
# active_layer_data observe =================================
    observe({
      req(rv_data$anndata,
          rv_selections$data_layer)
      layer <- rv_selections$data_layer
      if (layer=="X") {
        X_data <- rv_data$anndata$X
      } else if (layer == "raw") {
        X_data <- rv_data$anndata$raw$X
      } else { #} if (dat_loc == "layers") { #must be a layer
        is_layer <- any(layer %in% rv_data$anndata$layers$keys())
        #is_layer <- any(rv_data$anndata$layers$keys()==dat_loc)
        if (is_layer) {
          #X_data <- isolate(rv_data$anndata$layers[[dat_loc]]) #isolate
          X_data <- rv_data$anndata$layers$get(layer)
        } else {
          message("data not found")
          X_data <- NULL
        }
      }

      if(is.null(dimnames(X_data)[[1]]) & !is.null(X_data) ){
        dimnames(X_data) <- list(rv_data$anndata$obs_names,rv_data$anndata$var_names)
      }

      active_layer_data$layer <- layer
      active_layer_data$data <- X_data # is a matrix
      message("set -------> active_layer")

    })




# heat_data observe =================================
    observe({
      req(rv_data$anndata,
          rv_selections$data_layer,
          active_layer_data$data)

      # NOTES:
      #    this packs a reactive list of data for generating the heatmap.
      #    we need to return a
      #    - "filtered"/"subsetted data matrix
      #       - we need to "subset" the omics (10-2000ish omics)
      #       - subset the "samples" by meta-label
      #   - grouping variables
      #     - for samples (which will also be aggregated)
      #     - for omics
      #
      #  have already set a reactive active_layer_data$data -> X_data
      #
    message("in observer: heat_data packer")
      in_conf <- rv_data$config
      dat_source <- rv_selections$data_layer
      # not using the viz_now... for now
      X_data <- active_layer_data$data

      # this is all of the "active" omics (subsetting in side-selector)
      omics <- rv_selections$selected_omics$all_omics
      # # send signal back to side selector...
      # rv_selections$selected_omics$freeze <- FALSE
      omics_idx <- which(rv_data$anndata$var_names %in% omics)

      #   - subset omics... bug above :: heat_data observe =================================
      #omics_idx <- which(rv_data$anndata$var[[rv_selections$feat_subset]] %in% rv_selections$feat_subsel)
      #omics <- rv_data$anndata$var_names[omics_idx]

      # omic_grp_nm <- rv_selections$feat_group_by
      # # st var (omics)
      # if (!is.na( omic_grp_nm ) | !is.null(omic_grp_nm)) {
      #   omic_grp <- rv_data$anndata$var[[ omic_grp_nm ]][omics_idx]
      # } else { #subset to selected_omics since there was NO meta-category to subset against
      #   omic_grp <- (omics_idx>0) #hack a single group...
      # }

      # samples...
      # 1. subset
      # 2. grouping

      # FILTER the data matrix
      #

      #   - subset samples :: heat_data observe =================================
      samples_idx <- which(rv_data$anndata$obs[[rv_selections$observ_subset]] %in% rv_selections$observ_subsel)
      samples <- rv_data$anndata$obs_names[samples_idx]

      samp_grp_nm <- rv_selections$observ_group_by

      if (!is.na( samp_grp_nm ) | !is.null(samp_grp_nm)) {
        samp_grp <- rv_data$anndata$obs[[ samp_grp_nm ]][samples_idx]
      } else { #subset to all samles since there was NO meta-category to subset against (SHOULD NEVER HAPPEN)
        samp_grp <- (samples_idx>0) #hack a single group...
      }


      X_filtered <- X_data[samples,omics]

      X_tab <- as.data.table(X_filtered)

      # convert to table for bubbleplot?
      #always keep the sample_ID column
      X_ID = "sample_ID"
      # ---------- : prep data_matrix for heat_map/bubble
      X_fact <- rv_selections$observ_group_by

      tmp_meta <- as.data.table(rv_data$anndata$obs) # can probably just access anndata$obs directly since we don' tneed it to be a data_table?
      #row.names(tmp_meta) <- rv_data$anndata$obs_names
      tmp_meta <- tmp_meta[samples_idx,]
      #row.names(tmp_meta) <- samples
      #tmp_meta$sample_ID <- samples


      # # create the table
      # hm_data = data.table()
      # #TODO:  we already subsetted, so don't need to loop?  # can't enable "transpose" version until
      # # loop over sample_IDs / subset category
      # # this will be for bubbleplots?
      # keep_cols <- c( in_conf[UI == X_ID]$ID,
      #                 in_conf[UI == samp_grp_nm]$ID,
      #                 in_conf[UI == rv_selections$observ_subset]$ID,
      #                 in_conf[UI == samp_grp_nm]$ID)
      #
      #
      # grp_ys <- omic_grp
      # names(grp_ys) <- omics
      # for(omic_j in omics){
      #   tmp <- tmp_meta[, keep_cols, with = FALSE]
      #   colnames(tmp) <- c("X_ID", "X_nm", "sub","group_x")
      #   #tmp$X_ID = tmp_meta[[in_conf[UI == X_fact]$ID]]
      #   tmp$Y_nm <- omic_j
      #
      #   tmp$group_y <- grp_ys[omic_j]
      #   tmp$val = X_filtered[,omic_j ]
      #
      #   hm_data = rbindlist(list(hm_data, tmp))
      # }

# browser()
#       # TODO:  use data.table melt (or something rather than loop)
#       # melt(DT, id.vars = c("family_id", "age_mother"),
#       #      measure.vars = c("dob_child1", "dob_child2", "dob_child3"))
#       #
#       xtab <- as.data.table(t(X_filtered))
#       xtab$group_y <- omic_grp
#       xtab$omic <- omics
#
#       xtab2 <- as.data.table(X_filtered)
#       xtab2$group_x <- samp_grp
#

      message("---->  finishing: heat_data packer")


      heat_data$x_names <- samp_grp
      heat_data$y_names <- omics
      heat_data$x_group <- samp_grp_nm

       heat_data$type <- active_layer_data$layer #   dat_loc
       heat_data$data <- NULL #hm_data
       heat_data$mat <- X_filtered
       heat_data$obs_meta <- tmp_meta
       heat_data$ready <- TRUE
       heat_data$selected_omics <- rv_selections$selected_omics

       message("---->  DONE: heat_data packer")

    })




    #TODO:  split head_data into
    #
    #   agg_mat
    #   filt_mat
    #
    #       heat_data$aggregated
    #       heat_data$full



  })
}
    ##############################
    ##############################
    ##############################
    ##############################



## To be copied in the UI
# mod_playground_ui("playground_ui_1")

## To be copied in the server
# mod_playground_server("playground_ui_1")

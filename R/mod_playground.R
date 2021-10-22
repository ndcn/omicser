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

      )
      #,
      # # QC tab DEPRICATED FOR NOW
      # tabPanel(
      #   title = "QC",value = 'qc',
      #   mod_pg_vis_qc_ui(id=ns("pg_vis_qc_ui_1"))
      #
      # )
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
    #mod_pg_vis_qc_server("pg_vis_qc_ui_1",rv_data, rv_selections)

# REACTIVEVALUES =================================
    heat_data <- reactiveValues(
      x_names = NULL,
      y_names = NULL,
      x_group = NULL,
      y_group = NULL, #DEPRICATE?


      data = NULL,
      mat = NULL,
      obs_meta = NULL,
      ready = FALSE,

      x_annot = NULL,
      y_annot = NULL,
      selected_omics = NULL

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
          message("data layer not found")
          X_data <- NULL
        }
      }

      if(is.null(dimnames(X_data)[[1]]) & !is.null(X_data) ){
        dimnames(X_data) <- list(rv_data$anndata$obs_names,rv_data$anndata$var_names)
      }

      active_layer_data$layer <- layer
      active_layer_data$data <- X_data # is a matrix
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

      X_data <- active_layer_data$data

      # this is all of the "active" omics (subsetting in side-selector)
      omics <- rv_selections$selected_omics$all_omics
      omics_idx <- which(rv_data$anndata$var_names %in% omics)

      #   - subset samples :: heat_data observe =================================
      samples_idx <- which(rv_data$anndata$obs[[rv_selections$observ_subset]] %in% rv_selections$observ_subsel)
      samples <- rv_data$anndata$obs_names[samples_idx]

      samp_grp_nm <- rv_selections$observ_group_by

      if (!is.na( samp_grp_nm ) | !is.null(samp_grp_nm)) {
        samp_grp <- rv_data$anndata$obs[[ samp_grp_nm ]][samples_idx]
      } else { #subset to all samles since there was NO meta-category to subset against (SHOULD NEVER HAPPEN)
        samp_grp <- (samples_idx>0) #hack a single group...
        message(">>>>>>>>>>>>no sample group !?!")
      }


      X_filtered <- X_data[samples,omics]

      obs_meta <- as.data.table(rv_data$anndata$obs) # can probably just access anndata$obs directly since we don' tneed it to be a data_table?
      obs_meta <- obs_meta[samples_idx,]
      samp_annot <- obs_meta[ ,rv_data$shaddow_defs$exp_annot, with=FALSE]

      var_meta <- as.data.table(rv_data$anndata$var)
      feat_annot <- var_meta[omics_idx,rv_data$shaddow_defs$feat_annot,with=FALSE]


      message("---->  finishing: heat_data packer")

      heat_data$x_annot <- samp_annot
      heat_data$y_annot <- feat_annot
      heat_data$x_names <- samp_grp
      heat_data$y_names <- omics
      heat_data$x_group <- samp_grp_nm

       heat_data$type <- active_layer_data$layer #   dat_loc
       heat_data$data <- NULL #hm_data
       heat_data$mat <- X_filtered
       heat_data$obs_meta <- obs_meta  #might not need this any more
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

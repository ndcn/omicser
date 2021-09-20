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
        title = "Quantities", value = 'raw',
        mod_pg_vis_raw_ui(id=ns("pg_vis_raw_ui_1"))
      ),
      # volcano tab
      tabPanel(
        title = "Differential",value = 'comp',
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
    mod_pg_table_server("pg_pg_table_ui_1",rv_data, rv_selections) #active_layer_data ?
    mod_pg_vis_raw_server("pg_vis_raw_ui_1",rv_data, rv_selections,heat_data,box_data,varbox_data)
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

    box_data <- reactiveValues(
      x_name = NULL,
      y_name = NULL,
      colors = NULL,
      dat_source =NULL,
      data = NULL
    )


    varbox_data <- reactiveValues(
      x_name = NULL,
      y_name = NULL,
      colors = NULL,
      dat_source =NULL,
      data = NULL
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
      req(rv_data$ad,
          rv_selections$data_layer)

      layer <- rv_selections$data_layer
      if (layer=="X") {
        X_data <- isolate(rv_data$ad$X) #isolate
      } else if (layer == "raw") {
        X_data <- isolate(rv_data$ad$raw$X)#isolate
      } else { #} if (dat_loc == "layers") { #must be a layer
        is_layer <- any(layer %in% rv_data$ad$layers$keys())
        #is_layer <- any(rv_data$ad$layers$keys()==dat_loc)
        if (is_layer) {
          #X_data <- isolate(rv_data$ad$layers[[dat_loc]]) #isolate
          X_data <- isolate(rv_data$ad$layers$get(dat_loc) )
        } else {
          print("data not found")
          return(ret_vals)
        }
      }

      if(is.null(dimnames(X_data)[[1]])){
        dimnames(X_data) <- list(rv_data$ad$obs_names,rv_data$ad$var_names)
      }

      active_layer_data$layer <- layer
      active_layer_data$data <- X_data

    })

# heat_data observe =================================
    observe({
      req(rv_data$ad,
          rv_selections$data_source,
          active_layer_data$data)

      in_conf <- rv_data$config
      in_meta <- as.data.table(rv_data$ad$obs)
      row.names(in_meta) <- rv_data$ad$obs_names

      dat_source <- rv_selections$data_source

      group_action <- rv_selections$group_action


      #   if (heat_data$ready) {
      rv_selections$omics_list$viz_now <- FALSE
      # ---------- : prep data_matrix for heat_map/bubble
      # TODO: make an "x_is" oaram in side_selector
      x_is <- "obs" #haven't implimented visualizing by var yet...
      X_fact <- rv_selections$plot_x #in_fact <- rv_selections$observ_grpA


      X_data <- active_layer_data$data
      ## AGGREGATE along x-axis
      # if (x_is == "obs"){

        dat_j_grp <-rv_selections$feat_subset
        # subset var (omics)
        if (!is.na( dat_j_grp ) ) {  #maybe don't need this check?
          dat_js <- rv_selections$feat_subsel #index by number
          dat_js_set <- rv_data$ad$var[[ dat_j_grp ]]
          omic_js <- rv_data$ad$var_names[ dat_js_set %in% dat_js ]
          #X_data <- X_data[,dat_js_set %in% rv_selections$feat_subsel ]
        } else { #subset to omics_list
          dat_js <- rv_selections$omics_list$value
          dat_js_set <- rv_data$ad$var_names
          omic_js <-  dat_js  # dat_js_set[rv_data$ad$var_names %in% dat_js]
        }

        tmp_meta <- as.data.table(rv_data$ad$obs) # can probably just acces ad$obs directly since we don' tneed it to be a data_table?
        row.names(tmp_meta) <- rv_data$ad$obs_names

        if (is.na(rv_selections$observ_subset)) {
          in_subset <-  rv_selections$plot_x# don't try and subset at the end...
          in_subsel <- character(0)
        } else {
          in_subset <-  rv_selections$observ_subset
          in_subsel <- rv_selections$observ_subsel
        }

        if (group_action == "none" | rv_selections$observ_group_by=="NA") {
          # group along the plotting variable
          grp_x <- rv_selections$plot_x
        } else {
          grp_x <- rv_selections$observ_group_by
        }

        grp_y <- rv_selections$feat_group_by # NOTE: could be "omics" from selector
        if ((grp_y == "NA") | (grp_y == "")) {
          grp_ys <- dat_js
        } else {
          grp_ys <- rv_data$ad$var[[ grp_y ]][which( dat_js_set %in% dat_js )]
        }
        #group_ys <- rv_data$ad$var[[ grp_by$y ]] [which( dat_js_set %in% dat_js )]
        names(grp_ys) <- dat_js

        X_ID = "sample_ID"


      # create the table
      hm_data = data.table()
      #TODO:  we already subsetted, so don't need to loop?  # can't enable "transpose" version until
      # loop over sample_IDs / subset category

      browser()
      keep_cols <- c( in_conf[UI == X_ID]$ID, in_conf[UI == X_fact]$ID, in_conf[UI == in_subset]$ID, in_conf[UI == grp_x]$ID)
      for(omic_j in omic_js){
        tmp <- tmp_meta[, keep_cols, with = FALSE]
        colnames(tmp) <- c("X_ID", "X_nm", "sub","group_x")
        #tmp$X_ID = tmp_meta[[in_conf[UI == X_fact]$ID]]
        tmp$Y_nm <- omic_j

        tmp$group_y <- grp_ys[omic_j]
        tmp$val = X_data[,omic_j ]

        hm_data = rbindlist(list(hm_data, tmp))
      }

      if(length(in_subsel) != 0 & length(in_subsel) != nlevels(hm_data$sub)){
        hm_data <-  hm_data[hm_data$sub %in% in_subsel,]
      }


       heat_data$x_names <- unique(hm_data$X_nm)
       heat_data$y_names <- dat_js
       heat_data$x_group <- grp_x
       heat_data$y_group <- grp_y
       heat_data$x_source <- x_is
       heat_data$subsel <- in_subsel
       heat_data$subset <- in_subset
       heat_data$type <- active_layer_data$layer #   dat_loc
       heat_data$data <- hm_data
       heat_data$mat <- X_data
       heat_data$meta <- tmp_meta
       heat_data$ready <- TRUE

    }
    )


# box_data observe (req head_data) =================================
    observe({
        req(rv_data$ad,
            rv_selections$data_source)

      in_conf <- rv_data$config
      in_meta <- as.data.table(rv_data$ad$obs)
      row.names(in_meta) <- rv_data$ad$obs_names

      dat_source <- rv_selections$data_source
      if (dat_source == "obs") {
        in_fact <- rv_selections$plot_x #in_fact <- rv_selections$observ_grpA
        dat_key <- rv_selections$plot_y
        in_data <- isolate(rv_data$ad$obs[[dat_key]])
        names(in_data) <- rv_data$ad$obs_names
        grp_by <- rv_selections$observ_group_by

        if (is.na(rv_selections$observ_subset)) {
          in_subset <-  character(0) # NOTE: disabled set to NA
          in_subsel <- character(0)  # NOTE: disabled set to NA
        } else {
          in_subset <-  rv_selections$observ_subset
          in_subsel <- rv_selections$observ_subsel
        }

        tmp_meta <- as.data.table(isolate(rv_data$ad$obs))
        row.names(tmp_meta) <- rv_data$ad$obs_names

        x_is <- "obs"
        if (is.null(in_subset)) {
          in_subset <- in_fact
          in_subsel <- levels(factor(tmp_meta[[in_fact]]))
        }
        # now subset in_data

      } else { # is "X"

        if (is.null(heat_data$x_source)) {

          print("not data from hm_data yet")
          return()
        } else {
          dat_key <- heat_data$type
          # aggregate X  into data vector
          X_data <- heat_data$mat
          tmp_meta <- heat_data$meta
          ##TODO: is there a cleaner way to do these means than double coersion?
          ##
          in_data <- as.data.frame(Matrix::rowMeans(X_data[row.names(tmp_meta),],na.rm=TRUE))
          # force to positive???
          #add noise?
          #bx_data[val < 0]$val = 0
          set.seed(42)
          viz_noise <- rnorm(length(in_data)) * diff(range(in_data)) / 1000 # part per thousand  "visualization jitter"
          in_data <- in_data + viz_noise

          in_fact <- rv_selections$plot_x

          # copy these from heat_data
          in_subset <- heat_data$subset
          in_subsel <- heat_data$subsel
          grp_by <- heat_data$x_group
          x_is <- heat_data$x_source

        }
      }

      # do we need different ID and UI??
      bx_data <- tmp_meta[, c(in_conf[ID == in_fact]$ID, in_conf[ID == grp_by]$ID,     in_conf[UI == in_subset]$ID),  with = FALSE]
      colnames(bx_data) = c("X","grp","sub")
      bx_data$val <- in_data

      # # set the color of our boxes from config
      # ctmp <- strsplit(in_conf[UI == in_fact]$fCL, "\\|")[[1]]
      # names(ctmp) <- strsplit(in_conf[UI == in_fact]$fID, "\\|")[[1]]
      # bx_data$color <- ctmp[bx_data$X]
      # # map color onto the plot from our config
      ctmp = strsplit(in_conf[UI == in_fact]$fCL, "\\|")[[1]]  #pre-computed colors..
      names(ctmp) = levels(bx_data$X)
      #refactor X
      gg_lvl = levels(bx_data$X)[levels(bx_data$X) %in% unique(bx_data$X)]
      bx_data$X = factor(bx_data$X, levels = gg_lvl)
      grp_colors <- ctmp[gg_lvl]

      # SUBSET
      if(length(in_subsel) != 0 & length(in_subsel) != nlevels(bx_data$sub)){
        bx_data <-  bx_data[bx_data$sub %in% in_subsel,]
      }

      box_data$x_name <- in_fact
      box_data$y_name <- dat_key
      box_data$colors <- grp_colors
      box_data$dat_source  <- dat_source
      box_data$data <- bx_data

      }
    ) # observe box_data


# var_box_data observe ===================
      observe({
        req(rv_data$ad,
            rv_selections$data_source,
            rv_selections$plot_feats)



        group_action <- rv_selections$group_action

        in_conf <- rv_data$config
        in_meta <- as.data.table(rv_data$ad$obs)
        row.names(in_meta) <- rv_data$ad$obs_names

        in_fact <- rv_selections$plot_var_x #in_fact <- rv_selections$observ_grpA
        dat_key <- rv_selections$plot_var_y
        in_data <- isolate(rv_data$ad$var[[dat_key]])
        names(in_data) <- isolate(rv_data$ad$var_names)

        if (group_action == "none" | rv_selections$feat_group_by=="NA") {
          # group along the plotting variable
          grp_by <- in_fact
        } else {
          grp_by <- rv_selections$feat_group_by
        }

        if (is.na(rv_selections$feat_subset)) {
          in_subset <-  in_fact     # NOTE: disabled set to NA
          in_subsel <- character(0)  # NOTE: disabled set to NA
        } else {
          in_subset <-  rv_selections$feat_subset
          in_subsel <- rv_selections$feat_subsel
        }

        tmp_meta <- as.data.table(isolate(rv_data$ad$var))
        row.names(tmp_meta) <- rv_data$ad$var_names

        x_is <- "var"


        # do we need different ID and UI??
        bx_data <- tmp_meta[, c(in_conf[ID == in_fact]$ID, in_conf[ID == grp_by]$ID, in_conf[UI == in_subset]$ID),  with = FALSE]

        # we should never get here, but hack an early return if we don't have x&y available
        if (!dim(bx_data)[1]){
          return(varbox_data)
        }

        cols <- c("X","grp","sub")
        colnames(bx_data) = cols
        bx_data <- bx_data[,(cols) := lapply(.SD, as.factor), .SDcols = cols] # force the columns to be factors (i.e. if logical)

        bx_data$val <- in_data

        # # set the color of our boxes from config
        # ctmp <- strsplit(in_conf[UI == in_fact]$fCL, "\\|")[[1]]
        # names(ctmp) <- strsplit(in_conf[UI == in_fact]$fID, "\\|")[[1]]
        # bx_data$color <- ctmp[bx_data$X]
        # # map color onto the plot from our config
        ctmp = strsplit(in_conf[UI == in_fact]$fCL, "\\|")[[1]]  #pre-computed colors..
        names(ctmp) = levels(bx_data$X)
        #refactor X
        gg_lvl = levels(bx_data$X)[levels(bx_data$X) %in% unique(bx_data$X)]
        bx_data$X = factor(bx_data$X, levels = gg_lvl)
        grp_colors <- ctmp[gg_lvl]

        # SUBSET
        if(length(in_subsel) != 0 & length(in_subsel) != nlevels(bx_data$sub)){
          bx_data <-  bx_data[bx_data$sub %in% in_subsel,]
        }

        varbox_data$x_name <- in_fact
        varbox_data$y_name <- dat_key
        varbox_data$colors <- grp_colors
        varbox_data$dat_source  <- rv_selections$data_source
        varbox_data$data <- bx_data

      }
      )


  })
}

## To be copied in the UI
# mod_playground_ui("playground_ui_1")

## To be copied in the server
# mod_playground_server("playground_ui_1")

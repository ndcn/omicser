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
      tabPanel(
        title = "Table", value='table',
        mod_pg_table_ui(id=ns("pg_table_ui_2"))
      ),
      # ingest tab
      tabPanel(
        title = "Quantities", value = 'raw',
        mod_pg_vis_raw_ui(id=ns("pg_vis_raw_ui_1"))
      ),
      # table tab
      tabPanel(
        title = "Differential",value = 'comp',
        mod_pg_vis_comp_ui(id=ns("pg_vis_comp_ui_1"))

      ),
      # table tab
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
mod_playground_server <- function(id ,rv_in, p) {
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    mod_pg_table_server("pg_table_ui_2",rv_in, p)
    mod_pg_vis_raw_server("pg_vis_raw_ui_1",rv_in, p,heat_data,box_data)
    mod_pg_vis_comp_server("pg_vis_comp_ui_1",rv_in, p)
    mod_pg_vis_qc_server("pg_vis_qc_ui_1",rv_in, p)




    # p$data_source
    # p$plot_x
    # p$plot_y
    # p$observ_subset
    # p$observ_subsel
    # # p$feat_subset
    # # p$feat_subsel
    # # group (plotting)
    # p$feat_group_by
    # p$observ_group_by
    # # aggregate collapsing data matrix
    # p$feat_agg
    # p$observ_agg
    # p$group_action
    # p[["omics_list"]]

    observe({
      req(p$omics_list)
      print("---> viz?? ")
      if (heat_data$ready) {
        p$omics_list$viz_now <- FALSE  # reset it here...
      }
    })

    heat_data <- reactiveValues(
      x_names = NULL,
      y_names = NULL,
      x_source =NULL,
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


    # send this to the "raw" tab
    observe({
      req(rv_in$ad,
          p$data_source)

      # need: data, x_name,y_name,
      # data -> table with $variable
      #                    $value
      #                    $group (if group == variable then non-grouped...
      #                    else grouped along x_name)
      #
      #
print("in reactive: hm_data (playground)")
        # ret_vals <- list(
        #   x_names = NULL,
        #   y_names = NULL,
        #   x_source =NULL,
        #   type = NULL,
        #   data = NULL,
        #   mat = NULL,
        #   meta = NULL,
        #   ready = FALSE
        # )

      in_conf <- rv_in$config
      in_meta <- as.data.table(rv_in$ad$obs)
      row.names(in_meta) <- rv_in$ad$obs_names

      dat_source <- p$data_source

      group_action <- p$group_action

      # HEATMAP DATA ======================
      # ---------- : prep data_matrix for heat_map/bubble
      # TODO: make an "x_is" oaram in side_selector
      x_is <- "obs" #haven't implimented visualizing by var yet...
      X_fact <- p$plot_x #in_fact <- p$observ_grpA
      dat_loc <- p$plot_y
      dat_loc <- "X"
      if (dat_loc=="X") {
        X_data <- isolate(rv_in$ad$X) #isolate
      } else if (dat_loc == "raw") {
        X_data <- isolate(rv_in$ad$raw$X)#isolate
      } else if (dat_loc == "layers") { #must be a layer
        if (dat_loc %in% rv_in$ad$layers_keys()) {
          X_data <- isolate(rv_in$ad$layers[[dat_loc]])#isolate
        } else {
          print("data not found")
          return(ret_vals)
        }
      }

      # makesure i have rownames right for X_data?
      in_subset <- p$observ_subset
      in_subsel <- p$observ_subsel
      # subset by

      grp_by <- list( x = p$observ_group_by,
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

        grp_by <- list( y = p$observ_group_by,
                          x = p$feat_group_by)

        tmp_meta <- as.data.table(rv_in$ad$var) # can probably just acces ad$obs directly since we don' tneed it to be a data_table?
        row.names(tmp_meta) <- rv_in$ad$var_names
        # DO SUBSET

        X_ID <- in_conf[UI == X_fact]$ID

        if (grp_by$y != "NA") {
          group_ys <- rv_in$ad$obs[[ grp_by$y ]][which( dat_js_set %in% dat_js )]
        } else {
          group_ys <- dat_js
        }
        #group_ys <- rv_in$ad$obs[[ grp_by$y ]] [which( dat_js_set %in% dat_js )]
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

        if (grp_by$y != "NA") {
          group_ys <- rv_in$ad$var[[ grp_by$y ]][which( dat_js_set %in% dat_js )]
        } else {
          group_ys <- dat_js
        }
        #group_ys <- rv_in$ad$var[[ grp_by$y ]] [which( dat_js_set %in% dat_js )]
        names(group_ys) <- dat_js

        X_ID = "sample_ID"
      }
      # create the table
      hm_data = data.table()
      #TODO:  we already subsetted, so don't need to loop?
      # loop over sample_IDs / subset category
      keep_cols <- c( in_conf[UI == X_ID]$ID, in_conf[UI == X_fact]$ID, in_conf[UI == in_subset]$ID, in_conf[UI == grp_by$x]$ID)
      for(omic_j in omic_js){
        tmp <- tmp_meta[, keep_cols, with = FALSE]
        colnames(tmp) <- c("X_ID", "X_nm", "sub","group_x")
        #tmp$X_ID = tmp_meta[[in_conf[UI == X_fact]$ID]]
        tmp$Y_nm <- omic_j

        tmp$group_y <- group_ys[omic_j]
        #tmp$val = X_data[,which(dat_js_set==dat_j)]
        tmp$val = X_data[,omic_j ]

        hm_data = rbindlist(list(hm_data, tmp))
      }

      if(length(in_subsel) != 0 & length(in_subsel) != nlevels(hm_data$sub)){
        hm_data <-  hm_data[hm_data$sub %in% in_subsel,]
      }


       heat_data$x_names <- unique(hm_data$X_nm)
       heat_data$y_names <- dat_js
       heat_data$x_source <- x_is
       heat_data$type <- dat_loc
       heat_data$data <- hm_data
       heat_data$mat <- X_data
       heat_data$meta <- tmp_meta
       heat_data$ready <- TRUE
#
#      return(ret_vals)
    }
    )

      # send hm_data to heatmap_reactive
      # send X_data to violin (in case we need to collapse X)



    # box_data <- reactive({
    #   req(rv_in$ad,
    #         p$data_source)
      observe({
        req(rv_in$ad,
            p$data_source)

        # hm_data <- list(
        #   x_names = X_fact,
        #   y_names = grp_by$y,
        #   x_source = x_is,
        #   type = dat_loc,
        #   data = hm_data,
        #   mat = X_data,
        #   meta = tmp_meta,
        #   ready = TRUE
        # )
      print("in reactive: box_data (playground)")


      in_conf <- rv_in$config
      in_meta <- as.data.table(rv_in$ad$obs)
      row.names(in_meta) <- rv_in$ad$obs_names

      dat_source <- p$data_source
      # use X_data for in_data below??
      # Already ahve tmp_meta, in)fact, etc...
      # BOX/VIOLIN DATA ======================
      if (dat_source != "X") {
        if (dat_source == "obs") {
          in_fact <- p$plot_x #in_fact <- p$observ_grpA
          dat_key <- p$plot_y
          in_data <- isolate(rv_in$ad$obs[[dat_key]])
          names(in_data) <- rv_in$ad$obs_names
          grp_by <- p$observ_group_by
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
          grp_by <- p$feat_group_by # NOTE: could be individual "omics" from selector
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

        if ( !is.null(heat_data$x_source)){
          # hm_data <- list(
          #   x_names = X_fact,
          #   y_names = grp_by$y,
          #   x_source = x_is,
          #   type = dat_loc,
          #   data = hm_data,
          #   mat = X_data,
          #   meta = tmp_meta,
          #   ready = TRUE
          # )
          print("not data from hm_data yet")
        } else {

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

          in_fact <- p$plot_x
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
  )




  })
}

## To be copied in the UI
# mod_playground_ui("playground_ui_1")

## To be copied in the server
# mod_playground_server("playground_ui_1")

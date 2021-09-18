

#' side_selector UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_side_selector_ui <- function(id){
  ns <- NS(id)

  quant_plot_types <- c("Dot Plot"="dotplot","Heat Map"="hmap","Box/violin Plot"="box")
  comp_plot_types <- c("Volcano Plot"="volcano", "Dot Plot"="dotplot","Heat Map"="hmap","Box/violin Plot"="box")

  selector_tags <- tagList(
    # helpText(
    #    HTML("Load a Dataset first, and then selected the variables and quantiteis below.")
    #    ),
    #  hr(style = "border-top: 1px solid #000000;"),
    # #TODO:  change these to render prper HTML

    htmlOutput( ns( "ui_curr_database" )),
    uiOutput( ns( "ui_DIV_warn" )),
    htmlOutput( ns( "ui_db_type" )),

    fluidRow(
      hr(style = "border-top: 1px solid #000000;"),
      h4("Choose plotting variables-X & measure-Y"),
      column(
        width=2,
        "choose:"
      ),
      column(
        width = 5,
        offset = 0,
        shinyjs::disabled(radioButtons( ns("RB_obs_X"), "X and Y from:",
                                        choices = c("sample annot. (obs)"="obs",
                                                    # "omic annot. (var)"="var",
                                                    "data matrix (X)" = "X"),
                                        selected = "obs" , inline = TRUE)
        )
      ),
      column(
        width = 5,
        offset = 0,
        uiOutput(ns("ui_data_layer"))
      )
    ),
    hr(style = "border-top: 1px solid #000000;"),

    uiOutput(ns("ui_xy_select")),

    hr(style = "border-top: 1px solid #000000;"),
    mod_omic_selector_ui(ns("omic_selector_ui_1")),

    hr(style = "border-top: 1px solid #000000;"),
    uiOutput(ns("ui_obs_subset")),

    hr(style = "border-top: 1px dashed grey;"),

    fluidRow(
      column(
        width = 8,
        h4("Omic annotation summary:")
      ),
      column(
        width=3,
        offset=1,
        checkboxInput( ns("CB_xy_var_select"), "show?",
                       value = FALSE)
      )
    ),
    uiOutput(ns("ui_xy_var_select")),
    hr(style = "border-top: 1px dashed grey;"),

    fluidRow(
      uiOutput(ns("ui_text_warn"), width = "100%"),
        ), #fluidRow

    fluidRow(
      column(
        width = 10,
        offset = 1,
        uiOutput(ns("ui_groupB"))
      )
    )

  ) #taglist

  return(selector_tags)

}

#' side_selector Server Functions
#'
#' @noRd
mod_side_selector_server <- function(id, rv_data){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    ### Reactive expressions ============================================
    rv_selections <- reactiveValues(

      data_source = NULL,
      data_layer = NULL,
      #omics_names = NULL,
      omics_list = NULL,

      plot_x = NULL,
      plot_y = NULL,

      group_action = NULL,
      observ_group_by = NULL,
      observ_subset = NULL,
      observ_subsel = NULL,

      # aggregate obs
      feat_group_by = NULL,

      plot_var_x = NULL,
      plot_var_y = NULL,
      feat_subset = NULL,   # NOT ENABLEDj, using omics selector
      feat_subsel = NULL,  # NOT ENABLED
      # feature X & Y?
      plot_feats = NULL,

      # TODO:  make this conditional on what kind of data is loaded.. and just a single plot type
      raw_plot_type = NULL,
      comp_plot_type = NULL
    )
    ### OMICS  =========================================================

    all_omics <- reactive( names(rv_data$omics) )
    def_omics <- reactive( rv_data$default$omics)
    new_db_trig <- reactive( rv_data$trigger )
    omics_list <- mod_omic_selector_server("omic_selector_ui_1", all_omics ,def_omics,new_db_trig)

    ### Outputs =========================================================
    output$ui_curr_database <- renderUI({
      if (is.null(rv_data$database_name)) {
        out_text <- "No data loaded"
      } else {
          out_text <- paste("Current databse: ", rv_data$database_name)
      }
      out_text <- h4(out_text)
      return(out_text)
    })

    # Warning if no data loaded
    output$ui_DIV_warn <- renderUI( {
      if (is.null(rv_data$database_name)) {
        div(
          tags$br(),
          span(class = "warn", "No dataset loaded")
         )
        }
      })

    # show the factors that have been loaded
    output$ui_db_type <- renderUI({
      req(rv_data$omics_type)
      out_text <- paste("<h6>Data type: <i>", rv_data$omics_type, "</i>-omics</h6>")
      out_text <- HTML(out_text)
      return(out_text)
      })


    observeEvent(
      (rv_data$trigger>0),  #why does this happen twice?
      {
        req(rv_data$config,  # set when database is chosen
            rv_data$default)
        rv_data$trigger <- 0
        # update full omics list
        all_omics <- rv_data$omics
        shinyjs::enable("RB_obs_X")
      },
    ignoreNULL = TRUE,
    ignoreInit = TRUE
    ) #observe event



    ## dynamic subset UI group
    output$ui_data_layer <- renderUI({
      req(input$RB_obs_X,
          rv_data$config)
      choices <- rv_data$config[field=="layer"]$UI  # X, raw, or layers
      # default data_source is obs
      ret_tags <-  selectizeInput(ns("SI_data_layer"),
                           "matrix data:",
                           choices =  choices,
                           selected = choices[1])

      return(ret_tags)
    })


    ## dynamic subset UI group
    output$ui_obs_subset <- renderUI({
      req(input$RB_obs_X,
          rv_data$config)

        subs_label <- paste0("select ",isolate(input$SI_obs_subset),"s: ")
        choices_obs_x <- rv_data$config[grp == TRUE & field=="obs"]$UI

        ret_tags <- tagList(
          fluidRow(
          column(
            width = 5,
            selectizeInput(ns("SI_obs_subset"),
                             "subset samples:",
                           choices =  choices_obs_x,
                           selected = rv_data$default$obs_x),
            br(),
              actionButton(ns("CB_obs_sub_none"), "NONE", class = "btn btn-primary" ),
              actionButton(ns("CB_obs_sub_all"), "ALL", class = "btn btn-primary")
          ),
          column(width = 7,
                 offset = 0,
                 uiOutput(ns("ui_obs_subset_sel")),
          ),

        ) #fluidrow
        ) #taglist

        return(ret_tags)
    })


    ## dynamic subset selector checkboxes
    output$ui_obs_subset_sel <- renderUI({
      req(input$RB_obs_X,
          input$SI_obs_subset,
          rv_data$config)


      # check to see if we have the input$SI_obs_subset in our config
      # if not return andwait for the SI_obs_subset to be updated
      if ( !(input$SI_obs_subset %in% rv_data$config$UI) ) {
        print("input$SI_obs_subset no yet updated")
        return(NULL)
      }

      cfg <- isolate(rv_data$config)

      subs <- strsplit(cfg[UI == input$SI_obs_subset]$fID, "\\|")[[1]]
      # sorting hack
      subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
      if (all(!is.na(subs2))){
        names(subs2) <- seq_along(subs2)
        subs <- subs[as.numeric(names(sort(subs2)))]
      }
      subs_label <- paste0("select ",isolate(input$SI_obs_subset),"s: ")

      ret_tags <- tagList(
          checkboxGroupInput( ns("CB_obs_subsel"),
                              label = subs_label,
                              inline = TRUE,
                              choices = subs,
                              selected = subs)
      ) #taglist

      return(ret_tags)

    })

    # dynamic x and y selector
    output$ui_xy_select <- renderUI({
      req(input$RB_obs_X,
          rv_data$config) # do i need this?
      # TODO:  add "NA" to group and subset options...

      if (input$RB_obs_X == "obs") {
        #render X and Y from obs
        choices_x <- rv_data$config[grp == TRUE & field=="obs"]$UI
        choices_y <- rv_data$config[measure == TRUE & field=="obs"]$UI

        group_obs <- rv_data$config[grp == TRUE & field=="obs"]$UI # <- choices_x
        group_var <- NA                                         #


        def_x <- rv_data$default$obs_x
        def_y <- rv_data$default$obs_y
        def_grp_o <- rv_data$default$obs_subset
        def_grp_v <- group_var[1]

      } else { # X
        # choices X and group_obs same?
        choices_x <- rv_data$config[grp == TRUE & field=="obs"]$UI  # X, raw, or layers
        choices_y <- rv_data$config[measure == TRUE & field=="var"]$UI# X, raw, or layers
        #TODO: change this spot to toggle x and y?
        choices_y <- c("<omic selector>",choices_y)

        # add "omics" to group_var
        group_obs <- rv_data$config[grp == TRUE & field=="obs"]$UI # subset for aggregating?
        group_var <- rv_data$config[grp == TRUE & field=="var"]$UI # subset for aggregating

        #x_is_obs_or_var = rv_data$config[measure == TRUE & field=="var"]$UI # are we grouping by obs or var?
        def_x <- rv_data$default$obs_x
        def_y <- rv_data$default$obs_y
        def_grp_o <- rv_data$default$obs_subset
        def_grp_v <- rv_data$default$var_subset
      }


      to_return <-  tagList(
        fluidRow(
          column(
            width=2,
            "Choose:",br(),
            "X & Y",
            "to viz"
          ),
          column(
            width=5,
            offset=0,
            selectizeInput(ns("SI_x"),
                           label = "variable (x-axis)",
                           choices = choices_x,
                           selected = def_x)
          ),
          column(
            width = 5,
            offset = 0,
            selectizeInput(ns("SI_y"),
                           label = "measure (y-axis):",
                           choices = choices_y,
                           selected = def_y)

          )
        ),
        hr(style = "border-top: 1px dashed grey;"),
        fluidRow(
          column(
            width=2,
            offset=0,
            radioButtons(ns("RB_none_or_grp"),
                         label = "grouping",
                         choices = c("none","group by"),
                         selected = "none")
          ),
          column(
            width=5,
            offset=0,
            shinyjs::disabled(
              selectizeInput(ns("SI_group_obs"),
                             label = "observ group (samples)",
                             choices = group_obs,
                             selected = def_grp_o)
            )
          ),
          column(
            width = 5,
            offset = 0,
            shinyjs::disabled(
              selectizeInput(ns("SI_group_var"),
                             label = "feature group (omics)",
                             choices = group_var,
                             selected = def_grp_v)
            )
          )
        )#fluidRow
      ) #tagList

      return(to_return)

    })

    # dynamic x and y selector
    output$ui_xy_var_select <- renderUI({
      req(rv_data$config) # do i need this?



      if (input$CB_xy_var_select) {

        choices_x <- rv_data$config[grp == TRUE & field=="var"]$UI
        choices_y <- rv_data$config[measure == TRUE & field=="var"]$UI
        group_var <- c("<omic selector>",choices_x) # rv_data$config[grp == TRUE & field=="var"]$UI #


        to_return <-  tagList(
          fluidRow(
            column(
              width=2,
              "Choose:",br(),
              "group & meas"
            ),
            column(
              width=5,
              offset=0,
              selectizeInput(ns("SI_var_x"),
                             label = "omic group (x)",
                             choices = choices_x,
                             selected = choices_x[1])
            ),
            column(
              width = 5,
              offset = 0,
              selectizeInput(ns("SI_var_y"),
                             label = "measure (y):",
                             choices = choices_y,
                             selected = choices_y[1])

            )
          )
        )

        return(to_return)
      } else {

        return(NULL)

      }

    })


    ### observe s =========================================================
    observeEvent(input$CB_obs_sub_all, {
      req(rv_data$config)
      cfg <- isolate(rv_data$config)

      subs <- strsplit(cfg[UI == input$SI_obs_subset]$fID, "\\|")[[1]]
      subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
      if (all(!is.na(subs2))){
        names(subs2) <- seq_along(subs2)
        subs <- subs[as.numeric(names(sort(subs2)))]
      }
      subs_label <- paste0("select ",isolate(input$SI_obs_subset),"s: ")

      freezeReactiveValue(input, "CB_obs_subsel")
      updateCheckboxGroupInput(session,
                               inputId = "CB_obs_subsel",
                               label = subs_label,
                               choices = subs,
                               selected = subs,
                               inline = TRUE) # WARNING.  make sure subset is null is checked (length(0?))
    })


    observeEvent(input$CB_obs_sub_none, {
      req(rv_data$config)

      cfg <- isolate(rv_data$config)

      subs <- strsplit(cfg[UI == input$SI_obs_subset]$fID, "\\|")[[1]]
      subs2 <- as.numeric(gsub("[^[:digit:]]", "", subs))
      if (all(!is.na(subs2))){
        names(subs2) <- seq_along(subs2)
        subs <- subs[as.numeric(names(sort(subs2)))]
      }

      subs_label <- paste0("select ",isolate(input$SI_obs_subset),"s: ")

      freezeReactiveValue(input, "CB_obs_subsel")
      updateCheckboxGroupInput(session,
                               inputId = "CB_obs_subsel",
                               label = subs_label,
                               choices = subs,
                               selected = character(0), # character(0), #, "" , NULL
                               inline = TRUE) # WARNING.  make sure subset is null is checked (length(0?))
    })


    observe( {
      req(input$RB_none_or_grp)
      # toggle group selector
      if ( input$RB_none_or_grp=="none" ) {
        shinyjs::disable("SI_group_obs")
        #shinyjs::disable("SI_group_var")
      } else {
        shinyjs::enable("SI_group_obs")
        #shinyjs::enable("SI_group_var")
      }
    })

    observe({
      req(input$SI_y,
          input$RB_obs_X)

      if ( input$RB_obs_X=="X" ) {
        shinyjs::disable("SI_y")
      } else {
        shinyjs::enable("SI_y")
      }
    })


    observe({
      req(rv_data$config)

      choices_x <- rv_data$config[grp == TRUE & field=="var"]$UI
      choices_y <- rv_data$config[measure == TRUE & field=="var"]$UI

      if (isTruthy(choices_x) & isTruthy(choices_y) ) {
        shinyjs::enable("CB_xy_var_select")
      } else {
        shinyjs::disable("CB_xy_var_select")
      }
    })


    #TODO:  make the choices *named* paste0(rv_data$config$fID, fUI)
    # dat_source = rv_data$config$mat[fID == rv_data$aux_raw]$ID
    # paste(raw_choices, dat_sourcce)
    #{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }"

    #
      # updateSelectizeInput(session, "SI_y_info", server = TRUE,
      #                      choices = sel2_choices,
      #                      selected = sel2)
    # })
    # TODO:  wrap this in a module..
    #
    #

    observeEvent( omics_list$viz_now, {
      # route the chosen data type out...
      req(input$CB_obs_subsel,
          input$SI_obs_subset)
      if (omics_list$viz_now) {

        #TODO: add data_source to config values
        rv_selections$data_source <- input$RB_obs_X
        rv_selections$data_layer <- input$SI_data_layer

        rv_selections$plot_x <- input$SI_x
        rv_selections$plot_y <- input$SI_y

        rv_selections$observ_subset <- input$SI_obs_subset
        rv_selections$observ_subsel <- input$CB_obs_subsel

        # group (plotting)
        rv_selections$feat_group_by <- input$SI_group_var
        rv_selections$observ_group_by <- input$SI_group_obs

        rv_selections$group_action <- input$RB_none_or_grp

        # # Disabled
        rv_selections$feat_subset = NA # NOT ENABLEDj, using omics selector
        rv_selections$feat_subsel = NA # NOT ENABLED

        rv_selections$plot_var_x <- input$SI_var_x
        rv_selections$plot_var_y <- input$SI_var_y
        rv_selections$plot_feats <- input$CB_xy_var_select

        rv_selections$omics_list <- omics_list  # value & viz_now
        # do the subsetting here?  create a reactive "data blob"?

      } else {
        rv_selections$data_source <- NULL #blocks playground reactive from computing until we "submit" ($viz_now)
      }

    })


  ### RETURN =========================================================
  return(rv_selections)

  })
}

## To be copied in the UI
# mod_side_selector_ui("side_selector_ui_1")

## To be copied in the server
# mod_side_selector_server("side_selector_ui_1")

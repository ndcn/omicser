

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


  selector_tags <- tagList(
    # helpText(
    #    HTML("Load a Dataset first, and then selected the variables and quantiteis below.")
    #    ),
    #  hr(style = "border-top: 1px solid #000000;"),
    # #TODO:  change these to render prper HTML

    uiOutput( ns( "ui_DIV_warn" )),
    htmlOutput( ns( "ui_db_type" )),

    h1("Variables & *-Omics Selector"),
    hr(style = "border-top: 1px solid #000000;"),

    fluidRow(
      column(width = 3, offset = 0,style="border-right: 2px solid black",
             htmlOutput( ns( "ui_curr_database" ))
      ),
      column(width = 4, offset = 0,style="border-right: 2px solid black",
             uiOutput(ns("ui_data_layer"))
      ),
      column(width = 4, offset = 1,#style="border-right: 2px solid black",
             #style='padding-left:0px; padding-right:1px',
             actionButton(ns("AB_plot_now"), "Plot Now!", class = "btn btn-primary btn-large" ) #class="hidableDefault"
      )
    ),

    hr(style = "border-top: 1px solid #000000;"),
    h4("Experimental Factors"),
    fluidRow(
      column(width = 2,
             offset = 0,
             "Sample",br(),
             "meta-",
             "data"
      ),
      column(width = 10,
             offset = 0,
              uiOutput(ns("ui_meta_select")),
      )
    ),

    # omics selection --------------------
    hr(style = "border-top: 1px solid #000000;"),
    h4("*-Omics"),
    fluidRow(
      column(width = 2,
             offset = 0,
             "*-Omics",br(),
             "meta-",
             "data"),
      column(width = 10,
             offset = 0,

             uiOutput(ns("ui_omic_select")),
             hr(style = "border-top: 1px dashed grey;"),
             #TODO: CHANGE THIS "TARGET OMICS" OR SOMETHING
             mod_omic_selector_ui(ns("omic_selector_ui_1")),

      )
    ),

    hr(style = "border-top: 1px dashed grey;"),

    fluidRow(
      uiOutput(ns("ui_text_warn"), width = "100%")
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

      data_layer = NULL,
      #omics_names = NULL,
      selected_omics = NULL,

      observ_group_by = NULL,
      observ_group_by2 = NULL,
      observ_subset = NULL,
      observ_subsel = NULL,

      # aggregate obs
      feat_group_by = NULL,
      feat_subset = NULL,   # NOT ENABLEDj, using omics selector
      feat_subsel = NULL,  # NOT ENABLED

      #DISABLED.
      #TODO: depricate all these
      group_action = NULL,
      data_source = NULL,
      plot_x = NULL,
      plot_y = NULL,
      plot_var_x = NULL,
      plot_var_y = NULL,
      plot_feats = NULL,
      # TODO:  make this conditional on what kind of data is loaded.. and just a single plot type
      raw_plot_type = NULL,
      comp_plot_type = NULL,
      GO = FALSE

    )
    ### OMICS  =========================================================

    #new_db_trig <- reactive( rv_data$trigger )
    rv_config <- reactive({
      rv_data$config
    })

    obs_sub <- mod_subset_selector_server("subset_selector_ui_obs",rv_config,"obs")
    var_sub <- mod_subset_selector_server("subset_selector_ui_var",rv_config,"var")


    all_omics <- reactive( rv_data$ad$var_names )  #only changes when new database is loaded
    def_omics <- reactive( rv_data$default$omics )

    # filter omics from subsetting
    active_omics <- reactive({
      #this is the maybe subsetting
      # subset var (omics)
      if (!is.null( var_sub$set ) ) {
        if (!is.null( var_sub$select )) {
          if (length(var_sub$select)>0) {
            return (all_omics()[ rv_data$ad$var[[ var_sub$set ]]  %in% var_sub$select ])
          } else {
            print("everything unselected...")
          }
        }
        }

        return( all_omics() )

    })


    selected_omics <- mod_omic_selector_server("omic_selector_ui_1", active_omics, def_omics) #, new_db_trig)


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


    # observeEvent(
    #   (rv_data$trigger>0),  #why does this happen twice?
    #   {
    #     req(rv_data$config,  # set when database is chosen
    #         rv_data$default)
    #     rv_data$trigger <- 0
    #     # update full omics list
    #     all_omics <- rv_data$omics
    #     #shinyjs::enable("RB_obs_X")
    #   },
    # ignoreNULL = TRUE,
    # ignoreInit = TRUE
    # ) #observe event



    ## dynamic subset UI group
    output$ui_data_layer <- renderUI({
      req(rv_data$config)
      choices <- rv_data$config[field=="layer"]$UI  # X, raw, or layers
      # default data_source is obs
      ret_tags <-  selectizeInput(ns("SI_data_layer"),
                           "choose data layer:",
                           choices =  choices,
                           selected = choices[1])

      return(ret_tags)
    })


    # dynamic x and y selector
    output$ui_meta_select <- renderUI({
      req(rv_data$config)

      group_obs <- rv_data$config[grp == TRUE & field=="obs"]$UI # <- choices_x
      group_obs2 <- rv_data$config[grp == TRUE & field=="obs"]$UI
      def_grp_o <- rv_data$default$obs_subset

      to_return <-  tagList(
        fluidRow(
          column(
            width=5,
            offset=0,
            #shinyjs::disabled(
            selectizeInput(ns("SI_group_obs"),
                           label = "group by:",
                           choices = group_obs,
                           selected = def_grp_o)
            #)
          ),
          column(
            width=2,
            offset=0,
            checkboxInput(ns("CB_sub_grp"),
                         label = "sub-grouping:",
                         value = FALSE)
          ),
          column(
            width=5,
            offset=0,
            shinyjs::disabled(
            selectizeInput(ns("SI_group_obs2"),
                           label = "sub-group",
                           choices = group_obs2,
                           selected = def_grp_o)
            )
          )

        ),
        hr(style = "border-top: 1px dashed grey;"),
        #uiOutput(ns("ui_obs_subset")),
        mod_subset_selector_ui(ns("subset_selector_ui_obs")),
        )

      return(to_return)

    })


    # dynamic x and y selector
    output$ui_omic_select <- renderUI({
      req(rv_data$config)

      group_var <- rv_data$config[grp == TRUE & field=="var"]$UI                                         #
      def_grp_v <- group_var[1]

      to_return <-  tagList(

          fluidRow(
            column(
              width = 5,
              offset = 0,
              #shinyjs::disabled(
              selectizeInput(ns("SI_group_var"),
                             label = "group by:",
                             choices = group_var,
                             selected = def_grp_v)
              #)
            )
          ),#fluidRow
          #uiOutput(ns("ui_var_subset"))
          hr(style = "border-top: 1px dashed grey;"),

          mod_subset_selector_ui(ns("subset_selector_ui_var"))


          ) #tagList

      return(to_return)

    })

    ### observe s =========================================================


    observe({
      req(input$SI_group_obs2)
      if ( input$CB_sub_grp ) {
        shinyjs::enable("SI_group_obs2")
      } else {
        shinyjs::disable("SI_group_obs2")
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


      observeEvent(input$AB_plot_now, {
        # send signal to to update the selected omics
        selected_omics$freeze <- selected_omics$freeze+1
      })

      observeEvent(selected_omics$freeze, {
          # route the chosen data type out...
          # req(input$CB_obs_subsel,
          #     input$SI_obs_subset)
          #TODO: add data_source to config values
          #

        rv_selections$data_layer <- input$SI_data_layer

        rv_selections$observ_subset <- obs_sub$set #input$SI_obs_subset
        rv_selections$observ_subsel <- obs_sub$select #input$CB_obs_subsel
        # #
        rv_selections$feat_subset <- var_sub$set #input$SI_var_subset #NA # NOT ENABLEDj, using omics selector
        rv_selections$feat_subsel <- var_sub$select #input$CB_var_subsel #NA # NOT ENABLED

        # group (plotting)
        rv_selections$feat_group_by <- input$SI_group_var
        rv_selections$observ_group_by <- input$SI_group_obs
        rv_selections$observ_group_by2 <- input$SI_group_obs2  # could be null

        rv_selections$selected_omics <- selected_omics  # value & viz_now & all_active

        rv_selections$GO = TRUE  # use this to trigger the heatmap to update...


        # DEPRICATE
        rv_selections$data_source <- NA # "X" #DISABLED input$RB_obs_X
        rv_selections$plot_x <- NA # DISABLED input$SI_x
        rv_selections$plot_y <- NA # DISABLED....input$SI_y
        rv_selections$group_action <- NA #DISABLED input$RB_none_or_grp
        rv_selections$plot_var_x <- NA # NOT ENABLED#input$SI_var_x
        rv_selections$plot_var_y <- NA # NOT ENABLED input$SI_var_y
        rv_selections$plot_feats <- NA # NOT ENABLED input$CB_xy_var_select


    })


  ### RETURN =========================================================
  return(rv_selections)

  })
}

## To be copied in the UI
# mod_side_selector_ui("side_selector_ui_1")

## To be copied in the server
# mod_side_selector_server("side_selector_ui_1")

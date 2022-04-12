#' curation_lipid UI Function
#'
#' @description A shiny Module for the curation of lipid data sets.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @author Rico Derks
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#'
mod_curation_lipid_ui <- function(id){
  ns <- NS(id)
  tagList(
    h3("Lipidomics"),
    HTML("<p>
         For the curation of a new data set three files are needed. Please supply
         the following files:
         <ul>
         <li>File (csv) containing the data matrix as <b>n x m</b>: <b>n</b> = observations/samples,
         <b>m</b> = variables/lipids. File should include a first row containing
         lipid names, and a first column with sample id's. </li>
         <li>File (csv) containing variable information.</li>
         <li>File (csv) containing sample information.</li>
         </ul>
         </p>"),
    br(),
    fluidRow(
      column(width = 6,
             style = "border-right: 1px solid black",
             h4("File upload :"),
             ### upload files
             # data matrix
             fileInput(
               inputId = ns("fi_lipid_data"),
               label = "Lipid data matrix (csv) :",
               multiple = FALSE,
               accept = ".csv",
               buttonLabel = "Browse...",
               placeholder = "No file selected"
             ),
             htmlOutput(outputId = ns("fi_lipid_data_status")),
             # sample info
             fileInput(
               inputId = ns("fi_lipid_obs_info"),
               label = "Sample info (csv) :",
               multiple = FALSE,
               accept = ".csv",
               buttonLabel = "Browse...",
               placeholder = "No file selected"
             ),
             htmlOutput(outputId = ns("fi_lipid_obs_status")),
             # curate the data
             actionButton(inputId = ns("ab_lipid_curate"),
                          label = "Curate",
                          class = "btn btn-large btn-danger")
      ),
      column(width = 6,
             # show some parameters if all files are correctly loaded
             htmlOutput(outputId = ns("ui_lipid_curation_params"))
      )
    )
  )
}

#' curation_lipid Server Functions
#'
#' @importFrom tools file_ext
#' @importFrom waiter Waiter transparent spin_loaders
#'
#' @noRd
#'
mod_curation_lipid_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    #### get the database root and all names ####
    DB_ROOT <- golem::get_golem_options("db_root")
    DB_NAMES <- golem::get_golem_options("database_names")

    ### define some reactive values
    rv_data <- reactiveValues(data = list(data = NULL,
                                          status = NULL,
                                          message = NULL),
                              observations = list(data = NULL,
                                                  status = NULL,
                                                  message = NULL),
                              variables = list(data = NULL,
                                               status = NULL,
                                               message = NULL),
                              database = list(name = NULL,
                                              status = NULL,
                                              message = NULL),
                              test = list(tests = NULL,
                                          status = NULL,
                                          message = NULL))

    #### get the data matrix ####
    observeEvent(input$fi_lipid_data, {
      req(input$fi_lipid_data)

      # validate the extension, message is NOT shown!!
      ext <- tools::file_ext(input$fi_lipid_data$name)
      validate(need(ext == "csv", "Please upload a csv file"))

      # read the file
      res <- read_lipid_data(filename = input$fi_lipid_data$datapath,
                             data_type = "data_matrix")

      # get all everything
      rv_data$data$status <- res$status
      rv_data$data$message <- res$message
      rv_data$data$data <- res$data_df

      # get the lipid info
      res_var <- get_lipid_info(lipid_data = rv_data$data$data)
      # get all everything
      # rv_data$variables$status <- res_var$status
      # rv_data$variables$message <- res_var$message
      rv_data$variables$data <- res_var
    })

    # show if there is an error during upload
    output$fi_lipid_data_status <- renderUI({
      req(rv_data$data$status)

      # show if there is an error
      if(rv_data$data$status == "error"){
        HTML(paste("<text style='color:red; font-weight:bold'>Error:", rv_data$data$message, "</text>"))
      }
    })


    #### get sample info ####
    observeEvent(input$fi_lipid_obs_info, {
      req(input$fi_lipid_obs_info)

      # validate the extension, message is NOT shown!!
      ext <- tools::file_ext(input$fi_lipid_obs_info$name)
      validate(need(ext == "csv", "Please upload a csv file"))

      # read the file
      res <- read_lipid_data(filename = input$fi_lipid_obs_info$datapath,
                             data_type = "observations")

      # get all everything
      rv_data$observations$status <- res$status
      rv_data$observations$message <- res$message
      rv_data$observations$data <- res$data_df
    })

    # show if there is an error during upload
    output$fi_lipid_obs_status <- renderUI({
      req(rv_data$observations$status)

      # show if there is an error
      if(rv_data$observations$status == "error"){
        HTML(paste("<text style='color:red; font-weight:bold'>Error:", rv_data$observations$message, "</text>"))
      }
    })


    #### Curation parameters ####
    output$ui_lipid_curation_params <- renderUI({
      # if one of them is still NULL don't continue
      req(rv_data$data$status,
          rv_data$observations$status,
          rv_data$observations$data)

      # if there is any error do NOT continue
      status <- c(rv_data$data$status,
                  # rv_data$variables$status,
                  rv_data$observations$status)

      # if no error
      if(!("error" %in% status)){
        tagList(
          h4("Parameters :"),
          # input field for database name
          textInput(inputId = ns("ti_lipid_db_name"),
                    label = "Enter database name :"),
          # show any possible error messages
          htmlOutput(outputId = ns("lipid_db_name_status")),
          # select tests for differential expression
          checkboxGroupInput(inputId = ns("cbg_lipid_tests_de"),
                             label = "Select tests for differential expression :",
                             choices = c("t-test" = "ttest",
                                         "Mann-Whitney" = "mw"),
                             selected = "ttest"),
          # show any possible error messages
          htmlOutput(outputId = ns("lipid_test_status")),
          # select which group to use for the differential expression
          selectInput(inputId = ns("si_lipid_group_de"),
                      label = "Select a group for differential expression :",
                      choices = colnames(rv_data$observations$data)[-1],
                      multiple = FALSE),
          # select curation steps
          checkboxGroupInput(inputId = ns("cbg_lipid_curation_params"),
                             label = "Select curation steps :",
                             choices = c("Remove zero lipids" = "zero",
                                         "Remove lipids 2/3 NA" = "twothird"),
                             selected = c("zero", "twothird")),
          # show when curation is done
          htmlOutput(outputId = ns("lipid_curation_status"))
        )
      }
    })

    # show when curation is done
    output$lipid_curation_status <- renderUI({
      req(rv_data$curation_status)

      if(rv_data$curation_status == "Done"){
        HTML(paste("<text style='color:blue; font-weight:bold'>Curation done!</text>"))
      }
    })


    #### define waiter ####
    wtr <- waiter::Waiter$new(
      # set the waiter to the full app
      id = NULL,
      html = waiter::spin_whirly(),
      color = waiter::transparent(.5)
    )


    #### curation button clicked ####
    observeEvent(input$ab_lipid_curate, {
      # if one of them is still NULL don't continue
      req(rv_data$data$status,
          rv_data$observations$status,
          input$ti_lipid_db_name,
          input$si_lipid_group_de)

      # show the waiter
      wtr$show()

      ## initialize some variables here
      # initialize curation parameters
      remove_zero_lipids <- FALSE
      remove_twothird_lipids <- FALSE

      ## get going
      # get the parameters for the curation
      lipid_curation_params <- input$cbg_lipid_curation_params
      if("zero" %in% lipid_curation_params) {
        remove_zero_lipids <- TRUE
      }
      if("twothird" %in% lipid_curation_params) {
        remove_twothird_lipids <- TRUE
      }

      # get which tests to do for the differential expression table
      if(is.null(input$cbg_lipid_tests_de)) {
        # no test selected
        rv_data$test$status <- "error"
        rv_data$test$message <- "No test selected!!"
      } else {
        rv_data$test$tests <- input$cbg_lipid_tests_de
      }

      # get the group name selected for differential expression
      group_name <- input$si_lipid_group_de

      # get the new database name
      rv_data$database$name <- input$ti_lipid_db_name

      # check if the database name is valid
      status_db_name <- check_database_name(db_name = rv_data$database$name,
                                            current_db_names = DB_NAMES)
      # update everything
      rv_data$database$status <- status_db_name$status
      rv_data$database$message <- status_db_name$message

      ## check for errors
      # get all status'
      status <- c(rv_data$data$status,
                  rv_data$observations$status,
                  status_db_name$status,
                  rv_data$test$status)

      if("error" %in% status) {
        # do nothing

      } else {
        # continue
        # do the curation
        res_cur <- curate_lipidomics(data = list(data = rv_data$data$data,
                                                 obs = rv_data$observations$data,
                                                 var = rv_data$variables$data),
                                     db_name = rv_data$database$name,
                                     db_root = DB_ROOT,
                                     remove_zero_lipids = remove_zero_lipids,
                                     remove_twothird_lipids = remove_twothird_lipids,
                                     tests = rv_data$test$tests,
                                     test_group = group_name)

        rv_data$curation_status <- res_cur
      } # end of all fine

      # stop the waiter
      wtr$hide()
    })


    # show if there is an error in the database name
    output$lipid_db_name_status <- renderUI({
      req(rv_data$database$status)

      # show if there is an error
      if(rv_data$database$status == "error"){
        HTML(paste("<text style='color:red; font-weight:bold'>Error:", rv_data$database$message, "</text>"))
      }
    })


    # show error if no test is selected
    output$lipid_test_status <- renderUI({
      req(rv_data$test$status)

      # show if there is an error
      if(rv_data$test$status == "error"){
        HTML(paste("<text style='color:red; font-weight:bold'>Error:", rv_data$test$message, "</text>"))
      }
    })

  }) # end mod_curation_lipid_server
}

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
    # Variable / lipid information
    fileInput(
      inputId = ns("fi_lipid_var_info"),
      label = "Lipid info (csv) :",
      multiple = FALSE,
      accept = ".csv",
      buttonLabel = "Browse...",
      placeholder = "No file selected"
    ),
    htmlOutput(outputId = ns("fi_lipid_var_status")),
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

  )
}

#' curation_lipid Server Functions
#'
#' @importFrom tools file_ext
#'
#' @noRd
#'
mod_curation_lipid_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    ### define some reactive values
    rv_data <- reactiveValues(data = list(data = NULL,
                                          status = NULL,
                                          message = NULL),
                              observations = list(data = NULL,
                                                  status = NULL,
                                                  message = NULL),
                              variables = list(data = NULL,
                                               status = NULL,
                                               message = NULL))

    ### get the data matrix ###
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
    })

    # show if there is an error during upload
    output$fi_lipid_data_status <- renderUI({
      req(rv_data$data$status)

      if(rv_data$data$status == "error"){
        HTML(paste("<text style='color:red; font-weight:bold'>Error:", rv_data$data$message, "</text>"))
      }
    })


    ### get variable info ###
    observeEvent(input$fi_lipid_var_info, {
      req(input$fi_lipid_var_info)

      # validate the extension, message is NOT shown!!
      ext <- tools::file_ext(input$fi_lipid_var_info$name)
      validate(need(ext == "csv", "Please upload a csv file"))

      # read the file
      res <- read_lipid_data(filename = input$fi_lipid_var_info$datapath,
                             data_type = "variables")

      # get all everything
      rv_data$variables$status <- res$status
      rv_data$variables$message <- res$message
      rv_data$variables$data <- res$data_df
    })

    # show if there is an error during upload
    output$fi_lipid_var_status <- renderUI({
      req(rv_data$variables$status)

      if(rv_data$variables$status == "error"){
        HTML(paste("<text style='color:red; font-weight:bold'>Error:", rv_data$variables$message, "</text>"))
      }
    })


    ### get sample info ###
    observeEvent(input$fi_lipid_obs_info, {
      req(input$fi_lipid_obs_info)

      # validate the extension, message is NOT shown!!
      ext <- tools::file_ext(input$fi_lipid_sample_info$name)
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

      if(rv_data$observations$status == "error"){
        HTML(paste("<text style='color:red; font-weight:bold'>Error:", rv_data$observations$message, "</text>"))
      }
    })
  })
}

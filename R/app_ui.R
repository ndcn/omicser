#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),

    # golem::activate_js(), #already loaded in your golem by `bundle_resources()`

    # Your application UI logic
    fluidPage(
      shinyjs::useShinyjs(),

      titlePanel(
        fluidRow(
          col_4(
            h1("NDCN omicser"),
            h5("Browse and play for creative hypothesis generation")
          ),
          col_4(img(src = "www/logo.svg"))
        )
      ), # end titlePanel
      sidebarLayout(
        div(
          id = "cond_sidebar",
          sidebarPanel(
            width = 4,
            conditionalPanel(
              ' input.top_tab === "welcome" || input.top_tab === "add_info" || input.top_tab === "etc" ',
              mod_side_info_ui("side_info_ui_1")
            ),
            conditionalPanel(
              ' input.top_tab === "playground" || input.top_tab === "table" || input.top_tab === "ingest" || input.top_tab === "export" ',
              mod_side_selector_ui("side_selector_ui_1")
            )
          ) # sidebarpanel
        ), # div
        mainPanel(
          width = 8,
          tabsetPanel(
            type = "tabs", # pills look good
            id = "top_tab",
            tabPanel(
              title = "Welcome", value = "welcome",
              mod_welcome_ui(id = "welcome_ui_1")
            ),
            # ingest tab
            tabPanel(
              title = "Ingest", value = "ingest",
              mod_ingestor_ui(id = "ingestor_ui_1")
            ),
            # playground tab
            tabPanel(
              title = "Playground", value = "playground",
              mod_playground_ui(id = "playground_ui_1")
            ),
            # Additional Info tab
            tabPanel(
              title = "Additional Info", value = "add_info",
              mod_additional_info_ui(id = "additional_info_ui_1")
            ),
            # # Export tab
            # tabPanel(
            #   title = "Export", value = "export",
            #   # copy the landing module for now
            #   mod_export_ui(id = "export_ui_1")
            # ),
            # table tab
            tabPanel(
              title = "Data Table", value = "table",
              # DT::dataTableOutput("my_datatable_0")
              mod_table_ui(id = "pg_table_ui_1")
            ),
            # Etc tab
            tabPanel(
              title = "Etc", value = "etc",
              # copy the landing module for now
              mod_additional_info_ui(id = "additional_info_ui_1")
            )
          ) # tabsetpanel
        ) # mainpanel
      ), # end sidebarlayout
      # actionButton("alert", "xxx"),
      tags$footer(tags$div(
        class = "footer", checked = NA, HTML('
              <head>
              <style>
              .footer a:link {color: #008b42; background-color: transparent; text-decoration: none}
              .footer a:visited {color: #008b42; background-color: transparent; text-decoration: none}
              .footer a:hover {color: #008b42; background-color: transparent; text-decoration: underline;}
              .footer a:active {color: red; background-color: transparent; text-decoration: underline;}
              .footer div {padding: 0px 0px 10px; color: grey;}
              position:fixed;
              bottom:0;
              width:100%;
              height:50px;   /* Height of the footer */
              /*padding: 10px;*/
              </style>
              </head>

              <body>
              <hr>
              <div>
              <a href="https://chanzuckerberg.com/ndcn/" target="_blank">CZI NDCN</a>

              Shiny App credits: NDCN, DTI, Andy Henrie <a href="https://github.com/ergonyc/omicser" target="_blank">ergonyc/omicser@github</a>

              </div>
              </body>
              '),
        align = "left"
      ))
    ) # end fluidpage
  ) # end taglist
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www", app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "omicser"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}

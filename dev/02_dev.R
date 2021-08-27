# Building a Prod-Ready, Robust Shiny Application.
#
# README: each step of the dev files is optional, and you don't have to
# fill every dev scripts before getting started.
# 01_start.R should be filled at start.
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
#
#
###################################
#### CURRENT FILE: DEV SCRIPT #####
###################################

# Engineering

## Dependencies ----
## Add one line by package you want to add as dependency
usethis::use_package( "thinkr" )
# add more?
usethis::use_package("shiny")

#usethis::use_package("shinyWidgets")
usethis::use_package("ggplot2")
usethis::use_package("gridExtra")

usethis::use_package("ggdendro")
usethis::use_package("DT")
usethis::use_package("dplyr")

usethis::use_package("plotly")
usethis::use_package("Matrix")

usethis::use_package("anndata")
usethis::use_package("SingleCellExperiment")

usethis::use_package("tidyr") # called by dplyr as dependency?
# usethis::use_package("magrittr")
# usethis::use_package("purrr")
# usethis::use_package("tibble")

# usethis::use_package("xml2")
# usethis::use_package("readr")
# usethis::use_package("stringr")
# usethis::use_package("shinyWidgets")
# usethis::use_package("waiter")
# usethis::use_package("vroom")
# usethis::use_package("readxl")
# usethis::use_package("yaml")
# usethis::use_package("rclipboard")
# usethis::use_package("markdown")
## Add modules ----
## Create a module infrastructure in R/
golem::add_module( name = "ingestor" ) # Name of the module
golem::add_module( name = "export" ) # Name of the module
golem::add_module( name = "playground" ) # Name of the module
golem::add_module( name = "pg_vis_raw" ) # Name of the module
golem::add_module( name = "pg_vis_comp" ) # Name of the module
golem::add_module( name = "pg_table" ) # Name of the module
golem::add_module( name = "side_selector" ) # Name of the module
golem::add_module( name = "welcome" ) # Name of the module
golem::add_module( name = "side_info" ) # Name of the module
golem::add_module( name = "additional_info") # Name of the module
golem::add_module( name = "omic_selector") # Name of the module
golem::add_module( name = "pg_vis_qc") # Name of the module

## Add helper functions ----
## Creates fct_* and utils_*
# golem::add_fct( "module" ) # sever
# golem::add_utils( "module" ) #ui
#   e.g.  boxplot, heatmap, functions
golem::add_fct( "pg_vis_raw" ) # sever
#  e.g. volcano plot functions..
golem::add_fct( "pg_vis_comp" ) # sever

golem::add_fct( "ingestor" ) # sever


# ## External resources
# TODO: ADD THESE JS/CSS resources
# ## Creates .js and .css files at inst/app/www
# golem::add_js_file( "script" )
# golem::add_js_handler( "handlers" )
# golem::add_css_file( "custom" )

## Add internal datasets ----
## If you have data in your package
usethis::use_data_raw( name = "transcript_gold_standard", open = FALSE )
usethis::use_data_raw( name = "metabol_gold_standard", open = FALSE )
usethis::use_data_raw( name = "proteo_gold_standard", open = FALSE )
usethis::use_data_raw( name = "lipid_gold_standard", open = FALSE )
usethis::use_data_raw( name = "vilas_A", open = FALSE )
usethis::use_data_raw( name = "vilas_B", open = FALSE )
usethis::use_data_raw( name = "domenico_A", open = FALSE )
usethis::use_data_raw( name = "domenico_B", open = FALSE )
usethis::use_data_raw( name = "oscar_microglia", open = FALSE )
usethis::use_data_raw( name = "oscar_toy", open = FALSE )
usethis::use_data_raw( name = "yassene_example", open = FALSE )

## Tests ----
## Add one line by test you want to create
usethis::use_test( "app" )

# Documentation

## Vignette ----
usethis::use_vignette("omicser")
devtools::build_vignettes()

# ## Code Coverage----
# ## Set the code coverage service ("codecov" or "coveralls")
# usethis::use_coverage()
#
# # Create a summary readme for the testthat subdirectory
# covrpage::covrpage()
#
# ## CI ----
# ## Use this part of the script if you need to set up a CI
# ## service for your application
# ##
# ## (You'll need GitHub there)
# usethis::use_github()
#
# # # GitHub Actions
# usethis::use_github_action()
# # Chose one of the three
# # See https://usethis.r-lib.org/reference/use_github_action.html
# usethis::use_github_action_check_release()
# usethis::use_github_action_check_standard()
# usethis::use_github_action_check_full()
# # Add action for PR
# usethis::use_github_action_pr_commands()
#
# # Travis CI
# usethis::use_travis()
# usethis::use_travis_badge()
#
# # AppVeyor
# usethis::use_appveyor()
# usethis::use_appveyor_badge()
#
# # Circle CI
# usethis::use_circleci()
# usethis::use_circleci_badge()
#
# # Jenkins
# usethis::use_jenkins()
#
# # GitLab CI
# usethis::use_gitlab_ci()

# You're now set! ----
# go to dev/03_deploy.R
rstudioapi::navigateToFile("dev/03_deploy.R")


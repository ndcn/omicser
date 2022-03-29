#' curation_lipid
#'
#'
#' @title Read the data for lipid curation.
#'
#' @description A helper function to read the different data files for lipidomics
#'     data set curration.
#'
#' @param filename file name (csv) of the data matrix file
#' @param data_type what is the data type to read, i.e. 'data_marix', 'variables'
#'    or 'observations'.
#'
#' @return Returns a list with status, message and data as data.frame.
#'
#' @author Rico Derks
#'
#' @noRd
#'
read_lipid_data <- function(filename = NULL,
                            data_type = c("data_matrix", "variables", "observations")) {
  ### sanity check
  # is a filename supplied
  if(is.null(filename)) {
    stop("Please supply a file name!")
  }
  # does the file exist
  if(!file.exists(filename)) {
    stop("File doesn't exist!")
  }
  # check if data_type is correct
  if(!(data_type %in% c("data_matrix", "variables", "observations"))) {
    stop("Wrong 'data_type' supplied!")
  }

  # initialize list to return
  res <- list(status = NULL,
              message = NULL,
              data_df = NULL)

  # set check names parameter
  check_names <- ifelse(data_type == "data_matrix",
                        FALSE,
                        TRUE)

  # read the file
  data_df <- read.csv(file = filename,
                      header = TRUE,
                      na.strings = c(".", "", "NA"),
                      check.names = check_names)

  # what column(s) should be present.
  col_present <- switch(
    data_type,
    "data_matrix" = "sample_id",
    "variables" = "lipid_name",
    "observations" = "sample_id"
  )

  # check if the correct column exists
  if(!(col_present %in% colnames(data_df))) {
    res$status <- "error"
    res$message <- paste0("'", col_present ,"' column missing!")
  } else {
    res$status <- "ok"
    res$data_df <- data_df
  }

  # return the result
  return(res)
}


#' @title Check the database name
#'
#' @description Check if the database name doesn't contain any spaces or special
#'     characters.
#'
#' @param db_name character vector (length = 1) with the database name.
#' @param current_db_names named list with all the current database names.
#'
#' @return TRUE if database name is ok, FALSE if not
#'
#' @author Rico Derks
#'
#' @noRd
#'
check_database_name <- function(db_name = NULL,
                                current_db_names = NULL) {

  # initialize the list to return
  res <- list(status = NULL,
              message = NULL)

  # check if the database doesn't contain any special characters
  db_name_matches <- regexpr(text = db_name,
                             pattern = "[a-zA-Z//-_]+")

  if(nchar(db_name) == attributes(db_name_matches)$match.length) {
    # all ok
    res$status <- "ok"
  } else {
    # special character used
    res$status <- "error"
    res$message <- "Error in database name. Please only use 'a-z', 'A-Z' or '-_'!"
  }

  # check if database name exists
  # this is only done for database known to omicser
  # other folders are not detected
  if(db_name %in% unlist(current_db_names)) {
    res$status <- "error"
    res$message <- "Database name already exists!"
  }

  # return the result
  return(res)
}

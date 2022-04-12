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
#' @importFrom utils read.csv
#'
#' @noRd
#'
read_lipid_data <- function(filename = NULL,
                            data_type = c("data_matrix", "variables", "observations")) {
  ### sanity check
  # is a file name supplied
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


#' @title Extract all lipid information
#'
#' @description Extract which lipids are uploaded and combine with the data
#'     from the master lipids tables.
#'
#' @param lipid_data data.frame with the lipid data. Colnames contain the lipid
#'     names from the lipidyzer.
#'
#' @return data.frame containing only the lipids uploaded and with all extra
#'     information from master lipids table.
#'
#' @author Rico Derks
#'
#' @importFrom dplyr filter group_by summarise rename
#' @importFrom rlang .data
#'
#' @noRd
#'
get_lipid_info <- function(lipid_data = NULL) {
  # Sanity check
  if(is.null(lipid_data)) {
    stop("No data supplied!")
  }

  # get the lipid names (lipidyzer) from the uploaded data.
  lipidyzer_lipids <- colnames(lipid_data)

  # get all information about the lipids
  lipid_info <- master_table_lipids %>%
    filter(.data$lipidyzer_name %in% lipidyzer_lipids)

  # data.frame needs to be collapsed, for some lipids multiple entries
  lipid_info_clean <- lipid_info %>%
    group_by(.data$lipidyzer_name) %>%
    summarise(common_name = .data$common_name[1],
              species_shorthand = .data$species_shorthand[1],
              lipid_class = .data$lipid_class[1],
              abbreviation = .data$abbreviation[1],
              abbreviation_chains = .data$abbreviation_chains[1],
              regno = paste(.data$regno, collapse = ","),
              lm_id = paste(.data$lm_id, collapse = ","),
              name = paste(.data$name, collapse = "\t"),
              abbrev = .data$abbrev[1],
              abbrev_chains = .data$abbrev_chains[1],
              core = .data$core[1],
              main_class = .data$main_class[1],
              sub_class = .data$sub_class[1],
              class_level4 = .data$class_level4[1],
              kegg_id = paste(.data$kegg_id, collapse = ","),
              hmdb_id = paste(.data$hmdb_id, collapse = ","),
              chebi_id = paste(.data$chebi_id, collapse = ","),
              lipidbank_id = paste(.data$lipidbank_id, collapse = ","),
              pubchem_id = paste(.data$pubchem_cid, collapse = ","),
              class_lm_id = .data$class_lm_id[1],
              class_kegg_id = .data$class_kegg_id[1],
              class_hmdb_id = .data$class_hmdb_id[1],
              class_chebi_id = .data$class_chebi_id[1],
              class_lipidbank_id = .data$class_lipidbank_id[1],
              class_pubchem_id = .data$class_pubchem_cid[1],
              .groups = "drop") %>%
    rename(lipid_name = .data$lipidyzer_name) %>%
    as.data.frame()

  # set the row names
  rownames(lipid_info_clean) <- lipid_info_clean$lipid_name

  # return results
  return(lipid_info_clean)
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


#' @title Curate a lipidomics data set
#'
#' @description Curate a lipidomics data set generated with the Sciex Lipidzyer.
#'     The resulting {anndata} object is written to a file in hdf5 format.
#'
#' @param data list containing all data, i.e. data, observations info and
#'     variables info.
#' @param db_name name of the new database.
#' @param db_root root of all the databases.
#' @param steps which steps should be done.
#' @param remove_zero_lipids remove lipids which only contain NA's or zero's.
#'     Default is FALSE.
#' @param remove_twothird_lipids remove lipids which contain more than 2/3 NA's
#'     or zero's. Default is FALSE.
#' @param tests which tests to use for the differential expression table.
#'     t-test ('ttest') and or Mann-Whitney test ('mw').
#' @param test_group group on which to do the test(s).
#'
#' @return Nothing.
#'
#' @author Rico Derks
#'
#' @noRd
#'
curate_lipidomics <- function(data = NULL,
                              db_name = NULL,
                              db_root = NULL,
                              steps = NULL,
                              remove_zero_lipids = FALSE,
                              remove_twothird_lipids = FALSE,
                              tests = c("ttest", "mw"),
                              test_group = NULL) {
  #### sanity checks
  # is/are the correct tests selected
  tests <- match.arg(arg = tests,
                     several.ok = TRUE)

  # is a group selected
  if(is.null(test_group)) {
    stop("No group select for the test(s)!")
  }

  # check if the database folder exist, if not create it
  check <- create_db_folder(db_root = db_root,
                            db_name = db_name)

  # if creating db folder fails
  if(check == FALSE) {
    stop("New database folder doesn't exist and can NOT be created! Please check
         write permissions!")
  }

  #### Get going here
  # initialize list with all the data
  lipid_data <- list(data_matrix = NULL,
                     obs_info = NULL,
                     var_info = NULL)

  #### Get all data ####
  # lipid data: remove first column and convert to matrix
  lipid_data$data_matrix <- as.matrix(data$data[, -1])
  rownames(lipid_data$data_matrix) <- data$data$sample_id

  # sample info
  lipid_data$obs_info <- data$obs
  rownames(lipid_data$obs_info) <- lipid_data$obs_info$sample_id

  # variable info
  lipid_data$var_info <- data$var
  rownames(lipid_data$var_info) <- lipid_data$var_info$lipid_name

  #### Do some clean up ####
  # remove lipds which only contain NA's or zero's
  if (remove_zero_lipids == TRUE) {
    lipid_data <- clean_up(lipid_data = lipid_data)
  }

  #### Pre-processing ####
  lipid_data <- preproc_lipid_data(lipid_data = lipid_data,
                                   remove_twothird_lipids = remove_twothird_lipids)

  #### Pack into ann-data ####
  ad <- pack_lipid_anndata(data = lipid_data,
                           db_root = db_root,
                           db_name = db_name)

  #### Calculate differential expression table ####
  diff_exp_table <- diff_epxression(data = lipid_data,
                                    test_types = tests,
                                    test_group = test_group,
                                    db_root = db_root,
                                    db_name = db_name)

  #### write config file ####
  write_lipid_config(data = lipid_data,
                     diff_exp = diff_exp_table,
                     db_root = db_root,
                     db_name = db_name)

  # write to hdf5 file
  ad$write_h5ad(filename = file.path(db_root, db_name, "core_data.h5ad"))
  ad$write_h5ad(filename = file.path(db_root, db_name, "db_data.h5ad"))

  # return result
  return("Done")
}


#' @title Write the database config file
#'
#' @param data list with all the data
#' @param diff_exp the differential expressio data
#' @param db_root root path of all databases
#' @param db_name name of the new database
#'
#' @return Nothing
#'
#' @author Rico Derks
#'
#' @noRd
#'
write_lipid_config <- function(data = NULL,
                               diff_exp = NULL,
                               db_root = NULL,
                               db_name = NULL) {
  # define the type
  omic_type <- "lipid"
  aggregate_by_default <- ifelse(omic_type == "transcript", TRUE, FALSE) #e.g.  single cell

  # define the config list
  config_list <- list(
    ### grouping factors
    # if it needs to be in subset add here as well
    group_vars = colnames(data$var_info)[-1],
    group_obs = c("Group"),

    ### layer info
    layer_values = c("X"),
    # are the names of the layers used?
    layer_names = c("Conc."),

    # ANNOTATIONS / TARGETS
    # what adata$obs do we want to make default values for...
    # should just pack according to UI?

    ### observations
    default_obs = colnames(data$obs_info[1]),
    # heatmap default selected
    obs_annots = colnames(data$obs_info[-1]),

    ### variables
    default_var = c("feature_name"),
    # heatmap default selected
    var_annots = colnames(data$var_info)[-1],

    ### set the target features, looks like they are not loaded yet
    target_features = "",

    ### set the feature details when dot clicked in volcano plot
    # looks like this is not working, in Domenico script it works
    feature_deets = c(""),

    ### differential expression
    diffs = list(diff_exp_comps = levels(factor(diff_exp$versus)),
                 diff_exp_comp_type = levels(factor(diff_exp$comp_type)), #i don"t think we need this
                 diff_exp_obs_name = levels(factor(diff_exp$obs_name)),
                 diff_exp_tests = levels(factor(diff_exp$test_type))),

    ### meta info
    annotation_database =  NA,
    publication = "TBD",
    method = "bulk", # c("single-cell","bulk","other")
    omic_type = omic_type, # see above
    aggregate_by_default = aggregate_by_default, # see above
    organism = 'mmusculus',
    lab = "Giera",
    title = "Lipidomics",
    date = format(Sys.time(), "%a %b %d %X %Y")
  ) # end config list

  # write configuration to yaml file
  omicser::write_db_conf(config_list = config_list,
                         db_name = db_name,
                         db_root = db_root)

  # return nothing
  return(NULL)
}


#' @title Pack the lipidomics data into an anndata object
#'
#' @param data a list with all the data
#' @param db_root root path of all databases
#' @param db_name name of the database
#'
#' @return anndata object.
#'
#' @author Rico Derks
#'
#' @noRd
#'
pack_lipid_anndata <- function(data = NULL,
                               db_root = NULL,
                               db_name = NULL) {
  # prepare the data
  data_list <- list(
    data_mat = data$data_matrix,
    obs_meta = data$obs_info,
    var_annot = data$var_info,
    omics = NULL,
    sample_ID = data$obs_info$sample_id,
    etc = NULL,
    raw = NULL,
    uns = NULL
  )

  # make anndata object
  ad <- setup_database(database_name = db_name,
                       db_path = db_root,
                       data_in = data_list,
                       re_pack = TRUE)

  # return anndata object
  return(ad)
}


#' @title Calculate the differential expression table
#'
#' @param data list with all the data.
#' @param test_types what test do you want to do.
#' @param test_group group on which to do the test(s).
#' @param comp the comparissons.
#' @param db_root root path of all databases.
#' @param db_name name of the database.
#'
#' @return The differential expression results. The differential expression
#'     results are also save to a RDS file
#'
#' @author Rico Derks
#'
#' @importFrom reticulate import
#' @importFrom utils combn
#'
#' @noRd
#'
diff_epxression <- function(data = NULL,
                            test_types = c("ttest", "mw"),
                            test_group = NULL,
                            comp = NULL,
                            db_name = NULL,
                            db_root = NULL) {
  #### Sanity checks
  # are the correct tests selected
  test_types <- match.arg(arg = test_types,
                          several.ok = TRUE)

  if(is.null(test_group)) {
    stop("No group select for the test(s)!")
  }

  #### Get going
  # Create the comparisons
  # get the names
  comp_names <- sort(unique(data$obs_info[, test_group]))
  # create the combinations
  comp_types <- combn(x = comp_names,
                      m = 2)

  # initialize data.frame
  de_table <- data.frame()

  for(test_type in test_types) {
    # each test
    for(comp_type in seq(ncol(comp_types))) {
      # each comparison
      # extract the x and y of the comparison
      # comp_type <- comp_types[1]
      parts <- unlist(
        regmatches(x = comp_type,
                   m = regexec(text = comp_type,
                               pattern = "^\\{?([a-zA-Z0-9\\/\\+\\._]+)\\}?V\\{?([a-zA-Z0-9\\/\\+\\._]+)\\}?$"))
      )
      reference <- comp_types[2, comp_type]
      group <- comp_types[1, comp_type]

      # make a data.frame
      lipid_df <- cbind(data$obs_info[, c("sample_id", test_group)], as.data.frame(data$data_matrix))

      if(reference == "rest") {
        # if the reference is all other groups
      } else {
        # compute the de table
        tmp_de_table <- switch(
          test_type,
          "ttest" = {
            perform_ttest(data = lipid_df,
                          group = group,
                          reference = reference,
                          obs_name = test_group)
          },
          "mw" = {
            perform_mwtest(data = lipid_df,
                           group = group,
                           reference = reference,
                           obs_name = test_group)
          }
        )
        # combine the results
        de_table <- rbind(de_table, tmp_de_table)
      }
    } # end for-loop comp_types
  } # end for-loop test_types

  # make some column factors
  de_table$versus <- as.factor(de_table$versus)
  de_table$comp_type <- as.factor(de_table$comp_type)
  de_table$obs_name <- as.factor(de_table$obs_name)
  de_table$test_type <- as.factor(de_table$test_type)

  # # remove infinite values and NaN's
  # remove_idx <- c(which(is.infinite(diff_exp$logfoldchanges)),
  #                 which(is.nan(diff_exp$logfoldchanges)))
  # diff_exp <- diff_exp[-remove_idx, ]

  # save the table
  saveRDS(object = de_table,
          file = file.path(db_root, db_name, "db_de_table.rds"))

  # return de table
  return(de_table)
}


#' @title Perform statistical test for differential expression table
#'
#' @description Do a statistical test, i.e. t-test or Mann-Whitney test.
#'
#' @param data the data
#' @param group name of the group
#' @param reference name of the reference
#' @param obs_name column name which contains the group
#'
#' @return
#'
#' @author Rico Derks
#'
#' @importFrom dplyr filter mutate select rename rename_with
#' @importFrom tidyr pivot_longer nest unnest
#' @importFrom purrr map
#' @importFrom rlang !!
#' @importFrom broom tidy
#' @importFrom stats t.test p.adjust
#' @importFrom rlang .data
#'
#' @noRd
#'
perform_ttest <- function(data = NULL,
                          group = NULL,
                          reference = NULL,
                          obs_name = NULL) {
  # set column names in correct order
  # fold changes is group/reference
  col_names <- c(group_name = group, reference_name = reference)
  # sort for the column names
  col_names <- sort(col_names)

  # compute table
  de_table <- data %>%
    rename(stat_group = !!as.symbol(obs_name)) %>%
    # select the correct groups
    filter(.data$stat_group == reference |
             .data$stat_group == group) %>%
    # make long
    pivot_longer(cols = -c(.data$sample_id, .data$stat_group),
                 names_to = "lipid_name",
                 values_to = "value") %>%
    # get the data for each lipid
    nest(data = -.data$lipid_name) %>%
    # do the t-test
    mutate(ttest = map(.x = .data$data,
                       .f = ~ tidy(t.test(value ~ stat_group, data = .x)))) %>%
    # unfold all data
    unnest(cols = .data$ttest) %>%
    # remove unwanted columns
    select(-.data$data) %>%
    # rename some column names
    rename(names = .data$lipid_name,
           pvals = .data$p.value) %>%
    rename_with(~names(col_names), .cols = c(.data$estimate1, .data$estimate2)) %>%
    mutate(
      scores = -1,
      logfoldchanges = log2(.data$group_name / .data$reference_name),
      pvals_adj  = p.adjust(p = .data$pvals,
                            method = "fdr"),
      group = col_names["group_name"],
      comp_type = "grpVref",
      reference = col_names["reference_name"],
      test_type = "t-test",
      obs_name = obs_name,
      versus = paste(.data$group, "vs.", .data$reference)
    ) %>%
    select(.data$names, .data$scores, .data$logfoldchanges, .data$pvals,
           .data$pvals_adj, .data$group, .data$comp_type, .data$reference,
           .data$test_type, .data$obs_name, .data$versus)

  # return result
  return(de_table)
}


#' @title Perform a Mann-Whitney test on lipidyzer data
#'
#' @param data the data
#' @param group name of the group
#' @param reference name of the reference
#' @param obs_name column name which contains the group
#'
#' @return
#'
#' @author Rico Derks
#'
#' @importFrom dplyr filter mutate select rename rename_with summarise group_by pull
#' @importFrom tidyr pivot_longer pivot_wider nest unnest
#' @importFrom purrr map_dbl
#' @importFrom rlang !!
#' @importFrom broom tidy
#' @importFrom stats t.test p.adjust
#' @importFrom rlang .data
#'
#' @noRd
#'
perform_mwtest <- function(data = NULL,
                           group = NULL,
                           reference = NULL,
                           obs_name = NULL) {
  # set column names in correct order
  # fold changes is group/reference
  col_names <- c(group_name = group, reference_name = reference)
  # sort for the column names
  col_names <- sort(col_names)

  # compute table
  de_table <- data %>%
    rename(stat_group = !!as.symbol(obs_name)) %>%
    # select the correct groups
    filter(.data$stat_group == reference |
             .data$stat_group == group) %>%
    # make long
    pivot_longer(cols = -c(.data$sample_id, .data$stat_group),
                 names_to = "lipid_name",
                 values_to = "value") %>%
    # get the data for each lipid
    nest(data = -.data$lipid_name) %>%
    # do the Mann-Whitney-test
    mutate(pvals = map_dbl(.x = .data$data,
                           .f = ~ tidy(wilcox.test(value ~ stat_group, data = .x)) %>%
                             # remove unwanted columns
                             select(-.data$statistic, .data$method, .data$alternative) %>%
                             pull(.data$p.value))) %>%
    # unfold the data
    unnest(cols = .data$data) %>%
    # calculate the average of each group
    group_by(.data$lipid_name, .data$stat_group) %>%
    summarise(mean_value = mean(.data$value), .groups = "drop",
              pvals = .data$pvals[1]) %>%
    # make wide again
    pivot_wider(id_cols = -.data$stat_group,
                names_from = .data$stat_group,
                values_from = .data$mean_value) %>%
    # rename some columns
    rename(names = .data$lipid_name) %>%
    rename_with(~names(col_names), .cols = c(!!as.symbol(col_names[1]), !!as.symbol(col_names[2]))) %>%
    mutate(
      scores = -1,
      logfoldchanges = log2(.data$group_name / .data$reference_name),
      pvals_adj  = p.adjust(p = .data$pvals,
                            method = "fdr"),
      group = col_names["group_name"],
      comp_type = "grpVref",
      reference = col_names["reference_name"],
      test_type = "Mann-Whitney-test",
      obs_name = obs_name,
      versus = paste(.data$group, "vs.", .data$reference)
    ) %>%
    select(.data$names, .data$scores, .data$logfoldchanges, .data$pvals,
           .data$pvals_adj, .data$group, .data$comp_type, .data$reference,
           .data$test_type, .data$obs_name, .data$versus)

  # return the de table
  return(de_table)
}


#' @title Do some clean up of the lipidomics data
#'
#' @description Do some clean up of the lipidomics data. Remove all lipids which
#'     contain no data (all NA's).
#'
#' @param lipid_data list with all the data,
#'
#' @return list with the cleaned data.
#'
#' @author Rico Derks
#'
#' @noRd
#'
clean_up <- function(lipid_data = NULL) {
  ## remove columns with only NA's
  # get the id's
  remove_variable_idx <- which(apply(X = lipid_data$data_matrix,
                                     MARGIN = 2,
                                     FUN = function(x) {
                                       all(is.na(x))
                                     }))
  # remove from data matrix
  lipid_data$data_matrix <- lipid_data$data_matrix[, -remove_variable_idx]
  # remove from variable
  lipid_data$var_info <- lipid_data$var_info[-remove_variable_idx, ]

  # return cleaned data
  return(lipid_data)
}


#' @title Pre-processing of the lipid data
#'
#' @description Do several pre-processing steps on the lipid data.
#'
#' @param lipid_data list with all the lipid data
#' @param remove_twothird_lipids remove lipids which contain more than 2/3 NA's
#'     or zero's. Default is FALSE.
#'
#' @details If a remove_* parameter is set to FALSE removing of lipids is not
#'     done, but only saved in the variable info {var_info}. If set to TRUE the
#'     lipids are actually removed.
#'
#' @return list with pre-processed data
#'
#' @author Rico Derks
#'
#' @noRd
#'
preproc_lipid_data <- function(lipid_data = NULL,
                               remove_twothird_lipids = FALSE) {
  # Replace NA's by zero's: is this tricky because its sparse?
  lipid_data$data_matrix[which(is.na(lipid_data$data_matrix), arr.ind = TRUE)] <- 0

  ## remove lipids which contain 2/3 of zero's ##
  # get the index the lipids
  # check if 2/3 of a variable contains zero's
  excess_zero_conc <- colSums(lipid_data$data_matrix == 0) > 2/3 * dim(lipid_data$data_matrix)[1]
  # remove or record
  if(remove_twothird_lipids == TRUE) {
    # remove lipids
    lipid_data$data_matrix <- lipid_data$data_matrix[, -which(excess_zero_conc)]
    lipid_data$var_info <- lipid_data$var_info[-which(excess_zero_conc), ]
  } else {
    # record lipids
    lipid_data$var_info$excess_zero_conc <- excess_zero_conc
  }

  # return pre-processed data
  return(lipid_data)
}


#' @title Create the new database folder
#'
#' @description Check if the new database folder exists, if not create it.
#'
#' @param db_name name of the new database.
#' @param db_root path of the root of all databases.
#'
#' @return logical, indicating if the folder
#'
#' @author Rico Derks
#'
#' @noRd
#'
create_db_folder <- function(db_name = NULL,
                             db_root = NULL) {
  if(is.null(db_name) | is.null(db_root)) {
    stop("No database name or root supplied!")
  }

  # if directory doesn't exist create it
  if(!dir.exists(file.path(db_root, db_name))) {
    dir.create(path = file.path(db_root, db_name))
  }

  return(TRUE)
}


## code to prepare `yassene_example` dataset goes here
#######################################################################
#######################################################################
##
##  "WHAT IS THE CONTEXT?" lipidomics composition data
##  - data from yassenne / ricoo
##
#######################################################################
#######################################################################
#  formely `yassene_A` (conc and compos)
#
#
# Golem options (Inactive) ---------------------------
# Set options here
# options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode
# # Detach all loaded packages and clean your environment
# golem::detach_all_attached()
# # rm(list=ls(all.names = TRUE))
# # Document and reload your package
# golem::document_and_reload()


# ------------------------------------------
# 0. preamble/setup -------------------------
# ------------------------------------------
require(tidyverse)
require(reticulate)
reticulate::use_condaenv(required = TRUE, condaenv = 'sc39')
require(anndata)

require(glmnet)
require(readxl)
require(matrixStats)


# create the folder to contain the raw data
DB_NAME = "yassene_A_conc"
DB_DIR = file.path("data-raw",DB_NAME)
if (!dir.exists(DB_DIR)) {
  dir.create(DB_DIR)
}




# ------------------------------------------
# 1. documentation / provenance ------------
# ------------------------------------------


# TODO:  markdown file or html with some copy about the database
#  - lab, paper link/name
#  summarize results / data origin whatever


organism <- ""
lab <- ""
annotation_database <- ""



# ------------------------------------------
# 2. helper functions ----------------------
# ------------------------------------------



# ------------------------------------------
# 3. load data -----------------------------
# ------------------------------------------
RAW_DIR <- "ingest/Yassene_A"

# 2. load concentration data --------------------
csv_name <- "example_species_conc.csv"
file_name <- file.path(RAW_DIR,csv_name)
lipid_concentration <- read.csv(file=file_name, header=TRUE, as.is = TRUE)


conc_group <- lipid_concentration$Group
conc_group[conc_group==""] <- "CONTROL"  #i'm not sure this is true

meta_cols <- c("Name","Group")
obs_meta <- lipid_concentration[meta_cols]
obs_meta$Group <- conc_group

rownames(obs_meta) <- lipid_concentration$Name

conc_data <- lipid_concentration %>% dplyr::select(!matches(meta_cols))
rownames(conc_data) <- obs_meta$Name


conc_mat <- as.matrix( conc_data)
rownames(conc_mat) <- obs_meta$Name

class(conc_mat) <- "numeric"
#ZSCALE the matrix
zconc_raw <- scale(conc_mat) #including zeros

# zero out NA: is this tricky because its sparse?
conc_mat[which(is.na(conc_mat), arr.ind = TRUE)] <- 0

excess_zero_conc <- ( colSums(conc_mat==0) > 2/3*dim(conc_mat)[1] )
#poor_conc <- which(colSums(conc_mat==0) > 2/3*dim(conc_mat)[1])

# ZSCALE the matrix
zconc <- scale(conc_mat)


mu_conc <- colMeans(conc_mat)
var_conc <- colVars(conc_mat)

# obs_meta$mean <- rowMeans(conc_mat)
# obs_meta$var <- apply(conc_mat,1,var)

obs_meta$var_conc <- rowVars(conc_mat)
obs_meta$mean_conc <-rowMeans(conc_mat)


# experiment with the scaled (including)
regr_group <- obs_meta$Group
g <- unique(regr_group)
# g
#
test_conc <- zconc
regr_group <- as.numeric(factor(obs_meta$Group))

ind_rem_group <- which(regr_group == which(table(regr_group) < 3))

# remove groups with less than 3 observations..
if (length(ind_rem_group) > 0) {
  test_conc <- test_conc[-ind_rem_group, ]
  regr_group <- regr_group[-ind_rem_group]
}

lipids=colnames(test_conc)
var_ <- as.data.frame(lipids)

var_$mean_conc <- mu_conc
var_$var_conc <- var_conc

var_$excess_zero_conc <- excess_zero_conc

# choose the "significant" columns via lasso regression (glmnet)
set.seed(100)
cvfit <- cv.glmnet(test_conc, regr_group, nlambda = 100, alpha = .8, family = "multinomial", type.multinomial = "grouped")
coef <- coef(cvfit, s = "lambda.min")
tmp <- as.matrix(coef$"1")
tmp1 <- tmp[which(tmp != 0)]
coef_names <- rownames(tmp)[which(tmp != 0)][-1]
ind_coef <- which(colnames(test_conc) %in% coef_names)

uns <- list(comp_lasso_coef=coef)

var_$conc_sig_lasso_coef <- (colnames(test_conc) %in% coef_names)

rownames(var_) <- var_$lipids



# ------------------------------------------
# 4. pack into anndata                    --
# ------------------------------------------

X <- conc_mat
obs <- obs_meta

db_prefix = "core_data"
saveRDS(X, file = paste0(DB_DIR, "/", db_prefix, "_X.rds"))
saveRDS(obs, file = paste0(DB_DIR, "/", db_prefix, "_obs.rds"))
saveRDS(var_, file = paste0(DB_DIR, "/", db_prefix, "_var.rds"))
saveRDS(uns, file = paste0(DB_DIR, "/", db_prefix, "_uns.rds"))


ad <- AnnData(
  X = X,
  obs = obs,
  var = var_,
  uns = uns
)

ad



ad$write_h5ad(filename=file.path(DB_DIR,"core_data.h5ad"))



# ------------------------------------------
# 5. post processing                      --
# ------------------------------------------
# use scanpy to compute:
#     - differential measures (volcano)
#     - marginal stats
#     - dimension reductions (PCA,UMAP, tSNE)


#> AnnData object with n_obs × n_vars = 2 × 3
#>     obs: 'group'

require(anndata)
require(reticulate)
DB_NAME = "yassene_A_conc"
DB_DIR = file.path("data-raw",DB_NAME)
RAW_DIR <- "ingest/Yassene_A"

db_prefix = "core_data"
X = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_X.rds"))
obs = readRDS(  file = paste0(DB_DIR, "/", db_prefix, "_obs.rds"))
var_ = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_var.rds"))


ad <- read_h5ad(file.path(DB_DIR,"core_data.h5ad"))
ad


sc <- import("scanpy")


norm_adata <- sc$pp$normalize_total(ad,copy=TRUE)
#normpc_adata <- sc$pp$normalize_per_cell(ad,copy=TRUE)
z_adata <- sc$pp$scale(ad,copy=TRUE)


test_types <- c('wilcoxon','t-test_overestim_var')
comp_types <- c("grpVrest")


helper_function<-('data-raw/compute_de_table.R')
source(helper_function)

diff_exp <- compute_de_table(ad,comp_types, test_types, obs_names = 'Group')



# TODO:  add de_table into uns
#new_adata$uns$de_table <- diff_exp
ad$write_h5ad(filename=file.path(DB_DIR,"core_data_plus_de.h5ad"))

db_prefix = "de"
saveRDS(diff_exp, file = paste0(DB_DIR, "/", db_prefix, "_table.rds"))


# ------------------------------------------
# 6. dimension reduction - PCA / umap    --
# ------------------------------------------
require(anndata)
require(reticulate)

DB_NAME = "yassene_A_conc"
DB_DIR = file.path("data-raw",DB_NAME)
#RAW_DIR <- "ingest/Vilas_B"

ad <- read_h5ad(file.path(DB_DIR,"core_data_plus_de.h5ad"))




# ------------------------------------------
# 7 . create config and default files                   --
# ------------------------------------------
require(anndata)
require(reticulate)

DB_NAME = "yassene_A_conc"
DB_DIR = file.path("data-raw",DB_NAME)
#RAW_DIR <- "ingest/Vilas_B"

ad <- read_h5ad(file.path(DB_DIR,"core_data_plus_de.h5ad"))
db_prefix = "de"
diff_exp = readRDS( file = paste0(DB_DIR, "/", db_prefix, "_table.rds"))




###

# measures
#  This ordering is the "default"
measures <- list(obs = ad$obs_keys()[3:4],
                 var = ad$var_keys()[2:3])
# differentials  #if we care we need to explicitly state. defaults will be the order...
diffs <- list(diff_exp_groups =  levels(factor(diff_exp$group)),
              diff_exp_comp_type =  levels(factor(diff_exp$comp_type)),
              diff_exp_obs_name =  levels(factor(diff_exp$obs_name)),
              diff_exp_tests =  levels(factor(diff_exp$test_type)))

# Dimred
dimreds <- list(obsm = NA,
                varm = NA)

# what ad$obs do we want to make default values for...
# # should just pack according to UI?
default_factors <- c("Group")




helper_function<-('data-raw/create_config_table.R')

source(helper_function)

db_prefix = "omxr"
conf_and_def <- create_config_table(ad,
                                    measures,
                                    diffs,
                                    dimreds,
                                    default_factors,
                                    db_prefix= db_prefix,
                                    db_dir = DB_DIR)


yassene_A_conc_conf = readRDS( paste0(DB_DIR,"/",db_prefix,"_conf.rds") )
yassene_A_conc_def = readRDS( paste0(DB_DIR,"/",db_prefix,"_def.rds") )
yassene_A_conc_omics = readRDS( paste0(DB_DIR,"/",db_prefix,"_omics.rds") )
yassene_A_conc_meta = readRDS( paste0(DB_DIR,"/",db_prefix,"_meta.rds") )


usethis::use_data(yassene_A_conc_conf, yassene_A_conc_def, yassene_A_conc_omics, yassene_A_conc_meta, overwrite = TRUE)



















raw <- conc_mat
X <- zconc

layers <- list(comp = comp_mat)
layers$zcomp <- zcomp
layers$conc <- conc_mat
layers$zconc <- zconc

obsm <- NULL
varm <- NULL

## TODO:
##    1. calculate differential expression
##    2. calculate PCAs & UMAPS
##    3. incorporate lipid classes ... as feature_annotations (in var_)


###  here's some annotation
xlsx_name <- "Lipid_classification.xlsx"
file_name <- file.path(DATA_DIR,DB_NAME,xlsx_name)

lipid_class <- read_excel(file_name)

merge_data <- function(lipid_data, lipid_class = NULL, by = NULL) {
  # make the lipid data long and add some extra columns, same as in tidy_lipids
  lipid_data_long <- lipid_data %>%
    mutate(sample_type = factor(tolower(str_extract(string = .data$sample_name,
                                                    pattern = "([bB]lank|[qQ][cC]pool|[sS]ample)")))) %>%
    # join the meta data
    left_join(y = meta_data,
              by = c("sample_name" = by),
              suffix = c("", ".y"))


  return(lipid_data_long)
}

results_test <- do_stat_test(lipid_data = isolate(all_data$analysis_data),
                             group = input$test_select_group,
                             group1_name = input$test_group1,
                             group2_name = input$test_group2,
                             normalization = input$select_test_normalization,
                             transformation = input$select_test_transformation,
                             test = input$select_test)

volcano_plot(lipid_data = test_result(),
             pvalue_adjust = input$test_cor_pvalue,
             title = paste0(input$test_group1, " vs ", input$test_group2))


volcano_plot <- function(lipid_data, pvalue_adjust = FALSE, title = "") {
  # create y-axis title
  y_title <- ifelse(pvalue_adjust == FALSE,
                    "-log10(p value)",
                    "-log10(cor. p value)")

  # create the plot
  p <- lipid_data %>%
    mutate(show_p = case_when(
      pvalue_adjust == FALSE ~ .data$p_log10,
      pvalue_adjust == TRUE ~ .data$p_log10_adj
    )) %>%
    plot_ly(x = ~fc_log2,
            y = ~show_p,
            text = ~ShortLipidName,
            colors = rainbow(n = 100),
            customdata = lipid_data$ShortLipidName,
            source = "volcano_plot_click") %>%
    add_markers(color = ~LipidClass,
                size = 3) %>%
    layout(xaxis = list(zeroline = FALSE,
                        title = "log2(fold change)"),
           yaxis = list(title = y_title),
           shapes = list(vline(-1),
                         vline(1),
                         hline(-log10(0.05))),
           legend = list(orientation = "h"),
           title = list(text = title,
                        x = 0)) %>%
    event_register(event = "plotly_click")

  return(p)
}

vline <- function(x = 0, color = "blue") {
  list(
    type = "line",
    y0 = 0,
    y1 = 1,
    yref = "paper",
    x0 = x,
    x1 = x,
    line = list(color = color,
                width = 1,
                dash = "dash")
  )
}

hline <- function(y = 0, color = "blue") {
  list(
    type = "line",
    x0 = 0,
    x1 = 1,
    xref = "paper",
    y0 = y,
    y1 = y,
    line = list(color = color,
                width = 1,
                dash = "dash")
  )
}


# OBSERVABLES
#
observables <- list(obs = c("var_conc","mean_conc"),
                    var = c("var_conc", "var_comp","mean_conc","mean_comp"),
                    layers = c("comp", "zcomp","conc","zconc"),
                    raw = c("X"),
                    obsm = NA,
                    rawvar = NA)


# COMPARABLES
comparables <- list(varm = NA,
                    obsm = NA)
# Dimred
dimreds <- list(varm = NA,
                obsm = NA)



#######################################################################
#######################################################################
##
##  Create data tables / h5
##  1. config (XXX_conf.rds)  : structure of the database
##  2. defaults (XXX_def.rds) : the startup parameters
##  3. meta_data (XXX_meta.rds)  : data about the observations (obs)
##  4. omics  (XXX_omics.rds):  data about the  (var)
##  3. gene_expression (xxx_gexpr.h5):  count matrix of observ X var (X)
##
#######################################################################
#######################################################################

# X <- omicser::yassene_A_X
# obs <- omicser::yassene_A_obs
# var_ <- omicser::yassene_A_var
# obsm <- omicser::yyassene_A_obsm
# varm <- omicser::yassene_A_varm
# uns <-  omicser::yyassene_A_uns
# layers <-  omicser::yassene_A_layers


#TODO:  pack into a list or anndata structure for simplicity...

require("data.table")
## TODO:  check
##      - do we need default_1, 2 and multi?
##

helper_functions<-('data-raw/data_prep_functions.R')

source(helper_functions)

db_dir = "data-raw"
# ui_config <- make_ingest_file_primitives(X,obs,var_,obsm=NA,varm=NA,uns=NA,
#                                           gex.assay = NA, gex.slot = c("data", "scale.data", "counts"),
#                                           gene_mapping = FALSE, db_prefix = "Vilas_A", db_dir = db_dir,
#                                           default_omics1 = NA, default_omics2 = NA, default_multi = NA,
#                                           default_dimred = NA, chunk_size = 500, meta_to_include = NA, legend_cols = 4,
#                                           max_levels_ui = 50)

# make_ingest_files_primitive(X,obs,var_,obsm=obsm, varm=varm,
#                             uns=uns, layers = layers, gex.assay = NA, gex.slot = c("data", "scale.data", "counts"),
#                             gene_mapping = FALSE, db_prefix = db_prefix, db_dir = "data-raw",
#                             default_omics1 = NA, default_omics2 = NA, default_multi = NA,
#                             default_dimred = NA, chunk_size = 500, meta_to_include = NA, legend_cols = 4,
#                             max_levels_ui = 50)

db_prefix <- "yassene_A_"

make_ingest_file_primitives(X,obs,var_,obsm,varm,uns, layers,
                            observables, comparables, dimreds,
                            default_omic = NA, default_dimred = NA, meta_to_include = NA,
                            gex.assay = NA, gex.slot = c("data", "scale.data", "counts"),
                            gene_mapping = FALSE, db_prefix = db_prefix, db_dir = "data-raw",
                            chunk_size = 500,  legend_cols = 4,
                            max_levels_ui = 50)

#yassene_A_conf = readRDS(file.path(db_dir,"test1conf.rds"))
yassene_A_conf = readRDS( paste0(db_dir,"/",db_prefix,"conf.rds") )

# defaults:  list of meta1, meta2, omics1, omics2, omics (list of 10). dimred, grp1, grp2
yassene_A_def = readRDS( paste0(db_dir,"/",db_prefix,"def.rds") )

# list of vars )e/g/ 3000 genes with counts?
yassene_A_omics = readRDS( paste0(db_dir,"/",db_prefix,"omics.rds") )
# use this sorted one to resort everything before packing into anndata

yassene_A_meta = readRDS( paste0(db_dir,"/",db_prefix,"meta.rds") )
yassene_A_X = readRDS( paste0(db_dir,"/",db_prefix,"X.rds") )
yassene_A_obs = readRDS( paste0(db_dir,"/",db_prefix,"obs.rds") )
yassene_A_obsm = readRDS( paste0(db_dir,"/",db_prefix,"obsm.rds") )
yassene_A_var = readRDS( paste0(db_dir,"/",db_prefix,"var.rds") )
yassene_A_varm = readRDS( paste0(db_dir,"/",db_prefix,"varm.rds") )
yassene_A_uns = readRDS( paste0(db_dir,"/",db_prefix,"uns.rds") )
yassene_A_layers = readRDS( paste0(db_dir,"/",db_prefix,"layers.rds") )



# usethis::use_data(yassene_A_X,yassene_A_var,yassene_A_obs, overwrite = TRUE)
# usethis::use_data(yassene_A_obsm,yassene_A_varm,yassene_A_layers,yassene_A_uns, overwrite = TRUE)
usethis::use_data(yassene_A_conf, yassene_A_def, yassene_A_omics, yassene_A_meta, overwrite = TRUE)

#######################################################################
#######################################################################
##
##  Create some "differential tables" of type:
##  1. aggregated over "observations"
##  2. aggregated over "features"
##  3. aggretated over features and observations.
##
#######################################################################
#######################################################################

#######################################################################
#######################################################################
##
##  ANNDATA example
##
#######################################################################
#######################################################################

library(anndata)


# X <- omicser::yassene_A_X
# obs <- omicser::yassene_A_obs
# var_ <- omicser::yassene_A_var
# obsm <- omicser::yassene_A_obsm
# varm <- omicser::yassene_A_varm
# uns <- omicser::yassene_A_varm
# layers <- omicser::yassene_A_layers

X <- yassene_A_X
obs <- yassene_A_obs
var_ <- yassene_A_var
obsm <- yassene_A_obsm
varm <- yassene_A_varm
uns <-  yassene_A_uns
layers <-  yassene_A_layers


ad <- AnnData(
  X = X,
  obs = obs,
  var = var_,
  layers = layers,
  obsm = obsm,
  varm = varm,
  uns = uns
)

ad

# skip raw for now... all loaded into layers
# adraw <- AnnData(
#   X = raw,
#   var = var_
# )
#
# ad$raw <- adraw


#write_h5ad(anndata = ad, filename = file.path(db_dir,"data-raw/Vilas_A.h5ad"))
# anndata R wrapper is broken.. .invoke python
#
ad$write_h5ad(filename="data-raw/yassene_A.h5ad")
#> AnnData object with n_obs × n_vars = 2 × 3


## code to prepare `domenico_DIA`, `domenico_TMT``, and `domenico_transcriptomics`` dataset


##### DATASET PREPARATION ######
require(tidyverse)
require(data.table)

# Preprocess Transcriptomics Dataset -------------------------------------------------------------

DATA_DIR <- "/Users/ahenrie/Projects/NDCN_dev/dbrowse/ingest"
DATA_DIR <- "ingest"
DB_NAME <- "Domenico"


file_root <- file.path(DATA_DIR,DB_NAME)

DIA_csv <- "DIA/170823_aging_against_SC_merged_all_lib_2_all_candidates.csv"

file_path <- file.path(file_root,DIA_csv)



DIA_all_raw <- read.csv( file.path(file_root,DIA_csv), sep = ";", header = TRUE)
# subset transcriptomics dataset
DIA_proteins <- unique(DIA_all_raw %>% subset(select = c("UniProtIds", "Genes")))


load(file.path(file_root,"Transcriptomics/MuSC_proteome_transcriptome_merged.RData"))
#final.merge.with.union.RNAseq <-

#rename the CD34Hi and CD34Lo columns to include also the GSE-Number
final.merge.with.union.RNAseq <- rename(final.merge.with.union.RNAseq,
                                        "logFC_GSE155642.CD34Lo"="logFC_CD34Lo",
                                        "adj.P.Val_GSE155642.CD34Lo"="adj.P.Val_CD34Lo",
                                        "logFC_GSE155642.CD34Hi"="logFC_CD34Hi",
                                        "adj.P.Val_GSE155642.CD34Hi"="adj.P.Val_CD34Hi")
column_names <- c(names(DIA_proteins),
                  "logFC_GSE63860", "adj.P.Val_GSE63860",
                  "logFC_GSE47401", "adj.P.Val_GSE47401",
                  "logFC_GSE47177", "adj.P.Val_GSE47177",
                  "logFC_GSE155642.CD34Lo", "adj.P.Val_GSE155642.CD34Lo",
                  "logFC_GSE155642.CD34Hi", "adj.P.Val_GSE155642.CD34Hi")



# select only the rows that are also in Muscle-Stem-Cell Dataset
transcriptomics_data <- left_join(DIA_proteins, final.merge.with.union.RNAseq, by = c("UniProtIds")) %>% select(all_of(column_names))
# pivot the data into long-format
transcriptomics_long <- transcriptomics_data %>% pivot_longer(cols=3:12, names_to=c(".value", "label"), names_pattern="(.+)_(.+)")
# replace CD34Lo/Hi "." with "space", better for plotting name in heatmap later on (textwidth constraints)
transcriptomics_long$label <- gsub(".", " ", transcriptomics_long$label, fixed=TRUE)

#transcript_cols <- colnames(transcriptomics_data)[3:12] #used later to filter the dataset for csv-download
trdata_label <- unique(transcriptomics_long$label)      #used to create factor of label-column
transcriptomics_long$label <- factor(transcriptomics_long$label, levels = c("o / y", "g / y", "g / o", trdata_label))


#
# # save both dataframes into transcriptomics RData-File
# save(transcriptomics_long, file = "inputs/Transcriptomics/transcriptomics_data.RData")



# Preprocess DIA (Merge with transcriptomics) -----------------------------
# load("inputs/Transcriptomics/MuSC_proteome_transcriptome_merged.RData") #reload
# final.merge.with.union.RNAseq <- rename(final.merge.with.union.RNAseq,
#        "logFC_GSE155642.CD34Lo"="logFC_CD34Lo",
#        "adj.P.Val_GSE155642.CD34Lo"="adj.P.Val_CD34Lo",
#        "logFC_GSE155642.CD34Hi"="logFC_CD34Hi",
#        "adj.P.Val_GSE155642.CD34Hi"="adj.P.Val_CD34Hi")
# # reload DIA_all and reassert column_names
# DIA_all <- read.csv("inputs/DIA/170823_aging_against_SC_merged_all_lib_2_all_candidates.csv", sep = ";", header = TRUE)
DIA_all_subset <- DIA_all_raw %>% subset(select=c("Genes", "ProteinDescriptions", "UniProtIds", "AVG.Log2.Ratio", "Qvalue", "Comparison..group1.group2."))
column_names_subset <- c(colnames(DIA_all_subset),
                         "logFC_GSE63860", "adj.P.Val_GSE63860",
                         "logFC_GSE47401", "adj.P.Val_GSE47401",
                         "logFC_GSE47177", "adj.P.Val_GSE47177",
                         "logFC_GSE155642.CD34Lo", "adj.P.Val_GSE155642.CD34Lo",
                         "logFC_GSE155642.CD34Hi", "adj.P.Val_GSE155642.CD34Hi")
# Build merged "Trasncript-MuSc"-Dataset
# select only those rows that are also in Muscle-Stem-Cell Dataset
DIA_all_subset <- left_join(DIA_all_subset, final.merge.with.union.RNAseq, by = c("UniProtIds")) %>% select(all_of(column_names_subset))
#rename specified columns in order to change less parameters in several functions
DIA_all_subset <- rename(DIA_all_subset, logFC = AVG.Log2.Ratio, adj.P.Val = Qvalue, gene.name = Genes)
# Reorder Variable-levels for visulization
DIA_all_subset$Comparison..group1.group2. <- as.character(DIA_all_subset$Comparison..group1.group2.)
DIA_all_subset$Comparison..group1.group2. <- factor(DIA_all_subset$Comparison..group1.group2., levels = c("o / y", "g / y", "g / o"))

DIA_all <- DIA_all_subset

#save(DIA_all, file="inputs/DIA/Preprocessed_MuSc_Transcriptomics_merged.RData")
#write.table(DIA_all, "inputs/DIA/170823_aging_against_SC_merged_transcriptomics.csv", sep=";", dec = ".", row.names = FALSE, col.names = TRUE)



# Data Import TMT -------------------------------------------------------------
# Include Mass-Spec Datasets
list_names <- c("EDL_yg", "EDL_yo", "Gastr_yg", "Gastr_yo", "TA_yg", "TA_yo", "Sol_yg", "Sol_yo")
data_list <- setNames(replicate(length(list_names), data.frame()), list_names)
data_list[["EDL_yg"]] <- read.csv(file.path(file_root,"TMT/EDL_yg_new_annotated.csv"), sep = ",", header = TRUE, stringsAsFactors = FALSE)
data_list[["EDL_yo"]] <- read.csv(file.path(file_root,"TMT/EDL_yo_new_annotated.csv"), sep = ",", header = TRUE, stringsAsFactors = FALSE)
data_list[["Gastr_yg"]]  <- read.csv(file.path(file_root,"TMT/Gastr_yg_new_annotated.csv"), sep = ",", header = TRUE, stringsAsFactors = FALSE)
data_list[["Gastr_yo"]] <- read.csv(file.path(file_root,"TMT/Gastr_yo_new_annotated.csv"), sep = ",", header = TRUE, stringsAsFactors = FALSE)
data_list[["TA_yg"]] <- read.csv(file.path(file_root,"TMT/TA_yg_new_annotated.csv"), sep = ",", header = TRUE, stringsAsFactors = FALSE)
data_list[["TA_yo"]] <- read.csv(file.path(file_root,"TMT/TA_yo_new_annotated.csv"), sep = ",", header = TRUE, stringsAsFactors = FALSE)
data_list[["Sol_yg"]] <- read.csv(file.path(file_root,"TMT/Sol_yg_new_annotated.csv"), sep = ",", header = TRUE, stringsAsFactors = FALSE)
data_list[["Sol_yo"]] <- read.csv(file.path(file_root,"TMT/Sol_yo_new_annotated.csv"), sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Rename empty gene.name entries with the ID
for(i in list_names) {
  #data_list[[i]] <- data_list[[i]] %>% mutate(gene.name = as.character(gene.name))
  data_list[[i]]$gene.name[data_list[[i]]$gene.name == ""] <- NA
  data_list[[i]]$gene.name[is.na(data_list[[i]]$gene.name)] <- data_list[[i]]$ID[is.na(data_list[[i]]$gene.name)]
}

# Data Import TA 3rd Condition
TA_CTX <- read.csv(file.path(file_root,"TMT/TA_ctx_new_annotated.csv"), sep = ",", header = TRUE, stringsAsFactors = FALSE)
TA_CTX <- TA_CTX %>% dplyr::mutate(Group = factor("reg", levels = c("oy", "gy", "reg")))
# Rename empty gene.names with ID Column
TA_CTX$gene.name[TA_CTX$gene.name == ""] <- TA_CTX$ID[TA_CTX$gene.name == ""]

# Rename weird Proteins
rename_list <- list("P68433" = "H3c1", "P62806" = "H4c1", "P00761 SWISS-PROT:P00761" = "P00761")
for(i in list_names) {
  #print(paste0(i,":"))
  for(j in names(rename_list)) {
    data_list[[i]]$gene.name[data_list[[i]]$ID == j] <- rename_list[[j]]
    TA_CTX$gene.name[TA_CTX$ID == j] <- rename_list[[j]]
    #print(paste("renamed Protein", j,"to", rename_list[[j]]))
  }
}

# Save every Protein-Name from all TMT-datasets in one list
tmt_prot_list <- c()
for(i in list_names) {
  tmt_prot_list <- unique(c(tmt_prot_list, data_list[[i]]$gene.name))
}
tmt_prot_list <- unique(c(tmt_prot_list, TA_CTX$gene.name))



# Combine Sets ------------------------------------------------------------
# Function combines all conditions of each muscle-dataset
combineSets <- function(data_list, subset_string, naming) {
  for (i in 1:length(data_list)) {
    # Every second elem. in data_list, combine the two groups "gy" and "oy" into one dataset
    if(i %% 2 == 0) {
      data_list[[i]] <- data_list[[i]] %>% dplyr::mutate(Group = factor("oy", levels = c("oy", "gy", "reg")))
      data_list[[i-1]] <- data_list[[i-1]] %>% dplyr::mutate(Group = factor("gy", levels = c("oy", "gy", "reg")))
      a.name <- paste0(naming, sub(subset_string, '',names(data_list[i])), "_all")
      assign( a.name, rbind(data_list[[i-1]], data_list[[i]]), envir=parent.frame())
    }
  }
}
# Function call
combineSets(data_list, '_.*$', 'spec_')

# Add reg condition to TA_dataset
.spec_TA_all <- spec_TA_all
spec_TA_all <- rbind(spec_TA_all, TA_CTX)  #rbind adds TA_CTX asn new rows

# Select Muscle Dataset to work with --------------------------------------
muscles <- list(
  "Soleus" = spec_Sol_all,
  "Gastrocnemius" = spec_Gastr_all,
  "Tibialis Anterior" = spec_TA_all,
  "Extensor Digitorum Longus" = spec_EDL_all
)

# now pack into databases...



# Data Import DIA ---------------------------------------------------------
#DIA_all <- read.csv("inputs/DIA/170823_aging_against_SC_merged_all_lib_2_all_candidates.csv", sep = ";", header = TRUE)

save(DIA_all, file=file.path(file_root,"DIA/DIA_all.Rds"))

# save the data blobs into something easier to work with...

# transcriptomics_long
save(transcriptomics_long, file=file.path(file_root,"Transcriptomics/transcriptomics_long.Rds"))

# save TMT data
save(tmt_prot_list, file=file.path(file_root,"TMT/tmt_prot_list.Rds"))
save(muscles, file=file.path(file_root,"TMT/muscles.Rds"))



##########################################
##
##   set up databases for ingest
##
################################################


require(tidyverse)
DATA_DIR <- "/Users/ahenrie/Projects/NDCN_dev/dbrowse/ingest"
DATA_DIR <- "ingest"
DB_NAME <- "Domenico"


file_root <- file.path(DATA_DIR,DB_NAME)

# Data Import DIA ---------------------------------------------------------

#save(DIA_all, file=file.path(file_root,"DIA/DIA_all.Rds"))
load(file.path(file_root,"DIA/DIA_all.Rds"))
# pack into X, obs, vars, meta?
X <- NULL
obs <- DIA_all
vars <- NULL

# fixe some awkward names
obs <- obs %>% dplyr::rename(gene_name=gene.name)
obs <- obs %>% dplyr::rename(age_comparison=Comparison..group1.group2.)


domenico_DIA_X <- X
domenico_DIA_obs <- obs
domenico_DIA_vars <- vars

usethis::use_data(domenico_DIA_X, domenico_DIA_obs, domenico_DIA_vars,overwrite = TRUE)
usethis::use_data( domenico_DIA_obs,overwrite = TRUE)


# transcriptomics_long
#save(transcriptomics_long, file=file.path(file_root,"Transcriptomics/transcriptomics_long.Rds"))
load(file.path(file_root,"Transcriptomics/transcriptomics_long.Rds"))

# pack into X, obs, vars, meta?
X <- NULL
obs <- transcriptomics_long
vars <- NULL

domenico_A_X <- X
domenico_A_obs <- obs
domenico_A_vars <- vars

#usethis::use_data(domenico_TMT_X, domenico_TMT_obs, domenico_TMT_vars)
usethis::use_data(domenico_A_obs,overwrite = TRUE)






# save TMT data
#save(tmt_prot_list, file=file.path(file_root,"TMT/tmt_prot_list.Rds"))
#save(muscles, file=file.path(file_root,"TMT/muscles.Rds"))
load(file.path(file_root,"TMT/muscles.Rds"))

muscles1 <- muscles$Soleus
muscles$Soleus$muscle <- "Soleus"
muscles$Gastrocnemius$muscle <- "Gastrocnemius"
muscles$`Tibialis Anterior`$muscle <- "Tibialis Anterior"
muscles$`Extensor Digitorum Longus`$muscle <- "Extensor Digitorum Longus"

TMT_data_table <- dplyr::bind_rows(list(muscles$Soleus,muscles$Gastrocnemius,muscles$`Tibialis Anterior`,muscles$`Extensor Digitorum Longus`))
TMT_data_table <- dplyr::mutate(TMT_data_table,gene_x_muscle=paste0(gene.name,"_",muscle))


# pack into X, obs, vars, meta?
X <- NULL
obs <- TMT_data_table
vars <- NULL

# fixe some awkward names
obs <- obs %>% dplyr::rename(gene_name=gene.name)
obs$age_comparison <- obs$Group
#obs <- obs %>% dplyr::select(-c("GO","GO.ID")) # these are ugly...supress


domenico_TMT_X <- X
domenico_TMT_obs <- obs
domenico_TMT_vars <- vars

#usethis::use_data(domenico_TMT_X, domenico_TMT_obs, domenico_TMT_vars)
usethis::use_data(domenico_TMT_obs,overwrite = TRUE)



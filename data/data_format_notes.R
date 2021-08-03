# # AnnData vs SingleCellObject
# #The internals of an AnnData
# require(Matrix)
#
# # # DATA MATRIX - matrix
# # X = matrix(1:6, nrow = 2)
# # # annotation of OBSERVATIONS - data.frame (pandas) :: obs (obsm, obsp),
# # obs = data.frame(group = c("a", "b"), row.names = c("s1", "s2"))
# # obsm = list(
# #   ones = matrix(rep(1L, 10), nrow = 2),
# #   rand = matrix(rnorm(6), nrow = 2),
# #   zeros = matrix(rep(0L, 10), nrow = 2)
# # )
# #
# # # VARIABLES - :: var (varm, varp - named list())
# # var = data.frame(type = c(1L, 2L, 3L), row.names = c("var1", "var2", "var3"))
# # layers = list(
# #   spliced = matrix(4:9, nrow = 2),
# #   unspliced = matrix(8:13, nrow = 2)
# # )
# # varm = list(
# #   ones = matrix(rep(1L, 12), nrow = 3),
# #   rand = matrix(rnorm(6), nrow = 3),
# #   zeros =matrix(rep(0L, 12), nrow = 3)
# # )
# # # UNSTRUCTURED annotaion (dict)
# # uns = list(
# #   a = 1,
# #   b = data.frame(i = 1:3, j = 4:6, value = runif(3)),
# #   c = list(c.a = 3, c.b = 4)
# # )
#
#
# require(anndata)
#
# ad <- AnnData(
#   X=X,
#   obs = obs,
#   var = var,
#   layers = layers,
#   obsm = obsm,
#   varm = varm,
#   uns = uns
# )
# ad$var_names
# ad$n_vars
# ad$obs_names
# ad$n_obs0
#
# view <- ad[c("s1"),c("var1","var2")]
# view
#
# #
# # AnnData(
# #   X = NULL,
# #   obs = NULL,
# #   var = NULL,
# #   uns = NULL,
# #   obsm = NULL,
# #   varm = NULL,
# #   layers = NULL,
# #   raw = NULL,
# #   dtype = "float32",
# #   shape = NULL,
# #   filename = NULL,
# #   filemode = NULL,
# #   obsp = NULL,
# #   varp = NULL
# # )
# #
# # Raw(adata, X = NULL, var = NULL, varm = NULL)
#
#
# # SINGLE CELL
# library(SingleCellExperiment)
# counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
# sce <- SingleCellExperiment(counts)
# sce
#
# sce <- SingleCellExperiment(list(counts=counts))
# sce
#
# pretend.cell.labels <- sample(letters, ncol(counts), replace=TRUE)
# pretend.gene.lengths <- sample(10000, nrow(counts))
#
# sce <- SingleCellExperiment(list(counts=counts),
#                             colData=DataFrame(label=pretend.cell.labels),   # cell metadata
#                             rowData=DataFrame(length=pretend.gene.lengths), #feature metadata
#                             metadata=list(study="GSE111111")                # etc metadata
# )
# sce
#
# se <- SummarizedExperiment(list(counts=counts))
# as(se, "SingleCellExperiment")
#
# dim(assay(sce))
#
# colnames(colData(sce))
# colnames(rowData(sce))
#
# # if (!requireNamespace("BiocManager", quietly = TRUE))
# #   install.packages("BiocManager")
# #
# # BiocManager::install("scRNAseq")
# # library(scRNAseq)
# # sce <- ReprocessedAllenData("tophat_counts")
# # sce
#
# nrows <- 200
# ncols <- 6
# counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
# row_ranges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
#                      IRanges(floor(runif(200, 1e5, 1e6)), width=100),
#                      strand=sample(c("+", "-"), 200, TRUE),
#                      feature_id=sprintf("ID%03d", 1:200))
# col_data <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
#                      row.names=LETTERS[1:6])
#
# se1 <- SummarizedExperiment(assays=list(counts=counts),
#                      rowRanges=row_ranges, colData=col_data)
#
# # > se1
# # class: RangedSummarizedExperiment
# # dim: 200 6
# # metadata(0):                               # non row/col metadata - e.g. transcript lenth, gene symbol, genomic coordinates
# #   assays(1): counts                        # experiments... can be multiple e.g. raw count,  normalized data, etc
# # rownames: NULL
# # rowData names(1): feature_id               # feature metadat
# # colnames(6): A B ... E F
# # colData names(1): Treatment               # cell metadata ... e.g. batch of origin, treatment condition
#                                             # use for subsetting
# se2 <- SummarizedExperiment(assays=list(counts=counts), colData=col_data)
# se2
#
# x <- 1
# switch (x,
#          a = ,
#          b = 1,
#          c = 2,
#          stop("Unknown `x`", call. = FALSE)
# )
#
# switch(x, 1, 2, 3)
# x<- 3


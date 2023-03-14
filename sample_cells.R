suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(SingleCellExperiment)
  library(scran)
  library(scater)
  library(Matrix)
  library(stringr)
})

parser <- ArgumentParser()

parser$add_argument("--input_path",
  type = "character",
  help = "Input folder containing the MTX formatted data.."
)

args <- parser$parse_args()

create_sce <- function(path) {
  expression_matrix <- ReadMtx(mtx = file.path(path, "matrix.mtx"),
                               features = file.path(path, "bins.tsv"),
                               cells = file.path(path, "barcodes.tsv"))
  seurat_object <- CreateSeuratObject(counts = expression_matrix)
  return(as.SingleCellExperiment(seurat_object))
}

save_sce <- function(path, sce) {
	dir.create(path)
  writeMM(counts(sce), file.path(path, 'matrix.mtx'))
  write(colnames(sce), file.path(path, 'barcodes.tsv'))
  write.table(data.frame(rownames(sce))[c('rownames.sce.', 'rownames.sce.')],
              file=file.path(path, 'bins.tsv'),
              row.names=FALSE,
              col.names=FALSE,
              quote=FALSE,
              sep="\t")
}

format_output_path <- function(input_path, sampling) {
  to_parse <- str_split(input_path, '/')[[1]]
  return(file.path(
    '/home/gamazeps/data/matrices',
    to_parse[[6]], # dataset
    to_parse[[7]], # mark
    sampling,
    to_parse[[9]]  # feature engineeing
  ))
}

sce <- create_sce(args$input_path)
n_cells = ncol(sce)

# Here we reset the seed before each sampling in order to avoid each new sampling to be a superset
# of the previous one. See below the effect:
# > set.seed(0)
# > rnorm(2)
#   1.26295428488079 -0.326233360705649
#
# > set.seed(0)
# > rnorm(1)
#   1.26295428488079
#for (i in seq(500, n_cells, by=500)) {
#  # The seed is changed at each iteration because it would otherwise sample the same cells.
#  #
#  set.seed(i)
#  sampled_sce <- sce[, sample(ncol(sce), i, replace = FALSE)]
#  selection <- paste0('sampled_cell_n_', i, sep='')
#  path <- format_output_path(args$input_path, selection)
#  save_sce(path, sampled_sce)
#}

for (i in seq(20, 99, by=20)) {
  set.seed(i)
  sampled_sce <- sce[, sample(ncol(sce), as.integer(n_cells * i / 100), replace = FALSE)]
  selection <- paste0('sampled_cell_p_', i, sep='')
  path <- format_output_path(args$input_path, selection)
  save_sce(path, sampled_sce)
}

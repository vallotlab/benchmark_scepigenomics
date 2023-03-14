suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(Signac)
  library(SingleCellExperiment)
  library(scran)
  library(scater)
  library(Matrix)
#  library(ChromSCape)
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
cell_count <- colSums(counts(sce))

for (low in c(0, 10, 20, 30, 40, 50)) {
  bottom <- quantile(cell_count, low / 100, names=FALSE)
  top <- quantile(cell_count, (low + 50) / 100, names=FALSE)
  idx <- which((cell_count >= bottom) & (cell_count <= top))
  new_sce <- sce[, idx]
  selection <- paste0('filtered_cell_q_', low, '_', low + 50)
  path <- format_output_path(args$input_path, selection)
  save_sce(path, new_sce)
}

for (low in c(0, 10, 20, 30, 40, 50, 60)) {
  bottom <- quantile(cell_count, low / 100, names=FALSE)
  idx <- which(cell_count >= bottom)
  new_sce <- sce[, idx]
  selection <- paste0('filtered_cell_q_', low, '_', 100)
  path <- format_output_path(args$input_path, selection)
  save_sce(path, new_sce)
}

set.seed(0)
sampled_sce <- sce[, sample(ncol(sce), as.integer(ncol(sce) / 2), replace = FALSE)]
path <- format_output_path(args$input_path, 'filtered_cell_baseline')
save_sce(path, sampled_sce)

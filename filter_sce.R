suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(Signac)
  library(SingleCellExperiment)
  library(scran)
  library(scater)
  library(Matrix)
  library(ChromSCape)
})

parser <- ArgumentParser()

parser$add_argument("--input_path",
  type = "character",
  help = "Input folder containing the MTX formatted data.."
)
parser$add_argument("--output_path",
  type = "character",
  help = "Output basepath for the filtered matrices."
)
parser$add_argument("--percent_features",
  type = "double",
  help = "Number of regions to keep"
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
  writeMM(counts(sce), file.path(path, 'matrix.mtx'))
  write(colnames(sce), file.path(path, 'barcodes.tsv'))
  write.table(data.frame(rownames(sce))[c('rownames.sce.', 'rownames.sce.')],
              file=file.path(path, 'bins.tsv'),
              row.names=FALSE,
              col.names=FALSE,
              quote=FALSE,
              sep="\t")
}


sce <- create_sce(args$input_path)
k <- as.integer(nrow(sce) * args$percent_features)

#sce <- ChromSCape::find_top_features(sce, k)
dat <- FindVariableFeatures(as.Seurat(sce), nfeatures=k, selection.method='disp')
var.features <- dat@assays[[1]]@var.features
sce <- as.SingleCellExperiment(dat)
sce <- sce[match(var.features, rownames(sce)),]

save_sce(args$output_path, sce)

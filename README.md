# scHPTM benchmark

This repository contains the code for benchmarking the quality of the
representation of single-cell histone post translational modifications (scHPTM).

This is the companion to [Raimundo, F., Prompsy, P., Vert, J.P. and Vallot, C., 2022. Best practices for single-cell histone modification analysis. bioRxiv, pp.2022-09.](https://www.biorxiv.org/content/10.1101/2022.09.21.508811v1).

We evaluate the effect of the following factors:

-   Matrix construction.
-   Quality control (based on coverage).
-   Feature selection effect (using HVG or coverage).
-   Role of number of cells per experiment (by downsampling).
-   Role of coverage per cell (by downsampling cells based on coverage).
-   Role of normalization used (RPM or TF-IDF).
-   Role of dimension reduction algorithm (using 7 standard single-cell
    epigenetics pipelines).

The evaluation relies on having a robust co-assay (either transcriptomic or
surface proteins), and measuring how well the scHPTM representation conserves
the local geometry that can be observed in the reference co-assay.

## Building the original matrices

The original matrices are accessible in two formats:

- `.fragments.tsv`: for [Marsolier]() and [Zhang](). This is a raw format that contains the reads
                    per barcode.
- `matrix.mtx`: for [Zhu](), the data is already transformed into 5kbps bins.

### Folder structure

In this analysis we decided to follow a simple folder structure:

- All the data are stored under `$HOME/data`.
- The matrices are stored under `$HOME/data/matrices/$DATASET/$MARK/$CONDITION/$MATRIX_CONSTRUCTION`
- The embeddings are stored under `$HOME/data/matrices/$DATASET/$MARK/$CONDITION/$MATRIX_CONSTRUCTION`

Where:

- `$HOME`: is your home directory (`/home/username`)
- `$DATASET`: is one of `PairedTag_Zhu_2021`, `scChIP_Marsolier_2022`, `scCutTagPro_Zhang_2021`
- `$MARK`: is one of `H3K27ac`, `H3K27me3`, `H3K4me1`, `H3K4me3`, `H3K9me3`
- `$CONDITION`: is the in-silico experiment we did, unmodified data is the `raw` condition.
- `$MATRIX_CONDITION`: is either the binsize (e.g. `10k`) or the annotation used (e.g. `GeneTSS`).

### scChIP\_Marsolier\_2022

We downloaded the following H3K27me3 scChIP data from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164385):

- GSM5008855	`MM468_DMSO1_day60_H3K27me3 [scChIP-seq]`
- GSM5008856	`MM468_DMSO3_day77_H3K27me3 [scChIP-seq]`
- GSM5008857	`MM468_DMSO5_day131_H3K27me3 [scChIP-seq]`
- GSM5008858	`MM468_5FU1_day33_H3K27me3 [scChIP-seq]`
- GSM5008859	`MM468_5FU2_day67_H3K27me3 [scChIP-seq]`

The `DMSO` data are untreated MM468 cells, the `5FU` data are trested cells before resistance to
chemotherapy happens.

Using the web interface of ChromSCape we generated matrices for different binsizes (5kbps up to
1000kbps) in different analysis called `MarsollierK27.$BINSIZE`.

We also generated a matrix using the author annotation using the `_peaks.tsv.gz` files contained on
GEO and called it AuthorAnnotation.

These R objects are then converted into language agnostic 10x formated data using the
`Create Marsolier.ipynb` notebook, in that step we also filter regions known to have CNAs on this
cell line, as well as filter barcodes (more than 3k reads, less than 15k reads, and keep 10k most
covered regions).
You will need to replace `gamazeps` with your username in that notebook.

### PairedTag\_Zhu\_2021

The data is downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152020)

- `GSE152020_Paired-Tag_H3K27ac_DNA_filtered_matrix.tar.gz`
- `GSE152020_Paired-Tag_H3K27me3_DNA_filtered_matrix.tar.gz`
- `GSE152020_Paired-Tag_H3K4me1_DNA_filtered_matrix.tar.gz`
- `GSE152020_Paired-Tag_H3K4me3_DNA_filtered_matrix.tar.gz`
- `GSE152020_Paired-Tag_H3K9me3_DNA_filtered_matrix.tar.gz`

And each is moved to `$HOME/data/matrices/PairedTag_Zhu_2021/$MARK/raw/5k`.

Afterwards, the data is rebinned using `resize_bins.py`, here is an example command:

```
time python resize_bins.py \
--input_path=$HOME/data/matrices/PairedTag_Zhu_2021/H3K27me3/raw/5k \
--output_dir=$HOME/data/matrices/PairedTag_Zhu_2021/H3K27me3/raw/1000k \
--binsize=1000
```

The authors have already filtered the cells, so we keep them all.

### scCutTagPro\_Zhang\_2021

The data is downloaded from [zenodo]()

Matrices are built in the same fashion as for scChIP\_Marsolier\_2022

The cell filtering is done by keeping the same cells as the original authors (see the `` notebook).

## Building the preprocessing conditions

### Feature selection

Feature selection is done using the `filter_sce_features.R` script, here is the usage:

```
mkdir -p $HOME/data/matrices/PairedTag_Zhu_2021/H3K27me3/filtered_hvg_90/1000k
Rscript filter_sce_features.R \
--input_path=$HOME/data/matrices/PairedTag_Zhu_2021/H3K27me3/raw/1000k \
--output_path=$HOME/data/matrices/PairedTag_Zhu_2021/H3K27me3/filtered_hvg_90/1000k \
--mode=hvg
--fraction_features=0.9 \
```

- `fraction_features`: fraction of features to keep (0.9 for 90%)
- `mode`: either `hvg` for Seurat's HVG selection, or `coverage` for keeping the regions with the highest coverage.

### Cell selection by coverage

*Note*: you will need to replace `gamazeps` by your username.

Sample selection is done by using the `filter_cells_quality.R` script, here is the usage:

```
Rscript filter_sce_features.R \
--input_path=$HOME/data/matrices/PairedTag_Zhu_2021/H3K27me3/raw/1000k \
```

It will create the various `filtered_cell_q_x_y` (where the cells with coverage between x and y percentiles are kept).

### Downsamples cells

*Note*: you will need to replace `gamazeps` by your username.

Downsampling the cellls at random is done by using the `sample_cells.R` script, here is the usage:

```
Rscript sample_cells.R \
--input_path=$HOME/data/matrices/PairedTag_Zhu_2021/H3K27me3/raw/1000k \
```

It will create the various `sampled_cell_p_x`, where `x` is the percentage of cells kept.

## Running and evaluating the methods

### Building the embeddings

All the scripts assume that the data is stored in 10x mtx format.

```sh
time Rscript R_analysis.R \
--input_path=$HOME/data/matrices/scCutTagPro_Zhang_2021/H3K4me3/raw/50 \
--output_path=$HOME/data/embeddings/scCutTagPro_Zhang_2021/H3K4me3/raw/50 
time python3 nmf_process.py \
--input_path=$HOME/data/matrices/scCutTagPro_Zhang_2021/H3K4me3/raw/50 \
--output_path=$HOME/data/embeddings/scCutTagPro_Zhang_2021/H3K4me3/raw/50 
time python3 peakVI_process.py \
--input_path=$HOME/data/matrices/scCutTagPro_Zhang_2021/H3K4me3/raw/50 \
--output_path=$HOME/data/embeddings/scCutTagPro_Zhang_2021/H3K4me3/raw/50 
time python3 scale_process.py \
--input_path=$HOME/data/matrices/scCutTagPro_Zhang_2021/H3K4me3/raw/50 \
--output_path=$HOME/data/embeddings/scCutTagPro_Zhang_2021/H3K4me3/raw/50 
```

The embeddings will be stored in the `output_path`, each embeddin is stored in a `.csv` file named
after the method.

### Compute the scores

All the metrics will be stored in `$HOME/data/scores/DATASET/MARK/MATRIX_CONSTRUCTION`

The `compute_scores.py` script will attempt to compute the neighborhood score for all the different
matrix construction for a given mark, dataset and preporcessing.

- The `gt_path` is the `adt` extracted from the `.rds` object for scCUT&Tag, and the mode is `ADT`.
- The `gt_path` is the scRNA-seq matrix downloaded from GEO for the PairedTag, and the mode is `RNA`.

```sh
python3 compute_scores.py \
--output_path=$HOME/data/scores/scCutTagPro_Zhang_2021/H3K4me3/raw/  \
--embeddings_path=$HOME/data/embeddings/scCutTagPro_Zhang_2021/H3K4me3/raw/ \
--mode=ADT \
--gt_path=$HOME/data/input/scCutTagPro_Zhang_2021/H3K4me3.adt.csv 
```

The supervised metrics require an input 'ground truth' which is:
- `meta.tsv` from geo for Pairedtag
- The author annotation in the `.rds` files for scCut&TagPro
- Any matrix of the same mark for scChIP

```sh
time python3 compute_supervised_metrics.py '
--output_path=$HOME/data/scores/scCutTagPro_Zhang_2021/H3K4me1/raw/20k '
--embeddings_path=$HOME/data/embeddings/scCutTagPro_Zhang_2021/H3K4me1/raw/20k '
--gt_path=$HOME/data/input/scCutTagPro_Zhang_2021/H3K4me1')
```

### Analysis

- All the scores are stored under `data/scores`.
- All the runtime informations are stored in the `logs` folder.

The analysis colab is `scCUTTag_analysis.ipynb`.

# Original code

This code is based on a fork of the work originally done at Google, whose source
can be found [here](https://github.com/google-research/google-research/tree/master/schptm_benchmark).

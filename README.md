# scHPTM benchmark

This repository contains the code for benchmarking the quality of the
representation of single-cell histone post translational modifications (scHPTM).

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

# Running and evaluating the methods

All the packages assume that the data is stored in 10x mtx format.

Running the R methods on a matrix is done in the following way:

```sh
Rscript R_analysis.R  \
--input_path=/home/gamazeps/data/matrices/scCutTagPro_Zhang_2021/H3K4me3/raw/10k  \
--output_path=/home/gamazeps/data/embeddings/scCutTagPro_Zhang_2021/H3K4me3/raw/10k 
```

The embeddings will be stored in the `output_path`, each in a different csv file.


The `compute_scores.py` script will attempt to compute the neighborhood score for all the embeddings contained in its input folder

```sh
python3 compute_scores.py \
--output_path=/home/gamazeps/data/scores/scCutTagPro_Zhang_2021/H3K4me3/raw/10k  \
--embeddings_path=/home/gamazeps/data/embeddings/scCutTagPro_Zhang_2021/H3K4me3/raw/10k \
--mode=ADT \
--gt_path=/home/gamazeps/data/input/scCutTagPro_Zhang_2021/H3K4me3.adt.csv 
```


# Original code

This code is based on a fork of the work originally done at Google, whose source
can be found [here](https://github.com/google-research/google-research/tree/master/schptm_benchmark).

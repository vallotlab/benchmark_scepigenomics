import collections
import math
import os

from absl import logging
import anndata
import numpy as np
import scipy.io
import scipy.sparse
import tensorflow as tf


def grange(chrom, start, end):
  return f'{chrom}:{start}-{end}'


def bins_from_annotation(adata: anndata.AnnData,
                         annotation: os.PathLike) -> anndata.AnnData:
  """Generates a count matrix for the given annotation from a binned matrix.

  The binning size is infered from the first row of the matrix.

  Args:
    adata: Original binned matrix, all bins are assumed to be of the same size.
    annotation: Path to the annotation file.

  Returns:
    New count matrix
  """
  adata.var_names = adata.var_names.map(rename_coor)
  valid_bins = set(adata.var_names)
  start, end = adata.var_names[0].split(':')[1].split('-')
  binsize = int(end) - int(start)

  annot_index = []
  annot_rows = []
  with tf.io.gfile.GFile(annotation, 'r') as f:
    for line in f:
      splits = line.split(',')
      chrom, gene_start, gene_end = splits[0], int(splits[1]), int(splits[2])
      start = math.floor(gene_start / binsize) * binsize
      end = math.ceil(gene_end / binsize) * binsize
      acc = []
      for position in range(start, end, binsize):
        region = grange(chrom, position, position+binsize)
        if region in valid_bins:
          acc.append(region)
      if acc:
        annot_rows.append(adata[:, acc].X.sum(axis=1))
        annot_index.append(grange(chrom, gene_start, gene_end))

  new_adata = anndata.AnnData(scipy.sparse.csr_matrix(np.hstack(annot_rows)))
  new_adata.var_names = annot_index
  new_adata.obs = adata.obs
  return new_adata


def rename_coor(coor):
  if len(coor.split('_')) == 3:
    chrom, start, end = coor.split('_')
    return f'{chrom}:{start}-{end}'
  return coor


def merge_bins(adata, bin_size):
  orig_bins = collections.defaultdict(list)
  adata.var_names = adata.var_names.map(rename_coor)
  for coor in adata.var_names:
    chrom, start, end = coor.split(':')[0], int(coor.split(':')[1].split('-')[0]), int(coor.split(':')[1].split('-')[1])
    orig_bins[chrom].append((start, end))
  logging.info('Done with counting the bins')

  resized_bins_index = []
  resized_chrs = []
  resized_bins_counts = []
  for chrom, ranges in orig_bins.items():
    curr_bin = 0
    curr_acc = []
    for (start, end) in sorted(ranges):
      if start // bin_size == curr_bin:
        curr_acc.append(f'{chrom}:{start}-{end}')
      else:
        if curr_acc:
          # For the empty initialisation at the beginning of the chr.
          resized_bins_counts.append(adata[:, curr_acc].X.sum(axis=1))
          resized_bins_index.append(f'{chrom}:{curr_bin*bin_size}-{(curr_bin+1)*bin_size}')
        curr_acc = [f'{chrom}:{start}-{end}']
        curr_bin = start // bin_size
    resized_bins_counts.append(adata[:, curr_acc].X.sum(axis=1))
    resized_bins_index.append(f'{chrom}:{curr_bin*bin_size}-{(curr_bin+1)*bin_size}')
    resized_chrs.append(scipy.sparse.csr_matrix(np.hstack(resized_bins_counts)))
    resized_bins_counts = []
    logging.info(f'Done with {chrom}')

  new_adata = anndata.AnnData(scipy.sparse.csr_matrix(np.hstack([chrom.toarray() for chrom in resized_chrs])))
  new_adata.var_names = resized_bins_index
  new_adata.obs = adata.obs
  return new_adata

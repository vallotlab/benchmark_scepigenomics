"""Script for running PeakVI on a 10x formatted dataset."""

import os
from typing import Sequence

from absl import app
from absl import flags
import anndata
import pandas as pd
import scipy.io
import scipy.sparse
import scvi
import tensorflow as tf


FLAGS = flags.FLAGS
flags.DEFINE_string('input_path', None, 'Path to the 10x formatted folder.')
flags.DEFINE_string('output_path', None, 'Path to the output directory.')


def create_anndata(path: os.PathLike) -> anndata.AnnData:
  """Creates anndata object from raw data.

  Args:
    path: Path to the 10x formatted input files.

  Returns:
    anndata object for the experiment.
  """
  with tf.io.gfile.GFile(os.path.join(path, 'matrix.mtx'), mode='rb') as f:
    matrix = scipy.io.mmread(f)
  matrix = scipy.sparse.csr_matrix(matrix)
  adata = anndata.AnnData(matrix)
  adata = adata.transpose()
  with tf.io.gfile.GFile(os.path.join(path, 'barcodes.tsv'), mode='r') as f:
    barcodes = pd.read_csv(f, sep='\t', header=None)[0]
  adata.obs_names = barcodes
  with tf.io.gfile.GFile(os.path.join(path, 'bins.tsv'), mode='r') as f:
    bins = pd.read_csv(f, sep='\t', header=None)[0]
  adata.var_names = bins
  return adata


def main(argv: Sequence[str]) -> None:
  if len(argv) > 1:
    raise app.UsageError('Too many command-line arguments.')

  adata = create_anndata(FLAGS.input_path)
  scvi.model.PEAKVI.setup_anndata(adata)
  vae = scvi.model.PEAKVI(adata)
  vae.train()
  dr = pd.DataFrame(vae.get_latent_representation(), index=adata.obs_names)

  tf.io.gfile.makedirs(FLAGS.output_path)
  with tf.io.gfile.GFile(os.path.join(FLAGS.output_path, 'peakVI.csv'), 'w') as f:
    dr.to_csv(f)


if __name__ == '__main__':
  app.run(main)

import os
from typing import Sequence

from absl import app
from absl import flags
import anndata
import pandas as pd
import scipy.io
import scipy.sparse
import tensorflow as tf

import resize_bins_lib

FLAGS = flags.FLAGS
flags.DEFINE_string('input_path', None, 'Path to the 10x formatted folder.')
flags.DEFINE_string('output_dir', None, 'Path to the output directory.')
flags.DEFINE_integer('binsize', None, 'Number of bp per bin (in kbp).')
flags.DEFINE_enum('mode', 'bins', ['bins', 'annotation'],
                  'Number of bp per bin (in kbp)')
flags.DEFINE_string('annotation', None, 'Path to the annotation.')


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


def save_anndata(adata: anndata.AnnData, output_dir: os.PathLike,
                 input_path: os.PathLike):
  """Saves AnnData object in 10X format."""
  tf.io.gfile.makedirs(output_dir)
  with tf.io.gfile.GFile(os.path.join(output_dir, 'matrix.mtx'), mode='w') as f:
    scipy.io.mmwrite(f, adata.X.transpose())
  new_bins = pd.DataFrame(adata.var_names, columns=['var_names'])
  with tf.io.gfile.GFile(os.path.join(output_dir, 'bins.tsv'), mode='w') as f:
    new_bins.to_csv(
        f,
        sep='\t',
        index=False,
        header=False,
        columns=['var_names', 'var_names'])
  tf.io.gfile.copy(
      os.path.join(input_path, 'barcodes.tsv'),
      os.path.join(output_dir, 'barcodes.tsv'),
      overwrite=True)


def main(argv: Sequence[str]) -> None:
  del argv

  adata = create_anndata(FLAGS.input_path)
  if FLAGS.mode == 'bins':
    adata = resize_bins_lib.merge_bins(adata, FLAGS.binsize*(10**3))
  elif FLAGS.mode == 'annotation':
    adata = resize_bins_lib.bins_from_annotation(adata, FLAGS.annotation)

  save_anndata(adata, FLAGS.output_dir, FLAGS.input_path)


if __name__ == '__main__':
  flags.mark_flag_as_required('input_path')
  flags.mark_flag_as_required('output_dir')
  flags.mark_flag_as_required('binsize')
  app.run(main)

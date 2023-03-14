"""TODO(fraimundo): DO NOT SUBMIT without one-line documentation for compute_scores.

TODO(fraimundo): DO NOT SUBMIT without a detailed description of compute_scores.
"""

import os
from typing import Sequence

from absl import app
from absl import flags
from absl import logging
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse
from sklearn.neighbors import kneighbors_graph
import tensorflow as tf

FLAGS = flags.FLAGS
flags.DEFINE_string('gt_path', None, 'Path to the scRNA-seq ground truth.')
flags.DEFINE_string('embeddings_path', None, 'Path to the embeddings folder.')
flags.DEFINE_string('output_path', None, 'Directory where to write the scores.')
flags.DEFINE_string('mode', None, '')
flags.DEFINE_boolean('all', False, 'Whether to run all feature engineerings')

FILES = {
    'LSI.csv': 'Chromscape_LSI',
    'Chromscape_white_LSI.csv': 'Chromscape_white_LSI',
    'Chromscape_white_LSI_50.csv': 'Chromscape_LSI_50',
    'pca.csv': 'Chromscape_PCA',
    'Signac.csv': 'Signac',
    'Signac_10.csv': 'Signac_10',
    'SnapATAC.csv': 'SnapATAC',
    'cisTopic.csv': 'cisTopic',
    'peakVI.csv': 'PeakVI',
    'SCALE.csv': 'SCALE',
    'nmf.csv': 'NMF',
    'tfidf_nmf.csv': 'TFIDF-NMF',
}


_FRACTIONS = [.001, .003, .005, .01, 0.03, .05, .1]


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
  with tf.io.gfile.GFile(os.path.join(path, 'genes.tsv'), mode='r') as f:
    bins = pd.read_csv(f, sep='\t', header=None)[0]
  adata.var_names = bins
  with tf.io.gfile.GFile(os.path.join(path, 'meta.tsv'), mode='r') as f:
    meta = pd.read_csv(f, sep='\t', index_col='Cell_ID')
  adata.obs = meta
  return adata


def process_rna(rna):
  sc.pp.calculate_qc_metrics(rna, percent_top=None, log1p=False, inplace=True)
  sc.pp.normalize_total(rna, target_sum=1e4)
  sc.pp.log1p(rna)
  sc.pp.highly_variable_genes(rna, min_mean=0.0125, max_mean=3, min_disp=0.5)
  rna = rna[:, rna.var.highly_variable]
  sc.pp.regress_out(rna, ['total_counts'])
  sc.pp.scale(rna, max_value=10)
  sc.tl.pca(rna, svd_solver='arpack')
  return rna


def create_adt_h5a(path):
  with tf.io.gfile.GFile(path, mode='r') as f:
    adt = pd.read_csv(f, index_col=0).transpose()
  adata = anndata.AnnData(adt)
  adata.obs_names = adata.obs_names.map(lambda x: '-'.join(x.split('.')))
  return adata


def process_adt(adata):
  sc.tl.pca(adata, svd_solver='arpack')
  return adata


def compute_mutual_mixing(nng1, latent2, k):
  nng2 = kneighbors_graph(latent2, n_neighbors=k).todense().astype('bool')
  both = np.logical_and(nng1, nng2).sum()
  return both / (latent2.shape[0] * k)


def evaluate_methods(gt_path: os.PathLike,
                     embeddings_path: os.PathLike,
                     source: str,
                     mark: str,
                     feature_selection: str,
                     binsizes: Sequence[str],
                     out_paths: Sequence[str],
                     fractions: Sequence[float],
                     mode: str = 'RNA'):
  """Computes neighborhood score.

  Args:
    gt_path: Path to the 10X RNA formatted data and annotation.
    embeddings_path: Path to the csv embeddings.
    source: Dataset of origin.
    mark: Histone mark target.
    feature_selection: Feature selection used before DR.
    binsize: Feature engineering methods used to generate the matrix.
    out_paths: Paths to the scores.
    fractions: Percentage of the cells to use for kNN graph.

  Returns:
    Neighborhood score for the methods.

  """
  idx = None
  embeddings = dict()
  adata = None



  nng_rna = dict()

  for binsize, out_path in zip(binsizes, out_paths):
    for file, name in FILES.items():
      if tf.io.gfile.exists(os.path.join(embeddings_path, binsize, file)):
        with tf.io.gfile.GFile(
            os.path.join(embeddings_path, binsize, file), mode='r') as f:
          embeddings[name] = pd.read_csv(f, index_col=0)

    if binsize[-1] == 'k' and ((binsize != 'MacsPseudoBulk') and (binsize != 'SicerPseudoBulk')):
      binsize = int(binsize[:-1])

    if not embeddings:
      logging.info('There are no embeddings available.')
      continue

    if not adata:
      # We take the ids of the cells in the sample from an embeddings.
      # The fact that all embeddings are on the same cells is enforced when they are
      # added in the obsm field.
      for _, emb in embeddings.items():
        idx = emb.index
        break

      if mode == 'RNA':
        adata = create_anndata(gt_path)
        adata = adata[idx, :]
        logging.info('Starting to process the RNA matrix')
        adata = process_rna(adata)
        logging.info('Built the RNA matrix')
      elif mode == 'ADT':
        adata = create_adt_h5a(gt_path)
        adata = adata[idx, :]
        logging.info('Starting to process the ADT matrix')
        adata = process_adt(adata)
        logging.info('Built the ADT matrix')

    res = []
    for fraction in fractions:
      k = round(adata.n_obs * fraction)
      if k == 0:
        continue
      if k not in nng_rna:
        nng_rna[k] = kneighbors_graph(
            adata.obsm['X_pca'], n_neighbors=k).todense().astype('bool')
      for name, emb in embeddings.items():
        adata.obsm[name] = emb  # Enforces all embeddings being on the same cells.
        score = compute_mutual_mixing(nng_rna[k], adata.obsm[name], k)
        res.append({
            'Method': name,
            'K': k,
            'fraction': fraction,
            'frac. shared neighbors': score,
            'Mark': mark,
            'Dataset': source,
            'binsize': binsize,
            'feature_selection': feature_selection,
        })
        logging.info('Done with %s', name)

    if len(res) == 0:
      continue

    scores = pd.DataFrame(res)
    tf.io.gfile.makedirs(out_path)
    with tf.io.gfile.GFile(
        os.path.join(out_path, 'scores.csv'), mode='w') as f:
      scores.to_csv(f)
    logging.info('Done with %s', binsize)


def main(argv: Sequence[str]) -> None:
  if len(argv) > 1:
    raise app.UsageError('Too many command-line arguments.')

  to_parse = FLAGS.embeddings_path.split('/')
  if FLAGS.all:
    binsizes = os.listdir(FLAGS.embeddings_path)
    out_paths = [os.path.join(FLAGS.output_path, bs) for bs in binsizes]
    embeddings_path = FLAGS.embeddings_path
    feature_selection = to_parse[-1]
    mark = to_parse[-2]
    source = to_parse[-3]
    logging.info(f'{feature_selection}, {mark}, {source}')
    logging.info(binsizes)
    logging.info(out_paths)
  else:
    binsizes = [to_parse[-1]]
    out_paths = [FLAGS.output_path]
    embeddings_path = os.path.join(to_parse[:-1])
    feature_selection = to_parse[-2]
    mark = to_parse[-3]
    source = to_parse[-4]

  scores = evaluate_methods(
      gt_path=FLAGS.gt_path,
      embeddings_path=embeddings_path,
      source=source,
      mark=mark,
      feature_selection=feature_selection,
      binsizes=binsizes,
      out_paths=out_paths,
      fractions=_FRACTIONS,
      mode=FLAGS.mode,
  )


if __name__ == '__main__':
  app.run(main)

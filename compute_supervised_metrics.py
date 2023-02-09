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
import sklearn.cluster
import sklearn.metrics
from sklearn.neighbors import kneighbors_graph
import tensorflow as tf

FLAGS = flags.FLAGS
flags.DEFINE_string('gt_path', None, 'Path to the scRNA-seq ground truth.')
flags.DEFINE_string('embeddings_path', None, 'Path to the embeddings folder.')
flags.DEFINE_string('output_path', None, 'Directory where to write the scores.')

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


def create_anndata(path: os.PathLike, source: str) -> anndata.AnnData:
  """Creates anndata object from raw data.

  Args:
    path: Path to the 10x formatted input files.

  Returns:
    anndata object for the experiment.
  """
  if source == 'scCutTagPro_Zhang_2021':
    with tf.io.gfile.GFile(os.path.join(path, 'adt.csv'), mode='r') as f:
      adt = pd.read_csv(f, index_col=0).transpose()
    adata = anndata.AnnData(adt)
    with tf.io.gfile.GFile(os.path.join(path, 'l1.csv'), mode='r') as f:
      labels = pd.read_csv(f, sep=',', index_col=0)['x']
    labels.index = adata.obs_names
    adata.obs['Annotation'] = labels
    adata.obs_names = adata.obs_names.map(lambda x: '-'.join(x.split('.')))
    return adata


  with tf.io.gfile.GFile(os.path.join(path, 'matrix.mtx'), mode='rb') as f:
    matrix = scipy.io.mmread(f)
  matrix = scipy.sparse.csr_matrix(matrix)
  adata = anndata.AnnData(matrix)
  adata = adata.transpose()
  with tf.io.gfile.GFile(os.path.join(path, 'barcodes.tsv'), mode='r') as f:
    barcodes = pd.read_csv(f, sep='\t', header=None)[0]
  adata.obs_names = barcodes

  if source == 'PairedTag_Zhu_2021':
    features_fp = 'genes.tsv'
  if source == 'scChIP_Grosselin_2019':
    features_fp = 'bins.tsv'
  with tf.io.gfile.GFile(os.path.join(path, features_fp), mode='r') as f:
    bins = pd.read_csv(f, sep='\t', header=None)[0]
  adata.var_names = bins

  if source == 'PairedTag_Zhu_2021':
    with tf.io.gfile.GFile(os.path.join(path, 'meta.tsv'), mode='r') as f:
      meta = pd.read_csv(f, sep='\t', index_col='Cell_ID')
    adata.obs = meta
  if source == 'scChIP_Grosselin_2019':
    with tf.io.gfile.GFile(os.path.join(path, 'annot.txt'), mode='r') as f:
      labels = pd.read_csv(f, sep='\t', header=None)[0]
    labels.index = adata.obs_names
    adata.obs['Annotation'] = labels

  return adata


def evaluate_methods(gt_path: os.PathLike,
                     embeddings_path: os.PathLike,
                     source: str,
                     mark: str,
                     feature_selection: str,
                     binsize: str):
  """Computes supervised metrics.

  Args:
    gt_path: Path to the 10X RNA formatted data and annotation.
    embeddings_path: Path to the csv embeddings.
    source: Dataset of origin.
    mark: Histone mark target.
    feature_selection: Feature selection used before DR.
    binsize: Feature engineering method used to generate the matrix.

  Returns:
    Returns pd.DataFrame containing the supervised metrics.

  """
  idx = None
  embeddings = dict()

  for file, name in FILES.items():
    if tf.io.gfile.exists(os.path.join(embeddings_path, file)):
      with tf.io.gfile.GFile(
          os.path.join(embeddings_path, file), mode='r') as f:
        embeddings[name] = pd.read_csv(f, index_col=0)

  if not embeddings:
    logging.info('There are no embeddings available.')
    return None

  # We take the ids of the cells in the sample from an embeddings.
  # The fact that all embeddings are on the same cells is enforced when they are
  # added in the obsm field.
  for _, emb in embeddings.items():
    idx = emb.index
    break

  # get labels
  adata = create_anndata(gt_path, source)
  adata = adata[idx, :]
  logging.info('Built the gt matrix')

  if binsize[-1] == 'k' and binsize != 'PseudoBulk':
    binsize = int(binsize[:-1])

  res = []
  for name, emb in embeddings.items():
    adata.obsm[name] = emb  # Enforces all embeddings being on the same cells.

    silhouette = sklearn.metrics.silhouette_score(adata.obsm[name],
                                                  adata.obs['Annotation'])

    kmeans = sklearn.cluster.KMeans(
        n_clusters=len(adata.obs['Annotation'].unique()), random_state=0).fit(adata.obsm[name])
    adata.obs['predicted_clusters'] = kmeans.labels_

    kmeans_ari = sklearn.metrics.adjusted_rand_score(adata.obs['Annotation'],
                                              adata.obs['predicted_clusters'])
    kmeans_ami = sklearn.metrics.adjusted_mutual_info_score(
        adata.obs['Annotation'], adata.obs['predicted_clusters'])

    ward = sklearn.cluster.AgglomerativeClustering(
        n_clusters=len(adata.obs['Annotation'].unique())).fit(adata.obsm[name])
    adata.obs['predicted_clusters'] = ward.labels_

    ward_ari = sklearn.metrics.adjusted_rand_score(adata.obs['Annotation'],
                                              adata.obs['predicted_clusters'])
    ward_ami = sklearn.metrics.adjusted_mutual_info_score(
        adata.obs['Annotation'], adata.obs['predicted_clusters'])

    res.append({
        'Method': name,
        'Mark': mark,
        'Dataset': source,
        'binsize': binsize,
        'feature_selection': feature_selection,
        'kmeans_ami': kmeans_ami,
        'kmeans_ari': kmeans_ari,
        'ward_ami': ward_ami,
        'ward_ari': ward_ari,
        'silhouette': silhouette,
    })
    logging.info('Done with %s', name)
  return pd.DataFrame(res)


def main(argv: Sequence[str]) -> None:
  if len(argv) > 1:
    raise app.UsageError('Too many command-line arguments.')

  to_parse = FLAGS.embeddings_path.split('/')
  binsize = to_parse[-1]
  feature_selection = to_parse[-2]
  mark = to_parse[-3]
  source = to_parse[-4]

  scores = evaluate_methods(
      gt_path=FLAGS.gt_path,
      embeddings_path=FLAGS.embeddings_path,
      source=source,
      mark=mark,
      feature_selection=feature_selection,
      binsize=binsize,
  )

  if scores is None:
    return

  tf.io.gfile.makedirs(FLAGS.output_path)
  with tf.io.gfile.GFile(
      os.path.join(FLAGS.output_path, 'supervised_scores.csv'), mode='w') as f:
    scores.to_csv(f)


if __name__ == '__main__':
  app.run(main)

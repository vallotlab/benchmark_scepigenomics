import anndata
import numpy as np
import scipy.sparse
import tensorflow as tf

from google3.learning.brain.research.combini.single_cell.sc_cuttag import resize_bins_lib

def generate_dummy_data(binsize=1000):
  test_X = scipy.sparse.csr_matrix(np.eye(10))
  test_bins = [f'chr1:{i*binsize}-{(i+1)*binsize}' for i in range(0, 10)]
  #test_bins += [f'chr2:{i*binsize}-{(i+1)*binsize}' for i in range(0, 10)]
  test_adata = anndata.AnnData(test_X)
  test_adata.var_names = test_bins
  return test_adata

class MergeBinsTest(tf.test.TestCase):

  def test_proper_sum(self):
    test_adata = generate_dummy_data(1000)
    resized_adata = resize_bins_lib.merge_bins(test_adata, 2000)

    self.assertAllEqual(
        np.sum(test_adata.X, axis=1),
        np.sum(resized_adata.X, axis=1))


if __name__ == '__main__':
  tf.test.main()


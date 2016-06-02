from scipy.stats import norm
from dtcut import DTCUT_test;
import numpy as np;

class Test(DTCUT_test):

  recompute  = True;
  norm_dist  = None
  counter_i  = 0;
  clt_thresh = 10;

  def prepare_data(self):
    print "in prepare_data"
    self.norm_dist = (self.T.L.mean(), self.T.L.std());
    print self.norm_dist
  #edef

  #############################################################################

  def test(self, i_node):
    """ Do a CLT approximation to a permutation test
    """

    self.counter_i = self.counter_i + 1
    print "test %d" % self.counter_i;

    leaves  = self.T.get_tree_node_leaves(i_node);
    nleaves = len(leaves);

    cluster_score = np.sum(self.T.L[np.array(list(leaves))]);

    clt_mean = nleaves*self.norm_dist[0];
    clt_std  = nleaves*self.norm_dist[1];

    p_value = 1 - norm.cdf( (cluster_score - clt_mean) / clt_std);
    return p_value;

  #edef

  #############################################################################

#eclass

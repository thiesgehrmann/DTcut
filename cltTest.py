from scipy.stats import norm
from dtcut import DTCUT_test;
import numpy as np;

class Test(DTCUT_test):

  recompute    = True;
  norm_dist    = None
  clt_thresh   = 10;
  permutations = []

  #############################################################################

  def prepare_data(self):
    self.norm_dist    = (self.T.L.mean(), self.T.L.std());
    self.permutations = [ None for x in self.T.T ];
  #edef

  #############################################################################

  def test(self, i_node, **kwargs):

    if len(self.T.get_tree_node_leaves(i_node)) >= 10:
      return self.T.get_tree_node_p_value(i_node) if self.T.get_tree_node_visited(i_node) else self.cltTest(i_node);
    else:
      return self.permutationTest(i_node, kwargs['factor']);
    #fi

  #edef

  #############################################################################

  def permutationTest(self, i_node, factor):

    # First, determine how many permutations we need to get significant results
    alpha = self.p_thresh;
    nperms_to_sig = max(1000, ((factor+1) / alpha) * 1.1);

    # Check how many we have already done
    if self.permutations[i_node] == None:
      nperms, perms = self.permutations[i_node];
    #fi

    # And perform the necessary number of permutations to get it
    cluster_size = len(self.T.get_tree_node_leaves(i_node));
    new_perms    = np.concatenate([perms, np.array([self.performPermutation(cluster_size) for x in xrange(nperms_to_sig-nperms)])]);

    p_value = np.sum(new_perms > self.score(i_node)) / float(nperms_to_sig)

    self.permutations[i_node] = (nperms_to_sig, sorted(new_perms, descend=True)[0:int(self.p_thresh*nperms_to_sig)]);

    return p_value

  #edef

  #############################################################################

  def performPermutation(self, n):
    return np.sum(np.random.sample(self.T.L, n));
  #edef

  #############################################################################

  def score(self, i_node):
    leaves  = self.T.get_tree_node_leaves(i_node);
    nleaves = len(leaves);

    return np.sum(self.T.L[np.array(list(leaves))]);
  #edef

  #############################################################################

  def cltTest(self, i_node):
    """ Do a CLT approximation to a permutation test
    """

    nleaves       = len(self.T.get_tree_node_leaves(i_node));
    cluster_score = self.score(i_node);

    clt_mean = nleaves*self.norm_dist[0];
    clt_std  = nleaves*self.norm_dist[1];

    p_value = 1 - norm.cdf( (cluster_score - clt_mean) / clt_std);
    return p_value;

  #edef

  #############################################################################

#eclass

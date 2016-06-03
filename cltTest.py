from scipy.stats import norm
from dtcut import DTCUT_test;
import numpy as np;

class Test(DTCUT_test):

  recompute    = False;
  norm_dist    = None
  clt_thresh   = 10;
  permutations = []
  tests_done   = 0;

  #############################################################################

  def prepare_data(self):
    self.norm_dist    = (self.T.L.mean(), self.T.L.std());
    self.permutations = [ None for x in self.T.T ];
  #edef

  #############################################################################

  def test(self, i_node, **kwargs):

    print "entering test"
    #if len(self.T.get_tree_node_leaves(i_node)) >= 10:
    pval = self.T.get_tree_node_p_value(i_node) if self.T.get_tree_node_visited(i_node) else self.cltTest(i_node);
    #else:
      #pval = self.permutationTest(i_node, kwargs['factor']);
    #fi
    print "leaving test"
    self.tests_done = self.tests_done + 1 - int(self.T.get_tree_node_visited(i_node));
    return pval
  #edef

  #############################################################################

  def permutationTest(self, i_node, factor):

    print "entering permutationTest"

    # First, determine how many permutations we need to get significant results
    alpha = self.p_thresh;
    nperms_to_sig = max(1000, int(((factor+1) / alpha) * 1.1));

    # Check how many we have already done
    if self.permutations[i_node] == None:
      nperms, perms = (0, np.array([]));
    else:
      nperms, perms = self.permutations[i_node];
    #fi

    print "performing %d permutations" % (nperms_to_sig-nperms);

    # And perform the necessary number of permutations to get it
    cluster_size = len(self.T.get_tree_node_leaves(i_node));
    new_perms    = np.concatenate([perms, np.array([np.sum(np.random.choice(self.T.L, cluster_size, replace=False)) for x in xrange(nperms_to_sig-nperms)])]);

    p_value = np.sum(new_perms > self.score(i_node)) / float(nperms_to_sig)

    self.permutations[i_node] = (nperms_to_sig, sorted(new_perms, reverse=True)[0:int(self.p_thresh*nperms_to_sig)]);

    print "leaving permutationTest"
    return p_value

  #edef

  #############################################################################

  def score(self, i_node):
    leaves  = self.T.get_tree_node_leaves(i_node);
    return np.sum(self.T.L[np.array(list(leaves))]);
  #edef

  #############################################################################

  def cltTest(self, i_node):
    """ Do a CLT approximation to a permutation test
    """

    nleaves       = len(self.T.get_tree_node_leaves(i_node));
    cluster_score = self.score(i_node);

    print "cluster_score (%d, %d): %f" % (i_node, nleaves, cluster_score)

    clt_mean = nleaves*self.norm_dist[0];
    clt_std  = nleaves*self.norm_dist[1];

    p_value = 1 - norm.cdf( (cluster_score - clt_mean) / clt_std);
    return p_value;

  #edef

  #############################################################################

#eclass

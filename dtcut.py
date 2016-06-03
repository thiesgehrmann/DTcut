import numpy as np;
import copy;
import sys;


def prepare_tree(T):
  """
  The output of the linkage is confusing, and doesn't contain the leaf nodes.
  We add them here, and also do some pre-processing so that it is easy to keep
  track of the leaf nodes that descend from each internal node (an optimization).
  """
    
    # The prepared tree
    # Nodes have values (child_left_id, child_right_id, height, node_leaves, p_value, visited)
  P = [ [None, None, 0, set([i]), None, False] for i in xrange(len(T) + 1) ];
  
  for t in T:
    child_left_id, child_right_id, height, leaf_count = t;
    
    child_left  = P[int(child_left_id)];
    child_right = P[int(child_right_id)];
    
    node_leaves = child_right[3] | child_left[3];
    
    P.append( [int(child_left_id), int(child_right_id), height, node_leaves, None, False] );
  #efor
  
  return P;
#edef  

###############################################################################

class DTCUT_test:
  """
  This is the standard description of a test class.
  The test has access to *all* the raw data, and the tree structure in the DTCUT class.

  You can access any node in a tree with:
    self.T.T[i_node] = [ left_child(int), right_child(int), node_height(float), node_leaves(set), visited(bool) ]

  You can access the values for any leaf:
    self.X[leaf_id] = [ value ]*n (floats)

  You can access the labels for any leaf:
    self.labels[leaf_id] = label(bool)

  You can also access the ID of any leaf:
    self.IDS[leaf_id] = id (probably string)

  In one's own testing class, one only has to define the statistical test.
  Therefore, the structure of your test class code should be:

    from dtcut import DTCUT_test;
    class my_test(DTCUT_test):
      def test(self, i_node):
        p_value = # The magic statistical testing code
        return p_value;
  """

    # Tree
  T =  None;

    # Should we re-run the statistical test when the node is re-visited?
  recompute = False;
  p_thresh  = None;

  #############################################################################

  def __init__(self, T, p_thresh):
    self.T        = T;
    self.p_thresh = p_thresh
    self.prepare_data();
  #edef

  #############################################################################

  def prepare_data(self):
    # Make this function format the labels or data into whatever you need.

    return None;
  #edef

  #############################################################################

  def correct(self, pval, factor):
    """
    A wrapper function for pvalue correction, incase we want to do something more fancy
    """
    return pval * factor;
  #edef

  #############################################################################

  def test(self, i_node, **kwargs):
    print "You must define a test function in your statistical test class";
    return 1.0;
  #edef

  #############################################################################

  def is_significant(self, pval, factor):
    """
    A wrapper to test if a pvalue with a given correction factor is significant or not
    """
    return (self.correct(pval, factor) < self.p_thresh)
  #edef

  #############################################################################

  def descend_rule(self, i_node, p_thresh, factor, min_set_size, max_set_size):
    
    children       = [ self.T.get_tree_node_left_child_id(i_node), self.T.get_tree_node_right_child_id(i_node) ];
    valid_children = [ child for child in children if ( (child is not None) and (len(self.T.get_tree_node_leaves(child)) >= min_set_size) ) ]
    p_value        = self.T.get_tree_node_p_value(i_node);
    n_leaves       = len(self.T.get_tree_node_leaves(i_node))

    node_is_significant = self.is_significant(p_value, factor);

    if n_leaves > max_set_size:
      print "parent too big, return children"
      return valid_children;
    #fi

    if not(node_is_significant):
      print "Parent not sig. return children"
      return valid_children;
    #fi
    
    more_significant_children = [];
    for child_id in children:
      if (child_id is not None) and (self.T.get_tree_node_p_value(child_id) < p_value) and ( len(self.T.get_tree_node_leaves(child_id)) >= min_set_size):
        more_significant_children.append(child_id);
      #fi
    #edef
    
    if not(more_significant_children) and node_is_significant:
      print "Parent is significant, and there are no more sig children"
      return None;
    else:
      print "Children are more significant than parent"
      return valid_children;
    #fi
  #edef

#eclass

###############################################################################

class DTCUT:

  T            = None;
  X            = None;
  IDS          = None;
  L            = None;
  stat         = None;
  min_set_size = None;
  max_set_size = None;

  ###############################################################################

  def __init__(self, T, X, IDS, L):
    self.T   = copy.deepcopy(T);
    self.X   = np.matrix(X);
    self.IDS = np.array(IDS);
    self.L   = np.array(L);
  #edef

  ###############################################################################

  def get_clusters(self, S):
    """
    Get a list of length equal to the number of leaves.
    Each element describes the cluster number the leaf is in.
    """

    C = np.zeros(( (len(self.T)+1)/2, ), dtype=int) - 1;

    for (i, i_node) in enumerate(S):
      cluster_node = self.T[i_node];
      C[list(cluster_node[3])] = i + 1;
    #efor

    return C;
  #edef

  ###############################################################################

  def get_clusters_info(self, S):
    """
    Return a list of lists describing, for each significant node:
    * The p-value of the node,
    * The height that node exists at,
    * The IDS of the leaves at that node
    """

    I = [];
    
    for (i, i_node) in enumerate(S):
      cluster_node = self.T[i_node];
      I.append( ( cluster_node[2], cluster_node[4], [ self.IDS[leaf] for leaf in cluster_node[3] ] ) );
    #efor
    
    return I;
  #edef

  ###############################################################################

  def correct(self, pval, factor):
    """
    A wrapper function for pvalue correction, incase we want to do something more fancy
    """

    return pval * factor;
  #edef

  ###############################################################################

  def test_tree(self, stat, min_set_size=1, max_set_size=sys.maxint):
    """
    Loop through the tree and test each node in a BFS manner.

    We perform a dynamic correction of the nodes in the tree.
    Whenever we find a significant node, we set it aside, and test it later with an
    updated correction factor.
    This procedure is iterative, and corrects all significant nodes with updated correction factors until we don't perform any more tests.

    This is the dynamic correction algorithm of:
    "Proteny: Discovering and visualizing statistically significant syntenic clusters at the proteome level"

    Input:
     * stat:         The statistical test class
     * min_set_size: The smallest size a cluster may be
     * max_set_size: The largest size a cluster may be
    """

    self.stat         = stat;
    self.min_set_size = min_set_size;
    self.max_set_size = max_set_size;

    correction_factor = 0;
    testing           = [ -1 ];

    while True:
      significant = [];
      last_correction_factor = correction_factor;
      for i_node in testing:
          # Find significant nodes for a particular branch
        S, tests_done = self.traverse_node(i_node, correction_factor);
          # And keep track of them
        significant.extend(S);
          # And also the number of tests that were done additionally.
        correction_factor = correction_factor + tests_done;
      #efor

        # If we haven't done any more tests since last time, we can quit.
      if correction_factor == last_correction_factor:
        break;
      #fi
      testing = significant;
      last_correction_factor = correction_factor;
    #ewhile

    return significant;
  #edef

  ###############################################################################

  def traverse_node(self, root, correction_factor):
    """
    Do a BFS for a given branch.

    Initially the queue contains only the root node.
    We add a node to a 'significant' stack iff 
      * It is significant
      * Its children are /less/ significant

    Otherwise, we add its children to the queue.

    Returns significant nodes and the number of tests done, also the updated tree (I need to stop doing that.
    """


      # The tree traversal queue. Start with root node
    Q = [ root ];

      # Found significant nodes
    S = [];

    tests_done = correction_factor;

    while len(Q) > 0:

        # Traverse the queue
      i_node = Q.pop(0);

      visited = self.visit_node(i_node, tests_done);
      tests_done = tests_done + visited;


        # 
      next_tests = self.stat.descend_rule(i_node, self.stat.p_thresh, tests_done, self.min_set_size, self.max_set_size);
      if next_tests is None:
        S.append(i_node);
        print "SIGNIFICANT: Node %d, pvalue %.10f, correction factor %d" % (i_node, self.get_tree_node_p_value(i_node), tests_done);
      else:
        Q.extend(next_tests);
      #fi
    #ewhile

    return S, (tests_done - correction_factor);
  #edef

  ###############################################################################

  def visit_node(self, i_node, tests_done):
    """
    Visit a given node.
    This means, check its p-value, and the p-value of its children.
    It returns
     * The p-value of the node
     * Whether this node was first time we visited this node
    """

    more_significant_children = False;
    nodes_visited = 0;

    child_left_id, child_right_id, node_height, node_leaves, node_pvalue, node_visited = self.get_tree_node(i_node);

      # Do not evaluate nodes smaller than a certain number of leaves
    if len(node_leaves) < self.min_set_size or len(node_leaves) > self.max_set_size:
      self.set_tree_node_visited(i_node, False);
      self.set_tree_node_p_value(i_node, 1.0);
      return 0;
    #fi

      # Get the p-value of the node
    if node_visited == False or self.stat.recompute:
      node_pvalue = self.stat.test(i_node, factor=tests_done);

      self.set_tree_node_visited(i_node, True);
      self.set_tree_node_p_value(i_node, node_pvalue);
      nodes_visited = nodes_visited + (1 - int(node_visited));
    #fi

      # Get p-value of the LEFT and RIGHT child node
    child_nodes = [ child_left_id, child_right_id ];
    for child_id in child_nodes:
      if child_id != None:
        if self.get_tree_node_p_value(child_id) == None:
          child_node_pvalue = self.stat.test(child_id, factor=tests_done+nodes_visited);
          self.set_tree_node_visited(child_id, True);
          self.set_tree_node_p_value(child_id, child_node_pvalue);
          nodes_visited = nodes_visited + 1;
        #fi
      #fi
    #efor

    print "===\np-value of parent (%d) = %f\np-value of left (%d) = %f\np-value of right (%d) = %f\n\n" % (i_node,         self.get_tree_node_p_value(i_node), 
                                                                                                           child_left_id,  self.get_tree_node_p_value(child_left_id),
                                                                                                           child_right_id, self.get_tree_node_p_value(child_right_id))

    return nodes_visited;
  #edef

  #############################################################################

  def get_tree_node(self, i_node):
    return self.T[i_node];
  #edef

  #############################################################################

  def get_tree_node_left_child_id(self, i_node):
    return self.get_tree_node(i_node)[0];
  #edef

  #############################################################################

  def get_tree_node_right_child_id(self, i_node):
    return self.get_tree_node(i_node)[1];
  #edef

  #############################################################################

  def get_tree_node_height(self, i_node):
    return self.get_tree_node(i_node)[2];
  #edef

  #############################################################################

  def get_tree_node_leaves(self, i_node):
    return self.get_tree_node(i_node)[3];
  #edef

  #############################################################################

  def get_tree_node_p_value(self, i_node):
    return self.get_tree_node(i_node)[4];
  #edef

  #############################################################################

  def set_tree_node_p_value(self, i_node, pvalue):
    self.T[i_node][4] = pvalue;
  #edef

  #############################################################################

  def get_tree_node_visited(self, i_node):
    return self.get_tree_node(i_node)[5];
  #edef

  #############################################################################

  def set_tree_node_visited(self, i_node, visited):
    self.T[i_node][5] = visited;
  #edef

  #############################################################################

  def get_data_values(self, object_ids):
    if object_ids is None:
      return self.X;
    else:
      return self.X[object_ids];
    #fi
  #edef

  #############################################################################

  def set_data_values(self, object_ids, data):
    if object_ids is None:
      self.X = data;
    else:
      self.X[object_ids] = data;
    #fi
  #edef

  #############################################################################

  def get_labels(self, object_ids):
    if object_ids is None:
      return self.L;
    else:
      return self.L[object_ids];
    #fi
  #edef

  #############################################################################

  def set_labels(self, object_ids, labels):
    if object_ids is None:
      self.L = labels;
    else:
      self.L[object_ids] = labels;
    #fi

#eclass



def prepare_tree(T):
  """
  The output of the linkage is confusing, and doesn't contain the leaf nodes.
  We add them here, and also do some pre-processing so that it is easy to keep
  track of the leaf nodes that descend from each internal node (an optimization
  """
    
    # The prepared tree
    # Nodes have values (child_left_id, child_right_id, height, node_leaves, p_value)
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

class DT:

  T            = None;
  stat_func    = None;
  p_thresh     = None;
  min_set_size = None;

  ###############################################################################

  def __init__(self, T):
    self.T = T[:];
  #edef

  ###############################################################################

  def get_clusters_info(self, C, IDS):
    """
    For a given set of clusters, produce a table with the important information from these clusters
    """

    info = [];
    for c in C:
      cluster_node = self.T[c];
      info.append( ( cluster_node[2], cluster_node[4], [ IDS[leaf] for leaf in cluster_node[3] ] ));
    #efor
    return info;
  #efor

  ###############################################################################

  def correct(self, pval, factor):
    """
    A wrapper function for pvalue correction, incase we want to do something more fancy
    """

    return pval * factor;
  #edef

  ###############################################################################

  def test_tree(self, stat_func, p_thresh, min_set_size):
    """
    Loop through the tree and test each node in a BFS manner.

    We perform a dynamic correction of the nodes in the tree.
    Whenever we find a significant node, we set it aside, and test it later with an
    updated correction factor.
    This procedure is iterative, and corrects all significant nodes with updated correction factors until we don't perform any more tests.

    This is the dynamic correction algorithm of:
    "Proteny: Discovering and visualizing statistically significant syntenic clusters at the proteome level"
    """

    self.stat_func    = stat_func;
    self.p_thresh     = p_thresh;
    self.min_set_size = min_set_size;

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

      pval, children, children_more_significant, visited = self.visit_node(i_node);
      tests_done = tests_done + visited;

      if not(children_more_significant) and self.correct(pval, tests_done) < self.p_thresh:
        S.append(i_node);
        print "SIGNIFICANT: Node %d, pvalue %.10f, correction factor %d" % (i_node, pval, tests_done);
      else:
        Q.extend(children);
      #fi
    #ewhile

    return S, (tests_done - correction_factor);
  #edef

  ###############################################################################

  def visit_node(self, i_node):
    """
    Visit a given node.
    This means, check its p-value, and the p-value of its children.
    It returns
     * The p-value of the node
     * the list of children to test, 
     * Whether it's children are more significant
     * Whether this node was first time we visited this node
    """

    more_significant_children = False;
    children                  = [];

    child_left_id, child_right_id, node_height, node_leaves, node_pvalue, node_visited = self.T[i_node];

      # Do not evaluate nodes smaller than a certain number of leaves
    if len(node_leaves) < self.min_set_size:
      return 1.0, [], False, False;
    #fi

      # Get the p-value of the node
    if node_visited == False:
      node_pvalue        = self.stat_func(i_node, node_leaves);
      self.T[i_node][-1] = True;
      self.T[i_node][-2] = node_pvalue;
    #fi

      # Get p-value of the LEFT child node
    child_nodes = [ child_left_id, child_right_id ];
    for child_id in child_nodes:
      if child_id != None:
        children.append(child_id);
        child_child_left_id, child_child_right_id, child_height, child_node_leaves, child_node_pvalue, child_node_visited = self.T[child_id];
        if child_node_pvalue == None:
          child_node_pvalue = self.stat_func(child_id, child_node_leaves);
          self.T[child_id][-1] = child_node_visited;
          self.T[child_id][-2] = child_node_pvalue;
        #fi
        if node_pvalue > child_node_pvalue and len(child_node_leaves) >= self.min_set_size:
          more_significant_children = True;
        #fi
      #fi
    #efor

    return node_pvalue, children, more_significant_children, 1 - node_visited;
  #edef

#eclass


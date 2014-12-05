#!/bin/python

  # Standard python
import imp;
import sys;

  # Additional libraries
import matplotlib.pyplot as plt;
import numpy as np;
import fastcluster as fc;

  # Specific functions
from ibidas import *;
from scipy.cluster.hierarchy import dendrogram;
from scipy.spatial.distance import pdist;

###############################################################################

def get_data(cuffdiff_combine_file, labels_file):
  """
  Just some stuff for me, formatting the data
  """

  X = Load(cuffdiff_combine_file);
  P = Read(labels_file, format='tsv2')();

    # Just the individual conditions, and the slices that represent them (some are double)
  N = dict([ (n.split('_')[int(n.split('_', 3)[-1])-1],n) for n in X.Names if ('_value_' in n) ])

    # Get the expression values per test_id
  X = X.Get(_.test_id, *N.values()) / (('test_id',) + tuple(N.keys()));
  X = X.Cast(X.Type.getSubType(0), *(['real64']*len(N)))

    # Get rid of the proteins with all-zero expression values
  X = X[X.Get(tuple(X.Names[1:])).To(_.data, Do=_.Each(lambda x: sum(x))) != 0]

    # Get the IDS
  IDS = X.test_id();

  X_mat = zip(*X.Without(_.test_id)());
  L = [ True if prot in P else False for prot in IDS ];

  return IDS, X_mat, L;
#edef

def read_data(gene_ids, expression_file, labels_file):

  

  return IDS, X, L;
#edef
  

###############################################################################

def correct(pval, factor):
  """
  A wrapper function for pvalue correction, incase we want to do something more fancy
  """

  return pval * factor;
#edef

###############################################################################

def test_tree(T, stat_test, p_thresh, min_set_size=None):
  """
  Loop through the tree and test each node in a BFS manner.

  We perform a dynamic correction of the nodes in the tree.
  Whenever we find a significant node, we set it aside, and test it later with an
  updated correction factor.
  This procedure is iterative, and corrects all significant nodes with updated correction factors until we don't perform any more tests.

  This is the dynamic correction algorithm of:
  "Proteny: Discovering and visualizing statistically significant syntenic clusters at the proteome level"
  """
  correction_factor = 0;
  testing           = [ -1 ];

  while True:
    significant = [];
    last_correction_factor = correction_factor;
    for i_node in testing:
        # Find significant nodes for a particular branch
      S, tests_done = traverse_node(T, i_node, stat_test, correction_factor, p_thresh, min_set_size);
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
#efor
      

###############################################################################

def traverse_node(T, root, stat_test, correction_factor, p_thresh, min_set_size):
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

    pval, children, children_more_significant, visited = visit_node(T, i_node, stat_test, min_set_size);
    tests_done = tests_done + visited;

    if not(children_more_significant) and correct(pval, tests_done) < p_thresh:
      S.append(i_node);
      print "SIGNIFICANT: Node %d, pvalue %.10f, correction factor %d" % (i_node, pval, tests_done);
    else:
      Q.extend(children);
    #fi
  #ewhile

  return S, (tests_done - correction_factor);
#edef

###############################################################################

def prepare_tree(T):
  """
  The output of the linkage is confusing, and doesn't contain the leaf nodes.
  We add them here, and also do some pre-processing so that it is easy to keep
  track of the leaf nodes that descend from each internal node (an optimization
  """

    # The prepared tree
    # Nodes have values (child_left_id, child_right_id, height, node_leaves, p_value)
  P = [ (None, None, 0, set([i]), None, False) for i in xrange(len(T) + 1) ];

  for t in T:
    child_left_id, child_right_id, height, leaf_count = t;

    child_left  = P[int(child_left_id)];
    child_right = P[int(child_right_id)];

    node_leaves = child_right[3] | child_left[3];

    P.append( (int(child_left_id), int(child_right_id), height, node_leaves, None, False) );
  #efor

  return P;
#edef

###############################################################################

def visit_node(T, i_node, stat_test, min_set_size):
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

  child_left_id, child_right_id, node_height, node_leaves, node_pvalue, node_visited = T[i_node];

    # Do not evaluate nodes smaller than a certain number of leaves
  if len(node_leaves) < min_set_size:
    return 1, [], False, False;
  #fi

    # Get the p-value of the node
  if node_visited == False:
    node_pvalue = stat_test(i_node, node_leaves);
    T[i_node]   = ( child_left_id, child_right_id, node_height, node_leaves, node_pvalue, True );
  #fi

    # Get p-value of the LEFT child node
  if child_left_id != None:
    children.append(child_left_id);
    child_child_left_id, child_child_right_id, child_height, child_node_leaves, child_node_pvalue, child_node_visited = T[child_left_id];
    if child_node_pvalue == None:
      child_node_pvalue = stat_test(child_left_id, child_node_leaves);
      T[child_left_id] = ( child_child_left_id, child_child_right_id, child_height, child_node_leaves, child_node_pvalue, child_node_visited );
    #fi
    if node_pvalue > child_node_pvalue and len(child_node_leaves) >= min_set_size:
      more_significant_children = True;
    #fi
  #fi

    # Get p-value of the RIGHT child node
  if child_right_id != None:
    children.append(child_right_id);
    child_child_left_id, child_child_right_id, child_height, child_node_leaves, child_node_pvalue, child_node_visited = T[child_right_id];
    if child_node_pvalue == None:
      child_node_pvalue = stat_test(child_right_id, child_node_leaves);
      T[child_right_id] = ( child_child_left_id, child_child_right_id, child_height, child_node_leaves, child_node_pvalue, child_node_visited );
    #fi
    if node_pvalue > child_node_pvalue and len(child_node_leaves) >= min_set_size:
      more_significant_children = True
    #fi
  #fi

  return node_pvalue, children, more_significant_children, 1 - node_visited;

#edef

###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

  cuffdiff_combine = '/data/tmp/thiesgehrmann/rnaseq_schco3/cuffdiff_combine.dat';
  labels_file      = '/home/nfs/thiesgehrmann/groups/w/phd/tasks/cazyme_predicted/predicted_cazymes.txt';
  labels_file      = '/home/nfs/thiesgehrmann/groups/w/phd/data/schco3/go_annots/3043.txt';
  distance_metric  = 'correlation';
  linkage_method   = 'complete';
  

  expression_file = sys.argv[1];
  distance_metric = sys.argv[2];
  linkage_method  = sys.argv[3];
  stat_test_file  = sys.argv[4];
  output_dir      = sys.argv[5];
  label_files     = sys.argv[6:];

    # Import the statistical test
  TEST = imp.load_source('Test', stat_test_file);

    # Get the expression data matrix;
  print "Reading Data";
  IDS, X, L = ct.get_data(cuffdiff_combine, label_files);

    # Calculate the distances
  print "Calculating distances";
  Y = ct.pdist(X, distance_metric);
  
    # Calculate a dendrogram
  print "Constructing dendrogram";
  D = fc.linkage(Y, linkage_method, distance_metric);

    # Prepare the tree
  T = ct.prepare_tree(D);

  for i, label_set in enumerate(L):
    current_tree = T[:];
    print "TESTING label set: %s" % label_files[i];

      # Initialize the statistical test
    stat_test = Test(current_tree, IDS, X, label_set);
    stat_func = lambda i_node, subset: stat_test.test(i_node, subset);

      # Test the tree, getting significant nodes
    S = ct.test_tree(current_tree, stat_func, 0.05, min_set_size=20)

  #efor

#fi


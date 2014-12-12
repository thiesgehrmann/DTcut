#!/bin/env python

  # This package stuff
import dtcut as dtcut;
import cluster_figure as cf;

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


def read_data(ids_file, values_file, labels_files):
  """
  Read the input files.
    ids_file:    A list of object identifiers
    values_file: A vector describing each object
    labels_file: A list of files which contain lists of objects with a given annotation
  """

  IDS = Read(ids_file, format='tsv2').Get(0)();
  X   = zip(*Read(values_file).Cast('real64')());
  L   = [];

  for l_file in labels_files:
    labelled = Read(l_file, format='tsv2')();
    L.append([ True if id in labelled else False for id in IDS ]);
  #efor

  return IDS, X, L;
#edef


###############################################################################
###############################################################################
###############################################################################

def usage(arg0):

  print "%s: General purpose dynamic tree cutter based on a statistical test" % arg0;
  print " Usage: %s <ids_file> <values_file> <distance_metric> <linkage_method> <stat_test> <p-value_threshold> <min_set_size> <output_dir> <output_figs> <label_file_1> [label_file_2] ... [label_file_n]" % arg0;
  print "";
  print " Parameters";
  print "  ids_file:        The list of element IDS for your objects. For example, gene IDS";
  print "  values_file:     A description of each object, sorted identically to <ids_file>. For example, expression matrix.";
  print "  distance_metric: The metric used to calculate distances. For example, 'correlation' or 'euclidean'.";
  print "  linkage_method:  The linkage algorithm to use in clustering. E.g. 'complete'.";
  print "  stat_test:       The file describing the statistical test. See the example' EnrichmentTest.py'.";
  print "  p-value thresh:  The p-value threshold to use.";
  print "  min_set_size:    The minimum size a set must have in order to consider it.";
  print "  output_dir:      The output directory to put results.";
  print "  output_figs:     Output a figure describing the expression of the genes in each cluster. (True or False)";
  print "  label_file(s):   A file containing a list of genes that have a particular label";

#edef

if __name__ == '__main__':

  if len(sys.argv) < 11:
    usage(sys.argv[0]);
    sys.exit(1);
  #fi

  ids_file        = sys.argv[1];
  values_file     = sys.argv[2];
  distance_metric = sys.argv[3];
  linkage_method  = sys.argv[4];
  stat_test_file  = sys.argv[5];
  pvalue_thresh   = float(sys.argv[6]);
  min_set_size    = int(sys.argv[7]) if sys.argv[7] != 'None' else None;
  output_dir      = sys.argv[8];
  output_figs     = True if sys.argv[9] == 'True' else False;
  labels_files    = sys.argv[10:];

    # Import the statistical test
  TEST = imp.load_source('Test', stat_test_file);

    # Get the expression data matrix;
  print "Reading Data";
  IDS, X, L = read_data(ids_file, values_file, labels_files);

    # Calculate the distances
  print "Calculating distances";
  Y = pdist(X, distance_metric);
  
    # Calculate a dendrogram
  print "Constructing dendrogram";
  D = fc.linkage(Y, linkage_method, distance_metric);

    # Prepare the tree
  T = dtcut.prepare_tree(D);

  for i, label_set in enumerate(L):
    current_tree = dtcut.DTCUT(T);

    print "TESTING label set: %s" % labels_files[i];

      # Initialize the statistical test
    stat_test = TEST.Test(current_tree, IDS, X, label_set);
    stat_func = lambda i_node: stat_test.test(i_node);

      # Test the tree, getting significant nodes
    S = current_tree.test_tree(stat_func, pvalue_thresh, min_set_size);

      # Get information about these significant nodes
    clusters = current_tree.get_clusters(S);
    info     = current_tree.get_clusters_info(S, IDS);

      # Export the data
    output_prefix = '%s/%s' % (output_dir, labels_files[i].split('/')[-1]);
    Export(Rep(info) / ('height', 'p_value', 'members'), '%s.dtcut.tsv' % (output_prefix));
    if output_figs == True:
      cf.draw_clusters(X, D, clusters, info, '%s.dtcut.png' % (output_prefix));
    #fi
    
  #efor

#fi


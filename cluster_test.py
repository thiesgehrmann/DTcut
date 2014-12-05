#!/bin/env python

  # This package stuff
import dt as dt;

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

def read_data(gene_ids, expression_file, labels_files):

  X = Load(expression_file);
  P = Read(labels_files[0], format='tsv2')();
    
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

###############################################################################
###############################################################################
###############################################################################

def usage(arg0):

  print "%s: General purpose dynamic tree cutter based on a statistical test" % arg0;
  print " Usage: %s <expression_file> <distance_metric> <linkage_method> <stat_test> <output_dir> <label_file_1> [label_file_2] ... [label_file_n]" % arg0;
#edef

if __name__ == '__main__':

  cuffdiff_combine = '/data/tmp/thiesgehrmann/rnaseq_schco3/cuffdiff_combine.dat';
  labels_file      = '/home/nfs/thiesgehrmann/groups/w/phd/tasks/cazyme_predicted/predicted_cazymes.txt';
  labels_file      = '/home/nfs/thiesgehrmann/groups/w/phd/data/schco3/go_annots/3043.txt';
  distance_metric  = 'correlation';
  linkage_method   = 'complete';
  
  if len(sys.argv) < 9:
    usage(sys.argv[0]);
    sys.exit(1);
  #fi

  expression_file = sys.argv[1];
  distance_metric = sys.argv[2];
  linkage_method  = sys.argv[3];
  stat_test_file  = sys.argv[4];
  pvalue_thresh   = float(sys.argv[5]);
  min_set_size    = int(sys.argv[6]) if sys.argv[6] != 'None' else None;
  output_dir      = sys.argv[7];
  label_files     = sys.argv[8:];

    # Import the statistical test
  TEST = imp.load_source('Test', stat_test_file);

    # Get the expression data matrix;
  print "Reading Data";
  IDS, X, L = get_data(cuffdiff_combine, labels_file);
  L = [L];

    # Calculate the distances
  print "Calculating distances";
  Y = pdist(X, distance_metric);
  
    # Calculate a dendrogram
  print "Constructing dendrogram";
  D = fc.linkage(Y, linkage_method, distance_metric);

    # Prepare the tree
  T = dt.prepare_tree(D);

  for i, label_set in enumerate(L):
    current_tree = dt.DT(T);

    print "TESTING label set: %s" % label_files[i];

      # Initialize the statistical test
    stat_test = TEST.Test(current_tree, IDS, X, label_set);
    stat_func = lambda i_node, subset: stat_test.test(i_node, subset);

      # Test the tree, getting significant nodes
    S = current_tree.test_tree(stat_func, pvalue_thresh, min_set_size);
    print S;
    
  #efor

#fi


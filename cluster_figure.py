  #Additional libraries
import numpy as np;
import matplotlib.pyplot as plt;
import matplotlib.patches as mpatches

  # Specific functions
from scipy.cluster.hierarchy import dendrogram;

###############################################################################

def inverse_indexes(L):
  """
  Inverse an index list.
  So, rather than L[i] = j, you get L'[j] = i;
  """

  X = np.zeros((len(L),));

  for i, l in enumerate(L):
    X[l] = i;
  #efor

  return X;
#edef

###############################################################################

def draw_dendrogram(D, clusters, nclusters, ax):
  """
  Annotate the output of the dendrogram produced by scipy to contain the location of our discovered clusters.
  """

    # Loop through available colors
  colors = 'bgrcmy';
  used_colors = '';

    # Draw the deendrogram
  dend = dendrogram(D, color_threshold=0, no_labels=True)
  leaves = inverse_indexes(dend['leaves']);

    # Now add the annotations for each cluster
  for clust in xrange(1, nclusters+1):
      # Find the range of positions of leaves in the cluster
    leaf_locations = leaves[clusters == clust];
    color          = colors[clust % len(colors)];

      # Keep track of the colors used
    used_colors += color;

      # The leaves represent the locations of the cluster in the figure
    xmin = min(leaf_locations) * 10;
    xmax = max(leaf_locations) * 10;
    ymin = 0;
    ymax = 2;

      # Draw a rectangle around this cluster
    patch = mpatches.Rectangle( (xmin, ymin), xmax - xmin, ymax, alpha=0.4, color=color);
    ax.add_patch(patch);

      # Also add the cluster ID so that it can be linked to the graph
    plt.annotate('%d' % clust, ( (xmin+xmax)/2, ymin), xytext=(0, -4), textcoords='offset points', va='top', ha='center');

      # Add a line so that we can see how many leaves are in the cluster
    plt.plot((xmin, xmax), (ymin, ymin), '%c-' % color, linewidth=10);
  #efor

    # Return the colors used
  return used_colors;
#edef

###############################################################################

def draw_clusters(X, D, clusters, info, filename):
  """
    Draw a figure with a dendrogram with each significant cluster highlighted.
    Also draw the feature vectors for each leaf within a cluster.
  """
  plt.cla();


    # Count the number of clusters we found
  nclusters = max(clusters);

    # Knowing the number of clusters, design the figure size
  dendrogram_rows = 3; # The dendrogram needs more space
  figdim_x        = 4; # Two graphs per line
  figdim_y        = int(np.ceil(float(nclusters)/figdim_x)) + dendrogram_rows;
  fig             = plt.figure(figsize=(16*figdim_x, 3.0*figdim_y));

   # Make sure that the data is in numpy format.
  X = np.matrix(X);

    # Plot the dendrogram
  ax          = plt.subplot2grid((figdim_y, figdim_x), (0,0), colspan=figdim_x, rowspan=dendrogram_rows);
  used_colors = draw_dendrogram(D, clusters, nclusters, ax);

  plt.title('Significant clusters found in dendrogram: %s' % filename);

    # Now plot a graph showing the feature vector for each cluster
  xticks = np.arange(X.shape[1]);
  for cluster_k in xrange(1,nclusters+1):
      # Calculate the index of this subplot
    plot_id = cluster_k + (figdim_x * dendrogram_rows);
    plt.subplot(figdim_y, figdim_x, plot_id);

      # Get the additional information at this plot
    height, pvalue, leaves = info[cluster_k-1]; 

      # Find the objects in this cluster
    c_mem   = (clusters == cluster_k);
      # And select only these profiles
    ty = X[c_mem,:];
    tx = np.tile(xticks, ty.shape[0]).reshape(ty.shape[0],len(xticks));
    
      # Plot these with an alpha of 0.05
    plt.plot(tx.T, ty.T, 'k', alpha=0.05);
      # Then plot the average on top of it with the color used in the dendrogram
    plt.plot(tx.T, ty.mean(0).T, used_colors[cluster_k-1], linewidth=1);
    plt.xlim((0, ty.shape[1]));
    plt.title('Cluster %d (%d objects) (pvalue: %f) (height: %f)' % (cluster_k, ty.shape[0], pvalue, height), color=used_colors[cluster_k-1]);
  #efor

  plt.savefig(filename, dpi=200);

#edef


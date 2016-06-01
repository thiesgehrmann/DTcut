from scipy.stats import chi2_contingency; 
from dtcut import DTCUT_test;

class Test(DTCUT_test):

  def test(self, i_node):
    """
    Do a simple enrichment using the Chi^2 test.
    """
    
    subset        = self.get_tree_node_leaves(i_node)
    subset_labels = [ self.labels[sample] for sample in subset ];
    
    a = sum(self.labels);
    b = sum(subset_labels);
    
    c = len(self.labels)   - a;
    d = len(subset_labels) - b;
    
    a += 1; b += 1; c += 1; d += 1;
    
    chi2, pvalue, dof, expected = chi2_contingency([[a,b],[c,d]]);

    return pvalue;
  #edef

#eclass

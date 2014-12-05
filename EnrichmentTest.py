from scipy.stats import chi2_contingency; 

class Test:

  T      = None;
  IDS    = None;
  X      = None;
  labels = None;

  #############################################################################

  def __init__(self, T, IDS, X, labels):
    self.T      = T;
    self.IDS    = IDS;
    self.X      = X;
    self.labels = labels;
  #edef

  #############################################################################

  def test(self, i_node, sub_set):
    """
    Do a simple enrichment using the Chi^2 test.
    """
    
    subset_labels = [ self.labels[sample] for sample in sub_set ];
    
    a = sum(self.labels);
    b = sum(subset_labels);
    
    c = len(self.labels)   - a;
    d = len(subset_labels) - b;
    
    a += 1; b += 1; c += 1; d += 1;
    
    chi2, pvalue, dof, expected = chi2_contingency([[a,b],[c,d]]);

    return pvalue;
  #edef

#eclass

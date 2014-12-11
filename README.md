DT-cut
======

Dynamic Test Cut: A general purpose tool to dynamically cut a tree based on a statistical test at each node.

A general purpose implementation of the dynamic tree cutting algorithm in [Proteny](https://github.com/thiesgehrmann/proteny#proteny).

![An example output from DT-cut](https://raw.githubusercontent.com/thiesgehrmann/DTcut/master/example_data/example_output/GO:0008152.tsv.dtcut.png)

## Dependencies

For the core library functionality, you must have installed:
* Numpy: http://www.numpy.org/
* Scipy: http://www.scipy.org/
* Fastcluster: https://pypi.python.org/pypi/fastcluster

If you want to use the provided execution wrapper 'cluster_test.py', then you must also install:
* IBIDAS: https://pypi.python.org/pypi/Ibidas
* Matplotlib: http://matplotlib.org/

## Example usage

Example data in /example_data is prepared from:

* http://genome-www.stanford.edu/yeast_stress/data/rawdata/complete_dataset.txt [Gasch et al. (2000) Mol. Biol. Cell. 11(12) 4241-4257](http://genome-www.stanford.edu/yeast_stress/gasch.pdf)
* http://downloads.yeastgenome.org/curation/literature/gene_association.sgd.gz

The example can be run with the following command:
```
./cluster_test.py example_data/gene_ids.tsv example_data/expression_matrix.tsv correlation complete EnrichmentTest.py 0.05 20 ./ True example_data/GO\:0008152.tsv
```

It produces output like that seen in file [/example_data/example_output/readme/GO:0008152.tsv.dtcut.tsv](https://raw.githubusercontent.com/thiesgehrmann/DTcut/master/example_data/example_output/GO:0008152.tsv.dtcut.tsv).


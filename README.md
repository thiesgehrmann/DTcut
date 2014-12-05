DT-cut
======

Dynamic Test Cut: A general purpose tool to dynamically cut a tree based on a statistical test at each node.

A general purpose implementation of the dynamic tree cutting algorithm in [Proteny](https://github.com/thiesgehrmann/proteny#proteny).

## Dependencies

* Numpy: http://www.numpy.org/
* Scipy: http://www.scipy.org/
* Fastcluster: https://pypi.python.org/pypi/fastcluster
* IBIDAS: https://pypi.python.org/pypi/Ibidas

## Example usage

Example data in /example_data is prepared from:

* http://genome-www.stanford.edu/yeast_stress/data/rawdata/complete_dataset.txt
* http://downloads.yeastgenome.org/curation/literature/gene_association.sgd.gz

The example can be run with the following command:
```
./cluster_test.py example_data/gene_ids.tsv example_data/expression_matrix.tsv  correlation complete EnrichmentTest.py 0.05 20 ./ example_data/GO\:0008152.tsv
```

It produces output like:



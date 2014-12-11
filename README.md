DT-cut
======

Dynamic Test Cut: A general purpose tool to dynamically cut a tree based on a statistical test at each node.

A general purpose implementation of the dynamic tree cutting algorithm in [Proteny](https://github.com/thiesgehrmann/proteny#proteny).

![An example output from DT-cut](https://raw.githubusercontent.com/thiesgehrmann/DTcut/master/readme/GO:0008152.tsv.dtcut.png)

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

* http://genome-www.stanford.edu/yeast_stress/data/rawdata/complete_dataset.txt
* http://downloads.yeastgenome.org/curation/literature/gene_association.sgd.gz

The example can be run with the following command:
```
./cluster_test.py example_data/gene_ids.tsv example_data/expression_matrix.tsv correlation complete EnrichmentTest.py 0.05 20 ./ True example_data/GO\:0008152.tsv
```

It produces output like that seen in file [/readme/example_output.tsv](https://raw.githubusercontent.com/thiesgehrmann/DTcut/master/readme/GO:0008152.tsv.dtcut.tsv).

|Height|P-value|Members|
|------|-------|-------|
|0.406364473724 | 8.89184904984e-18 | ['YJR026W', 'YJR027W', 'YJR028W', 'YJR029W', 'YER160C', 'YAR009C', 'YER172C', 'YML039W', 'YML040W', 'YML045W', 'YHR214C-B', 'YMR045C', 'YMR046C', 'YBR012W-A', 'YMR050C', 'YMR051C', 'YBL005W', 'YBL005W-A', 'YER138C', 'YBR012W-B'] |
|0.92738340635  | 1.4995446452e-06  | ['YML019W', 'YFL026W', 'YMR011W', 'YAL007C', 'YLR404W', 'YAL012W', 'YFL037W', 'YKL077W', 'YPL244C', 'YLR380W', 'YDL096C', 'YEL071W', 'YGR172C', 'YPR074C', 'YHR084W', 'YDR372C', 'YER003C', 'YIL123W', 'YDL145C', 'YPL012W', 'YMR296C', 'YGL202W', 'YLR153C', 'YDR393W', 'YDR261C', 'YMR307W', 'YPL028W', 'YPL030W', 'YOR175C', 'YKR090W', 'YGL225W', 'YJL192C', 'YAR003W', 'YLR179C', 'YDR158W', 'YGL039W', 'YPL048W', 'YMR305C', 'YLR192C', 'YAR002C-A', 'YLR195C', 'YLR452C', 'YDR432W', 'YDR433W', 'YIL015W', 'YBR133C', 'YLL018C', 'YNL109W', 'YDR189W', 'YGL008C', 'YGL012W', 'YMR123W', 'YGR014W', 'YGR124W', 'YPL094C', 'YDL095W', 'YGL026C', 'YER090W', 'YLR036C', 'YPL246C', 'YBR166C', 'YKL128C', 'YBL032W', 'YBR171W', 'YJR040W', 'YNL066W', 'YOR254C', 'YMR149W', 'YJL145W', 'YHL003C', 'YML046W', 'YNL078W', 'YML085C', 'YER119C-A', 'YML124C', 'YDL236W', 'YER122C', 'YCL045C', 'YGL200C', 'YNL291C', 'YBR205W', 'YKL165C', 'YJR139C', 'YBR210W', 'YNR021W', 'YER146W', 'YGR078C', 'YDR410C', 'YPL256C', 'YPL163C', 'YLR300W', 'YGL097W', 'YLR379W', 'YJL196C', 'YMR199W', 'YPR145W', 'YOR311C', 'YJR105W', 'YLR059C', 'YBR243C', 'YKL209C', 'YKL211C', 'YBR252W', 'YMR091C', 'YIL053W', 'YLR328W', 'YHR024C', 'YDR309C', 'YNR067C', 'YEL001C', 'YDL052C', 'YDL055C', 'YMR215W', 'YLR342W', 'YLR088W', 'YML126C', 'YIL039W', 'YOL007C', 'YJR143C', 'YGL148W', 'YBR023C'] |
|1.07143779724  | 3.61612238635e-22 | ['YLR359W', 'YKL120W', 'YBR291C', 'YBR294W', 'YER069W', 'YIL094C', 'YKL001C', 'YOR375C', 'YDR354W', 'YCL009C', 'YLR364W', 'YIL116W', 'YJR137C', 'YCL030C', 'YNL220W', 'YOL064C', 'YKR069W', 'YMR300C', 'YEL063C', 'YMR062C', 'YGR204W', 'YCL018W', 'YJL060W', 'YLR302C', 'YLR303W', 'YFR030W', 'YOR184W', 'YAR015W', 'YDR034C', 'YGL234W', 'YNR050C', 'YJR109C', 'YER052C', 'YBR248C', 'YHR018C', 'YOR200W', 'YOR201C', 'YOR202W', 'YOR203W', 'YPR167C', 'YDR035W', 'YDL182W', 'YGL009C', 'YDL059C', 'YJR010W', 'YOR222W', 'YIL074C', 'YOL058W', 'YER081W', 'YOL140W', 'YDL198C', 'YIR034C'] |



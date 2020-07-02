# vembrane: variant filtering using python expressions

## Authors

* Jan Forster (@jafors)
* Till Hartmann (@tedil)
* Elias Kuthe (@eqt)
* Johannes Köster (@johanneskoester)
* Christopher Schröder (@christopher-schroeder)

## Examples

```vembrane test.bcf 'ANN["Gene_Name"] == "CDH2" and ANN["Annotation_Impact"] == "HIGH"```

# vembrane: variant filtering using python expressions

## Authors

* Jan Forster (@jafors)
* Till Hartmann (@tedil)
* Elias Kuthe (@eqt)
* Johannes Köster (@johanneskoester)
* Christopher Schröder (@christopher-schroeder)

## Examples

Remove annotations where gene equals "CHD2" and its impact is "HIGH"
```vembrane test.bcf 'ANN["Gene_Name"] == "CDH2" and ANN["Annotation_Impact"] == "HIGH"```

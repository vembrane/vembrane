# vembrane: variant filtering using python expressions

## Authors

* Jan Forster (@jafors)
* Till Hartmann (@tedil)
* Elias Kuthe (@eqt)
* Johannes Köster (@johanneskoester)
* Christopher Schröder (@christopher-schroeder)

## Examples

Vembrane allows to simultaneously filter variants based on any `INFO` field, `CHROM`, `POS`, `REF`, `ALT`, `QUAL`, and the annotation field `ANN`. When filtering based on `ANN`, annotation entries are filtered first. If no annotation entry remains, the entire variant is deleted.

* Only keep annotations and variants where gene equals "CHD2" and its impact is "HIGH": 
  ```
  vembrane variants.bcf 'ANN["Gene_Name"] == "CDH2" and ANN["Annotation_Impact"] == "HIGH"'
  ```
* Remove variants with quality less than 30:
  ```
  vembrane variants.vcf 'QUAL >= 30'
  ```
* Only keep annotations and variants where feature (transcript) is ENST00000307301:
  ```
  vembrane variants.bcf 'ANN["Feature"] == "ENST00000307301"'
  ```
* Only keep annotations and variants where protein position is less than 10:
  ```
  vembrane variants.bcf 'ANN["Protein_position"] < 10'
  ```
* Only keep annotations and variants where consequence contains the word stream (mathing upstream and downstream):
  ```
  vembrane variants.vcf 're.search("stream", ANN["Consequence"])'
  ```

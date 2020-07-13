# vembrane: variant filtering using python expressions

## Authors

* Jan Forster (@jafors)
* Till Hartmann (@tedil)
* Elias Kuthe (@eqt)
* Johannes Köster (@johanneskoester)
* Christopher Schröder (@christopher-schroeder)
* Felix Mölder (@felixmoelder)

## Examples

Vembrane allows to simultaneously filter variants based on any `INFO` field, `CHROM`, `POS`, `REF`, `ALT`, `QUAL`, and the annotation field `ANN`. When filtering based on `ANN`, annotation entries are filtered first. If no annotation entry remains, the entire variant is deleted.

* Only keep annotations and variants where gene equals "CDH2" and its impact is "HIGH": 
  ```
  vembrane variants.bcf 'ANN["Gene_Name"] == "CDH2" and ANN["Annotation_Impact"] == "HIGH"'
  ```
* Only keep variants with quality at least 30:
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
* Only keep annotations and variants where consequence contains the word "stream" (matching "upstream" and "downstream"):
  ```
  vembrane variants.vcf 're.search("stream", ANN["Consequence"])'
  ```

## Custom ANN types

vembrane parses the following annotation fields to a custom type:
* (snpeff) `cDNA.pos / cDNA.length`, `CDS.pos / CDS.length` and `AA.pos / AA.length` are re-exposed as `cDNA`, `CDS` and `AA` respectively with attributes `pos` and `length`, e.g. can be accessed like this: `ANN["cDNA"].pos`
* `CLIN_SIG` is split at `'&'` into a list of entries

Any unknown annotation field will be left as is.

## Missing values in annotations

If a certain annotation field lacks a value, it will be replaced with the special value of `NA`. Comparing with this value will always result in `False`, e.g.
`ANN["cDNA"].pos > 0` will always evaluate to `False` *if* there was no value in the "cDNA.pos / cDNA.length" field of ANN (otherwise the comparison will be carried out with the usual semantics).
One way to handle optional values is by asserting that the field is not None, e.g `ID and "foo" in ID`.

## Development
### pre-commit hooks
Since we enforce code formatting with `black` by checking for that in CI, we can avoid "fmt" commits by ensuring formatting is done upon comitting changes:
1. make sure `pre-commit` is installed on your machine / in your env (should be available in pip, conda, archlinux repos, ...)
2. run `pre-commit install`. This will activate pre-commit hooks to your _local_ .git

Now when calling `git commit`, your changed code will be formatted with `black`, checked with`flake8`, get trailing whitespace removed and trailing newlines added (if needed)

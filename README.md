# vembrane: variant filtering using python expressions

Vembrane allows to simultaneously filter variants based on any `INFO` field, `CHROM`, `POS`, `REF`, `ALT`, `QUAL`, and the annotation field `ANN`. When filtering based on `ANN`, annotation entries are filtered first. If no annotation entry remains, the entire variant is deleted.


## Filter expression
The filter expression can be any valid python expression that evaluates to `bool`. However, functions and symbols available have been restricted to the following:

 * `any`, `all`
 * `min`, `max`, `sum`
 * `list`, `dict`, `set`, `tuple`
 * `zip`, `map`
 * `str`, `float`, `int`
 * Any function or symbol from [`math`](https://docs.python.org/3/library/math.html)
 * Regular expressions via [`re`](https://docs.python.org/3/library/re.html)

## Available fields
The following VCF fields can be accessed in the filter expression:

|Name|Type|Interpretation|Example expression|
|---|---|---|---|
|`INFO`|`Dict[str, Any¹]`| `INFO field -> Value`  | `INFO["DP"] > 0`|
|`ANN`| `Dict[str, Any²]`| `ANN field -> Value` | `ANN["Gene_Name"] == "CDH2"`|
|`CHROM`| `str` | Chromosome Name  |  `CHROM == "chr2"` |
|`POS`| `int` | Chromosomal position  | `24 < POS < 42`|
|`ID`| `str`  | Variant ID |  `ID == "rs11725853"` |
|`REF`| `str` |  Reference allele  | `REF == "A"` |
|`ALT`| `List[str]` |  Alternative alleles  | `"C" in ALT or ALT[0] == "G"`|
|`QUAL`| `float`  | Quality |  `QUAL >= 60` |
|`FILTER`|  |   |  |
|`FORMAT`|`Dict[str, Dict[str, Any¹]]`| `Format -> (Sample -> Value)` | `FORMAT["DP"][SAMPLES[0]] > 0` |
|`SAMPLES`|`List[str]`| `[Sample]`  |  `"Tumor" in SAMPLES` |

 ¹ depends on type specified in VCF header

 ² for the usual snpeff and vep annotations, custom types have been specified; any unknown ANN field will simply be of type `str`. If something lacks a custom parser/type, please consider filing an issue in the [issue tracker](https://github.com/vembrane/vembrane/issues).


## Examples

* Only keep annotations and variants where gene equals "CDH2" and its impact is "HIGH":
  ```
  vembrane 'ANN["Gene_Name"] == "CDH2" and ANN["Annotation_Impact"] == "HIGH"' variants.bcf
  ```
* Only keep variants with quality at least 30:
  ```
  vembrane 'QUAL >= 30' variants.vcf
  ```
* Only keep annotations and variants where feature (transcript) is ENST00000307301:
  ```
  vembrane 'ANN["Feature"] == "ENST00000307301"' variants.bcf
  ```
* Only keep annotations and variants where protein position is less than 10:
  ```
  vembrane 'ANN["Protein"].start < 10' variants.bcf
  ```
* Only keep variants where mapping quality is exactly 60:
  ```
  vembrane 'INFO["MQ"] == 60' variants.bcf
  ```
* Only keep annotations and variants where consequence contains the word "stream" (matching "upstream" and "downstream"):
  ```
  vembrane 're.search("stream", ANN["Consequence"])' variants.vcf
  ```
* Only keep annotations and variants where CLIN_SIG contains "pathogenic", "likely_pathogenic" or "drug_response":
  ```
  vembrane 'any(entry in ANN["CLIN_SIG"] for entry in ("pathogenic", "likely_pathogenic", "drug_response"))' variants.vcf
  ```

## Custom ANN types

vembrane parses the following annotation fields to a custom type:
* (snpeff) `cDNA.pos / cDNA.length`, `CDS.pos / CDS.length` and `AA.pos / AA.length` are re-exposed as `cDNA`, `CDS` and `AA` respectively with properties `start`, `end` and `length`, e.g. can be accessed like this: `ANN["cDNA"].start`
* (vep) `cDNA_position`, `CDS_position` and `Protein_position` are re-exposed as `cDNA`, `CDS` and `Protein` respectively with properties `start`, `end` and `length`, e.g. can be accessed like this: `ANN["cDNA"].start`
* `CLIN_SIG` is split at `'&'` into a list of entries

Any unknown annotation field will be left as is.

## Missing values in annotations

If a certain annotation field lacks a value, it will be replaced with the special value of `NA`. Comparing with this value will always result in `False`, e.g.
`ANN["cDNA"].start > 0` will always evaluate to `False` *if* there was no value in the "cDNA.pos / cDNA.length" field of (snpeff) ANN (otherwise the comparison will be carried out with the usual semantics).
One way to handle optional values is by asserting that the field is not None, e.g `ID and "foo" in ID`.

## Development
### pre-commit hooks
Since we enforce code formatting with `black` by checking for that in CI, we can avoid "fmt" commits by ensuring formatting is done upon comitting changes:
1. make sure `pre-commit` is installed on your machine / in your env (should be available in pip, conda, archlinux repos, ...)
2. run `pre-commit install`. This will activate pre-commit hooks to your _local_ .git

Now when calling `git commit`, your changed code will be formatted with `black`, checked with`flake8`, get trailing whitespace removed and trailing newlines added (if needed)

## Authors

* Marcel Bargull (@mbargull)
* Jan Forster (@jafors)
* Till Hartmann (@tedil)
* Johannes Köster (@johanneskoester)
* Elias Kuthe (@eqt)
* Felix Mölder (@felixmoelder)
* Christopher Schröder (@christopher-schroeder)

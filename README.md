# vembrane: variant filtering using python expressions

vembrane allows to simultaneously filter variants based on any `INFO` field, `CHROM`, `POS`, `REF`, `ALT`, `QUAL`, and the annotation field `ANN`. When filtering based on `ANN`, annotation entries are filtered first. If no annotation entry remains, the entire variant is deleted.

## `vembrane filter`

### Filter expression
The filter expression can be any valid python expression that evaluates to `bool`. However, functions and symbols available have been restricted to the following:

 * `all`, `any`
 * `abs`, `len`, `max`, `min`, `round`, `sum`
 * `enumerate`, `filter`, `iter`, `map`, `next`, `range`, `reversed`, `sorted`, `zip`
 * `dict`, `list`, `set`, `tuple`
 * `bool`, `chr`, `float`, `int`, `ord`, `str`
 * Any function or symbol from [`math`](https://docs.python.org/3/library/math.html)
 * Regular expressions via [`re`](https://docs.python.org/3/library/re.html)

### Available fields
The following VCF fields can be accessed in the filter expression:

|Name|Type|Interpretation|Example expression|
|---|---|---|---|
|`INFO`|`Dict[str, Any¹]`| `INFO field -> Value`  | `INFO["DP"] > 0`|
|`ANN`| `Dict[str, Any²]`| `ANN field -> Value` | `ANN["Gene_Name"] == "CDH2"`|
|`CHROM`| `str` | Chromosome Name  |  `CHROM == "chr2"` |
|`POS`| `int` | Chromosomal position  | `24 < POS < 42`|
|`ID`| `str`  | Variant ID |  `ID == "rs11725853"` |
|`REF`| `str` |  Reference allele  | `REF == "A"` |
|`ALT`| `str` |  Alternative allele³  | `ALT == "C"`|
|`QUAL`| `float`  | Quality |  `QUAL >= 60` |
|`FILTER`| `List[str]` | Filter tags | `"PASS" in FILTER` |
|`FORMAT`|`Dict[str, Dict[str, Any¹]]`| `Format -> (Sample -> Value)` | `FORMAT["DP"][SAMPLES[0]] > 0` |
|`SAMPLES`|`List[str]`| `[Sample]`  |  `"Tumor" in SAMPLES` |
|`INDEX`|`int`| `Index of variant in the file`  |  `INDEX < 10` |

 ¹ depends on type specified in VCF header

 ² for the usual snpeff and vep annotations, custom types have been specified; any unknown ANN field will simply be of type `str`. If something lacks a custom parser/type, please consider filing an issue in the [issue tracker](https://github.com/vembrane/vembrane/issues).

 ³ vembrane does not handle multi-allelic records itself. Instead, such files should be
 preprocessed by either of the following tools (preferably even before annotation):
 - [`bcftools norm -m-any […]`](http://samtools.github.io/bcftools/bcftools.html#norm)
 - [`gatk LeftAlignAndTrimVariants […] --split-multi-allelics`](https://gatk.broadinstitute.org/hc/en-us/articles/360037225872-LeftAlignAndTrimVariants)
 - [`vcfmulti2oneallele […]`](http://lindenb.github.io/jvarkit/VcfMultiToOneAllele.html)


### Examples

* Only keep annotations and variants where gene equals "CDH2" and its impact is "HIGH":
  ```
  vembrane filter 'ANN["Gene_Name"] == "CDH2" and ANN["Annotation_Impact"] == "HIGH"' variants.bcf
  ```
* Only keep variants with quality at least 30:
  ```
  vembrane filter 'QUAL >= 30' variants.vcf
  ```
* Only keep annotations and variants where feature (transcript) is ENST00000307301:
  ```
  vembrane filter 'ANN["Feature"] == "ENST00000307301"' variants.bcf
  ```
* Only keep annotations and variants where protein position is less than 10:
  ```
  vembrane filter 'ANN["Protein"].start < 10' variants.bcf
  ```
* Only keep variants where mapping quality is exactly 60:
  ```
  vembrane filter 'INFO["MQ"] == 60' variants.bcf
  ```
* Only keep annotations and variants where consequence contains the word "stream" (matching "upstream" and "downstream"):
  ```
  vembrane filter 're.search("(up|down)stream", ANN["Consequence"])' variants.vcf
  ```
* Only keep annotations and variants where CLIN_SIG contains "pathogenic", "likely_pathogenic" or "drug_response":
  ```
  vembrane filter 'any(entry in ANN["CLIN_SIG"] for entry in ("pathogenic", "likely_pathogenic", "drug_response"))' variants.vcf
  ```

### Custom `ANN` types
`vembrane` parses entries in the annotation field as outlined in [Types.md](Types.md)

### Missing values in annotations

If a certain annotation field lacks a value, it will be replaced with the special value of `NA`. Comparing with this value will always result in `False`, e.g.
`ANN["MOTIF_POS"] > 0` will always evaluate to `False` *if* there was no value in the "MOTIF_POS" field of ANN (otherwise the comparison will be carried out with the usual semantics).

*Explicitly* handling missing/optional values in INFO or FORMAT fields can be done by checking for NA, e.g.: `INFO["DP"] is NA`.

Handling missing/optional values in fields other than INFO or FORMAT can be done by checking for None, e.g `ID is not None`.

## `vembrane table`
In addition to the `filter` subcommand, vembrane (`≥ 0.5`) also supports writing tabular data with the `table` subcommand.
In this case, an expression which evaluates to `tuple` is expected, for example:
```
vembrane table 'CHROM, POS, 10**(-QUAL/10)', ANN["CLIN_SIG"] > table.tsv`.
```

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

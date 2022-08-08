[![CI](https://github.com/vembrane/vembrane/actions/workflows/main.yml/badge.svg)](https://github.com/vembrane/vembrane/actions/workflows/main.yml) [![DOI](https://zenodo.org/badge/276383670.svg)](https://zenodo.org/badge/latestdoi/276383670)

# vembrane: variant filtering using python expressions

vembrane allows to simultaneously filter variants based on any `INFO` or `FORMAT` field, `CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`, and the annotation field `ANN`. When filtering based on `ANN`, annotation entries are filtered first. If no annotation entry remains, the entire variant is deleted.

vembrane relies on [pysam](https://pysam.readthedocs.io/en/latest/) for reading/writing VCF/BCF files.

## Installation
vembrane is available in [bioconda](https://bioconda.github.io/) and can either be installed into an existing conda environment with `mamba install -c bioconda vembrane` or into a new named environment `mamba create -n environment_name -c bioconda vembrane`.
Alternatively, if you are familiar with git and [poetry](https://python-poetry.org/), clone this repository and run `poetry install`.

## `vembrane filter`

### Usage
vembrane takes two positional arguments: The filter expression and the input file; the latter may be omitted to read from `stdin` instead, making it easy to use vembrane in pipe chains.
```
usage: vembrane filter [options] expression [input vcf/bcf]

options:
  -h, --help            show this help message and exit
  --output OUTPUT, -o OUTPUT
                        Output file. If not specified, output is written to STDOUT.
  --output-fmt {vcf,bcf,uncompressed-bcf}, -O {vcf,bcf,uncompressed-bcf}
                        Output format.
  --annotation-key FIELDNAME, -k FIELDNAME
                        The INFO key for the annotation field. Defaults to "ANN".
  --aux NAME PATH, -a NAME PATH
                        Path to an auxiliary file containing a set of symbols.
  --keep-unmatched      Keep all annotations of a variant if at least one of them
                        passes the expression (mimics SnpSift behaviour).
  --preserve-order      Ensures that the order of the output matches that of the input. 
                        This is only useful if the input contains breakends (BNDs)
                        since the order of all other variants is preserved anyway.
```


### Filter expression
The filter expression can be any valid python expression that evaluates to `bool`. However, functions and symbols available have been restricted to the following:

 * `all`, `any`
 * `abs`, `len`, `max`, `min`, `round`, `sum`
 * `enumerate`, `filter`, `iter`, `map`, `next`, `range`, `reversed`, `sorted`, `zip`
 * `dict`, `list`, `set`, `tuple`
 * `bool`, `chr`, `float`, `int`, `ord`, `str`
 * Any function or symbol from [`math`](https://docs.python.org/3/library/math.html)
 * Any function from [`statistics`](https://docs.python.org/3/library/statistics.html)
 * Regular expressions via [`re`](https://docs.python.org/3/library/re.html)
 * custom functions:
   * `without_na(values: Iterable[Any]) -> Iterable[Any]` (keep only values that are not `NA`)
   * `replace_na(values: Iterable[Any], replacement: Any) -> Iterable[Any]` (replace values that are `NA` with some other fixed value)
   * genotype related:
     * `count_hom`, `count_het` , `count_any_ref`, `count_any_var`, `count_hom_ref`, `count_hom_var`
     * `is_hom`, `is_het`, `is_hom_ref` , `is_hom_var`
     * `has_ref`, `has_var`

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
|`INDEX`|`int`| Index of variant in the file  |  `INDEX < 10` |

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
  Using set operations, the same may also be expressed as:
  ```
  vembrane filter 'not {"pathogenic", "likely_pathogenic", "drug_response"}.isdisjoint(ANN["CLIN_SIG"])' variants.vcf
  ```

### Custom `ANN` types
`vembrane` parses entries in the annotation field as outlined in [Types.md](Types.md)

### Missing values in annotations

If a certain annotation field lacks a value, it will be replaced with the special value of `NA`. Comparing with this value will always result in `False`, e.g.
`ANN["MOTIF_POS"] > 0` will always evaluate to `False` *if* there was no value in the "MOTIF_POS" field of ANN (otherwise the comparison will be carried out with the usual semantics).

Since you may want to use the regex module to search for matches, `NA` also acts as an empty `str`, such that `re.search("nothing", NA)` returns nothing instead of raising an exception.

*Explicitly* handling missing/optional values in INFO or FORMAT fields can be done by checking for NA, e.g.: `INFO["DP"] is NA`.

Handling missing/optional values in fields other than INFO or FORMAT can be done by checking for None, e.g `ID is not None`.

Sometimes, multi-valued fields may contain missing values; in this case, the `without_na` function can be convenient, for example: `mean(without_na(FORMAT['DP'][s] for s in SAMPLES)) > 2.3`. It is also possible to replace `NA` with some constant value with the `replace_na` function: `mean(replace_na((FORMAT['DP'][s] for s in SAMPLES), 0.0)) > 2.3`

### Auxiliary files
`vembrane` supports additional files, such as lists of genes or ids with the `--aux NAME path/to/file` option. The file should contain one item per line and is parsed as a set. For example `vembrane filter --aux genes genes.txt "ANN['SYMBOL'] in AUX['genes']" variants.vcf` will keep only records where the annotated symbol is in the set specified in `genes.txt`.

## `vembrane table`

In addition to the `filter` subcommand, vembrane (`≥ 0.5`) also supports writing tabular data with the `table` subcommand.
In this case, an expression which evaluates to `tuple` is expected, for example:
```sh
vembrane table 'CHROM, POS, 10**(-QUAL/10), ANN["CLIN_SIG"]' input.vcf > table.tsv
```

When handling **multi-sample VCFs**, you often want to iterate over all samples in a record by looking at a `FORMAT` field for all of them.
However, if you use a standard Python list comprehension (something like `[FORMAT['DP'][sample] for sample in SAMPLES]`), this would yield a single column with a list containing one entry per sample (something like `[25, 32, 22]` for three samples with the respective depths).

In order to have a separate column for each sample, you can use the **`for_each_sample()`** function in both the main `vembrane table` expression and the `--header` expression.
It should contain one [lambda expression](https://docs.python.org/3/reference/expressions.html#lambda) with exactly one argument, which will be substituted by the sample names in the lambda expression.

For example, you could specifiy expressions for the `--header` and the main VCF record evaluation like this:
```sh
vembrane table --header 'CHROM, POS, for_each_sample(lambda sample: f"{sample}_depth")' 'CHROM, POS, for_each_sample(lambda s: FORMAT["DP"][s])' input.vcf > table.tsv
```
Given a VCF file with samples `Sample_1`, `Sample_2` and `Sample_3`, the header would expand to be printed as:
```
CHROM  POS   Sample_1_depth   Sample_2_depth   Sample_3_depth
```
and the expression to evaluate on each VCF record would become:
```python
(CHROM, POS, FORMAT['DP']['Sample_1'], FORMAT['DP']['Sample_2'], FORMAT['DP']['Sample_3'])
```

When not supplying a `--header` expression, the entries of the expanded main expression become the column names in the header.
When supplying a header via `--header`,  its `for_each_sample()` expects an expression which can be evaluated to `str` and must have the same number of fields as the main expression.

Please note that, as anywhere in vembrane, you can use arbitrary Python expressions in `for_each_sample()` lambda expressions.
So you can for example perform computations on fields or combine multiple fields into one value:
```sh
vembrane table 'CHROM, POS, for_each_sample(lambda sample: FORMAT["AD"][sample] / FORMAT["DP"][sample] * QUAL)' input.vcf > table.tsv
```

## `vembrane annotate`

vembrane is able to annotate vcf files with a given table-like file. In addition to the vcf and annotation file, the user has to provide a configuration file.

Configuration (Example):

```yaml
## example.yaml
annotation:
    file: "example.tsv" # the table-like annotation file column with header
    columns:
      chrom: "chrom" # column name of the annotation file refering to the chromosome
      start: "chromStart" # column name of the annotation file refering to the chromosome start
      stop: "chromEnd" # column name of the annotation file refering to the chromosome end
    delimiter: "\t" # delimiter of the columns
    values:
    - value: # a new annotation entry in the info field of the vcf
        vcf_name: "genehancer_score" # the name of annotation entry
        number: "1" # number of values for each entry
        description: "Score from genehancer." # description of this entry in the header
        type: "Float" # type of the values
        expression: "DATA['score'][0]" # any python expression to calculate the value(s)
                                       # DATA['score'] refers to the 'score' column of the annotation field
    - value: # a second annotation entry to annotate
        vcf_name: "genehancer_score2"
        number: "1"
        description: "Score from genehancer."
        type: "Float"
        expression: "log(max(DATA['score']) * 2)"
```

example.tsv (Example):
```
chrom	chromStart	chromEnd	name	score
chr10	76001	77000	HJSDHKD	463
chr10	120054	130024	HJSJHKD	463
chr10	432627	492679	IDASJLD	327
chr10	540227	872071	SZAGHSD	435
chr10	654480	1000200	HSJKJSD	12
```

Exemplary invocation: `vembrane annotate example.yaml example.bcf > annotated.vcf`.

Internally for each vcf record the overlapping regions of the annotation file are determined and stored in `DATA`. The expression may then access the `DATA` object and its columns by the columns names to generate a single or multiple values of cardinality `number` of type `type`. These values are stored in the new annotation entry under the name `vcf_name` and with header description `description`.

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
* David Lähnemann (@dlaehnemann)
* Felix Mölder (@felixmoelder)
* Christopher Schröder (@christopher-schroeder)

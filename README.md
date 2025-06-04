[![CI](https://github.com/vembrane/vembrane/actions/workflows/main.yml/badge.svg)](https://github.com/vembrane/vembrane/actions/workflows/main.yml)
[![Zenodo DOI](https://zenodo.org/badge/276383670.svg)](https://zenodo.org/badge/latestdoi/276383670)
[![Paper DOI:10.1093/bioinformatics/btac810](http://img.shields.io/badge/DOI-10.1093/bioinformatics/btac810-3c799f.svg)](https://doi.org/10.1093/bioinformatics/btac810)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/vembrane/README.html)

# vembrane: variant filtering using python expressions

vembrane allows to simultaneously filter variants based on any `INFO` or `FORMAT` field, `CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`, and the annotation field `ANN`. When filtering based on `ANN`, annotation entries are filtered first. If no annotation entry remains, the entire variant is deleted.

vembrane relies on [pysam](https://pysam.readthedocs.io/en/latest/) for reading/writing VCF/BCF files.

For a comparison with similar tools have a look at the [vembrane benchmarks](https://github.com/vembrane/vembrane-benchmark).

## Installation
vembrane is available in [bioconda](https://bioconda.github.io/) and can either be installed into an existing conda environment with `mamba install -c conda-forge -c bioconda vembrane` or into a new named environment `mamba create -n environment_name -c conda-forge -c bioconda vembrane`.
Alternatively, if you are familiar with git and [uv](https://docs.astral.sh/uv/), clone this repository and run `uv sync`.
See [docs/develop.md](docs/develop.md) for further details.

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
  --aux NAME=PATH, -a NAME=PATH
                        Path to an auxiliary file containing a set of symbols.
  --ontology PATH       Path to an ontology in OBO format. 
                        The ontology is loaded into memory and can be used in expressions via the SO symbol.
                        May be compressed with gzip, bzip2 or xz.
                        Defaults to built-in ontology (from sequenceontology.org).
  --keep-unmatched      Keep all annotations of a variant if at least one of them
                        passes the expression (mimics SnpSift behaviour).
  --preserve-order      Ensures that the order of the output matches that of the input.
                        This is only useful if the input contains breakends (BNDs)
                        since the order of all other variants is preserved anyway.
```


### Filter expression
The filter expression can be any valid python expression that evaluates to a value of *type* `bool`.
If you want to use truthy values, you need to wrap the expression in `bool()`, or aggregate multiple values via `any()` or `all()`.

However, functions and symbols available have been restricted to the following:

 * `all`, `any`
 * `abs`, `len`, `max`, `min`, `round`, `sum`
 * `enumerate`, `filter`, `iter`, `map`, `next`, `range`, `reversed`, `sorted`, `zip`
 * `dict`, `list`, `set`, `tuple`
 * `bool`, `chr`, `float`, `int`, `ord`, `str`
 * Any function or symbol from [`math`](https://docs.python.org/3/library/math.html)
 * Any function from [`statistics`](https://docs.python.org/3/library/statistics.html)
 * Regular expressions via [`re`](https://docs.python.org/3/library/re.html)
 * custom functions:
   * `without_na(values: Iterable[T]) -> Iterable[T]` (keep only values that are not `NA`)
   * `replace_na(values: Iterable[T], replacement: T) -> Iterable[T]` (replace values that are `NA` with some other fixed value)
   * genotype related:
     * `count_hom`, `count_het` , `count_any_ref`, `count_any_var`, `count_hom_ref`, `count_hom_var`
     * `is_hom`, `is_het`, `is_hom_ref` , `is_hom_var`
     * `has_ref`, `has_var`

### Available fields
The following VCF fields can be accessed in the filter expression:

| Name      | Type                         | Interpretation                                                                                     | Example expression             |
| --------- | ---------------------------- | -------------------------------------------------------------------------------------------------- | ------------------------------ |
| `INFO`    | `Dict[str, Any¹]`            | `INFO field -> Value`                                                                              | `INFO["DP"] > 0`               |
| `ANN`     | `Dict[str, Any²]`            | `ANN field -> Value`                                                                               | `ANN["Gene_Name"] == "CDH2"`   |
| `CHROM`   | `str`                        | Chromosome Name                                                                                    | `CHROM == "chr2"`              |
| `POS`     | `int`                        | Chromosomal position (1-based)                                                                     | `24 < POS < 42`                |
| `END`     | `int`                        | Chromosomal end position (1-based, inclusive, NA for breakends); also accessible via `INFO["END"]` | `24 < END < 42`                |
| `ID`      | `str`                        | Variant ID                                                                                         | `ID == "rs11725853"`           |
| `REF`     | `str`                        | Reference allele                                                                                   | `REF == "A"`                   |
| `ALT`     | `str`                        | Alternative allele³                                                                                | `ALT == "C"`                   |
| `QUAL`    | `float`                      | Quality                                                                                            | `QUAL >= 60`                   |
| `FILTER`  | `List[str]`                  | Filter tags                                                                                        | `"PASS" in FILTER`             |
| `FORMAT`  | `Dict[str, Dict[str, Any¹]]` | `Format -> (Sample -> Value)`                                                                      | `FORMAT["DP"][SAMPLES[0]] > 0` |
| `SAMPLES` | `List[str]`                  | `[Sample]`                                                                                         | `"Tumor" in SAMPLES`           |
| `INDEX`   | `int`                        | Index of variant in the file                                                                       | `INDEX < 10`                   |

 ¹ depends on type specified in VCF header

 ² for the usual snpeff and vep annotations, custom types have been specified; any unknown ANN field will simply be of type `str`. If something lacks a custom parser/type, please consider filing an issue in the [issue tracker](https://github.com/vembrane/vembrane/issues).

 ³ vembrane does not handle multi-allelic records itself. Instead, such files should be
 preprocessed by either of the following tools (preferably even before annotation):
 - [`bcftools norm -m-any […]`](http://samtools.github.io/bcftools/bcftools.html#norm)
 - [`gatk LeftAlignAndTrimVariants […] --split-multi-allelics`](https://gatk.broadinstitute.org/hc/en-us/articles/360037225872-LeftAlignAndTrimVariants)
 - [`vcfmulti2oneallele […]`](http://lindenb.github.io/jvarkit/VcfMultiToOneAllele.html)


### Examples

* Only keep annotations and variants where gene equals "CDH2" and its impact is "HIGH":
  ```sh
  vembrane filter 'ANN["Gene_Name"] == "CDH2" and ANN["Annotation_Impact"] == "HIGH"' variants.bcf
  ```
* Only keep variants with quality at least 30:
  ```sh
  vembrane filter 'QUAL >= 30' variants.vcf
  ```
* Only keep annotations and variants where feature (transcript) is ENST00000307301:
  ```sh
  vembrane filter 'ANN["Feature"] == "ENST00000307301"' variants.bcf
  ```
* Only keep annotations and variants where protein position is less than 10:
  ```sh
  vembrane filter 'ANN["Protein_position"].start < 10' variants.bcf
  ```
* Only keep variants where the ID matches the regex pattern `^rs[0-9]+`:
  ```sh
  vembrane filter 'bool(re.search("^rs[0-9]+", ID or ""))' variants.vcf
  ```
* Only keep variants where mapping quality is exactly 60:
  ```sh
  vembrane filter 'INFO["MQ"] == 60' variants.bcf
  ```
* Only keep annotations and variants where CLIN_SIG contains "pathogenic", "likely_pathogenic" or "drug_response":
  ```sh
  vembrane filter \
    'any(entry in ANN["CLIN_SIG"]
         for entry in ("pathogenic", "likely_pathogenic", "drug_response"))' \
    variants.vcf
  ```
  Using set operations, the same may also be expressed as:
  ```sh
  vembrane filter \
    'not {"pathogenic", "likely_pathogenic", "drug_response"}.isdisjoint(ANN["CLIN_SIG"])' \
    variants.vcf
  ```
* Filter on sample specific values:
  * by sample name:
    ```sh
    vembrane filter 'FORMAT["DP"]["specific_sample_name"] > 0' variants.vcf
    ```
  * by sample index:
    ```sh
    vembrane filter 'FORMAT["DP"][0] > 0' variants.vcf
    ```
  * by sample name based on the index in the list of `SAMPLES`:
    ```sh
    vembrane filter 'FORMAT["DP"][SAMPLES[0]] > 0' variants.vcf
    ```
  * using all or a subset of `SAMPLES`:
      ```sh
      vembrane filter 'mean(FORMAT["DP"][s] for s in SAMPLES) > 10' variants.vcf
      ```

* Filter on genotypes for specific samples (named "kid", "mom", "dad"):
  ```sh
  vembrane filter \
    'is_het("kid") and is_hom_ref("mom") and is_hom_ref("dad") and \
     all(FORMAT["DP"][s] > 10 for s in ["kid", "mom", "dad"])' \
    variants.vcf
  ```
* Explicitly access the `GT` field for the first sample in the file:
  ```sh
  vembrane filter 'FORMAT["GT"][0] == (1, 1)' variants.vcf
  ```

### Custom `ANN` types
`vembrane` parses entries in the annotation field as outlined in [docs/ann_types.md](docs/ann_types.md).

### Missing values in annotations

If a certain annotation field lacks a value, it will be replaced with the special value of `NA`. Comparing with this value will always result in `False`, e.g.
`ANN["MOTIF_POS"] > 0` will always evaluate to `False` *if* there was no value in the "MOTIF_POS" field of ANN (otherwise the comparison will be carried out with the usual semantics).

For fields with custom types, such as `ANN["Protein_position"]` which is of type `PosRange` with attributes `start`, `end` and `length`, trying to access `ANN["Protein_position"].start` will result in `NA` if there was no value for `"Protein_position"` in the annotation of the respective record, i.e. the access will return `NA` instead of raising an `AttributeError`.
In general, any attribute access on `NA` will result in `NA` (and issue a warning to stderr).

Since you may want to use the regex module to search for matches, `NA` also acts as an empty `str`, such that `re.search("nothing", NA)` returns nothing instead of raising an exception.

*Explicitly* handling missing/optional values in INFO or FORMAT fields can be done by checking for NA, e.g.: `INFO["DP"] is NA`.

Handling missing/optional values in fields other than INFO or FORMAT can be done by checking for None, e.g `ID is not None`.

Sometimes, multi-valued fields may contain missing values; in this case, the `without_na` function can be convenient, for example: `mean(without_na(FORMAT['DP'][s] for s in SAMPLES)) > 2.3`. It is also possible to replace `NA` with some constant value with the `replace_na` function: `mean(replace_na((FORMAT['DP'][s] for s in SAMPLES), 0.0)) > 2.3`

### Auxiliary files
`vembrane` supports additional files, such as lists of genes or ids with the `--aux NAME=path/to/file` option. The file should contain one item per line and is parsed as a set. For example `vembrane filter --aux genes=genes.txt "ANN['SYMBOL'] in AUX['genes']" variants.vcf` will keep only records where the annotated symbol is in the set specified in `genes.txt`.

### Ontologies
`vembrane` supports ontologies in OBO format. The ontology is loaded into memory and can be accessed in the filter expression via the `SO` symbol. This enables filtering based on relationships between ontology terms. 
For example, `vembrane filter --ontology so.obo 'ANN["Consequence"].any_is_a("intron_variant")'` will keep only records where at least one of the consequences is an intron variant *or a subtype thereof*.
If no ontology is provided, the built-in ontology from sequenceontology.org (date: 2024-06-06) is loaded automatically if the `SO` symbol is accessed.

There are three relevant classes/types:
- `Term`: Represents a term in the ontology. It inherits from `str` and can be used as such.
- `Consequences`: Represents a list of terms. It inherits from `list` and can be used as such.
- `SO`: Represents the ontology itself. It is a singleton and can be used to access the ontology.

The following functions are available for ontologies, where `term` is a single `Term` and `terms` is a `Consequences` object:
- `SO.get_id(term: Term) -> str`: Convert from term name (e.g. `stop_gained`) to accession (e.g. `SO:0001587`).
- `SO.get_term(id_: str) -> Term`: Convert from accession (e.g. `SO:0001587`) to term name (e.g. `stop_gained`).
- `terms.most_specific_terms() -> Consequences`: Narrow down the list of terms to the most specific ones, e.g. `frameshift_variant&frameshift_truncation&intron_variant&splice_site_variant&splice_donor_5th_base_variant` will lead to `frameshift_truncation&splice_donor_5th_base_variant`.
- `term.ancestors() -> Consequences`: Get *all* ancestral levels of a term, all the way to the ontology's root node.
- `term.descendants() -> Consequences`: Get *all* descendant levels of a term, all the way to the ontology's respective leave nodes.
- `term.parents() -> Consequences`: Get immediate parents of a term.
- `term.children() -> Consequences`: Get immediate children of a term.
- `term.is_a(parent: Term) -> bool`: Check if there is a path from `term` to `parent`, i.e. whether `term` is the `parent` type or a subtype of it.
- `terms.any_is_a(parent: Term) -> bool`: Check if any of the terms is a subtype of `parent`.
- `term.is_ancestor(other: Term) -> bool`: Check if `term` is an ancestor of `other`.
- `term.is_descendant(other: Term) -> bool`: Check if `term` is a descendant of `other`. (Same as `is_a`)
- `term.path_length(target: Term) -> int | None`: Get the shortest path length from `term` to `target` *or vice versa*. Returns `None` if no path exists.

## `vembrane tag`
While `vembrane filter` removes/skips records which do not pass the supplied expression,
`vembrane tag` applies tags to records which do pass the expression.
In other words, `tag` is a non-destructive version of `filter`, which only adds tags to records (hence never removes information from the input).
As such, multiple tags can be applied to a single record in the same run.

Note that the VCF specification explicitly defines the `FILTER` field as a "list of codes for filters that *fail*" (emphasis by us).
As such, `PASS` *should* indicate that all filters have passed, and `.` (missing) *should* indicate that no filters have been applied.
Custom tags should therefore indicate whether a record failed a certain filter, and not whether it passed.

However, the default behaviour is to tag records which pass the expression, and not those which fail.
To toggle between applying tags in a positive (tagging records that *pass*) or negative (tagging records that *fail*) sense, use the `--tag-mode [pass|fail]`/`-m [pass|fail]` flag.

Tag names must neither contain whitespace nor semicolons. Additionally, the tag name must also not be `"0"`.

### Examples
* Tag records with quality at least 30 with a tag named `quality_at_least_30`:
  ```sh
  vembrane tag --tag quality_at_least_30="QUAL >= 30" variants.vcf
  ```
* Tag records with quality at least 30 with a tag named `quality_at_least_30`
  and records which have "Illumina" in their list of platforms with a tag named `illumina`:
  ```sh
  vembrane tag -t quality_at_least_30="QUAL >= 30" -t illumina "'Illumina' in INFO['platformnames']" variants.vcf
  ```
* Tag records with quality *less than* 30 with a tag named `q_below_30` using the `--tag-mode fail` setting:
    ```sh
    vembrane tag --tag-mode fail --tag q_below_30="QUAL >= 30" variants.vcf
    ```
* Tag records with quality *less than* 30 with a tag named `q_below_30` by negating the expression itself:
    * ```sh
      vembrane tag --tag q_below_30="not (QUAL >= 30)" variants.vcf
      ```
    * ```sh
      vembrane tag --tag q_below_30="QUAL < 30" variants.vcf
      ```


## `vembrane table`

In addition to the `filter` subcommand, vembrane (`≥ 0.5`) also supports writing tabular data with the `table` subcommand.
In this case, an expression which evaluates to `tuple` is expected, for example:
```sh
vembrane table 'CHROM, POS, 10**(-QUAL/10), ANN["CLIN_SIG"]' input.vcf > table.tsv
```

When handling **multi-sample VCFs**, you often want to iterate over all samples in a record by looking at a `FORMAT` field for all of them.
Therefore, `vembrane table` defaults to a long table format:
In this case, the first column will always be called `SAMPLE` and there's an additional variable of the same name available for the expressions.
For example:
```sh
vembrane table 'CHROM, POS, FORMAT["AD"][SAMPLE] / FORMAT["DP"][SAMPLE] * QUAL' input.vcf > long_table.tsv
```
will yield a table with the columns `'SAMPLE'`, `'CHROM'`, `'POS'`, and `'FORMAT["AD"][SAMPLE] / FORMAT["DP"][SAMPLE] * QUAL'`.

If you instead want a wide table format, where each sample has its own column, you can toggle this behaviour with the `--wide` flag:
```sh
vembrane table --wide --header 'CHROM, POS, for_each_sample(lambda sample: f"{sample}_depth")' 'CHROM, POS, for_each_sample(lambda s: FORMAT["DP"][s])' input.vcf > table.tsv
```

This makes use of the **`for_each_sample()`** function in both the main `vembrane table` expression and the `--header` expression.
It should contain one [lambda expression](https://docs.python.org/3/reference/expressions.html#lambda) with exactly one argument, which will be substituted by the sample names in the lambda expression.

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
vembrane table --wide 'CHROM, POS, for_each_sample(lambda sample: FORMAT["AD"][sample] / FORMAT["DP"][sample] * QUAL)' input.vcf > table.tsv
```



## `vembrane table ALL`
If you want to extract all information from a VCF file, including every single `INFO`, `FORMAT` and annotation `ANN`/`INFO["ANN"]` field that is defined in the header, you can use the `table` subcommand with the pseudo-expression `ALL`:
```sh
vembrane table 'ALL' input.vcf > table.tsv
```
To control the naming convention of the columns, you can use the `--naming-convention` option with the following allowed values:
  - `dictionary`: The column names are rendered as a python dictionary acces, e.g. `INFO["DP"]`.
  - `underscore`: The column names are rendered with underscores, e.g. `INFO_DP`.
  - `slash`: The column names are rendered with slashes, e.g. `INFO/DP` (`bcftools` style).
The default is `dictionary`.

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

## `vembrane structured`

The `structured` subcommand allows you to convert VCF records into structured data formats such as JSON, JSONL, or YAML based on a [YTE template](https://yte-template-engine.github.io).

### Usage
```
usage: vembrane structured [options] template [input vcf]

options:
  -h, --help            show this help message and exit
  --annotation-key FIELDNAME, -k FIELDNAME
                        The INFO key for the annotation field. Defaults to "ANN".
  --output OUTPUT, -o OUTPUT
                        Output file. If not specified, output is written to STDOUT.
  --output-fmt {json,jsonl,yaml}
                        Output format. If not specified, can be automatically determined from the --output file extension.
```

### Examples

* Convert VCF records to JSON format using a YTE template:
  ```sh
  vembrane structured template.yml input.vcf --output output.json
  ```

* Convert VCF records to YAML format and write to STDOUT:
  ```sh
  vembrane structured template.yml input.vcf --output-fmt yaml
  ```

* Convert VCF records to JSONL format and write to a file:
  ```sh
  vembrane structured template.yml input.vcf --output output.jsonl
  ```

In the template file, you can define the desired structure and expressions to retrieve data from the VCF record.
The YTE template thereby models the desired structure into which each VCF record shall be converted.
Inside of the template, VCF record specific variable are accessible analogous to expressions in other vembrane commands, for example:

```yaml
variant:
  chromosome: ?CHROM
  position: ?POS
  reference_allele: ?REF
  alternative_allele: ?ALT
  ?if ANN:
    ?if ANN["GENE"]:
      gene: ?ANN["GENE"]
    impact: ?ANN["IMPACT"]
```

As can be seen, YTE supports the specification of Python expressions for templating.
This works by prefixing strings with `?`.
More YTE details and examples can be found in the [YTE documentation](https://yte-template-engine.github.io).

A more complex example, leveraging most capabilities of YTE, is the following:

```yaml
__variables__:
  samples_with_af: "?[sample for sample in SAMPLES if FORMAT['AF'][sample] is not NA]"

variant:
  chrom: ?CHROM
  pos: ?POS
  ref: ?REF
  alt: ?ALT
  qual: ?QUAL
  ?if ID is not None:
    id: ?ID
  ?if INFO["SVLEN"] is not NA:
    svlen: ?INFO["SVLEN"]
  ?if ANN:
    ?if ANN["SYMBOL"]:
      gene: ?ANN["SYMBOL"]
    impact: ?ANN["IMPACT"]
  ?if samples_with_af:
    samples:
      ?for sample in samples_with_af:
        ?sample:
          allelic_fraction: ?f"{FORMAT['AF'][sample]:.0%}"
```

* We define a variable at the top, collecting all samples having a value in the AF format field.
* If the variant record has a value for ID, this is included in the output.
* If the variant record has a value for INFO/SVLEN, this is included in the output. Note that unlike all the primary optional fields like ID, QUAL etc., missing values in INFO and FORMAT are represented as `NA` instead of `None`.
* If the record has annotation, we show gene symbol (if present) and impact.
* If there is at least one sample with allele frequency (`AF`) information, we show this in a substructure with an entry for each such sample.


## Citation
Check the "Cite this repository" entry in the sidebar for citation options.

Also, please read [should-I-cite-this-software](https://github.com/mr-c/shouldacite/blob/main/should-I-cite-this-software.md) for background.

## Authors

* Marcel Bargull (@mbargull)
* Jan Forster (@jafors)
* Till Hartmann (@tedil)
* Johannes Köster (@johanneskoester)
* Elias Kuthe (@eqt)
* David Lähnemann (@dlaehnemann)
* Felix Mölder (@felixmoelder)
* Christopher Schröder (@christopher-schroeder)

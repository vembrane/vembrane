[![CI](https://github.com/vembrane/vembrane/actions/workflows/main.yml/badge.svg)](https://github.com/vembrane/vembrane/actions/workflows/main.yml)
[![Zenodo DOI](https://zenodo.org/badge/276383670.svg)](https://zenodo.org/badge/latestdoi/276383670)
[![Paper DOI:10.1093/bioinformatics/btac810](http://img.shields.io/badge/DOI-10.1093/bioinformatics/btac810-3c799f.svg)](https://doi.org/10.1093/bioinformatics/btac810)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/vembrane/README.html)

![logo dark](./logo_dark.svg#gh-dark-mode-only)
![logo light](./logo.svg#gh-light-mode-only)

# VCF/BCF transformation using python expressions

vembrane allows to transform VCF/BCF files by specifying flexible Python expressions over `INFO` or `FORMAT` field, `CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`, and the annotation field `ANN`/`CSQ`.
Supported transformations include filtering, tagging, annotation, sorting, and conversion to tabular (CSV/TSV), structured formats (JSON/JSONL/YAML), and [FHIR](https://www.hl7.org/fhir).
vembrane relies on [pysam](https://pysam.readthedocs.io/en/latest/) or [cyvcf2](https://github.com/brentp/cyvcf2) for reading/writing VCF/BCF files.

See **https://vembrane.github.io** for a high-level overview.
For a comparison with similar tools have a look at the [vembrane benchmarks](https://github.com/vembrane/vembrane-benchmark).

## Installation
vembrane is available in [bioconda](https://bioconda.github.io/) and can either be installed into an existing conda environment with `mamba install -c conda-forge -c bioconda vembrane` or into a new named environment `mamba create -n environment_name -c conda-forge -c bioconda vembrane`.
Alternatively, if you are familiar with git and [uv](https://docs.astral.sh/uv/), clone this repository and run `uv sync`.
See [docs/develop.md](docs/develop.md) for further details.

## Subcommands

vembrane provides several subcommands for different tasks:

*   `filter`: Filters VCF/BCF files based on flexible Python expressions. See the detailed documentation below.
    ```sh
     vembrane filter 'CHROM == "chr3" and ANN["Consequence"].any_is_a("frameshift_variant")' variants.bcf
     ```

*   `tag`: A non-destructive version of `filter`. Instead of removing records, it adds a user-defined tag to the `FILTER` field of records that match an expression. For more details, see [`docs/tag.md`](docs/tag.md).
    ```sh
    vembrane tag --tag quality_at_least_30="QUAL >= 30" variants.vcf
    ```

*   `table`: Creates tabular (TSV) files from VCF/BCF data. The columns can be defined with flexible Python expressions. For more details, see [`docs/table.md`](docs/table.md).
    ```sh
    vembrane table 'CHROM, POS, 10**(-QUAL/10), ANN["CLIN_SIG"]' input.vcf > table.tsv
    ```

*   `annotate`: Annotates VCF/BCF files with data from an external table-like file (e.g., TSV, CSV), based on genomic coordinates. For more details, see [`docs/annotate.md`](docs/annotate.md).
    ```sh
    vembrane annotate example.yaml example.bcf > annotated.vcf
    ```

*   `sort`: Sorts VCF/BCF files. The sort order can be defined via flexible Python expressions and thus be based on any field and annotation that occurs in the VCF/BCF (e.g. impact or clinical significance). The primary use case is variant prioritization. For more details, see [`docs/sort.md`](docs/sort.md).
    ```sh
    vembrane sort input.vcf 'round(ANN["gnomad_AF"], 1), -ANN["REVEL"]' > prioritized.vcf
    ```

*   `structured`: Converts VCF records into structured formats like JSON, JSONL, or YAML using a flexible YTE template. For more details, see [`docs/structured.md`](docs/structured.md).
    ```sh
    vembrane structured template.yml input.vcf --output output.json
    ```

*   `fhir`: Converts VCF records into FHIR observations. For more details, see [`docs/fhir.md`](docs/fhir.md).
    ```sh
    vembrane fhir tumor GRCh38 --profile mii_molgen_v2025.0.0 --output-fmt json --annotation-key ANN < sample.vcf > sample-tumor.fhir.json
    ```

## `vembrane filter`
Full documentation is also available in [`docs/filter.md`](docs/filter.md).

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
                        The INFO key for the annotation field. This defaults to "ANN", but tools might
                        use other field names. For example, default VEP annotations can be parsed by
                        setting "CSQ" here.
  --aux NAME=PATH, -a NAME=PATH
                        Path to an auxiliary file containing a set of symbols.
  --context CONTEXT     Python statement defining a context for given Python expressions.
                        Extends eventual definitions given via --context-file.
  --context-file CONTEXT_FILE
                        Path to Python script defining a context for given Python expressions.
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
| `ANN`²    | `Dict[str, Any³]`            | `ANN field -> Value`²                                                                              | `ANN["SYMBOL"] == "CDH2"`²     |
| `CHROM`   | `str`                        | Chromosome Name                                                                                    | `CHROM == "chr2"`              |
| `POS`     | `int`                        | Chromosomal position (1-based)                                                                     | `24 < POS < 42`                |
| `END`     | `int`                        | Chromosomal end position (1-based, inclusive, NA for breakends); also accessible via `INFO["END"]` | `24 < END < 42`                |
| `ID`      | `str`                        | Variant ID                                                                                         | `ID == "rs11725853"`           |
| `REF`     | `str`                        | Reference allele                                                                                   | `REF == "A"`                   |
| `ALT`     | `str`                        | Alternative allele⁴                                                                                | `ALT == "C"`                   |
| `QUAL`    | `float`                      | Quality                                                                                            | `QUAL >= 60`                   |
| `FILTER`  | `List[str]`                  | Filter tags                                                                                        | `"PASS" in FILTER`             |
| `FORMAT`  | `Dict[str, Dict[str, Any¹]]` | `Format -> (Sample -> Value)`                                                                      | `FORMAT["DP"][SAMPLES[0]] > 0` |
| `SAMPLES` | `List[str]`                  | `[Sample]`                                                                                         | `"Tumor" in SAMPLES`           |
| `INDEX`   | `int`                        | Index of variant in the file                                                                       | `INDEX < 10`                   |

 ¹ depends on type specified in VCF header

 ² if your VCF defines annotations under a key other than `ANN` (e.g. VEP's `CSQ`) you have to specify this via the `--annotation-key` flag (e.g. `--annotation-key CSQ`). You can (and should, for portability) still use `ANN` in your expressions then (although the given annotation key works as well, e.g. `CSQ["SYMBOL"]`).

 ³ for the usual snpeff and vep annotations, custom types have been specified; any unknown ANN field will simply be of type `str`. If something lacks a custom parser/type, please consider filing an issue in the [issue tracker](https://github.com/vembrane/vembrane/issues).

 ⁴ vembrane does not handle multi-allelic records itself. Instead, such files should be
 preprocessed by either of the following tools (preferably even before annotation):
 - [`bcftools norm -m-any […]`](http://samtools.github.io/bcftools/bcftools.html#norm)
 - [`gatk LeftAlignAndTrimVariants […] --split-multi-allelics`](https://gatk.broadinstitute.org/hc/en-us/articles/360037225872-LeftAlignAndTrimVariants)
 - [`vcfmulti2oneallele […]`](http://lindenb.github.io/jvarkit/VcfMultiToOneAllele.html)


### Examples

* Only keep annotations and variants where gene equals "CDH2" and its impact is "HIGH":
  ```sh
  vembrane filter 'ANN["SYMBOL"] == "CDH2" and ANN["Annotation_Impact"] == "HIGH"' variants.bcf
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

`vembrane` parses entries in the annotation field (`ANN` or whatever you specify under `--annotation-key`) as outlined in [docs/ann_types.md](docs/ann_types.md).

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

### Custom context
In addition to the default context for expressions given to vembrane, consisting of variables containing the VCF/BCF record information, helper functions and useful Python builtins (like `any()`), it is possible to define additional custom context via the flags

* `--context STATEMENT`: providing a Python statement,
* `--context-file PATH`: providing a Python script to evaluate (considering the parent directory of the script as search path for eventual additional imports).

In both cases, the provided input is evaluated with Python. Any global variables (or functions) become available in the Python expressions that you provide to vembrane for filtering or the other subcommands. Note that the code you provide as context is not sandboxed and should be trusted. Carefully review any code you get from the internet or AI.

An example use case is the import of additional Python modules. For example, you could use this to randomly subsample 10% of the VCF/BCF records given in a file `input.vcf`:

```bash
vembrane filter --context "import random; random.seed(42)" "random.random() < 0.1" input.vcf > subsampled.vcf
```

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

# `vembrane sort`

The `sort` subcommand allows to sort VCF/BCF files via keys defined as Python expressions (in ascending order).
The Python expressions are analogous to expressions in other vembrane commands.
This feature loads the entire VCF file into memory in order to maximize performance.
It is thus meant to sort small, already filtered VCF files, e.g. for **prioritizing records for the human eye**.
For large VCF files, the only relevant sorting is usually by position, which is better done with e.g. bcftools (and usually the sorting that variant callers output).

### Usage

```
usage: vembrane sort [-h] [--output OUTPUT] [--output-fmt {vcf,bcf,uncompressed-bcf}] [--preserve-annotation-order]
                     [--chunk-size CHUNK_SIZE] [--annotation-key FIELDNAME] [--aux NAME=PATH] [--ontology PATH]
                     [--overwrite-number-info FIELD=NUMBER] [--overwrite-number-format FIELD=NUMBER]
                     [--backend {cyvcf2,pysam}]
                     [vcf] expression

Sort VCF records by one or multiple Python expressions that encode keys for the desired order. This feature loads the entire
VCF file into memory in order to maximize performance. It is thus meant to sort small, already filtered VCF files, e.g. for
prioritizing records for the human eye. For large VCF files, the only relevant sorting is usually by position, which is better
done with e.g. bcftools (and usually the sorting that variant callers output).

positional arguments:
  vcf                   The VCF/BCF file containing the variants. If not specified, reads from STDIN.
  expression            Python expression (or tuple of expressions) returning orderable values to sort the VCF records by
                        (ascending, smallest values coming first). If multiple expressions are provided as a tuple, they are
                        prioritized from left to right with lowest priority on the right. NA/NaN values are sorted to the end.

options:
  -h, --help            show this help message and exit
  --output OUTPUT, -o OUTPUT
                        Output file, if not specified, output is written to STDOUT.
  --output-fmt {vcf,bcf,uncompressed-bcf}, -O {vcf,bcf,uncompressed-bcf}
                        Output format.
  --preserve-annotation-order
                        If set, annotations are not sorted within the records, but kept in the same order as in the input VCF
                        file. If not set (default), annotations are sorted within the record according to the given keys if
                        any of the sort keys given in the python expression refers to an annotation.
  --chunk-size CHUNK_SIZE
                        Number of VCF records to sort in memory. If the VCF file exceeds this number of records, external
                        sorting is used.
  --annotation-key FIELDNAME, -k FIELDNAME
                        The INFO key for the annotation field. This defaults to 'ANN', but tools might use other field names.
                        For example, default VEP annotations can be parsed by setting 'CSQ' here.
  --aux NAME=PATH, -a NAME=PATH
                        Path to an auxiliary file containing a set of symbols
  --ontology PATH       Path to an ontology in OBO format. May be compressed with gzip, bzip2 and xz. Defaults to built-in
                        ontology (from sequenceontology.org).
  --overwrite-number-info FIELD=NUMBER
                        Overwrite the number specification for INFO fields given in the VCF header. Example: `--overwrite-
                        number cosmic_CNT=.`
  --overwrite-number-format FIELD=NUMBER
                        Overwrite the number specification for FORMAT fields given in the VCF header. Example: `--overwrite-
                        number-format DP=2`
  --backend {cyvcf2,pysam}, -b {cyvcf2,pysam}
                        Set the backend library.
```

### Example

The following command sorts records first by `gnomad_AF` (binned and ascending), and then by `REVEL` score (descending).
The descending sort is achieved by negating the `REVEL` value (`-ANN['REVEL']`) in the key expression.

```bash
vembrane sort input.vcf 'round(ANN["gnomad_AF"], 1), -ANN["REVEL"]' > prioritized.vcf
```
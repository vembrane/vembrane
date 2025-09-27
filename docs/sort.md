# `vembrane sort`

The `sort` subcommand allows to sort VCF/BCF files via keys defined as Python expressions (in ascending order).
The Python expressions are analogous to expressions in other vembrane commands.
This feature is primarily meant to sort 
already filtered VCF files, e.g. for **prioritizing records for the human eye**. 
For large VCF files, the only relevant sorting is usually by position, 
which is better done with e.g. bcftools (and usually the standard sorting 
that variant callers output).

### Usage

```
usage: vembrane sort [-h] [--output OUTPUT] [--output-fmt {vcf,bcf,uncompressed-bcf}] [--preserve-annotation-order]
                     [--max-in-mem-records MAX_IN_MEM_RECORDS] [--annotation-key FIELDNAME] [--aux NAME=PATH]
                     [--ontology PATH] [--overwrite-number-info FIELD=NUMBER] [--overwrite-number-format FIELD=NUMBER]
                     [--backend {cyvcf2,pysam}]
                     expression [vcf]

positional arguments:
  expression            Python expression (or tuple of expressions) returning orderable values (keys) to sort the VCF records
                        by. By default keys are considered in ascending order. To sort by descending order, use
                        `desc(<expression>)` on the entire expression or on individual items of the tuple. If multiple
                        expressions are provided as a tuple, they are prioritized from left to right with lowest priority on
                        the right. NA/NaN values are always sorted to the end. Expressions on annotation entries will cause
                        the annotation with the minimum key value (or maximum if descending) to be considered to sort the
                        record.
  vcf                   The VCF/BCF file containing the variants. If not specified, reads from STDIN. (default: -)

options:
  -h, --help            show this help message and exit
  --output OUTPUT, -o OUTPUT
                        Output file, if not specified, output is written to STDOUT. (default: -)
  --output-fmt {vcf,bcf,uncompressed-bcf}, -O {vcf,bcf,uncompressed-bcf}
                        Output format. (default: vcf)
  --preserve-annotation-order
                        If set, annotations are not sorted within the records, but kept in the same order as in the input VCF
                        file. If not set (default), annotations are sorted within the record according to the given keys if
                        any of the sort keys given in the python expression refers to an annotation.
  --max-in-mem-records MAX_IN_MEM_RECORDS
                        Number of VCF records to sort in memory. If the VCF file exceeds this number of records, external
                        sorting is used. (default: 100000)
  --annotation-key FIELDNAME, -k FIELDNAME
                        The INFO key for the annotation field. This defaults to 'ANN', but tools might use other field names.
                        For example, default VEP annotations can be parsed by setting 'CSQ' here. (default: ANN)
  --aux NAME=PATH, -a NAME=PATH
                        Path to an auxiliary file containing a set of symbols
  --context CONTEXT     Python statement defining a context for given Python expressions. Extends eventual definitions given via --context-file. Any global variables (or
                        functions) become available in the Python expressions.
                        Note that the code you pass here is not sandboxed and should
                        be trusted. Carefully review any code you get from the internet or AI.
  --context-file CONTEXT_FILE
                        Path to Python script defining a context for given Python expressions. Any global variables (or functions) become available in the Python expressions.
                        Note that the code you pass here is not sandboxed and should
                        be trusted. Carefully review any code you get from the internet or AI.
  --ontology PATH       Path to an ontology in OBO format. May be compressed with gzip, bzip2 and xz. Defaults to built-in
                        ontology (from sequenceontology.org).
  --overwrite-number-info FIELD=NUMBER
                        Overwrite the number specification for INFO fields given in the VCF header. Example: `--overwrite-
                        number cosmic_CNT=.`
  --overwrite-number-format FIELD=NUMBER
                        Overwrite the number specification for FORMAT fields given in the VCF header. Example: `--overwrite-
                        number-format DP=2`
  --backend {cyvcf2,pysam}, -b {cyvcf2,pysam}
                        Set the backend library. (default: cyvcf2)
```

### Examples

The following command sorts records first by `gnomad_AF` (binned into orders of magnitude and ascending), and then by `REVEL` score (descending).
The descending sort is achieved by marking the `REVEL` value  as descending via `desc()` in the key expression.

```bash
vembrane sort 'int(log10(ANN["gnomad_AF"])), desc(ANN["REVEL"])' input.vcf > prioritized.vcf
```

In case of non-numeric values, an order can be defined ad-hoc via an inline dictionary.
For example, in order to get variants with high impact in one of their annotations first, we can define the following.

```bash
vembrane sort '{"HIGH": 0, "MODERATE: 1, "LOW": 2, "MODIFIER" 3}[ANN["IMPACT"]]' input.vcf > prioritized.vcf
```
Since ascending sort is the default, variants with at least one `HIGH` in their annotations will come first.
Moreover, `vembrane` will also sort the annotation entries in the corresponding order, with higher impacts coming first.
This behavior can be disabled via the flag `--preserve-annotation-order`.
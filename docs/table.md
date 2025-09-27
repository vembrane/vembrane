# `vembrane table`


In addition to the `filter` subcommand, vembrane (`≥ 0.5`) also supports writing tabular data with the `table` subcommand.

### Usage

```
usage: vembrane table [-h] [--separator CHAR] [--header TEXT] [--naming-convention CONVENTION] [--wide] [--long LONG] [--output OUTPUT] [--annotation-key FIELDNAME]
                      [--aux NAME=PATH] [--context CONTEXT] [--context-file CONTEXT_FILE] [--ontology PATH] [--overwrite-number-info FIELD=NUMBER]
                      [--overwrite-number-format FIELD=NUMBER] [--backend {cyvcf2,pysam}]
                      expression [vcf]

Convert VCF/BCF records to tabular format.

positional arguments:
  expression            A comma-separated tuple of expressions that define the table column contents. Use ALL to output all fields.
  vcf                   Path to the VCF/BCF file to be filtered. Defaults to '-' for stdin. (default: -)

options:
  -h, --help            show this help message and exit
  --separator CHAR, -s CHAR
                        Define the field separator (default: \t). (default: )
  --header TEXT         Override the automatically generated header. Provide "auto" (default) to automatically generate the header from the expression. Provide a comma
                        separated string to manually set the header. Provide "none" to disable any header output. (default: auto)
  --naming-convention CONVENTION
                        The naming convention to use for column names when generating the header for the ALL expression. (default: dictionary)
  --wide                Instead of using long format with a special SAMPLE column, generate multiple columns per sample with the `for_each_sample` utility function.
  --long LONG           Long format is now the default. For wide format, use `--wide` instead.
  --output OUTPUT, -o OUTPUT
                        Output file, if not specified, output is written to STDOUT. (default: -)
  --annotation-key FIELDNAME, -k FIELDNAME
                        The INFO key for the annotation field. This defaults to 'ANN', but tools might use other field names. For example, default VEP annotations can be
                        parsed by setting 'CSQ' here. (default: ANN)
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
  --ontology PATH       Path to an ontology in OBO format. May be compressed with gzip, bzip2 and xz. Defaults to built-in ontology (from sequenceontology.org).
  --overwrite-number-info FIELD=NUMBER
                        Overwrite the number specification for INFO fields given in the VCF header. Example: `--overwrite-number cosmic_CNT=.`
  --overwrite-number-format FIELD=NUMBER
                        Overwrite the number specification for FORMAT fields given in the VCF header. Example: `--overwrite-number-format DP=2`
  --backend {cyvcf2,pysam}, -b {cyvcf2,pysam}
                        Set the backend library. (default: cyvcf2)
```

### Examples

An expression which evaluates to `tuple` is expected, for example:
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
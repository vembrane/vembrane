# `vembrane table`

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
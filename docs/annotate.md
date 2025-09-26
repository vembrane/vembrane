# `vembrane annotate`

vembrane is able to annotate vcf files with a given table-like file. In addition to the vcf and annotation file, the user has to provide a configuration file.

### Usage

```
usage: vembrane annotate [-h] [--output OUTPUT] [--output-fmt {vcf,bcf,uncompressed-bcf}] [--annotation-key FIELDNAME] [--aux NAME=PATH] [--context CONTEXT]
                         [--context-file CONTEXT_FILE] [--ontology PATH] [--overwrite-number-info FIELD=NUMBER] [--overwrite-number-format FIELD=NUMBER]
                         [--backend {cyvcf2,pysam}]
                         config [vcf]

Add new INFO field annotations to a VCF/BCF from other data sources, using a configuration file.

positional arguments:
  config                The configuration file.
  vcf                   Path to the VCF/BCF file to be filtered. Defaults to '-' for stdin. (default: -)

options:
  -h, --help            show this help message and exit
  --output OUTPUT, -o OUTPUT
                        Output file, if not specified, output is written to STDOUT. (default: -)
  --output-fmt {vcf,bcf,uncompressed-bcf}, -O {vcf,bcf,uncompressed-bcf}
                        Output format. (default: vcf)
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
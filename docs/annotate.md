# `vembrane annotate`

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
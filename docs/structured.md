# `vembrane structured`

The `structured` subcommand allows you to convert VCF records into structured data formats such as JSON, JSONL, or YAML based on a [YTE template](https://yte-template-engine.github.io).

### Usage
```
usage: vembrane structured [-h] [--output-fmt {json,jsonl,yaml}] [--output OUTPUT] [--annotation-key FIELDNAME] [--aux NAME=PATH] [--context CONTEXT]
                           [--context-file CONTEXT_FILE] [--ontology PATH] [--overwrite-number-info FIELD=NUMBER] [--overwrite-number-format FIELD=NUMBER]
                           [--backend {cyvcf2,pysam}]
                           template [vcf]

Create structured output from a VCF/BCF and a YTE template.

positional arguments:
  template              File containing a YTE template with the desired structure per record and expressions that retrieve data from the VCF/BCF record.
  vcf                   Path to the VCF/BCF file to be filtered. Defaults to '-' for stdin. (default: -)

options:
  -h, --help            show this help message and exit
  --output-fmt {json,jsonl,yaml}
                        Output format. If not specified, can be automatically determined from the --output file extension.
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
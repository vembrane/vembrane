# `vembrane structured`

The `structured` subcommand allows you to convert VCF records into structured data formats such as JSON, JSONL, or YAML based on a [YTE template](https://yte-template-engine.github.io).

### Usage
```
usage: vembrane structured [options] template [input vcf]

options:
  -h, --help            show this help message and exit
  --annotation-key FIELDNAME, -k FIELDNAME
                        The INFO key for the annotation field. This defaults to "ANN", but tools might
                        use other field names. For example, default VEP annotations can be parsed by
                        setting "CSQ" here.
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
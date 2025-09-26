# `vembrane fhir`

The `fhir` subcommand allows to convert VCF records into [HL7 FHIR](https://fhir.org/) observations.

### Usage

The subcommand takes three positional arguments, the VCF file, the sample to be used (must be the name of a sample column in the VCF file), and the genome assembly to consider.
In addition, various optional Python expressions can be provided in order to allow more information to be encoded as FHIR observations (see `options`).
The subcommand expects the VCF to be annotated with [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html), using the options `--vcf_info_field ANN --numbers --symbol --hgvs --hgvsg` (`--numbers`, `--symbol`, and `--hgvs` are part of `--everything` as well).

```
usage: vembrane fhir [-h] [--url URL] [--status STATUS] --profile {mii_molgen_v2025.0.0} [--id-source ID_SOURCE] [--genomic-source-class GENOMIC_SOURCE_CLASS] [--sample-allelic-frequency SAMPLE_ALLELIC_FREQUENCY] [--sample-allelic-read-depth SAMPLE_ALLELIC_READ_DEPTH] [--confidence-status CONFIDENCE_STATUS]
                     [--detection-limit DETECTION_LIMIT] [--output-fmt {json,jsonl,yaml}] [--output OUTPUT] [--annotation-key FIELDNAME] [--aux NAME=PATH] [--ontology PATH] [--overwrite-number-info FIELD=NUMBER] [--overwrite-number-format FIELD=NUMBER] [--backend {cyvcf2,pysam}]
                     [vcf] sample {GRCh37,GRCh38}

positional arguments:
  vcf                   Path to the VCF/BCF file to be filtered. Defaults to '-' for stdin.
  sample                The sample to use for generating FHIR output.
  {GRCh37,GRCh38}       The reference assembly used for read mapping.

options:
  -h, --help            show this help message and exit
  --url URL, -u URL     Generic url used as identifier by FHIR e.g. http://<institute>/<department>/VCF
  --status STATUS, -s STATUS
                        Status of findings. E.g. final, preliminary, ...
  --profile {mii_molgen_v2025.0.0}
                        The FHIR profile to use for generating the output, see https://github.com/vembrane/vembrane/tree/main/vembrane/modules/assets/fhir/profiles for available profiles and the degree of support.
  --id-source ID_SOURCE
                        URL to the source of IDs found in the ID column of the VCF file. IDs are only used if this is given.
  --genomic-source-class GENOMIC_SOURCE_CLASS
                        The genomic source class of the given variants as defined by LOINC: https://loinc.org/48002-0. Either provide the name as a string or a Python expression that evaluates to the name, e.g., for Varlociraptor '"Somatic" if INFO["PROB_SOMATIC"] > 0.95 else ...'.
  --sample-allelic-frequency SAMPLE_ALLELIC_FREQUENCY
                        Python expression calculating the the samples allelic frequencyas percentage. E.g. "FORMAT['AF'][sample][0] * 100"
  --sample-allelic-read-depth SAMPLE_ALLELIC_READ_DEPTH
                        Python expression accessing the the samples allelic read depth.Default is: "FORMAT['AD'][sample][1]"
  --confidence-status CONFIDENCE_STATUS
                        Python expression for calculating the variants confidence status being High, Intermediate or Low. E.g. "'High' if QUAL >= 20 else ('Intermediate' if QUAL >= 10 else 'Low')"
  --detection-limit DETECTION_LIMIT
                        Detection limit / sensitivity of the analysis in percent (e.g. 95).
  --output-fmt {json,jsonl,yaml}
                        Output format. If not specified, can be automatically determined from the --output file extension.
  --output OUTPUT, -o OUTPUT
                        Output file, if not specified, output is written to STDOUT.
  --annotation-key FIELDNAME, -k FIELDNAME
                        The INFO key for the annotation field.
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
                        Set the backend library.
```

### Examples

For example, running

```bash
vembrane fhir tumor GRCh38 --profile mii_molgen_v2025.0.0 --output-fmt json --annotation-key ANN < sample.vcf > sample-tumor.fhir.json
```

would create FHIR observations for all VCF records in `sample.vcf`, assuming the GRCh38 genome assembly and using the [MII molgen FHIR profile 2025.0.0](https://www.medizininformatik-initiative.de/Kerndatensatz/KDS_Molekulargenetischer_Befundbericht_V2025/implementation-guides-ImplementationGuide-2025.x-DE-MIIIGModulMolGenDE-TechnischeImplementierung-Variante-Observation.html).
The annotation key usually has to be set to either `ANN` or `CSQ`, depending on where the annotations are stored in the VCF records (see option `--vcf_info_field` of VEP).
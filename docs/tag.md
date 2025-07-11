# `vembrane tag`
While `vembrane filter` removes/skips records which do not pass the supplied expression,
`vembrane tag` applies tags to records which do pass the expression.
In other words, `tag` is a non-destructive version of `filter`, which only adds tags to records (hence never removes information from the input).
As such, multiple tags can be applied to a single record in the same run.

Note that the VCF specification explicitly defines the `FILTER` field as a "list of codes for filters that *fail*" (emphasis by us).
As such, `PASS` *should* indicate that all filters have passed, and `.` (missing) *should* indicate that no filters have been applied.
Custom tags should therefore indicate whether a record failed a certain filter, and not whether it passed.

However, the default behaviour is to tag records which pass the expression, and not those which fail.
To toggle between applying tags in a positive (tagging records that *pass*) or negative (tagging records that *fail*) sense, use the `--tag-mode [pass|fail]`/`-m [pass|fail]` flag.

Tag names must neither contain whitespace nor semicolons. Additionally, the tag name must also not be `"0"`.

### Examples
* Tag records with quality at least 30 with a tag named `quality_at_least_30`:
  ```sh
  vembrane tag --tag quality_at_least_30="QUAL >= 30" variants.vcf
  ```
* Tag records with quality at least 30 with a tag named `quality_at_least_30`
  and records which have "Illumina" in their list of platforms with a tag named `illumina`:
  ```sh
  vembrane tag -t quality_at_least_30="QUAL >= 30" -t illumina "'Illumina' in INFO['platformnames']" variants.vcf
  ```
* Tag records with quality *less than* 30 with a tag named `q_below_30` using the `--tag-mode fail` setting:
    ```sh
    vembrane tag --tag-mode fail --tag q_below_30="QUAL >= 30" variants.vcf
    ```
* Tag records with quality *less than* 30 with a tag named `q_below_30` by negating the expression itself:
    * ```sh
      vembrane tag --tag q_below_30="not (QUAL >= 30)" variants.vcf
      ```
    * ```sh
      vembrane tag --tag q_below_30="QUAL < 30" variants.vcf
      ```
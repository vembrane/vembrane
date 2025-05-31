# Changelog

## [2.0.0](https://github.com/vembrane/vembrane/compare/v1.0.7...v2.0.0) (2025-05-31)


### ⚠ BREAKING CHANGES

* Sequence ontology support ([#173](https://github.com/vembrane/vembrane/issues/173))
* `table 'ALL'` pseudo-expression for converting all VCF information to long table formats; table now defaults to --long format, and --long has been replaced by --wide ([#189](https://github.com/vembrane/vembrane/issues/189))

### Features

* `table 'ALL'` pseudo-expression for converting all VCF information to long table formats; table now defaults to --long format, and --long has been replaced by --wide ([#189](https://github.com/vembrane/vembrane/issues/189)) ([d93b817](https://github.com/vembrane/vembrane/commit/d93b8170fc89c6ec1c35b9dcb572eea075720727))
* new subcommand for generating structured output from vcf using a yte template ([#184](https://github.com/vembrane/vembrane/issues/184)) ([a07bb85](https://github.com/vembrane/vembrane/commit/a07bb8583898c8e760e56b89f48e15edb498fcab))
* Sequence ontology support ([#173](https://github.com/vembrane/vembrane/issues/173)) ([a307316](https://github.com/vembrane/vembrane/commit/a307316d15de032b30457f0de16121eb704a3edc))


### Bug Fixes

* for long table output, deal with VCF files without SAMPLES by including an empty SAMPLE column ([#191](https://github.com/vembrane/vembrane/issues/191)) ([7cf1e2d](https://github.com/vembrane/vembrane/commit/7cf1e2d21ef450729b2d91fba880488af9fe2c5b))

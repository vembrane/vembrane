# Changelog

## [2.4.0](https://github.com/vembrane/vembrane/compare/v2.3.4...v2.4.0) (2025-09-27)


### Features

* add ability to define a custom context for expressions ([#228](https://github.com/vembrane/vembrane/issues/228)) ([4be7631](https://github.com/vembrane/vembrane/commit/4be7631c02562668d3e09874351df0fa5ae6258c))

## [2.3.4](https://github.com/vembrane/vembrane/compare/v2.3.3...v2.3.4) (2025-09-18)


### Bug Fixes

* vembrane structured failing when VCF has no ANN annotations (issue [#224](https://github.com/vembrane/vembrane/issues/224)) ([#225](https://github.com/vembrane/vembrane/issues/225)) ([73a780d](https://github.com/vembrane/vembrane/commit/73a780daabe805a4887a5a1bc49f76aa98b4d081))


### Documentation

* generalize intro and mention homepage ([#226](https://github.com/vembrane/vembrane/issues/226)) ([170cab2](https://github.com/vembrane/vembrane/commit/170cab26d9140e9c22698e52dc30147de0874102))

## [2.3.3](https://github.com/vembrane/vembrane/compare/v2.3.2...v2.3.3) (2025-09-15)


### Bug Fixes

* specify artifact unpack directory ([#221](https://github.com/vembrane/vembrane/issues/221)) ([dbec559](https://github.com/vembrane/vembrane/commit/dbec55983d667056f38d72ca54cb419bf9e86e3f))

## [2.3.2](https://github.com/vembrane/vembrane/compare/v2.3.1...v2.3.2) (2025-09-15)


### Bug Fixes

* upload correct build artifacts during CI release process ([#219](https://github.com/vembrane/vembrane/issues/219)) ([3d37220](https://github.com/vembrane/vembrane/commit/3d37220e7e3549f5c37f40f5ecfc74724acae336))

## [2.3.1](https://github.com/vembrane/vembrane/compare/v2.3.0...v2.3.1) (2025-09-15)


### Bug Fixes

* improve CLI help (consistently mention BCF, clarify sort behavior) ([#217](https://github.com/vembrane/vembrane/issues/217)) ([b69ff3b](https://github.com/vembrane/vembrane/commit/b69ff3bc681d28a9bad6a7d8344f2187eb477a0b))

## [2.3.0](https://github.com/vembrane/vembrane/compare/v2.2.0...v2.3.0) (2025-09-10)


### Features

* always allow to use ANN in expressions, even if the VCF/BCF uses e.g. CSQ for annotations ([#208](https://github.com/vembrane/vembrane/issues/208)) ([713d993](https://github.com/vembrane/vembrane/commit/713d99310c558c1839a513a49ee5882571738d09))


### Bug Fixes

* more modern error base class constructor syntax ([#212](https://github.com/vembrane/vembrane/issues/212)) ([584c3d4](https://github.com/vembrane/vembrane/commit/584c3d4b90e9a3d77431a021c333f950ac4633aa))


### Documentation

* in example, stratify by order of magnitude when sorting by gnomAD allele frequency ([#209](https://github.com/vembrane/vembrane/issues/209)) ([b08597e](https://github.com/vembrane/vembrane/commit/b08597eb19d05f80f4b0dd100d6dbc3422d195d3))
* mention new logic of always being able to use ANN regardless of the used annotation key in the detailed docs for filter ([#211](https://github.com/vembrane/vembrane/issues/211)) ([8872250](https://github.com/vembrane/vembrane/commit/8872250535ef016ba83d3114324cba6cf2467262))

## [2.2.0](https://github.com/vembrane/vembrane/compare/v2.1.0...v2.2.0) (2025-07-17)


### Features

* new subcommand for sorting VCF files via one or multiple Python expressions ([#197](https://github.com/vembrane/vembrane/issues/197)) ([528b41c](https://github.com/vembrane/vembrane/commit/528b41c54dd9220c65a1978b13559519a7fc4b3c))


### Bug Fixes

* handle SV record for vembrane fhir ([21c543b](https://github.com/vembrane/vembrane/commit/21c543b3350bcb7ff2869a473d2ff92fc4637b20))
* improve error handling ([6eec426](https://github.com/vembrane/vembrane/commit/6eec42692bf57ec5ee4d521f91f07ee56d31b863))
* improve error handling ([9725e06](https://github.com/vembrane/vembrane/commit/9725e06397ffa474d320f833b795029947e9bc76))


### Documentation

* add missing help text to commands ([#203](https://github.com/vembrane/vembrane/issues/203)) ([bac54e9](https://github.com/vembrane/vembrane/commit/bac54e95585b17f151e48e26057ef968960f295e))
* complete vep options for vembrane fhir ([8148a43](https://github.com/vembrane/vembrane/commit/8148a438c1db040e7118fa775981f4ae6e655938))
* separate documentation into smaller parts ([#206](https://github.com/vembrane/vembrane/issues/206)) ([5808c36](https://github.com/vembrane/vembrane/commit/5808c36e50b6fc3d03fd688c09912e2dd0cbe42a))

## [2.1.0](https://github.com/vembrane/vembrane/compare/v2.0.0...v2.1.0) (2025-06-13)


### Features

* subcommand for generating FHIR from VCF ([#187](https://github.com/vembrane/vembrane/issues/187)) ([7ebd37f](https://github.com/vembrane/vembrane/commit/7ebd37f666d7a92c8da733ae0cf2d684bd411318))


### Documentation

* improve explanation of specifying (and using) alternative `--annotation-key`s ([#193](https://github.com/vembrane/vembrane/issues/193)) ([17506ff](https://github.com/vembrane/vembrane/commit/17506ff5d81748edfe59cbc32200a3bab6241d83))

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

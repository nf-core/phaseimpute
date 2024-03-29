# nf-core/phaseimpute: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of nf-core/phaseimpute, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Changed`

- [#18](https://github.com/nf-core/phaseimpute/pull/18)
  - Maps and region by chromosome
  - update tests config files
  - correct meta map propagation
  - Test impute and test sim works
- [#19](https://github.com/nf-core/phaseimpute/pull/19) - Changed reference panel to accept a csv, update modules and subworkflows (glimpse1/2 and shapeit5)
- [#20](https://github.com/nf-core/phaseimpute/pull/20) - Added automatic detection of vcf contigs for the reference panel and automatic renaming available
- [#22](https://github.com/nf-core/phaseimpute/pull/20) - Add validation step for concordance analysis. Input channels changed to
  match inputs steps.

### `Fixed`

- [#15](https://github.com/nf-core/phaseimpute/pull/15) - Changed test csv files to point to nf-core repository
- [#16](https://github.com/nf-core/phaseimpute/pull/16) - Removed outdir from test config files

### `Dependencies`

### `Deprecated`

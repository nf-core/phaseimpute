# nf-core/phaseimpute: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of nf-core/phaseimpute, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- [#20](https://github.com/nf-core/phaseimpute/pull/20) - Added automatic detection of vcf contigs for the reference panel and automatic renaming available
- [#22](https://github.com/nf-core/phaseimpute/pull/20) - Add validation step for concordance analysis. Input channels changed to match inputs steps. Outdir folder organised by steps. Modules config by subworkflows.
- [#26](https://github.com/nf-core/phaseimpute/pull/26) - Added QUILT method
- [#47](https://github.com/nf-core/phaseimpute/pull/47) - Add possibility to remove samples from reference panel. Add glimpse2 chunking method. Add full-size test parameters.
- [#58](https://github.com/nf-core/phaseimpute/pull/58) - Add external params posfile and chunks. Add glimpse2 phasing and imputation.
- [#67](https://github.com/nf-core/phaseimpute/pull/67) - Export CSVs from each step.
- [#71](https://github.com/nf-core/phaseimpute/pull/71) - Allow external panel to be used in step impute.
- [#97](https://github.com/nf-core/phaseimpute/pull/97) - Add dog reference panel and config to test pipeline with other species.
- [#102](https://github.com/nf-core/phaseimpute/pull/102) - Add dog panel test.
- [#119](https://github.com/nf-core/phaseimpute/pull/119) - Add dog test with panelprep and imputation.

### `Changed`

- [#18](https://github.com/nf-core/phaseimpute/pull/18)
  - Maps and region by chromosome
  - update tests config files
  - correct meta map propagation
  - Test impute and test sim works
- [#19](https://github.com/nf-core/phaseimpute/pull/19) - Changed reference panel to accept a csv, update modules and subworkflows (glimpse1/2 and shapeit5)
- [#40](https://github.com/nf-core/phaseimpute/pull/40) - Add STITCH method. Reorganize panelprep subworkflows.
- [#51](https://github.com/nf-core/phaseimpute/pull/51) - Update all process and fix linting errors. Remove fastqc added by the template.
- [#56](https://github.com/nf-core/phaseimpute/pull/56) - Move to nf-test to check the output files names generated. Fix validation and concatenation by chromosomes missing. Add dedicated GLIMPSE1 subworkflow. Fix posfile generation to be done once for glimpse and stitch.
- [#68](https://github.com/nf-core/phaseimpute/pull/68) - QUILT can handle external params chunks and hap-legend files.
- [#78](https://github.com/nf-core/phaseimpute/pull/78) - Separate validate step from panel preparation.
- [#84](https://github.com/nf-core/phaseimpute/pull/84) - Change depth computation to use SAMTOOLS_DEPTH and make separation by chromosome only if regions are specified.
- [#85](https://github.com/nf-core/phaseimpute/pull/85) - Use external params in individual tests for tools.
- [#86](https://github.com/nf-core/phaseimpute/pull/86) - Move `bcftools_convert` to `vcf_sites_extract_bcftools`.
- [#88](https://github.com/nf-core/phaseimpute/pull/88) - Improve multiqc report with more information.
- [#91](https://github.com/nf-core/phaseimpute/pull/91) - Update metro map with all steps and remove deprecated ones.
- [#93](https://github.com/nf-core/phaseimpute/pull/93) - Add support for CRAM file.
- [#93](https://github.com/nf-core/phaseimpute/pull/93) - Check contigs name at workflow level for BAM and VCF.
- [#93](https://github.com/nf-core/phaseimpute/pull/93) - Samples remove with multiallelics records.
- [#93](https://github.com/nf-core/phaseimpute/pull/93) - Samtools merge in BAM_REGION sbwf.
- [#93](https://github.com/nf-core/phaseimpute/pull/93) - Fix glimpse2_phase output file names.
- [#93](https://github.com/nf-core/phaseimpute/pull/93) - Fix fai combination to fasta.
- [#96](https://github.com/nf-core/phaseimpute/pull/96) - Simplify csv export
- [#96](https://github.com/nf-core/phaseimpute/pull/96) - Use only legend file as posfile for all imputation workflow.
- [#100](https://github.com/nf-core/phaseimpute/pull/100) - Update bcftools, samtools, ... nf-core modules. All indexing is now done with the file creation for most of them.
- [#101](https://github.com/nf-core/phaseimpute/pull/101) - Set `--compute_freq` as `false` by default.
- [#102](https://github.com/nf-core/phaseimpute/pull/102) - Compute chr name from whole vcf.
- [#102](https://github.com/nf-core/phaseimpute/pull/102) - Only warn the user if some contigs are absent from files, the regions to compute is now the intersection of regions, panel, posfile, chunks, map.
- [#102](https://github.com/nf-core/phaseimpute/pull/102) - Update all test and recompute snapshot to match new version of the phaseimpute test dataset.
- [#103](https://github.com/nf-core/phaseimpute/pull/103) - Update Glimpse2 phase, gunzip and multiqc

### `Fixed`

- [#15](https://github.com/nf-core/phaseimpute/pull/15) - Changed test csv files to point to nf-core repository
- [#16](https://github.com/nf-core/phaseimpute/pull/16) - Removed outdir from test config files
- [#65](https://github.com/nf-core/phaseimpute/pull/65) - Separate stitch output by individuals
- [#75](https://github.com/nf-core/phaseimpute/pull/75) - Set frequency computation with VCFFIXUP process as optional with --compute_freq. Use Glimpse_chunk on panel vcf to compute the chunk and not makewindows on fasta.

### `Dependencies`

### `Deprecated`

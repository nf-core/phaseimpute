# nf-core/phaseimpute: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev

### `Added`

- [#175](https://github.com/nf-core/phaseimpute/pull/175) - Add support for all input files in `.json` or `.yaml` format.
- [#181](https://github.com/nf-core/phaseimpute/pull/181) - Add nf-co2footprint plugin to the config file.
- [#184](https://github.com/nf-core/phaseimpute/pull/184) - Add support `.csi` index for `.bam` files.
- [#188](https://github.com/nf-core/phaseimpute/pull/188) - Add documentation for all subworkflows.
- [#210](https://github.com/nf-core/phaseimpute/pull/200) - Add BEAGLE5 support for genotype imputation.
- [#211](https://github.com/nf-core/phaseimpute/pull/211) - Add MINIMAC4 support for genotype imputation.

### `Changed`

- [#166](https://github.com/nf-core/phaseimpute/pull/166) - Bump version to 1.1.0dev and update `CHANGELOG.md`.
- [#170](https://github.com/nf-core/phaseimpute/pull/170) - Update TEMPLATE to nf-core tools version 3.1.2.
- [#175](https://github.com/nf-core/phaseimpute/pull/175) - Update TEMPLATE to nf-core tools version 3.2.0. Move `CHRCHECK` functions to the workflow directory.
- [#182](https://github.com/nf-core/phaseimpute/pull/182) - Add dark version of the metromap and dynamically change it in the README.
- [#185](https://github.com/nf-core/phaseimpute/pull/185) - Add `--sampleNames_file` option for `STICH` and `QUILT`.
- [#187](https://github.com/nf-core/phaseimpute/pull/187) - Update modules and subworkflows.
- [#197](https://github.com/nf-core/phaseimpute/pull/187) - Update to nf-core/tools version 3.3.1 and update nf-test.
- [#188](https://github.com/nf-core/phaseimpute/pull/188) - Change name of `BAM_IMPUTE_GLIMPSE2` to `BAM\_VCF_IMPUTE_GLIMPSE2`.
- [#199](https://github.com/nf-core/phaseimpute/pull/199) - Update all modules to latest nf-core versions.
- [#201](https://github.com/nf-core/phaseimpute/pull/201) - Update TEMPLATE to nf-core tools version 3.3.2.
- [#209](https://github.com/nf-core/phaseimpute/pull/209) - Update TEMPLATE to nf-core tools version 3.4.1.

### `Fixed`

- [#166](https://github.com/nf-core/phaseimpute/pull/166) - Fix depth type to `number` to enable float.
- [#179](https://github.com/nf-core/phaseimpute/pull/179) - Fix VCF usage in `GLIMPSE2`.
- [#183](https://github.com/nf-core/phaseimpute/pull/183) - Remove wrongfully added files in `BAM_EXTRACT_REGION_SAMTOOLS`.
- [#185](https://github.com/nf-core/phaseimpute/pull/185) - Fix CSV generation and check that all mentioned path files exist.
- [#189](https://github.com/nf-core/phaseimpute/pull/189) - Set meta map id as string to avoid error when using numbers in csv files.

### `Dependencies`

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `bcftools` | 1.20        | 1.21        |
| `gunzip`   | 1.10        | 1.13        |
| `lbzip2`   |             | 2.5         |
| `multiqc`  | 1.27        | 1.29        |
| `r-stitch` | 1.6.10      | 1.7.3       |
| `shapeit5` | 1.0.0       | 5.1.1       |
| `vcflib`   | 1.0.3       | 1.0.14      |
| `beagle5`  |             | 5.2         |
| `minimac4` |             | 4.1.6       |

### New contributors

[Gaspard Ichas](https://github.com/gichas)

## v1.0.0 - Black Labrador [2024-12-09]

Initial release of nf-core/phaseimpute, created with the [nf-core](https://nf-co.re/) template.
Special thanks to [Matthias Hörtenhuber](https://github.com/mashehu), [Mazzalab](https://github.com/mazzalab) and [Sofia Stamouli](https://github.com/sofstam) for the review of this release.

### `Added`

- [#20](https://github.com/nf-core/phaseimpute/pull/20) - Added automatic detection of vcf contigs for the reference panel and automatic renaming available.
- [#22](https://github.com/nf-core/phaseimpute/pull/20) - Add validation step for concordance analysis. Input channels changed to match inputs steps. Outdir folder organised by steps. Modules config by subworkflows.
- [#26](https://github.com/nf-core/phaseimpute/pull/26) - Added QUILT method.
- [#47](https://github.com/nf-core/phaseimpute/pull/47) - Add possibility to remove samples from reference panel. Add glimpse2 chunking method. Add full-size test parameters.
- [#58](https://github.com/nf-core/phaseimpute/pull/58) - Add external params posfile and chunks. Add glimpse2 phasing and imputation.
- [#67](https://github.com/nf-core/phaseimpute/pull/67) - Export CSVs from each step.
- [#71](https://github.com/nf-core/phaseimpute/pull/71) - Allow external panel to be used in step impute.
- [#97](https://github.com/nf-core/phaseimpute/pull/97) - Add dog reference panel and config to test pipeline with other species.
- [#102](https://github.com/nf-core/phaseimpute/pull/102) - Add dog panel test.
- [#119](https://github.com/nf-core/phaseimpute/pull/119) - Add dog test with panelprep and imputation.
- [#118](https://github.com/nf-core/phaseimpute/pull/118) - Explain how to customize arguments in the pipeline.
- [#111](https://github.com/nf-core/phaseimpute/pull/111) - Add nf-test for all subworkflow, workflow, modules and functions.
- [#131](https://github.com/nf-core/phaseimpute/pull/131) - Set normalisation as optional. Fix extension detection function. Add support for validation with vcf files. Concatenate vcf only if more than one file. Change `--phased` to `--phase` for consistency.
- [#143](https://github.com/nf-core/phaseimpute/pull/143) - Improve contigs warning and error logging. The number of chromosomes contigs is summarized if above `max_chr_names`.
- [#146](https://github.com/nf-core/phaseimpute/pull/146) - Add `seed` parameter for `QUILT`.
- [#164](https://github.com/nf-core/phaseimpute/pull/164) - Add additional requirement on input schema `"uniqueEntries": ["panel", "chr"]` and `end` should be greater than `start` in regions.

### `Changed`

- [#18](https://github.com/nf-core/phaseimpute/pull/18) - Maps and region by chromosome. Update tests config files. Correct meta map propagation. `test_impute` and `test_sim` works.
- [#19](https://github.com/nf-core/phaseimpute/pull/19) - Changed reference panel to accept a csv, update modules and subworkflows (glimpse1/2 and shapeit5)
- [#40](https://github.com/nf-core/phaseimpute/pull/40) - Add `STITCH` method. Reorganize panelprep subworkflows.
- [#51](https://github.com/nf-core/phaseimpute/pull/51) - Update all process and fix linting errors. Remove `FASTQC` added by the template.
- [#56](https://github.com/nf-core/phaseimpute/pull/56) - Move to nf-test to check the output files names generated. Fix validation and concatenation by chromosomes missing. Add dedicated GLIMPSE1 subworkflow. Fix posfile generation to be done once for glimpse and stitch.
- [#68](https://github.com/nf-core/phaseimpute/pull/68) - `QUILT` can handle external params chunks and hap-legend files.
- [#78](https://github.com/nf-core/phaseimpute/pull/78) - Separate validate step from panel preparation.
- [#84](https://github.com/nf-core/phaseimpute/pull/84) - Change depth computation to use `SAMTOOLS_DEPTH` and make separation by chromosome only if regions are specified.
- [#85](https://github.com/nf-core/phaseimpute/pull/85) - Use external params in individual tests for tools.
- [#86](https://github.com/nf-core/phaseimpute/pull/86) - Move `BCFTOOLS_CONVERT` to `VCF_SITES_EXTRACT_BCFTOOLS`.
- [#88](https://github.com/nf-core/phaseimpute/pull/88) - Improve multiQC report with more information.
- [#91](https://github.com/nf-core/phaseimpute/pull/91) - Update metro map with all steps and remove deprecated ones.
- [#93](https://github.com/nf-core/phaseimpute/pull/93) - Add support for CRAM file.
- [#93](https://github.com/nf-core/phaseimpute/pull/93) - Check contigs name at workflow level for BAM and VCF.
- [#93](https://github.com/nf-core/phaseimpute/pull/93) - Samples remove with multi-allelics records.
- [#93](https://github.com/nf-core/phaseimpute/pull/93) - Samtools merge in `BAM_REGION` subworkflow.
- [#93](https://github.com/nf-core/phaseimpute/pull/93) - Fix glimpse2_phase output file names.
- [#93](https://github.com/nf-core/phaseimpute/pull/93) - Fix fai combination to fasta.
- [#96](https://github.com/nf-core/phaseimpute/pull/96) - Simplify csv export
- [#96](https://github.com/nf-core/phaseimpute/pull/96) - Use only legend file as posfile for all imputation workflow.
- [#100](https://github.com/nf-core/phaseimpute/pull/100) - Update bcftools, samtools, ... nf-core modules. All indexing is now done with the file creation for most of them.
- [#101](https://github.com/nf-core/phaseimpute/pull/101) - Set `--compute_freq` as `false` by default.
- [#102](https://github.com/nf-core/phaseimpute/pull/102) - Compute chr name from whole vcf.
- [#102](https://github.com/nf-core/phaseimpute/pull/102) - Only warn the user if some contigs are absent from files, the regions to compute is now the intersection of regions, panel, posfile, chunks, map.
- [#102](https://github.com/nf-core/phaseimpute/pull/102) - Update all test and recompute snapshot to match new version of the phaseimpute test dataset.
- [#103](https://github.com/nf-core/phaseimpute/pull/103) - Update `GLIMPSE2_PHASE`, `GUNZIP` and `MULTIQC`
- [#135](https://github.com/nf-core/phaseimpute/pull/135) - Impute by batch of 100 individuals by default using `--batch_size` parameter. All individuals BAM files are gathered and VCF are allowed for `GLIMPSE1` and `GLIMPSE2`. Channel preprocessing of stitch is done in stitch subworkflow. Genotype likelihood computation for `GLIMPSE1` is now done outside of the subworkflow and merge the resulting vcf with all the samples. New test added to check batch separation. Improve `usage.md` documentation. Add validation to initialization of the pipeline to ensure compatibility between tools, steps and the files provided by the user.
- [#139](https://github.com/nf-core/phaseimpute/pull/139) - Update all nf-core modules.
- [#146](https://github.com/nf-core/phaseimpute/pull/146) - Remove conda CI check for PR due to Nextflow error.
- [#144](https://github.com/nf-core/phaseimpute/pull/144) - Documentation updates.
- [#148](https://github.com/nf-core/phaseimpute/pull/148) - Fix AWS fulltest github action for manual dispatch.
- [#149](https://github.com/nf-core/phaseimpute/pull/149) - Remove the map file from the AWS fulltest.
- [#152](https://github.com/nf-core/phaseimpute/pull/152) - Fix URLs in the documentation and remove tools citation in the README, use a white background for all images in the documentation.
- [#153](https://github.com/nf-core/phaseimpute/pull/153) - Update and simplify subworkflows snapshot and check only for files names (no md5sum for bam and vcf files due to timestamp).
- [#157](https://github.com/nf-core/phaseimpute/pull/157) - Add `chunk_model` as parameter for better control over `GLIMPSE2_CHUNK` and set window size in `GLIMPSE1_CHUNK` and `GLIMPSE2_chunk` to 4mb to reduce number of chunks (empirical).
- [#160](https://github.com/nf-core/phaseimpute/pull/160) - Improve `CHANGELOG.md` and add details to `usage.md`
- [#158](https://github.com/nf-core/phaseimpute/pull/158) - Remove frequency computation and phasing from full test to reduce cost and computational time.
- [#164](https://github.com/nf-core/phaseimpute/pull/164) - Rename `BAM_REGION_SAMTOOLS` to `BAM_EXTRACT_REGION_SAMTOOLS`. Remove `GLIMPSE2_SPLITREFERENCE` as it is not used. Add more steps to `test_all` profile for more exhaustiveness.
- [#163](https://github.com/nf-core/phaseimpute/pull/163) - Improve configuration for demanding processes. Use Genome in a Bottle VCF benchmarking file for AWS full test. Moved from `glimpse1` to `glimpse2` for the full test profile.
- [#165](https://github.com/nf-core/phaseimpute/pull/165) - Update metro map and add logo to the documentation.

### `Fixed`

- [#15](https://github.com/nf-core/phaseimpute/pull/15) - Changed test csv files to point to nf-core repository.
- [#16](https://github.com/nf-core/phaseimpute/pull/16) - Removed `outdir` from test config files.
- [#65](https://github.com/nf-core/phaseimpute/pull/65) - Separate stitch output by individuals.
- [#75](https://github.com/nf-core/phaseimpute/pull/75) - Set frequency computation with `VCFFIXUP` process as optional with `--compute_freq`. Use `GLIMPSE_CHUNK` on panel vcf to compute the chunk and not makewindows on fasta.
- [#117](https://github.com/nf-core/phaseimpute/pull/117) - Fix directories in CSV.
- [#151](https://github.com/nf-core/phaseimpute/pull/151) - Fix `Type not supported: class org.codehaus.groovy.runtime.GStringImpl` error due to `String` test in `getFileExtension()`.
- [#158](https://github.com/nf-core/phaseimpute/pull/158) - Fix contigs usage when regions is only a subset of the given contigs (e.g. if panel file has the 22 chr and the region file only 2 then only the 2 common will be processed). Fix `multiQC` samples names for better comprehension. Fix `-resume` errors when `ch_fasta` is use by adding `cache = 'lenient'` in necessary processes. Fix `--window-size` of `GLIMPSE_CHUNK` from `4` to `4000000`.
- [#153](https://github.com/nf-core/phaseimpute/pull/153) - Fix getFileExtension function. Fix image in `usage.md`. Fix small warnings and errors with updated language server. `def` has been added when necessary, `:` use instead of `,` in assertions, `_` added to variables not used in closures, `for` loop replaced by `.each{}`, remove unused code / input.
- [#161](https://github.com/nf-core/phaseimpute/pull/161) - Fix `VCF_SPLIT_BCFTOOLS` when only one sample present by updating `BCFTOOLS_PLUGINSPLIT` and adding `BCFTOOLS_QUERY` to get truth samples names for renaming the resulting files.
- [#162](https://github.com/nf-core/phaseimpute/pull/162) - Fix `fai` usage when provided by `genomes` parameter.
- [#164](https://github.com/nf-core/phaseimpute/pull/164) - Improve documentation writing
- [#163](https://github.com/nf-core/phaseimpute/pull/163) - Fix MULTIQC samples names (add post-processing for clean up `FILTER_CHR_DWN`, `FILTER_CHR_INP`, `GAWK_ERROR_SPL`, `GAWK_RSQUARE_SPL`). Fix output panel `publishDir`. Fix java version to `17` in `ci.yml` due to new nextflow version.

### `Dependencies`

| Dependency    | New version |
| ------------- | ----------- |
| `bcftools`    | 1.20        |
| `bedtools`    | 2.31.1      |
| `gawk`        | 5.3.0       |
| `glimpse-bio` | 1.1.1       |
| `glimpse-bio` | 2.0.1       |
| `gunzip`      | 1.10        |
| `htslib`      | 1.21        |
| `multiqc`     | 1.25.1      |
| `r-quilt`     | 1.0.5       |
| `r-stitch`    | 1.6.10      |
| `samtools`    | 1.21        |
| `shapeit5`    | 1.0.0       |
| `tabix`       | 1.11        |
| `vcflib`      | 1.0.3       |

### `Deprecated`

### `Contributors`

[Louis Le Nézet](https://github.com/LouisLeNezet)
[Anabella Trigila](https://github.com/atrigila)
[Eugenia Fontecha](https://github.com/eugeniafontecha)
[Maxime U Garcia](https://github.com/maxulysse)
[Matias Romero Victorica](https://github.com/mrvictorica)
[Nicolas Schcolnicov](https://github.com/nschcolnicov)
[Hemanoel Passarelli](https://github.com/hemanoel)
[Matthias Hörtenhuber](https://github.com/mashehu)
[Sofia Stamouli](https://github.com/sofstam)

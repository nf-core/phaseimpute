<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-phaseimpute_logo_dark.png">
    <img alt="nf-core/phaseimpute" src="docs/images/nf-core-phaseimpute_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/phaseimpute/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/phaseimpute/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/phaseimpute/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/phaseimpute/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/phaseimpute/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/phaseimpute)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23phaseimpute-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/phaseimpute)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/phaseimpute** is a bioinformatics pipeline to phase and impute genetic data. Different steps are available each corresponding to a dedicated modes.

### Main steps of the pipeline

The **phaseimpute** pipeline is constituted of 5 main steps:

| Metro map                                                              | Modes                                                                                                                                                                                                                                                                                                                                             |
| ---------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="docs/images/metro/MetroMap.png" alt="metromap" width="800"/> | - **Panel preparation**: Phasing, QC, variant filtering, variant annotation of the reference panel <br> - **Imputation**: Impute the target dataset on the reference panel <br> - **Simulate**: Simulation of the target dataset from high quality target data <br> - **Concordance**: Concordance between the target dataset and a truth dataset |

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

The basic usage of this pipeline is to impute a target dataset based on a phased panel.
First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,file,index
1_BAM_1X,/path/to/.bam,/path/to/.bai
1_CRAM_1X,/path/to/.cram,/path/to/.crai
```

Each row represents a bam file with its index file. For some tools and steps, you will also need to submit a samplesheet with the reference panel.

A final samplesheet file for the reference panel may look something like the one below. This is for 3 chromosomes.

```csv
panel,chr,vcf,index
Phase3,1,ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz,ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.csi
Phase3,2,ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz,ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.csi
Phase3,3,ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz,ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.csi
```

Now, you can run the pipeline using:

```bash
nextflow run nf-core/phaseimpute \
   -profile <docker/singularity/.../institute> \
   --input <samplesheet.csv>  \
   --genome "GRCh38" \
   --panel <phased_reference_panel.csv> \
   --steps "panelprep,impute" \
   --tools "glimpse1" \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/phaseimpute/usage) and the [parameter documentation](https://nf-co.re/phaseimpute/parameters).

## Description of the different steps of the pipeline

Here is a short description of the different steps of the pipeline.
For more information please refer to the [documentation](https://nf-core.github.io/phaseimpute/usage/).

| steps           | Flow chart                                                                       | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| --------------- | -------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **--panelprep** | <img src="docs/images/metro/PanelPrep.png" alt="Panel preparation" width="600"/> | The preprocessing mode is responsible to the preparation of the multiple input file that will be used by the phasing process. <br> The main processes are : <br> - **Haplotypes phasing** of the reference panel using [**Shapeit5**](https://odelaneau.github.io/shapeit5/). <br> - **Normalize** the reference panel to select only the necessary variants. <br> - **Chunking the reference panel** in a subset of region for all the chromosomes. <br> - **Extract** the positions where to perform the imputation.                                                                                                                                                                                                                                                                                                                                                                           |
| **--impute**    | <img src="docs/images/metro/Impute.png" alt="Impute target" width="600"/>        | The imputation mode is the core mode of this pipeline. <br> It is constituted of 3 main steps: <br> - **Imputation**: Impute the target dataset on the reference panel using either: <br> &emsp; - [**Glimpse1**](https://odelaneau.github.io/GLIMPSE/glimpse1/index.html): It's come with the necessety to compute the genotype likelihoods of the target dataset (done using [BCFTOOLS_mpileup](https://samtools.github.io/bcftools/bcftools.html#mpileup)). <br> &emsp; - [**Glimpse2**](https://odelaneau.github.io/GLIMPSE/glimpse2/index.html) <br> &emsp; - [**Stitch**](https://github.com/rwdavies/stitch) This steps does not require a reference panel but needs to merge the samples. <br> &emsp; - [**Quilt**](https://github.com/rwdavies/QUILT) <br> - **Ligation**: all the different chunks are merged together then all chromosomes are reunited to output one VCF per sample. |
| **--simulate**  | <img src="docs/images/metro/Simulate.png" alt="simulate_metro" width="600"/>     | The simulation mode is used to create artificial low informative genetic information from high density data. This allow to compare the imputed result to a _truth_ and therefore evaluate the quality of the imputation. <br> For the moment it is possible to simulate: <br> - Low-pass data by **downsample** BAM or CRAM using [SAMTOOLS_VIEW -s](https://www.htslib.org/doc/samtools-view.html) at different depth.                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| **--validate**  | <img src="docs/images/metro/Validate.png" alt="concordance_metro" width="600"/>  | This mode compare two vcf together to compute a summary of the differences between them. <br> This step use [**Glimpse2**](https://odelaneau.github.io/GLIMPSE/glimpse2/index.html) concordance process.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/phaseimpute/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/phaseimpute/output).

## Credits

nf-core/phaseimpute was originally written by Louis Le Nézet.

We thank the following people for their extensive assistance in the development of this pipeline:

- Anabella Trigila
- Saul Pierotti
- Eugenia Fontecha
- Matias Romero Victorica

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#phaseimpute` channel](https://nfcore.slack.com/channels/phaseimpute) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/phaseimpute for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

You can cite one of the main imputation methods ([`QUILT`](https://github.com/rwdavies/QUILT)) as follows:

> **Rapid genotype imputation from sequence with reference panels.**
>
> Davies, R. W., Kucka, M., Su, D., Shi, S., Flanagan, M., Cunniff, C. M., Chan, Y. F., & Myers, S.
>
> _Nature genetics_ 2021 June 03. doi: [10.1038/s41588-021-00877-0](https://doi.org/10.1038/s41588-021-00877-0)

You can cite one of the main imputation methods ([`GLIMPSE`](https://github.com/odelaneau/GLIMPSE)) as follows:

> **Efficient phasing and imputation of low-coverage sequencing data using large reference panels.**
>
> Rubinacci, S., Ribeiro, D. M., Hofmeister, R. J., & Delaneau, O.
>
> _Nature genetics_ 2021. doi:[]()

> **Imputation of low-coverage sequencing data from 150,119 UK Biobank genomes**
>
> Rubinacci, S., Hofmeister, R. J., Sousa da Mota, B., & Delaneau, O.
>
> _Nature genetics_ 2023. doi:[]()

You can cite one of the main imputation methods ([`STITCH`](https://github.com/rwdavies/STITCH)) as follows:

> **Rapid genotype imputation from sequence without reference panels.**
>
> Davies, R. W., Flint, J., Myers, S., & Mott, R.
>
> _Nature genetics_ 2016 . doi: []().

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

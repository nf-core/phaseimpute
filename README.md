<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-phaseimpute_logo_dark.png">
    <img alt="nf-core/phaseimpute" src="docs/images/nf-core-phaseimpute_logo_light.png">
  </picture>
</h1>

**Multi-step pipeline dedicated to genetic imputation from simulation to validation**

[![GitHub Actions CI Status](https://github.com/nf-core/phaseimpute/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/phaseimpute/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/phaseimpute/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/phaseimpute/actions/workflows/linting.yml) [![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/phaseimpute/results) [![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX) [![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/phaseimpute)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23phaseimpute-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/phaseimpute) [![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core) [![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core) [![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/phaseimpute** is a bioinformatics pipeline to phase and impute genetic data.

<img src="docs/images/metro/MetroMap_animated.svg" alt="metromap" width="800"/>

The pipeline is constituted of five main steps:

1. **Check chromosomes names**: Validates the presence of the different contigs in all variants and alignment files, ensuring data compatibility for further processing.

2. **Panel preparation**: Performs phasing, QC, variant filtering, and variant annotation of the reference panel.

3. **Imputation**: Imputes genotypes in the target dataset using the reference panel.

4. **Simulation**: Generates simulated datasets from high-quality target data for testing and validation purposes.

5. **Concordance**: Evaluates the accuracy of imputation by comparing the imputed data against a truth dataset.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

The primary function of this pipeline is to impute a target dataset based on a phased panel. Begin by preparing a samplesheet with your input data, formatted as follows:

```csv title="samplesheet.csv"
sample,file,index
SAMPLE_1X,/path/to/.<bam/cram>,/path/to/.<bai,crai>
```

Each row represents either a bam or a cram file along with its corresponding index file. Ensure that all input files have consistent file extensions.

For certain tools and steps within the pipeline, you will also need to provide a samplesheet for the reference panel. Here's an example of what a final samplesheet for a reference panel might look like, covering three chromosomes:

```csv title="panel.csv"
panel,chr,vcf,index
Phase3,1,ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz,ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.csi
Phase3,2,ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz,ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.csi
Phase3,3,ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz,ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.csi
```

## Running the pipeline

Execute the pipeline with the following command:

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
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/phaseimpute/usage) and the [parameter documentation](https://nf-co.re/phaseimpute/parameters).

## Description of the different steps of the pipeline

Here is a short description of the different steps of the pipeline.
For more information please refer to the [usage documentation](https://nf-co.re/phaseimpute/usage).

1. **`--panelprep`**: The preprocessing mode is responsible for preparing multiple input files that will be used by the phasing and imputation process.
   - The main processes are:
     - **Haplotypes phasing** of the reference panel using [**Shapeit5**](https://odelaneau.github.io/shapeit5/).
     - **Normalize** the reference panel to select only the necessary variants.
     - **Chunking the reference panel** into a subset of regions for all the chromosomes.
     - **Extract** the positions where to perform the imputation.

<img src="docs/images/metro/PanelPrep.png" alt="Panel preparation" width="600"/>

2. **`--impute`**: The imputation mode is the core mode of this pipeline.
   - It consists of 3 main steps:
     - **Imputation**: Impute the target dataset on the reference panel using either:
       - [**Glimpse1**](https://odelaneau.github.io/GLIMPSE/glimpse1/index.html).
       - [**Glimpse2**](https://odelaneau.github.io/GLIMPSE/).
       - [**Stitch**](https://github.com/rwdavies/stitch).
       - [**Quilt**](https://github.com/rwdavies/QUILT)
     - **Ligation**: Merges the different chunks together and reunites all chromosomes to output one VCF per sample.

<img src="docs/images/metro/Impute.png" alt="Impute target" width="600"/>

3. **`--simulate`**: The simulation mode is used to create artificial low informative genetic information from high-density data. This allows for the comparison of the imputed result to a _truth_ dataset, evaluating the quality of the imputation.
   - Currently, it is possible to simulate:
     - Low-pass data by **downsampling** BAM or CRAM using [`samtools view -s`](https://www.htslib.org/doc/samtools-view.html) at different depths.

<img src="docs/images/metro/Simulate.png" alt="simulate_metro" width="600"/>

4. **`--validate`**: This mode compares two VCF files to compute a summary of the differences between them.
   - This step uses the [**Glimpse2**](https://odelaneau.github.io/GLIMPSE/) concordance process.

<img src="docs/images/metro/Validate.png" alt="concordance_metro" width="600"/>

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/phaseimpute/results) tab on the nf-core website pipeline page.
For more details on the output files and reports, please refer to the [output documentation](https://nf-co.re/phaseimpute/output).

## Credits

nf-core/phaseimpute was originally written by Louis Le Nézet & Anabella Trigila.

We thank the following people for their extensive assistance in the development of this pipeline:

- Saul Pierotti
- Eugenia Fontecha
- Matias Romero Victorica
- Hemanoel Passarelli

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#phaseimpute` channel](https://nfcore.slack.com/channels/phaseimpute) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->

If you use nf-core/phaseimpute for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX)

An extensive list of references for the tools used by the pipeline, including QUILT, GLIMPSE, and STITCH, can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

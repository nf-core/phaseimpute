<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/logo/nf-core-phaseimpute_logo_dark.png">
    <img alt="nf-core/phaseimpute" src="docs/images/logo/nf-core-phaseimpute_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new/nf-core/phaseimpute)
[![GitHub Actions CI Status](https://github.com/nf-core/phaseimpute/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/phaseimpute/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/phaseimpute/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/phaseimpute/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/phaseimpute/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.14329225-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.14329225)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.4.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.4.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/phaseimpute)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23phaseimpute-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/phaseimpute)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/phaseimpute** is a bioinformatics pipeline to phase and impute genetic data.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/images/metro/MetroMap_animated_dark.svg">
  <img src="docs/images/metro/MetroMap_animated.svg" alt="metromap"/>
</picture>

The whole pipeline consists of five main steps, each of which can be run separately and independently. Users are not required to run all steps sequentially and can select specific steps based on their needs:

1. **QC: Chromosome Name Check**: Ensures compatibility by validating that all expected contigs are present in the variant and alignment files.

2. **Simulation (`--simulate`)**: Generates artificial datasets by downsampling high-density data to simulate low-pass genetic information. This enables the comparison of imputation results against a high-quality dataset (truth set). Simulations may include:
   - **Low-pass data generation** by downsampling BAM or CRAM files with [`samtools view -s`](https://www.htslib.org/doc/samtools-view.html) at different depths.

3. **Panel Preparation (`--panelprep`)**: Prepares the reference panel through phasing, quality control, variant filtering, and annotation. Key processes include:
   - **Normalization** of the reference panel to retain essential variants.
   - **Phasing** of haplotypes in the reference panel using [Shapeit5](https://odelaneau.github.io/shapeit5/).
   - **Chunking** of the reference panel into specific regions across chromosomes.
   - **Position Extraction** for targeted imputation sites.

4. **Imputation (`--impute`)**: This is the primary step, where genotypes in the target dataset are imputed using the prepared reference panel. The main steps are:
   - **Imputation** of the target dataset using tools like [Glimpse1](https://odelaneau.github.io/GLIMPSE/glimpse1/index.html), [Glimpse2](https://odelaneau.github.io/GLIMPSE/), [Stitch](https://github.com/rwdavies/stitch), [Quilt](https://github.com/rwdavies/QUILT), [Beagle5](https://faculty.washington.edu/browning/beagle/beagle.html) or [Minimac4](https://github.com/statgen/Minimac4).
   - **Ligation** of imputed chunks to produce a final VCF file per sample, with all chromosomes unified.

5. **Validation (`--validate`)**: Assesses imputation accuracy by comparing the imputed dataset to a truth dataset. This step leverages the [Glimpse2](https://odelaneau.github.io/GLIMPSE/) concordance process to summarize differences between two VCF files.

For more detailed instructions, please refer to the [usage documentation](https://nf-co.re/phaseimpute/usage).

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

Run one of the steps of the pipeline (imputation with glimpse1) using the following command and test profile:

```bash
nextflow run nf-core/phaseimpute \
   -profile test, <docker/singularity/.../institute> \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/phaseimpute/usage) and the [parameter documentation](https://nf-co.re/phaseimpute/parameters).

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

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md). Further development tips can be found in the [development documentation](docs/development.md).

For further information or help, don't hesitate to get in touch on the [Slack `#phaseimpute` channel](https://nfcore.slack.com/channels/phaseimpute) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/phaseimpute for your analysis, please cite it using the following doi: [10.5281/zenodo.14329225](https://doi.org/10.5281/zenodo.14329225)

An extensive list of references for the tools used by the pipeline, including QUILT, GLIMPSE, and STITCH, can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

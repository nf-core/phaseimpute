# nf-core/phaseimpute: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

## Panel preparation outputs `--step panelprep`

This step of the pipeline performs a QC of the reference panel data and produces the necessary files for imputation (`--step impute`). It has two optional modes: reference panel phasing with SHAPEIT5 and removal of specified samples from reference panel.

- [Remove Multiallelics](#multiallelics) - Remove multiallelic sites from the reference panel
- [Convert](#convert) - Convert reference panel to .hap and .legend files
- [Posfile](#posfile) - Produce a TSV with the list of positions to genotype (for STITCH/QUILT)
- [Sites](#sites) - Produce a TSV with the list of positions to genotype (for GLIMPSE1)
- [Glimpse Chunk](#glimpse) - Create chunks of the reference panel

### Convert

- `prep_panel/haplegend/`
  - `*.hap`: a .hap file for the reference panel.
  - `*.legend*`: a .legend file for the reference panel.

[bcftools](https://samtools.github.io/bcftools/bcftools.html) aids in the conversion of vcf files to .hap and .legend files. A .samples file is also generated. Once that you have generated the hap and legend files for your reference panel, you can skip the reference preparation step and directly submit these files for imputation (to be developed). The hap and legend files are input files used with `--tools quilt`.

### Posfile

- `prep_panel/posfile/`
  - `*.hap`: a .txt file with the list of position to genotype.

[bcftools query](https://samtools.github.io/bcftools/bcftools.html) produces tab-delimited files per chromosome that can be gathered into a samplesheet and directly submitted for imputation with `--tools stitch` using the parameter `--posfile`.

### Sites

- `prep_panel/sites/`
  - `vcf/`
    - `*.vcf.gz`: VCF with biallelic SNPs only.
    - `*.csi`: Index file for VCF.
  - `tsv/`
    - `*.txt.gz`: TXT file for biallelic SNPs.
    - `*.tbi`: Index file for TSV.

[bcftools query](https://samtools.github.io/bcftools/bcftools.html) produces VCF (`*.vcf.gz`) files per chromosome. These QCed VCFs can be gathered into a csv and used with all the tools in `--step impute` using the flag `--panel`.

In addition, [bcftools query](https://samtools.github.io/bcftools/bcftools.html) produces tab-delimited files (`*_tsv.txt`) and, together with the VCFs, they can be gathered into a samplesheet and directly submitted for imputation with `--tools glimpse1` and `--posfile` (not yet implemented).

### Glimpse Chunk

- `prep_panel/chunks/`
  - `*.txt`: TXT file containing the chunks obtained from running Glimpse chunks.

[Glimpse1 chunk](https://odelaneau.github.io/GLIMPSE/) defines chunks where to run imputation. For further reading and documentation see the [Glimpse1 documentation](https://odelaneau.github.io/GLIMPSE/glimpse1/commands.html). Once that you have generated the chunks for your reference panel, you can skip the reference preparation step and directly submit this file for imputation.

## QUILT imputation mode

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [QUILT](#quilt) - Perform imputation
- [Concatenate](#concatenate) - Concatenate all imputed chunks into a single VCF.

### QUILT

- `imputation/quilt/`
- `quilt.*.vcf.gz`: Imputed VCF for a specific chunk.
- `quilt.*.vcf.gz.tbi`: TBI for the Imputed VCF for a specific chunk.

[quilt](https://github.com/rwdavies/QUILT) performs the imputation. This step will contain the VCF for each of the chunks.

### Concat

- `imputation/quilt/bcftools/concat`
- `.*.vcf.gz`: Imputed and ligated VCF for all the input samples.

[bcftools concat](https://samtools.github.io/bcftools/bcftools.html) will produce a single VCF from a list of imputed VCFs in chunks.

## STITCH imputation mode

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [STITCH](#stitch) - Perform imputation
- [Concatenate](#concatenate) - Concatenate all imputed chunks into a single VCF

### STITCH

- `imputation/stitch/`
- `stitch.*.vcf.gz`: Imputed VCF for a specific chunk.
- `stitch.*.vcf.gz.tbi`: TBI for the Imputed VCF for a specific chunk.

[STITCH](https://github.com/rwdavies/STITCH) performs the imputation. This step will contain the VCF for each of the chunks.

### Concat

- `imputation/stitch/bcftools/concat`
- `.*.vcf.gz`: Imputed and concatenated VCF for all the input samples.

[bcftools concat](https://samtools.github.io/bcftools/bcftools.html) will produce a single VCF from a list of imputed VCFs.

## GLIMPSE2 imputation mode

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [GLIMPSE2](#glimpse2) - Perform imputation
- [Concatenate](#concatenate) - Concatenate all imputed chunks into a single VCF

### GLIMPSE2 output files

- `imputation/glimpse2/concat`
- `.*.vcf.gz`: Imputed and concatenated VCF for all the input samples.

## Reports

Reports contain useful metrics and pipeline information for the different modes.

- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.
[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

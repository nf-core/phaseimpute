# nf-core/phaseimpute: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

## Panel preparation outputs `--steps panelprep`

This step of the pipeline performs a QC of the reference panel data and produces the necessary files for imputation (`--steps impute`).

It has the following optional modes:

- `--normalize` - Normalize the reference panel with `bcftools norm` and remove multiallelic sites. It also allow to remove samples using `--remove_samples`.
- `--compute_freq` - Compute allele frequencies with `vcffixup`.
- `--phase` - Phase the reference panel with `SHAPEIT5`.

The pipeline will produce the following outputs:

- [Normalize reference panel](#panel-directory) - Remove multiallelic sites from the reference panel and compute allele frequencies if needed.
- [Convert](#haplegend-directory) - Convert reference panel to `.hap` and `.legend` files.
- [Posfile](#sites-directory) - Produce a `.tsv` file with the list of positions to genotype for the different tools.
- [Chromosomes chunks](#chunks-directory) - Create chunks of the reference panel.
- [CSV](#csv-directory) - Obtained `.csv` files from this step.

The directory structure from `--steps panelprep` is:

```tree
├── panel
├── haplegend
├── sites
├── chunks
│   ├── glimpse1
│   └── glimpse2
├── csv
```

### Panel directory

<details markdown="1">
<summary>Output files</summary>

- `prep_panel/panel/`
  - `*.vcf.gz`: The reference panel VCF files after all the preprocessing steps are completed.
  - `*.tbi`: The index file for the prepared reference panel.

</details>

A directory containing the reference panel per chromosome after preprocessing.
The files will be normalized if the flag `--normalize` is used (with `_normalized` suffix). The files will have their allele frequency computed if the flaq `--compute_freq` is used (with `_fixup` suffix).
The files will be phased if the flag `--phase` is used (with `_phased` suffix).

### Haplegend directory

<details markdown="1">
<summary>Output files</summary>

- `prep_panel/haplegend/`
  - `*.hap`: a `.hap` file for the reference panel containing the genotype.
  - `*.legend*`: a `.legend` file for the reference panel containing the variants informations.
  - `*.samples`: a `.samples` file for the reference panel containing the samples informations.

</details>

[`bcftools convert`](https://samtools.github.io/bcftools/bcftools.html#convert) aids in the conversion of VCF files to `.hap` and `.legend` files. A `.samples` file is also generated. Once that you have generated the hap and legend files for your reference panel, you can skip the reference preparation steps and directly submit these files for imputation. The hap and legend files can be used as input files with the `--tools quilt` option.

### Sites directory

<details markdown="1">
<summary>Output files</summary>

- `prep_panel/sites/`
  - `*.vcf.gz`: A VCF file with biallelic SNPs only.
  - `*.csi`: Index file of the VCF file.

</details>

[`bcftools query`](https://samtools.github.io/bcftools/bcftools.html#query) produces VCF (`*.vcf.gz`) files per chromosome. These QCed VCF files can be gathered into a CSV file and used with all the tools in `--steps impute` using the flag `--panel`.

### Chunks directory

<details markdown="1">
<summary>Output files</summary>

- `prep_panel/chunks/`
  - `*.txt`: Text file containing the chunks obtained after running `GLIMPSE1_CHUNK`.

</details>

[Glimpse1 chunk](https://odelaneau.github.io/GLIMPSE/glimpse1/) defines the chunks where imputation will be performed. For further reading and documentation see the [Glimpse1 documentation](https://odelaneau.github.io/GLIMPSE/glimpse1/commands.html). Once you have generated the chunks for your reference panel, you can skip the reference preparation steps and directly submit this file for imputation.

### CSV directory

<details markdown="1">
<summary>Output files</summary>

- `prep_panel/csv/`
  - `chunks_glimpse1.csv`: A CSV file containing the list of chunks obtained for each chromosome and panel.
  - `panel.csv`: A CSV file containing the final phased and prepared for each chromosome and input panel.
  - `posfile.csv`: A CSV file containing the final list of panel positions, in VCF and TSV files, for each chromosome and input panel.

</details>

## Imputation outputs `--steps impute`

The results from `--steps impute` will have the following directory structure:

```tree
├── batch
├── csv
├── <glimpse1|glimpse2|quilt|stitch|beagle5|minimac4>
│   ├── concat/
│   └── samples/
├── stats
```

<details markdown="1">
<summary>Output files</summary>

- `imputation/batch/all.batchi.id.txt`: List of samples names processed in the i^th^ batch.
- `imputation/csv/`
  - `impute.csv`: A single CSV file containing the path to a VCF file and its index, of each imputed sample with their corresponding tool.
- `imputation/[glimpse1,glimpse2,quilt,stitch]/`
  - `concat/all.batch*.vcf.gz`: The concatenated VCF files of all imputed samples by batches.
  - `concat/all.batch*.vcf.gz.tbi`: The index file for the concatenated imputed VCF files of the samples.
  - `samples/*.vcf.gz`: A VCF file of each imputed sample.
  - `samples/*.vcf.gz.tbi`: The index file of the imputed VCF files.
- `imputation/*.<tool>.bcftools_stats.txt`: The statistics of the imputed VCF target file produced by [`BCFTOOLS_STATS`](https://samtools.github.io/bcftools/bcftools.html#stats.)

</details>

[`bcftools concat`](https://samtools.github.io/bcftools/bcftools.html#concat) will produce a single VCF file from a list of imputed VCF files in chunks.

## Simulation outputs `--steps simulate`

The results from `--steps simulate` will have the following directory structure:

```tree
├── csv
├── samples
```

<details markdown="1">
<summary>Output files</summary>

- `simulation/`
  - `csv`:
    - `simulate.csv`: Samplesheet listing all downsampled target alignment files.
  - `*.depth_*x.bam`: An alignment file from the target file downsampled at the desired depth.
  - `*.bam.csi`: The corresponding index of the alignment file.

</details>

## Validation outputs `--steps validate`

The results from `--steps validate` will have the following directory structure:

```tree
├── concat
├── samples
├── stats
```

<details markdown="1">
<summary>Output files</summary>

- `validation/`
  - `concat/all.truth.vcf.gz`: The concatenated VCF file of all truth sample.
  - `concat/all.truth.vcf.gz.tbi`: The index file of the concatenated truth VCF file of the samples.
  - `samples/*.vcf.gz`: A VCF file of each truth sample.
  - `samples/*.vcf.gz.tbi`: The index file of the truth VCF file.
  - `stats/`:
    - `*.truth.bcftools_stats.txt`: The statistics of the truth VCF target file produced by [`BCFTOOLS_STATS`](https://samtools.github.io/bcftools/bcftools.html#stats.)
    - `*.P<panel name>_T<imputation tool>_SNP.txt`: Concordance metrics of the SNPs variants obtained with [`GLIMPSE2_CONCORDANCE`](https://odelaneau.github.io/GLIMPSE/docs/documentation/concordance/).
    - `AllSamples.txt`: Aggregation of the above `GLIMPSE_CONCORDANCE` output across samples and tools.

</details>

## Reports

Reports contain useful metrics and pipeline information for the different modes.

- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline.
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Parameters used by the pipeline run: `params.json`.
  - Report generated by [`nf-co2footprint`](https://nextflow-io.github.io/nf-co2footprint/): `co2footprint_trace.txt`, `co2footprint_summary.txt` and `co2footprint_report`.html.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

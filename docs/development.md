# Development

## Style

Names of releases are composed of a color + a dog breed.

## Features and tasks

- [x] Add automatic detection of chromosome name to create a renaming file for the vcf files
- [] Add support of `anyOf()` or `oneOf()` in the nf-core schema for the map, panel and region files
- [] Add nf-test for all modules and subworkflows
- [] Remove all TODOs
- [] Check if panel is necessary depending on the tool selected
- [x] Set modules configuration as full path workflow:subworkflow:module
- [] Add support for the map files
- [x] Add support for imputation by individuals or by groups of individuals

## Run tests

### Launch with Nextflow

```bash
nextflow run main.nf -profile singularity,test --outdir results -resume
nextflow run main.nf -profile singularity,test_sim --outdir results -resume
nextflow run main.nf -profile singularity,test_validate --outdir results -resume
nextflow run main.nf -profile singularity,test_all --outdir results -resume
nextflow run main.nf -profile singularity,test_quilt --outdir results -resume
```

### Launch with nf-test

```bash
nf-test test --verbose --profile singularity --tag test_all
nf-test test --verbose --profile singularity --tag test_all --update-snapshot #To update the snaps of a given test
```

## Problematic

### Channel management and combination

If only one specie at a time, then only one fasta file and only one map file (normally ?)
Do we want to be able to compute multiple panel at the same time ?
If so we need to correctly combine the different channel depending on their meta map.

All channel need to be identified by a meta map as follow:

- I : individual id
- P : panel id
- R : region used
- M : map used
- T : tool used
- G : reference genome used (is it needed ?)
- S : simulation (depth or genotype array)

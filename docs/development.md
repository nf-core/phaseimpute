# Development

## Features and tasks

- [x] Add automatic detection of chromosome name to create a renaming file for the vcf files
- [] Add automatic detection of chromosome name to create a renaming file for the bam files
- [] Make the different tests workflows work
  - [x] Simulation
  - [x] Validation
  - [] Preprocessing
  - [x] Imputation
  - [] Validation
  - [] Postprocessing
- [] Add support of `anyOf()` or `oneOf()` in the nf-core schema for the map, panel and region files
- [] Add nf-test for all modules and subworkflows
- [] Remove all TODOs
- [] Check if panel is necessary depending on the tool selected
- [x] Set modules configuration as full path workflow:subworkflow:module
- [] Where should the map file go (separate csv or in panel csv)
- [] Add support for imputation by individuals or by groups of individuals

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
nf-test test --verbose --profile singularity --tag test_all --update-snap #To update the snaps of a given test
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

## Open questions

How to use different schema ?

- Use nf-validation
  For the moment use different input / step.
  In the futur, if/else logic will be added in the yml nf-core schema.

What's the use of dumpcustomsoftware ?
Will be deleted

How to add to multiQC ?
Take exemple on Sarek.
All report file are in a dedicated channel.

How to add nf-test ?
Add in `tests` folder and run with tag.
Add tags.yml

How to run stub tests ?
Use nf-test

How to run the tests ?
nf-test option tag

What's the use of the template branch ?
TEMPLATE branch have the skeleton for all common part of the pipeline.
Will be asked to be merged to dev from time to time.

When is it necessary to merge to master / main ?
First release, create a false PR to first commit that will be checked by whole community + 2 reviewers approval.

What should be the Github action ?
All GA come from the TEMPLATE branch.

# Development

To contribute to this pipeline you will need to install the development environment:
This is possible only on linux or MacOs machine as Nextflow only work on these platform.

```bash
conda env create -f environment.yml
conda activate nf-core-phaseimpute-1.0dev
```

## Add new module

```bash
nf-core modules install
```

## Run tests

```bash
nextflow run main.nf -profile singularity,test --outdir results -resume
```

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

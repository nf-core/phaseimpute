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

## Open questions

How to use different schema ?
- Use nf-validation
- Use If Else statement

What's the use of dumpcustomsoftware ?

How to add to multiQC ?

How to add nf-test ?

How to run stub tests ?

How to run the tests ?

What's the use of the template branch ?

When is it necessary to merge to master / main ?

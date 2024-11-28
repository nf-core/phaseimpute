# Tips for development

## Channel management and combination

All channels need to be identified by a meta map. To follow which information is available, the `meta` argument
is suffixed with a combination of the following capital letters:

- I : individual id
- P : panel id
- R : region used
- M : map used
- T : tool used
- G : reference genome used (is it needed ?)
- S : simulation (depth or genotype array)

Therefore, the following channel operation example includes a meta map containing the panel id with the region and tool used:

```nextflow
ch_panel_for_impute.map {
    metaPRT, vcf, index -> ...
}
```

## Release names

The names of releases are composed of a color and a dog breed.

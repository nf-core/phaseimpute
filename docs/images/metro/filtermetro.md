# How to generate and filter metro map.

## What we want

- Easily updatable diagram
- Expandable / filter at different depth
- Customizable
- Non used channel not visible (cleaning)
- Correct names
- Easy to read and identify, workflow, subworkflow, function, process, files
- Easy to use and configurable

## Available solutions with pro and cons

- `-with-dag` generate a mermaid diagram mostly complete
  - Pros : Easy to use, already implemented, dependant of the workflow
  - Cons : Difficult to read, messy (many not used channels), not exhaustive of the pipeline

## Run pipeline with dag

```bash
nextflow run main.nf -profile test_sim,singularity --outdir results -with-dag -preview
```

## Extract mermaid from html

```bash
python docs/images/metro/metro.py -f results/pipeline_info/pipeline_dag_2024-05-25_23-27-14.html
```

````bash
#Get last html
htmlfile=$(find results/pipeline_info -name "pipeline_dag_*.html" -printf "%T@ %p\n" | sort -n | tail -1 | awk '{print $2}')
#Extract mermaid
content=$(sed -n '/<pre class="mermaid"/,/<\/pre>/p' $htmlfile \
    | sed 's/<pre class="mermaid" style="text-align: center;">//g' \
    | sed 's/<\/pre>//g')
#Register into markdown
mdfile="docs/images/metro/mermaid.md"
touch $mdfile
echo '```mermaid' > $mdfile
echo "$content" >> $mdfile
echo '```' >> $mdfile
````

## Filter

### Extract all empty

```bash
grep -oP 'v\d+\[" "\]' $mdfile | sed 's/\[" "\]//g' > docs/images/metro/empty.txt
```

### Extract all relation

```bash
grep -oP 'v\d+ --> v\d+' $mdfile | sed 's/ --> /\t/g' > docs/images/metro/relationships.txt
```

### Structure

- Flowchart
  - Subgraph1
    - Name
    - Type
    - Content
      - Process1
        - Name
        - Type
      - Process2
        - Name
        - Type

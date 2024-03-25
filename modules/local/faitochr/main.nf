process FAITOCHR {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(fai)

    output:
    tuple val(meta), path("*.txt"), emit: annot_chr
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Check if chr prefix is present in the chromosome names
    col1="chr"
    col2=""
    if [ \$(awk 'NR==1 {print \$1}' ${fai} | grep -c '^chr') -eq 1 ]; then
        col1=""
        col2="chr"
    fi

    # Take the fai file and add/remove the chr prefix to the chromosome names
    # Keep only first column, remove chr prefix if present, add chr prefix if needed
    # chr prefix is added only on number only chromosome names or XYMT
    awk -F'\t' '{print \$1}'  ${fai} | \
        sed 's/^chr//g' | \
        awk -v col1=\${col1} -v col2=\${col2} \
            'BEGIN {OFS=" "} {if (\$1 ~ /^[0-9]+|[XYMT]\$/) print col1\$1, col2\$1; else print \$1, \$1}' \
        > ${prefix}.txt

    # We should have a file with the chromosome names in the second column corresponding to the fai format

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | grep -o 'GNU Awk [0-9.]*' | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | grep -o 'GNU Awk [0-9.]*' | cut -d ' ' -f 3)
    END_VERSIONS
    """
}

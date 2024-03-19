process FAITOCHR {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(fai), val(addchr)

    output:
    tuple val(meta), path("*.txt"), emit: annot_chr
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Take the fai file and add the chr prefix to the chromosome names
    if [ "${addchr}" = true ]; then
        col1=""
        col2="chr"
    else
        col1="chr"
        col2=""
    fi
    awk -F'\t' '{print \$1}'  ${fai} | \
        sed 's/chr//g' | \
        awk -v col1=\${col1} -v col2=\${col2} 'BEGIN {OFS=" "} {print col1\$1, col2\$1}' > ${prefix}.txt

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

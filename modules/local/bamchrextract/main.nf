process BAMCHREXTRACT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.txt"), emit: chr
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        head \\
        $input \| \\
        grep '^@SQ' | cut -d$'\t' -f2 \\
        > ${prefix}.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version |& sed '1!d; s/^.*samtools //' )
        grep: \$( grep --version |& grep -o -E '[0-9]+\\.[0-9]+' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version |& sed '1!d; s/^.*samtools //' )
        grep: \$( grep --help |& grep -o -E '[0-9]+\\.[0-9]+\\.[0-9]+' )
    END_VERSIONS
    """
}

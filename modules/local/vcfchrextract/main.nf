process VCFCHREXTRACT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.18--h8b25389_0':
        'biocontainers/bcftools:1.18--h8b25389_0' }"

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
    bcftools \\
        head \\
        $input \\
        \| grep -o -E '^##contig=<ID=([^,>]*)' | cut -d'=' -f3 \\
        > ${prefix}.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
        grep: \$( grep --version |& sed '1!d; s/^.*grep (GNU grep) //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
        grep: \$( grep --version |& sed '1!d; s/^.*grep (GNU grep) //' )
    END_VERSIONS
    """
}

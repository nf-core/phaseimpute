process CONCATENATE {
    label 'process_single'

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path('*.txt'), emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk '(NR == 1) || (FNR > 1)' $input > ${prefix}.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -1 | grep -o -E '([0-9]+.){1,2}[0-9]')
    END_VERSIONS
    """
}

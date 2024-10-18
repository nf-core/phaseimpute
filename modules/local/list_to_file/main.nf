process LIST_TO_FILE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(input, arity: '0..*'), val(id)

    output:
    tuple val(meta), path('*.id.txt'), path('*.noid.txt'), emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Take all files of the input and list them in a file
    # and add as second column the id
    awk 'BEGIN {
        split("${input}", f, " ");
        ids = "${id}";
        gsub(/[\\[\\]]/, "", ids);
        split(ids, i, ", ");
        for (j in f) print f[j], i[j]
    }' > ${prefix}.id.txt

    # Take all files of the input and list them in a file
    # without the id

    awk 'BEGIN {
        split("${input}", f, " ");
        for (j in f) print f[j]
    }' > ${prefix}.noid.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}

process ADD_COLUMNS {
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
    awk '(NR>=2) && (NR<=10)' $input | \\
    awk 'NR==1{\$(NF+1)="ID"} NR>1{\$(NF+1)="${meta.id}"}1' | \\
    awk 'NR==1{\$(NF+1)="Region"} NR>1{\$(NF+1)="${meta.region}"}1' | \\
    awk 'NR==1{\$(NF+1)="Depth"} NR>1{\$(NF+1)="${meta.depth}"}1' | \\
    awk 'NR==1{\$(NF+1)="GPArray"} NR>1{\$(NF+1)="${meta.gparray}"}1' | \\
    awk 'NR==1{\$(NF+1)="Tools"} NR>1{\$(NF+1)="${meta.tools}"}1' | \\
    awk 'NR==1{\$(NF+1)="Panel"} NR>1{\$(NF+1)="${meta.panel}"}1' > \\
    ${prefix}.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -1 | grep -o -E '([0-9]+.){1,2}[0-9]')
    END_VERSIONS
    """
}

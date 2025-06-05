process IMPUTE5_CONVERTREF {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/lindonkambule/impute5"

    input:
    tuple val(meta),
          path(ref_bcf),
          path(ref_index),
          val(region)

    output:
    tuple val(meta),
          path("*_xcf.bcf"), 
          path("*_xcf.bcf.csi"), 
          path("*_xcf.bin"),
          path("*_xcf.fam"),     emit: xcf_file
    path "versions.yml",       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maf = task.ext.maf ?: 0.03125  
    def chrom = region.split(':')[0] 

    """
    xcftools_static view \\
        --i ${ref_bcf} \\
        --o ${prefix}_xcf.bcf \\
        -O sh \\
        --r ${chrom} \\
        -m ${maf} \\
        -T ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xcftools: \$(xcftools_static --version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_xcf.bcf \\
          ${prefix}_xcf.bcf.csi \\
          ${prefix}_xcf.bin \\
          ${prefix}_xcf.fam
          
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xcftools: "stub"
    END_VERSIONS
    """
}
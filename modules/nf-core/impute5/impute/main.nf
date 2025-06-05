process IMPUTE5_IMPUTE {
    tag "$meta.id"
    label 'process_high'

    container "lindonkambule/impute5:v1.2.0"

    input:
    tuple val(meta),
          path(target_bcf),
          path(target_index),
          path(ref_xcf),
          path(ref_index),
          path(ref_bin),
          path(ref_fam),
          path(gen_map),
          val(region)

    output:
    tuple val(meta), path("*_imputed.bcf"), emit: bcf
    tuple val(meta), path("*_imputed.log"), emit: log
    path 'versions.yml',                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = "_imputed"
    def out_prefix = "${prefix}${suffix}"
    
    // Gestion sécurisée de la région
    def regionArg = ""
    def bufferArg = ""
    
    if (region && region != "" && region != "[]") {
        try {
            def (chr, coords) = region.split(':')
            def (start, end) = coords.split('-').collect{ it as Integer }
            def bufStart = Math.max(1, start - (task.ext.buffer_size ?: 250000))
            def bufEnd = end + (task.ext.buffer_size ?: 250000)
            def bufferRegion = "${chr}:${bufStart}-${bufEnd}"
            
            regionArg = "--r ${region}"
            bufferArg = "--buffer-region ${bufferRegion}"
        } catch (Exception e) {
            log.warn "Invalid region format: ${region}. Proceeding without region specification."
        }
    }
    
    """
    impute5_v1.2.0_static \\
        --g ${target_bcf} \\
        --h ${ref_xcf} \\
        --m ${gen_map} \\
        ${regionArg} \\
        ${bufferArg} \\
        --o ${out_prefix}.bcf \\
        --l ${out_prefix}.log \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        impute5: \$(impute5_v1.2.0_static --help | grep -E "^IMPUTE5" | head -n1 | sed 's/IMPUTE5 v//' | sed 's/ .*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = "_imputed"
    def out_prefix = "${prefix}${suffix}"
    """
    touch ${out_prefix}.bcf ${out_prefix}.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        impute5: "stub"
    END_VERSIONS
    """
}
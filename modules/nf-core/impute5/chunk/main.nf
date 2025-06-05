process IMPUTE5_CHUNK {
    tag "$meta.id"
    label 'process_medium'
    container "lindonkambule/impute5:v1.2.0"

    input:
    tuple val(meta),
          path(ref_xcf),
          path(ref_index),
          path(ref_bin),
          path(ref_fam),
          path(target_bcf),
          path(target_index),
          val(region)

    output:
    tuple val(meta),
          path("chunks_*.txt"), emit: chunks
    tuple val(meta),
          path(ref_xcf),     emit: ref_xcf
    tuple val(meta),
          path(ref_index),   emit: ref_index
    tuple val(meta),
          path(ref_bin),     emit: ref_bin
    tuple val(meta),
          path(ref_fam),     emit: ref_fam
    tuple val(meta),
          path(target_bcf),  emit: target
    tuple val(meta),
          path(target_index),emit: tbi
    path 'versions.yml',              emit: versions

    script:
    def args   = task.ext.args   ?: ''
    def chrom  = region.split(':')[0] 
    """
    imp5Chunker_v1.2.0_static \\
      --h ${ref_xcf} \\
      --g ${target_bcf} \\
      --r ${chrom} \\
      ${args} \\
      --o chunks_${chrom}.txt

    cat <<-END_VERSIONS > versions.yml
    \"${task.process}\":
      impute5_chunker: \$(imp5Chunker_v1.2.0_static --version)
    END_VERSIONS
    """

    stub:
    def chrom = region ? region.split(':')[0] : 'test'
    """
    # create empty outputs
    touch chunks_${chrom}.txt
    touch stub_ref.xcf
    touch stub_ref.xcf.csi
    touch stub_ref.bin
    touch stub_ref.fam
    touch stub_target.bcf
    touch stub_target.bcf.csi
    
    cat <<-END_VERSIONS > versions.yml
    \"${task.process}\":
      impute5_chunker: \"stub\"
    END_VERSIONS
    """
}
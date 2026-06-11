process VCFLIB_VCFFIXUP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/169e4e28f26469eb05baf60eab777bccadd747ac75038c6bb22149cd40c2ff38/data':
        'community.wave.seqera.io/library/bcftools_vcflib:0b47030679d1eff1' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.{vcf,bcf}{,.gz}"), emit: vcf
    tuple val("${task.process}"), val('vcflib'), val("1.0.14"), topic: versions, emit: versions_vcflib
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    def args2   = task.ext.args2  ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}.fixed"

    if ( "$vcf" == "${prefix}.vcf.gz" ) {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    """
    vcffixup \\
        $args \\
        $vcf \\
        | bcftools view -Ob -o ${prefix}.bcf
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}.fixed"

    if ( "$vcf" == "${prefix}.vcf.gz" ) {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    """
    echo | gzip > ${prefix}.vcf.gz
    """
}

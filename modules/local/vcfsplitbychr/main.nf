process VCFSPLITBYCHR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(index), val(chromosomes)  
    
    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.csi"), emit: vcf_per_chr
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chr_list = chromosomes.join(' ')  
    """
    for chr in ${chr_list}; do
        echo "Processing chromosome: \$chr"
        bcftools view \\
            --regions \$chr \\
            --output-type z \\
            $args \\
            --output ${prefix}_\${chr}.vcf.gz \\
            ${vcf}
        
        bcftools index --csi ${prefix}_\${chr}.vcf.gz
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_chr1.vcf.gz
    touch ${prefix}_chr1.vcf.gz.csi
    touch ${prefix}_chr2.vcf.gz
    touch ${prefix}_chr2.vcf.gz.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
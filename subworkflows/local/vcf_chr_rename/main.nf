include { BCFTOOLS_ANNOTATE           } from '../../../modules/nf-core/bcftools/annotate/main.nf'
include { BCFTOOLS_INDEX as VCF_INDEX } from '../../../modules/nf-core/bcftools/index/main.nf'

workflow VCF_CHR_RENAME {
    take:
    ch_vcf          // channel: [ [id, ref], vcf, csi ]
    file_chr_rename // file

    main:

    ch_versions = Channel.empty()

    // Rename the chromosome without prefix
    BCFTOOLS_ANNOTATE(ch_vcf
        .combine(Channel.of([[], [], []]))
        .combine(Channel.of(file_chr_rename))
    )
    
    VCF_INDEX(BCFTOOLS_ANNOTATE.out.vcf)
    ch_versions = ch_versions.mix(VCF_INDEX.out.versions.first())

    ch_vcf_rename = BCFTOOLS_ANNOTATE.out.vcf
        .combine(VCF_INDEX.out.csi)

    emit:
    vcf_rename     = ch_vcf_rename         // [ meta, vcf, csi ]
    versions       = ch_versions           // channel: [ versions.yml ]
}
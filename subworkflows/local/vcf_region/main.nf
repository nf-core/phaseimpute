include { BCFTOOLS_VIEW as VIEW_VCF_REGION } from '../../../modules/nf-core/bcftools/view/main.nf'
include { BCFTOOLS_INDEX                   } from '../../../modules/nf-core/bcftools/index/main.nf'


workflow VCF_REGION {
    take:
    ch_vcf          // channel: [ [id], vcf ]
    ch_region       // channel: [ [region], val(region) ]
    ch_fasta        // channel: [ fasta ]

    main:

    ch_versions = Channel.empty()

    // Filter the region of interest of the panel file
    ch_input_region = ch_vcf
        .combine(ch_fasta)
        .combine(ch_region)
        .map{ metaI, vcf, index, fasta, metaR, region ->
            [metaI + metaR, vcf, index, region+",chr"+region]}

    VIEW_VCF_REGION(ch_input_region, [], [], [])
    ch_versions = ch_versions.mix(VIEW_VCF_REGION.out.versions.first())

    BCFTOOLS_INDEX(VIEW_VCF_REGION.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

    ch_vcf_region = VIEW_VCF_REGION.out.vcf
        .combine(BCFTOOLS_INDEX.out.csi)

    emit:
    vcf_region    = ch_vcf_region   // channel: [ metaIR, vcf, index ]
    versions      = ch_versions     // channel: [ versions.yml ]

}

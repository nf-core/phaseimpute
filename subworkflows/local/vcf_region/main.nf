include { BCFTOOLS_VIEW as VIEW_VCF_REGION } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_INDEX                   } from '../../../modules/nf-core/bcftools/index'


workflow VCF_REGION {
    take:
    ch_vcf          // channel: [ [id], vcf ]
    ch_region       // channel: [ [chr, region], region ]

    main:

    ch_versions = Channel.empty()

    // Filter the region of interest of the vcf file
    ch_input_region = ch_vcf
        .combine(ch_region)
        .map{
            metaI, vcf, index, metaCR, region ->
            [metaI + metaCR, vcf, index, region+",chr"+region]
        }

    VIEW_VCF_REGION(ch_input_region, [], [], [])
    ch_versions = ch_versions.mix(VIEW_VCF_REGION.out.versions.first())

    BCFTOOLS_INDEX(VIEW_VCF_REGION.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

    ch_vcf_region = VIEW_VCF_REGION.out.vcf
        .combine(BCFTOOLS_INDEX.out.csi)

    emit:
    vcf_region    = ch_vcf_region   // channel: [ [id, chr, region], vcf, index ]
    versions      = ch_versions     // channel: [ versions.yml ]

}

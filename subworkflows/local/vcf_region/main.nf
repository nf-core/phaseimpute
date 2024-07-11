include { BCFTOOLS_VIEW } from '../../../modules/nf-core/bcftools/view'

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

    BCFTOOLS_VIEW(ch_input_region, [], [], [])
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())

    ch_vcf_region = BCFTOOLS_VIEW.out.vcf
        .combine(BCFTOOLS_VIEW.out.tbi)

    emit:
    vcf_region    = ch_vcf_region   // channel: [ [id, chr, region], vcf, tbi ]
    versions      = ch_versions     // channel: [ versions.yml ]

}

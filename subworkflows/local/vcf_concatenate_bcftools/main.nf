include { BCFTOOLS_CONCAT } from '../../../modules/nf-core/bcftools/concat'

workflow VCF_CONCATENATE_BCFTOOLS {

    take:
    ch_vcf_tbi // channel: [ [id, panel, chr, tools], vcf, tbi ]

    main:

    ch_versions = Channel.empty()

    // Keep only id from meta
    ch_vcf_tbi_grouped = ch_vcf_tbi
        .map{ metaIPTC, vcf, tbi -> [metaIPTC.subMap("id", "tools", "panel") + ["chr": "all"], vcf, tbi] }
        .groupTuple( by:0 )

    // Ligate and concatenate chunks
    BCFTOOLS_CONCAT(ch_vcf_tbi_grouped)
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions.first())

    // Join VCFs and TBIs
    ch_vcf_tbi_join = BCFTOOLS_CONCAT.out.vcf
        .join(BCFTOOLS_CONCAT.out.tbi)

    emit:
    vcf_tbi      = ch_vcf_tbi_join // channel: [ [id], vcf, tbi ]
    versions     = ch_versions     // channel: [ versions.yml ]
}

include { BCFTOOLS_CONCAT } from '../../../modules/nf-core/bcftools/concat'

workflow VCF_CONCATENATE_BCFTOOLS {

    take:
    ch_vcf_index // channel: [ [id, panel_id, chr, tools], vcf, index ]

    main:

    ch_versions = channel.empty()

    // Keep only id from meta
    ch_vcf_index_grouped = ch_vcf_index
        .map{ metaIPTC, vcf, index -> [metaIPTC.subMap("id", "tools", "panel_id", "batch") + ["chr": "all"], vcf, index] }
        .groupTuple( by:0 )
        .map{ metaIPTC, vcf, index -> [metaIPTC, vcf, index, vcf.size() ] } // Compute number of records
        .branch{ it ->
            one: it[3] == 1
            more: it[3] > 1
        }

    // Ligate and concatenate chunks
    BCFTOOLS_CONCAT(ch_vcf_index_grouped.more.map{ it -> [it[0], it[1], it[2]] })
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions.first())

    // Join VCFs and TBIs
    ch_vcf_index_concat = BCFTOOLS_CONCAT.out.vcf
        .join(
            BCFTOOLS_CONCAT.out.tbi.mix(
                BCFTOOLS_CONCAT.out.csi
            )
        )

    ch_vcf_index_join = ch_vcf_index_grouped.one
        .map{ it -> [it[0], it[1][0], it[2][0]] }
        .mix(ch_vcf_index_concat)

    emit:
    vcf_index    = ch_vcf_index_join // channel: [ [id], vcf, index ]
    versions     = ch_versions     // channel: [ versions.yml ]
}

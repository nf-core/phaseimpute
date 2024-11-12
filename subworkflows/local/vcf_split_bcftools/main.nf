include { BCFTOOLS_PLUGINSPLIT  } from '../../../modules/nf-core/bcftools/pluginsplit'
include { BCFTOOLS_QUERY        } from '../../../modules/nf-core/bcftools/query/main'

workflow VCF_SPLIT_BCFTOOLS {
    take:
    ch_vcf          // channel: [ [id, chr, tools], vcf, index ]

    main:

    ch_versions = Channel.empty()

    BCFTOOLS_QUERY(ch_vcf, [], [], []) // List samples

    ch_samples = ch_vcf
        .join(
            BCFTOOLS_QUERY.out.output.splitText().groupTuple()
                .ifEmpty{
                    error "VCF_SPLIT_BCFTOOLS: File does not have any samples information"
                }
        )
        .branch {
            one : it[3].size() == 1
            multiple : it[3].size() > 1
        }

    BCFTOOLS_PLUGINSPLIT(ch_samples.multiple.map{it[0..2]}, [], [], [], [])
    ch_versions = ch_versions.mix(BCFTOOLS_PLUGINSPLIT.out.versions.first())

    ch_vcf_samples = BCFTOOLS_PLUGINSPLIT.out.vcf
        .transpose()
        .map{metaITC, vcf -> [metaITC + [id: vcf.getBaseName().tokenize(".")[0]], vcf]}

    ch_tbi_samples = BCFTOOLS_PLUGINSPLIT.out.tbi
        .transpose()
        .map{metaITC, tbi -> [metaITC + [id: tbi.getBaseName().tokenize(".")[0]], tbi]}

    ch_vcf_tbi_samples = ch_samples.one
        .map{ it[0..2] }
        .mix(ch_vcf_samples
            .join(ch_tbi_samples)
        )

    emit:
    vcf_tbi        = ch_vcf_tbi_samples   // channel: [ [id, chr, tools], vcf, index ]
    versions       = ch_versions          // channel: [ versions.yml ]

}

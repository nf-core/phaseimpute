include { BCFTOOLS_PLUGINSPLIT  } from '../../../modules/nf-core/bcftools/pluginsplit'
include { BCFTOOLS_INDEX        } from '../../../modules/nf-core/bcftools/index'

workflow VCF_SAMPLES_BCFTOOLS {
    take:
    ch_vcf          // channel: [ [id, chr, tools], vcf ]

    main:

    ch_versions = Channel.empty()

    BCFTOOLS_PLUGINSPLIT(ch_vcf, [], [], [], [])
    ch_versions = ch_versions.mix(BCFTOOLS_PLUGINSPLIT.out.versions.first())

    ch_vcf_samples = BCFTOOLS_PLUGINSPLIT.out.vcf
        .transpose()
        .map{metaITC, vcf -> [metaITC + [id: vcf.getBaseName().tokenize(".")[0]], vcf]}

    BCFTOOLS_INDEX(ch_vcf_samples)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

    ch_vcf_tbi_samples = ch_vcf_samples
        .join(BCFTOOLS_INDEX.out.tbi)

    emit:
    vcf_tbi_join   = ch_vcf_tbi_samples   // channel: [ [id, chr, tools], vcf, index ]
    versions       = ch_versions          // channel: [ versions.yml ]

}

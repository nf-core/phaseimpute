include { BCFTOOLS_ANNOTATE           } from '../../../modules/nf-core/bcftools/annotate/main.nf'
include { BCFTOOLS_INDEX              } from '../../../modules/nf-core/bcftools/index/main.nf'
include { FAITOCHR                    } from '../../../modules/local/faitochr/main.nf'

workflow VCF_CHR_RENAME {
    take:
    ch_vcf          // channel: [ [id], vcf, index ]
    ch_fasta        // channel: [ [id], fasta, fai ]

    main:

    ch_versions = Channel.empty()

    // Generate the chromosome renaming file
    FAITOCHR(ch_fasta.map{ metaG, fasta, fai -> [metaG, fai] })
    ch_versions = ch_versions.mix(FAITOCHR.out.versions)

    // Rename the chromosome without prefix
    BCFTOOLS_ANNOTATE(
        ch_vcf // channel: [ [id], vcf, index ]
            .combine(Channel.of([[],[],[]]))
            .combine(FAITOCHR.out.annot_chr.map{it[1]})
    )
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())

    BCFTOOLS_INDEX(BCFTOOLS_ANNOTATE.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

    ch_vcf_renamed = BCFTOOLS_ANNOTATE.out.vcf
        .combine(BCFTOOLS_INDEX.out.csi, by:0)

    emit:
    vcf_renamed    = ch_vcf_renamed        // [ meta, vcf, csi ]
    versions       = ch_versions           // channel: [ versions.yml ]
}

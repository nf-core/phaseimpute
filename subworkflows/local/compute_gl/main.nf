include { BCFTOOLS_MPILEUP          } from '../../../modules/nf-core/bcftools/mpileup/main.nf'
include { BCFTOOLS_INDEX            } from '../../../modules/nf-core/bcftools/index/main.nf'


workflow COMPUTE_GL {

    take:
    ch_input   // channel: [ [id, ref], bam, bai ]
    ch_target  // channel: [ [panel], sites, tsv]
    ch_fasta   // channel: [ [ref], fasta, fai]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_mpileup = ch_input
        .combine(ch_target)
        .map{metaI, bam, bai, metaP, sites, tsv ->
                [metaI + metaP, bam, sites, tsv]}

    BCFTOOLS_MPILEUP(
        ch_mpileup,
        ch_fasta,
        false
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions.first())

    ch_output = BCFTOOLS_MPILEUP.out.vcf
        .combine(BCFTOOLS_MPILEUP.out.tbi, by:0)

    ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_MPILEUP.out.stats.map{ it[1] })

    emit:
    vcf           = ch_output                    // channel: [ [id, panel], vcf, tbi ]
    versions      = ch_versions                  // channel: [ versions.yml ]
    multiqc_files = ch_multiqc_files
}
